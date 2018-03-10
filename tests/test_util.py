#!/usr/bin/env python

import datetime
import os.path
import shutil
import tempfile

from nose.tools import raises
import numpy as np

import clawpack.geoclaw.util as util
from clawpack.geoclaw.util import NOAA_API_URL, fetch_noaa_tide_data


class TestFetchNoaaTideData:
    station = '1234567'
    begin_date = datetime.datetime(2000, 10, 30, 12, 0)
    end_date = datetime.datetime(2000, 10, 30, 12, 24)

    water_level_url = '{}?product={}&station={}'.format(
        NOAA_API_URL, 'water_level', station)
    predictions_url = '{}?product={}&station={}'.format(
        NOAA_API_URL, 'predictions', station)

    @classmethod
    def setup_class(cls):
        cls.temp_dir = tempfile.mkdtemp()

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.temp_dir)

    def test_retrieve_and_cache(self):
        cache_dir = os.path.join(self.temp_dir,
                                 self.test_retrieve_and_cache.__name__)

        water_level_response = \
            ('Date Time, Water Level, Sigma, O, F, R, L, Quality\n'
             '2000-10-30 12:00,1.001,0.001,0,0,0,0,v\n'
             '2000-10-30 12:06,1.002,0.002,0,0,0,0,v\n'
             '2000-10-30 12:12,1.003,0.003,0,0,0,0,v\n'
             '2000-10-30 12:18,1.004,0.004,0,0,0,0,v\n'
             '2000-10-30 12:24,1.005,0.005,0,0,0,0,v\n')

        predictions_response = ('Date Time, Prediction\n'
                                '2000-10-30 12:00,1.101\n'
                                '2000-10-30 12:06,1.102\n'
                                '2000-10-30 12:12,1.103\n'
                                '2000-10-30 12:18,1.104\n'
                                '2000-10-30 12:24,1.105\n')

        date_time_range = ['2000-10-30T12:00',
                           '2000-10-30T12:06',
                           '2000-10-30T12:12',
                           '2000-10-30T12:18',
                           '2000-10-30T12:24']

        date_time_expected = np.array(date_time_range, dtype='datetime64[m]')
        water_level_expected = np.array([1.001, 1.002, 1.003, 1.004, 1.005])
        prediction_expected = np.array([1.101, 1.102, 1.103, 1.104, 1.105])

        # monkey patch urllib to return mock data
        def mock_read_response(url):
            if 'product=water_level' in url:
                return water_level_response
            if 'product=predictions' in url:
                return predictions_response
            raise AssertionError
        self._monkey_patch_urlopen(mock_read_response)

        # first time, should fetch data and save in cache
        self._fetch_and_assert(date_time_expected, water_level_expected,
                               prediction_expected, cache_dir)

        # configure endpoint to not return data anymore
        self._monkey_patch_urlopen(lambda url: 'Should read from cache')

        # second time, should read data from cache
        self._fetch_and_assert(date_time_expected, water_level_expected,
                               prediction_expected, cache_dir)

    @raises(ValueError)
    def test_api_error(self):
        cache_dir = os.path.join(self.temp_dir, self.test_api_error.__name__)

        # configure endpoint to return an error
        self._monkey_patch_urlopen(lambda url: 'Something went wrong')

        # should raise ValueError
        fetch_noaa_tide_data(self.station, self.begin_date, self.end_date,
                             cache_dir=cache_dir)

    @raises(ValueError)
    def test_date_time_range_mismatch(self):
        cache_dir = os.path.join(self.temp_dir,
                                 self.test_date_time_range_mismatch.__name__)

        # missing first two entries
        water_level_response = \
            ('Date Time, Water Level, Sigma, O, F, R, L, Quality\n'
             '2000-10-30 12:12,1.003,0.003,0,0,0,0,v\n'
             '2000-10-30 12:18,1.004,0.004,0,0,0,0,v\n'
             '2000-10-30 12:24,1.005,0.005,0,0,0,0,v\n')

        # missing last two entries
        predictions_response = ('Date Time, Prediction\n'
                                '2000-10-30 12:00,1.101\n'
                                '2000-10-30 12:06,1.102\n'
                                '2000-10-30 12:12,1.103\n')

        # monkey patch urllib to return mock data
        def mock_read_response(url):
            if 'product=water_level' in url:
                return water_level_response
            if 'product=predictions' in url:
                return predictions_response
            raise AssertionError
        self._monkey_patch_urlopen(mock_read_response)

        # should raise ValueError
        fetch_noaa_tide_data(self.station, self.begin_date, self.end_date,
                             cache_dir=cache_dir)

    def test_missing_values(self):
        cache_dir = os.path.join(self.temp_dir,
                                 self.test_missing_values.__name__)

        # missing fist two water level values
        water_level_response = \
            ('Date Time, Water Level, Sigma, O, F, R, L, Quality\n'
             '2000-10-30 12:00,,0.001,0,0,0,0,v\n'
             '2000-10-30 12:06,,0.002,0,0,0,0,v\n'
             '2000-10-30 12:12,1.003,0.003,0,0,0,0,v\n'
             '2000-10-30 12:18,1.004,0.004,0,0,0,0,v\n'
             '2000-10-30 12:24,1.005,0.005,0,0,0,0,v\n')

        # missing a prediction value
        predictions_response = ('Date Time, Prediction\n'
                                '2000-10-30 12:00,1.101\n'
                                '2000-10-30 12:06,1.102\n'
                                '2000-10-30 12:12,1.103\n'
                                '2000-10-30 12:18,\n'
                                '2000-10-30 12:24,1.105\n')

        date_time_range = ['2000-10-30T12:00',
                           '2000-10-30T12:06',
                           '2000-10-30T12:12',
                           '2000-10-30T12:18',
                           '2000-10-30T12:24']

        date_time_expected = np.array(date_time_range, dtype='datetime64[m]')
        water_level_expected = np.array([np.nan, np.nan, 1.003, 1.004, 1.005])
        prediction_expected = np.array([1.101, 1.102, 1.103, np.nan, 1.105])

        # monkey patch urllib to return mock data
        def mock_read_response(url):
            if 'product=water_level' in url:
                return water_level_response
            if 'product=predictions' in url:
                return predictions_response
            raise AssertionError
        self._monkey_patch_urlopen(mock_read_response)

        self._fetch_and_assert(date_time_expected, water_level_expected,
                               prediction_expected, cache_dir)

    def _fetch_and_assert(self, date_time_expected, water_level_expected,
                          prediction_expected, cache_dir):
        # fetch data
        date_time, water_level, prediction = fetch_noaa_tide_data(
            self.station, self.begin_date, self.end_date, cache_dir=cache_dir)

        # make sure data was correctly retrieved
        np.testing.assert_equal(date_time, date_time_expected)
        np.testing.assert_equal(water_level, water_level_expected)
        np.testing.assert_equal(prediction, prediction_expected)

    @staticmethod
    def _monkey_patch_urlopen(mock_read_response):
        """Modify 'util.urlopen' to return a custom response."""
        def mock_urlopen(url):
            class MockHttpResponse:
                def __enter__(self):
                    return self

                def __exit__(self, exc_type, exc_value, traceback):
                    pass

                def read(self):
                    # encode mock response
                    text = mock_read_response(url)
                    return text.encode('utf-8')

            return MockHttpResponse()

        # replace 'util.urlopen' with mock implementation
        util.urlopen = mock_urlopen

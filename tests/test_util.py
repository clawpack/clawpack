#!/usr/bin/env python

import datetime

from nose import SkipTest
from nose.tools import raises
import numpy as np

try:
    import mockfs
    import requests
    import requests_mock
except ImportError:
    raise SkipTest('Required libraries not available')

from clawpack.geoclaw.util import NOAA_API_URL, fetch_noaa_tide_data


class TestFetchNoaaTideData:
    station = '1234567'
    begin_date = datetime.datetime(2000, 10, 30, 12, 0)
    end_date = datetime.datetime(2000, 10, 30, 12, 24)

    water_level_url = '{}?product={}&station={}'.format(
        NOAA_API_URL, 'water_level', station)
    predictions_url = '{}?product={}&station={}'.format(
        NOAA_API_URL, 'predictions', station)

    def setup(self):
        mockfs.replace_builtins()

    def teardown(self):
        mockfs.restore_builtins()

    @requests_mock.Mocker()
    def test_retrieve_and_cache(self, mock):
        water_level_response = \
            ('Date Time, Water Level, Sigma, O, F, R, L, Quality\n'
             '2000-10-30 12:00,1.001,0.001,0,0,0,0,v\n'
             '2000-10-30 12:06,1.002,0.002,0,0,0,0,v\n'
             '2000-10-30 12:12,1.003,0.003,0,0,0,0,v\n'
             '2000-10-30 12:18,1.004,0.004,0,0,0,0,v\n'
             '2000-10-30 12:24,1.005,0.005,0,0,0,0,v\n')

        predictions_response = \
            ('Date Time, Prediction\n'
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

        # rewire requests to return mock data
        mock.get(self.water_level_url, text=water_level_response)
        mock.get(self.predictions_url, text=predictions_response)

        # first time, should fetch data and save in cache
        self._fetch_and_assert(date_time_expected, water_level_expected,
                               prediction_expected)

        # configure endpoint to not return data anymore
        mock.get(self.water_level_url, text='Should read from cache')
        mock.get(self.predictions_url, text='Should read from cache')

        # second time, should read data from cache
        self._fetch_and_assert(date_time_expected, water_level_expected,
                               prediction_expected)

    @raises(ValueError)
    @requests_mock.Mocker()
    def test_api_error(self, mock):
        # configure endpoint to return an error
        mock.get(NOAA_API_URL, text='Something went wrong')

        # should raise ValueError
        fetch_noaa_tide_data(self.station, self.begin_date, self.end_date)

    @raises(ValueError)
    @requests_mock.Mocker()
    def test_date_time_range_mismatch(self, mock):
        # missing first two entries
        water_level_response = \
            ('Date Time, Water Level, Sigma, O, F, R, L, Quality\n'
             '2000-10-30 12:12,1.003,0.003,0,0,0,0,v\n'
             '2000-10-30 12:18,1.004,0.004,0,0,0,0,v\n'
             '2000-10-30 12:24,1.005,0.005,0,0,0,0,v\n')

        # missing last two entries
        predictions_response = \
            ('Date Time, Prediction\n'
             '2000-10-30 12:00,1.101\n'
             '2000-10-30 12:06,1.102\n'
             '2000-10-30 12:12,1.103\n')

        # rewire requests to return mock data
        mock.get(self.water_level_url, text=water_level_response)
        mock.get(self.predictions_url, text=predictions_response)

        # should raise ValueError
        fetch_noaa_tide_data(self.station, self.begin_date, self.end_date)

    @requests_mock.Mocker()
    def test_missing_values(self, mock):
        # missing fist two water level values
        water_level_response = \
            ('Date Time, Water Level, Sigma, O, F, R, L, Quality\n'
             '2000-10-30 12:00,,0.001,0,0,0,0,v\n'
             '2000-10-30 12:06,,0.002,0,0,0,0,v\n'
             '2000-10-30 12:12,1.003,0.003,0,0,0,0,v\n'
             '2000-10-30 12:18,1.004,0.004,0,0,0,0,v\n'
             '2000-10-30 12:24,1.005,0.005,0,0,0,0,v\n')

        # missing a prediction value
        predictions_response = \
            ('Date Time, Prediction\n'
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

        # rewire requests to return mock data
        mock.get(self.water_level_url, text=water_level_response)
        mock.get(self.predictions_url, text=predictions_response)

        self._fetch_and_assert(date_time_expected, water_level_expected,
                               prediction_expected)

    def _fetch_and_assert(self, date_time_expected, water_level_expected,
                          prediction_expected):
        # fetch data
        date_time, water_level, prediction = fetch_noaa_tide_data(
            self.station, self.begin_date, self.end_date)

        # make sure data was correctly retrieved
        np.testing.assert_equal(date_time, date_time_expected)
        np.testing.assert_equal(water_level, water_level_expected)
        np.testing.assert_equal(prediction, prediction_expected)

#!/usr/bin/env python

r"""Utility functions for GeoClaw"""

import os
import urllib2
import tarfile

def strip_archive_extensions(path, extensions=["tar", "tgz", "bz2", "gz"]):
    r"""
    Strip off archive extensions defined in *extensions* list.

    Return stripped path calling this function recursively until all splitext
    does not provide an extension in the *extensions* list.

    """

    if os.path.splitext(path)[-1][1:] in extensions:
        return strip_archive_extensions(os.path.splitext(path)[0])
    else:
        return path


def get_remote_file(url, output_dir=None, file_name=None, force=False,  
                         verbose=False, ask_user=False):
    r"""Fetch file located at *url* and store at *output_dir*.

    :Input:
    
     - *url* (path) - URL to file to be downloaded.
     - *output_dir* (path) - Directory that the remote file will be downloaded
       to.  Defaults to the current working directory returned by *os.getcwd()*.
     - *file_name* (string) - Name of local file.  This defaults to the name of
       the remote file.
     - *force* (bool) - Force downloading of remote file regardless of whether
       it exists locally or not.  Default is *False*
     - *verbose* (bool) - Print out status information.  Default is *False*
     - *ask_user* (bool) - Whether to ask the user if it is ok to download the
       file before proceeding.  Default is *False*

    :Raises:
     
    Exceptions are raised from the *urllib2* module having to do with errors
    fetching the remote file.  Please see its documentation for more details of
    the exceptions that can be raised.

    returns nothing
    """

    if output_dir is None:
        output_dir = os.getcwd()

    if file_name is None:
        file_name = os.path.basename(url)
        
    output_path = os.path.join(output_dir, file_name)
    unarchived_output_path = strip_archive_extensions(output_path)

    if not os.path.exists(unarchived_output_path) or force:

        if ask_user:
            ans = raw_input("  Ok to download topo file and save as %s?  \n"
                            % unarchived_output_path,
                            "     Type y[es], n[o].")
            if ans.lower() in ['y', 'yes']:
                if verbose:
                    print "*** Aborting download."
                return None
            
        if not os.path.exists(output_path):
            # Fetch remote file, will raise a variety of exceptions depending on
            # the retrieval problem if it happens
            if verbose:
                print "Downloading %s to %s..." % (url, output_path)
            with open(output_path, "w") as output_file:
                remote_file = urllib2.urlopen(url)
                output_file.write(remote_file.read())
            if verbose:
                print "Done downloading."

        if tarfile.is_tarfile(output_path):
            if verbose:
                print "Un-archiving %s to %s..." % (output_path, 
                                                    unarchived_output_path)
            with tarfile.open(output_path, mode="r:*") as tar_file:
                tar_file.extractall(path=output_dir)
            if verbose:
                print "Done un-archiving."
        # TODO: Should check here if a file is a bare compressed file (no tar)
    else:
        if verbose:
            print "Skipping %s because it already exists locally." % url
        return None

    return unarchived_output_path
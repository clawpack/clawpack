#!/usr/bin/python
"""
iotools
=========
    Provides several useful functions for easy input and output from data files

    Contains:
        convertd2e
        datafile2array
        array2datafile


"""

from __future__ import absolute_import
import string
import re
import numpy
from six.moves import range


#====================================================================
def convertd2e (numberstring=" "):
    """
    def convertd2e (numberstring):

        takes a string and replaces all d's or D's with e's.
        usefule for reading data files with doubles output by 
        Fortran, which uses d, or D instead of e, for the exponent.
    """
    Dd=re.compile("[Dd]")
    newstring=Dd.sub("e",numberstring)
    return newstring
    # end convertd2e ==================================================


#======================================================================
def datafile2array (datafile=" ",sep=None, dtype="float",skiplines=0, \
    skipfirstcols=0, skiplastcols=0):
    """
        open and read data from a ascii text file into a
        numpy array. The number of rows and columns in the data file will match
        the size of the array.

        sep is the character seperation between fields in the data file
        dtype is the how the file data is to be interpretted.
        skiplines = n skips the first n line(s). 
        skipfirstcols = n skips the first n column(s)
        skiplastcols = n skips the last n column(s).
    """
    fid=open(datafile)
    data=fid.readlines()
    fid.close()

    dataarray=[]
    for row in range(skiplines,len(data)):
        data[row]=convertd2e(data[row])
        data[row]=string.split(data[row],sep)
        if data[row]!=[]:
            if dtype!=" ":
                for col in range(skipfirstcols,len(data[row])-skiplastcols) :
                    if dtype=="float":
                        data[row][col]=float(data[row][col])
                    elif dtype=="int":
                        data[row][col]=int(data[row][col])
            if dataarray!=[]: 
                if len(data[row])-skipfirstcols-skiplastcols==len(dataarray[0]):
                    dataarray.append(data[row][skipfirstcols:len(data[row])-skiplastcols])
            else:
                dataarray.append(data[row][skipfirstcols:len(data[row])-skiplastcols])
    
    dataarray=numpy.array(dataarray)
    return dataarray
    # end loaddatafile ==================================================

def array2datafile (dataarray,datafile=" ",sep=""):
    """
    def array2datafile (array,datafile=" "):

        output a tuple, list or numpy array into a ascii text file. 
        The number of rows and columns in the data file will match
        the size of the array.
   
    """
    fid=open(datafile,'w')
    shp=numpy.shape(dataarray)
    if len(shp)>1:
        for row in range(shp[0]):
            for col in range(shp[1]):
                fid.write("%s %s" % (dataarray[row][col],sep))

            fid.write("\n")
    else:
        for col in range(shp[0]) :
            fid.write("%s %s" % (dataarray[col],sep))

    fid.close()
    return
    # end array2datafile =====================================================
    


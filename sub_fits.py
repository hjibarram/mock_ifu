#!/usr/bin/python
import sys
import numpy as np
from numpy import float_
from numpy import absolute as abs
from numpy import random as ran
import cosmolopy.distance as cd
import time
#import my_auto_ssp_elines_rnd as ssp
import pyfits as pyf
import cosmolopy.distance as cd
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wfit
from scipy.interpolate.interpolate import interp1d
#import my as my
import os.path as ptt
import matplotlib

def sycall(comand):
    from subprocess import call
    line=comand.split(" ")
    fcomand=[]
    fcomand.extend(line)
    call(fcomand)

def wfits(name, data, hdr):
    if ptt.exists(name) == False:
        wfit(name,data,hdr)
    else:
        sycall("rm "+name)
        wfit(name,data,hdr)
        
sys.argv=filter(None,sys.argv)
file1=sys.argv[1]
file2=sys.argv[2]
file3=sys.argv[3]

[data1, head1]=gdata(file1, 0, header=True)
[data2, head2]=gdata(file2, 0, header=True)

data0=data1-data2
crpix1=head1["CRPIX1"]
crval1=head1["CRVAL1"]
cd1_1=head1["CD1_1"]
cd1_2=head1["CD1_2"]
crpix2=head1["CRPIX2"]
crval2=head1["CRVAL2"]
cd2_1=head1["CD2_1"]
cd2_2=head1["CD2_2"]

hd=pyf.PrimaryHDU().header
hd["CRPIX1"]=crpix1
hd["CRVAL1"]=crval1
hd["CD1_1"]=cd1_1
hd["CD1_2"]=cd1_2
hd["CRPIX2"]=crpix2
hd["CRVAL2"]=crval2
hd["CD2_1"]=cd2_1
hd["CD2_2"]=cd2_2
hd["NAXIS"]=head1["NAXIS"]
hd["NAXIS1"]=head1["NAXIS1"]
hd["NAXIS2"]=head1["NAXIS2"]
wfits(file3, data0, hd)
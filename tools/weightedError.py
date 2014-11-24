#!/usr/bin/python3
############################################################################################
#
# Filename : "weightedError.py"
# Author   : A. Ngo
# 
# 2014
#
############################################################################################
#
#
# This Python script is used to plot the weighted error distribution as a normed histogram and 
# compare it against the standard normal distribution, in order to show the unbiasedness of 
# the estimated field.
# Three HDF5-files are required in order to generate this plot:
# 1.) Yfield.h5
# 2.) Y_old.h5
# 3.) Estimation_Variance.h5
#
# USAGE:
# python3 weightedError.py <Path to the three HDF5-files> <zf>
#
# Meaningful values of the zoom factor zf lie between 0 anf 0.32.
# The PNG output can be found inside the sub-directory 'histo/'.
############################################################################################

import os
import sys
from os.path import basename
import numpy as np
#import tables as pt
import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from mayavi import mlab

def loadHDF5( hdf5file, hdf5group ):
    arr = []
    hdf5filebase = hdf5file.rsplit(".",1)[0]
    basename = os.path.basename( hdf5filebase )
    # Load HDF5 data to array
    a = h5py.File( hdf5file )
    dataset = a[ hdf5group ]
    arr = dataset[...]
    a.close()
    # Some statistics about the data:
    aMax = arr.max()
    aMin = arr.min()
    # print( basename + ": min=" + str(aMin) + ": max=" + str(aMax) )
    Nx=0
    Ny=0
    Nz=0
    if( arr.ndim == 3 ):
        Nx=arr.shape[2]
        Ny=arr.shape[1]
        Nz=arr.shape[0]
        # print( basename + ": 3D: " + str(Nx) + " x " + str(Ny) + " x " + str(Nz) )
    else:
        Nx=arr.shape[1]
        Ny=arr.shape[0]
        # print( basename + ": 2D: " + str(Nx) + " x " + str(Ny) )
    return arr
 

def loadMeasureZone( hdf5file, hdf5group, fraction ):
    arr = []
    zone = []
    hdf5filebase = hdf5file.rsplit(".",1)[0]
    basename = os.path.basename( hdf5filebase )
    # Load HDF5 data to array
    a = h5py.File( hdf5file )
    dataset = a[ hdf5group ]
    arr = dataset[...]
    a.close()
    # Some statistics about the data:
    aMax = arr.max()
    aMin = arr.min()
    # print( basename + ": min=" + str(aMin) + ": max=" + str(aMax) )
    Nx=0
    Ny=0
    Nz=0
    if( arr.ndim == 3 ):
        Nx=arr.shape[2]
        Ny=arr.shape[1]
        Nz=arr.shape[0]
        # print( basename + ": 3D: " + str(Nx) + " x " + str(Ny) + " x " + str(Nz) )
    else:
        Nx=arr.shape[1]
        Ny=arr.shape[0]
        # print( basename + ": 2D: " + str(Nx) + " x " + str(Ny) )
    #
    NxStart = int(fraction*Nx)
    NxEnd = Nx - NxStart
    NyStart = int(fraction*Ny)
    NyEnd = Ny - NyStart
    if( arr.ndim == 3 ):
        zone = arr[:,NyStart:NyEnd,NxStart:NxEnd]
    else:
        zone = arr[NyStart:NyEnd,NxStart:NxEnd]
    return zone



if __name__ == '__main__' :

    matplotlib.rcParams.update({'font.size': 18})
    matplotlib.rcParams.update({'font.family': 'serif'})

    sPathToHDF5 = sys.argv[1]
    xStart = float(sys.argv[2])

    Y_orig = []
    Y_est = []
    V_est = []

    if(xStart == 0):
        Y_orig = loadHDF5( sPathToHDF5 + "/Yfield.h5", "YField" )
    else:
        Y_orig = loadMeasureZone( sPathToHDF5 + "/Yfield.h5", "YField", xStart )
    Y_orig_mean = np.average( Y_orig )
    Y_orig_std = np.std( Y_orig )
    Y_orig_var = np.var( Y_orig )
    print( "Y_orig: (mean,sigma,variance) = " + str(Y_orig_mean) + ", " + str(Y_orig_std) + ", " + str(Y_orig_var) )

    if(xStart == 0):
        Y_est = loadHDF5( sPathToHDF5 + "/Y_old.h5", "Y_old" )
    else:
        Y_est = loadMeasureZone( sPathToHDF5 + "/Y_old.h5", "Y_old", xStart )
    Y_est_mean = np.average( Y_est )
    Y_est_std = np.std( Y_est )
    Y_est_var = np.var( Y_est )
    print( "Y_est : (mean,sigma,variance) " + str(Y_est_mean) + ", " + str(Y_est_std) + ", " + str(Y_est_var) )

    if(xStart == 0):
        V_est = loadHDF5( sPathToHDF5 + "/Estimation_Variance.h5", "sigma2" )
    else:
        V_est = loadMeasureZone( sPathToHDF5 + "/Estimation_Variance.h5", "sigma2", xStart )


    if not os.path.exists('histo'):
        os.makedirs('histo')

    fig = plt.figure(figsize=(8,8))
    plt.axis('off')
    im0 = []
    if( V_est.ndim == 3 ):
        im0 = plt.contourf( V_est[int(0.5*V_est.shape[0]),:,:], extend = 'both', cmap = "jet")
    else:
        im0 = plt.contourf( V_est, extend = 'both', cmap = "jet")
    cb = plt.colorbar(im0,orientation='vertical',pad=0.05,shrink=0.99)
    pngfile = "histo/V_est.png"
    plt.savefig( pngfile )
    os.system("convert -trim "+pngfile+" "+pngfile)
    plt.clf()



    W_err = (Y_orig - Y_orig_mean - Y_est + Y_est_mean) / V_est
    #Nx=W_err.shape[1]
    #Ny=W_err.shape[0]
    #print( "2D: " + str(Nx) + " x " + str(Ny) )

#    plt.axis('off')
#    if( W_err.ndim == 3 ):
#        im0 = plt.contourf( W_err[int(0.5*V_est.shape[0]),:,:], extend = 'both', cmap = "jet")
#    else:
#        im0 = plt.contourf( W_err, extend = 'both', cmap = "jet")
#    cb = plt.colorbar(im0,orientation='vertical',pad=0.05,shrink=0.99)
#    pngfile = "histo/wErr.png"
#    plt.savefig( pngfile )
#    os.system("convert -trim "+pngfile+" "+pngfile)
#    plt.clf()

    hist, bins = np.histogram( W_err, bins=100, normed=True )
    width = 0.6 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    plt.axis([-5, 5, -0.01, 0.41])
    # frame1 = plt.gca()
    # frame1.axes.get_xaxis().set_visible(False)
    # frame1.axes.get_yaxis().set_visible(False)
    #plt.axis('off')
    plt.bar(center, hist, align='center', width=width)
    t = np.arange(-5,5,0.05)
    plt.plot(t,1/np.sqrt(2*np.pi)*np.exp(-0.5*t*t),linewidth='2',color='red')
    pngfile = "histo/hist.png"
    plt.savefig( pngfile )
    os.system("convert -trim "+pngfile+" "+pngfile)
    plt.clf()

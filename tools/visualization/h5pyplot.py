#!/usr/bin/python3
#################################################################################################
#
# Filename: "h5pyplot.py"
# Authors : A. Ngo and O. Klein
# 
# 2014
#
#################################################################################################
#
# This Python script can be used to create 2-D plots from a 2-D array strored in HDF5
# or 2-D cross-sectional plots from a 3-D array strored in HDF5.
# Examples of such HDF5 files (.h5) can be found in the DATA and BUFFER_DG or BUFFER_FEM
# subdirectories created by the dune-gesis programs. These files store the parameter field
# or the finite element solutions on a structured mesh with rectangular or hexahedral elements.
#
#
# Usage:
# python3 h5pyplot.py <file> <group> <dim> <min> <max> <L> <zf> <plot style> <color style>
#
# file  : name of the HDF5 file
# group : name of the dataset inside the HDF5 file
# dim   : dimension of the multi-dim. array, the only options are dim=2 or dim=3
# min   : minimal value of the colorbar
# max   : maximal value of the colorbar
# L     : colorbar divided into L discrete levels, meaningful values are 16 or 32
# zf    : zoom factor, meaningful values are within the interval [0.0..0.32] with 0.0 meaning the full view without zoom
# ps    : plot style, the only options are C (contour) or some other character which is not C
# cs    : color style, the only options are J (Jet scheme) or some other character which activates the customizable colormap 'cmapRdYlBu'.
#
# If maxValue > minValue, the colorbar is scaled over the specified range [minValue..maxValue].
# If minValue > maxValue, the colorbar is scaled over the full range of ALL values of the HDF5 array.
#
# Examples:
# python3 h5pyplot.py Yfield.h5 YField 3 -7.2 -3.4 32 0.2 C J
# python3 h5pyplot.py Y_estimated.h5 Y_est 2 1 0 32 0.2 C J
#
#
#
#################################################################################################

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

def cmapRdYlBu (name = 'ownMap', N = 18) :
    #colList = ['#000000','#313695','#4575b4','#74add1','#abd9e9','#e0f3f8','#fee090','#fdae61','#f46d43','#d73027','#a50026','#FFFFFF']
    colList = ['#303090','#313695','#4070b4','#4575b4','#7075b4','#74add1','#abd9e9','#dbe9e9','#f0f0f0','#f6ee91','#f9d090','#fdae61','#f48d43','#f46d43','#d73027','#d93027','#a50026','#a70026']
    #colList = ['#67001f','#b2182b','#d6604d','#e4a582','#eddbc7','#efdfcf','#f7f7f7','#e1c5f0','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061']
    #colList = colList[::-1]
    cmap = col.LinearSegmentedColormap.from_list(name,colList,N)
    #cmap.set_under(color = '#000000')
    #cmap.set_over(color = '#cc8888')
    cmap.set_under(color = '#3030F0')
    cmap.set_over(color = '#ff1493')
    #return cm.jet
    return cmap

def plot2d (hdf5file, hdf5group, minVal, maxVal, steps, xStart, style, color) :
    hdf5filebase = hdf5file.rsplit(".",1)[0]
    basename = os.path.basename( hdf5filebase )

    bUseJet = False;
    bContourF = False;
    if( style == "C" ):
        bContourF = True;
    if( color == "J" ):
        bUseJet = True;

    a = h5py.File( hdf5file )
    dataset = a[ hdf5group ]
    arr = dataset[...]
    a.close()

    aMax = arr.max()
    aMin = arr.min()
    print( basename + ": min=" + str(aMin) + ": max=" + str(aMax) )

    arr_mean = np.average( arr )
    arr_std = np.std( arr )
    arr_var = np.var( arr )
    print( basename + ": (mean,sigma,var) = " + str(arr_mean) + ", " + str(arr_std) + ", " + str(arr_var) )    

    if( maxVal < minVal ):
        maxVal = aMax
        minVal = aMin

    Nx=arr.shape[1]
    Ny=arr.shape[0]

    NxStart = int(xStart*Nx) #- 3
    NxEnd = Nx - NxStart #- 6
    print("NxStart="+str(NxStart))
    print("NxEnd="+str(NxEnd))

    NyStart = int(xStart*Ny) #- 3
    NyEnd = Ny - NyStart #- 6
    print("NyStart="+str(NyStart))
    print("NyEnd="+str(NyEnd))

    arr2oben = arr[NyStart:NyEnd,NxStart:NxEnd]
    #arr2oben = arr2oben[::-1,::]
    print(arr2oben.shape)

    cmap = cmapRdYlBu('ownMap',steps)
    cm.register_cmap(cmap=cmap)

    levels = [minVal + i * (maxVal - minVal)/steps for i in range(steps+1)]

    fig = plt.figure(figsize=(8,8))
    plt.axis('off')

    if bContourF == True :
        if bUseJet == False:
            im0 = plt.contourf(arr2oben, extend = 'both', levels = levels, cmap = 'ownMap')
        else:
            im0 = plt.contourf(arr2oben, extend = 'both', levels = levels, cmap = "jet")
    else :
        im0 = plt.imshow(arr2oben, cmap = "jet", vmax=maxVal, vmin=minVal)


    if not os.path.exists('PNG'):
        os.makedirs('PNG')

    if not os.path.exists('PNG/Top'):
        os.makedirs('PNG/Top')
    pngfile = "PNG/Top/" + basename + ".png"
    plt.savefig( pngfile )
    os.system("convert -trim "+pngfile+" "+pngfile)
    plt.clf()
    
    fig = plt.figure(figsize=(8,8))
    plt.axis('off')
    cb = plt.colorbar(im0,orientation='vertical',pad=0.05,shrink=0.99,ticks=[aMin,aMax,minVal,maxVal])
    labels = [minVal + i * (maxVal - minVal)/4 for i in range(5)]
    cb.set_ticks( labels )
    cb.ax.tick_params(labelsize=20)
    if maxVal > 999:
        ticklabels = ["${:0.2}$".format(d) for d in labels]
    else:
        ticklabels = ["${:0.2f}$".format(d) for d in labels]
    cb.set_ticklabels( ticklabels )
    if not os.path.exists('PNG/Legend'):
        os.makedirs('PNG/Legend')
    pngfile = "PNG/Legend/" + basename + ".png"
    plt.savefig( pngfile )
    os.system("convert -trim "+pngfile+" "+pngfile)
    plt.clf()

def plot3d (hdf5file, hdf5group, minVal, maxVal, steps, xStart, style, color ) :
    print("")
    print("Plotting " + str(hdf5file))
    hdf5filebase = hdf5file.rsplit(".",1)[0]
    basename = os.path.basename( hdf5filebase )

    bUseJet = False;
    bContourF = False;
    if( style == "C" ):
        bContourF = True;
    if( color == "J" ):
        bUseJet = True;

#    a = pt.openFile( hdf5file )
#    a.root
#    arrObject = a.getNode("/",hdf5group)
#    arr = arrObject.read()
#    a.close()

    a = h5py.File( hdf5file )
    dataset = a[ hdf5group ]
    arr = dataset[...]
    a.close()

    aMax = arr.max()
    aMin = arr.min()
    print( basename + ": min=" + str(aMin) + ": max=" + str(aMax) )

    arr_mean = np.average( arr )
    arr_std = np.std( arr )
    arr_var = np.var( arr )
    print( basename + ": (mean,sigma,var) = " + str(arr_mean) + ", " + str(arr_std) + ", " + str(arr_var) )    

    if( maxVal < minVal ):
        maxVal = aMax
        minVal = aMin

    Nx=arr.shape[2]
    Ny=arr.shape[1]
    Nz=arr.shape[0]


    mx=int(Nx*25.0/50.0)
    my=int(Ny*25.0/50.0)
    mz=int(Nz*2.5/5)

    NxStart = int(xStart*Nx) #- 3
    NxEnd = Nx - NxStart #- 6
    print("NxStart="+str(NxStart))
    print("NxEnd="+str(NxEnd))

    NyStart = int(xStart*Ny) #- 3
    NyEnd = Ny - NyStart #- 6
    print("NyStart="+str(NyStart))
    print("NyEnd="+str(NyEnd))

    arr2oben = arr[mz,NyStart:NyEnd,NxStart:NxEnd]
    arr2oben = arr2oben[::-1,::]
    print(arr2oben.shape)
    arr2seite = arr[::-1,my,NxStart:NxEnd]
    print(arr2seite.shape)
    arr2vorne = arr[::-1,NyStart:NyEnd,mx]
    arr2vorne = arr2vorne[::,::-1]
    print(arr2vorne.shape)
    arr2vorne = arr2vorne.transpose()

    if bContourF == True :
        arr2oben = arr2oben[::-1,::]
        arr2seite = arr2seite[::-1,::]
        arr2vorne = arr2vorne[::-1,::]

    cmap = cmapRdYlBu('ownMap',steps)
    cm.register_cmap(cmap=cmap)

    levels = [minVal + i * (maxVal - minVal)/steps for i in range(steps+1)]
    
    fig = plt.figure(figsize=(8,8))
    plt.axis('off')

    im0 = 0;
    if bContourF == True :
        if bUseJet == False:
            im0 = plt.contourf(arr2oben, extend = 'both', levels = levels, cmap = 'ownMap')
        else:
            im0 = plt.contourf(arr2oben, extend = 'both', levels = levels, cmap = "jet")
    else :
        im0 = plt.imshow(arr2oben, cmap = "jet", vmax=maxVal, vmin=minVal)
    #plt.contour(arr2oben, linewidths=2, colors='w')

    if not os.path.exists('PNG'):
        os.makedirs('PNG')

    if not os.path.exists('PNG/Top'):
        os.makedirs('PNG/Top')
    pngfile = "PNG/Top/" + basename + ".png"
    plt.savefig( pngfile )
    os.system("convert -trim "+pngfile+" "+pngfile)
    plt.clf()
    
    fig = plt.figure(figsize=(8,4))
    plt.axis('off')
    if bContourF == True :
        if bUseJet == False:
            plt.contourf(arr2seite, extend = 'both', levels = levels, cmap = 'ownMap')
        else:
            plt.contourf(arr2seite, extend = 'both', levels = levels, cmap = "jet")
    else :
        plt.imshow(arr2seite, cmap = "jet", vmax=maxVal, vmin=minVal)
    if not os.path.exists('PNG/Side'):
        os.makedirs('PNG/Side')
    pngfile = "PNG/Side/" + basename + ".png"
    plt.savefig( pngfile )
    os.system("convert -trim "+pngfile+" "+pngfile)
    plt.clf()

    fig = plt.figure(figsize=(4,8))
    plt.axis('off')
    if bContourF == True :
        if bUseJet == False:
            plt.contourf(arr2vorne, extend = 'both', levels = levels, cmap = 'ownMap')
        else:
            plt.contourf(arr2vorne, extend = 'both', levels = levels, cmap = "jet")            
    else :
        plt.imshow(arr2vorne, cmap="jet", vmax=maxVal, vmin=minVal)
    if not os.path.exists('PNG/Front'):
        os.makedirs('PNG/Front')
    pngfile = "PNG/Front/" + basename + ".png"
    plt.savefig( pngfile )
    os.system("convert -trim "+pngfile+" "+pngfile)
    plt.clf()

    fig = plt.figure(figsize=(8,8))
    plt.axis('off')
    cb = plt.colorbar(im0,orientation='vertical',pad=0.05,shrink=0.99)

    labels = [minVal + i * (maxVal - minVal)/4 for i in range(5)]
    cb.set_ticks( labels )
    cb.ax.tick_params(labelsize=20)
    if maxVal > 999:
        ticklabels = ["${:0.2}$".format(d) for d in labels]
    else:
        ticklabels = ["${:0.2f}$".format(d) for d in labels]
    cb.set_ticklabels( ticklabels )
    if not os.path.exists('PNG/Legend'):
        os.makedirs('PNG/Legend')
    pngfile = "PNG/Legend/" + basename + ".png"
    plt.savefig( pngfile )
    os.system("convert -trim "+pngfile+" "+pngfile)
    plt.clf()

if __name__ == '__main__' :
    matplotlib.rcParams.update({'font.size': 22})
    matplotlib.rcParams.update({'font.family': 'serif'})
    if int(sys.argv[3]) == 3 :
        plot3d (hdf5file = sys.argv[1], hdf5group = sys.argv[2], minVal = float(sys.argv[4]), maxVal = float(sys.argv[5]), steps = int(sys.argv[6]), xStart = float(sys.argv[7]), style = sys.argv[8], color = sys.argv[9] )
    else :
        if int(sys.argv[3]) == 2 :
            plot2d (hdf5file = sys.argv[1], hdf5group = sys.argv[2], minVal = float(sys.argv[4]), maxVal = float(sys.argv[5]), steps = int(sys.argv[6]), xStart = float(sys.argv[7]), style = sys.argv[8], color = sys.argv[9])
        else :
            raise (IndexError, 'dimension must be 2 or 3')

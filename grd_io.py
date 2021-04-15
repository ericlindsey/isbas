# -*- coding: utf-8 -*-
"""
Routines to read and write a GMT-format netCDF4 file

Created on Tue May  3 20:49:09 2016
@author: elindsey
"""

from scipy.io import netcdf
import numpy as np
import sys,subprocess

def read_grd(fname):
    ncfile = open_ncfile(fname)
    # read the coordinates, and data in variable named 'z'.
    naming=get_naming(ncfile)
    if naming == 'xy':
        kx, ky = 'x', 'y'
    elif naming == 'lonlat':
        kx, ky = 'lon', 'lat'  
    x = ncfile.variables[kx][:].copy()
    y = ncfile.variables[ky][:].copy()
    z = ncfile.variables['z'][:].copy()
    ncfile.close()
    return x,y,z

def read_naming(fname):
    # for a file that is not already open
    ncfile = open_ncfile(fname)
    naming=get_naming(ncfile)
    ncfile.close()
    return naming

def get_naming(ncfile):
    # for a file that is already open
    if 'x' in ncfile.variables:
        naming = 'xy'
    elif 'lon' in ncfile.variables:
        naming = 'lonlat'
    else:
        keylist=[key for key in ncfile.variables.keys()]
        print('Error: netcdf coordinate system not (x,y) or (lon,lat). Found variables: %s'%keylist)
        sys.exit(1)
    return naming

def open_ncfile(fname):
    try:
        ncfile = netcdf.netcdf_file(fname,'r')
    except TypeError:
        print('Got TypeError while reading GRD file %s. Creating a temporary NetCDF version 3 copy.\nTo speed this process, convert the file beforehand (use nccopy -k classic <infile> <outfile>).'%fname)
        subprocess.call('nccopy -k classic %s %s.nc3'%(fname,fname), shell=True)
        ncfile = netcdf.netcdf_file('%s.nc3'%fname,'r')
    return ncfile


def grd_shape(fname):
    ncfile = netcdf.netcdf_file(fname,'r')
    # read the data in variable named 'z'.
    data = ncfile.variables['z'][:]
    shp = data.shape
    data = None
    act=ncfile.variables['z'].actual_range
    print(act)
    ncfile.close()
    return shp

def write_grd(x, y, z, filename, title='Created by write_grd using scipy.io.netcdf', naming='lonlat'):
    '''Write COARDS compliant netcdf (grd) file.
    Modified from the GmtPy library, c2009 Sebastian Heimann.
    http://emolch.github.com/gmtpy'''
    
    assert y.size, x.size == z.shape
    ny, nx = z.shape
    nc = netcdf.netcdf_file(filename, 'w')
    assert naming in ('xy', 'lonlat')

    if naming == 'xy':
        kx, ky = 'x', 'y'
    else:
        kx, ky = 'lon', 'lat'

    nc.node_offset = 0
    if title is not None:
        nc.title = title

    nc.Conventions = 'COARDS/CF-1.0'
    nc.createDimension(kx, nx)
    nc.createDimension(ky, ny)

    xvar = nc.createVariable(kx, np.float32, (kx,))
    yvar = nc.createVariable(ky, np.float32, (ky,))
    if naming == 'xy':
        xvar.long_name = kx
        yvar.long_name = ky
    else:
        xvar.long_name = 'longitude'
        xvar.units = 'degrees_east'
        yvar.long_name = 'latitude'
        yvar.units = 'degrees_north'

    zvar = nc.createVariable('z', np.float32, (ky, kx))
    
    xvar[:] = x.astype(np.float32)
    yvar[:] = y.astype(np.float32)
    zvar[:] = z.astype(np.float32)
    
    xvar.actual_range=np.float32([np.min(x),np.max(x)])
    yvar.actual_range=np.float32([np.min(y),np.max(y)])
    zvar.actual_range=np.float32([np.min(z),np.max(z)])

    nc.close()
    
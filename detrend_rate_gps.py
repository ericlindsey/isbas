#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Use a GPS file to fit and replace the trend in the InSAR file
GPS file format:
    lat lon num_pts rate
num_pts is the number of pixels around the GPS point to average in the InSAR .grd

Created on Fri Jan 19 10:26:59 2018

@author: elindsey
"""

import grd_io
import detrend_ts

import os, argparse
import numpy as np


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Remove trend estimated from GPS from an InSAR dataset')
    parser.add_argument('gpsfile',type=str,help='GPS los velocity file. Columns lon,lat,num_pts,rate. Required.')
    parser.add_argument('insarfile',type=str,help='InSAR .grd file. Required.')
    args = parser.parse_args()
    
    
    xvec,yvec,insardat = grd_io.read_grd(args.insarfile)
    print('read %s, dimensions (%d,%d)'%(args.insarfile,len(xvec),len(yvec)))

    X,Y = np.meshgrid(xvec,yvec)
    gpsdat = np.loadtxt(args.gpsfile)
    print('read %s, length %d'%(args.gpsfile,len(gpsdat)))    
    
    G,d=detrend_ts.construct_gpsconstraint(0,insardat,gpsdat,X,Y)
        
    #select the specified number of trend parameters
    G=G[:,0:3]

    #invert for a trend that matches the value around given GPS points
    m=np.dot(np.linalg.pinv(G),d)
    print(m)
    #reconstruct the model - a bit ugly here dealing with unknown number of trendparams
    onemat=np.ones(np.shape(X))
    model = detrend_ts.reconstruct_model_nparams(3,m,X,Y,onemat)  
        
    insar_detrend = insardat-model
    
    grdnaming=grd_io.read_naming(args.insarfile)
    outfile='%s_fitgps.grd'%os.path.splitext(args.insarfile)[0]
    
    print('writing %s'%outfile)
    
    grd_io.write_grd(xvec,yvec,insar_detrend,outfile,naming=grdnaming)
    
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

SBAS/ISBAS algorithm components

Created on Wed Jun  1 13:39:09 2016
@author: elindsey
"""

#local imports
import os,sys,subprocess
import numpy as np
import itertools

import grd_io
import detrend_ts

def read_baselines(fname):
    # read the baseline table and store in date-sorted order
    strdat=np.genfromtxt(fname,str,usecols=0)
    numdat=np.genfromtxt(fname,usecols=(1,2,4))
    sortorder=np.argsort(numdat[:,1])

    ids=strdat[sortorder]
    jdates=[str(int(np.floor(i))) for i in numdat[sortorder,0]]    
    dates=numdat[sortorder,1]
    bperp=numdat[sortorder,2]
    
    return ids,jdates,dates,bperp
    
def get_igram_ids(sat,fname,ids):
    # read intf.in and convert to a list of tuples of IDs
    if sat=='ALOS':
        istart=13
        iend=18
    elif sat=='ALOS2':
        istart=12
        iend=17
    elif sat=='S1':
        istart=3
        iend=11
        #use only the date part of the id
        for i in range(len(ids)):
            ids[i]=ids[i][15:23]
    else:
        print('Error: satellite %s not yet implemented. Edit the code to provide location of ID number in granule name.'%sat)
        sys.exit(1)
    igrams=np.genfromtxt(fname,dtype=str)
    igram_ids=[]
    for igram in igrams:
        igramsplit=igram.split(':')
        strid0=igramsplit[0][istart:iend]
        strid1=igramsplit[1][istart:iend]
        id0=np.where(ids==strid0)[0][0]
        id1=np.where(ids==strid1)[0][0]
        igram_ids.append( (id0,id1) )
    return igram_ids
    
def get_g(igram_ids,dates):
    # very "pythonic" but hard to read: iterates over each interferogram (for k in range(len(igram_ids)))
    # and produces a row of values which are either zero or the timespan between the i+1th and ith dates, depending whether the igram covers that time
    idarr=np.arange(1,len(dates))
    G=[[int(boolval)*(dates[i+1]-dates[i]) for i,boolval in enumerate(np.logical_and(igram_ids[k][0]<idarr, igram_ids[k][1]>=idarr))] for k in range(len(igram_ids))]
    return np.array(G)
    
def load_igrams(trackdir, unwrap, igram_ids, jdates):
    data=None
    i=0
    for igram in igram_ids:
        fname='%s/intf/%s_%s/%s'%(trackdir,jdates[igram[0]],jdates[igram[1]],unwrap)
        xvec,yvec,zdata=grd_io.read_grd(fname)
        if data is None:
            #allocate array now that we know the x,y size
            data=np.zeros((len(igram_ids),len(yvec),len(xvec)))
        data[i,:,:]=zdata
        i=i+1
    return xvec,yvec,data

def ref_igrams(X,Y,data,reflon,reflat,refradius):
    for i in range(len(data)):
        xn,yn,meannear,mediannear=detrend_ts.get_near_data(X,Y,data[i,:,:],reflon,reflat,refradius)
        data[i,:,:] -= mediannear
    return data

def check_data(xvec,yvec,data,nt,igram_ids,unw_check_threshold,grdnaming):
    # function to check all loops and identify pixels where the loop value is close to zero. Mark these pixels as valid for each igram in the triplet
    # create a NaN matrix the same size as data (a 3D matrix containing all of the interferograms)
    validdata=np.zeros(np.shape(data))
    
    # output directory
    os.makedirs('unw_check', exist_ok=True)
    #iterate over loop list
    for loop in itertools.combinations(range(nt),3):
        #get interferogram IDs for each loop
        try:
            i0=igram_ids.index((loop[0],loop[1]))
            i1=igram_ids.index((loop[1],loop[2]))
            i2=igram_ids.index((loop[0],loop[2]))
        except ValueError:
            #one of the interferograms does not exist, skip this loop
            continue
        #compute the loop residual
        cycle=data[i0,:,:] + data[i1,:,:] - data[i2,:,:]
        #get zero if the pixel has a large residual, otherwise 1
        
        #cyclevalid checks if the values are less than the threshold. Any NaN values in the loop will just result in zero vote
        mask=~np.isnan(cycle)
        mask[mask]=np.abs(cycle[mask])<unw_check_threshold
        cyclevalid=1*mask
        #note, testing 'less than' on NaN throws a RuntimeWarning, but still gives the correct answer:
        #cyclevalid=1*(np.abs(cycle)<threshold)

        #validdata tracks the number of total 'votes' for a good pixel for each interferogram.
        validdata[i0,:,:] += cyclevalid
        validdata[i1,:,:] += cyclevalid
        validdata[i2,:,:] += cyclevalid
        
        loopname='%d_%d_%d'%(loop[0],loop[1],loop[2])
        grd_io.write_grd(xvec, yvec, cycle, 'unw_check/loop_%s.grd'%loopname, naming=grdnaming)
        grd_io.write_grd(xvec, yvec, cyclevalid, 'unw_check/cyclevalid_%s.grd'%loopname, naming=grdnaming)

    # mask pixels with too few votes.
    # note, it's possible for a bad pixel to still get 1 vote, if there is an opposing unwrapping error in another ifg.
    minvalid=2
    data[validdata<minvalid]=np.nan
        
    # TODO: get the interferogram directory and write out valid / invalid / (number of votes)
    
    for i in range(len(validdata)):
        grd_io.write_grd(xvec, yvec, validdata[i,:,:], 'unw_check/valid_%d.grd'%i, naming=grdnaming)

    return data


def timeseries_invert(G,data,dates,ts_type,mingood):
    if (ts_type == 'SBAS'):
        # do the standard SBAS inversion
        print('Running standard SBAS...')
        ts_dirname='sbas'
        model = sbas_invert(G, data) 
    elif (ts_type == 'ISBAS'):
        # do the ISBAS inversion
        nt=len(dates)
        if (mingood > 0):
            # use pixels with min. number of good interferograms
            print('Running improved ISBAS for pixels with minimum %d good interferograms...'%mingood)
            ts_dirname='nsbas_%d'%mingood
        else:
            # use pixels with a non-zero determinant of G
            print('Running improved ISBAS for pixels which are fully constrained (full-rank G)...')
            ts_dirname='nsbas_full'
        model = isbas_invert(G, data, nt, mingood)
    else:
        print('Error: timeseries type %s not recognized'%ts_type)
        sys.exit(1)
    print('Done.')
    # reconstruct the timeseries
    print('Reconstructing the timeseries')
    ts = reconstruct_ts(model, dates)
    return ts_dirname,ts
    
def sbas_invert(G,data):
    Ginv=np.linalg.pinv(G)
    model=np.tensordot(Ginv,data,axes=1)
    return model
    
def isbas_invert(G,data,nt,mingood=0):
    # version with for loops
    ni,nx,ny=np.shape(data)
    model=np.zeros((nt-1,nx,ny))*np.nan
    progressbar=10
    numgood=0
    Gstore={}
    for i in range(nx):
        for j in range(ny):
            # get a subset of G which operates on just the ifgs for which the pixel is correlated
            Igood = np.where(~np.isnan(data[:,i,j]))[0]            
            if (len(Igood) > 0 ):
                Ggood=G[Igood,:]
                #option 1: mingood > 0: check that the pixel has more good interferograms than 'mingood'
                #option 2: mingood <= 0: check if Ggood is full-rank
                if ( (mingood > 0 and len(Igood) > mingood) or (mingood <= 0 and np.linalg.matrix_rank(Ggood) == nt-1) ):
                    # pseudoinverse of G - doing this inside the loop is slow, so cache results
                    Ikey=hash(Igood.tobytes())
                    if Ikey in Gstore:
                        Ginv=Gstore[Ikey]
                    else:
                        Ginv=np.linalg.pinv(Ggood)
                        Gstore[Ikey]=Ginv
                    # compute the model ts components
                    model[:,i,j]=np.dot(Ginv,data[Igood,i,j])
                    numgood+=1
            if(i*j/(nx*ny)*100 > progressbar):
                
                print('%d%% finished, %d good pixels so far'%(progressbar,numgood))
                progressbar+=10
    print('100%% finished, %d good pixels total'%(numgood))
    return model
    
def reconstruct_ts(model,dates):
    # the timeseries calculation returns the displacements between each epoch, divided by the timespan.
    # here we sum and scale them to return the final timeseries.
    nt=len(dates)
    ntmodel,nx,ny=np.shape(model)   
    ts=np.zeros((nt,nx,ny))
    for k in range(1,nt):
        ts[k]=ts[k-1]+model[k-1]*(dates[k]-dates[k-1])
    return ts

def check_ts(ts,data,igram_ids,xvec,yvec,unw_check_threshold,grdnaming):
    #predict each interferogram from the timeseries, and compare it to the actual
    os.makedirs('ts_check', exist_ok=True)
    for i,pair in enumerate(igram_ids):
        synth_igram=ts[pair[1],:,:]-ts[pair[0],:,:]
        resid=synth_igram-data[i,:,:]
        grd_io.write_grd(xvec, yvec, resid, 'ts_check/resid_%d.grd'%i, naming=grdnaming)
        #mask any data that exceed the threshold
        mask=~np.isnan(resid)
        mask[mask]=np.abs(resid[mask])>unw_check_threshold
        #print(i)
        #print('non-nan before masking: %d'%np.count_nonzero(~np.isnan(data[i])))
        data[i,mask]=np.nan
        #print('non-nan after masking:  %d'%np.count_nonzero(~np.isnan(data[i])))
        
    return data
    
def convert_rad2mm(ts,prmfile):
    with open(prmfile) as f:
        for line in f:
            if 'wavelen' in line:
                wavelen=float(line.split()[2])
    print('Converting from radians to mm and flipping the sign, using radar wavelength %f'%wavelen)
    tsmm=ts*(-1000)*wavelen/(4*np.pi)
    return tsmm
    
def calc_rate(ts,dates):
    nt,nx,ny=np.shape(ts)
    rate=np.zeros((nx,ny))*np.nan
    #progressbar=10
    for i in range(nx):
        for j in range(ny):
            # this is also a slow step, for a large number of pixels
            idx = np.isfinite(ts[:,i,j])
            if(len(dates[idx])==nt):
                pfit=np.polyfit(dates[idx],ts[idx,i,j],1)
                rate[i,j]=pfit[0]
            #if(i*j/(nx*ny)*100 > progressbar):
            #    print('%d%% finished'%progressbar)
            #    progressbar+=10
    #print('100% finished')
    #units of dates are in days. convert to years
    rate=rate*365.25
    return rate

def save_ts(outdir,xvec,yvec,dates,ts,grdnaming):
    os.makedirs(outdir, exist_ok=True)
    for i in range(len(dates)):
        grd_io.write_grd(xvec, yvec, ts[i], '%s/ts_mm_%04d.grd'%(outdir,dates[i]), naming=grdnaming)

def save_rate(outdir,xvec,yvec,rate,grdnaming):
    os.makedirs(outdir, exist_ok=True)
    grd_io.write_grd(xvec, yvec, rate, '%s/rate_mm_yr.grd'%outdir, naming=grdnaming)


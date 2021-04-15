# -*- coding: utf-8 -*-
"""
 SBAS/ISBAS algorithm driver program

Created on Wed Jun  1 13:39:09 2016
@author: elindsey
"""
#system-wide imports
import os, argparse, configparser, time
import numpy as np
import detrend_ts
import ts_func

if __name__ == '__main__':
    globalstart=time.time()
    parser = argparse.ArgumentParser(description='Run ISBAS/SBAS on a GMTSAR-formatted dataset')
    parser.add_argument('config',type=str,help='supply name of config file to setup processing options. Required.')
    args = parser.parse_args()
    
    ############################################################
    ##### read config file
    
    config=configparser.ConfigParser()
    config.optionxform = str #make the config file case-sensitive
    config.read(args.config)
    
    # timeseries options
    ts_type   = config.get('timeseries-config','ts_type')
    mingood   = config.getint('timeseries-config','nsbas_min_intfs')
    
    # ref. pixel options
    reflat    = config.getfloat('timeseries-config','reflat')
    reflon    = config.getfloat('timeseries-config','reflon')
    refnum    = config.getint('timeseries-config','refnum')
   
    # file input options
    SAT=config.get('timeseries-config','sat_name')
    unwrapfile= config.get('timeseries-config','unwrap_file')
    grdnaming = config.get('timeseries-config','grdnaming') # this is really unnecessary but we keep it for now
    prmfile   = config.get('timeseries-config','prm_file')
    intfs     = config.get('timeseries-config','intf_file')
    baselines = config.get('timeseries-config','baseline_file')
    
    #unwrapping check options
    unw_check_threshold = config.getfloat('timeseries-config','unw_check_threshold')
    
    # rate (velocity) options
    calcrate  = config.getboolean('timeseries-config','calcrate')
   
    # de-trending options
    detrend   = config.getboolean('timeseries-config','detrend')
    trendparams = config.getint('timeseries-config','trendparams')
    gpsfile   = config.get('timeseries-config','gps_file')
    constrainedtrend = config.getboolean('timeseries-config','constrainedtrend')
    
    # plotting options
    makeplots = config.getboolean('timeseries-config','makeplots')
    
    ############################################################
    ##### read and setup input data
    
    # read baseline table
    print('reading baselines from %s'%baselines)
    currdir=os.getcwd()
    ids,jdates,dates,bperp = ts_func.read_baselines(os.path.join(currdir,baselines))
    
    # get dates of all interferogram pairs and convert to integer indexes
    print('reading list of interferograms from %s' %intfs)
    igram_ids = ts_func.get_igram_ids(SAT, os.path.join(currdir,intfs), ids)
 
    # load data
    print('Reading unwrapped interferograms matching filename %s from directory %s'%(unwrapfile,currdir) )
    start=time.time()
    xvec,yvec,data = ts_func.load_igrams(currdir, unwrapfile, igram_ids, jdates)
    X,Y = np.meshgrid(xvec,yvec)
    print('Completed in %.1f seconds.'%(time.time()-start))
    
    ############################################################
    ##### do computations
    
    # construct the G matrix
    # print('Constructing the design matrix')
    G = ts_func.get_g(igram_ids, dates)
   
    # subtract the median value of the data around the reference point from each interferogram
    print('Referencing interferograms to median value of nearest %d pixels to point (%f,%f) '%(refnum,reflon,reflat))
    start=time.time()
    data = ts_func.ref_igrams(X,Y,data,reflon,reflat,refnum)
    print('Completed in %.1f seconds.'%(time.time()-start))
    
    # check for unwrapping errors:
    # if the sum of values around a cycle is not close to zero, that pixel probably has
    # an unwrapping error in one of the 3 interferograms. If a pixel never appears valid
    # in any cycle of 3 interferograms, we mask it with NaN.
    if unw_check_threshold>0:
        print('Checking for unwrapping errors: cycle voting method')
        start=time.time()
        data = ts_func.check_data(xvec,yvec,data,len(ids),igram_ids,unw_check_threshold,grdnaming)
        print('Completed in %.1f seconds.'%(time.time()-start))
    
    # invert igrams for timeseries using either sbas or nsbas    
    start=time.time()
    ts_dirname,ts = ts_func.timeseries_invert(G,data,dates,ts_type,mingood)
    print('Completed in %.1f seconds.'%(time.time()-start))

    # convert from radians to mm of displacement. Note we also flip the sign.
    tsmm = ts_func.convert_rad2mm(ts,os.path.join(currdir,prmfile))
    
    if detrend:
        outdir = os.path.join(currdir,ts_dirname+'_detrend')
        # estimate trend and remove it
        if(constrainedtrend or gpsfile==''):
            if (gpsfile==''):
                print('Detrending using all data, without gps')
                gpsdat=np.array([[reflon,reflat,refnum,0]])
            else:
                print('Detrending using all data, with gps as constraint')
                #format: x0,y0,numpix,[v0, or vector of displacements at each satellite acquisition time]
                gpsdat=np.loadtxt(gpsfile)
            tsmm = detrend_ts.detrend_constraints(X,Y,tsmm,gpsdat,trendparams)
        else:
            print('Detrending using only data around gps points')
            #read GPS file
            #format: x0,y0,numpix,[v0, or vector of displacements at each satellite acquisition time]
            gpsdat=np.loadtxt(gpsfile)
            tsmm = detrend_ts.detrend_fitpoints(X,Y,tsmm,gpsdat,trendparams)
    else:
        outdir = os.path.join(currdir,ts_dirname)
    
    # compute velocities
    if calcrate:
        print('Computing velocity...')
        start=time.time()
        rate = ts_func.calc_rate(tsmm,dates)
        print('Completed in %.1f seconds.'%(time.time()-start))
    
    ############################################################
    ##### write out results
    
    # save output as grd files for each time step
    print('Saving timeseries steps as .grd files to directory %s'%outdir)
    ts_func.save_ts(outdir,xvec,yvec,dates,tsmm,grdnaming)   
    
    # save rates
    if calcrate:
        ts_func.save_rate(outdir,xvec,yvec,rate,grdnaming)
    
    ############################################################
    ##### make figures
    
    if makeplots:
        print('Making figure showing detrended velocity and sample timeseries')
        import matplotlib.pyplot as plt

        #plot velocity
        plt.figure(1)
        plt.clf()
        if detrend:
            plt.imshow(rate,interpolation='nearest')
        else:
            plt.imshow(rate,interpolation='nearest')
        plt.colorbar()
        plt.title('LOS rate (mm/yr)')
        plt.savefig('%s/rate_plot.png'%outdir)

        # in the future:
        # plot_ts.interactiveplot(xvec, yvec, dates, ts)
        # plot_ts.movie(xvec, yvec, dates, ts)
        #
        # working version for now:
        # plot timeseries as a consecutive series of images
        nt=len(dates)
        ny=int(np.floor(np.sqrt(nt)))
        nx=int(np.ceil(nt/ny))
        tsmin=np.min(tsmm)
        tsmax=np.max(tsmm)
        tsfig = plt.figure(figsize=(10, 10), constrained_layout=False)
        ax_grid = tsfig.add_gridspec(ny,nx, wspace=0, hspace=0)
        axs=ax_grid.subplots()
        for i,ax in enumerate(axs.flat):
            if i<nt:
                im=ax.imshow(tsmm[i,:,:],vmin=tsmin, vmax=tsmax)
                ax.xaxis.set_visible([])
                ax.yaxis.set_visible([])
                ax.set_title(jdates[i], x=0.1+nt/800, y=0.8-nt/700)
            else:
                ax.axis('off')
        cax = tsfig.add_axes([0.25,0, .5, 0.5])
        cax.axis('off')
        tsfig.colorbar(im, orientation="horizontal",label="LOS displacement (mm)")
        plt.savefig('%s/ts_plot.png'%outdir)

    print('Finished! Total time: %.1f seconds.\n'%(time.time()-globalstart))


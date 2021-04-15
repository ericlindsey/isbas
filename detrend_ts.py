# -*- coding: utf-8 -*-
"""
De-trending components for SBAS algorithm

Created on Wed Jul  6 18:06:01 2016
@author: elindsey
"""
import numpy as np

def detrend_fitpoints(X,Y,ts,gpsdat,trendparams=3):
    # for each date, fit a trend to velocities near the GPS points
    # (note, is there some better network-based way to do this?)
    #
    onemat=np.ones(np.shape(X))
    
    #estimate trends for each time series component
    for i in range(np.size(ts,0)):
        #extract data around the GPS points for use in fit equations
        G,d=construct_gpsconstraint(i,ts[i],gpsdat,X,Y)
        
        #select the specified number of trend parameters
        G=G[:,0:trendparams]
        
        #invert for a trend that matches the value around given GPS points
        m=np.dot(np.linalg.inv(G),d)
        
        #reconstruct the model - a bit ugly here dealing with unknown number of trendparams
        model = reconstruct_model_nparams(trendparams,m,X,Y,onemat)  
            
        ts[i]=ts[i]-model

    return ts

def detrend_constraints(X,Y,ts,gpsdat,trendparams=3):
    # for each date, remove a trend by fitting to the entire image,
    # but using a small number of GPS as constraints
    #
    onemat=np.ones(np.shape(X))
        
    #determine which method to use
    if(trendparams == range(np.size(gpsdat,0))):
        invmethod='inv'
    else:
        invmethod='pinv'
    
    #estimate trends for each time series component
    for i in range(np.size(ts,0)):
        #extract data around the GPS points for use in constraint equations
        C,z=construct_gpsconstraint(i,ts[i],gpsdat,X,Y)
        
        #remove nan points from the data and form G matrix
        d=ts[i].ravel()
        indx=~np.isnan(d)
        x=X.ravel()[indx]
        y=Y.ravel()[indx]
        onevec=onemat.ravel()[indx]
        d=d[indx]
        #Create the trend -fitting matrix
        #assumes 6 parameters here (this is the maximum for now), then remove extra parameters before fitting     
        G=np.column_stack([onevec,x,y,x**2,y**2,x*y])
        
        #select the specified number of trend parameters
        G=G[:,0:trendparams]
        C=C[:,0:trendparams]
        
        #invert for a trend subject to constraint that it matches the value around given GPS points
        m=constrained_lsq(G,d,C,z,method=invmethod)
        
        #reconstruct the model - a bit ugly here dealing with unknown number of trendparams
        model = reconstruct_model_nparams(trendparams,m,X,Y,onemat)      
            
        ts[i]=ts[i]-model

    return ts
    
def construct_gpsconstraint(imagenum,image,gpsdat,X,Y):
    # create the matrices used to fit a trend to values near the gps points
    #
    for j in range(np.size(gpsdat,0)):
        x0=gpsdat[j,0]
        y0=gpsdat[j,1]
        numpix=int(gpsdat[j,2])
        if (np.size(gpsdat,axis=1) == 4):
            v0=gpsdat[j,3] #use average velocity
        else:
            v0=gpsdat[j,3+imagenum] #use displacement specified at a particular epoch
        
             
        xn,yn,meannear,mediannear=get_near_data(X,Y,image,x0,y0,numpix)
        
        #constraint equation is that the sum of trend values at all the selected points is equal to the sum of their mean
        #assumes 6 parameters here (this is the maximum for now), remove extra parameters later before fitting
        Gpoint=np.array([numpix, np.sum(xn), np.sum(yn), np.sum(xn**2), np.sum(yn**2), np.sum(xn*yn)])
        dpoint=numpix*(meannear-v0)
        if (j==0):
            G=np.array(Gpoint,ndmin=2)
            d=dpoint
        else:
            G=np.row_stack(( G,Gpoint ))
            d=np.row_stack(( d,dpoint )) 
            
    return G,d
    
def get_near_data(X,Y,image,x0,y0,numpix):
    # get the mean and median of numpix points near x0,y0, and also return list of
    # the numpix nearest non-nan X,Y coordinates
    distarr=np.sqrt((X-x0)**2+(Y-y0)**2)
    distmask=np.ma.array(distarr,mask=np.isnan(image))
    nearindx=np.ma.argsort(distmask.ravel())[0:numpix]
    meannear=np.mean(image.ravel()[nearindx])
    mediannear=np.median(image.ravel()[nearindx])
    xn=X.ravel()[nearindx]
    yn=Y.ravel()[nearindx]
    return xn,yn,meannear,mediannear
    
def reconstruct_model_nparams(trendparams,m,X,Y,onemat):
    # given a model and a number of trend parameters, reconstruct the trend
    #
    model = {
    1: lambda m,X,Y,onemat: m[0]*onemat,
    2: lambda m,X,Y,onemat: m[0]*onemat+m[1]*X,
    3: lambda m,X,Y,onemat: m[0]*onemat+m[1]*X+m[2]*Y,
    4: lambda m,X,Y,onemat: m[0]*onemat+m[1]*X+m[2]*Y+m[3]*X**2,
    5: lambda m,X,Y,onemat: m[0]*onemat+m[1]*X+m[2]*Y+m[3]*X**2+m[4]*Y**2,
    6: lambda m,X,Y,onemat: m[0]*onemat+m[1]*X+m[2]*Y+m[3]*X**2+m[4]*Y**2+m[5]*X*Y
    }[trendparams](m,X,Y,onemat) 
    
    return model
    
def constrained_lsq(G,d,C,z=0,method='inv'):
    # minimize:   Gm = d 
    # subject to: Cm = z
    # This function works best when the KKT matrix is invertible in the classic sense, and you use the standard 'inv' method
    # Underdetermined/overconstrained problems may fail inelegantly - accurate results with pinv() are not guaranteed!
    #
    # we introduce dummy Lagrange multipliers l, and solve the KKT equations:
    # (from http://stanford.edu/class/ee103/lectures/constrained-least-squares/constrained-least-squares_slides.pdf)
    # [2*G.T*G  C.T] [m] = [2*G.T*d]
    # [C         0 ] [l] = [z]
    #
    #double check the shape of d
    d=np.array(d,ndmin=2).T
    if (np.size(d,0)==1):
        d=d.T
    #double check the shape of C - avoids problem if only one constraint is passed
    C=np.array(C,ndmin=2)
    #form the big matrix -- sorry, these double-parentheses are hard to read
    K=np.column_stack(( np.row_stack(( 2*np.dot(G.T,G),C )),np.row_stack(( C.T,np.zeros(( np.size(C,0),np.size(C,0) )) )) ))
    #form the data column
    if (np.size(z)!=np.size(C,0)):
        z=np.zeros((np.size(C,0),1)) #if z was given incorrectly (or it was left at the default zero)
    D=np.row_stack(( 2*np.dot(G.T,d),z ))
    #compute the model vector
    if(method=='inv'):
        m=np.dot(np.linalg.inv(K),D)
    elif(method=='pinv'):
        m=np.dot(np.linalg.pinv(K),D)
    #separate the model from lagrange multipliers
    mm=m[0:np.size(G,1)]
    return mm    
    #ml=m[np.size(G,1),-1]
    #return mm,ml

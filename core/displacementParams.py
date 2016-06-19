# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 16:30:27 2016

@author: benedikt
modified by TCO
"""


import numpy as np
from tracking3d.core  import utilsForProcessing
#import re

class Speed:#class to measure displacements parameters: mean speed, average velocity, and tumbling rate 
    def __init__(self,stackDistanceClass,numPart):
        self.MO = stackDistanceClass.MO #magnification of the microscope objective
        self.pxSize = stackDistanceClass.pxSize #reconstructed pixel sitz from stackDistance class
        self.index_medium = stackDistanceClass.index_medium #refractive index of the medium surrounding particles
        self.index_chamber = stackDistanceClass.index_chamber #refractive index of the material of the chamber (coverslip for example)
        self.directoryPath = stackDistanceClass.directoryPath
        self.numPart = numPart
        self.meanSpeed = None
        self.AverageVelocity = None
        self.TumblingRate = None

        
    def ComputeParams(self):
        z_init = 0
        tt = 0# tumbling time
        # arrays                                            
        tab1 = []
        tab2 = []
        slist = []
        #tstep = []
        ## Initiating data
        timestampsPath = self.directoryPath+'\\timestamps.txt'
        px_py_d_vol_meanampPath = self.directoryPath+'\\SingleParticuleMethod\\Part'+str(self.numPart).zfill(3)+'\\px_py_d_vol_meanamp.txt'
        tab1 = np.genfromtxt(timestampsPath,dtype=str)
        tab2 = np.genfromtxt(px_py_d_vol_meanampPath, dtype=str, skip_header=1)
        nb_frames = np.shape(tab2)[0] #number of frames for the tracked particule

        # hologram numbers
        holonr = tab2[:,0]
        holonr = holonr.astype(np.int)
        first_holo = holonr[0] #first hologram of the sequence (can be different than 0)

        # time
        ttab = tab1[first_holo:first_holo+nb_frames,3]
        t = ttab.astype(np.float)*1e-3                 # from millisecond to s
        
                
        # x-value
        xtab = tab2[:,1]
        x_px = xtab.astype(np.float)             # in pixels
        x = x_px*self.pxSize*1e-6                # in m
        
        # y-value
        ytab = tab2[:,2] 
        y_px = ytab.astype(np.float)             # in pixels
        y = y_px*self.pxSize*1e-6                # in m  
        
        
        # z-value
        dist_rec = tab2[:,3]            
        dist_rec = dist_rec.astype(np.float)            # in cm
        zconvert = self.index_medium*1e4/(self.MO**2*self.index_chamber)                # in um
        z = zconvert*(dist_rec-dist_rec[0])*1e-6+z_init        # conversion factor dist_rec to z in m
        
        voltab = tab2[:,4] #volume computed
        vol = voltab.astype(np.float)

        
        meanamptab = tab2[:,5] #mean amplitude computed
        meanamp = meanamptab.astype(np.float)
        
        ## Mean speed computation
        for i in range (0, len(holonr)-1):
            #print i
            dx = x[i+1]-x[i]
            dy = y[i+1]-y[i]
            dz = z[i+1]-z[i]
            ds = np.sqrt((dx)**2+(dy)**2+(dz)**2)
            
            slist.append(ds)
            
        #    dt = t[i+1]-t[i]                           # specific time step from i to (i+1)
        #    tstep.append(dt)
        
        s_tot = sum(slist)                              # total distance travelled in m
        t_tot = t[-1]-t[0]                                  # total time passed in s
        mv = s_tot/t_tot                                # mean speed in m/s
        
        ## Average velocity
        s_a = np.sqrt((x[-1]-x[0])**2+(y[-1]-y[0])**2+(z[-1]-z[0])**2)  # in m
        
        v_a = s_a/t_tot                                 # average velocity in m/s
        
        ## tumbling rate
        for j in range(1, len(holonr)-1):
               
            dx1 = x[j]-x[j-1]
            dy1 = y[j]-y[j-1]
            dz1 = z[j]-z[j-1]
            #dz1=0
            d1 = np.sqrt(dx1**2+dy1**2+dz1**2)        
           
            dx2 = x[j+1]-x[j]
            dy2 = y[j+1]-y[j]
            dz2 = z[j+1]-z[j]
            #dz2=0
            d2 = np.sqrt(dx2**2+dy2**2+dz2**2)
            theta = np.arccos((dx1*dx2+dy1*dy2+dz1*dz2)/(d1*d2))
            theta *= 180./np.pi #angle in degrees
            #theta = float(np.arccos((x[j]-x[j-1])*(x[j+1]-x[j])+(y[j]-y[j-1])*(y[j+1]-y[j])+(z[j]-z[j-1])*(z[j+1]-z[j])/(d1*d2)))
            #print theta
            
#            if theta <= float((180.-12.)/360.*2.*np.pi): 
#                tt = t[j+1]-t[j]+tt
            if np.abs(theta) > 12.:
                tt += t[j+1]-t[j]
                
        tr = tt/(t[-1]-t[1])
        self.SaveTrackingParams(holonr,ttab,x,y,z,vol,meanamp)
        self.meanSpeed = mv
        self.AverageVelocity = v_a
        self.TumblingRate = tr
        self.SaveSpeedParams()

    def SaveTrackingParams(self,holonr,time,x,y,z,volume,meanamp):
        #Save coordinates for the single particule tracking
        f = open(self.directoryPath+'\\SingleParticuleMethod\\Part'+str(self.numPart).zfill(3)+'\\Tracking_values.txt','w')
        #f=open(self.stackDistance.directoryPath+'\\Part'+str(numPart).zfill(3)+'_coord.txt','w')
        #f.writelines('holoNumber'+' '+'x[px]'+' '+'y[px]'+' '+'z[um]'+' '+'volume[um3]'+' '+'meanamp'+'\n')
        f.writelines('holoNumber\ttime[ms]\tx[m]\ty[m]\tz[m]\tvolume[um3]\tmeanamp\n')
        
        for k in range(np.shape(holonr)[0]):
            holo = str(holonr[k]).zfill(5)
            t = str(time[k])
            px = str(x[k])
            py = str(y[k])
            d = str(z[k])
            vol = str(volume[k])
            amp = str(meanamp[k])
            f.writelines(holo+'\t'+t+'\t'+px+'\t'+py+'\t'+d+'\t'+vol+'\t'+amp+'\n')
        f.close()
    def GetSpeedParams(self):
        return [self.meanSpeed,self.AverageVelocity,self.TumblingRate]
    def SaveSpeedParams(self):
        ColumnTitles = ['Mean Speed [m/s]','Average velocity [m/s]','Tumbling Rate']
        Values = [self.meanSpeed,self.AverageVelocity,self.TumblingRate]
        fname = self.directoryPath+'\\SingleParticuleMethod\\Part'+str(self.numPart).zfill(3)+'\\SpeedParams.txt'
        utilsForProcessing.SaveParamsFile(ColumnTitles,Values,fname)     
    
            





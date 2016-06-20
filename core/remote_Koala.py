# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 11:49:13 2016

@author: zerodeku
"""

from tracking3d.core import utilsForProcessing as utils
from tracking3d.core import remoteControls

import cv2
import PIL.Image
import os

import time
import numpy as np

class RemoteKoala:
    def __init__(self, config, ip_Koala, configType, directoryPath, 
                 first_holo_seq):
        self.config = config #config Number associated to MO in Koala
        self.ip_Koala = ip_Koala #ip of remoted Koala computer
        self.wavelength = None #source wavelength
        self.configType = configType #configtype 1=single wavelength,2=dual 
        #wavelength
        self.directoryPath = directoryPath #path of sequence directory
        self.remoteCommands = remoteControls #class remote Controls
        self.first_holo_seq = first_holo_seq #first hologram of the sequence 
        #we want to analyze
        self.password = None #password for the remote access
        self.connected = False #remote status
        self.pxSize = 1 #reconstructed pixel size
        self.width = 1024 #reconstructed image width
        self.height = 1024 #reconstructed image  height
        self.lastconfig = 0 #last configuration opened
        self.lastfirst_holo_seq = -1 #last hologram opened
        self.amp = None #amplitude image
        self.ph = None#phase image
        self.seqFormat = 1 #sequence format 0=avi (avi has to be saved through 
        #ImageJ in tif sequence holo_0000i.tif),1=raw,2=tif sequence
        self.holoWidth = 1024 #hologram width
        self.holoHeight = 1024 #hologram height
        self.bpp = None #bit per pixel
        self.framesNumber = None #frames Number of the sequence
        
    @classmethod
    def from_dataset(cls, db, remoteParams): #initialize remoteKoala class 
        k = cls(db.ProductionNo, db.IPKoala, db.configType, 
                remoteParams.directoryPath, remoteParams.first_holo_seq) 
        return k
    def DefineSeqFormat(self): #Define the format of the sequence (avi,raw,tif)
        holodir = self.directoryPath+'\\Holograms'
        if os.path.isfile(holodir+'\\holo.avi'):
            #WARNING, not possible to acquire holo directly from holo.avi, the 
            #sequence has to be saved in tif sequence with name holo_00000.tif,
            #holo_00001.tif,... this can be done with ImageJ: drag and drop 
            #holo.avi, File/Save As/Image Sequence
            #Format=TIFF, Name=holo_,Strat at 0, Digit=5
            self.seqFormat = 0
            cap = cv2.VideoCapture(holodir+'\\holo.avi')
            self.holoWidth = (int)(cap.get(3))
            self.holoHeight = (int)(cap.get(4))
            #print cap.get(7)
            self.framesNumber = (int)(cap.get(7))
            #print self.framesNumber
            self.bpp = 8
            cap.release()
            
        elif os.path.isfile(holodir+'\\holo.raw'):#Raw format
            self.seqFormat = 1
            kbin_header_dtype = np.dtype([("width", "i4"), ("height", "i4"), ("bpp", "i4"), 
                                          ("frames_number", "i4")])
            f = open(holodir+'\\holo.raw', 'rb')
            kbin_header = np.fromfile(f, dtype=kbin_header_dtype, count=1) 
            self.holoWidth = kbin_header['width']
            self.holoHeight = kbin_header['height']   
            self.bpp = kbin_header['bpp']
            self.framesNumber = kbin_header['frames_number']
        else:
            self.seqFormat = 2
            self.framesNumber = self.frames_Number(holodir, '.tif')#holo numbers
            self.bpp = 8
    def frames_Number(self, path, extension):
        list_dir = []
        count = 0
        try:
            list_dir = os.listdir(path)
            for files in list_dir:
                if files.endswith(extension):
                    count += 1
        except:
            count = -1
          
        return count
    def path_xth_holo(self, x): #compute path of the xth hologram
        fname = self.directoryPath+'\\Holograms'
        holopath = self.directoryPath+'\\Holograms\\holo.tif'
        if self.seqFormat == 0:#avi format
            #holo.avi has to be saved in tiff format (see DefineSeqFormat)
            holopath = self.directoryPath+'\\Holograms\\holo_'+\
            str(x).zfill(5)+'.tif'
            
        elif self.seqFormat == 1:
            f = open(fname+'\\holo.raw', 'rb')           
            beginat = x*(self.holoHeight*self.holoWidth)
            size = self.holoHeight*self.holoWidth
            f.seek(beginat+16, 1)
            tmp = np.fromfile(f, dtype='uint8', count=size)
            f.close
            Z = np.reshape(tmp, (self.holoHeight, self.holoWidth))          
            holopath = self.directoryPath+'\\Holograms\\holo.tif'
            img = PIL.Image.fromarray(Z)
            img.save(holopath)
        
        elif self.seqFormat == 2:
            holopath = self.directoryPath+'\\Holograms\\'+str(x).zfill(5)+\
            '_holo.tif'
            self.frames_number = self.frames_Number(self.directoryPath+\
            '\\Holograms', '.tif')
            self.bpp = 8
  
        return holopath
        
    def ConnectKoala(self):
        #connect to Koala
        [ret, user] = self.remoteCommands.Connect(self.ip_Koala, 'user')
        if self.password == None:
             #ask for password 
            passwd = utils.getPassword("Password for " 
                                                    +user, 'admin')        
            self.password = passwd
        self.remoteCommands.Login(self.password)#login to remote with password
        self.connected = True #connected to Koala

    def OpenWindows(self):
        try:
            self.remoteCommands.OpenWindows() #try to open windows until ok
        except:
            time.sleep(1)
            self.OpenWindows()
    def OpenConfig(self):
        self.remoteCommands.OpenConfig(self.config)#open configuration
        self.OpenWindows() #open windows (holo,intensity, phase)
        self.lastconfig = self.config #last config open is the actual config
        
    def OpenHologram(self):
        path = self.path_xth_holo(self.first_holo_seq)#path of the xth hologram
        self.remoteCommands.OpenHologram(path, self.configType)
        time.sleep(.01)
        #get the reconstructed pixel size from Koala
        self.pxSize = self.remoteCommands.GetPxSizeUm()
        #get phase image width from Koala
        self.width = self.remoteCommands.GetPhaseWidth()
        #get phase image height from Koala
        self.height = self.remoteCommands.GetPhaseHeight()
        #last first Holo of the sequence if the actual
        self.lastfirst_holo_seq = self.first_holo_seq
    def OpenHoloFromPath(self, path):
        #remote open hologram from path and configType
        self.remoteCommands.OpenHologram(path, self.configType)
    
    def CloseConnection(self):
        self.remoteCommands.Logout() #remote logout
        self.connected = False #connection is False
    def AcquireWavefront(self):
        #Get intensity image
        self.amp = self.remoteCommands.getIntensity32fImage('lambda1', 
                                                            self.width,
                                                            self.height)
        #get phase image
        self.ph = self.remoteCommands.getPhase32fImage('lambda1', 
                                                       self.width, 
                                                       self.height)
    def AcquireWavefrontBcg(self):
        #Get intensity image
        amp = self.remoteCommands.getIntensity32fImage('lambda1', 
                                                       self.holoWidth, 
                                                       self.holoHeight)
        #get phase image
        ph = self.remoteCommands.getPhase32fImage('lambda1', self.holoWidth, 
                                                  self.holoHeight)
        return amp, ph
        
    def GetLambdaNm(self):
        #get the wavelength from Koala
        self.wavelength = self.remoteCommands.GetLambdaNm(1)
    def new_hv_seg(self, hv_seg):
        #add segments in Koala for tilt correction
        self.remoteCommands.AddCorrSegmentAndComputePhaseCorrection(hv_seg, 2)
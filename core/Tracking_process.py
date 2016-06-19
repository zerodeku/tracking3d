# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 14:36:15 2014

@author: tcolomb
"""


from scipy import ndimage
import scipy.stats


from tracking3d.core import binkoala
from tracking3d.core import utilsForProcessing as utils
from tracking3d.core import remoteControls
from tracking3d.core import ThresholdMethod
from tracking3d.core import UnwrapK
from tracking3d.core import displacementParams
from tracking3d.core import preprocessing
from tracking3d.core import cplxprocessing
import cv2
import PIL.Image
import os
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from threading import Thread
from itertools import islice
import sys, time
from PyQt4 import QtGui, QtCore
import numpy as np
import shutil
            
class remoteKoala:
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
        self.seqFormat = 2 #sequence format 0=avi (avi has to be saved through 
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
#            cap=cv2.VideoCapture(fname+'\\holo.avi')
#            #ret, frame =cap.read()
#            frameCounter=0
#            play=True
#            while(play):
#                ret,frame=cap.read()
#                if frame is not None:
#                    if frameCounter==x:
#                        img=PIL.Image.fromarray(frame)
#                        img.save(holopath)
#                        play=False
#                    frameCounter+=1
#                else:
#                    play=False
#            cap.release()
            
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
        time.sleep(.1)
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
class aberrationCorrection:
    #class for aberr. compensation using filtering (!not profile)
    def __init__(self, filter_size_aberration):
        self.filter_size_aberration = filter_size_aberration#uniform filter size
        self.ph = None #phase image
        self.complex = None #complex image =amp*exp(i*ph)
        self.ph_filter = None #filtered phase
        self.amp_filter = None#filtere amplitude
    @classmethod
    def from_dataset(cls, recParams): 
        #init class from dataset recParams.filter.size_aberration
        k = cls(recParams.filter_size_aberration)
        return k   
    def PhaseUniformFilter(self): #Uniform filter only on the phase
        tmp_out = np.exp(complex(0., 1.)*self.ph) #conversion phase to complex
        tmp_real = np.real(tmp_out) #real part of tmp_out
        tmp_img = np.imag(tmp_out) # imaginariy part of tmp_out
        #uniform filter on real and imaginaray part
        tmp_real = ndimage.uniform_filter(tmp_real, 
                                          size=self.filter_size_aberration)
        tmp_img = ndimage.uniform_filter(tmp_img, 
                                         size=self.filter_size_aberration)
        #recombine real and imaginary part in a complex image
        tmp_out = tmp_out*np.exp(-complex(0., 1.)*\
            np.angle(tmp_real+tmp_img*complex(0., 1.)))
        #extract filter phase from the complex image
        self.ph_filter = np.angle(tmp_out)
    def ComplexUniformFilter(self): #uniform filter on the complex wavefront
        #same procedure thant PhaseUniformFilter but using  amplitude and phase
        tmp_out = self.complex
        tmp_real = np.real(tmp_out)
        tmp_img = np.imag(tmp_out)
        tmp_real = ndimage.uniform_filter(tmp_real, 
                                          size=self.filter_size_aberration)
        tmp_img = ndimage.uniform_filter(tmp_img, 
                                         size=self.filter_size_aberration)
        tmp = tmp_real+tmp_img*complex(0., 1.)
        self.ph_filter = np.angle(tmp_out*np.exp(-complex(0., 1.)*\
                        np.angle(tmp)))
        self.amp_filter = np.abs(tmp_out)-np.abs(tmp)
        
class stackDistance:
    #class to reconstruct amplitude and phase stack for different 
    #reconstruction distances and some methods
    def __init__(self, CCDPixelSizeUM, NA, ROI, transmission, index_medium, 
                 index_chamber, index_sample, roi_sample, central_rec_dist, 
                 range_rec_dist, ratio_DOF, first_holo_seq, unwrap, 
                 use_filter_size_aberration, filter_size_aberration,  
                 use_filter_amplitude, bool_deriv, threshold, 
                 bool_ThresholdMethod, particlecontrast, 
                 it_eros, it_dil, configType, directoryPath, saveAllStack):
        #self.pixel_size_um=pixel_size_um
        self.CCDPixelSizeUM = CCDPixelSizeUM
        self.MO = None #magnification of the microscope objective
        self.NA = NA #numerical aperture of the microscope objective
        self.ROI = ROI# ROI in Koala obtained from Database
        self.transmission = transmission #DHM model transmission or reflection
        self.index_medium = index_medium #refractive index of the medium 
                                        #surrounding particles
        self.index_chamber = index_chamber #refractive index of the material 
                                    #of the chamber (coverslip for example)
        self.index_sample = index_sample #refractive index of the particle
        self.roi_sample = roi_sample #Region of interest in pixel around the 
                                        #xy position used to evaludate focus
        self.central_rec_dist = central_rec_dist #central rec distance
        self.range_rec_dist = range_rec_dist #range of the rec distance 
                                    #d=central_rec_dist+-range_rec_dist
        self.ratio_DOF = ratio_DOF #Depth of field ratio, larger is the ratio, 
                        #larger are the stack image number for a given range 
        self.first_holo_seq = first_holo_seq #first holo of the seq. to analyze
        self.unwrap = unwrap
        self.use_filter_size_aberration = use_filter_size_aberration #booleean 
                                                    #True=use uniform filtering
        self.filter_size_aberration = filter_size_aberration #filter size of 
                #the uniform filtering used for aberraction compensation class
        self.use_filter_amplitude = use_filter_amplitude #boolean 
                        #True=ComplexUniformFilter, False=PhaseUniformFilter
        self.bool_deriv = bool_deriv #use derivation to try eliminating non 
                                #moving particle between holo i+1 and holo i
        self.threshold = threshold #threshold value for particle detection 
                                    #(0<threshold<255) (images are converted 
                                        #in 8bit between min and max values)
        self.bool_ThresholdMethod = bool_ThresholdMethod #Threshold method: 
                                                #Manual, Huang, otsu,....
        self.particlecontrast = particlecontrast #particle contrast, 
                                #positive=particle brighter than background, 
                                    #negative=particle darker that background
        self.it_eros = it_eros #number of iteration for erosion on the mask 
                                            #image defined by threshold method
        self.it_dil = it_dil #number of dilatation for erosion on the mask 
                                            #image defined by threshold method
        self.configType = configType #config Type 1=single, 2=dual-wavelength
        self.directoryPath = directoryPath #sequence directoryPath
        self.saveAllStack = saveAllStack #True=save stack for each hologram, 
                                    #False, save only two stacks i and i+1
        self.phaseProjection = None #phase Projection (maximal value on the 
                                                        #stack for each pixel)
        self.phaseProjection0 = None #phase Projection of holo i-1
        self.phaseProjection1 = None #phase Profjection of holo i
        self.holoNumber = None #hologram number
        self.step_rec_dist = None #rec. dist step (defined by NA and ratio_DOF)
        self.remoteKoala = None #class remote Koala
        self.number_stack_image = None #number of images in stack
        self.threshold_filter = None #mask defined by threshold method, 
                                        #threshold value, it_eros and it_dil
        self.xy_position = None #xy position of particles
        self.x = None #x position
        self.y = None #y position
        self.z = None #z position
        self.t = None #time position (hologram number)
        self.volume = None #volume measurement
        self.mean_amp = None #mean amplitude
        self.pxSize = None #reconstructed pixel size
        self.projectedWithPoints = None #projection image with draw xy points
        self.focusMethod = 0 #Autofocus method
        self.UnwrapK = UnwrapK.UnwrapK()
        self.singleParticuleXYPositionDetectionMethod = 0 #0=maxValue, 
                                                            #1=with threshold
        self.BcgWavefront = None#Background wavefront at given d
        self.BcgCorrection = False#True=Background remove
        self.prop = None#propagation method
        self.apo = None#apodization method
        self.bcg_rec_dist = 0.#rec. dist to acquire background
        self.RecordedBcgWavefront = None#Recorded background wavefront
    @classmethod
    def from_dataset(cls, db, dhmParams, sampleParams, holderParams, recParams, 
                     remoteParams): 
        #init class from dataset and database ObjectiveDB.db3
        k = cls(db.CCDPixelSizeUM, db.NA, db.ConfigROI, db.TransmissionDHM, 
                sampleParams.index_medium, sampleParams.index_chamber, 
                sampleParams.index_sample, holderParams.roi_sample, 
                holderParams.central_rec_dist, holderParams.range_rec_dist, 
                holderParams.ratio_DOF, remoteParams.first_holo_seq, 
                recParams.unwrap, recParams.use_filter_size_aberration, 
                recParams.filter_size_aberration, 
                recParams.use_filter_amplitude, recParams.bool_deriv, 
                recParams.threshold, recParams.bool_ThresholdMethod, 
                recParams.particlecontrast, recParams.it_eros, 
                recParams.it_dil, db.configType, remoteParams.directoryPath, 
                remoteParams.saveAllStack)
        return k
        
    def ResetPhaseProjection(self):
        #Reset all phase Projection images
        self.phaseProjection = None
        self.phaseProjection1 = None
        self.phaseProjection0 = None
    
    def SamplePlaneStepUM(self):
        #convert the step of rec dist in real displacement in the sample plane
        dof = self.remoteKoala.wavelength/(self.NA**2)*1e-7
        self.step_rec_dist = self.MO**2*dof/self.index_medium/self.ratio_DOF*\
            self.index_chamber
        zconvert = 1e4/self.MO**2*self.index_medium/self.index_chamber
        stepUm = self.step_rec_dist*zconvert
        stackHeightUM = self.range_rec_dist*2*zconvert
        return stepUm, stackHeightUM
    def ComputeMO(self):
        #Compute the effective magnification and not the one from MO value
        self.MO = self.CCDPixelSizeUM/self.pxSize
    def InitBcgPropagation(self):
        width = self.remoteKoala.holoWidth
        height = self.remoteKoala.holoHeight
        self.apo = preprocessing.apodization(width, height)
        if os.path.exists(self.directoryPath+'\\Background'):
            fname = utils.FindFileInDirectory(self.directoryPath+\
                    '\\Background', '.bin')
            if np.shape(fname)[0] != 0:
                fname = fname[0]
                fname = str(fname)
                rec_dist_text = fname.translate(None, 'BcgWavefront_bin')
                rec_dist_text = rec_dist_text[0:-1]
                rec_dist = (float)(rec_dist_text)      
                self.bcg_rec_dist = rec_dist
                self.RecordedBcgWavefront = \
                                binkoala.read_mat_cplx_bin(self.directoryPath\
                                +'\\Background\\'+fname)
                size = np.shape(self.RecordedBcgWavefront)
                wavelength = self.remoteKoala.wavelength*1e-9
                CCDPixelSize = self.CCDPixelSizeUM*1e-6
                self.prop = cplxprocessing.fresnel_propagation(size[0], size[1], 0, 
                                                               CCDPixelSize, wavelength)
        
    def AcquireAndSaveStack(self, adjustTiltEachTime=False, 
                            onlyAmplitudeComputation=False):
        #Acquire and save Stack for a given hologram
        #adjustTilteEachTime: can be modified if user want to compensate for 
        #tilt for each hologram
        #onlyAmplitudeComputation if only the amplitude want to be analyzed, 
        #by default is False
        
        #Test if objects are not None
        if self.remoteKoala is None:
            utils.ErrorMessage('AcquireAndSaveStack: remoteKoala is None')
        if self.holoNumber is None:
            utils.ErrorMessage('AcquireAndSaveStack: holoNumber is wrong')
            return
        if self.remoteKoala.wavelength is None:
            utils.ErrorMessage('AcquireAndSaveStack: wavelength is not defined')
            return
        dof = self.remoteKoala.wavelength/(self.NA**2)*1e-7 #define the depth 
        #of focus of the MO converstion of the depth of focus to reconstruction
        #distance step, have to take in account the magnification and the 
        #ref index of the medium and chamber and the user defined ratio_DOF
        self.step_rec_dist = self.MO**2*dof/self.index_medium/self.ratio_DOF*\
                            self.index_chamber
        #number of images in the stack is two times the range of reconstruction
                            #distance divided by step
        self.number_stack_image = (int)(np.ceil(2*self.range_rec_dist/\
                                self.step_rec_dist))+1
        #define the name of the directory where the stack images are saved
        if self.saveAllStack:
            directoryNumber = self.holoNumber#the holo number for all stack
        else:
            directoryNumber = self.holoNumber%2 #modulo 2 if not all stack
        #create directory if they do not exist, and delete previous data 
        utils.CreateDirectory(self.directoryPath+\
                            '\\Intensity', str(directoryNumber).zfill(5))
        utils.DeleteAllFilesInDirectory(self.directoryPath+\
                            '\\Intensity\\'+str(directoryNumber).zfill(5))
        utils.CreateDirectory(self.directoryPath+'\\Phase', 
                                           str(directoryNumber).zfill(5))
        utils.DeleteAllFilesInDirectory(self.directoryPath\
                                +'\\Phase\\'+str(directoryNumber).zfill(5))
        
        #utils.DeleteAllFilesInDirectory(self.directoryPath+'\\Projection')
        path = self.remoteKoala.path_xth_holo(self.holoNumber)#define holo path
                                                            #from holoNumber
        self.remoteKoala.OpenHoloFromPath(path)#open the hologram from path 
                                                            #through remote
        if adjustTiltEachTime: #compute tilt correction if asked to be done for
                                #all holograms
            self.remoteKoala.remoteCommands.ComputePhaseCorrection(1, 1)
        #time.sleep(.1)
        #Get the reconstructed pixel size from Koala
        self.pxSize = self.remoteKoala.pxSize
        recdist = []
        #iteration from 0 to the number of stack images to save
        for k in range(0, self.number_stack_image):
            #define reconstruction distance
            d = -self.range_rec_dist+k*self.step_rec_dist+self.central_rec_dist
            #set rec dist in Koala
            self.remoteKoala.remoteCommands.SetRecDistCM(d)
            recdist.append(d)
            pathIntensity = self.directoryPath+'\\Intensity\\'\
                            + str(directoryNumber).zfill(5)+'\\intensity_'+\
                            str(k).zfill(3)+'.bin'
            pathPhase = self.directoryPath+'\\Phase\\'+\
            str(directoryNumber).zfill(5)+'\\phase_'+str(k).zfill(3)+'.bin'
            #save amplitude
            self.remoteKoala.remoteCommands.SaveAmpFloatToFile('lambda1', 
                                                               pathIntensity)
            #save phase
            self.remoteKoala.remoteCommands.SavePhaseFloatToFile('lambda1', 
                                                                 pathPhase)
            if not self.use_filter_size_aberration and self.unwrap:
                self.UnwrapK.PerformUnwrapFromFile(pathPhase)
        #if aberraction compensation with uniform filter is  True
        if self.use_filter_size_aberration:
            #iterate for all previous saved data
            for k in range(0, self.number_stack_image):
                pathPhase = self.directoryPath+'\\Phase\\'+\
                str(directoryNumber).zfill(5)+'\\phase_'+str(k).zfill(3)+'.bin'            
                #read phase from .bin file
                [ph, header] = binkoala.read_mat_bin(pathPhase)
                #object correction from aberrationCorrection class
                correction = aberrationCorrection(self.filter_size_aberration)
                if self.use_filter_amplitude:#perform complex uniform filter
                    pathAmplitude = self.directoryPath+'\\Intensity\\'+\
                    str(directoryNumber).zfill(5)+'\\intensity_'+\
                    str(k).zfill(3)+'.bin'
                    #read amplitude image
                    [amp, headeramp] = binkoala.read_mat_bin(pathAmplitude)
                    #define complex child of object correction
                    correction.complex = amp*np.exp(complex(0., 1.)*ph)
                    correction.ComplexUniformFilter() #peform filtering
                    amp = correction.amp_filter #get the amp_filter
                
                else:#perform correction only on phase
                    correction.ph = ph
                    correction.PhaseUniformFilter()
                ph = correction.ph_filter #get the phase filter
                
                w = header[0][3] #width from bin file
                h = header[0][4] #height from bin file
                pxSizeM = header[0][5] #pixel size in meter from bin file
                self.pxSize = pxSizeM*1e6 #conversion in micron
                if not onlyAmplitudeComputation:
                    hconv = header[0][6] 
                    binkoala.write_mat_bin(pathPhase, ph, w, h, pxSizeM, 
                                           hconv, 1)#rewrite phase bin
                    if self.unwrap:
                        self.UnwrapK.PerformUnwrapFromFile(pathPhase)
                if self.use_filter_amplitude:
                    amp = correction.amp_filter
                    binkoala.write_mat_bin(pathAmplitude, amp, w, h, pxSizeM, 
                                           -1, 0) #rewrite amplitude bin
        if self.BcgCorrection:
#            if self.BcgWavefront is None:
#                print 'Bcg Wavefront is None'
#                return
            #iterate for all previous saved data  
            for k in range(0, self.number_stack_image):
                pathPhase = self.directoryPath+'\\Phase\\'+\
                str(directoryNumber).zfill(5)+'\\phase_'+str(k).zfill(3)+'.bin'
                #read phase from .bin file
                [ph, header] = binkoala.read_mat_bin(pathPhase)
                pathAmplitude = self.directoryPath+'\\Intensity\\'+\
                str(directoryNumber).zfill(5)+'\\intensity_'+\
                str(k).zfill(3)+'.bin' 
                [amp, headeramp] = binkoala.read_mat_bin(pathAmplitude)
                
                wavefront = amp*np.exp(complex(0., 1.)*ph)
                #self.BcgWavefront=self.BcgFromDistance(recdist[k])
                self.BcgFromDistance(recdist[k])
#               size=np.shape(self.BcgWavefront)
#                apo=preprocessing.apodization(size[0],size[1])  
#                wavefront=apo(wavefront)
                wavefront /= self.BcgWavefront
                
                amp = np.abs(wavefront)
                ph = np.angle(wavefront)
                
                w = header[0][3] #width from bin file
                h = header[0][4] #height from bin file
                pxSizeM = header[0][5] #pixel size in meter from bin file
                self.pxSize = pxSizeM*1e6 #conversion in micron
                hconv = header[0][6] 
                #rewrite phase bin
                binkoala.write_mat_bin(pathPhase, ph, w, h, pxSizeM, hconv, 1)
                if self.unwrap:
                    self.UnwrapK.PerformUnwrapFromFile(pathPhase)
                #rewrite amplitude bin
                binkoala.write_mat_bin(pathAmplitude, amp, w, h, pxSizeM, -1, 0)
             
         
    def BcgFromDistance(self, d):#compute the background for position d
        if self.ROI is None:
            'print self.ROI is None'
            return None
        else:
            self.prop.set_dist((d-self.bcg_rec_dist)*1e-2)#adjust prop distance
            tmp = self.prop(self.RecordedBcgWavefront)#perform propagation
            tmp = tmp[self.ROI[0]:self.ROI[0]+self.ROI[2], 
                      self.ROI[1]:self.ROI[1]+self.ROI[3]]
            self.BcgWavefront = tmp
            #return tmp
    def OpenPhFrombin(self, holo_number, dist_number):
        #Open phase image from bin according to holo number and distance
        if self.saveAllStack:
            path = self.directoryPath+'\\Phase\\'+str(holo_number).zfill(5)+\
            '\\phase_'+str(dist_number).zfill(3)+'.bin'
        else:
            path = self.directoryPath+'\\Phase\\'+str(holo_number%2).zfill(5)+\
            '\\phase_'+str(dist_number).zfill(3)+'.bin'
        [ph, header] = binkoala.read_mat_bin(path) 
        return[ph, header]
    def OpenAmpFrombin(self, holo_number, dist_number):
        #Open amplitude image from bin according to holo number and distance
        if self.saveAllStack:
            path = self.directoryPath+'\\Intensity\\'+\
            str(holo_number).zfill(5)+'\\intensity_'+\
            str(dist_number).zfill(3)+'.bin'
        else:
            path = self.directoryPath+'\\Intensity\\'+\
            str(holo_number%2).zfill(5)+'\\intensity_'+\
            str(dist_number).zfill(3)+'.bin'
        [amp, header] = binkoala.read_mat_bin(path) 
        return[amp, header]  
    def ZProjection(self):
        #Compute the phase projection as the maximal (or minimal if 
        #particuleconstrast is negative) value in z direction
        for k in range(0, self.number_stack_image):
            [ph0, header] = self.OpenPhFrombin(self.holoNumber, k)
            if k == 0:
                ph = ph0
            else:
                if self.particlecontrast == 0:
                    ph = np.maximum(ph, ph0)
                else:
                    ph = np.minimum(ph, ph0)
        #phase projection is the 
        if self.unwrap:
            self.phaseProjection = ph-np.nanmean(ph)
        else:
            self.phaseProjection = np.angle(np.exp(complex(0., 1.)*\
                                (ph-np.nanmean(ph))))            
        
    def Threshold(self):#Compute the threshold mask to be able to detect xy position
        #bool_ThresholdMethod= different methods available
            #-1=manual
            #0=Huang
            #1=Intermodes
            #2=MaxEntropy
            #3=ReniEntropy
            #4=Triangle
            #5=Otsu
            #6=Yen
            #7=Moments
            #8=IsoData
            #9=peakDetection
        #testthresholdin2steps=False
        sign = 1.
        if self.particlecontrast == 1:#if particle constrast is negative
            sign = -1.
        #if derivation is used (to limit impact of non moving particles)
        if self.bool_deriv:
            if self.unwrap:
                imtothreshold = self.phaseProjection1-self.phaseProjection0
                imtothreshold = imtothreshold-np.nanmean(imtothreshold)
            else:
                imtothreshold = np.angle(np.exp(complex(0., 1.)*\
                (self.phaseProjection1-self.phaseProjection0)))
                imtothreshold = np.angle(np.exp(complex(0., 1.)*\
                (imtothreshold-np.nanmean(imtothreshold))))  
            
            imtothreshold *= sign
        else:
            imtothreshold = sign*self.phaseProjection1
        #convert image to threshold in 8bit between min and max value
        minValue = np.min(imtothreshold)
        maxValue = np.max(imtothreshold)
        imtothreshold = (imtothreshold-minValue)/(maxValue-minValue)*255
        if self.bool_ThresholdMethod == -1:
            #Manual method, self.threshold is defined by the user in main panel
            imtothreshold = np.asarray(imtothreshold)
            imtothreshold = np.asarray(imtothreshold, dtype=np.uint8)
            mask = np.where(imtothreshold <= self.threshold, 0, 1)
        if self.bool_ThresholdMethod == 5:
            #Otsu method (implemented in cv2 library)
            #imtothreshold=np.asarray(imtothreshold)
            imtothreshold = np.asarray(imtothreshold, dtype=np.uint8)
            ret, mask = cv2.threshold(imtothreshold, 0, 255, cv2.THRESH_OTSU)
            self.threshold = (int)(ret) #computed threshold
        if self.bool_ThresholdMethod != -1 and\
            self.bool_ThresholdMethod != 5 and\
            self.bool_ThresholdMethod != 9:
            #initialize object of class ThresholdMethod
            ThresholdClass = \
            ThresholdMethod.ThresholdMethod(self.bool_ThresholdMethod, 
                                            imtothreshold)
            #compute the threshold according to method and image
            ThresholdClass.ComputeThreshold()
            #get the threshold from object
            self.threshold = ThresholdClass.threshold
            #construct mask
            mask = np.where(imtothreshold <= self.threshold, 0, 1)
        if self.bool_ThresholdMethod == 9:#peak detect method from ndimage
            imtothreshold = np.asarray(imtothreshold)
            imtothreshold = np.asarray(imtothreshold, dtype=np.uint8)
            filt2 = np.copy(imtothreshold)
            filt2[filt2 < self.threshold] = 0 #apply threshold
            #detect peaks from ndimage library
            labeled_image, number_of_objects = ndimage.label(filt2)
            peak_slices = ndimage.find_objects(labeled_image)
            mask = np.zeros(np.shape(imtothreshold)).astype(np.float32)
            for peak_slice in peak_slices:
                dx, dy = peak_slice
                mask[dx.start:dx.stop, dy.start:dy.stop] = 1
        if self.threshold == -1:
            utils.ErrorMessage("Error on Threshold evaluation")
        if np.sum(mask) == 0: #Warning if threshold is too high and mask=0
            print 'threshold is too high, no particle are detected'
#        if self.bool_ThresholdMethod==0:
#            thresh=(self.threshold+180.)/360.*255
#            
#            img=255.*np.asarray((imtothreshold+np.pi)/2/np.pi)
#            img=np.asarray(img,dtype=np.uint8)
#            img=cv2.GaussianBlur(img,(5,5),0)
#            ret,filt=cv2.threshold(img,thresh,255,cv2.THRESH_BINARY)
#            if filt.max()==0:
#                print 'threshold is too high, no particle are detected'
#            else:
#                filt=filt/filt.max()
        if not self.it_eros == 0: #erosion of the mask defined by self.it_eros
            mask = ndimage.binary_erosion(mask, 
                                          iterations=self.it_eros).astype(np.float32)     
        if not self.it_dil == 0:#dilatation of the new mask
            #dilatation to avoid several points very close
            mask = ndimage.binary_dilation(mask, 
                                           iterations=self.it_dil).astype(np.float32)   
        if np.sum(mask) == 0:#Warning if threshold is too high and mask=0
            print 'threshold is too high, no particle are detected'
# Test code: 
#        if testthresholdin2steps and self.bool_deriv:
#            filt=ndimage.binary_dilation(filt,iterations=self.roi_sample).astype(np.float32)
#            img=self.phaseProjection1*filt
#            labeled_image, number_of_objects = ndimage.label(img)
#            peak_slices = ndimage.find_objects(labeled_image)
#            filt=np.zeros(np.shape(img)).astype(np.float32)
#            for peak_slice in peak_slices:
#                area=np.zeros(np.shape(img)).astype(np.float32)
#                dx,dy=peak_slice
#                area[dx.start:dx.stop,dy.start:dy.stop]=1.
#                area=area*self.phaseProjection1
#                index_max=area.argmax()
#                idx=np.unravel_index(index_max,np.shape(img))
#                filt[idx[0],idx[1]]=1.

        self.threshold_filter = mask #define self.threshold_filter from computed mask
    
    def ThresholdSingleParticule(self, imtothreshold):
        minValue = np.min(imtothreshold)
        maxValue = np.max(imtothreshold)
        imtothreshold = (imtothreshold-minValue)/(maxValue-minValue)*255
        if self.bool_ThresholdMethod == -1:
            #Manual method, self.threshold is defined by the user in main panel
            imtothreshold = np.asarray(imtothreshold)
            imtothreshold = np.asarray(imtothreshold, dtype=np.uint8)
            mask = np.where(imtothreshold <= self.threshold, 0, 1)
        if self.bool_ThresholdMethod == 5:
            #Otsu method (implemented in cv2 library)
            imtothreshold = np.asarray(imtothreshold)
            imtothreshold = np.asarray(imtothreshold, dtype=np.uint8)
            ret, mask = cv2.threshold(imtothreshold, 0, 255, cv2.THRESH_OTSU)
            self.threshold = (int)(ret) #computed threshold
        if self.bool_ThresholdMethod != -1 and\
            self.bool_ThresholdMethod != 5 and\
            self.bool_ThresholdMethod != 9:
            #initialize object of class ThresholdMethod
            ThresholdClass = \
            ThresholdMethod.ThresholdMethod(self.bool_ThresholdMethod, 
                                            imtothreshold)
            #compute the threshold according to method and image
            ThresholdClass.ComputeThreshold()
            #get the threshold from object
            self.threshold = ThresholdClass.threshold
            mask = np.where(imtothreshold <= self.threshold, 0, 1)#construct mask
        if self.bool_ThresholdMethod == 9:#peak detection method from ndimage
            imtothreshold = np.asarray(imtothreshold)
            imtothreshold = np.asarray(imtothreshold, dtype=np.uint8)
            filt2 = np.copy(imtothreshold)
            filt2[filt2 < self.threshold] = 0 #apply threshold
            #detect peaks from ndimage library
            labeled_image, number_of_objects = ndimage.label(filt2)
            peak_slices = ndimage.find_objects(labeled_image)
            mask = np.zeros(np.shape(imtothreshold)).astype(np.float32)
            for peak_slice in peak_slices:
                dx, dy = peak_slice
                mask[dx.start:dx.stop, dy.start:dy.stop] = 1
        if np.sum(mask) == 0: #Warning if threshold is too high and mask=0
            print 'threshold is too high, no particle are detected'
#        if self.bool_ThresholdMethod==0:
#            thresh=(self.threshold+180.)/360.*255
#            
#            img=255.*np.asarray((imtothreshold+np.pi)/2/np.pi)
#            img=np.asarray(img,dtype=np.uint8)
#            img=cv2.GaussianBlur(img,(5,5),0)
#            ret,filt=cv2.threshold(img,thresh,255,cv2.THRESH_BINARY)
#            if filt.max()==0:
#                print 'threshold is too high, no particle are detected'
#            else:
#                filt=filt/filt.max()
        if not self.it_eros == 0: #erosion of the mask defined by self.it_eros
            mask = ndimage.binary_erosion(mask, 
                                          iterations=self.it_eros).astype(np.float32)     
        if not self.it_dil == 0:#dilatation of the new mask
            #dilatation to avoid several points very close
            mask = ndimage.binary_dilation(mask, 
                                           iterations=self.it_dil).astype(np.float32)  
        if np.sum(mask) == 0:#Warning if threshold is too high and mask=0
            print 'threshold is too high, no particle are detected'
            #if mask is only zeros, the program will stops, so we replace
            #by only ones values
            mask = np.ones(np.shape(mask))
# Test code: 
#        if testthresholdin2steps and self.bool_deriv:
#            filt=ndimage.binary_dilation(filt,iterations=self.roi_sample).astype(np.float32)
#            img=self.phaseProjection1*filt
#            labeled_image, number_of_objects = ndimage.label(img)
#            peak_slices = ndimage.find_objects(labeled_image)
#            filt=np.zeros(np.shape(img)).astype(np.float32)
#            for peak_slice in peak_slices:
#                area=np.zeros(np.shape(img)).astype(np.float32)
#                dx,dy=peak_slice
#                area[dx.start:dx.stop,dy.start:dy.stop]=1.
#                area=area*self.phaseProjection1
#                index_max=area.argmax()
#                idx=np.unravel_index(index_max,np.shape(img))
#                filt[idx[0],idx[1]]=1.

        return mask #define self.threshold_filter from computed mask
    def ExtractXYPositionSingleParticuleFromHolo(self, initialPosition, 
                                                 maxdisplacement, 
                                                 plotImage=False):
        method = self.singleParticuleXYPositionDetectionMethod
        #method=0:max value in region of research, 1 threshold method
        #Extract XY Position for Single particule tracking in the actual hologram
        #Initial position is the position xy for the previous hologram
        #maxdisplacement is the maximal displacement in pixel between holo i-1 and i
        #defined in the main panel by Max Displacement [pixels]
        
        #maximal distance to look for the particle=Roi size used for the 
        #focus+2*max displacement(2 times because it can be positive or 
        #negative direction)
        #h=max(self.roi_sample,2*maxdisplacement)
        h = 2*maxdisplacement
        if h == 0:
            h = self.roi_sample
        #define region of interest around the initial position according to h
        left = initialPosition[0]-h/2 
        top = initialPosition[1]-h/2
        left = max(0, left) #left and top cannot be smaller than 0
        top = max(0, top)
        sign = 1.
        if self.particlecontrast == 1:#change sign if particle contrast<0
            sign = -1.
        if self.bool_deriv:# define image to threshold by derivation if True
            if self.unwrap:
                imtothreshold = self.phaseProjection1-self.phaseProjection0
                imtothreshold -= np.nanmean(imtothreshold)
            else:
                imtothreshold = np.angle(np.exp(complex(0., 1.)\
                    *(self.phaseProjection1-self.phaseProjection0)))
                imtothreshold = np.angle(np.exp(complex(0., 1.)\
                    *(imtothreshold-np.nanmean(imtothreshold))))
            imtothreshold *= sign
        else:
            imtothreshold = sign*self.phaseProjection1
        #crop image to ROI:
        imtothreshold = imtothreshold[left:left+h, top:top+h]
        
        if method == 0:#take the maximum of image
            xy_position = np.unravel_index(imtothreshold.argmax(), 
                                           imtothreshold.shape)
            xy_img = np.copy(xy_position)
            xy_position = np.sum([xy_position, [left, top]], axis=0)
            imgToshow = imtothreshold
        if method == 1:#compute the center of mass of projection*mask
            #imtothreshold=imtothreshold[left:left+h,top:top+h]
            threshold_filter = self.ThresholdSingleParticule(imtothreshold)
#            plt.figure(48)
#            plt.imshow(threshold_filter)
            imgToshow = threshold_filter*imtothreshold
            centroid = ndimage.measurements.center_of_mass(threshold_filter*\
                        imtothreshold)
            xy_img = np.copy(centroid)
            xy_position = np.sum([centroid, [left, top]], axis=0)
        #print xy_img
        if plotImage:
            plt.figure(49)
            plt.gcf().clear()
            plt.imshow(imgToshow)
            plt.autoscale(False)
            plt.axhline(y=xy_img[0])
            plt.axvline(x=xy_img[1])
            plt.plot(xy_img[0], xy_img[1])
            plt.show()
#        plt.figure(50)
#        plt.imshow(self.threshold_filter)
        #compute image of projection with red point for measured xy position
        self.ProjectionSingleParticuleWithDetectedPoints(xy_position)
        return xy_position
        
    def ExtractXYPositionFromHolo(self):
        #Extract all the XY position in the entire image for Automatic tracking
        lab, num_features = ndimage.measurements.label(self.threshold_filter)
        #find center of mass = particle detection
        xy_position = ndimage.measurements.center_of_mass(self.threshold_filter, 
                                                          labels=lab, index=range(num_features+1))
        del xy_position[0]#first xy_position value is alway (nan,nan)
        self.xy_position = xy_position
        #compute images of projection with detected points
        self.ProjectionWithDetectedPoints()
    def zFromxy(self, plot=False):
        #compute z position from all the xy position
        if self.phaseProjection1 is None:
            utils.ErrorMessage("The phase projection is null, because you change the projection process\n Do Adjust Detection again.")
            return
        #mask=np.zeros((np.shape(self.phaseProjection1)[0],np.shape(self.phaseProjection1)[1]))
        if plot: #if user want to polt, values have to be intialized to empty array
            if np.shape(self.xy_position)[0] == 0:
                utils.ErrorMessage("You have to select a point before adjusting focus")
                return
            self.x = []
            self.y = []
            self.z = []
            self.t = []
            self.volume = []
            self.mean_amp = []
        for xyPoints in self.xy_position: #for all detected xypsosition
            self.x = np.append(self.x, xyPoints[0]) #array x of xypoints
            self.y = np.append(self.y, xyPoints[1]) #array y of xypoints
            #perform measurement of focus reconstruction distance
            [d_focus, vol, m] = self.focusAtxyposition(xyPoints, plot)
            self.z = np.append(self.z, d_focus) #array z with measured 
                                                #reconstruction distance
            self.t = np.append(self.t, self.holoNumber) #array time using 
                            #holoNumber (conversion as to be done with 
                            #timestamps.txt saved by koala)
            self.volume = np.append(self.volume, vol) #array of volume
            self.mean_amp = np.append(self.mean_amp, m) #array of mean amp
            #mask[xyPoints[0],xyPoints[1]]=1 #mask that defined the pixel 
            #x,y of position of particle
        #init rgbarray to convert to rgb image    
        rgbArray = np.zeros((np.shape(self.phaseProjection1)[0], 
                             np.shape(self.phaseProjection1)[1], 3), 'uint8')
        #gray level between 0 and 2pi for the projection iamge
        im = (self.phaseProjection1+np.pi)/2/np.pi
        rgbArray[:, :, 0] = im*255
        rgbArray[:, :, 1] = im*255
        rgbArray[:, :, 2] = im*255
        #save projection image in the right directory 
        #(create directory if not exists)
        img = PIL.Image.fromarray(rgbArray)
        utils.CreateDirectory(self.directoryPath, 'Projection')
        img.save(self.directoryPath+'\\Projection\\Proj_'+\
                        str(self.holoNumber).zfill(5)+'.jpg')

    
    def focusAtxyposition(self, xyPoints, plot, *args):  #arg=partNumber
        #compute the focus reconstruction distance for a single xy position
        focusCrit = np.zeros(self.number_stack_image)#init focus criteria
        #select the region of interest around xyposition to evaluate focus
        h = self.roi_sample
        left = xyPoints[0]-h/2
        top = xyPoints[1]-h/2
        left = max(0, left)
        top = max(0, top)
        for k in range(0, self.number_stack_image):
            #compute criteria for each reconstruction distance
            focusCrit[k] = self.FocusFromMethod(xyPoints, k)
        
        #Depending of focus criteria, min, max or other value defined the focus
        if self.focusMethod == 0:
            index_min = focusCrit.argmin()
        if self.focusMethod == 1:
            index_min = focusCrit.argmin()
        if self.focusMethod == 2:
            index_min = focusCrit.argmax()
        if self.focusMethod == 3:
            index_min = focusCrit.argmax()
        if self.focusMethod == 4:
            index_min = focusCrit.argmax()
        if self.focusMethod == 5:
            index_min = focusCrit.argmax()
        if self.focusMethod == 6:
            index_min = focusCrit.argmax()
        if self.focusMethod == 7:
            index_min = focusCrit.argmax() 
        if self.focusMethod == 8:
            i1 = focusCrit.argmin()
            i2 = focusCrit.argmax()
            index_min = (i1+i2)/2.
        #plot the focus criteria and the resulting amplitude and phase image
        if plot: 
            g = plt.figure(1)
            g.clear()
            plt.subplot(221)
            rec_dist = np.linspace(0, self.number_stack_image-1, self.number_stack_image)
            rec_dist *= self.step_rec_dist
            rec_dist += -self.range_rec_dist+self.central_rec_dist
            d = -self.range_rec_dist+index_min*self.step_rec_dist+self.central_rec_dist
            plt.plot(rec_dist, focusCrit, 'bo')
            plt.plot([d, d], [np.min(focusCrit), np.max(focusCrit)], 'k-', lw=2)
            plt.title('focusCrit')
            plt.xlabel('rec dist [cm]')

            [imgfocusamp, header] = self.OpenAmpFrombin(self.holoNumber, (int)(index_min))
            left_display = xyPoints[0]-150/2
            top_display = xyPoints[1]-150/2
            left_display = max(0, left_display)
            top_display = max(0, top_display)
            imgfocusamp = imgfocusamp[left_display:left_display+150, top_display:top_display+150]
            plt.subplot(223)
            plt.imshow(imgfocusamp, cmap='gray')
            plt.axis('off')
            [imgfocusph, header] = self.OpenPhFrombin(self.holoNumber, (int)(index_min))
            imgfocusph = imgfocusph[left_display:left_display+150, top_display:top_display+150]
            plt.subplot(224)
            plt.imshow(imgfocusph, cmap='gray')
            plt.axis('off')
            plt.show()
        
        #compute focus reconstruction distance from the computed index
        d_focus = -self.range_rec_dist+index_min*self.step_rec_dist
        
        [amp, header] = self.OpenAmpFrombin(self.holoNumber, (int)(index_min))
        amp = amp[left:left+self.roi_sample, top:top+self.roi_sample]
        m = np.mean(amp)
        #compute the volume of particle at position xy
        [ph, header] = self.OpenPhFrombin(self.holoNumber, (int)(index_min))
        ph = ph[left:left+self.roi_sample, top:top+self.roi_sample]
        mask = np.where(ph >= 0, 1., 0)#mask to take in account only 
        #the particle and not the surrounding medium noise
        delta_index = self.index_sample-self.index_medium  
        #volume measurement
        volume = np.sum(mask*ph)*delta_index*self.pxSize**2
        if len(args) != 0:
            #args=number of particule, is one if single particule tracking 
            #method, focused intensity and phase image are saved
            partNumber = (int)(args[0])
            
            #original phase file
            if self.saveAllStack:
                originalphasepath = self.directoryPath+'\\Phase\\'+\
                    str(self.holoNumber).zfill(5)+\
                    '\\phase_'+str((int)(index_min)).zfill(3)+'.bin'
                originalamppath = self.directoryPath+'\\Intensity\\'+\
                    str(self.holoNumber).zfill(5)+'\\intensity_'+\
                    str((int)(index_min)).zfill(3)+'.bin'
            else:
                originalphasepath = self.directoryPath+'\\Phase\\'+\
                    str(self.holoNumber%2).zfill(5)+'\\phase_'+\
                    str((int)(index_min)).zfill(3)+'.bin'
                originalamppath = self.directoryPath+'\\Intensity\\'+\
                    str(self.holoNumber%2).zfill(5)+'\\intensity_'+\
                    str((int)(index_min)).zfill(3)+'.bin'
            newphasepath = self.directoryPath+'\\SingleParticuleMethod\\Part'+\
                str(partNumber).zfill(3)+'\\Phase\\'+\
                str(self.holoNumber).zfill(5)+'_phase.bin'
            newampphase = self.directoryPath+'\\SingleParticuleMethod\\Part'+\
                str(partNumber).zfill(3)+'\\Intensity\\'+\
                str(self.holoNumber).zfill(5)+'_intensity.bin'
            
            shutil.copy2(originalphasepath, newphasepath)
            shutil.copy2(originalamppath, newampphase)
            
        return [d_focus, volume, m]
        
    def FocusFromMethod(self, xyPoints, positionStack):
        #Definition of the focus method
        #xyPoints is the xy position, position stack is the position in the 
        #stack (related to rec distance)
        #region of interest arround xy position
        h = self.roi_sample
        left = xyPoints[0]-h/2
        top = xyPoints[1]-h/2
        left = max(0, left) #left cannot be <0
        top = max(0, top) #top cannot be <0
        if self.focusMethod == 0 or self.focusMethod == 5: 
            #method min amp std or or max std deviation
            [img, header] = self.OpenAmpFrombin(self.holoNumber, positionStack)
            img = img[left:left+self.roi_sample, top:top+self.roi_sample]
            m = np.mean(img)
            fx = np.std(img)/m
        if self.focusMethod == 1: #absolute min amplitude value
            [img, header] = self.OpenAmpFrombin(self.holoNumber, positionStack)
            img = img[left:left+self.roi_sample, top:top+self.roi_sample]
            fx = np.min(img)
        if self.focusMethod == 2: #absolute max phase value
            [img, header] = self.OpenPhFrombin(self.holoNumber, positionStack)
            img = img[left:left+self.roi_sample, top:top+self.roi_sample]
            fx = np.max(img)
            #Test on amplitude
#            [img,header]=self.OpenAmpFrombin(self.holoNumber,positionStack)
#            img=img[left:left+self.roi_sample,top:top+self.roi_sample]
#            m=np.mean(img)
        if self.focusMethod == 3:
            [img, header] = self.OpenPhFrombin(self.holoNumber, positionStack)
            img = img[left:left+self.roi_sample, top:top+self.roi_sample]
            fx = np.sum(img)
            #Test on amplitude
#            [img,header]=self.OpenAmpFrombin(self.holoNumber,positionStack)
#            img=img[left:left+self.roi_sample,top:top+self.roi_sample]
#            m=np.mean(img)
        if self.focusMethod == 4:#Max phase(valuePos-mean)
            [img_entire, header] = self.OpenPhFrombin(self.holoNumber, positionStack)
            img = img_entire[left:left+self.roi_sample, top:top+self.roi_sample]
            m = np.mean(img)           
            value_xy = img_entire[xyPoints[0], xyPoints[1]]
            fx = np.abs(value_xy-m)
        if self.focusMethod == 6:#'Skewness'
            [img, header] = self.OpenAmpFrombin(self.holoNumber, positionStack)
            img = img[left:left+self.roi_sample, top:top+self.roi_sample]
            fx = np.abs(scipy.stats.skew(np.hstack(img), bias=False))
        if self.focusMethod == 7:#'Max std from Sobel'
            [img, header] = self.OpenAmpFrombin(self.holoNumber, positionStack)
            img = img[left:left+self.roi_sample, top:top+self.roi_sample]
            sx = ndimage.sobel(img, axis=0, mode='constant')
            sy = ndimage.sobel(img, axis=1, mode='constant')
            sob = np.hypot(sx, sy)
            m = np.mean(sob)
            fx = np.std(sob)/m
        if self.focusMethod == 8:#'Max std from Sobel on phase image'
            modelsobel = 'reflect'
            [img, header] = self.OpenPhFrombin(self.holoNumber, positionStack)
            img = img[left:left+self.roi_sample, top:top+self.roi_sample]
            sx = ndimage.sobel(img, axis=0, mode=modelsobel)
            sy = ndimage.sobel(img, axis=1, mode=modelsobel)
            sob = np.hypot(sx, sy)
            m = np.mean(sob)
            fx = np.std(sob)/m
#            [img,header]=self.OpenAmpFrombin(self.holoNumber,positionStack)
#            img=img[left:left+self.roi_sample,top:top+self.roi_sample]
#            m=np.mean(img)
            #fx=np.sum(img)
        return fx
        
    def FocusCriteria(self, amp):#certainly useless
        m = np.mean(amp)
        fx = np.std(amp)/m
        return [fx, m]
        
    def ProjectionSingleParticuleWithDetectedPoints(self, xypoint):
        #Computation of Projection image with xy position in red for single particule tracking
        #Mask = 1 only on xy position pixel
        mask = np.zeros(np.shape(self.phaseProjection))
        mask[xypoint[0], xypoint[1]] = 1
        #Create rectangle equal to roi_sample around the xy position
        x = min((np.shape(self.phaseProjection)[0])-self.roi_sample-1, 
                max(0, (int)(xypoint[0]-self.roi_sample/2)))
        y = min((np.shape(self.phaseProjection)[1])-self.roi_sample-1, 
                max(0, (int)(xypoint[1]-self.roi_sample/2)))
        for a in range(x, x+self.roi_sample):
            mask[a, y] = 1.
            mask[a, y+self.roi_sample] = 1.
        for b in range(y, y+self.roi_sample):
            mask[x, b] = 1.
            mask[x+self.roi_sample, b] = 1.
        #init rgbarray
        rgbArray = np.zeros((np.shape(self.phaseProjection1)[0], 
                             np.shape(self.phaseProjection1)[1], 3), 'uint8')
        im = (self.phaseProjection1-np.min(self.phaseProjection1))/\
                (np.max(self.phaseProjection1)-np.min(self.phaseProjection1))
        rgbArray[:, :, 0] = im*255
        rgbArray[:, :, 1] = im*255
        rgbArray[:, :, 2] = im*255
        #put red value where mask=1
        rgbArray[:, :, 1] = np.where(mask == 1, 0, rgbArray[:, :, 0])
        rgbArray[:, :, 2] = np.where(mask == 1, 0, rgbArray[:, :, 2])
        rgbArray[:, :, 0] = np.where(mask == 1, 255, rgbArray[:, :, 1])
        self.projectedWithPoints = rgbArray
    def ProjectionWithDetectedPoints(self):
        #computation of projection with detected points for automatic tracking
        for xyPoints in self.xy_position:
            x = min((np.shape(self.phaseProjection)[0])-self.roi_sample-1, 
                    max(0, (int)(xyPoints[0]-self.roi_sample/2)))
            y = min((np.shape(self.phaseProjection)[1])-self.roi_sample-1, 
                    max(0, (int)(xyPoints[1]-self.roi_sample/2)))
            for a in range(x, x+self.roi_sample):
                self.threshold_filter[a, y] = 1.
                self.threshold_filter[a, y+self.roi_sample] = 1.
            for b in range(y, y+self.roi_sample):
                self.threshold_filter[x, b] = 1.
                self.threshold_filter[x+self.roi_sample, b] = 1.

        rgbArray = np.zeros((np.shape(self.phaseProjection1)[0], 
                             np.shape(self.phaseProjection1)[1], 3), 'uint8')
        im = (self.phaseProjection1-np.min(self.phaseProjection1))/\
                (np.max(self.phaseProjection1)-np.min(self.phaseProjection1))
        rgbArray[:, :, 0] = im*255
        rgbArray[:, :, 1] = im*255
        rgbArray[:, :, 2] = im*255
        rgbArray[:, :, 1] = np.where(self.threshold_filter == 1, 0, rgbArray[:, :, 0])
        rgbArray[:, :, 2] = np.where(self.threshold_filter == 1, 0, rgbArray[:, :, 2])
        rgbArray[:, :, 0] = np.where(self.threshold_filter == 1, 255, rgbArray[:, :, 1])
        self.projectedWithPoints = rgbArray

        
            
class tracking(Thread):
    #class for tracking along time (along hologram sequence)
    def __init__(self, roi_sample, max_displacement):
        Thread.__init__(self)#define a thread to be able to interact with gui      
        self.roi_sample = roi_sample #region of interest around xy position
        self.max_displacement = max_displacement#max displacement of particle 
                                            #between holo i and i+1 in pixels
        self.stackDistance = None #object of class stackDistance
        self.pxSize = None #reconstructed pixel size
        self.last_holo_seq = None #last holo of sequence that has to be 
                                    #analyzed (defined by user in main panel)
        self.importXYPosition = False #import XY position from txt file 
                    #(implemented for a given customer, has to be coded for 
                                                    #other format of txt file)
        self.XYImportFilePath = None #path of the xy position text file
        self.pxpyd_particule = None #Position X Position Y Distance of single 
                                        #particule (single particule tracking)
    def run(self):
        #run thread ExtractXYZPositionFromHolo()
        self.ExtractXYZPositionFromHolo()
    @classmethod
    def from_dataset(cls, holderParams, sampleParams):
        #init object tracking from dataset
        k = cls(holderParams.roi_sample, sampleParams.max_displacement)
        return k
    def SaveCoordinates(self):
        #Save Coordinates xyz for all holograms (there is no tracking here, 
        #just a list of coordinates!)
        f = open(self.stackDistance.directoryPath+'\\coordinates.txt', 'w')
        f.writelines('holoNumber\tx[px]\ty[px]\tdelta_dist_rec[cm]\tvolume[um3]\tmeanamp\n')
        #number of coordinates for all holograms 
        n_coord = np.shape(self.stackDistance.x)[0]  
        for k in range(n_coord):
            time = str((int)(self.stackDistance.t[k])).zfill(5)#hologram number
            xpos = str(self.stackDistance.x[k])#x position
            ypos = str(self.stackDistance.y[k])#y position
            zpos = str(self.stackDistance.z[k])#z position
            volumepos = str(self.stackDistance.volume[k])#volume
            amp = str(self.stackDistance.mean_amp[k])#mean amplitude
            f.writelines(time+'\t'+xpos+'\t'+ypos+'\t'+zpos+'\t'+volumepos+'\t'+amp+'\n')
        f.close()
        if self.importXYPosition:
            #if xy position is imported from a txt file, add missing value 
            #(z in particular)
        #method for a given customer xy position file
            self.ExportModifiedDonner()
    def SaveSingleParticuleCoordinates(self, numPart):
        #Save coordinates for the single particule tracking
        f = open(self.stackDistance.directoryPath+\
                '\\SingleParticuleMethod\\Part'+str(numPart).zfill(3)+\
                '\\px_py_d_vol_meanamp.txt', 'w')
        f.writelines('holoNumber\tx[px]\ty[px]\tdist_rec[cm]\tvolume[um3]\tmeanamp\n')
        n_coord = np.shape(self.pxpyd_particule)[0]
        
        for k in range(n_coord):
            #♥holo=str((int)(self.stackDistance.t[k])).zfill(5)
            holo = str((int)(k+self.stackDistance.first_holo_seq)).zfill(5)
            px = str(self.pxpyd_particule[k][0])
            py = str(self.pxpyd_particule[k][1])
            d = str(self.pxpyd_particule[k][2])
            vol = str(self.pxpyd_particule[k][3])
            amp = str(self.pxpyd_particule[k][4])
            f.writelines(holo+'\t'+px+'\t'+py+'\t'+d+'\t'+vol+'\t'+amp+'\n')
        f.close()
    
    def files_number(self, path, extension):
        #count the number of file with the extenstion in a given pah
        list_dir = []
        count = 0
        try:
            list_dir = os.listdir(path)
            for file in list_dir:
                if file.endswith(extension):
                    count += 1
        except:
            count = -1
          
        return count
    
    def BcgCompute(self):
        value = self.last_holo_seq-self.stackDistance.first_holo_seq
        progressbar = ProgressBar(value) # create Progress Barr
        progressbar.show() #Show progress bar
        tic = time.clock() #absolute time to evaluate entire process time
        backgroundstate = self.stackDistance.BcgCorrection
        self.stackDistance.BcgCorrection = False
        holo_number = self.stackDistance.first_holo_seq
        last_holo = self.last_holo_seq-1
        #define holo path from holoNumber
        path = self.stackDistance.remoteKoala.path_xth_holo(holo_number)
        self.stackDistance.remoteKoala.OpenHoloFromPath(path)
        
        utils.CreateDirectory(self.stackDistance.directoryPath, '\\Background')
        pathIntensity = self.stackDistance.directoryPath+'\\Background\\intensity.bin'
        pathPhase = self.stackDistance.directoryPath+'\\Background\\phase.bin'
        #save amplitude
        self.stackDistance.remoteKoala.remoteCommands.SaveAmpFloatToFile('lambda1', pathIntensity)
        #save phase
        self.stackDistance.remoteKoala.remoteCommands.SavePhaseFloatToFile('lambda1', pathPhase)         
        [ph, header] = binkoala.read_mat_bin(pathPhase) #read phase from .bin file
        [amp, headeramp] = binkoala.read_mat_bin(pathIntensity)
                
        
        #amp,ph = self.stackDistance.remoteKoala.AcquireWavefrontBcg()
        rec_distance = self.stackDistance.remoteKoala.remoteCommands.GetRecDistCM()
        #amp=self.stackDistance.remoteKoala.amp
        #ph=self.stackDistance.remoteKoala.ph
        wavefront = amp*np.exp(complex(0., 1.)*ph)
        toc = time.clock() #absolut time
        timeremaining = (toc-tic)*self.last_holo_seq-1 #time remaining
        progressbar.SetTimeRemaining(timeremaining) #set time remaining in progress bar
        holo_number += 1
        while holo_number < last_holo+1:
            tic = time.clock()
            #define holo path from holoNumber
            path = self.stackDistance.remoteKoala.path_xth_holo(holo_number)
            self.stackDistance.remoteKoala.OpenHoloFromPath(path)
            self.stackDistance.remoteKoala.remoteCommands.SaveAmpFloatToFile('lambda1', 
                                                                             pathIntensity)
            self.stackDistance.remoteKoala.remoteCommands.SavePhaseFloatToFile('lambda1', 
                                                                               pathPhase)     
            [ph, header] = binkoala.read_mat_bin(pathPhase) #read phase from .bin file
            [amp, headeramp] = binkoala.read_mat_bin(pathIntensity)
            #amp,ph=self.stackDistance.remoteKoala.AcquireWavefrontBcg()
#            amp=self.stackDistance.remoteKoala.amp
#            ph=self.stackDistance.remoteKoala.ph
            wavefront += amp*np.exp(complex(0., 1.)*ph)
            value = holo_number-self.stackDistance.first_holo_seq
            progressbar.progressbar.setValue(value)
            QtGui.qApp.processEvents()
            #stop if "Stop button" is pressed on progressBar
            if not progressbar._active:
                break
            toc = time.clock()
            timeremaining = (toc-tic)*(self.last_holo_seq-holo_number)
            progressbar.SetTimeRemaining(timeremaining)
            holo_number += 1
        
        wavefront /= holo_number-1
        wavefront = self.stackDistance.apo(wavefront)
        wavefront = np.fft.fft2(wavefront)
        self.stackDistance.BcgWavefront = wavefront
        self.stackDistance.BcgCorrection = backgroundstate
        size = np.shape(wavefront)
        files = utils.FindFileInDirectory(self.stackDistance.directoryPath+\
                        '\\Background', '.bin')
        for fname in files:
            os.remove(self.stackDistance.directoryPath+'\\Background\\'+\
                    str(fname))
        binkoala.write_mat_cplx_bin(self.stackDistance.directoryPath+\
                    '\\Background\\BcgWavefront_'+str(rec_distance)+\
                    '.bin', wavefront, size[0], size[1])
        self.stackDistance.holoNumber = self.stackDistance.first_holo_seq
        progressbar._active = False
        
    def TrackSingleParticule(self, partNumber):
        #Method to track a single particule (Loop on particule will use this 
        #method successively)
        #create directory if they do not exists
        utils.CreateDirectory(self.stackDistance.directoryPath, '\\SingleParticuleMethod\\Part'+str(partNumber).zfill(3)+'\\Intensity')
        utils.CreateDirectory(self.stackDistance.directoryPath, '\\SingleParticuleMethod\\Part'+str(partNumber).zfill(3)+'\\Phase')
        #initialize position xy d of particule
        self.pxpyd_particule = []
        #initial position (partNumber-1 because partNumber=1,2,....). 
        #The initial position is defined by the user by selecting a point in 
        #the image
        initialPosition = self.stackDistance.xy_position[partNumber-1]
        #initial position should be x,y,z,volume and mean, so we add 0 value 
        #for volume and mean amplitude
        initialPosition += (0., )
        initialPosition += (0., )
        #central rec distance is define by the previous position rec distance
        self.central_rec_dist = initialPosition[2]
        #set central rec dist to object stackDistance
        self.stackDistance.central_rec_dist = initialPosition[2]
        #set holoNumber to first holo
        self.stackDistance.holoNumber = self.stackDistance.first_holo_seq
        # create Progress Barr
        value = self.last_holo_seq-self.stackDistance.first_holo_seq
        progressbar = ProgressBar(value)
        progressbar.show()#Show progress bar
        #absolute time at the beginning, in order to evaluate the time for the 
        #entire process
        tic = time.clock()
        #Acquire Stack for the first holo
        self.stackDistance.AcquireAndSaveStack()
        #Compute the Z Projection of the first holo
        self.stackDistance.ZProjection()
        #Define phaseProjection1
        self.stackDistance.phaseProjection1 = self.stackDistance.phaseProjection
        #if derivation, perform stack reconstruction for the holo "first holo+1"
        if self.stackDistance.bool_deriv:
            #phase projection=previous phase projection (=first holo)
            self.stackDistance.phaseProjection0 = self.stackDistance.phaseProjection
            self.stackDistance.holoNumber += 1#change the holoNumber to first +1
            self.stackDistance.AcquireAndSaveStack() #Acquire stack
            self.stackDistance.ZProjection()
            #phaseProjection1= computed phase Projection
            self.stackDistance.phaseProjection1 = self.stackDistance.phaseProjection
        #Extract xy position in the case of single particule tracking
        xypoint = self.stackDistance.ExtractXYPositionSingleParticuleFromHolo(initialPosition, 0)
        plt.ion()
        plt.figure(101)
        plt.imshow(self.stackDistance.projectedWithPoints)
        #replace intial XY Position by the computed xy position
        initialPosition = self.replace_position(initialPosition, xypoint)
        #compute z position for measured xy position
        [d_focus, volume, m] = self.stackDistance.focusAtxyposition(initialPosition, False, *[partNumber])
        d_focus += self.central_rec_dist #we must add central reconstruction 
                                            #distance to the computed d_focus 
        #replace reconstruction distance with the computed one
        initialPosition = self.replace_rec_dist(initialPosition, d_focus)
        # replace volumen and mean amplitude with the computed ones
        initialPosition = self.replace_vol_and_mean(initialPosition, volume, m)
        #add computed initial position to array position x,y,d,volume,mean amp
        self.pxpyd_particule.append(initialPosition)
        #for the next holo,phaseprojection0 becomes previous phaseprojection1
        self.stackDistance.phaseProjection0 = self.stackDistance.phaseProjection1
        toc = time.clock() #absolut time
        timeremaining = (toc-tic)*self.last_holo_seq-1 #time remaining
        progressbar.SetTimeRemaining(timeremaining) #set time remaining 
        k = 0
        #iterative computation along the number of hologram to last holo
        while self.stackDistance.holoNumber < self.last_holo_seq-1:
            k += 1
            tic = time.clock()
            time.sleep(0.05)
            progressbar.progressbar.setValue(self.stackDistance.holoNumber+1-\
                    self.stackDistance.first_holo_seq)
            QtGui.qApp.processEvents()
            if not progressbar._active:
                #stop if "Stop button" is pressed on progressBar
                plt.ioff()
                break
            #set the stackDistance hologram number to next one
            self.stackDistance.holoNumber += 1
            #set the central rec dist from previous computed rec dist
            self.central_rec_dist = initialPosition[2] 
            self.stackDistance.central_rec_dist = initialPosition[2]
            #Perform Stack acquis., Zprojection and detection of xyz,vol,mean
            self.stackDistance.AcquireAndSaveStack()
            self.stackDistance.ZProjection()
            self.stackDistance.phaseProjection1 = self.stackDistance.phaseProjection
            xypoint = self.stackDistance.ExtractXYPositionSingleParticuleFromHolo(initialPosition, 
                                                                                  self.max_displacement)
            plt.imshow(self.stackDistance.projectedWithPoints)
            initialPosition = self.replace_position(initialPosition, xypoint)
            [d_focus, volume, m] = self.stackDistance.focusAtxyposition(initialPosition, False, *[partNumber])
            d_focus += self.central_rec_dist
            initialPosition = self.replace_rec_dist(initialPosition, d_focus)
            initialPosition = self.replace_vol_and_mean(initialPosition, volume, m)
            
            #add the new position to array
            self.pxpyd_particule.append(initialPosition)
            self.stackDistance.phaseProjection0 = self.stackDistance.phaseProjection1
            toc = time.clock()
            timeremaining = (toc-tic)*(self.last_holo_seq-1-\
                            self.stackDistance.holoNumber)
            progressbar.SetTimeRemaining(timeremaining)
        
        plt.ioff()
        #process is finished or stop
        progressbar._active = False
        #Save single part coordinates for the particuler defined by partNumber
        self.SaveSingleParticuleCoordinates(partNumber)
        #Compute speed parameters

        ParamsSpeed = displacementParams.Speed(self.stackDistance, partNumber)        
        ParamsSpeed.ComputeParams()
        ParamsSpeed.SaveSpeedParams()
#        [meanSpeed,AverageVelocity,TumblingRate]=ParamsSpeed.GetSpeedParams()
#        print "Particule "+str(partNumber)+" speed parameters are:"
#        print "Mean speed: ", meanSpeed,"[m/s]"
#        print "Average velocity: ", AverageVelocity,"[m/s]"
#        print "Tumbling rate: ", TumblingRate  
        
        #Create 3D graph of the tracked particule
        self.ComputexyzpositionSingleParticuleMethod(partNumber)
        #init stackDistance holoNumber to the first holo seq
        self.stackDistance.holoNumber = self.stackDistance.first_holo_seq
        #reset phase projection
        self.stackDistance.phaseProjection = None
        self.stackDistance.phaseProjection0 = None
        self.stackDistance.phaseProjection1 = None

    def replace_rec_dist(self, tup, d):
        #replace reconstruction distance in the tuple
        lst = list(tup)
        lst[2] = d
        return tuple(lst)
    def replace_position(self, tup, pos):
        #replace position x,y in the tuple
        lst = list(tup)
        lst[0] = pos[0]
        lst[1] = pos[1]
        return tuple(lst)
    def replace_vol_and_mean(self, tup, vol, mean):
        #replace volume and mean amplitude in tuple
        lst = list(tup)
        lst[3] = vol
        lst[4] = mean
        return tuple(lst)
    def ExtractXYZPositionFromHolo(self, SaveOnlyStack):
        #Extract XYZ Position From Holo in automatic tracking
    
        #init x,y,z,t, volumen and mean amplitude
        self.x = []
        self.y = []
        self.z = []
        self.t = []
        self.volume = []
        self.mean_amp = []
        #init x,y,z,t, volumen and mean amplitude in object stackDistance
        self.stackDistance.x = np.array([])
        self.stackDistance.y = np.array([])
        self.stackDistance.z = np.array([])
        self.stackDistance.t = np.array([])
        self.stackDistance.volume = np.array([])
        self.stackDistance.mean_amp = np.array([])
        #init first hologram to analyze
        self.stackDistance.holoNumber = self.stackDistance.first_holo_seq
        progressbar = ProgressBar(self.last_holo_seq-\
                        self.stackDistance.first_holo_seq)
        progressbar.show()
        tic = time.clock()
        #Compute stack zprojection,.... for the first holo
        self.stackDistance.AcquireAndSaveStack()
        if not SaveOnlyStack:
            self.stackDistance.ZProjection()
            self.stackDistance.phaseProjection1 = self.stackDistance.phaseProjection
            if self.stackDistance.bool_deriv:
                self.stackDistance.phaseProjection0 = self.stackDistance.phaseProjection
                self.stackDistance.holoNumber += 1
                self.stackDistance.AcquireAndSaveStack()
                self.stackDistance.ZProjection()
                self.stackDistance.phaseProjection1 = self.stackDistance.phaseProjection
                
            
            if self.importXYPosition: #if xy position file is imported, do not perform thresholding
                #threshold_filter is zeros    
                self.stackDistance.threshold_filter =\
                                            np.zeros(np.shape(self.stackDistance.phaseProjection1))
                #xy position for each hologram is read from the text file
                self.ReadXYPositionFromDonner(self.stackDistance.holoNumber)
            else:#if xy position is not imported, we have to compute if from threshold
                self.stackDistance.Threshold()#computed threshold mask depending of the method
                #Extract xy positions of particules from thresholded image
                self.stackDistance.ExtractXYPositionFromHolo()
            self.stackDistance.zFromxy()#compute z position for each xy detected position
            #if derivate is used, phaseProjection1 becomes phaseProjection0
            self.stackDistance.phaseProjection0 = self.stackDistance.phaseProjection1
        toc = time.clock()
        timeremaining = (toc-tic)*self.last_holo_seq-self.stackDistance.first_holo_seq
        progressbar.SetTimeRemaining(timeremaining)
        #successive reconstruction of holograms to determine x,y,z position for all particules
        while self.stackDistance.holoNumber < self.last_holo_seq-1:
            tic = time.clock()
            time.sleep(0.05)
            progressbar.progressbar.setValue(self.stackDistance.holoNumber-\
                                        self.stackDistance.first_holo_seq)
            QtGui.qApp.processEvents()
            if not progressbar._active:#stop if stop button pressed in progress bar
                break
            self.stackDistance.holoNumber += 1 #increment the hologram number
            self.stackDistance.AcquireAndSaveStack()
            if not SaveOnlyStack:
                self.stackDistance.ZProjection()
                self.stackDistance.phaseProjection1 = self.stackDistance.phaseProjection
                
                if self.importXYPosition:
                    self.stackDistance.threshold_filter = \
                                        np.zeros(np.shape(self.stackDistance.phaseProjection0))
                    self.ReadXYPositionFromDonner(self.stackDistance.holoNumber)
                else:
                    self.stackDistance.Threshold()
                    self.stackDistance.ExtractXYPositionFromHolo()
                self.stackDistance.zFromxy()
                self.stackDistance.phaseProjection0 = self.stackDistance.phaseProjection1
            toc = time.clock()
            timeremaining = (toc-tic)*(self.last_holo_seq-1-self.stackDistance.holoNumber)
            progressbar.SetTimeRemaining(timeremaining)

        progressbar._active = False   
        self.SaveStackParams()#save the stack parameters 
        if not SaveOnlyStack:
            #save all coordinates x,y,z volume mean amplitude WITHOUT tracking consideration!
            self.SaveCoordinates()
            if not self.importXYPosition:
                #track particule from first holo  partby minimizing dist. to part. on the next holo
                self.Computexyzposition()
            #init to state before begin tracking
        self.stackDistance.holoNumber = self.stackDistance.first_holo_seq
        self.stackDistance.phaseProjection = None
        self.stackDistance.phaseProjection0 = None
        self.stackDistance.phaseProjection1 = None
        
    def XYFileAndNumberHoloCompatible(self):
        #verifiy if XY text file is compatible with the number of holo (works only for a customer)
        NumberHolosSeq = self.last_holo_seq-self.stackDistance.first_holo_seq
        NumberHolosDonner = self.NumberHolosSeqDonner()
        compatible = True
        if NumberHolosSeq != NumberHolosDonner:
            utils.ErrorMessage('First holo and last holo are not compatbile with the XY file\nXY File contains '+str(NumberHolosDonner)+' frames and your sequence has '+ str(NumberHolosSeq)+' frames.')
            compatible = False
        return compatible
        
    def NumberHolosSeqDonner(self):
        #Define number of holo of customer xy text file
        with open(self.XYImportFilePath, 'rU') as f:
            reader = csv.reader(f, delimiter='\t')
            first_row = next(reader)
            identical_text = first_row[0]
            identical_text = identical_text[:6]
            num_cols = len(first_row)
        return num_cols/3
        
    def ReadXYPositionFromDonner(self, hologramNumber):
        #Read xy position from customer text file 
        xyPosition = []
        with open(self.XYImportFilePath, 'rU') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in islice(reader, 1, None):
                x = row[3*hologramNumber]
                y = row[3*hologramNumber+1]
                if x != '' and y != '':
                    xpos = (float)(x)
                    ypos = (float)(y)
                    coord = (ypos, xpos)
                    xyPosition.append(coord)
        self.stackDistance.xy_position = xyPosition
        self.stackDistance.ProjectionWithDetectedPoints()
    def SaveCoordinatesParticule(self, firstImage, particulenumber, zconvert, x, y, z, vol, amp):
        #Save Coordinates for a single particule from automatic tracking
        f = open(self.stackDistance.directoryPath+'\\coordinates_'+str(particulenumber).zfill(3)\
                            +'.txt', 'w')
        f.writelines('holoNumber\txpos[um]\typos[um]\tzpos[um]\tdist_rec[cm]\tvolume[um3]\tmeanamp\n')
        n_coord = np.shape(x)[0]    
        for k in range(n_coord):
            timeholo = str(firstImage+k).zfill(5)
            xpos = str(x[k])
            ypos = str(y[k])
            zpos = str(z[k])
            dist_rec = str(z[k]/zconvert+self.stackDistance.central_rec_dist)
            volume = str(vol[k])
            meanamp = str(amp[k])
            f.writelines(timeholo+'\t'+xpos+'\t'+ypos+'\t'+zpos+'\t'+dist_rec+'\t'+volume+\
                            '\t'+meanamp+'\n')
        f.close()
         
    def ReadCoordinatesParticule(self, partNumber):
        #Read coordinated for a given particule for single particule tracking (or loop on particule)
        with open(self.stackDistance.directoryPath+'\\SingleParticuleMethod\\Part'+\
                        str(partNumber).zfill(3)+'\\px_py_d_vol_meanamp.txt', 'r') as f:
            reader = csv.reader(f)
            Coordinates = []
            for row in islice(reader, 1, None):
                line = row[0].split()
                X1 = [(float)(line[0]), (float)(line[1]), (float)(line[2]), (float)(line[3]), 
                      (float)(line[4]), (float)(line[5])]
                Coordinates.append(X1)    
        return np.array(Coordinates)
    
    def ReadCoordinates(self):
        #Read Coordinates of all particule (multiple particles for all holograms)
        with open(self.stackDistance.directoryPath+'\\coordinates.txt', 'r') as f:
            reader = csv.reader(f)
            Coordinates = []
            for row in islice(reader, 1, None):
                line = row[0].split()
                X1 = [(float)(line[0]), (float)(line[1]), (float)(line[2]), 
                      (float)(line[3]), (float)(line[4]), (float)(line[5])]
                Coordinates.append(X1)
    
        return Coordinates
    
    def ExportModifiedDonner(self):
        #Export Modified text customer xy file
        coordinatesPath = self.stackDistance.directoryPath+'\\coordinates.txt'
        DonnerPath = self.XYImportFilePath
        newFilePath = self.XYImportFilePath[:-4]+'_z.txt'
        coord = []
        zconvert = 1e4/self.stackDistance.MO**2*\
                    self.stackDistance.index_medium/self.stackDistance.index_chamber
        with open(coordinatesPath, 'r') as c:
            reader = csv.reader(c, delimiter='\t')
            for row in islice(reader, 1, None):
                holoNumber = (int)(row[0])
                z = (float)(row[3])
                z *= zconvert
                coord.append((holoNumber, z))
          
        newcoordDonner = []
        with open(DonnerPath, 'rU') as f:
            reader = csv.reader(f, delimiter='\t')
            first_row = next(reader)
            identical_text = first_row[0]
            identical_text = identical_text[:6]
            num_cols = len(first_row)
        numberFrames = num_cols/3
        with open(DonnerPath, 'rU') as d:
            reader = csv.reader(d, delimiter='\t')
            p = 0
            for k in range(numberFrames):
                first_row.insert(4*k+3, identical_text+str(k).zfill(3)+'z')
            newcoordDonner.append(first_row)
            for row in islice(reader, 1, None):          
                for k in range(numberFrames):
                    coord_frame = [i for i in coord if i[0] == k]
                    #print np.shape(coord_frame)
                    if p < np.shape(coord_frame)[0]:
                        value = str(coord_frame[p][1])
                    else:
                        value = ''
                    row.insert(4*k+3, value)
                p += 1
                newcoordDonner.append(row)
       
        f = open(newFilePath, 'w')
        for row in newcoordDonner:
            line = ''
            for value in row:
                line += value+'\t'
            line = line[:-1]+'\n'
            f.writelines(line)
        f.close()

    def ComputexyzpositionSingleParticuleMethod(self, partNumber):
        #plot 3d graph of single particule from single particule tracking (or loop on particuler)
        
        #z conversion according to magnification and refractive index 
        zconvert = 1e4/self.stackDistance.MO**2*\
            self.stackDistance.index_medium/self.stackDistance.index_chamber
        Coordinates = self.ReadCoordinatesParticule(partNumber) #Read the coordinates
        xpos = Coordinates[:, 1]*self.stackDistance.pxSize #convert px coord to x position in micron
        ypos = Coordinates[:, 2]*self.stackDistance.pxSize#convert px coords to y position in micron
        zpos = Coordinates[:, 3]*zconvert #convert rec distance coordinates to z position in micron
        
        #plot figure
        hf = plt.figure(30)
        ha = hf.add_subplot(111, projection='3d')
        ha.plot(xpos, ypos, zpos)
        ha.set_xlabel('X Position [um]')
        ha.set_ylabel('Y Position [um]')
        ha.set_zlabel('Z Position [um]')
        
        #axis can be modified if necessary
        #ha.set_zlim([-self.stackDistance.range_rec_dist*zconvert,self.stackDistance.range_rec_dist*zconvert])
        #ha.set_xlim([0.,1024.*self.stackDistance.pxSize])
        #ha.set_ylim([0.,1024.*self.stackDistance.pxSize])
        
        #Save Plot
        plt.savefig(self.stackDistance.directoryPath+'\\SingleParticuleMethod\\Part'+\
                        str(partNumber).zfill(3)+'\\3DGraph.png')
        plt.show() 
        #2D image
        directorySave = self.stackDistance.directoryPath+'\\SingleParticuleMethod\\Part'+\
                        str(partNumber).zfill(3)+'\\xy_trackImage\\'
        directoryOpen = self.stackDistance.directoryPath+'\\SingleParticuleMethod\\Part'+\
                        str(partNumber).zfill(3)+'\\Phase\\'
        utils.CreateDirectory(self.stackDistance.directoryPath+\
                    '\\SingleParticuleMethod\\Part'+str(partNumber).zfill(3), 'xy_trackImage')
        k = 0    
        #Save focused phase image with xy detected position with a red point
        for holo in Coordinates[:, 0]:
            phpath = directoryOpen+str((int)(holo)).zfill(5)+'_phase.bin'
            ph, header = binkoala.read_mat_bin(phpath)
            ph = (ph-np.min(ph))/(np.max(ph)-np.min(ph))*255
            rgbArray = np.zeros((np.shape(ph)[0], np.shape(ph)[1], 3), 'uint8')
            rgbArray[:, :, 0] = ph
            rgbArray[:, :, 1] = ph
            rgbArray[:, :, 2] = ph
            x = (int)(Coordinates[k, 1])
            y = (int)(Coordinates[k, 2])
            rgbArray[x-2:x+2, y-2:y+2, 0] = 255
            rgbArray[x-2:x+2, y-2:y+2, 1] = 0
            rgbArray[x-2:x+2, y-2:y+2, 2] = 0
            im = PIL.Image.fromarray(rgbArray)
            im.save(directorySave+str((int)(holo)).zfill(5)+'_phase_xy.jpg')
            k += 1
            
    
    def Computexyzposition(self):
        #compute trajectories from automatic tracking from coordinates.txt (all coordinates)
    
        #z conversion according to MO and refractive index
        zconvert = 1e4/self.stackDistance.MO**2*\
            self.stackDistance.index_medium/self.stackDistance.index_chamber
        Coordinates = self.ReadCoordinates()#Read all coordinates
        firstImage = (int)(Coordinates[0][0])#Define first image (=first holo)
        lastImage = (int)(Coordinates[-1][0])#Define last image (=lastholo)
        #Create directory to save 3D Graph
        utils.CreateDirectory(self.stackDistance.directoryPath, '3DGraph')
        ## follow particles
        #first xy positionS (several particule)
        particulesFrame0 = [row for row in Coordinates if row[0] == firstImage]
        #define number of particules of the first frame (first holo)
        numberOfParticules = np.shape(particulesFrame0)[0]
        colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
        gf = plt.figure(20)
        ga = gf.add_subplot(111, projection='3d')
        #track each particule defined on the first image
        for p in range(numberOfParticules):
            part_k_minus_1 = particulesFrame0[p] #position of particule p in the first frame
            #initialize coordinates for particule with number of hologram size
            tpart = np.zeros(lastImage-firstImage+1)#time (in term of holo number)
            xpart = np.zeros(lastImage-firstImage+1)#x position
            ypart = np.zeros(lastImage-firstImage+1)#y position
            zpart = np.zeros(lastImage-firstImage+1)#z position
            volpart = np.zeros(lastImage-firstImage+1)#volume
            amppart = np.zeros(lastImage-firstImage+1)#mean amplitude
            #write coordinates for the initial holo
            xpart[0] = part_k_minus_1[1]
            ypart[0] = part_k_minus_1[2]
            zpart[0] = part_k_minus_1[3]
            volpart[0] = part_k_minus_1[4]
            amppart[0] = part_k_minus_1[5]
            #previous position is defined by Pt1
            Pt1 = part_k_minus_1 
            #vol for part at position Pt1 (vol can be a criteria to select next position)
            volume_Pt = Pt1[4]
            #track the particule iteratively from the second holo to the last one
            #by searching for the closest detected particule on the next holo
            for k in range(1, lastImage-firstImage+1):
                #Extract all particule coordinates there are in the holo number=firstImage+k
                particulesFramek = [row for row in Coordinates if row[0] == (float)(firstImage+k)]
                #Find the Nearest point of Pt1 (holo i) in the list of particules of holo i+1
                #and the distance (dist) between Pt1 and the closest point
                [index_nearest, dist] = self.FindNearestPoints(Pt1, particulesFramek, volume_Pt)
                
                #if the nearest detected particule is smaller than the maximal 
                #displacement defined by the user, the new position if defined
                #by coordinatess of the closest particule
                if dist < self.max_displacement*self.stackDistance.pxSize and index_nearest != -1:
                    #certainly particule not detected
                    Pt_k = particulesFramek[index_nearest]
                #if the distance is larger, we assume that the particule was 
                #not detected (for example if derivation was used and the particule 
                #did not move) or if threshold fails. In this case, we assume
                #the particule did not move and next position is the same than 
                #previous
                else:
                    Pt_k = Pt1
        
                Pt1 = Pt_k #new previous position is the detected one
                #write the new coordinates for the tracked particule
                tpart[k] = k
                xpart[k] = Pt_k[1]
                ypart[k] = Pt_k[2]
                zpart[k] = Pt_k[3]
                volpart[k] = Pt_k[4]
                amppart[k] = Pt_k[5]
            ##computation of the physical position
            xpos = xpart*self.stackDistance.pxSize #conversion from piyel to x position in um
            ypos = ypart*self.stackDistance.pxSize #conversion from piyel to y position in um
            zpos = zpart*zconvert #conversion of reconstruction distance in um with zconvert
            
            #plot trajectories for single particule
            hf = plt.figure(30)
            ha = hf.add_subplot(111, projection='3d')
            
            ha.plot(xpos, ypos, zpos)
            ha.set_xlabel('X Position [um]')
            ha.set_ylabel('Y Position [um]')
            ha.set_zlabel('Z Position [um]')
            #force zlimit to be between +-range rec distance converted in um
            ha.set_zlim([-self.stackDistance.range_rec_dist*zconvert, 
                         self.stackDistance.range_rec_dist*zconvert])
            #force xlimit and ylimit to be between 1024 pixels*pixel size in um
            ha.set_xlim([0., 1024.*self.stackDistance.pxSize])
            ha.set_ylim([0., 1024.*self.stackDistance.pxSize])
            #save figure
            plt.savefig(self.stackDistance.directoryPath+'\\3DGraph\\3d_'+str(p).zfill(5)+'.png')
            plt.show()  
            
            #plot 3D trajectories for each tracked particule in the same graph
            ga.plot(xpos, ypos, zpos, color=colors[p%7])
            if p == numberOfParticules-1:
                ga.set_xlabel('X Position [um]')
                ga.set_ylabel('Y Position [um]')
                ga.set_zlabel('Z Position [um]')
                ga.set_zlim([-self.stackDistance.range_rec_dist*zconvert, 
                             self.stackDistance.range_rec_dist*zconvert])
                ga.set_xlim([0., 1024.*self.stackDistance.pxSize])
                ga.set_ylim([0., 1024.*self.stackDistance.pxSize])
                plt.figure(20)
                plt.savefig(self.stackDistance.directoryPath+'\\3DGraph\\3d_all.png')
        
            ##save coordinates particles in a text file
            self.SaveCoordinatesParticule(firstImage, p, zconvert, xpos, ypos, zpos, volpart, 
                                          amppart)
            
            #save 2D image of the Projection phase with red point on the measured position
            utils.CreateDirectory(self.stackDistance.directoryPath, 
                                               '\\Follow_part_proj\\Part'+str(p).zfill(3))
            for k in range(0, lastImage-firstImage+1):
                path = self.stackDistance.directoryPath+'\\Projection\\Proj_'+\
                        str(firstImage+k).zfill(5)+'.jpg'
                img = np.array(PIL.Image.open(path))
                x = (int)(xpart[k])
                y = (int)(ypart[k])
                img[x-2:x+2, y-2:y+2, 0] = 255
                img[x-2:x+2, y-2:y+2, 1] = 0
                img[x-2:x+2, y-2:y+2, 2] = 0
                im = PIL.Image.fromarray(img)
                im.save(self.stackDistance.directoryPath+'\\Follow_part_proj\\Part'+\
                            str(p).zfill(3)+'\\Proj_'+str(firstImage+k).zfill(5)+'.jpg')
                
    def FindNearestPoints(self, Pt, ListPoint, vol_Pt, K=1, limit_dist_2_part=0, 
                          distanceWithZ=False):
#        method to determine the closest point to Pt (previous position) in the list of point of actual holo
#        K is the number of closest particule to detect, by default K=1, but K>1
#        can be used, but another criteria has to be used to choose between the 
#        closest points if these points are closer thant limit_dist_2_part
#        distanceWithZ=True if we want to take in accont the zpostion for the computation
#        of distance. By default is false because error on z position evaluation could
#        induce error on particle tracking
        pxsize = self.stackDistance.pxSize #reconstructed pixel size
        #zconvert used to convert reconstruction distance to position distance
        zconvert = 1e4/self.stackDistance.MO**2*\
                    self.stackDistance.index_medium/self.stackDistance.index_chamber
        #mask used to convert x[pixel],y[pixels],z[rec dist] to x[um],y[um],z[um]
        mask = np.array([pxsize, pxsize, zconvert])
        ndata = np.shape(ListPoint)[0] #number of particules
        if ndata != 0: #if the number of part is not zero (if 0, no particule found, and dist=1e10)
            #compute the distance between Pt1 and all particules
            D = []
            Point = []
            for k in range(ndata):
                if not distanceWithZ:
                    D.append(ListPoint[k][1:3]*mask[0:2]) #convert points position in um
                    Point.append(Pt[1:3]*mask[0:2])  #convert previous point in um
                else:
                    D.append(ListPoint[k][1:4]*mask[0:3])
                    Point.append(Pt[1:4]*mask[0:3])
            #compute the distance between points
            a = np.asarray(D)
            b = np.asarray(Point)
            ndata = np.shape(D)[0]
            K = K if K < ndata else ndata
            sqd = np.sqrt(((a-b)**2).sum(axis=1))
            
            #sort the distance
            idx = np.argsort(sqd)             
            index_part = idx[0]
            dist = np.abs(sqd[idx[:K]][0])
            
            #if two particules are detected, filter with volume
            #not activated if K=1 (detect closest particule)
            if K > 1 and np.shape(idx)[0] > 1:
                dist_two_part = np.abs(sqd[idx[:K]][0]-sqd[idx[:K]][1])
                if dist_two_part < limit_dist_2_part:
                    print 'two particiles are close'
                    vol_Pt = Pt[4]
                    vol = []  
                    for index in idx[:K]:
                        vol.append(ListPoint[index][4])
                    diff_vol = np.abs((np.asarray(vol)-vol_Pt)/vol_Pt)
                    index_part = diff_vol.argmin()
                    index_part = idx[index_part]
        else:
            index_part = -1 #no particle found 
            dist = 1e10
        return [index_part, dist]
        
    def SaveStackParams(self):
        #Save stack parameters used for the computation (to perform tracking from coordinates)
        f = open(self.stackDistance.directoryPath+'\\rec_dist_params.txt', 'w')
        f.writelines('central rec dist_[cm]\trange rec dist [cm]\tstep rec dist [cm]\tpxSize [um]\tMO\tindex medium\tindex_chamber\n')
        f.writelines(str(self.stackDistance.central_rec_dist)+\
        '\t'+str(self.stackDistance.range_rec_dist)+'\t'+str(self.stackDistance.step_rec_dist)+\
        '\t'+str(self.stackDistance.pxSize)+'\t'+str(self.stackDistance.MO)+'\t'+\
        str(self.stackDistance.index_medium)+'\t'+str(self.stackDistance.index_chamber)+'\n')
        f.close()
    def OpenStackZparams(self): 
        #Open stack parameters in cas of Par. Tracking from Coordinates
        range_rec_dist = 0.
        step_rec_dist = 0.
        central_rec_dist = 0.
        pxSize = 1.
        with open(self.stackDistance.directoryPath+'\\rec_dist_params.txt', 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            a = 0
            for row in reader:
                if a != 0:
                    #line=row[0].split()
                    central_rec_dist = float(row[0])
                    range_rec_dist = float(row[1])
                    step_rec_dist = float(row[2])
                    pxSize = float(row[3])
                a += 1
        
#        print central_rec_dist
#        print range_rec_dist
#        print step_rec_dist
        self.stackDistance.central_rec_dist = central_rec_dist
        self.stackDistance.range_rec_dist = range_rec_dist
        self.stackDistance.step_rec_dist = step_rec_dist
        self.stackDistance.pxSize = pxSize

            
class ProgressBar(QtGui.QWidget):
    #class for the Progress bar
    def __init__(self, total, parent=None):
        super(ProgressBar, self).__init__(parent)
        self.progressbar = QtGui.QProgressBar()
        self.progressbar.setMinimum(0)
        self.progressbar.setMaximum(total)
        self.button = QtGui.QPushButton('Stop')
        self.button.clicked.connect(self.handleButton)
        self.timeremaining = 0.
        self.text = QtGui.QLabel('')
        main_layout = QtGui.QGridLayout()
        main_layout.addWidget(self.button, 0, 0)
        main_layout.addWidget(self.progressbar, 0, 1)
        main_layout.addWidget(self.text, 1, 0)
        self.setLayout(main_layout)
        self.setWindowTitle('Progress')
        self._active = True

    def SetTimeRemaining(self, timeremaining):
        #set the time remainin in the progressbar
        self.timeremaining = timeremaining
        self.updatetimeremaining()
        
    def updatetimeremaining(self):
        #update the time remaininng
        text = time.strftime("%H:%M:%S", time.gmtime(self.timeremaining))
        self.text.setText(text)
    def handleButton(self):
        if self._active:
            self._active = False
            
    def closeEvent(self, event):
        self._active = False
        
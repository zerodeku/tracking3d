# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 10:10:46 2014

@author: tcolomb
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Aug 05 14:27:06 2014

@author: tcolomb
"""

import os
import time
import numpy as np
import clr
import System
from System import Array
from Tkinter import *
from guidata.qt.QtGui import *
import PIL.Image as Image
import pylab
import matplotlib.pyplot as plt
import cv2


## to use external .Net library for remote control
pythonLoc = os.path.dirname(os.path.abspath(__file__))
clr.AddReference(pythonLoc+"\\TCPClient")
import KoalaClient
obj = KoalaClient.KoalaTCPClient()

#def obj_MsgBoxEvent(source,args):
#    print 'Error Koala'
#    return true
#obj.MsgBoxEvent+=KoalaClient.KoalaTCPClient.MsgBoxEventHandler(obj_MsgBoxEvent)

#%% Remote functions

def getPhaseImage(contrast, width, height):#Get the phase image 8bit from contrast 
    if contrast == 'lambda1':
        obj.SelectDisplayWL(8192)#select display
    elif contrast == 'lambda2':
        obj.SelectDisplayWL(16384)
    elif contrast == 'synth':
        obj.SelectDisplayWL(32768)
    elif contrast == 'synth_short':
        obj.SelectDisplayWL(65536)
    elif contrast == 'mapping':
        obj.SelectDisplayWL(131072)
    else:
        obj.SelectDisplayWL(8192)
                   
    w = width
    h = height
    stride = (w/4)*4
    if w % 4 != 0:
        stride += 4
    size = stride*h
    
    bufP = Array.CreateInstance(System.Byte, size)
    obj.GetPhaseImage(bufP)
    
    floatbufP = np.zeros(size, dtype=int)
    for x in range(0, size-1):
        floatbufP[x] = int(bufP[x])
    out = np.reshape(floatbufP, (h, stride))
    return out
    
def getIntensityImage(contrast, width, height):#Get the intensity 8bit image
    if contrast == 'lambda1':
        obj.SelectDisplayWL(8192)
    elif contrast == 'lambda2':
        obj.SelectDisplayWL(16384)
    elif contrast == 'synth':
        obj.SelectDisplayWL(32768)
    elif contrast == 'synth_short':
        obj.SelectDisplayWL(65536)
    elif contrast == 'mapping':
        obj.SelectDisplayWL(131072)
    else:
        obj.SelectDisplayWL(8192)
            
        
    w = width
    h = height
    stride = (w/4)*4
    if w % 4 != 0:
        stride += 4
    size = stride*h
    
    bufP = Array.CreateInstance(System.Byte, size)
    obj.GetIntensityImage(bufP)
    
    floatbufP = np.zeros(size, dtype=int)
    for x in range(0, size-1):
        floatbufP[x] = int(bufP[x])
    out = np.reshape(floatbufP, (h, stride))
    return out
    
def getIntensity32fImage(contrast, width, height): #get float intensity image
    if contrast == 'lambda1':
        obj.SelectDisplayWL(2048)
    elif contrast == 'lambda2':
        obj.SelectDisplayWL(4096)
    else:
        obj.SelectDisplayWL(2048)
    
    size = width*height
    
    bufP = Array.CreateInstance(System.Single, size)
    obj.GetIntensity32fImage(bufP)
    floatbufP = np.zeros(size, dtype=float)
    for x in range(0, size-1):
        floatbufP[x] = float(bufP[x])
    out = np.reshape(floatbufP, (height, width))
    return out
    
def getPhase32fImage(contrast, width, height):#get float phase image
    if contrast == 'lambda1':
        obj.SelectDisplayWL(8192)
    elif contrast == 'lambda2':
        obj.SelectDisplayWL(16384)
    elif contrast == 'synth':
        obj.SelectDisplayWL(32768)
    elif contrast == 'synth_short':
        obj.SelectDisplayWL(65536)
    elif contrast == 'mapping':
        obj.SelectDisplayWL(131072)
    else:
        obj.SelectDisplayWL(8192)
    
    size = width*height
    
    bufP = Array.CreateInstance(System.Single, size)
    obj.GetPhase32fImage(bufP)
    floatbufP = np.zeros(size, dtype=float)
    for x in range(0, size-1):
        floatbufP[x] = float(bufP[x])
    out = np.reshape(floatbufP, (height, width))
    return out

def SavePhaseFloatToFile(contrast, fname): #save float phase image in bin format
    if contrast == 'lambda1':
        obj.SelectDisplayWL(8192)
    elif contrast == 'lambda2':
        obj.SelectDisplayWL(16384)
    elif contrast == 'synth':
        obj.SelectDisplayWL(32768)
    elif contrast == 'synth_short':
        obj.SelectDisplayWL(65536)
    elif contrast == 'mapping':
        obj.SelectDisplayWL(131072)
    else:
        obj.SelectDisplayWL(8192)
  
    winID = 4
    obj.SaveImageFloatToFile(winID, fname, True)#True mean .bin format
    
def SaveAmpFloatToFile(contrast, fname):#save float intensity image in bin format
    if contrast == 'lambda1':
        obj.SelectDisplayWL(2048)
    elif contrast == 'lambda2':
        obj.SelectDisplayWL(4096)
    else:
        obj.SelectDisplayWL(2048) 
    winID = 2
    obj.SaveImageFloatToFile(winID, fname, True)#True mean .bin format
    

def OpenWindows(): #Open all windows
    obj.OpenPhaseWin()
    obj.OpenIntensityWin()
    obj.OpenHoloWin()
    
def OpenHologram(*args):
    #open hologram from path defined, if path not defined, user is asked to choose a file
    if len(args) == 1:
        holopath = unicode(QFileDialog.getOpenFileName())
        obj.LoadHolo(holopath, args[0])
    elif len(args) == 2:
        obj.LoadHolo(args[0], args[1])
    
def OpenConfig(configNumber):#Open configuration
    obj.OpenConfig(configNumber)
    time.sleep(0.2)
    
def Connect(IP, user):#Connect the remote according to IP
    return obj.Connect(IP, user, True)  
    
def Login(password): #Login with the password
    obj.Login(password)
    
def Logout():
    obj.Logout()
    
def GetPxSizeUm(): #Get reconstructed pixel size in um
    return obj.GetPxSizeUm()
    
def GetPhaseWidth():
    return obj.GetPhaseWidth()
    
def GetPhaseHeight(): 
    return obj.GetPhaseHeight()

def GetLambdaNm(wavelength): #return the wavelength of source 1,2 or 3
    return obj.GetLambdaNm(wavelength, False)
def AddCorrSegmentAndComputePhaseCorrection(hv_seg, degree):
    #Add correction segment and compute phase correction
    obj.AddCorrSegment(hv_seg.top, hv_seg.left, hv_seg.length, hv_seg.orient)
    obj.ComputePhaseCorrection(1, degree)
    
def GetRecDistCM():
    return obj.GetRecDistCM()
    
def SetRecDistCM(d):
    obj.SetRecDistCM(d)
    obj.OnDistanceChange()
    
def SetUnwrap2DState(state):
    obj.SetUnwrap2DState(state)
    
def SetUnwrap2DMethod(method):
    obj.SetUnwrap2DMethod(method)
    
def Unwrap2D(method):
    obj.SetUnwrap2DMethod(1) #method for unwrapping
    obj.SetUnwrap2DState(True)
    
def ResetCorrSegments():
    obj.ResetCorrSegment()

def AlgoResetPhaseMask():
    #does not work
    obj.AlgoResetPhaseMask()
        
def SaveImageToTiff(winID, path, width, height, *args):
    w = width
    h = height
    bufP = Array.CreateInstance(System.Single, h*w)
    if winID == 2:
        obj.GetIntensity32fImage(bufP)
    elif winID == 4:
        obj.GetPhase32fImage(bufP)
    
    img = np.zeros(h*w, dtype=float)    
    for x in range(0, h*w-1):
        img[x] = float(bufP[x])
    img = np.reshape(img, (h, w))
    if len(args) == 0:
        minVal = np.min(img)
        maxVal = np.max(img)
    img = (img-minVal)/(maxVal-minVal)
    c = np.around(img*255.)
    d = np.maximum(0., c)
    c = np.minimum(255., d)
    d = c.astype(np.uint8)
    p = Image.fromarray(d)
    p.save(path)
def SaveImageToTiffThroughKoala(winID, path):
    obj.SaveImageToFile(winID, path)
    
    

#Method not used in 3D Tracking!!
def CorrectionTilt(contrast): #Perform tilt correction on a given contrast (not used in 3D tracking)
    for k in range(0, 2):
        [x_line1, y_line1] = pylab.ginput(2)
        r = (x_line1[0], x_line1[1], y_line1[0], y_line1[1])
        dx = abs(r[0] - r[2])
        dy = abs(r[1] - r[3])
        if dx > dy: # horizontal
            orient = 0
            left = int(min(r[0], r[2]))
            top = int(r[1])
            length = int(abs(r[0] - r[2]))
            #x_draw=(r[0],r[2])
            #y_draw=(r[1],r[1])
        else: # vertical
            orient = 1
            top = int(min(r[1], r[3]))
            left = int(r[0])
            length = int(abs(r[1] - r[3]))
            #x_draw=(r[0],r[0])
            #y_draw=(r[1],r[3])
        obj.AddCorrSegment(top, left, length, orient)
        #plt.plot(x_draw, y_draw)# Plot the first line
    obj.ComputePhaseCorrection(1, 1)
    #PhImg=getPhImg(contrast)
    #plt.imshow(PhImg)

def CorrectionAberration(a, b):
    #Perform higher order correction on a given contrast (not used in 3D tracking)
    for k in range(0, 2):
        [x_line1, y_line1] = pylab.ginput(2)
        r = (x_line1[0], x_line1[1], y_line1[0], y_line1[1])
        dx = abs(r[0] - r[2])
        dy = abs(r[1] - r[3])
        if dx > dy: # horizontal
            orient = 0
            left = int(min(r[0], r[2]))
            top = int(r[1])
            length = int(abs(r[0] - r[2]))
#            x_draw=(r[0],r[2])
#            y_draw=(r[1],r[1])
        else: # vertical
            orient = 1
            top = int(min(r[1], r[3]))
            left = int(r[0])
            length = int(abs(r[1] - r[3]))
            #x_draw=(r[0],r[0])
            #y_draw=(r[1],r[3])
        obj.AddCorrSegment(top, left, length, orient)
        #plt.plot(x_draw, y_draw)# Plot the first line
    obj.ComputePhaseCorrection(a, b)

def ComputePhaseCorrection(method, degree): #compute phase correction according to method and degree
    obj.ComputePhaseCorrection(method, degree)

def ExtractProfile(): #extract profile from image: not used in 3D Tracking
    [x, y] = ginput(2)
    x_draw = (x[0], y[0])
    y_draw = (x[1], y[1])
    plt.plot(x_draw, y_draw)
    obj.SetPhaseProfileState(True)
    obj.ExtractPhaseProfile(x[0], x[1], y[0], y[1])
    length = obj.GetPhaseProfileLength()
    bufP = Array.CreateInstance(System.Double, length)
    obj.GetPhaseProfile(bufP)
    floatbufP = np.zeros(length, dtype=float)
    for x in range(0, length-1):
        floatbufP[x] = float(bufP[x])
    return floatbufP  
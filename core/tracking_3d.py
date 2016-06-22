# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 13:15:23 2016

@author: zerodeku
"""

from tracking3d.core import remote_Koala
from tracking3d.core import SqliteDB
import time, os
from tracking3d.core import utilsForProcessing as utils
from guidata.qt import QtGui
import numpy as np
from PIL import Image
import cv2

class Tracking3D:
    tracking = None
    def __init__(self):
        self.remoteKoala = None
        self.db = None
        # remote parameters
        self.contrast = ['lambda1', 'mapping']
        self.centRecDist = -1.7 # initial reconstruction distance
        self.width = 750 # image width
        self.height = 750 # image height
        self.pixSize = 1 # real value is updated when holograph is loaded
        self.KoalaViewSize = 748 # Koala software phase image view size
        self.dof = None # depth of field
        self.ratioDOF = None
        self.NA = None
        self.MO = None
        self.CCDPixelSizeUM = None
        self.zConvert = None
        self.stepUM = None
        self.stackHeightUM = None
        self.rangeRecDist = None
        self.tiltAjust = True
        self.adjacentSubtraction = True
        # sample related parameters
        self.chamberIndex = None
        self.mediumIndex = None
        
        # tracking parameters and data
        self.phaseBcg = None
        self.meanPhaseIntensity = None
        self.numHolos = 0
        self.reconsDist = 0
        self.threshold = 40
        self.directory = None
        self.bcgRemove = True
        self.tempDir = None
        self.tempPhaseDir = None
        self.tempPhaseBcgRemovedDir = None 
        self.tempPhaseStackDir = None
        self.zSearchDX = 20
        self.zsearchDY = 20
        self.nextSearchDX = 0
        self.nextSearchDY = 0
        
        
    @classmethod # RemoteKoala singleton method
    def getInstance(cls):
        if cls.tracking == None:
            cls.tracking = Tracking3D()
        return cls.tracking
    
    # initialize remoteKoala
    def initRemoteKoala(self, data):
        db = SqliteDB.DBReader()
        db.ProductionNo = data.MODB
        db.ParamFromChoice()
        self.db = db
        self.directory = db.DefaultDirectory
        self.remoteKoala = remote_Koala.RemoteKoala.from_dataset(self.db, data)
        self.remoteKoala.DefineSeqFormat()
        self.numHolos = self.remoteKoala.framesNumber
        self.tempDir = utils.CreateDirectory(self.directory, 'tmp')
        
    # connect Koala and open windows
    def setupKoala(self):
        self.remoteKoala.ConnectKoala()
        self.remoteKoala.OpenConfig()
        self.remoteKoala.OpenHologram()
        time.sleep(0.1)
        self.reconsAdjust()
        self.remoteKoala.GetLambdaNm()
        self.updatePhaseImage()
    
    def initSampleParams(self, data):
        self.mediumIndex = data.index_medium
        self.chamberIndex = data.index_chamber
        self.sampleIndex = data.index_sample
        self.particleSize = data.particle_size
        self.particleMaxSpeed = data.max_speed
    
    def initTrackingParams(self, data):
        self.rangeRecDist = data.range_rec_dist
        self.ratioDOF = data.ratio_DOF
        self.bcgRemove = data.bcgRemove
        self.centRecDist = data.central_rec_dist
        
        self.NA = self.db.NA
        self.CCDPixelSizeUM = self.db.CCDPixelSizeUM
        self.pxSize = self.remoteKoala.pxSize
        self.dof = self.remoteKoala.wavelength / (self.NA**2) * 1e-7
        self.MO = self.CCDPixelSizeUM / self.pxSize
        self.stepRecDist = self.MO**2 * self.dof / self.mediumIndex \
                               / self.ratioDOF * self.chamberIndex
        
        self.zConvert = 1e4 / self.MO**2 * self.mediumIndex / self.chamberIndex
#        self.stepUM = self.stepRecDist * self.zConvert
#        self.stackHeightUM = self.rangeRecDist * 2 * self.zConvert
        self.stepUM = 1
        self.stackHeightUM = 70
        
        
    # basic ajustments for reconstruction
    def reconsAdjust(self):
        self.remoteKoala.remoteCommands.SetRecDistCM(self.centRecDist)
        self.remoteKoala.remoteCommands.SetUnwrap2DState(True)
        
        # define some line segments for tilt correction
        if self.tiltAjust:
            hvSegs = []
            width, height = self.KoalaViewSize, self.KoalaViewSize
            margin = 50
            hvSegs.append(HVSeg(0, margin, height, 1))
            hvSegs.append(HVSeg(0, width - margin, height, 1))
            hvSegs.append(HVSeg(margin, 0, width, 0))
            hvSegs.append(HVSeg(height - margin, 0, width, 0))
            hvSegs.append(HVSeg(height / 2, 0, width, 0))
            hvSegs.append(HVSeg(0, width / 2, height, 1))
            for hvSeg in hvSegs:
                self.remoteKoala.new_hv_seg(hvSeg)

    def getPhaseImageByHoloNum(self, holoNum):
        self.remoteKoala.first_holo_seq = holoNum
        self.remoteKoala.OpenHologram()
        phasePath = self.tempDir + '\\phase.tif'
        self.remoteKoala.remoteCommands.SaveImageToTiffThroughKoala(4, phasePath)
        img = Image.open(phasePath)
        self.curPhase = np.array(img)
    
    # save all phase images at a certain rec dist to disk
    def saveAllPhaseImages(self):
        progressBar = ProgressBar(self.numHolos)
        progressBar.setWindowTitle('Saving phase')
        progressBar.show()
        
        self.tempPhaseDir = utils.CreateDirectory(self.tempDir, 'phase')
        for f in xrange(self.numHolos):
            phasePath = self.tempPhaseDir + '\\phase_' +  str(f).zfill(5) + '.tif'
            self.remoteKoala.first_holo_seq = f
            self.remoteKoala.OpenHologram()
            self.remoteKoala.remoteCommands.SaveImageToTiffThroughKoala(4, phasePath)
            progressBar.progressbar.setValue(f + 1)
            QtGui.qApp.processEvents()
        self.remoteKoala.first_holo_seq = 0
        
        progressBar._active = False
        progressBar.close()
    
    # save all phase images at a series of rec dist to disk
    def saveAllPhaseStackImages(self):
        progressBar = ProgressBar(self.numHolos)
        progressBar.setWindowTitle('Saving phase stacks')
        progressBar.show()
        
        self.tempPhaseStackDir = utils.CreateDirectory(self.tempDir, 'phaseStacks')
        for f in xrange(self.numHolos):
            self.remoteKoala.first_holo_seq = f
            self.remoteKoala.OpenHologram()
            stackPath = utils.CreateDirectory(self.tempPhaseStackDir, 
                                              '\\phaseStack_' +  str(f).zfill(5))
            for stackNum in xrange(self.stackHeightUM + 1):
                d = (stackNum - self.stackHeightUM / 2) / self.zConvert
                self.remoteKoala.remoteCommands.SetRecDistCM(d)
                imgPath = stackPath + '\\phase_' +  str(stackNum).zfill(3) + '.tif'
                self.remoteKoala.remoteCommands.SaveImageToTiffThroughKoala(4, imgPath)
                
            progressBar.progressbar.setValue(f + 1)
            QtGui.qApp.processEvents()
        self.remoteKoala.first_holo_seq = 0
        
        progressBar._active = False
        progressBar.close()
    
    # compute the background as the median of all the phase images
    def computeBcg(self):
        progressBar = ProgressBar(self.numHolos)
        progressBar.setWindowTitle('Computing bcg')
        progressBar.show()
        phaseAll = []
        step = 5
        for f in xrange(0, self.numHolos, step):
            phasePath = self.tempPhaseDir + '\\phase_' +  str(f).zfill(5) + '.tif'
            img = Image.open(phasePath)
            phaseAll.append(np.array(img))
            progressBar.progressbar.setValue(f + step)
            QtGui.qApp.processEvents()
        
        phaseAll = np.asanyarray(phaseAll)
        self.meanPhaseIntensity = np.sum(np.mean(phaseAll, 0))
        self.phaseBcg = np.median(phaseAll, 0).astype('uint8')
        path = self.tempDir + '\\phaseBcg.tif'
        img = Image.fromarray(self.phaseBcg)
        img.save(path)
        
        progressBar._active = False
        progressBar.close()    
    
    def computeStackBcg(self):
        progressBar = ProgressBar(self.stackHeightUM + 1)
        progressBar.setWindowTitle('Computing bcg')
        progressBar.show()
        step = 5
        stackBcgDir = utils.CreateDirectory(self.tempDir, 'stackBcg')
        meanIntensityAll = []
        
        if self.tempPhaseStackDir is None:
            self.tempPhaseStackDir = self.tempDir + '\\phaseStacks'
        for stackNum in xrange(self.stackHeightUM + 1):
            phaseAll = []
            for f in xrange(0, self.numHolos, step):
                imgPath = self.tempPhaseStackDir + '\\phaseStack_' + \
                    str(f).zfill(5) + '\\phase_' + str(stackNum).zfill(3) + '.tif'
                img = Image.open(imgPath)
                phaseAll.append(np.array(img))
        
            phaseAll = np.asanyarray(phaseAll)
            meanIntensityAll.append(np.sum(np.mean(phaseAll, 0)))
            phaseBcg = np.median(phaseAll, 0).astype('uint8')
            path = stackBcgDir + '\\bcg_' + str(stackNum).zfill(3) + '.tif'
            img = Image.fromarray(phaseBcg)
            img.save(path)       
            
            progressBar.progressbar.setValue(stackNum + 1)
            QtGui.qApp.processEvents()
        
        np.savetxt(self.tempDir + '\\meanIntensity.txt', np.asanyarray(meanIntensityAll))
        progressBar._active = False
        progressBar.close()  
    
    def subtractBcg(self):        
        progressBar = ProgressBar(self.numHolos)
        progressBar.setWindowTitle('Subtracting bcg')
        progressBar.show()
        phase_removedBcg = utils.CreateDirectory(self.tempDir, 'phase_removedBcg')
        self.tempPhaseBcgRemovedDir = phase_removedBcg
        for f in xrange(self.numHolos):
            phasePath = self.tempPhaseDir + '\\phase_' +  str(f).zfill(5) + '.tif'
            img = Image.open(phasePath)
            imgArr = np.array(img).astype('int')
            
            factor = self.meanPhaseIntensity / imgArr.sum()
            imgArr = imgArr * factor
            
            imgArr = imgArr - self.phaseBcg
            imgArr[imgArr <=  self.threshold] = 0
            imgArr[imgArr > self.threshold] = 255
            img = Image.fromarray(imgArr.astype('uint8'))
            phasePath = phase_removedBcg + '\\phase_' +  str(f).zfill(5) + '.tif'
            img.save(phasePath)
            
            progressBar.progressbar.setValue(f + 1)
            QtGui.qApp.processEvents()
        
        progressBar._active = False
        progressBar.close()
    
    # to be corrected
    def subtractStackBcg(self):
        progressBar = ProgressBar(self.numHolos)
        progressBar.setWindowTitle('Subtracting bcg stack')
        progressBar.show()
        phase_removedBcg = utils.CreateDirectory(self.tempDir, 'phaseStack_removedBcg')
        
        meanIntensityAll = np.loadtxt(self.tempDir + '\\meanIntensity.txt')
        for f in xrange(self.numHolos):
            stackDir = self.tempDir + '\\phaseStacks\\phaseStack_' +  str(f).zfill(5)
            stackDirOut = utils.CreateDirectory(phase_removedBcg, 'phaseStack_' + str(f).zfill(5))
            for stackNum in xrange(self.stackHeightUM + 1):
                phasePath = stackDir + '\\phase_' + str(stackNum).zfill(3) + '.tif'
                img = Image.open(phasePath)
                bcg = Image.open(self.tempDir + '\\stackBcg' + '\\bcg_' +
                                    str(stackNum).zfill(3) + '.tif')
                imgArr = np.array(img).astype('int')
                
                factor = meanIntensityAll[stackNum] / imgArr.sum()
                imgArr = imgArr * factor
                
                imgArr = imgArr - np.array(bcg)
                imgArr[imgArr <=  self.threshold] = 0
                imgArr[imgArr > self.threshold] = 255
                img = Image.fromarray(imgArr.astype('uint8'))
                phasePathOut = stackDirOut + '\\phase_' +  str(stackNum).zfill(3) + '.tif'
                img.save(phasePathOut)
            
            progressBar.progressbar.setValue(f + 1)
            QtGui.qApp.processEvents()
        
        progressBar._active = False
        progressBar.close()
        
    # not used in direct 3d tracking
    def imgPreprocessing(self):
        progressBar = ProgressBar(self.numHolos)
        progressBar.setWindowTitle('Subtracting bcg')
        progressBar.show()        
        
        srcDir = self.tempPhaseBcgRemovedDir
        if srcDir is None:
            srcDir = self.tempDir + '\\phase_removedBcg'
        outDir = utils.CreateDirectory(self.tempDir, 'phase_postProcessed')
        prevImage = np.array(Image.open(srcDir + '\\phase_' + str(0).zfill(5) + '.tif'))
        for f in xrange(1, self.numHolos):
            img = Image.open(srcDir + '\\phase_' + str(f).zfill(5) + '.tif')
            imgArr = np.array(img)
            # subtract adjacent image
            if self.adjacentSubtraction:
                tempArr = np.copy(imgArr)
                imgArr[prevImage == 0] = 0
                prevImage = tempArr
            
            # erosion
            kernel = np.ones((5, 5), np.uint8)
            imgArr = cv2.erode(imgArr, kernel, iterations = 1)
            
            # dialate and then close
            imgArr = cv2.dilate(imgArr, kernel, iterations = 2)
            imgArr = cv2.morphologyEx(imgArr, cv2.MORPH_CLOSE, kernel)
            
            # save
            img = Image.fromarray(imgArr.astype('uint8'))
            phasePath = outDir + '\\phase_' +  str(f).zfill(5) + '.tif'
            img.save(phasePath)
        
            progressBar.progressbar.setValue(f + 1)
            QtGui.qApp.processEvents()
        
        progressBar._active = False
        progressBar.close()  
    # deprecated
#    def updatePhaseImage(self):
#        self.prevPhase = self.curPhase
#        self.curPhase = self.remoteKoala.remoteCommands.getPhase32fImage(
#            self.contrast[0], self.width, self.height) 
       
    def detectAndTrack3d(self):
        utils.Log("Start detection and tracking")
        openTracks = []
        closedTracks = []
        curPoses = []
        
       
    def updatePhaseImage(self):
        self.prevPhase = self.curPhase
        self.getPhaseImageByHoloNum(self.holoNum)
    
    def resetPhaseImage(self):
        self.remoteKoala.first_holo_seq = 0
        self.getPhaseImageByHoloNum(0)
        self.prevPhase = None
        
    def saveCurPhase(self):
        path = self.directory + "\\test_phase.tif"
        self.remoteKoala.remoteCommands.SaveImageToTiffThroughKoala(4, path) 
        
    def performTracking(self):
        utils.Log("Start tracking ...")
#        self.tempPhaseDir = self.tempDir + '\\phase'
#        if not os.path.exists(self.tempPhaseDir):
#            self.saveAllPhaseImages()
##        self.saveAllPhaseImages()
#        if self.bcgRemove:
#            self.computeBcg()
#            self.subtractBcg()
#        self.xyFromPhase()
#        self.testChangeRecDist()
#        self.saveAllPhaseStackImages()
#        self.computeStackBcg()
#        self.subtractStackBcg()
    # for testing change step rec dist
    
    def testChangeRecDist(self):
        for d in xrange(-15, 15, 1):
            self.remoteKoala.remoteCommands.SetRecDistCM(d)
            time.sleep(0.2)

    # for testing the progress bar
    def testProgressBar(self):
        progressBar = ProgressBar(self.numHolos)
        progressBar.show()
        time.sleep(1)
        progressBar.progressbar.setValue(self.numHolos / 2)
        time.sleep(1)
        progressBar._active = False
        progressBar.close()        
    
class HVSeg:
    def __init__(self, top, left, length, orient):
        self.top = top
        self.left = left
        self.length = length
        self.orient = orient
    
class Track:
    def __init__(self):
        self.startHoloNum = -1
        self.endHoloNum = -1
        self.pos = []

class Position:
    def __init__(self, pos):
        self.isOpen = True
        self.pos = None
        self.curPos = pos
        self.zStack = []
        self.holoNum = -1
    def closePosition(self):
        self.isOpen = False
        
        
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
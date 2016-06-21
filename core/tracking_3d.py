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

class Tracking3D:
    tracking = None
    def __init__(self):
        self.remoteKoala = None
        self.db = None
        # remote parameters
        self.contrast = ['lambda1', 'mapping']
        self.centRecDist = -2.5 # initial reconstruction distance
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
        
        # sample related parameters
        self.chamberIndex = None
        self.mediumIndex = None
        
        # tracking parameters and data
        self.phaseBcg = None
        self.meanPhaseIntensity = None
        self.numHolos = 0
        self.holoNum = 0
        self.prevPhase = None
        self.curPhase = None
        self.closedtracks = []
        self.activeTracks = []
        self.reconsDist = 0
        self.threshMethod = None
        self.zfocusMethod = None
        self.directory = None
        self.bcgRemove = True
        self.tempDir = None
        self.tempPhaseDir = None
        
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
        
        # some testing tasks
#        self.getCurPhaseImage()
#        self.iterOverAllFrames()
#        self.computeBcg()
#        self.testProgressBar()
#        self.saveAllPhaseImages()
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
        
        self.stepUM = self.stepRecDist * self.zConvert
        self.stackHeightUM = self.rangeRecDist * 2 * self.zConvert
        
        
    # basic ajustments for reconstruction
    def reconsAdjust(self):
        self.remoteKoala.remoteCommands.SetRecDistCM(self.centRecDist)
        self.remoteKoala.remoteCommands.SetUnwrap2DState(True)
        
        # define some line segments for tilt correction
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
    
    # save all phase images at current rec dist to disk
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
    
    def subtractBcg(self):        
        progressBar = ProgressBar(self.numHolos)
        progressBar.setWindowTitle('Subtracting bcg')
        progressBar.show()
        
        phase_removedBcg = utils.CreateDirectory(self.tempDir, 'phase_removedBcg')
        for f in xrange(self.numHolos):
            phasePath = self.tempPhaseDir + '\\phase_' +  str(f).zfill(5) + '.tif'
            img = Image.open(phasePath)
            img_arr = np.array(img).astype('int')
            
            factor = self.meanPhaseIntensity / img_arr.sum()
            img_arr = img_arr * factor
            
            img_arr = img_arr - self.phaseBcg
            img_arr[img_arr < 0] = 0
            img_arr[img_arr > 10] = 255

#            img_arr = cv2.equalizeHist(img_arr)
            img = Image.fromarray(img_arr.astype('uint8'))        
            phasePath = phase_removedBcg + '\\phase_' +  str(f).zfill(5) + '.tif'
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
        self.tempPhaseDir = self.tempDir + '\\phase'
        if not os.path.exists(self.tempPhaseDir):
            self.saveAllPhaseImages()
#        self.saveAllPhaseImages()
        if self.bcgRemove:
            self.computeBcg()
            self.subtractBcg()
        


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
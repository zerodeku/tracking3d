# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 13:15:23 2016

@author: zerodeku
"""

from tracking3d.core import remote_Koala
from tracking3d.core import SqliteDB
import time
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
        self.recDist = -2.5 # initial reconstruction distance
        self.width = 750 # image width
        self.height = 750 # image height
        self.pixSize = 1 # real value is updated when holograph is loaded
        self.KoalaViewSize = 748 # Koala software phase image view size
        # tracking parameters and data
        self.phaseBcg = None
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
        self.removeBcg = False
        self.progressBar = None
        self.tempDir = None
        
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
        self.pixSize = self.remoteKoala.pxSize
        self.reconsAdjust()
        self.updatePhaseImage()
        
        # some testing tasks
#        self.getCurPhaseImage()
        self.iterOverAllFrames()
    
    # basic ajustments for reconstruction
    def reconsAdjust(self):
        self.remoteKoala.remoteCommands.SetRecDistCM(self.recDist)
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
            
    # go through all holos one by one
    def iterOverAllFrames(self):
        for f in xrange(self.remoteKoala.framesNumber):
            tic = time.clock()
            self.getPhaseImageByHoloNum(f)
            toc = time.clock()
            print "Elapsed time: " + str(toc - tic)
        self.remoteKoala.first_holo_seq = 0
        
    def saveCurPhase(self):
        path = self.directory + "\\test_phase.tif"
        self.remoteKoala.remoteCommands.SaveImageToTiffThroughKoala(4, path) 
        
    # compute the background as the median of all the phase images
    def computeBcg(self):
        phaseAll = []
        self.progressBar = ProgressBar(self.numHolos)
        self.progressBar.show()
        for f in xrange(self.remoteKoala.framesNumber):
            self.progressBar.progressbar.setValue(self.numHolos - 1 - f)
            self.getPhaseImageByHoloNum(f)
            phaseAll.append(self.curPhase)
        self.progressBar._active = False
        self.resetPhaseImage()
        self.phaseBcg = np.median(phaseAll, 0)
        utils.DisplayImage(self.phaseBcg)
        path = utils.CreateDirectory(self.directory, 'phaseBcg.tif')
        bcg = np.array(self.phaseBcg)
        img = Image.fromarray(bcg)
        img.save(path)

    def performTracking(self):
        if self.removeBcg:
            self.computeBcg()

    
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
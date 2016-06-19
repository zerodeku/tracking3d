# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 13:15:23 2016

@author: zerodeku
"""

from tracking3d.core import remote_Koala
from tracking3d.core import SqliteDB
import time

class Tracking3D:
    tracking = None
    def __init__(self):
        self.remoteKoala = None
        self.db = None
        # remote parameters        
        self.recDist = -2.5 # initial reconstruction distance
        self.width = 1024 # image width
        self.height = 1024 # image height
        self.pixSize = 1 # real value is updated when holograph is loaded
        self.KoalaViewSize = 748 # Koala software phase image view size
        # tracking parameters and data
        self.numHolos = 0
        self.holoNumber = 0
        self.prevPhase = None
        self.curPhase = None
        self.closedtracks = []
        self.activeTracks = []
        self.reconsDist = 0
        self.threshMethod = None
        self.zfocusMethod = None
        
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
        self.remoteKoala = remote_Koala.RemoteKoala.from_dataset(self.db, data)
        
    # connect Koala and open windows
    def setupKoala(self):
        self.remoteKoala.ConnectKoala()
        self.remoteKoala.OpenConfig()
        self.remoteKoala.OpenHologram()
        time.sleep(0.1)
        self.pixSize = self.remoteKoala.remoteCommands.GetPxSizeUm()
        self.reconsAdjust()
        
    
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
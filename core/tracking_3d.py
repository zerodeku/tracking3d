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
from enum import Enum
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.io
import pickle

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
        self.threshold = 50
        self.directory = None
        self.bcgRemove = True
        self.tempDir = None
        self.tempPhaseDir = None
        self.tempPhaseBcgRemovedDir = None 
        self.tempPhaseStackDir = None
        self.zSearchDX = 5
        self.zSearchDY = 5
        self.nextSearchDX = 0
        self.nextSearchDY = 0
        self.allFramePoints = None
        self.allFramePoses = None
        self.allPosObj = None
        self.curPhase = None
        self.prevPhase = None
        self.holoNum = 0
        self.fps = -1
        self.allTracks = None
        self.trackDX = 3
        self.trackDY = 3
        self.trackDZ = 5
        
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
        self.getFPS()
        
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
        print self.pxSize
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
        step = 1
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
       
    def detect3DPoses(self):
        utils.Log("--Start detection and tracking")
        progressBar = ProgressBar(self.numHolos)
        progressBar.setWindowTitle('Detecting 3d')
        progressBar.show()           
        
        detector = self.getDetector()
        allFramePoses = []
        allFramePoints = []
        allPosObj = []
        numFrames = self.numHolos
        for f in xrange(numFrames):
            stackDir = self.tempDir + '\\phaseStack_removedBcg\\phaseStack_' +\
                str(f).zfill(5)
            framePoses = [[], [], []]
            framePoints = [[], [], []]
            openPoses = []
            framePosObj = []
            for stackNum in xrange(self.stackHeightUM + 1):
                imgPath = stackDir + '\\phase_' + str(stackNum).zfill(3) + '.tif'
                # store nonzero points
                img = Image.open(imgPath)
                imgArr = np.array(img)
                nind = np.nonzero(imgArr)
                xs = nind[1] * self.pxSize
                ys = nind[0] * self.pxSize
                zs = [stackNum] * len(xs)                
                framePoints[0] += xs.tolist()
                framePoints[1] += ys.tolist()
                framePoints[2] += zs
                # extract positions
                img = cv2.imread(imgPath, cv2.IMREAD_GRAYSCALE)
                blobs = detector.detect(img)
                stackPoses = []
                for blob in blobs:
#                    print str(stackNum) + ': '
#                    print blob.pt
                    stackPoses.append(Position(f, (blob.pt[0] * self.pxSize, 
                                        blob.pt[1] * self.pxSize, stackNum)))
                    
                nextOpenPoses = []
                for pos0 in openPoses:
                    pos0.isOpen = False
                    for pos in stackPoses:
                        if not pos.isMerged and self.isPosInROI(pos0, pos):
                            pos0.isOpen = True
                            pos0.xyPoses.append(pos.xyPoses[-1])
                            pos0.zPoses.append(stackNum)
                            nextOpenPoses.append(pos0)
                            pos.isMerged = True
                            break
                    
                    if not pos0.isOpen or stackNum == self.stackHeightUM:
                        pos0.closePosition()
                        if len(pos0.zPoses) >= 4:
                            framePosObj.append(pos0)
                            framePoses[0].append(pos0.pos[0])
                            framePoses[1].append(pos0.pos[1])
                            framePoses[2].append(pos0.pos[2])
                        
                for pos in stackPoses:
                    if not pos.isMerged:
                        nextOpenPoses.append(pos)
                openPoses = nextOpenPoses
            
            allPosObj.append(framePosObj)
            allFramePoints.append(framePoints)
            allFramePoses.append(framePoses)
            
            progressBar.progressbar.setValue(f + 1)
            QtGui.qApp.processEvents()
        
        self.allPosObj = allPosObj
        self.allFramePoints = allFramePoints
        self.allFramePoses = allFramePoses
        
        progressBar._active = False
        progressBar.close()         
        
    def track3d(self):
        utils.Log('--Perform tracking in 3d')
        
        progressBar = ProgressBar(self.numHolos)
        progressBar.setWindowTitle('Tracking in 3d')
        progressBar.show()           
        if self.allPosObj is None:
            self.constructPosObj()
            
        allTracks = []
        numFrames = self.numHolos
        openTracks = []
        for f in xrange(numFrames):
#            print "Frame: " + str(f)
            framePoses = self.allPosObj[f]
#            print "No frame poses: " + str(len(framePoses))
            nextOpenTracks = []
#            print "tracks : poses " + str(len(openTracks)) + ":" + str(len(framePoses))
            for track in openTracks:
                track.isOpen = False
                for pos in framePoses:
                    if not pos.isTracked and self.isPosInTrackROI(track, pos):
                        track.isOpen = True
                        track.poses.append(pos.pos)
                        nextOpenTracks.append(track)
                        pos.isTracked = True
#                        print "Tracked pos"
                        break
                
                # Check the next frame for possible connection
                if not track.isOpen and f + 1 < numFrames:
                    nextFramePoses = self.allPosObj[f + 1]
                    for pos in nextFramePoses:
                        if self.isPosInTrackROI(track, pos):
                            track.isOpen = True
                            track.poses.append(track.poses[-1])
                            nextOpenTracks.append(track)
                            break
            
                if not track.isOpen or f == numFrames - 1:
#                    print "Closed a track with length: " + str(len(track.poses))
                    track.endHolo = f
                    track.isOpen = False
                    allTracks.append(track)
                    
            for pos in framePoses:
                if not pos.isTracked:
#                    print "New track found!"
                    newTrack = Track(f, pos.pos)
                    nextOpenTracks.append(newTrack)

            openTracks = nextOpenTracks
            progressBar.progressbar.setValue(f + 1)
            QtGui.qApp.processEvents()
            
        self.allTracks = allTracks
        self.exportTracks()
#        for f in xrange(numFrames):
#            print "Frame " + str(f) + ':' + str(len(allTracks[f]))
        progressBar._active = False
        progressBar.close()
    
    def constructPosObj(self):
        path = self.tempDir + '\\results\\allFramePoses.mat'
        tempDict = scipy.io.loadmat(path)
        allFramePoses = tempDict.values()[1]
        allPosObj = []

        for f in xrange(np.size(allFramePoses, 0)):
            framePosObj = []
            framePoses = allFramePoses[f]
            for i in xrange(np.size(framePoses[0])):
                posObj = Position(f, (framePoses[0][0][i], framePoses[1][0][i], 
                                      framePoses[2][0][i]))
                posObj.closePosition()
                framePosObj.append(posObj)
            allPosObj.append(framePosObj)
        self.allPosObj = allPosObj
    def exportPositions(self):
        utils.Log('--Exporting position results')
        self.resultDir = utils.CreateDirectory(self.tempDir, 'results')
        holoNum = xrange(self.numHolos)
        res = dict(holoNum=holoNum, framePoints=np.asarray(self.allFramePoints))
        scipy.io.savemat(self.resultDir + '\\allFramePoints.mat', res)
        res = dict(holoNum=holoNum, framePoses=np.asarray(self.allFramePoses))
        scipy.io.savemat(self.resultDir + '\\allFramePoses.mat', res)
        # drawFrames
        self.drawFrames(self.numHolos)        
        
    def exportTracks(self):
        utils.Log('--Exporting tracking results')
        self.resultDir = utils.CreateDirectory(self.tempDir, 'results')
        holoNum = xrange(self.numHolos)
        res = dict(holoNum=holoNum, tracks=np.asarray(self.allTracks))
        scipy.io.savemat(self.resultDir + '\\allTracks.mat', res)
        
    def isPosInROI(self, pos0, pos):
        (x0, y0) = pos0.xyPoses[-1]
        (x, y) = pos.xyPoses[-1]
        if abs(x - x0) <= self.zSearchDX and abs(y - y0) <= self.zSearchDY:
            return True
        return False
    
    def isPosInTrackROI(self, track, pos):
        (x0, y0, z0) = track.poses[-1]
        (x, y, z) = pos.pos
        if abs(x - x0) <= self.trackDX and abs(y - y0) <= self.trackDY and \
            abs(z - z0) <= self.trackDZ:
            return True
        return False
    
    def getDetector(self):
        # Setup SimpleBlobDetector parameters.
        params = cv2.SimpleBlobDetector_Params()
        # Change thresholds and distance
        params.minThreshold = 20;
        #params.maxThreshold = 256;
        params.minDistBetweenBlobs = 75     
        # Filter by Color.
        params.filterByColor = True
        params.blobColor = 255       
        # Filter by Area.
        params.filterByArea = True
        params.minArea = 20
        # Filter by Circularity
        params.filterByCircularity = False
        params.minCircularity = 0.1       
        # Filter by Convexity
        params.filterByConvexity = False
        params.minConvexity = 0.87       
        # Filter by Inertia
        params.filterByInertia = False
        params.minInertiaRatio = 0.01       
        # Create a detector with the parameters
        detector = cv2.SimpleBlobDetector(params)
        return detector
        
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

    def drawFrames(self, numFrames):
        fig = plt.figure(2)
        frameDir = utils.CreateDirectory(self.resultDir, 'frames3d')
        for f in xrange(numFrames):
            ax = Axes3D(fig)
            ax.set_aspect('equal')
            
            framePoints = self.allFramePoints[f]
            X = np.asarray(framePoints[0])
            Y = np.asarray(framePoints[1])
            Z = np.asarray(framePoints[2])
            points = ax.scatter(X, Y, Z, s=0.1, marker='.', color='g', depthshade=False)
            
            framePoses = self.allFramePoses[f]
    #        print framePoses
            X = np.asarray(framePoses[0])
            Y = np.asarray(framePoses[1])
            Z = np.asarray(framePoses[2])
            poses = ax.scatter(X, Y, Z, s=80, color='k', depthshade=False)
            
    #        points.remove()
            ax.set_xlim(0, 70)
            ax.set_ylim(0, 70)
            ax.set_zlim(0, 70)
            ax.set_xlabel('x (um)')
            ax.set_ylabel('y (um)')
            ax.set_zlabel('z (um)')
            plt.grid()
            plt.show()
            plt.savefig(frameDir + '\\frame_' + str(f).zfill(5) + '.png')
            plt.cla()
        
        plt.close()
        
    def getFPS(self):
        timePath = self.db.DefaultDirectory + '\\timestamps.txt'
        tab1 = np.genfromtxt(timePath, dtype=str)
        times = tab1[:, 3]
        times = times.astype(np.float)
        steps = np.diff(times, 1)
        self.fps = 1 / (np.mean(steps) * 1e-3)
    
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
        
        
    def performTracking(self):
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
        self.numHolos = 120
#        self.computeStackBcg()
#        self.subtractStackBcg()
        
        self.detect3DPoses()
        self.exportPositions()
        self.track3d()

        
# line segment class        
class HVSeg:
    def __init__(self, top, left, length, orient):
        self.top = top
        self.left = left
        self.length = length
        self.orient = orient

# Trajectory class    
class Track:
    def __init__(self, startHolo, pos):
        self.startHolo = startHolo
        self.endHolo = -1
        self.poses = [pos]
        self.isOpen = True

# Ending type of tracjectory    
class EndType(Enum):
    OutOfField = 1
    LoseTrack = 2

class Position:
    def __init__(self, holoNum, (x, y, z)):
        self.isOpen = True
        self.isMerged = False
        self.isTracked = False
        self.xyPoses = [(x, y)]
        self.zPoses = [z]
        self.holoNum = holoNum
        self.pos = (-1, -1, -1)
    def closePosition(self):
        h = len(self.zPoses)
        if h % 2 != 0:
            self.pos = self.xyPoses[h / 2] + (self.zPoses[h / 2], )
        else:
            (x0, y0) = self.xyPoses[h / 2 - 1]
            (x1, y1) = self.xyPoses[h / 2]
            (x, y) = ((x0 + x1) / 2.0, (y0 + y1) / 2.0)
            self.pos = (x, y, np.mean(self.zPoses))
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
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 13:56:55 2014

@author: tcolomb
"""


import numpy as np
import os
import sys
import time

from tracking3d.gui import gui_construction
from tracking3d.core import Tracking_process
from tracking3d.core import utilsForProcessing as utils
from tracking3d.core import UnwrapK
from guidata.dataset.qtwidgets import DataSetEditGroupBox
from guidata.qt.QtGui import QCheckBox, QVBoxLayout, QGroupBox
from guidata.qt.QtGui import QMenuBar, QAction, QWidget, QMainWindow, QSplitter
from guidata.qt.QtCore import SIGNAL
from guidata.qt import QtCore
import matplotlib.pyplot as plt
from guiqwt.builder import make
from threading import Thread
from tracking3d.core import SqliteDB


class KMainWindow(Thread, QMainWindow):
    def __init__(self):
        #initialize the main panel
        QMainWindow.__init__(self)
        Thread.__init__(self)
        widget = QWidget()
        self.setCentralWidget(widget)
        #Define menu bar and associated event
        self.menu_bar = QMenuBar(self)
        menu = self.menu_bar.addMenu('File')
        exit_action = QAction('Exit', self)
        exit_action.triggered.connect(self.close)
        defaultdir_action = QAction('Default directory', self)
        defaultdir_action.triggered.connect(self.defineDefaultDirectory)
        menu.addAction(exit_action)
        menu.addAction(defaultdir_action)
        self.setMenuBar(self.menu_bar)
        #Define main object
        self.remoteKoala = None#class remoteKoala to send commands to Koala
        self.aberrationCorrection = None#class to compensate for aberration 
                                        #(uniform filtering)
        self.stackDistance = None #class to acquire stack from holo and perform
                                #xy particle detection
        self.tracking = None #class to perform the tracking (follows particule 
                            #from an hologram to another)
        self.win_aberration = None #window to display Filtering correction
        self.win_projection = None #image window of projected phas with mask of
                                    #particule obtained from threshold
        self.win_derivate = None #image window of the derivation 
                                #(projection1-projection0)
        self.win_polynomial = None #image window of the phase where tilt 
                                #correction, selection of particule can be done
        self.db = None #class of the database methods (read and write from 
                        #objectiveDB.db3)
        self.focusPoint = None #focus Point array (x,y,rec dist)
        self.trackingMode = 1#mode of tracking (0=Automatic, 1=single particule
                            #2=loop on single particule)
        self.UnwrapK = None

        self.setWindowTitle('3D Tracking')
        #init dataset from gui_construction
        self.gb_dhm = DataSetEditGroupBox("DHM", 
                                          gui_construction.DHMParameters, 
                                          comment='')
        self.gb_sample = DataSetEditGroupBox("Sample", 
                                             gui_construction.SampleParameters, 
                                             comment='')
        self.gb_holder = DataSetEditGroupBox("Focus", 
                                             gui_construction.HolderParameters, 
                                             comment='')
        self.gb_remote = DataSetEditGroupBox("Sequence", 
                                             gui_construction.RemoteParameters, 
                                             comment='')
        self.gb_rec = DataSetEditGroupBox("Preprocessing Parameters", 
                                          gui_construction.ReconstructionParameters, 
                                          comment='')
        self.gb_proc = DataSetEditGroupBox("Process", 
                                           gui_construction.ProcessParameters, 
                                           comment='')
        self.gb_partchoice = DataSetEditGroupBox("Particle Choice", 
                                                 gui_construction.ParticuleChoice, 
                                                 comment='')  
        
        #create checkbox button to display the different windows
        self.gb_btn = QGroupBox(self)
        self.gb_btn.setTitle('Display Results for')
        self.vbox = QVBoxLayout(self.gb_btn)
       
        self.cb_polynomial = QCheckBox('Polynomial Correction', self)
        self.cb_polynomial.setTristate(False)
        self.cb_polynomial.stateChanged.connect(self.show_polynomial)
        
        self.cb_aberration = QCheckBox('Filtering Correction', self)
        self.cb_aberration.setTristate(False)
        self.cb_aberration.stateChanged.connect(self.show_aberration)

        self.cb_projection = QCheckBox('Projection with particle detection', self)
        self.cb_projection.setTristate(False)
        self.cb_projection.stateChanged.connect(self.show_projection)
 
        self.vbox.addWidget(self.cb_polynomial)
        self.vbox.addWidget(self.cb_aberration)
        self.vbox.addWidget(self.cb_projection)

        self.gb_btn.setLayout(self.vbox)
        self.gb_btn.setEnabled(True)
        
        #intialize data object for the different windows
        self.aberration_data = None #image window self.win_aberration
        self.projection_data = None #image window self.win_projection
        self.derivate_data = None #image window self.win_derivate 
        self.phase_data = None #image window self.win_polynomial

        #associate events to dataset apply buttons
        self.connect(self.gb_dhm, SIGNAL("apply_button_clicked()"), 
                     self.update_dhm_params)
        self.connect(self.gb_sample, SIGNAL("apply_button_clicked()"), 
                     self.update_sample_params)
        self.connect(self.gb_holder, SIGNAL("apply_button_clicked()"), 
                     self.update_holder_params)
        self.connect(self.gb_rec, SIGNAL("apply_button_clicked()"), 
                     self.update_rec_params)
        self.connect(self.gb_remote, SIGNAL("apply_button_clicked()"), 
                     self.update_remote_params)
        
        #organize dataset in the main panel
        splitter1 = QSplitter(QtCore.Qt.Vertical)
        splitter1.addWidget(self.gb_dhm)
        splitter1.addWidget(self.gb_remote)
        
        splitter2 = QSplitter(QtCore.Qt.Vertical)
        splitter2.addWidget(self.gb_sample)
        splitter2.addWidget(self.gb_holder)
        splitter = QSplitter(self)

        splitter.addWidget(splitter1)
        splitter.addWidget(splitter2)

        splitter.addWidget(self.gb_rec)

        splitter3 = QSplitter(QtCore.Qt.Vertical)
        splitter3.addWidget(self.gb_proc)
        splitter3.addWidget(self.gb_partchoice)
        splitter.addWidget(splitter3)

        splitter.addWidget(self.gb_btn)
        self.setCentralWidget(splitter)
        self.setContentsMargins(10, 5, 10, 5)   
        
        #Get all parameters from the dataset
        self.gb_dhm.get()
        self.gb_sample.get()
        self.gb_holder.get()
        self.gb_rec.get()
        self.gb_remote.get()
        self.gb_partchoice.get()
#    
    def closeEvent(self, event):
        #event when main panel is closed
        if self.remoteKoala.connected:
            self.remoteKoala.CloseConnection()
        self.show_aberration(False)
        self.show_projection(False)
        self.show_polynomial(False)

    def defineDefaultDirectory(self, event):
        #write the new default directory in the database from menu selection
        db = SqliteDB.DBReader()
        db.DefineDefaultDirectory()
        
    def saveParamsToFile_change(self):
        #save dataset parameters to text file to be able to reload them
        self.update_rec_params()
        Data = []
        #Data.append(('pixel_size_um',(float)(self.gb_dhm.dataset.pixel_size_um),'float'))
        Data.append(('MODB', (float)(self.gb_dhm.dataset.MODB), 'float'))
        #Data.append(('NA',(float)(self.gb_dhm.dataset.NA),'float'))
        #Data.append(('transmission',(int)(self.gb_dhm.dataset.transmission),'bool'))
        Data.append(('index_medium', (float)(self.gb_sample.dataset.index_medium), 'float'))
        Data.append(('index_chamber', (float)(self.gb_sample.dataset.index_chamber), 'float'))        
        Data.append(('index_sample', (float)(self.gb_sample.dataset.index_sample), 'float'))
        Data.append(('max_displacement', (float)(self.gb_sample.dataset.max_displacement), 'float'))
        Data.append(('roi_sample', (int)(self.gb_holder.dataset.roi_sample), 'int'))
        Data.append(('central_rec_dist', (float)(self.gb_holder.dataset.central_rec_dist), 'float'))
        Data.append(('range_rec_dist', (float)(self.gb_holder.dataset.range_rec_dist), 'float'))
        Data.append(('ratio_DOF', (float)(self.gb_holder.dataset.ratio_DOF), 'float'))
        Data.append(('focusMethod', (int)(self.gb_holder.dataset.focusMethod), 'int'))
        Data.append(('first_holo_seq', (int)(self.gb_remote.dataset.first_holo_seq), 'int'))
        Data.append(('last_holo_seq', (int)(self.gb_remote.dataset.last_holo_seq), 'int'))
        Data.append(('use_filter_size_aberration', (int)(self.gb_rec.dataset.use_filter_size_aberration), 'bool'))
        Data.append(('filter_size_aberration', (int)(self.gb_rec.dataset.filter_size_aberration), 'int'))
        Data.append(('use_filter_amplitude', (int)(self.gb_rec.dataset.use_filter_amplitude), 'bool'))
        Data.append(('XYImport', (int)(self.gb_rec.dataset.XYImport), 'bool'))
        Data.append(('XYImportFilePath', self.gb_rec.dataset.XYImportFilePath, 'string'))
        Data.append(('bool_deriv', (int)(self.gb_rec.dataset.bool_deriv), 'bool'))
        Data.append(('particlecontrast', (int)(self.gb_rec.dataset.particlecontrast), 'int'))
        Data.append(('threshold', (float)(self.gb_rec.dataset.threshold), 'float'))
        Data.append(('bool_ThresholdMethod', (int)(self.gb_rec.dataset.bool_ThresholdMethod), 'int'))
        Data.append(('it_eros', (int)(self.gb_rec.dataset.it_eros), 'int'))
        Data.append(('it_dil', (int)(self.gb_rec.dataset.it_dil), 'int'))
        Data.append(('saveAllStack', (int)(self.gb_remote.dataset.saveAllStack), 'bool'))
        utils.SaveParamsFileByLine(Data, self.gb_remote.dataset.directoryPath+'\\Tracking_parameters.txt')
    
    def loadParamsFromFile_change(self):
        #Load saved dataset parameters from text file
        fname = utils.OpenTxtFile('Open Parameters File', 
                                               self.gb_remote.dataset.directoryPath)
        Data = utils.ReadParamsFromLine(fname)
        for k in range(2):
            for info in Data:
                if info[2] == 'float':
                    value = (float)(info[1])
                elif info[2] == 'int':
                    value = (int)(info[1])
                elif info[2] == 'bool':
                    value = (int)(info[1])
                elif info[2] == 'string':
                    value = info[1]
                else:
                    utils.ErrorMessage('info data in the file is unknown')

                setattr(self.gb_dhm.dataset, info[0], value)
                setattr(self.gb_sample.dataset, info[0], value)
                setattr(self.gb_holder.dataset, info[0], value)      
                setattr(self.gb_remote.dataset, info[0], value)
                setattr(self.gb_rec.dataset, info[0], value)
            
            #get the new parameters in the main panel
            self.gb_dhm.get()
            self.gb_sample.get()
            self.gb_holder.get()
            self.gb_rec.get()
            self.gb_remote.get()
        #as parameters change, change children of the different class
        self.MODB_change() #change class parameters related to MO (config,configType,NA, ...)
        self.update_remote_params()
        self.update_sample_params()
        self.update_rec_params()
        
    def update_dhm_params(self):
        #update the dhm_params dataset and modify parameters of corresponding class accordingly       
        self.process() #initialize  object or modify them if they are not None        
        self.stackDistance.from_dataset(self.db, self.gb_dhm.dataset, 
                                        self.gb_sample.dataset, 
                                        self.gb_holder.dataset, 
                                        self.gb_rec.dataset, 
                                        self.gb_remote.dataset)
    def update_sample_params(self): 
        #update sample_params
        self.process()
        #modify class parameters related to sample_params
        self.stackDistance.index_medium = self.gb_sample.dataset.index_medium
        self.stackDistance.index_sample = self.gb_sample.dataset.index_sample
        self.stackDistance.index_chamber = self.gb_sample.dataset.index_chamber
        self.tracking.max_displacement = self.gb_sample.dataset.max_displacement
        self.tracking.roi_sample = self.gb_holder.dataset.roi_sample
        self.UpdateSamplePlaneStepUM()

    def update_holder_params(self):
        #updade holder_params (see below the parameters that have to be modified in objects)
        self.process()
        self.stackDistance.central_rec_dist = self.gb_holder.dataset.central_rec_dist
        self.stackDistance.range_rec_dist = self.gb_holder.dataset.range_rec_dist
        self.stackDistance.ratio_DOF = self.gb_holder.dataset.ratio_DOF
        self.stackDistance.roi_sample = self.gb_holder.dataset.roi_sample
        self.UpdateSamplePlaneStepUM()
        #Reset the phaseProjection data of stack Distance
        self.stackDistance.ResetPhaseProjection()

    def update_remote_params(self): 
        #update remote parameters depending of what has changed
        self.process()
        directorychanged = not self.gb_remote.dataset.directoryPath == self.remoteKoala.directoryPath
        if directorychanged:
            self.remoteKoala.directoryPath = self.gb_remote.dataset.directoryPath#change the sequence directory path
            self.remoteKoala.DefineSeqFormat() #define the sequence format for the new directory
            self.gb_remote.dataset.last_holo_seq = (int)(self.remoteKoala.framesNumber)-1 #define the last_holo of the sequence
            self.gb_remote.get()#write the last holo seq in panel
            self.tracking.last_holo_seq = self.gb_remote.dataset.last_holo_seq+1 #write new last holo value in tracking object

            self.stackDistance.directoryPath = self.gb_remote.dataset.directoryPath #change directory in stackDistance object
            self.tracking.stackDistance = self.stackDistance #modify the stackDistance object in class tracking
                                                                           
        #Test if the directoryPath is correct (exists and have Holograms sequence)
        if self.remoteKoala.directoryPath == None:
            utils.ErrorMessage("Please provide a sequence directory")
            return
        else:
            if not os.path.exists(self.remoteKoala.directoryPath+'\\Holograms'):
                utils.ErrorMessage("Please provide a valid sequence directory")
                return
        #connect to koala only if it is not the case
        if not self.remoteKoala.connected:
            self.remoteKoala.ConnectKoala()
        #test if the MO configuration has changed (not necessary to reload it if it is the same)
        configchanged = not self.remoteKoala.lastconfig == self.db.ProductionNo
       
        
        if configchanged or directorychanged: #change configuration also if directory changed
            self.remoteKoala.config = self.db.ProductionNo #read the config from the database
            self.stackDistance.phaseProjection1 = None #reset phase projection 1
            self.remoteKoala.lastfirst_holo_seq = -1 #init first holo seq to -1 to be sure to compute again this value 
            self.remoteKoala.OpenConfig() #remote koala command to open the configuration
            self.remoteKoala.GetLambdaNm()
            self.gb_holder.dataset.central_rec_dist = self.remoteKoala.remoteCommands.GetRecDistCM() #get the default reconstruction distance from remoteCommand class (=object of remoteKoala class)
            self.gb_holder.get() #write parameters in the main panel
        if not self.remoteKoala.lastfirst_holo_seq == self.gb_remote.dataset.first_holo_seq: #if first_holo change (always true if config changed or directory changed as self.remoteKoala.lastfirst_holo_seq=-1)
            self.remoteKoala.amp = None #reset amplitude image 
            self.stackDistance.phaseProjection1 = None # reset phase image
            self.remoteKoala.first_holo_seq = self.gb_remote.dataset.first_holo_seq       
            self.remoteKoala.OpenHologram() #Open hologram from remote 
            self.drawPolynomialImage() #draw the polynomial correction image =self.win_polynomial
        self.tracking.last_holo_seq = self.gb_remote.dataset.last_holo_seq+1
        
        #when user change the first image to analyze
        firstholochanged = not self.gb_remote.dataset.first_holo_seq == self.remoteKoala.first_holo_seq
        if firstholochanged:
            self.remoteKoala.first_holo_seq = self.gb_remote.dataset.first_holo_seq
            self.remoteKoala.amp = None
            self.stackDistance.ResetPhaseProjection()
            self.remoteKoala.first_holo_seq = self.gb_remote.dataset.first_holo_seq
            self.stackDistance.first_holo_seq = self.gb_remote.dataset.first_holo_seq
            self.stackDistance.use_filter_size_aberration = self.gb_rec.dataset.use_filter_size_aberration
            self.remoteKoala.OpenHologram()
            self.drawPolynomialImage()
        self.stackDistance.pxSize = self.remoteKoala.pxSize
        self.UnwrapK.pxSize = self.remoteKoala.pxSize
        #self.UnwrapK.hconv=self.remoteKoala.wavelength/4./np.pi*1e-9
        self.UnwrapK.hconv = 1e-9
        self.UpdateSamplePlaneStepUM()
        self.stackDistance.InitBcgPropagation()
        self.tracking.stackDistance = self.stackDistance
    def update_rec_params(self): 
        self.process()
        aberrationchanged = not self.gb_rec.dataset.use_filter_size_aberration == self.stackDistance.use_filter_size_aberration
        self.remoteKoala.from_dataset(self.db, self.gb_remote.dataset)   
        self.stackDistance.from_dataset(self.db, self.gb_dhm.dataset, self.gb_sample.dataset, self.gb_holder.dataset, 
                                        self.gb_rec.dataset, self.gb_remote.dataset)
        self.stackDistance.use_filter_size_aberration = self.gb_rec.dataset.use_filter_size_aberration
        self.aberrationCorrection.from_dataset(self.gb_rec.dataset)
        if aberrationchanged: #reconstruct the compensated polynomial iamge
            self.remoteKoala.first_holo_seq = self.gb_remote.dataset.first_holo_seq
            self.remoteKoala.amp = None
            self.stackDistance.ResetPhaseProjection()
            self.remoteKoala.first_holo_seq = self.gb_remote.dataset.first_holo_seq
            self.stackDistance.first_holo_seq = self.gb_remote.dataset.first_holo_seq
            self.stackDistance.use_filter_size_aberration = self.gb_rec.dataset.use_filter_size_aberration
            self.remoteKoala.OpenHologram()
            self.drawPolynomialImage()
       
    def adjust_aberration_change(self):
        #event when parametesr of aberraction compensation changed
        self.aberrationCorrection.filter_size_aberration = self.gb_rec.dataset.filter_size_aberration        
        if self.remoteKoala.amp == None or self.remoteKoala.ph == None:
            self.remoteKoala.AcquireWavefront()#acquire amp and phase data if at least one is None
        if self.gb_rec.dataset.use_filter_amplitude:#if use_filter_amplitude->perform complex uniform filter
            amp = self.remoteKoala.amp #get amplitude from remote class
            ph = self.remoteKoala.ph #get phase from remote class
            self.aberrationCorrection.complex = amp*np.exp(complex(0., 1.)*ph) #set complex in  aberractionCorrection object
            self.aberrationCorrection.ComplexUniformFilter() #compute filter
        else:
            self.aberrationCorrection.ph = self.remoteKoala.ph
            self.aberrationCorrection.PhaseUniformFilter()
        #display result in Filtering correction image 
        self.gb_btn.setEnabled(True)
        if self.aberration_data is None:
            self.aberration_data = make.image(data=self.aberrationCorrection.ph_filter, colormap="gray")
        else:
            self.aberration_data.set_data(self.aberrationCorrection.ph_filter)
        if self.win_aberration is not None:
            self.win_aberration.get_plot().replot()
        
    def adjust_detection_change(self):
        
        #event when user change the xy detection (from adjust detection button, but also if focus detection is changed, or treshold method changed)
        if not self.stackDistance.holoNumber == self.gb_remote.dataset.first_holo_seq:#if first holo changed
            self.stackDistance.holoNumber = self.gb_remote.dataset.first_holo_seq
            #reset phase projection
            self.stackDistance.ResetPhaseProjection()
        
        #Get the new parametesrs of xy detection
        self.stackDistance.threshold = self.gb_rec.dataset.threshold
        self.stackDistance.particleconstrast = self.gb_rec.dataset.particlecontrast
        self.stackDistance.remoteKoala = self.remoteKoala
        #Get the new parameters of aberraction compensation
        self.stackDistance.filter_size_aberration = self.gb_rec.dataset.filter_size_aberration
        self.stackDistance.use_filter_size_aberration = self.gb_rec.dataset.use_filter_size_aberration
        self.stackDistance.use_filter_amplitude = self.gb_rec.dataset.use_filter_amplitude
        
        #if phase projection is not (also when it was reset if first holo changed)
        if self.stackDistance.phaseProjection1 is None:
            self.stackDistance.AcquireAndSaveStack() #peform new stack
            self.stackDistance.ZProjection() #compute zprojection
            self.stackDistance.phaseProjection1 = self.stackDistance.phaseProjection
                  
        #if derivation is used compute stack and zprojection for next holo
        if self.stackDistance.bool_deriv and self.stackDistance.phaseProjection0 is None: 
            self.stackDistance.phaseProjection0 = self.stackDistance.phaseProjection
            self.stackDistance.holoNumber += 1
            self.stackDistance.AcquireAndSaveStack()
            self.stackDistance.ZProjection()
            self.stackDistance.phaseProjection1 = self.stackDistance.phaseProjection
        #Perform threshold for xy detection according to the method defined in the class
        if self.gb_rec.dataset.bool_ThresholdMethod != -2:
            self.stackDistance.Threshold()
            self.stackDistance.ExtractXYPositionFromHolo() #Extract XY Position 
        
        #Display results in Projection image = self.derivate_data send to win_projection
        self.gb_btn.setEnabled(True)
        if self.projection_data is None:
            self.projection_data = make.image(data=self.stackDistance.projectedWithPoints)
            if self.stackDistance.bool_deriv:
                imtothreshold = np.angle(np.exp(complex(0., 1.)*(self.stackDistance.phaseProjection1-self.stackDistance.phaseProjection0)))
                imtothreshold = imtothreshold-np.mean(imtothreshold)
                self.derivate_data = make.image(data=imtothreshold)
        else:
            self.projection_data.set_data(self.stackDistance.projectedWithPoints)
            if self.stackDistance.bool_deriv:
                imtothreshold = np.angle(np.exp(complex(0., 1.)*(self.stackDistance.phaseProjection1-self.stackDistance.phaseProjection0)))
                imtothreshold = imtothreshold-np.mean(imtothreshold)
                self.derivate_data = make.image(data=imtothreshold)
        if self.win_projection is not None:
            self.win_projection.get_plot().replot()
        if self.win_derivate is not None and self.stackDistance.bool_deriv:
            self.win_derivate.get_plot().replot()
    
        self.stackDistance.holoNumber = self.gb_remote.dataset.first_holo_seq
        self.gb_rec.dataset.threshold = self.stackDistance.threshold
        self.gb_rec.get()

    def adjust_focus_change(self):
        #event if the adjust focus detection is pressed
        if self.focusPoint is None: #display error message if no point is selected
            utils.ErrorMessage('You have to select a point')
        else:
            self.update_holder_params() #update the holder params (central rec dist, range, ratio,...)
            self.adjust_detection_change() #perform new xy detection with the new parameters
            
            #define the focus parameters (method, region of interest around particule position,...)
            self.stackDistance.focusMethod = self.gb_holder.dataset.focusMethod
            self.stackDistance.roi_sample = self.gb_holder.dataset.roi_sample
            #get which particule want to be focus adjusted
            particle = self.gb_partchoice.dataset.particulechoice
            #find focus ONLY for this particule, so user has to adjust for all of them independently
            self.stackDistance.xy_position = [self.focusPoint[particle]] #focusPoint is the xy position selected by user
            self.stackDistance.zFromxy(True) #compute reconstruction distance from xy position, and plot result ("True")
            d_focus = self.stackDistance.z[0]+self.stackDistance.central_rec_dist #compute focus reconstruction distance
            
            #write computed reconstruction distance in the main panel
            self.gb_holder.dataset.central_rec_dist = d_focus
            self.gb_holder.get()
            self.central_rec_dist_change()
            #replace the focus reconstruction distance for the point (Particule choice) that is computed
            if self.focusPoint is not None:
                particle = self.gb_partchoice.dataset.particulechoice
                self.focusPoint[particle] = self.replace_rec_dist(self.focusPoint[particle], self.stackDistance.z[0]+self.stackDistance.central_rec_dist)
#            for k in range(np.shape(self.focusPoint)[0]):
#                self.focusPoint[k] = self.replace_rec_dist(self.focusPoint[k],self.stackDistance.z[k]+self.stackDistance.central_rec_dist)
        

    def replace_rec_dist(self, tup, d):
        #replace reconstruction distance in tuple self.focusPoint
        lst = list(tup)
        lst[2] = d
        return tuple(lst)
    def roi_sample_change(self):
        #event when roi sample is changed->update  class object
        self.stackDistance.roi_sample = self.gb_holder.dataset.roi_sample
        self.tracking.roi_sample = self.gb_holder.dataset.roi_sample
    def focusMethod_change(self):
        #event when focus method change->modify class object
        self.stackDistance.focusMethod = self.gb_holder.dataset.focusMethod
    
    def construct_stack_change(self):
                #event when tracking is performed by the user (Perform tracking button)
        #initialize all class
        self.process()
        self.stackDistance = self.stackDistance.from_dataset(self.db, 
                                                             self.gb_dhm.dataset, 
                                                             self.gb_sample.dataset,
                                                             self.gb_holder.dataset, 
                                                             self.gb_rec.dataset, 
                                                             self.gb_remote.dataset)       
        self.stackDistance.use_filter_size_aberration = self.gb_rec.dataset.use_filter_size_aberration
        self.stackDistance.filter_size_aberration = self.gb_rec.dataset.filter_size_aberration
        self.stackDistance.first_holo_seq = self.gb_remote.dataset.first_holo_seq
        self.stackDistance.holoNumber = self.gb_remote.dataset.first_holo_seq
        self.stackDistance.focusMethod = self.gb_holder.dataset.focusMethod
        self.stackDistance.BcgCorrection = self.gb_rec.dataset.bcgRemove
        self.stackDistance.ROI = self.db.ConfigROI 
        self.stackDistance.InitBcgPropagation()                                              
        if self.stackDistance.remoteKoala is None:        
            self.stackDistance.remoteKoala = self.remoteKoala
        self.stackDistance.pxSize = self.remoteKoala.pxSize
        self.stackDistance.ComputeMO()
        self.tracking.stackDistance = self.stackDistance
        saveAllStackValue = self.gb_remote.dataset.saveAllStack
        self.stackDistance.saveAllStack = True

        t1 = Thread(target=self.tracking.ExtractXYZPositionFromHolo(True), args=[])
        t1.start()
        t1.join()
        del t1
        self.stackDistance.saveAllStack = saveAllStackValue

    def tracking_change(self):
        #event when tracking is performed by the user (Perform tracking button)
        #initialize all class
        self.process()
        self.stackDistance = self.stackDistance.from_dataset(self.db, 
                                                             self.gb_dhm.dataset, 
                                                             self.gb_sample.dataset,
                                                             self.gb_holder.dataset, 
                                                             self.gb_rec.dataset, 
                                                             self.gb_remote.dataset)       
        self.stackDistance.use_filter_size_aberration = self.gb_rec.dataset.use_filter_size_aberration
        self.stackDistance.filter_size_aberration = self.gb_rec.dataset.filter_size_aberration
        self.stackDistance.first_holo_seq = self.gb_remote.dataset.first_holo_seq
        self.stackDistance.holoNumber = self.gb_remote.dataset.first_holo_seq
        self.stackDistance.focusMethod = self.gb_holder.dataset.focusMethod
        self.stackDistance.BcgCorrection = self.gb_rec.dataset.bcgRemove                                                       
        if self.stackDistance.remoteKoala is None:        
            self.stackDistance.remoteKoala = self.remoteKoala
        self.stackDistance.pxSize = self.remoteKoala.pxSize
        self.stackDistance.ComputeMO()
        self.stackDistance.InitBcgPropagation()
        self.tracking.stackDistance = self.stackDistance
        
        #import xy position if available and asked for
        self.tracking.importXYPosition = self.gb_rec.dataset.XYImport
        fileExists = True
        if self.tracking.importXYPosition:
            self.tracking.XYImportFilePath = self.gb_rec.dataset.XYImportFilePath
            fileExists = utils.FileExists(self.tracking.XYImportFilePath, '.txt')
            if not self.tracking.XYFileAndNumberHoloCompatible():
                return
        if fileExists: #true if no xy importation or if txt file really exists
            if self.trackingMode == 0: #tracking mode = automatic
                if self.gb_rec.dataset.bool_ThresholdMethod == -2:
                    print 'No threshold is not compatible with this tracking mode'
                    return
                #perform a new thread with extraction of XYZ position with tracking object
                t1 = Thread(target=self.tracking.ExtractXYZPositionFromHolo(False), args=[])
                t1.start()
                t1.join()
            else:# tracking mode = single particle or loop on single particule
                if self.focusPoint is not None:#test that some focus point are defined
                    #set xy initial position to stackDistance object
                    self.tracking.stackDistance.xy_position = self.focusPoint
                    if self.gb_rec.dataset.bool_ThresholdMethod == -2:
                        self.stackDistance.singleParticuleXYPositionDetectionMethod = 0
                    else:
                        self.stackDistance.singleParticuleXYPositionDetectionMethod = 1
                    #iteration for all selected single particle track
                    for k in range(np.shape(self.focusPoint)[0]):
                        #new thread to track every single particule
                        t1 = Thread(target=self.tracking.TrackSingleParticule(k+1), args=[])
                        t1.start()
                        t1.join()
                else:#no single particule exists
                    utils.ErrorMessage('Single Particle initial Position is not defined')
        else:
            utils.ErrorMessage('Importation of xy position cannot be imported as the file is not a txt file')
    def trackingMode_change(self):
        #event whe tracking mode is changed->update corresponding class objects
        self.trackingMode = self.gb_proc.dataset.trackingMode
        #if new tracking mode is automatic or single particule->reset all points
        #a new point has to be defined !!
        if self.trackingMode == 0 or self.trackingMode == 1:
            self.win_polynomial.singlepointMethod = True
            self.win_polynomial.reset_pts()
        else:
            #if loop on single particule, modify win_polynomial paramater to be able to draw several points
            self.win_polynomial.singlepointMethod = False
        
    #method that seems to be useless !!   
    def run(self):
        for k in range(10):
            if self.projection_data is None:
                self.projection_data = make.image(data=self.stackDistance.projectedWithPoints)
            else:
                self.projection_data.set_data(self.stackDistance.projectedWithPoints)
            if self.win_projection is not None:
                self.win_projection.get_plot().replot()
            time.sleep(1)
        
    def particle_tracking_change(self):
        #event when Part. tracking from coord is pressed (only available for automatic tracking)
        self.tracking.stackDistance = self.stackDistance
        self.tracking.OpenStackZparams() #Open Stack parameters to be able to do tracking from saved coordinates
        self.tracking.Computexyzposition()
    def threshold_change(self):
        #event when treshold value is changed
        self.stackDistance.threshold = self.gb_rec.dataset.threshold
        self.adjust_detection_change()#perform again adjust detection as threshold changed
    
    def XYImport_change(self):
        #event when Import XY Position is pressed
        self.tracking.importXYPosition = self.gb_rec.dataset.XYImport
    
    def XYImportFilePath_change(self):
        #event when the path of the xy import file is changed
        self.gb_rec.dataset.get()
        self.tracking.XYImportFilePath = self.gb_rec.dataset.XYImportFilePath
    def use_filter_size_aberration_change(self):
        #event when use uniform filter for aberration compensation is changed
        self.stackDistance.ResetPhaseProjection()
        self.stackDistance.use_filter_size_aberration = self.gb_rec.dataset.use_filter_size_aberration
        self.remoteKoala.OpenHologram()

    def use_filter_amplitude_change(self):
        #event when use amplitude filter is changed
        self.stackDistance.ResetPhaseProjection()
        self.stackDistance.use_filter_amplitude = self.gb_rec.dataset.use_filter_amplitude
        self.remoteKoala.OpenHologram()
        self.drawPolynomialImage() #computed again the polynomial image according to new value
        
    def bool_deriv_change(self):
        #event if derive is changed
        self.stackDistance.bool_deriv = self.gb_rec.dataset.bool_deriv
        self.stackDistance.ResetPhaseProjection()
    
    def bool_ThresholdMethod_change(self):
        #event when treshold method is changed->peform new detection
        self.stackDistance.bool_ThresholdMethod = self.gb_rec.dataset.bool_ThresholdMethod
        if self.stackDistance.bool_ThresholdMethod != -2:
            self.adjust_detection_change()
    def particlecontrast_change(self):   
        #event when particule contrast changed->perfom new detectioin
        self.stackDistance.particlecontrast = self.gb_rec.dataset.particlecontrast  
        self.adjust_detection_change()
    def it_eros_change(self):
        #event when erosion changed
        self.stackDistance.it_eros = self.gb_rec.dataset.it_eros
        self.adjust_detection_change()
    def it_dil_change(self):
        #event when dilatation change
        self.stackDistance.it_dil = self.gb_rec.dataset.it_dil
        self.adjust_detection_change()
        
    def filter_size_aberration_change(self):
        #event when filter size aberration change->change class object
        self.aberrationCorrection.filter_size_aberration = self.gb_rec.dataset.filter_size_aberration
        self.stackDistance.ResetPhaseProjection()
    
    def MODB_change(self):
        #event when MO is changed =Koala configuration is changed
        if self.db is None: #test if object self.db is None, read it if True
            self.db = SqliteDB.DBReader()
        self.db.ProductionNo = self.gb_dhm.dataset.MODB #set config number to database class object
        self.db.ParamFromChoice() #Define the database parameters from the config choice (configType, NA, Magnification,...)
#        #Write the parameters associated to the configuration to the object clss
#        self.stackDistance.NA=self.db.NA
#        self.stackDistance.CCDPixelSizeUM=self.db.CCDPixelSizeUM
#        self.stackDistance.ROI=self.db.ConfigROI
#        self.stackDistance.InitBcgPropagation() 
#        #self.stackDistance.MO=self.db.Magnification
#        self.tracking.stackDistance=self.stackDistance
        if self.remoteKoala.connected: #send Koala remote commands if remote is connected
            self.remoteKoala.config = self.db.ProductionNo #set config number to remoteKoala class object
            self.remoteKoala.OpenConfig() #open configuration
            self.remoteKoala.OpenHologram() #open hologram
            self.drawPolynomialImage()  #acquire phase image and display it in win_polynomial
        self.stackDistance.pxSize = self.remoteKoala.pxSize #set reconstruction pixel size
        #Write the parameters associated to the configuration to the object clss
        self.stackDistance.NA = self.db.NA
        self.stackDistance.CCDPixelSizeUM = self.db.CCDPixelSizeUM
        self.stackDistance.ROI = self.db.ConfigROI
        self.stackDistance.InitBcgPropagation()
        #self.stackDistance.MO = self.db.Magnification
        self.tracking.stackDistance = self.stackDistance
        self.UpdateSamplePlaneStepUM()
        #write in main panel the default reconstruction distance of the configuration
        self.gb_holder.dataset.central_rec_dist = self.remoteKoala.remoteCommands.GetRecDistCM()
        self.gb_holder.get()

    def UpdateSamplePlaneStepUM(self):
        self.stackDistance.ratio_DOF = self.gb_holder.dataset.ratio_DOF
        self.stackDistance.ComputeMO()
        samplePlanestepUM, stackHeightUM = self.stackDistance.SamplePlaneStepUM()
        self.gb_holder.dataset.samplePlaneStepUM = samplePlanestepUM
        self.gb_holder.dataset.stackHeightUM = stackHeightUM
        self.gb_holder.get()
    
    def central_rec_dist_change(self):
        #event when the central reconstruction distance is changed
        if self.remoteKoala is not None: #only if remoteKoala is not None
            d = self.gb_holder.dataset.central_rec_dist #get central rec dist.
            self.remoteKoala.remoteCommands.SetRecDistCM(d)#set the rec dist in Koala through remote
            #for tracking mode loop on particle, change the initial rec dist depending of particule choice
            if self.gb_proc.dataset.trackingMode == 2 and self.focusPoint is not None:
                particle = self.gb_partchoice.dataset.particulechoice
                self.focusPoint[particle] = self.replace_rec_dist(self.focusPoint[particle], d)
            if self.gb_rec.dataset.bcgRemove:
                self.tracking.stackDistance.BcgFromDistance(d)
            self.drawPolynomialImage()#display the result in the Image window
    
    def range_rec_dist_change(self):
        self.stackDistance.range_rec_dist = self.gb_holder.dataset.range_rec_dist
        self.UpdateSamplePlaneStepUM()
    
    def unwrap_change(self):#event when the unwrap button is pressed
        self.stackDistance.unwrap = self.gb_rec.dataset.unwrap
        self.drawPolynomialImage()
    def bcgRemove_change(self):#event when background remove button changed
        if self.tracking.stackDistance is None:
            self.process()
            self.tracking.stackDistance = self.stackDistance
        self.tracking.stackDistance.BcgCorrection = self.gb_rec.dataset.bcgRemove
        if self.tracking.stackDistance.BcgCorrection:
            d = self.gb_holder.dataset.central_rec_dist
            self.tracking.stackDistance.BcgFromDistance(d)#compute the Bcg from d
        self.drawPolynomialImage()
    
    def bcgCompute_change(self):#compute the background from first holo to last holo
        utils.ErrorMessage("You have to reset the ROI in Koala before to continue")
        self.process()
        self.stackDistance = self.stackDistance.from_dataset(self.db, 
                                                             self.gb_dhm.dataset, 
                                                             self.gb_sample.dataset, 
                                                             self.gb_holder.dataset, 
                                                             self.gb_rec.dataset, 
                                                             self.gb_remote.dataset)       
        self.stackDistance.use_filter_size_aberration = self.gb_rec.dataset.use_filter_size_aberration
        self.stackDistance.filter_size_aberration = self.gb_rec.dataset.filter_size_aberration
        self.stackDistance.first_holo_seq = self.gb_remote.dataset.first_holo_seq
        self.stackDistance.holoNumber = self.gb_remote.dataset.first_holo_seq
        self.stackDistance.focusMethod = self.gb_holder.dataset.focusMethod                                        
        if self.stackDistance.remoteKoala is None:
            self.stackDistance.remoteKoala = self.remoteKoala
        self.stackDistance.pxSize = self.remoteKoala.pxSize
        self.stackDistance.ComputeMO()
        self.stackDistance.InitBcgPropagation()#initialize the propagation tools
        self.tracking.stackDistance = self.stackDistance
        self.tracking.BcgCompute() #compute the background from reconstruction
        tmp = self.tracking.stackDistance.BcgWavefront
        del tmp
        self.remoteKoala.first_holo_seq = -1
        self.remoteKoala.OpenConfig()
        self.update_remote_params()
    
    def init_buttons(self):
        #init buttons event
        self.gb_dhm.dataset.set_loadParamsFromFile_change_cb(self.loadParamsFromFile_change)
        
        self.gb_dhm.dataset.set_saveParamsToFile_change_cb(self.saveParamsToFile_change)
        
        self.gb_rec.dataset.set_adjust_detection_change_cb(self.adjust_detection_change)
        self.gb_rec.dataset.set_adjust_aberration_change_cb(self.adjust_aberration_change)
        self.gb_holder.dataset.set_adjust_focus_change_cb(self.adjust_focus_change)
        
        self.gb_proc.dataset.set_construct_stack_change_cb(self.construct_stack_change)
        self.gb_proc.dataset.set_trackingMode_change_cb(self.trackingMode_change)
        self.gb_proc.dataset.set_tracking_change_cb(self.tracking_change)
        self.gb_proc.dataset.set_particle_tracking_change_cb(self.particle_tracking_change)
        self.gb_rec.dataset.set_threshold_change_cb(self.threshold_change)
        
        self.gb_rec.dataset.set_unwrap_change_cb(self.unwrap_change)
        self.gb_rec.dataset.set_use_filter_size_aberration_change_cb(self.use_filter_size_aberration_change)
        self.gb_rec.dataset.set_filter_size_aberration_change_cb(self.filter_size_aberration_change)
        
        self.gb_rec.dataset.set_bcgRemove_change_cb(self.bcgRemove_change)
        self.gb_rec.dataset.set_bcgCompute_change_cb(self.bcgCompute_change)
        
        self.gb_rec.dataset.set_bool_deriv_change_cb(self.bool_deriv_change)
        self.gb_rec.dataset.set_bool_ThresholdMethod_change_cb(self.bool_ThresholdMethod_change)
        self.gb_rec.dataset.set_particlecontrast_change_cb(self.particlecontrast_change)
        self.gb_rec.dataset.set_it_eros_change_cb(self.it_eros_change)
        self.gb_rec.dataset.set_it_dil_change_cb(self.it_dil_change)
        self.gb_rec.dataset.set_use_filter_amplitude_change_cb(self.use_filter_amplitude_change)
        self.gb_rec.dataset.set_XYImport_change_cb(self.XYImport_change)
        
        self.gb_holder.dataset.set_focusMethod_change_cb(self.focusMethod_change)
        self.gb_holder.dataset.set_range_rec_dist_change_cb(self.range_rec_dist_change)
        self.gb_holder.dataset.set_central_rec_dist_change_cb(self.central_rec_dist_change)
        self.gb_holder.dataset.set_ratio_DOF_change_cb(self.UpdateSamplePlaneStepUM)
        self.gb_holder.dataset.set_roi_sample_change_cb(self.roi_sample_change)
        self.gb_rec.dataset.set_XYImportFilePath_change_cb(self.XYImportFilePath_change)
        self.gb_dhm.dataset.set_MODB_change_cb(self.MODB_change)
        self.gb_partchoice.dataset.set_particulechoice_change_cb(self.particulechoice_change)
    def process(self): 
        #init the class object if there are None
        if self.db is None:
            self.db = SqliteDB.DBReader()
            self.db.ProductionNo = self.gb_dhm.dataset.MODB
            self.db.ParamFromChoice()
        if self.aberrationCorrection is None and self.remoteKoala is None and self.stackDistance is None and self.tracking is None:
            self.init_buttons()
        
        
        if self.aberrationCorrection is None:
            self.aberrationCorrection = Tracking_process.aberrationCorrection.from_dataset(self.gb_rec.dataset)
            #self.aberrationCorrection.from_dataset(self.gb_rec.dataset)
        
        if self.tracking is None:
            self.tracking = Tracking_process.tracking.from_dataset(self.gb_holder.dataset, 
                                                                   self.gb_sample.dataset)
            #self.tracking.from_dataset(self.gb_sample.dataset)
        
        if self.remoteKoala is None:
            self.remoteKoala = Tracking_process.remoteKoala.from_dataset(self.db, 
                                                                         self.gb_remote.dataset)
            #self.remoteKoala.from_dataset(self.gb_remote.dataset,self.gb_rec.dataset)
            self.remoteKoala.DefineSeqFormat()
            self.gb_remote.dataset.last_holo_seq = (int)(self.remoteKoala.framesNumber)-1
            self.gb_remote.get()
            self.tracking.last_holo_seq = self.gb_remote.dataset.last_holo_seq+1
        if self.stackDistance is None:
            
            self.stackDistance = Tracking_process.stackDistance.from_dataset(self.db, 
                                                                             self.gb_dhm.dataset, 
                                                                             self.gb_sample.dataset,
                                                                             self.gb_holder.dataset, 
                                                                             self.gb_rec.dataset,
                                                                             self.gb_remote.dataset)

            self.stackDistance.remoteKoala = self.remoteKoala
            #self.stackDistance.InitBcgPropagation() 
        if self.UnwrapK is None:
            self.UnwrapK = UnwrapK.UnwrapK()
   
    def particulechoice_change(self):
        #event when user choose another particule (only when tracking mode = Loop on single particule)
        particule = self.gb_partchoice.dataset.particulechoice
        d_focus = list(self.focusPoint[particule])[2]
        self.gb_holder.dataset.central_rec_dist = d_focus
        self.gb_holder.get()
        self.central_rec_dist_change()
   
    def _show_win(self, state, win, data, title):
        #show window for all win_xxx
        if state != 0: # show window
            if win is None:
                win = gui_construction.PhaseDialog(title,
                                                   point_selected_cb=self.drawpoint)#apply self.drapwpoint when point_selected_cb event occurs
                win.show()
                plot = win.get_plot()
                plot.add_item(data)
        else:
            if win is not None:
                win.close()
                win = None
        return win
    
    
    def drawpoint(self, point):
        #event when a point is drawn in a image window (selection of particules)
        self.stackDistance.xy_position = None #reset xy_position is stackDistance
        point += (self.gb_holder.dataset.central_rec_dist,) #add reconstruction distance to draw point
        if self.gb_proc.dataset.trackingMode != 2 or self.focusPoint is None:#reset point if not Loop on single particule
            self.resetPoint()
        self.focusPoint.append(point) #add new point (only one point for automatic or single particule tracking)
        self.stackDistance.xy_position = self.focusPoint #set xy position to stackDistance object
        newchoice = [] #init the choice of particules for Particule choice list to empty
        #create new choices for particule choice with the number of drawn particules
        for k in range(np.shape(self.focusPoint)[0]):
            choice = (k, 'Part '+str(k+1))
            newchoice.append(choice)
        #set the new list of particule choice to the list buttons Particule choice
        self.gb_partchoice.dataset.choices.set(newchoice)
        self.gb_partchoice.set()
        self.gb_partchoice.dataset.particulechoice = np.shape(newchoice)[0]-1
        self.gb_partchoice.get()
        
            
    def resetPoint(self):
        #reset all selected points automatically if tracking mode =0 or 1, or if done manually with button in image
        self.focusPoint = []
        self.stackDistance.xy_position = self.focusPoint
    def _show_winpolynomial(self, state, win, title):
        #shwo win_polynomial image and define events
        if state != 0: # show window
            if win is None:
                win = gui_construction.PolynomialDialog(title, 
                                                        hvseg_selected_cb=self.applyPolynomialCorrection, 
                                                        reset_seg_cb=self.resetSegment, 
                                                        reset_pts_cb=self.resetPoint, 
                                                        point_selected_cb=self.drawpoint, 
                                                        area_selected_cb=self.area_selected)
                win.show()
                plot = win.get_plot()
                plot.add_item(self.phase_data)
        else:
            if win is not None:
                win.close()
                win = None
        return win
    
    def area_selected(self, rect):
        #plot image of selected area, no interest now, but could be used for offset or other methods
        plt.figure(1)
        width = self.remoteKoala.remoteCommands.GetPhaseWidth()
        height = self.remoteKoala.remoteCommands.GetPhaseHeight()
        phase = self.remoteKoala.remoteCommands.getPhaseImage('lambda1', width, height)
        plt.imshow(phase[rect[0]:rect[1], rect[2]:rect[3]])
        #plt.imshow(self.phase_data[rect[0]:rect[1],rect[2]:rect[3]])
    def resetSegment(self):
        #reset Seggments in Koala when reset buttons is pressed
        self.remoteKoala.remoteCommands.ResetCorrSegments()
    def applyPolynomialCorrection(self, hv_seg):
        #apply Polynomial correction in Koala from profile drawn on win_polynomial
        self.remoteKoala.new_hv_seg(hv_seg)#add new segments in Koala phase image and perform correction
        self.drawPolynomialImage() #draw polynomial image
    def drawPolynomialImage(self):
        #display corrected phase image
        width = self.remoteKoala.remoteCommands.GetPhaseWidth()
        height = self.remoteKoala.remoteCommands.GetPhaseHeight()
        if not self.stackDistance.unwrap:#get Koala Image
            phase = self.remoteKoala.remoteCommands.getPhase32fImage('lambda1', width, height)#get phase image
        else:
            phase = self.remoteKoala.remoteCommands.getPhase32fImage('lambda1', width, height)#Get the float image to be able to peform unwrapping
            phase = self.UnwrapK.PerformUnwrap(phase)#perforum the unwrapping
        if self.gb_rec.dataset.bcgRemove: #if bck is removed
            if self.stackDistance.BcgWavefront is not None:                              
                phase -= np.angle(self.stackDistance.BcgWavefront)
                phase = np.angle(np.exp(complex(0., 1.)*phase))
            else:
                print 'bcg wavefront is none'
                return
        if self.phase_data is None:
            self.phase_data = make.image(data=phase, colormap="gray")
        else:
            self.phase_data.set_data(phase)
        if self.win_polynomial is not None:
            self.win_polynomial.get_plot().replot()
        if not self.cb_polynomial.isChecked():
            self.cb_polynomial.nextCheckState()
        
    def show_aberration(self, state):
        #show win_aberration image with aberration_data
        self.win_aberration = self._show_win(state, 
                                             self.win_aberration, 
                                             self.aberration_data, 
                                             'Aberration compensation')
        
    def show_projection(self, state):
        #show and update win_projection image
        self.win_projection = self._show_win(state, 
                                             self.win_projection, 
                                             self.projection_data, 'Projection')
        if self.stackDistance.bool_deriv:#if derive, show also win_derivate
            self.win_derivate = self._show_win(state, 
                                               self.win_derivate, 
                                               self.derivate_data,'Derivate')
    def show_polynomial(self, state):
        #Show win_polynomial
        self.win_polynomial = self._show_winpolynomial(state, 
                                                       self.win_polynomial, 
                                                       'Polynomial Correction')
            
if __name__ == '__main__':
    # Create QApplication
    from guidata.qt.QtGui import QApplication
    #import sys
    
    app = QApplication.instance()  # checks if QApplication already exists
    if not app:    # create QApplication if it doesnt exist
        app = QApplication(sys.argv)
    
    window = KMainWindow()
    window.show()
    #sys.exit(app.exec_())
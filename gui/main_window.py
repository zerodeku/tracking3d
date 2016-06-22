# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 18:19:25 2016

@author: zerodeku
"""

import sys
from threading import Thread

from guidata.qt.QtGui import QMainWindow, QWidget, QSplitter, QPushButton
from guidata.qt import QtCore
from guidata.qt.QtCore import SIGNAL
from guidata.dataset.qtwidgets import DataSetEditGroupBox

from tracking3d.gui import gui_construction
from tracking3d.core import utilsForProcessing as utils
from tracking3d.core import tracking_3d

class MainWindow(Thread, QMainWindow):
    def __init__(self):
        # init main panel
        QMainWindow.__init__(self)
        Thread.__init__(self)
        widget = QWidget()
        self.setCentralWidget(widget)
        self.setWindowTitle('3D tracking')

        # define main objects
        self.tracking = tracking_3d.Tracking3D.getInstance()
        
        # init dataset from gui_construction
        self.gb_sequence = DataSetEditGroupBox("Sequence", 
                    gui_construction.RemoteParameters, comment='')
        self.gb_sample = DataSetEditGroupBox("Sample", 
                    gui_construction.SampleParameters, comment='')
        self.gb_tracking = DataSetEditGroupBox("Tracking", 
                    gui_construction.TrackingParameters, comment='')
            
        self.btn_process = QPushButton("Start tracking", self)
        self.btn_process.clicked.connect(self.start_tracking)

        # associate events to dataset apply buttons
        self.connect(self.gb_sequence, SIGNAL("apply_button_clicked()"), 
                     self.update_remote_params)
        self.connect(self.gb_tracking, SIGNAL("apply_button_clicked()"), 
                     self.update_tracking_params)     
        self.connect(self.gb_sample, SIGNAL("apply_button_clicked()"), 
                     self.update_sample_parameters)
                     
        # organize the app panels
        splitter1 = QSplitter(QtCore.Qt.Vertical)
        splitter1.addWidget(self.gb_sequence)
        splitter1.addWidget(self.gb_sample)
        splitter1.addWidget(self.gb_tracking)
        splitter1.addWidget(self.btn_process)
        
        splitter = QSplitter(self)
        splitter.addWidget(splitter1)
        self.setCentralWidget(splitter)
        self.setContentsMargins(10, 5, 10, 5)   
        
        # get all params from datasets
        self.gb_sequence.get()
        self.gb_sample.get()
        self.gb_tracking.get()
        
    # when the main program is closed
    def closeEvent(self, event):
        remoteKoala = self.tracking.remoteKoala
        if remoteKoala is not None and remoteKoala.connected:
            self.tracking.remoteKoala.remoteCommands.ResetCorrSegments()
            self.tracking.remoteKoala.CloseConnection()

    # update the remote parameters
    def update_remote_params(self):
        utils.Log("Sequence parameters initialized")
        self.tracking.initRemoteKoala(self.gb_sequence.dataset)
        self.gb_sequence.dataset.last_holo_seq = \
                    (int)(self.tracking.remoteKoala.framesNumber) - 1
        self.gb_sequence.get()
#        self.gb_sequence.setEnabled(False)
        self.tracking.setupKoala()
    
    def update_sample_parameters(self):
        utils.Log("Sample paramters initialized")
        self.tracking.initSampleParams(self.gb_sample.dataset)
      
    # update tracking parameters 
    def update_tracking_params(self):
        utils.Log("Tracking paramters initialized")
        self.tracking.initTrackingParams(self.gb_tracking.dataset)
        self.gb_tracking.dataset.samplePlaneStepUM = self.tracking.stepUM
        self.gb_tracking.dataset.stackHeightUM = self.tracking.stackHeightUM
        self.gb_tracking.get()
        
    # initiate tracking process
    def start_tracking(self):
        utils.Log("Start tracking process")
        self.tracking.performTracking()
        
        
# main to run the application
if __name__ == '__main__':
    # create QApplication
    from guidata.qt.QtGui import QApplication
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    
    window = MainWindow()
    window.show()
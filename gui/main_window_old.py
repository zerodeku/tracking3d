# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 18:19:25 2016

@author: zerodeku
"""

import sys
from threading import Thread

from guidata.qt.QtGui import QMainWindow, QWidget, QSplitter
from guidata.qt import QtCore
from guidata.qt.QtCore import SIGNAL
from guidata.dataset.qtwidgets import DataSetEditGroupBox

from tracking3d.gui import gui_construction
from tracking3d.core import SqliteDB
from tracking3d.core import Tracking_process
from tracking3d.core import utilsForProcessing as utils

class MainWindow(Thread, QMainWindow):
    def __init__(self):
        # init main panel
        QMainWindow.__init__(self)
        Thread.__init__(self)
        widget = QWidget()
        self.setCentralWidget(widget)

        # define main objects
        self.remoteKoala = None
        self.aberrationCorrection = None
        self.stackDistance = None
        self.tracking = None
        self.db = None
        self.focusPoint = None
        self.trackingMode = 0 # 0 - Automatic, 1 - Single
        self.unwrapK = None
        self.setWindowTitle('3D tracking')
        
        # init dataset from gui_construction
        self.gb_dhm = DataSetEditGroupBox("DHM", 
                    gui_construction.DHMParameters, comment='')
        self.gb_sample = DataSetEditGroupBox("Sample", 
                    gui_construction.SampleParameters, comment='')
        self.gb_remote = DataSetEditGroupBox("Sequence", 
                    gui_construction.RemoteParameters, comment='')       
        self.gb_holder = DataSetEditGroupBox("Focus", 
                    gui_construction.HolderParameters, comment='')
        self.gb_rec = DataSetEditGroupBox("Preprocessing Parameters", 
                    gui_construction.ReconstructionParameters, comment='')
#        self.gb_proc = DataSetEditGroupBox("Process", 
#                    gui_construction.ProcessParameters, comment='')
#        self.gb_partchoice = DataSetEditGroupBox("Particle Choice", 
#                    gui_construction.ParticuleChoice, comment='')
                    
        # associate events to dataset apply buttons
        self.connect(self.gb_dhm, SIGNAL("apply_button_clicked()"), 
                     self.update_dhm_params)
        self.connect(self.gb_sample, SIGNAL("apply_button_clicked()"), 
                     self.update_sample_params)
        self.connect(self.gb_remote, SIGNAL("apply_button_clicked()"), 
                     self.update_remote_params)
        self.connect(self.gb_holder, SIGNAL("apply_button_clicked()"), 
                     self.update_holder_params)
        
        # organize the subpanels
        splitter1 = QSplitter(QtCore.Qt.Vertical)
        splitter1.addWidget(self.gb_dhm)
        splitter1.addWidget(self.gb_sample)
        splitter1.addWidget(self.gb_remote)
        splitter1.addWidget(self.gb_holder)
        
        splitter2 = QSplitter(QtCore.Qt.Vertical)
        
        splitter = QSplitter(self)
        splitter.addWidget(splitter1)
        splitter.addWidget(splitter2)
        self.setCentralWidget(splitter)
        self.setContentsMargins(10, 5, 10, 5)   
        
        # get all params from datasets
        self.gb_dhm.get()

    def update_dhm_params(self):
        utils.Log("DHM parameters updated")
        self.process()
        
#        self.stackDistance.from_dataset(self.db, self.gb_dhm.dataset)
        
    def update_sample_params(self):
        utils.Log("sample parameters updated")
        self.process()
        
    def update_remote_params(self):
        utils.Log("sequence parameters updated")
        self.process()
      
    def update_holder_params(self):
        utils.Log("focus parameter updated")
        self.process()
        
        
        
    def process(self):
        if self.db is None:
            self.db = SqliteDB.DBReader()
            self.db.ProductionNo = self.gb_dhm.dataset.MODB
            self.db.ParamFromChoice()
        
#        if self.stackDistance is None:
#            self.stackDistance = Tracking_process.stackDistance.from_dataset(
#            self.db, self.gb_dhm.dataset, None, None, None, None)

# main to run the application
if __name__ == '__main__':
    # create QApplication
    from guidata.qt.QtGui import QApplication
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    
    window = MainWindow()
    window.show()
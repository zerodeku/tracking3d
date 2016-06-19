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
        self.gb_tracking = DataSetEditGroupBox("Tracking", 
                    gui_construction.TrackingParameters, comment='')
           
        # associate events to dataset apply buttons
        self.connect(self.gb_sequence, SIGNAL("apply_button_clicked()"), 
                     self.update_remote_params)
        self.connect(self.gb_tracking, SIGNAL("apply_button_clicked()"), 
                     self.update_tracking_params)        
        
        # organize the app panels
        splitter1 = QSplitter(QtCore.Qt.Vertical)
        splitter1.addWidget(self.gb_sequence)
        splitter1.addWidget(self.gb_tracking)
        
        splitter = QSplitter(self)
        splitter.addWidget(splitter1)
        self.setCentralWidget(splitter)
        self.setContentsMargins(10, 5, 10, 5)   
        
        # get all params from datasets
#        self.gb_sequence.get()
        
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
        self.gb_sequence.setEnabled(False)
        self.tracking.setupKoala()
    
    # update tracking parameters 
    def update_tracking_params(self):
        utils.Log("Tracking paramters initialized")
                
        
        
# main to run the application
if __name__ == '__main__':
    # create QApplication
    from guidata.qt.QtGui import QApplication
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    
    window = MainWindow()
    window.show()
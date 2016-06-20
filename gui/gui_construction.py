# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 20:01:10 2016

@author: zerodeku
"""
from guiqwt.config import _
from guidata.dataset.datatypes import DataSet, ActivableDataSet, BeginGroup, EndGroup
from guidata.dataset.dataitems import ChoiceItem, FloatItem, DirectoryItem, \
IntItem, ButtonItem

from tracking3d.core import SqliteDB

class RemoteParameters(ActivableDataSet):#Sequence parameters
    # read database to know default directory
    db = SqliteDB.DBReader()
    directoryPath = DirectoryItem("Seq. dir.", default=db.DefaultDirectory)
    first_holo_seq = IntItem(_("First hologram"), default=0, min=0)
    # last holo seq is defined with framesNumberProp that is modified when a 
    # new directory is choose
    last_holo_seq = IntItem(_("Last hologram"), default=0, min=0)
    saveAllStack = ChoiceItem(_("Save All stacks"), [(True, 'Yes'), 
                              (False, 'False')], default=False)
    #create list of choice from database   
    db = SqliteDB.DBReader()
    choice = db.choiceItemDB()
    MODB = ChoiceItem(_("Magnification"), choice).set_prop("display", callback=None)               
                              
class TrackingParameters(DataSet):
    # particle contrast
    particlecontrast = ChoiceItem(_("Particle Contrast"), [(0, 'Positive'), 
                                  (1, 'Negative')]).set_prop("display", 
                                    callback=None)
    # whether to remove background
    bcgRemove = ChoiceItem(_("Background remove"), [(True, 'Yes'), (False, 'No')], 
                         default=False).set_prop("display", callback=None)
    # set threshold for xy detection
#    threshold = FloatItem(("Threshold"),
#                          default=10., slider=False, min=0., 
#                          max=255, step=1,).set_prop("display", callback=None)
    # choose threshold method
    bool_ThresholdMethod = ChoiceItem(_("Threshold Method"), [(-2, 'No threshold/max value'), 
                                      (-1, 'Manual'), (0, 'Huang'), (2, 'Intermodes'), 
                                (3, 'ReniEntropy'), (4, 'Triangle'), (5, 'Otsu'), (6, 'Yen'), 
                                (7, 'Moments'), (8, 'IsoData'), (9, 'Peak Detection')], 
                            default=0).set_prop("display", callback=None)
            #Threshold method
            #-1=manual
            #0=Huang
            #1=Intermodes
            #2=MaxEntropy
            #3=ReniEntropy
            #4=Triangle
            #5=Otsu
            #6=Yen
            #7=Moments
            #8=IsoData
            #9=Peak Detection
    #erosion and dilatation iterations
    it_eros = IntItem(_("Iter. erosion"),
                          default=1, slider=False, unit='pixels', 
                          min=0).set_prop("display", callback=None)
    it_dil = IntItem(_("Iter. dilatation"),
                          default=4, slider=False, unit='pixels', 
                          min=0).set_prop("display", callback=None)

class ProcessCommand():
    def _tracking_change(self, item, value, parent):
        if self.tracking_change_cb is not None:
            self.tracking_change_cb()
    def set_tracking_change_cb(self, cb):
        self.tracking_change_cb = cb
    processCommand = ButtonItem("Start tracking", callback=_tracking_change)
    tracking_change_cb = None


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
#    saveAllStack = ChoiceItem(_("Save All stacks"), [(True, 'Yes'), 
#                              (False, 'False')], default=False)
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
                         default=True).set_prop("display", callback=None)
    # set threshold for xy detection
#    threshold = FloatItem(("Threshold"),
#                          default=10., slider=False, min=0., 
#                          max=255, step=1,).set_prop("display", callback=None)
#    # choose threshold method
#    bool_ThresholdMethod = ChoiceItem(_("Threshold Method"), [(-2, 'No threshold/max value'), 
#                                      (-1, 'Manual'), (0, 'Huang'), (2, 'Intermodes'), 
#                                (3, 'ReniEntropy'), (4, 'Triangle'), (5, 'Otsu'), (6, 'Yen'), 
#                                (7, 'Moments'), (8, 'IsoData'), (9, 'Peak Detection')], 
#                            default=0).set_prop("display", callback=None)
#            #Threshold method
#            #-1=manual
#            #0=Huang
#            #1=Intermodes
#            #2=MaxEntropy
#            #3=ReniEntropy
#            #4=Triangle
#            #5=Otsu
#            #6=Yen
#            #7=Moments
#            #8=IsoData
#            #9=Peak Detection
#    #erosion and dilatation iterations
#    it_eros = IntItem(_("Iter. erosion"),
#                          default=1, slider=False, unit='pixels', 
#                          min=0).set_prop("display", callback=None)
#    it_dil = IntItem(_("Iter. dilatation"),
#                          default=4, slider=False, unit='pixels', 
#                          min=0).set_prop("display", callback=None)
    central_rec_dist = FloatItem(_("Central Rec. dist."), default=-1.3, 
            min=-100, max=100, unit='cm').set_prop("display", callback=None)
    range_rec_dist = FloatItem(_("Range Rec. dist."), default=1, min=-100, max=100, 
                           unit='cm').set_prop("display", callback=None)
    ratio_DOF = FloatItem(_("Ratio of DOF"), default=1, 
                      min=0.1).set_prop("display", callback=None)
    samplePlaneStepUM = FloatItem(_("Step size along z"), default=0, 
                                  unit='um').set_prop("display", active=False)
    stackHeightUM = FloatItem(_("Stack Height"), default=0, 
                              unit='um').set_prop("display", active=False)

class SampleParameters(DataSet):#Focus parameters
    index_medium = FloatItem(_("Medium index"), default=1.33, min=1.)
    index_chamber = FloatItem(_("Chamber index"), default=1.52, min=1.)
    index_sample = FloatItem(_("Sample index"), default=1.39, min=1.)
    particle_size = FloatItem(_("Particle diameter (um)"), default=1., min=1.)
    max_speed = FloatItem(_("Max speed (um/s)"), default=30., min=1.)
 

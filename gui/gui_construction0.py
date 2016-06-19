# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 12:56:39 2014

@author: tcolomb
"""

from guiqwt.config import _
#from guiqwt.tools import RectangleTool, CircleTool, ColormapTool, SelectTool, SegmentTool, DisplayCoordsTool
from guiqwt.tools import SegmentTool, PointTool, RectangleTool
from guiqwt.plot import ImageDialog
from guiqwt.panels import ID_XCS, ID_YCS, ID_OCS, ID_ITEMLIST, ID_CONTRAST
from guidata.qt.QtGui import QPushButton, QDialogButtonBox
#from guiqwt.interfaces import IShapeItemType
#from guidata.dataset.datatypes import DataSet, ActivableDataSet, ValueProp, FuncProp, BeginGroup, EndGroup
from guidata.dataset.datatypes import DataSet, ActivableDataSet, BeginGroup, EndGroup, ValueProp, GetAttrProp, FuncProp
from guidata.dataset.dataitems import ChoiceItem, FloatItem, IntItem, ButtonItem, DirectoryItem, FileOpenItem
from tracking3d.core import SqliteDB
import numpy as np
detectionprop = False
framesNumberprop = 0
propImport = ValueProp(False)

#class to define polynomial correction segments
class hv_segment:
    def __init__(self, left, top, length, orient):
        self.left = left
        self.top = top
        self.length = length
        self.orient = orient # true -> vertical, false -> horizontal

#class to be able to modify the list of buttons dynamically (used for Particule choice)
class ChoicesVariable(object):
    def __init__(self):
        self.choices = [(0, 'Part1', None)] #default value is a single particule
    def set(self, choices):
        normalized = []
        for idx, c in enumerate(choices):
            normalized.append(self._normalize_choice(idx, c))
        self.choices = normalized
    def __call__(self, *args):
        return self.choices
    def _normalize_choice(self, idx, choice_tuple):
        img = None
        if isinstance(choice_tuple, tuple):
            if len(choice_tuple) == 2:
                key, value = choice_tuple
            else:
                key, value, img = choice_tuple
        else:
            key = idx
            value = choice_tuple
        if isinstance(value, str):
            value = str(value)
        return(key, value, img)
           
class DHMParameters(DataSet):
    #Magnification button
    def _MODB_change(self, item, value):
        if self.MODB_change_cb is not None:
            self.MODB_change_cb()
    def set_MODB_change_cb(self, cb):
        self.MODB_change_cb = cb   
    #create list of choice from database   
    db = SqliteDB.DBReader()
    choice = db.choiceItemDB()
    MODB = ChoiceItem(_("Magnification"), choice).set_prop("display", callback=_MODB_change)
    MODB_change_cb = None
    
    #Save Parameters buttons definition and associated events
    def _saveParamsToFile_change(self, item, value, parent):
        if self.saveParamsToFile_change_cb is not None:
            self.saveParamsToFile_change_cb()
    def set_saveParamsToFile_change_cb(self, cb):
        self.saveParamsToFile_change_cb = cb 
    saveParamsToFile = ButtonItem("Save Parameters to file", callback=_saveParamsToFile_change)
    saveParamsToFile_change_cb = None
    
    #Load Parameters buttons definition and associated events
    def _loadParamsFromFile_change(self, item, value, parent):
        if self.loadParamsFromFile_change_cb is not None:
            self.loadParamsFromFile_change_cb()
    def set_loadParamsFromFile_change_cb(self, cb):
        self.loadParamsFromFile_change_cb = cb 
    loadParamsFromFile = ButtonItem("Load Saved Parameters", callback=_loadParamsFromFile_change)
    loadParamsFromFile_change_cb = None

class RemoteParameters(ActivableDataSet):#Sequence parameters
    #read database to know default directory
    db = SqliteDB.DBReader()
    directoryPath = DirectoryItem("Seq. dir.", default=db.DefaultDirectory)
    first_holo_seq = IntItem(_("First hologram"), default=0, min=0)
    #last holo seq is defined with framesNumberProp that is modified when a new directory is choose
    last_holo_seq = IntItem(_("Last hologram"), default=framesNumberprop)
    saveAllStack = ChoiceItem(_("Save All stacks"), [(True, 'Yes'), (False, 'False')], default=False)


class SampleParameters(DataSet):#Sample parameters
    index_medium = FloatItem(_("Medium index"), default=1.33, min=1.)
    index_chamber = FloatItem(_("Chamber index closing the chamber"), default=1.52, min=1.)
    index_sample = FloatItem(_("Sample index"), default=1., min=1.)
    max_displacement = FloatItem(_("Max Displacement"), default=20, min=1., unit='pixels')

class HolderParameters(DataSet):#Focus parameters
    #Define central rec dis button and associated event when changed
    def _central_rec_dist_change(self, item, value):
        if self.central_rec_dist_change_cb is not None:
            self.central_rec_dist_change_cb()
    def set_central_rec_dist_change_cb(self, cb):
        self.central_rec_dist_change_cb = cb  
    central_rec_dist = FloatItem(_("Central Rec. dist."), default=-1.3, min=-100., max=100, 
                                 unit='cm').set_prop("display", callback=_central_rec_dist_change)
    central_rec_dist_change_cb = None
    
    def _range_rec_dist_change(self, item, value):
        if self.range_rec_dist_change_cb is not None:
            self.range_rec_dist_change_cb()
    def set_range_rec_dist_change_cb(self, cb):
        self.range_rec_dist_change_cb = cb  
    range_rec_dist = FloatItem(_("Range Rec. dist."), default=1, min=-100., max=100, 
                               unit='cm').set_prop("display", callback=_range_rec_dist_change)
    range_rec_dist_change_cb = None
    
    def _ratio_DOF_change(self, item, value):
        if self.ratio_DOF_change_cb is not None:
            self.ratio_DOF_change_cb()
    def set_ratio_DOF_change_cb(self, cb):
        self.ratio_DOF_change_cb = cb  
    ratio_DOF = FloatItem(_("Ratio of DOF"), default=1, 
                          min=0.1).set_prop("display", callback=_ratio_DOF_change)
    ratio_DOF_change_cb = None
    
    samplePlaneStepUM = FloatItem(_("Step size in sample plane"), default=0, 
                                  unit='um').set_prop("display", active=False)
    stackHeightUM = FloatItem(_("Stack Height"), default=0, 
                              unit='um').set_prop("display", active=False)
    #Define roi sample button and associated event
    def _roi_sample_change(self, item, value):
        if self.roi_sample_change_cb is not None:
            self.roi_sample_change_cb()
    def set_roi_sample_change_cb(self, cb):
        self.roi_sample_change_cb = cb   
    roi_sample = IntItem(_("ROI size (to find focus)"), default=30, min=1, 
                         unit='pixels').set_prop("display", callback=_roi_sample_change)
    roi_sample_change_cb = None
    
    #Define default focus method and associated event when changed
    def _focusMethod_change(self, item, value):
        if self.focusMethod_change_cb is not None:
            self.focusMethod_change_cb()
    def set_focusMethod_change_cb(self, cb):
        self.focusMethod_change_cb = cb   
    
    focusMethod = ChoiceItem(_("Focus Method"), [(0, 'Min amp std'), (1, 'Absolute Min Amp'), 
                            (2, 'Absolute Max Phase'), (3, 'Integral Density'), 
                            (4, 'Max phase(valuePos-mean)'), (5, 'Max std dev'), (6, 'Skewness'), 
                            (7, 'Max std from Sobel'), 
                            (8, 'Max std from Sobel on phase image')]).set_prop("display", callback=_focusMethod_change)
    focusMethod_change_cb = None
    
    #Event to adjust the focus detection an init of button
    def _adjust_focus_change(self, item, value, parent):
        if self.adjust_focus_change_cb is not None:
            self.adjust_focus_change_cb()
    def set_adjust_focus_change_cb(self, cb):
        self.adjust_focus_change_cb = cb   
    focus_perform = ButtonItem("Adjust Focus Detection", callback=_adjust_focus_change)
    adjust_focus_change_cb = None
    
class ReconstructionParameters(ActivableDataSet):#Preprocessing Parameters
         
    #events associated to different parameters
    def _threshold_change(self, item, value):
        if self.threshold_change_cb is not None:
            self.threshold_change_cb()
    def set_threshold_change_cb(self, cb):
        self.threshold_change_cb = cb  
        
    def _unwrap_change(self, item, value):
        if self.unwrap_change_cb is not None:
            self.unwrap_change_cb()
    def set_unwrap_change_cb(self, cb):
        self.unwrap_change_cb = cb 
    
    def _use_filter_size_aberration_change(self, item, value):
        if self.use_filter_size_aberration_change_cb is not None:
            self.use_filter_size_aberration_change_cb()
    def set_use_filter_size_aberration_change_cb(self, cb):
        self.use_filter_size_aberration_change_cb = cb 
        
    def _use_filter_amplitude_change(self, item, value):
        if self.use_filter_amplitude_change_cb is not None:
            self.use_filter_amplitude_change_cb()
    def set_use_filter_amplitude_change_cb(self, cb):
        self.use_filter_amplitude_change_cb = cb 
        
    def _filter_size_aberration_change(self, item, value):
        if self.filter_size_aberration_change_cb is not None:
            self.filter_size_aberration_change_cb()
    def set_filter_size_aberration_change_cb(self, cb):
        self.filter_size_aberration_change_cb = cb  
        
    def _bool_deriv_change(self, item, value):
        if self.bool_deriv_change_cb is not None:
            self.bool_deriv_change_cb()
    def set_bool_deriv_change_cb(self, cb):
        self.bool_deriv_change_cb = cb
        
    def _XYImport_change(self, item, value):
        if self.XYImport_change_cb is not None:
            self.XYImport_change_cb()
    def set_XYImport_change_cb(self, cb):
        self.XYImport_change_cb = cb   
        
    def _bool_ThresholdMethod_change(self, item, value):
        if self.bool_ThresholdMethod_change_cb is not None:
            self.bool_ThresholdMethod_change_cb()
    def set_bool_ThresholdMethod_change_cb(self, cb):
        self.bool_ThresholdMethod_change_cb = cb
    
    def _particlecontrast_change(self, item, value):
        if self.particlecontrast_change_cb is not None:
            self.particlecontrast_change_cb()
    def set_particlecontrast_change_cb(self, cb):
        self.particlecontrast_change_cb = cb
        
    def _it_eros_change(self, item, value):
        if self.it_eros_change_cb is not None:
            self.it_eros_change_cb()
    def set_it_eros_change_cb(self, cb):
        self.it_eros_change_cb = cb 
        
    def _it_dil_change(self, item, value):
        if self.it_dil_change_cb is not None:
            self.it_dil_change_cb()
    def set_it_dil_change_cb(self, cb):
        self.it_dil_change_cb = cb
        
    def _XYImportFilePath_change(self, item, value):
        if self.XYImportFilePath_change_cb is not None:
            self.XYImportFilePath_change_cb()
    def set_XYImportFilePath_change_cb(self, cb):
        self.XYImportFilePath_change_cb = cb
    
    #Initialization of buttons in different groups
    unwrap = ChoiceItem(_("Unwrap"), [(True, 'Yes'), (False, 'False')], 
                        default=False).set_prop("display", callback=_unwrap_change)    
    unwrap_change_cb = None
    #Group Filtering
    Filtering_group = BeginGroup("Filtering")     
    use_filter_size_aberration = ChoiceItem(_("Use uniform filter for aberration compensation"), 
                                            [(True, 'Yes'), (False, 'False')], 
                                             default=False).set_prop("display", 
                                            callback=_use_filter_size_aberration_change)    
    use_filter_size_aberration_change_cb = None
    filter_size_aberration = IntItem(_("Filter size (aberr. compens.)"),
                          default=100, slider=False, unit='pixels', 
                          min=1, max=300).set_prop("display", 
                            callback=_filter_size_aberration_change)
    filter_size_aberration_change_cb = None
    
    use_filter_amplitude = ChoiceItem(_("Use amplitude Filter"), [(False, 'No'), 
                                      (True, 'Yes')], default=False).set_prop("display", 
                                        callback=_use_filter_amplitude_change) 
    use_filter_amplitude_change_cb = None
    
    
    def _adjust_aberration_change(self, item, value, parent):
        if self.adjust_aberration_change_cb is not None:
            self.adjust_aberration_change_cb()
    def set_adjust_aberration_change_cb(self, cb):
        self.adjust_aberration_change_cb = cb    

    aberration_perform = ButtonItem("Adjust Filtering", callback=_adjust_aberration_change)
    adjust_aberration_change_cb = None
    
    
    def _bcgRemove_change(self, item, value):
        if self.bcgRemove_change_cb is not None:
            self.bcgRemove_change_cb()
    def set_bcgRemove_change_cb(self, cb):
        self.bcgRemove_change_cb = cb
        
    def _bcgCompute_change(self, item, value, parent):
        if self.bcgCompute_change_cb is not None:
            self.bcgCompute_change_cb()
    def set_bcgCompute_change_cb(self, cb):
        self.bcgCompute_change_cb = cb
        
    _Filtering_group = EndGroup("Filtering")

    Filtering_group = BeginGroup("Background")
    bcgRemove = ChoiceItem(_("Background remove"), [(True, 'Yes'), (False, 'False')], 
                         default=False).set_prop("display", callback=_bcgRemove_change)    
    bcgRemove_change_cb = None
    
    bcgCompute = ButtonItem("Compute Background", callback=_bcgCompute_change)
    bcgCompute_change_cb = None
    
    _Filtering_group = EndGroup("Background") 
    
    #Import xy (enable/disable importation of file or detection parameters)
    propImport = GetAttrProp("XYImport")# to enable/disable Import XY File path if True or False

    XYImport = ChoiceItem(_("Import XY Position"), [(False, 'No'), (True, 'Yes')], 
                          default=False).set_prop("display", 
                            callback=_XYImport_change, store=propImport)    
    XYImport_change_cb = None
    
    XYImportFilePath = FileOpenItem("File Path", formats='txt', 
                                    default=r'D:\Mesures\3DTracking_Donner\15.05.03_46.68\Sequence_import_xy').set_prop("display", active=FuncProp(propImport, lambda x: x == True), callback=_XYImportFilePath_change)
    XYImportFilePath_change_cb = None
    
    #Group Detection parameters
    Detection_group = BeginGroup("Detection XY").set_prop("display", active=FuncProp(propImport, lambda x: x == False))
    bool_deriv = ChoiceItem(_("Derive"), [(True, 'Yes'), (False, 'No')], 
                            default=False).set_prop("display", callback=_bool_deriv_change)
    bool_deriv_change_cb = None
    particlecontrast = ChoiceItem(_("Particle Contrast"), [(0, 'Positive'), 
                                  (1, 'Negative')]).set_prop("display", 
                                    callback=_particlecontrast_change)
    particlecontrast_change_cb = None
    threshold = FloatItem(("Threshold"),
                          default=10., slider=False, min=0., 
                          max=255, step=1,).set_prop("display", callback=_threshold_change)
    threshold_change_cb = None
    bool_ThresholdMethod = ChoiceItem(_("Threshold Method"), [(-2, 'No threshold/max value'), 
                                      (-1, 'Manual'), (0, 'Huang'), (2, 'Intermodes'), 
                                (3, 'ReniEntropy'), (4, 'Triangle'), (5, 'Otsu'), (6, 'Yen'), 
                                (7, 'Moments'), (8, 'IsoData'), (9, 'Peak Detection')], 
                            default=0).set_prop("display", callback=_bool_ThresholdMethod_change)
    bool_ThresholdMethod_change_cb = None
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
                          min=0).set_prop("display", callback=_it_eros_change)
    it_eros_change_cb = None
    it_dil = IntItem(_("Iter. dilatation"),
                          default=4, slider=False, unit='pixels', 
                          min=0).set_prop("display", callback=_it_dil_change)
    it_dil_change_cb = None
    
    def _adjust_detection_change(self, item, value, parent):
        if self.adjust_detection_change_cb is not None:
            self.adjust_detection_change_cb()
    def set_adjust_detection_change_cb(self, cb):
        self.adjust_detection_change_cb = cb   
    detection_perform = ButtonItem("Adjust Detection", callback=_adjust_detection_change)
    adjust_detection_change_cb = None
    _Detection_group = EndGroup("Detection XY")

class ParticuleChoice(DataSet):#Particule choice dataset
    #Has to be independent of any other dataset for unknown error when we want to change dynamically the choices
    choices = ChoicesVariable() #defined from Class ChoiceVariable used to modify dynamically the list
    def _particulechoice_change(self, item, value):
        if self.particulechoice_change_cb is not None:
            self.particulechoice_change_cb()
    def set_particulechoice_change_cb(self, cb):
        self.particulechoice_change_cb = cb
    particulechoice = ChoiceItem(_("Particle Choice"), choices).set_prop("display", callback=_particulechoice_change)   
    particulechoice_change_cb = None 
    
class ProcessParameters(DataSet):#Process buttons          
    #reconstruct only the stack without tracking procedure
    def _construct_stack_change(self, item, value, parent):
        if self.construct_stack_change_cb is not None:
            self.construct_stack_change_cb()
    def set_construct_stack_change_cb(self, cb):
        self.construct_stack_change_cb = cb    

    construct_stack = ButtonItem("Compute only stacks", callback=_construct_stack_change)
    construct_stack_change_cb = None
    
    
    #Tracking mode choice
    def _trackingMode_change(self, item, value):
        if self.trackingMode_change_cb is not None:
            self.trackingMode_change_cb()
    def set_trackingMode_change_cb(self, cb):
        self.trackingMode_change_cb = cb
    
    trackingMode = ChoiceItem(_("Tracking Mode"), [(0, 'Automatic'), (1, 'Single Particle'), 
                              (2, 'Loop on Particles')], 
                                default=1).set_prop("display", callback=_trackingMode_change)
    trackingMode_change_cb = None 
   
    #perform tracking according to tracking mode
    def _tracking_change(self, item, value, parent):
        if self.tracking_change_cb is not None:
            self.tracking_change_cb()
    def set_tracking_change_cb(self, cb):
        self.tracking_change_cb = cb    

    tracking_perform = ButtonItem("Perform tracking", callback=_tracking_change)
    tracking_change_cb = None
    
    #tracking from text file (only in automatic)
    def _particle_tracking_change(self, item, value, parent):
        if self.particle_tracking_change_cb is not None:
            self.particle_tracking_change_cb()
    def set_particle_tracking_change_cb(self, cb):
        self.particle_tracking_change_cb = cb    

    particle_tracking_perform = ButtonItem("Part. tracking from coord.", callback=_particle_tracking_change)
    particle_tracking_change_cb = None

class PhaseDialog(ImageDialog):#Class to define Image windows
    def __init__(self, title, point_selected_cb=None):
        super(PhaseDialog, self).__init__(wintitle=title, toolbar=True, edit=True)
        self.point_selected_cb = point_selected_cb
    def register_tools(self):
        ImageDialog.register_tools(self)
        self.add_tool(PointTool, handle_final_shape_cb=self.draw_point)
        #self.add_tool(CrossSectionTool)
    def draw_point(self, shape):
        r = shape.get_pos()
        point = (r[1], r[0])
        if self.point_selected_cb is not None:
            self.point_selected_cb(point)
    def create_plot(self, options):
        ImageDialog.create_plot(self, options)

class PolynomialDialog(ImageDialog): #class to define polynomial window
    def __init__(self, title, hvseg_selected_cb=None, reset_seg_cb=None, reset_pts_cb=None, 
                 reset_rec_cb=None, point_selected_cb=None, area_selected_cb=None):
        super(PolynomialDialog, self).__init__(wintitle=title, toolbar=True, edit=True)
        #define event
        self.hvseg_selected_cb = hvseg_selected_cb
        self.reset_seg_cb = reset_seg_cb
        self.reset_pts_cb = reset_pts_cb
        self.reset_rec_cb = reset_rec_cb
        self.point_selected_cb = point_selected_cb
        self.area_selected_cb = area_selected_cb
        self.itemsListToDelete = None
        self.plot = self.panels[ID_ITEMLIST].listwidget.plot
        self.singlepointMethod = True
    
        
    def install_button_layout(self): 
        #install new button on the ImageDialog
        ResetButton = QPushButton("Reset Seg.")
        ResetButton.setDefault(False)
        
        ResetPointsButton = QPushButton("Delete all Points")
        ResetPointsButton.setDefault(False)
        
        ResetRectangleButton = QPushButton("Delete all Rectangle")
        ResetRectangleButton.setDefault(False)
        
        buttonBox = QDialogButtonBox()
        buttonBox.addButton(ResetButton, QDialogButtonBox.ActionRole)
        buttonBox.addButton(ResetPointsButton, QDialogButtonBox.ActionRole)
        buttonBox.addButton(ResetRectangleButton, QDialogButtonBox.ActionRole)
        
        #Define event according to button pressed
        ResetButton.clicked.connect(self.reset_seg)
        ResetPointsButton.clicked.connect(self.reset_pts)
        ResetRectangleButton.clicked.connect(self.reset_rec)
        self.button_layout.addWidget(buttonBox)
        
        
    def register_tools(self):
        ImageDialog.register_tools(self)
        self.add_tool(SegmentTool, handle_final_shape_cb=self.draw_seg)
        self.DrawPoint = PointTool
        self.add_tool(self.DrawPoint, handle_final_shape_cb=self.draw_point)
        self.add_tool(RectangleTool, handle_final_shape_cb=self.draw_rect)

    def draw_point(self, shape):
        #draw point event
        self.reset_last_pt() #reset last point drawn (do nothing if loop on single part.)
        r = shape.get_pos()
        point = (r[1], r[0])
        if self.point_selected_cb is not None:
            self.point_selected_cb(point)
    
    def reset_last_pt(self):
        if self.singlepointMethod:
            self.listItems('PointShape')
            if len(self.itemsListToDelete) == 2:
                point = [self.itemsListToDelete[0]]
                self.plot.del_items(point)
                self.plot.replot()
    
    def reset_pts(self): #reset point event
        self.listItems('PointShape')
        self.plot.del_items(self.itemsListToDelete) #delete points in the image
        self.plot.replot()
        if self.reset_pts_cb is not None:
            self.reset_pts_cb()
    
    def reset_seg(self): #delete all profiles in the image
        self.listItems('SegmentShape')
        self.plot.del_items(self.itemsListToDelete)
        self.plot.replot()
        if self.reset_seg_cb is not None:
            self.reset_seg_cb()
    
    def reset_rec(self): #delete all rectangles in the image
        self.listItems('RectangleShape')
        self.plot.del_items(self.itemsListToDelete)
        self.plot.replot()
        if self.reset_rec_cb is not None:
            self.reset_rec_cb()
    
    def listItems(self, Type): #delete all Items of the given Type
        items = self.panels[ID_ITEMLIST].listwidget.plot.items[:]
        self.itemsListToDelete = []
        for item in items:
            name = type(item).__name__
            if name == Type:
                self.itemsListToDelete.append(item)      
    
    def create_plot(self, options):
        ImageDialog.create_plot(self, options)
        
    def draw_rect(self, shape): #event when a rectangle is drawn
        r = shape.get_rect()
        left = min(r[1], r[3])
        up = min(r[0], r[2])
        height = np.abs(r[0]-r[2])
        width = np.abs(r[1]-r[3])
        rectModify = [left, left+width, up ,up+height]
        #allows to send event to main code
        if self.area_selected_cb is not None:
            self.area_selected_cb(rectModify)
            
    def draw_seg(self, shape):#event when profile is drawn
        #Convert python segment format to remote Koala  segment format
        r = shape.get_bounding_rect_coords()
        # layout of r is x0 y0 x1 y1
        dx = abs(r[0] - r[2])
        dy = abs(r[1] - r[3])
        if dx > dy: # horizontal
            orient = False
            left = int(min(r[0], r[2]))
            top = int(r[1])
            length = int(abs(r[0] - r[2]))
            shape.set_rect(r[0], r[1], r[2], r[1])
        else: # vertical
            orient = True
            top = int(min(r[1], r[3]))
            left = int(r[0])
            length = int(abs(r[1] - r[3]))
            shape.set_rect(r[0], r[1], r[0], r[3])
            
        if self.hvseg_selected_cb is not None:
            hv = hv_segment(left, top, length, orient)
            self.hvseg_selected_cb(hv)

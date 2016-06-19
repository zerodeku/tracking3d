# -*- coding: utf-8 -*-
"""
Created on Tue Feb 02 11:09:25 2016

@author: tcolomb
"""

import subprocess
import numpy as np
from tracking3d.core import binkoala
import os

#class to perform 2D unwrapping using UnwrapK.exe
class UnwrapK:
    def __init__(self):
        self.directoryPath = os.path.abspath(os.path.join(os.path.dirname(__file__), 
                                                          os.pardir))+'\\Unwrap'
        self.NMPMask = None #Definition of a mask for non measured point, not use here
        self.pxSize = 1.#needed to save bin file, but value not important for the unwrap
        self.hconv = 1.#needed to save bin file, but value not important for the unwrap
        #define the path of UnwrapK.exe
        self.UnwrapKPath = os.path.abspath(os.path.join(os.path.dirname(__file__), 
                                                        os.pardir))+'\\Unwrap\\UnwrapK.exe'
        
        
    def PerformUnwrap(self, phase):
        if self.pxSize is not None and self.hconv is not None and phase is not None:
            width = np.shape(phase)[0]
            height = np.shape(phase)[1]
            phName = self.directoryPath+'\\temp.bin'
            if self.NMPMask is not None:
                phase *= self.NMPMask
            binkoala.write_mat_bin(phName, phase, width, height, self.pxSize*1e-6, self.hconv, 1)
            subprocess.Popen([self.UnwrapKPath, phName, phName], universal_newlines=True, 
                             stdout=subprocess.PIPE).communicate()
            unwrapped_phase, header = binkoala.read_mat_bin(phName)
            return unwrapped_phase
        else:
            print 'pxSize, hconv or phase is None'
    
    def PerformUnwrapFromFile(self, phName):
        subprocess.Popen([self.UnwrapKPath, phName, phName], universal_newlines=True, 
                                stdout=subprocess.PIPE).communicate()
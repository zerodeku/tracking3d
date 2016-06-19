# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 14:33:33 2014

@author: naspert
"""

import numpy
from scipy import ndimage
#import matplotlib

class apodization:
    # key points in the profile (normalized for L=512)
    X1 = 23
    X2 = 69

    # coefficients of the polynomial defining the profile
    A3 = -4136.7344456
    A2 = 433.60302457

    B3 = 344.72787047
    B2 = -123.88657845
    B1 = 13.913043478
    B0 = 0.5


    def profile(self, len):
        x2 = self.X2*len/512
        
        x1 = self.X1*len/512
        
        idx = numpy.arange(0, float(len))/(len - 1)
        x = numpy.ones((1, len))
        x[:, 0:x1] = idx[0:x1]*idx[0:x1]*(self.A2 + self.A3*idx[0:x1])
        
        idx[x1:x2] = idx[x1:x2] - self.X1/512.
        
        x[:, x1:x2] = self.B0 + self.B1*idx[x1:x2] + self.B2*idx[x1:x2]*idx[x1:x2] + self.B3*idx[x1:x2]*idx[x1:x2]*idx[x1:x2]
        x[:, len - x2 - 1:len] = x[:, x2::-1]
        return x
        
    def __init__(self, width, height):
        self.width = width
        self.height = height
        xp = self.profile(width)
        yp = self.profile(height)
        self.a = numpy.outer(yp[0, :], xp[0, :])
        #matplotlib.pyplot.imshow(self.a)
        
    def __call__(self, input_holo):        
        return numpy.multiply(self.a, input_holo)
    
class zerosup:
    def __init__(self, width, height):
        self.width = width
        self.height = height
    def __call__(self, input_img):
        return input_img - ndimage.uniform_filter(input_img, 5)
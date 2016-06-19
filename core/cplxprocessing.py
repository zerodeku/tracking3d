# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 16:21:05 2014

@author: naspert
"""

import numpy
from numpy.fft import fft2, fftshift, ifft2
from scipy import ndimage

from guiqwt.builder import make
#import pykoala.gui.gui_helpers
from tracking3d.core import binkoala
#import pdb

    
class fourier_filter_descriptor:
    shape = -1
    cutband = False
    intersection = True
    def cutfilter(self, mask):
        raise NotImplementedError('Method not overridden in derived class...')
        
    def fromstring(self, str_desc):
        str_split = str_desc.split(' ')
        item_split = map(int, str_split)
        return self.fromlist(item_split)
        
    def fromlist(self, list_desc):
        if list_desc[0] == 0:
            return fourier_filter_rect(list_desc[1], list_desc[2], list_desc[3], list_desc[4], 
                                       list_desc[5], list_desc[6])
        elif list_desc[0] == 1:
            return fourier_filter_ell(list_desc[1], list_desc[2], list_desc[3], list_desc[4], 
                                      list_desc[5], list_desc[6])
        else:
            raise ValueError('Invalid shape specified %d', list_desc[0])
    def __init__(self):
        return
                                               
class fourier_filter_rect(fourier_filter_descriptor):
    def __init__(self, left, top, width, height, cutband, intersection):
        #pdb.set_trace()
        fourier_filter_descriptor.__init__(self)
        self.shape = 0
        self.left = left
        self.top = top
        self.width = width
        self.height = height
        self.cutband = cutband
        self.intersection = intersection
    def __str__(self):
        return '%d %d %d %d %d %d %d' % (self.shape, self.left, self.top, self.width, 
                                         self.height, self.cutband, self.intersection)
    def cutfilter(self, mask):
        if self.cutband:
            tmp = numpy.ones(numpy.shape(mask))
            tmp[self.top:self.top+self.height, self.left:self.left+self.width] = 0.
        else:
            tmp = numpy.zeros(numpy.shape(mask))
            tmp[self.top:self.top+self.height, self.left:self.left+self.width] = 1.
        if self.intersection:
            mask = mask*tmp
        else:
            mask = mask + tmp
            idx = mask > 1.
            mask[idx] = 1.
        return mask  
        
class fourier_filter_ell(fourier_filter_descriptor):
    def __init__(self, centerx, centery, radiusx, radiusy, cutband, intersection):
        #pdb.set_trace()
        fourier_filter_descriptor.__init__(self)
        self.shape = 1
        self.centerx = centerx
        self.centery = centery
        self.radiusx = radiusx
        self.radiusy = radiusy
        self.cutband = cutband
        self.intersection = intersection
    def __str__(self):
        return '%d %d %d %d %d %d %d' % (self.shape, self.centerx, self.centery, self.radiusx, self.radiusy, self.cutband, self.intersection)
    def cutfilter(self, mask):
        top = self.centery - self.radiusy
        left = self.centerx - self.radiusx
        width = 2*self.radiusx
        height = 2*self.radiusy
        if self.cutband:
            tmp = numpy.ones(numpy.shape(mask))
            fg = 0.
        else:
            tmp = numpy.zeros(numpy.shape(mask))
            fg = 1.
            
        for h in range(top, top+height):
            for k in range(left, left+width):
                x = float(k - self.centerx)
                y = float(h - self.centery)
                d = x*x/float(self.radiusx*self.radiusx) + y*y/float(self.radiusy*self.radiusy)
                if d <= 1.:
                    tmp[h, k] = fg
                    
        if self.intersection:
            mask = mask*tmp
        else:
            mask = mask + tmp
            idx = mask > 1.
            mask[idx] = 1.
        return mask
        
class fourier_filter:
    def do_mask(self):
        self.mask = numpy.ones((self.height, self.width))
        for f in self.filters:
            self.mask = f.cutfilter(self.mask)
        
        self.mask = ndimage.uniform_filter(self.mask, 5)
        self.shifted_mask = fftshift(self.mask)
        
    def __init__(self, width, height, str_desc):
        self.width = width
        self.height = height
        self.mask = numpy.ones((height, width))
        self.filters = []
        self.current_spectrum = None
        self.mask_update_cb = None
        lines = str_desc.splitlines()
        lines.pop(0) #remove header
        for l in lines:
            self.filters.append(fourier_filter_descriptor().fromstring(l))
        self.do_mask()    
    
    def edit_mask(self, mask_update_cb=None):
        self.win = pykoala.gui.gui_helpers.MaskTool(filter_cb=self.mask_select, reset_filter_cb=self.clear_mask)
        self.mask_update_cb = mask_update_cb
        plot = self.win.get_plot()
        if self.current_spectrum is None:
            self.mask_image = make.image(data=self.mask, colormap="hot")
        else:
            self.mask_image = make.image(data=self.current_spectrum*self.mask, colormap="hot")
        plot.add_item(self.mask_image)
        self.win.show()
        return self.win
    
    def mask_select(self, shape_descr):
        #print "mask_select" + str(shape_descr)
        self.filters.append(fourier_filter_descriptor().fromlist(shape_descr))
        self._refresh_mask_plot()        
        
    def _refresh_mask_plot(self):
        self.do_mask()
        if self.current_spectrum is None:
            self.mask_image.set_data(self.mask)
        else:
            self.mask_image.set_data(self.current_spectrum*self.mask)
        self.win.get_plot().replot()
        if self.mask_update_cb is not None:
            self.mask_update_cb()
        
    def clear_mask(self):
        self.filters = []
        self.do_mask()
        self._refresh_mask_plot()
        
        
    def __call__(self, input_cplx):
        self.current_spectrum = fftshift(numpy.log(abs(input_cplx) + 1))
        return input_cplx*self.shifted_mask

    def __str__(self):
        tmp = ''
        for f in self.filters:
            tmp += f.__str__()
        return tmp
class spectrum_center:
    def __init__(self, kx, ky, width, height):
        self.kx = kx
        self.ky = ky
        self.width = width
        self.height = height
    
    def set_kx_ky(self, kx, ky):
        self.kx = kx
        self.ky = ky
    
    def __call__(self, input_cplx):
        #kx_s = (self.kx + self.width/2) % self.width
        #ky_s = (self.ky + self.height/2) % self.height
        kx_s = self.kx
        ky_s = self.ky
        res = numpy.zeros(numpy.shape(input_cplx)).astype(complex)
        res[0:self.height-ky_s, 0:self.width-kx_s] = input_cplx[ky_s:self.height, kx_s:self.width]
        res[0:self.height-ky_s, self.width - kx_s:self.width] = input_cplx[ky_s:self.height, 0:kx_s]
        res[self.height-ky_s:self.height, self.width - kx_s:self.width] = input_cplx[0:ky_s, 0:kx_s]
        res[self.height-ky_s:self.height, 0:self.width - kx_s] = input_cplx[0:ky_s, kx_s:self.width]
        return res

class refholo:
    def __init__(self):
        self.rh = None
        self.origin_file = None
        self.ref_holo_path = None
        self.ref_holo_filter_size = None
        
    def gen_rh(self, input_spec, ampl_div=True):
        ift = ifft2(input_spec)
        tmp_ph = numpy.exp(-numpy.angle(ift)*complex(0., 1.))
        if ampl_div:
            tmp_a = 1./numpy.abs(ift)
        else:
            tmp_a = numpy.ones(numpy.shape(input_spec))
    
        self.rh = tmp_a*tmp_ph
        self.origin_file = False        
        self.ref_holo_path = None
        self.ref_holo_filter_size = None        
     
    def save_ref_holo(self, path):
        s = numpy.shape(self.rh)
        binkoala.write_mat_cplx_bin(path, self.rh, s[1], s[0])
    def load_ref_holo(self, path, holoref_filter_size):   
        if path is not None:
            rh = binkoala.read_mat_cplx_bin(path) #open recorded holoref
            self.rh = self.ref_holo_filter(rh, holoref_filter_size)            
            self.origin_file = True
            self.ref_holo_path = path
            self.ref_holo_filter_size = holoref_filter_size
            
    def ref_holo_filter(self, rh, holoref_filter_size):
        #rh=binkoala.read_mat_cplx_bin(path) #open recorded holoref        
        if holoref_filter_size != 1:
            amp = numpy.absolute(rh) #keep original amplitude (not filtered)
            rh = numpy.exp(numpy.angle(rh)*complex(0., 1.)) #amp(rh)=1 to filter
            #%% Apply a polynomial fit of 2nd order to suppress most of the
            #phase jumps before being able to perform the filter on real and img
            #part of rh
            s = numpy.shape(rh)
            height = s[0]
            width = s[1]
            #As only the fit would to be done, don't care about pxsize and lambda
            #The propagation is done with d=0! height and width are the rh size              
            #k=pykoala.core.swl.swl_koala(height,width,pxsize_m,lambda_m,0)
            k = pykoala.core.phaseprocessing
            ph_corr = k.phase_correction(width, height, 0, 0, width, height)
                     
            ##define 2 horizontal and 2 vertical profiles on rh border
            # let 40 pixel from the border of the wavefront
            horiz_seg = k.hv_segment(40, 40, width-80, False) 
            vert_seg = k.hv_segment(40, 40, height-80, True)
            ph_corr.new_hv_seg(horiz_seg)
            ph_corr.new_hv_seg(vert_seg)               
            horiz_seg = k.hv_segment(40, height-80, width-80, False)
            vert_seg = k.hv_segment(width-80, 40, height-80, True)              
            ph_corr.new_hv_seg(horiz_seg)
            ph_corr.new_hv_seg(vert_seg)   
            #perform Fit
            ph_corr.perform_fit(numpy.angle(rh), 2)
            #Extract the computed mask
            mask = ph_corr.phase_mask_img
            #Apply mask to rh
            rh = rh*numpy.exp(mask*complex(0., 1.))
            #%%
            #uniform filter apply separatedly on real and img part or rh
            tmp_real = numpy.real(rh) 
            tmp_img = numpy.imag(rh)
            tmp_real = ndimage.uniform_filter(tmp_real, size=holoref_filter_size)
            tmp_img = ndimage.uniform_filter(tmp_img, size=holoref_filter_size)
            #Compute the filtered complex rh from real and img part
            rh = (tmp_real+tmp_img*complex(0., 1.))
            #Introduce again the phase jumps
            rh = rh*numpy.exp(-mask*complex(0., 1.))
            #again the phase jumps using -mask
            #Create final rh by using the original amplitude and the filtered phase
            rh = amp*numpy.exp(numpy.angle(rh)*complex(0., 1.))

        return rh
            
    def __call__(self, input_spec):
        return fft2(ifft2(input_spec)*self.rh)
        
class fresnel_propagation:
    def __init__(self, width, height, dist_m, pxsize_m, lambda_m):
        self.pxsize_m = pxsize_m
        self.lambda_m = lambda_m
        self.dist_m = dist_m
        self.propagation_mask = numpy.zeros((height, width))
        nx = numpy.arange(0, float(width))/(width - 1) - 0.5
        self.ct = -numpy.pi*lambda_m*dist_m/(pxsize_m*pxsize_m)
        ny = numpy.arange(0, float(height))/(height - 1) - 0.5
        # nx^2 + ny^2
        self.propagation_base = fftshift((numpy.tile(ny*ny, (width, 1)).transpose() +\
            numpy.tile(nx*nx, (height, 1)))*complex(0., 1.))
        self.propagation_mask = numpy.exp(self.ct*self.propagation_base)
    def set_dist(self, dist_m):
        if abs(dist_m - self.dist_m) < 1e-6:
            return
        self.dist_m = dist_m
        self.ct = -numpy.pi*self.lambda_m*self.dist_m/(self.pxsize_m*self.pxsize_m)
        self.propagation_mask = numpy.exp(self.ct*self.propagation_base)
    def __call__(self, input_cplx):
        return ifft2(self.propagation_mask*input_cplx)
        
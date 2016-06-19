# -*- coding: utf-8 -*-
"""
Created on Wed Jan 06 13:33:39 2016

@author: tcolomb
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 14:29:38 2015

@author: tcolomb
"""

#import cv2
import numpy as np
#from matplotlib import pyplot as plt
#from tracking3d.core import binkoala
#from scipy import ndimage
#import os
import sys
#from skimage.restoration import unwrap_phase


#class for different Threshold Method (python implementation of existing methods)
class ThresholdMethod:
    def __init__(self, method, img):#initialize the class with method and image
        #img has to be a 8bit image
        self.method = method
        #0 = Huang
        #1 = Intermodes
        #2 = MaxEntropy
        #3 = ReniEntropy
        #4 = Triangle
        #5 = Otsu (not implemented here use cv2)
        #6 = Yen
        #7 = Moments
        #8 = IsoData
        self.img = img
        self.threshold = -1
    def ComputeThreshold(self):#compute threshold from method
        histogram = np.histogram(self.img, bins=256)#compute histogram of 8bit image
        #choose the method and define the self.threshold
        if self.method == 0:
            self.threshold = self.Huang(histogram[0])
        if self.method == 1:
            self.threshold = self.Intermodes(histogram[0])
        if self.method == 2:
            self.threshold = self.MaxEntropy(histogram[0])
        if self.method == 3:
            self.threshold = self.ReniEntropy(histogram[0])
        if self.method == 4:
            self.threshold = self.Triangle(histogram[0])
        if self.method == 6:
            self.threshold = self.Yen(histogram[0])
        if self.method == 7:
            self.threshold = self.Moments(histogram[0])
        if self.method == 8:
            self.threshold = self.IsoData(histogram[0])
    def Huang(self, data):
        #From http://rsb.info.nih.gov/ij/plugins/download/AutoThresholder.java
        #Implements Huang's fuzzy thresholding method 
        #Uses Shannon's entropy function (one can also use Yager's entropy function) 
        #Huang L.-K. and Wang M.-J.J. (1995) "Image Thresholding by Minimizing  
        #the Measures of Fuzziness" Pattern Recognition, 28(1): 41-51
        #M. Emre Celebi  06.15.2007
        #Ported to ImageJ plugin by G. Landini from E Celebi's fourier_0.8 routines
        threshold = -1
    #Determine the first non-zero bin */
        first_bin = 0
        for ih in range(0, 256):
            if data[ih] != 0:
                first_bin = ih
                break
    #Determine the last non-zero bin */
        last_bin = 255
        for ih in list(reversed(range(first_bin, 256))):
            if data[ih] != 0:
                last_bin = ih
                break
    
        term = 1.0/(last_bin-first_bin)
        mu_0 = np.zeros(256)
        sum_pix = 0
        num_pix = 0
        for ih in range(first_bin, 256):
            sum_pix += ih*data[ih]
            num_pix += data[ih]
            #NUM_PIX cannot be zero ! */
            mu_0[ih] = sum_pix/num_pix
            
        mu_1 = np.zeros(256)
        sum_pix = 0
        num_pix = 0
        for ih in list(reversed(range(0, last_bin+1))):
            sum_pix += ih*data[ih]
            num_pix += data[ih]
            #NUM_PIX cannot be zero ! */
            mu_1[ih-1] = sum_pix/num_pix
    
    #Determine the threshold that minimizes the fuzzy entropy */
        threshold = -1
        min_ent = sys.float_info.max
        for it in range(0, 256):
            ent = 0.0
            for ih in range(0, it+1):
                #Equation (4) in Ref. 1 */
                mu_x = 1.0/(1.0+term*np.abs(ih-mu_0[it]))
                if not ((mu_x < 1e-6) or (mu_x > 0.999999)):
                    #Equation (6) & (8) in Ref. 1 */
                    ent += data[ih]*(-mu_x*np.log(mu_x)-(1.0-mu_x)*np.log(1.0-mu_x))
            for ih in range(it+1, 256): 
                mu_x = 1.0/(1.0+term*np.abs(ih-mu_1[it]))
                if not ((mu_x < 1e-6) or (mu_x > 0.999999)):
                    #Equation (6) & (8) in Ref. 1 */
                    ent += data[ih]*(-mu_x*np.log(mu_x)-(1.0-mu_x)*np.log(1.0-mu_x))
    
            #No need to divide by NUM_ROWS * NUM_COLS * LOG(2) ! */
            if ent < min_ent:
                min_ent = ent
                threshold = it
        return threshold
    
    def bimodalTest(self, y):
        length = np.size(y)
        b = False
        modes = 0
        for k in range(1, length-1):
            if y[k-1] < y[k] and y[k+1] < y[k]:
                modes += 1
                if modes > 2:
                    return False
        if modes == 2:
            b = True
        return b
    def Intermodes(self, data):
#        	// J. M. S. Prewitt and M. L. Mendelsohn, "The analysis of cell images," in
#		// Annals of the New York Academy of Sciences, vol. 128, pp. 1035-1053, 1966.
#		// ported to ImageJ plugin by G.Landini from Antti Niemisto's Matlab code (GPL)
#		// Original Matlab code Copyright (C) 2004 Antti Niemisto
#		// See http://www.cs.tut.fi/~ant/histthresh/ for an excellent slide presentation
#		// and the original Matlab code.
#		//
#		// Assumes a bimodal histogram. The histogram needs is smoothed (using a
#		// running average of size 3, iteratively) until there are only two local maxima.
#		// j and k
#		// Threshold t is (j+k)/2.
#		// Images with histograms having extremely unequal peaks or a broad and
#		// ﬂat valley are unsuitable for this method.
        iHisto = np.zeros(256)
        it = 0
        threshold = -1
        for i in range(0, 256):
            iHisto[i] = data[i]
        tHisto = iHisto
        while not self.bimodalTest(iHisto):
#            //smooth with a 3 point running mean filter
            for i in range(1, 255):
                tHisto[i] = (iHisto[i-1]+ iHisto[i]+iHisto[i+1])/3.
            tHisto[0] = (iHisto[0]+iHisto[1])/3 #//0 outside
            tHisto[255] = (iHisto[254]+iHisto[255])/3 #0 outside
            iHisto = tHisto
            it += 1
            if it > 10000:
                threshold = -1
                print "Intermodes Threshold not found after 10000 iterations."
                return threshold
#The threshold is the mean between the two peaks.
        tt = 0
        for i in range(1, 255):
            if iHisto[i-1] < iHisto[i] and iHisto[i+1] < iHisto[i]:
                tt += i
        threshold = (int)(np.floor(tt/2.))
        return threshold
        
    def MaxEntropy(self, data):
#        // Implements Kapur-Sahoo-Wong (Maximum Entropy) thresholding method
#        // Kapur J.N., Sahoo P.K., and Wong A.K.C. (1985) "A New Method for
#        // Gray-Level Picture Thresholding Using the Entropy of the Histogram"
#        // Graphical Models and Image Processing, 29(3): 273-285
#        // M. Emre Celebi
#        // 06.15.2007
#        // Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines
        threshold = -1
#        int ih, it;
#        int first_bin;
#        int last_bin;
#        tot_ent;  /* total entropy */
#        max_ent;  /* max entropy */
#        ent_back; /* entropy of the background pixels at a given threshold */
#        ent_obj;  /* entropy of the object pixels at a given threshold */
        norm_histo = np.zeros(256) #/* normalized histogram */
        P1 = np.zeros(256) #/* cumulative normalized histogram */
        P2 = np.zeros(256)
        
        total = 0.0
        for ih in range(0, 256):
            total += data[ih]
#        for ih in range(0, 256):
#            norm_histo[ih]=data[ih]/total
        norm_histo = data/total
        P1[0] = norm_histo[0]
        P2[0] = 1.0-P1[0]
        for ih in range(1, 256):
            P1[ih] = P1[ih-1] + norm_histo[ih]
            P2[ih] = 1.0 - P1[ih]
        #/* Determine the first non-zero bin */
        first_bin = 0
        for ih in range(0, 256):
            if not np.abs(P1[ih]) < 2.220446049250313E-16:
                first_bin = ih
                break
        #/* Determine the last non-zero bin */
        last_bin = 255
        for ih in list(reversed(range(first_bin, 256))):
            if not np.abs(P2[ih]) < 2.220446049250313E-16:
                last_bin = ih
                break
#        // Calculate the total entropy each gray-level
#        // and find the threshold that maximizes it 
        max_ent = sys.float_info.min
        for it in range(first_bin, last_bin+1):
            #/* Entropy of the background pixels */
            ent_back = 0.0
            for ih in range(0, it+1):
                if data[ih] != 0:
                    ent_back -= (norm_histo[ih]/P1[it])*np.log(norm_histo[ih]/P1[it])
            #/* Entropy of the object pixels */
            ent_obj = 0.0
            for ih in range(it+1, 256):
                if data[ih] != 0:
                    ent_obj -= (norm_histo[ih]/P2[it])*np.log(norm_histo[ih]/P2[it])
            
                    
            #/* Total entropy */
            tot_ent = ent_back + ent_obj
            #// IJ.log(""+max_ent+"  "+tot_ent);
            if max_ent < tot_ent:
                max_ent = tot_ent
                threshold = it
                
        return threshold
    
    def partialSum(self, y, j):
        x = 0
        for i in range(0, j+1):
            x += y[i]
        return x
    def ReniEntropy(self, data):
#        // Kapur J.N., Sahoo P.K., and Wong A.K.C. (1985) "A New Method for
#        // Gray-Level Picture Thresholding Using the Entropy of the Histogram"
#        // Graphical Models and Image Processing, 29(3): 273-285
#        // M. Emre Celebi
#        // 06.15.2007
#        // Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines
#        int threshold; 
#        int opt_threshold;
#        int ih, it;
#        int first_bin;
#        int last_bin;
#        int tmp_var;
#        int t_star1, t_star2, t_star3;
#        int beta1, beta2, beta3;
#        double alpha;/* alpha parameter of the method */
#        double term;
#        double tot_ent;  /* total entropy */
#        double max_ent;  /* max entropy */
#        double ent_back; /* entropy of the background pixels at a given threshold */
#        double ent_obj;  /* entropy of the object pixels at a given threshold */
#        double omega;
        norm_histo = np.zeros(256) #/* normalized histogram */
        P1 = np.zeros(256) #/* cumulative normalized histogram */
        P2 = np.zeros(256)
        
        total = 0.0
        for ih in range(0, 256):
            total += data[ih]
        norm_histo = data/total
        P1[0] = norm_histo[0]
        P2[0] = 1.0-P1[0]
        for ih in range(1, 256):
            P1[ih] = P1[ih-1] + norm_histo[ih]
            P2[ih] = 1.0 - P1[ih]
        
        first_bin = 0
        for ih in range(0, 256):
            if not np.abs(P1[ih]) < 2.220446049250313E-16:
                first_bin = ih
                break
        #/* Determine the last non-zero bin */
        last_bin = 255
        for ih in list(reversed(range(first_bin, 256))):
            if not np.abs(P2[ih]) < 2.220446049250313E-16:
                last_bin = ih
                break
#        /* Maximum Entropy Thresholding - BEGIN */
#        /* ALPHA = 1.0 */
#        /* Calculate the total entropy each gray-level
#        and find the threshold that maximizes it 
#        */
        threshold = 0# was MIN_INT in original code, but if an empty image is processed it gives 
            #an error later on.
        max_ent = 0.0
        for it in range(first_bin, last_bin+1):
            #/* Entropy of the background pixels */
            ent_back = 0.0
            for ih in range(0, it+1):
                if data[ih] != 0:
                    ent_back -= (norm_histo[ih]/P1[it])*np.log(norm_histo[ih]/P1[it])
            #/* Entropy of the object pixels */
            ent_obj = 0.0
            for ih in range(it+1, 256):
                if data[ih] != 0:
                    ent_obj -= (norm_histo[ih]/P2[it])*np.log(norm_histo[ih]/P2[it])
            #/* Total entropy */
            tot_ent = ent_back + ent_obj
            #// IJ.log(""+max_ent+"  "+tot_ent);
            if max_ent < tot_ent:
                max_ent = tot_ent
                threshold = it
        t_star2 = threshold
        #/* Maximum Entropy Thresholding - END */
        threshold = 0#//was MIN_INT in original code, but if an empty image is processed it gives
            #an error later on.
        max_ent = 0.0
        alpha = 0.5
        term = 1.0/(1.0 - alpha)
        for it in range(first_bin, last_bin+1):
            #/* Entropy of the background pixels */
            ent_back = 0.0
            for ih in range(0, it+1):
                ent_back += np.sqrt(norm_histo[ih]/P1[it])
            #/* Entropy of the object pixels */
            ent_obj = 0.0
            for ih in range(it+1, 256):
                ent_obj += np.sqrt(norm_histo[ih]/P2[it])
            #/* Total entropy */
            if ent_back*ent_obj > 0.0:
                tot_ent = term*np.log(ent_back*ent_obj)
            else:
                tot_ent = 0
#            tot_ent = term*((ent_back*ent_obj)>0.0 ? Math.log (ent_back*ent_obj): 0.0);
            if tot_ent > max_ent:
                max_ent = tot_ent
                threshold = it
        t_star1 = threshold
        threshold = 0 #//was MIN_INT in original code, but if an empty image is processed it gives
            #an error later on.
        max_ent = 0.0
        alpha = 2.0
        term = 1.0/(1.0-alpha)
        for it in range(first_bin, last_bin+1):
            #/* Entropy of the background pixels */
            ent_back = 0.0
            for ih in range(0, it+1):
                ent_back += (norm_histo[ih]*norm_histo[ih])/(P1[it]*P1[it])
            #/* Entropy of the object pixels */
            ent_obj = 0.0
            for ih in range(it+1, 256):
                ent_obj += (norm_histo[ih]*norm_histo[ih])/(P2[it]*P2[it])
            #/* Total entropy */
            if ent_back*ent_obj > 0.0:
                tot_ent = term*np.log(ent_back*ent_obj)
            else:
                tot_ent = 0
#            tot_ent = term *( ( ent_back * ent_obj ) > 0.0 ? Math.log(ent_back * ent_obj ): 0.0 );
            if tot_ent > max_ent:
                max_ent = tot_ent
                threshold = it
        t_star3 = threshold
        #/* Sort t_star values */
        if t_star2 < t_star1:
            tmp_var = t_star1
            t_star1 = t_star2
            t_star2 = tmp_var
        if t_star3 < t_star2:
            tmp_var = t_star2
            t_star2 = t_star3
            t_star3 = tmp_var
        if t_star2 < t_star1:
            tmp_var = t_star1
            t_star1 = t_star2
            t_star2 = tmp_var
        #/* Adjust beta values */
        if np.abs(t_star1-t_star2) <= 5:
            if np.abs(t_star2-t_star3) <= 5:
                beta1 = 1
                beta2 = 2
                beta3 = 1
            else:
                beta1 = 0
                beta2 = 1
                beta3 = 3
        else:
            if np.abs(t_star2-t_star3) <= 5:
                beta1 = 3
                beta2 = 1
                beta3 = 0
            else:
                beta1 = 1
                beta2 = 2
                beta3 = 1
#        //IJ.log(""+t_star1+" "+t_star2+" "+t_star3);
#        /* Determine the optimal threshold value */
        omega = P1[t_star3] - P1[t_star1]
        opt_threshold = (int)(t_star1*(P1[t_star1]+0.25*omega*beta1)+\
        0.25*t_star2*omega*beta2+t_star3*(P2[t_star3]+0.25*omega*beta3))
        return opt_threshold
    
    def Triangle(self, data):
#        //  Zack, G. W., Rogers, W. E. and Latt, S. A., 1977,
#        //  Automatic Measurement of Sister Chromatid Exchange Frequency,
#        // Journal of Histochemistry and Cytochemistry 25 (7), pp. 741-753
#        //
#        //  modified from Johannes Schindelin plugin
#        // 
#        // find min and max
        minimum = 0
        dmax = 0
        maximum = 0
        min2 = 0
        #for i in range(0,np.size(data)):
        for i in range(0, 256):
            if data[i] > 0:
                minimum = i
                break
        if minimum > 0:
            minimum -= 1#// line to the (p==0) point, not to data[min]
#        // The Triangle algorithm cannot tell whether the data is skewed to one side or another.
#        // This causes a problem as there are 2 possible thresholds between the max and the 2 
            #extremes of the histogram.
#        // Here I propose to find out to which side of the max point the data is furthest, and 
            #use that as the other extreme.
        for i in list(reversed(range(1, 256))):
            if data[i] > 0:
                min2 = i
                break
        if min2 < 255:
            min2 += 1 #// line to the (p==0) point, not to data[min]
        for i in range(0, 256):
            if data[i] > dmax:
                maximum = i
                dmax = data[i]
#        // find which is the furthest side
#        //IJ.log(""+min+" "+max+" "+min2);
        inverted = False
        if (maximum-minimum) < (min2-maximum):
#            // reverse the histogram
#            //IJ.log("Reversing histogram.");
            inverted = True
            left = 0         #// index of leftmost element
            right = 255 #// index of rightmost element
            while left < right:
                #// exchange the left and right elements
                temp = data[left] 
                data[left] = data[right] 
                data[right] = temp
                #// move the bounds toward the center
                left += 1
                right -= 1
            minimum = 255-min2
            maximum = 255-maximum

        if minimum == maximum:
            #//IJ.log("Triangle:  min == max.");
            return minimum
        #// describe line by nx * x + ny * y - d = 0
        #double nx, ny, d;
        #// nx is just the max frequency as the other point has freq=0
        nx = np.float32(data[maximum]) #  -min; data[min]; lowest value bmin = (p=0)% in the image
        ny = minimum - maximum
        d = np.sqrt(nx*nx + ny*ny)
        nx /= d
        ny /= d
        d = nx * minimum + ny * data[minimum]
        #// find split point
        split = minimum

        splitDistance = 0.0
        for i in range(minimum+1, maximum+1):
            newDistance = nx*i+ny*data[i]-d
            if newDistance > splitDistance:
                split = i
                splitDistance = newDistance
        split -= 1 
        if inverted:
            #// The histogram might be used for something else, so let's reverse it back
            left = 0
            right = 255
            while left < right:
                temp = data[left]
                data[left] = data[right]
                data[right] = temp
                left += 1
                right -= 1
            return 255-split
        else:
            return split
            
    def Yen(self, data):
#        // Implements Yen  thresholding method
#        // 1) Yen J.C., Chang F.J., and Chang S. (1995) "A New Criterion 
#        //    for Automatic Multilevel Thresholding" IEEE Trans. on Image 
#        //    Processing, 4(3): 370-378
#        // 2) Sezgin M. and Sankur B. (2004) "Survey over Image Thresholding 
#        //    Techniques and Quantitative Performance Evaluation" Journal of 
#        //    Electronic Imaging, 13(1): 146-165
#        //    http://citeseer.ist.psu.edu/sezgin04survey.html
#        //
#        // M. Emre Celebi
#        // 06.15.2007
#        // Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines
#        int threshold;
#        int ih, it;
#        double crit;
#        double max_crit;
        norm_histo = np.zeros(256) #/* normalized histogram */
        P1 = np.zeros(256) # /* cumulative normalized histogram */
        P1_sq = np.zeros(256) # 
        P2_sq = np.zeros(256) #

        total = 0.0
        for ih in range(0, 256):
            total += data[ih]
        norm_histo = data/total
        P1[0] = norm_histo[0]
        for ih in range(1, 256):
            P1[ih] = P1[ih-1] + norm_histo[ih]

        P1_sq[0] = norm_histo[0]*norm_histo[0]
        for ih in range(1, 256):
            P1_sq[ih] = P1_sq[ih-1] + norm_histo[ih] * norm_histo[ih]

        P2_sq[255] = 0.0
        for ih in list(reversed(range(0, 255))):
            P2_sq[ih] = P2_sq[ih + 1] + norm_histo[ih + 1] * norm_histo[ih + 1]
        #/* Find the threshold that maximizes the criterion */
        threshold = -1
        max_crit = sys.float_info.min
        for it in range(0, 256):
            temp1 = 0.0
            temp2 = 0.0
            if P1_sq[it] * P2_sq[it] > 0.0:
                temp1 = -1.0*np.log(P1_sq[it]*P2_sq[it])
            if P1[it]*(1.0-P1[it]) > 0.0:
                temp2 = 2*np.log(P1[it]*(1.0-P1[it]))
            crit = temp1+temp2
            if crit > max_crit:
                max_crit = crit
                threshold = it

        return threshold
        
    def Moments(self, data):
#W. Tsai, "Moment-preserving thresholding: a new approach," Computer Vision,
#Graphics, and Image Processing, vol. 29, pp. 377-393, 1985.
#Ported to ImageJ plugin by G.Landini from the the open source project FOURIER 0.8
#by  M. Emre Celebi , Department of Computer Science,  Louisiana State University in Shreveport
#Shreveport, LA 71115, USA
#http://sourceforge.net/projects/fourier-ipal
# http://www.lsus.edu/faculty/~ecelebi/fourier.htm
        total = 0.0
        m0 = 1.0
        m1 = 0.0
        m2 = 0.0
        m3 = 0.0
        summ = 0.0
        p0 = 0.0
#        double cd, c0, c1, z0, z1;	/* auxiliary variables */
        threshold = -1
        
        histo = np.ones(256)
        #total=np.sum(data)
        for i in range(0, 256):
            total += data[i]
        histo = data/total    #//normalised histogram
        
        #/* Calculate the first, second, and third order moments */
        for i in range(0, 256):
            m1 += i * histo[i]
            m2 += i * i * histo[i]
            m3 += i * i * i * histo[i]
#        First 4 moments of the gray-level image should match the first 4 moments
#        of the target binary image. This leads to 4 equalities whose solutions
#        are given in the Appendix of Ref. 1 
        cd = m0*m2-m1*m1
        c0 = (-m2*m2+m1*m3)/cd
        c1 = (m0*-m3+m2*m1)/cd
        z0 = 0.5*(-c1-np.sqrt(c1*c1-4.0*c0))
        z1 = 0.5*(-c1+np.sqrt(c1*c1-4.0*c0))
        p0 = (z1-m1)/(z1-z0)  #/* Fraction of the object pixels in the target binary image */
        
#        // The threshold is the gray-level closest  
#        // to the p0-tile of the normalized histogram 
        summ = 0
        for i in range(0, 256):
            summ += histo[i]
            if summ > p0:
                threshold = i
                break

        return threshold
        
    def IsoData(self, data):
#Also called intermeans
#Iterative procedure based on the isodata algorithm [T.W. Ridler, S. Calvard, Picture 
#thresholding using an iterative selection method, IEEE Trans. System, Man and 
#Cybernetics, SMC-8 (1978) 630-632.] 
#The procedure divides the image into objects and background by taking an initial threshold,
#then the averages of the pixels at or below the threshold and pixels above are computed. 
#The averages of those two values are computed, the threshold is incremented and the 
#process is repeated until the threshold is larger than the composite average. That is,
# threshold = (average background + average objects)/2
#The code in ImageJ that implements this function is the getAutoThreshold() method in the
#ImageProcessor class. 
#        //
#From: Tim Morris (dtm@ap.co.umist.ac.uk)
#Subject: Re: Thresholding method?
#posted to sci.image.processing on 1996/06/24
#The algorithm implemented in NIH Image sets the threshold as that grey
#value, G, for which the average of the averages of the grey values
#below and above G is equal to G. It does this by initialising G to the
#lowest sensible value and iterating:
#            
#L = the average grey value of pixels with intensities < G
#H = the average grey value of pixels with intensities > G
#is G = (L + H)/2?
#yes => exit
#no => increment G and repeat
#        //
#However the code commented below have the occasional discrepancy with IJ code.
#so I am reproducing the original IJ implementation for compatibility.
#        int level;
        maxValue = np.size(data) - 1
#        double result, sum1, sum2, sum3, sum4
        
        minimum = 0
        while data[minimum] == 0 and minimum < maxValue:
            minimum += 1

        maximum = maxValue
        while data[maximum] == 0 and maximum > 0:
            maximum -= 1
        if minimum >= maximum:
            level = np.size(data)/2
            return level
        
        movingIndex = minimum
#        inc = max(maximum/40, 1)
        condition = True
        while condition:
            sum1 = sum2 = sum3 = sum4 = 0.0
            for i in range(minimum, movingIndex+1):
                sum1 += i*data[i]
                sum2 += data[i]
            for i in range(movingIndex+1, maximum+1):
                sum3 += i*data[i]
                sum4 += data[i]
            result = (sum1/sum2 + sum3/sum4)/2.0
            movingIndex += 1
            condition = ((movingIndex+1) <= result and movingIndex < maximum-1)
        
#        //.showProgress(1.0);
        level = (int)(np.round(result))
        return level
#/*
#		// this should work straight away, but for some images there is a discrepancy with IJ
#		int i, l, toth, totl, h, g=0;
#		for (i = 1; i < 256; i++){
#			if (data[i] > 0){
#				g = i + 1;
#				break;
#			}
#		}
#		while (true){
#			l = 0;
#			totl = 0;
#			for (i = 0; i < g; i++) {
#				 totl = totl + data[i];
#				 l = l + (data[i] * i);
#			}
#			h = 0;
#			toth = 0;
#			for (i = g + 1; i < 256; i++){
#				toth += data[i];
#				h += (data[i]*i);
#			}
#			if (totl > 0 && toth > 0){
#				l /= totl;
#				h /= toth;
#				if (g == (int) Math.round((l + h) / 2.0))
#					break;
#			}
#			g++;
#			if (g > 254)
#				return -1;
#		}
#		return g+1; // +1 to match ImageJ output, but some images are off by +1 or -1   

        
            
        
          
        
        

# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 17:51:17 2014

Some handy tools based on numpy and ROOT package, wrapped into handy form.
Include:
    1. Determing background : background
    2. Smoothing.
    3. Peak searching.
@author: eastwoodknight
"""

import ROOT as rt
import numpy as np
from scipy import signal
#from read_files import get_front_spectrs



def GetHist(hist,bins):
    """get data from TH1F root hist in a form of np.array"""
    new_hist = []
    for i in range(bins):
        a = hist.GetBinContent(i)
        new_hist.append(a)
    return np.array(new_hist)
  
  
  
def SetHist(hist,np_hist):
    """set data from TH1F root hist in a form of np.array"""
    for i in range(len(np_hist)):
        hist.SetBinContent(i,np_hist[i])



def background(y,niter=30,parameters=''):
    """
    return the background of the given spectrum y;
    look into ROOT's TSpectrum manual to configure parameters"""
    lenY = len(y)
    hist = rt.TH1F('hist','Spectrum',lenY,1,lenY+1)
    SetHist(hist,y)
    
    s = rt.TSpectrum()
    sm_hist = s.Background(hist,niter,parameters)
    return GetHist(sm_hist,lenY)



def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
        
     This method is based on the convolution of a scaled window with the signal.
     The signal is prepared by introducing reflected copies of the signal 
     (with the window size) in both ends so that transient parts are minimized
     in the begining and end part of the output signal.
      
      input:
           x: the input signal 
           window_len: the dimension of the smoothing window; should be an odd integer
           window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
               flat window will produce a moving average smoothing.
   
      output:
           1D np.array - the smoothed signal
           
      example:
   
       t=linspace(-2,2,0.1)
       x=sin(t)+randn(len(t))*0.1
       y=smooth(x)
       
     see also: 
       
     numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
     scipy.signal.lfilter  
    """
    if x.ndim != 1:
         raise ValueError("spectrum.tools smooth: Only 1-dimension arrays are valid.")
    if x.size < window_len:
         raise ValueError("spectrum.tools smooth: Input vector should be bigger than window size.")
    if window_len<3:
         raise ValueError( "spectrum.tools smooth: Window length is too small .")
 
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
         raise ValueError( "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
   
   
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': #moving average
         w=np.ones(window_len,'d')
    else:
         w=eval('np.'+window+'(window_len)') 
    y=np.convolve(w/w.sum(),s,mode='valid')
    y = y[(window_len/2) :-(window_len/2)]# -1]
    
    return y
  
  
def search_peaks(data_x,data_y,parameters): #bartlett hanning
    """
    Based on scipy.signal.find_peaks_cwt peak-finding tool.
    threshold parameter allow to cut spectrum on the given level.
    Return 2 arrays, containing x and y coordinated of found peaks. 'ricker'
    """
    data_x, data_y = np.array(data_x), np.array(data_y)
    f_wavelet = eval('signal.'+parameters.wavelet)
    peak_id = signal.find_peaks_cwt(data_y,widths=parameters.widths,wavelet=f_wavelet,min_length=parameters.min_length,\
                                    min_snr=parameters.min_snr, noise_perc=parameters.noise_perc)
    if len(peak_id)==0:
        raise Exception('No peak was found.')
    peak_id = np.array(peak_id) - 1
    xpeaks,ypeaks = data_x[peak_id], data_y[peak_id]
    numbers = ypeaks>(ypeaks.max()*parameters.noise_perc)
    ypeaks = ypeaks[numbers]
    xpeaks = xpeaks[numbers]#-1
    return xpeaks, ypeaks



class Search_peak_error(Exception): pass
    
def search_peaks1(data_x,data_y,parameters): #bartlett hanning
    """
    Based on scipy.signal.find_peaks_cwt peak-finding tool.
    threshold parameter allow to cut spectrum on the given level.
    Return 2 arrays, containing x and y coordinated of found peaks. 'ricker'
    """
    data_x, data_y = np.array(data_x), np.array(data_y)
    f_wavelet = eval('signal.'+parameters.wavelet)
    peak_id = signal.find_peaks_cwt(data_y,widths=parameters.widths,wavelet=f_wavelet,min_length=parameters.min_length,\
                                    min_snr=parameters.min_snr, noise_perc=parameters.noise_perc)
    if len(peak_id)==0:
        raise Search_peak_error('No peak was found')
    peak_id = np.array(peak_id) - 1
    xpeaks,ypeaks = data_x[peak_id], data_y[peak_id]
    numbers = ypeaks>(ypeaks.max()*parameters.noise_perc)
    ypeaks = ypeaks[numbers]
    xpeaks = xpeaks[numbers]#-1
    return xpeaks, ypeaks



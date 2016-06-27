# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 15:08:07 2014

Code for processing automatic line-calibration of alpha-scale in our experiments.
Based on spectrum_tools.py 
 
@author: eastwoodknight
"""

from Pyspectr.spectrum_tools import *
#from read_files import * #get_front_spectrs
from scipy.optimize import leastsq
import numpy as np

#Alpha_energies = [6040,6143.8,6264,6899.2,7137,7922,8699,9261]

class Calibration_Error(Exception): pass 
 
def show_spectrs(xsample0,sample0,xsample,xmin,xmax,xpeaks,ypeaks,sample,energies,title=None,solution=None):

    if title:
        fig.suptitle(title)
    ax[0].set_title('Raw spectrum')
    ax[0].plot(xsample0,sample0,linestyle='steps') 
    ax[0].set_ylim([0,max(ypeaks)*2.15])
    
    ax[1].set_title('Processed spectrum with marked calibration peaks')
    ax[1].set_xlim([xmin,xmax])					
    ax[1].set_ylim([0,max(ypeaks)*1.25])
    ax[1].plot(xsample,sample,linestyle='steps',color='b',linewidth=3) 
    ax[1].plot(xpeaks,ypeaks,'ro',linewidth=4)
    
    try:
        ax[2].set_xlim([xmin,xmax])
        ax[2].set_title('Calibration function')
        ax[2].plot(xpeaks,energies,'ro',linewidth=4,label='Data points')
        ax[2].plot(xpeaks,solution[0]*np.ones(len(xpeaks))+solution[1]*xpeaks,linewidth=2,label='Fitting line')
    except:
        pass
    #print solution[0]*np.ones(len(xpeaks))+solution[1]*xpeaks
    ax[2].legend(loc='lower right')

    
    
    
def filter_spectr(xmin,xmax,xsample,sample,properties): #window_smooth=7,smooth_wavelet='hanning',background_options='BACK1_ORDER4'):
#    sample = np.array(sample)
    x_ind = (xsample >= xmin) & (xsample <= xmax)
    sample = sample[x_ind]
    xsample = xsample[x_ind] #).tolist() 
    if (sample == 0).sum() == len(sample):
        raise Calibration_Error('No data in given sample')
    if properties.window_smooth: 
        sample = smooth(sample,properties.window_smooth,properties.smooth_wavelet)
        sample_background = background(sample,parameters=properties.background_options) #,BACK1_INCLUDE_COMPTON
        sample -= sample_background
        sample[ sample< properties.threshold ] = 0  
    return xsample,sample
  


def calibrate_area(xsample,sample,xmin,xmax,calibration_properties,filter_properties,search_properties): #threshold=0.25,sigma=3,visualize=True,energies=Alpha_energies):  
    """ 
    Example of input parameter's classes:

    class record: pass
    calibration_properties = record()
    calibration_properties.visualize= False
    calibration_properties.sigma = 15 # параметр для уточнения середины пика методом среднего взвешенного, должен соотв. примерно полуширине пика, может быть None 
    calibration_properties.dlt = 32 # параметр для определения области соответствия расчетного положения пиков (соотв. пропорции калибровочных энергий) и реальных найденных пиков (метод spectrum_tools.search_peaks), если реальные пики не находятся вблизи расчетных 
                                    # в области dlt, то они заменяются на расчетные, в противном случае выбирается пик наиболее близкий к расчетному. see calibration.calibrate_area
    calibration_properties.energies = [6040,6143,6264,6899.2,7137,7922,8699,9261]#[7137,7922,8699,9261]#[6040,6143,6264,8699,9261]# 

    search_properties = record() #for noisy spectrum
    search_properties.widths=np.arange(1,5)
    search_properties.wavelet= 'ricker'
    search_properties.min_length=1.
    search_properties.min_snr=0.9
    search_properties.noise_perc=0.3 
    
    filter_properties = record()
    filter_properties.window_smooth=7
    filter_properties.smooth_wavelet='blackman'
    filter_properties.background_options='BACK1_ORDER8,BACK1_INCLUDE_COMPTON' 
    filter_properties.threshold= 3    
    
    Return xpeaks, solution (result of linear fitting of spectrum)
    """
#    sample, xsample = np.array(sample), np.array(xsample)
    xsample0 = xsample #make copies to visualize it later
    sample0 = sample
    
    #initiatian of properties
    try:
        weighted_average_sigma = calibration_properties.sigma
        dlt = calibration_properties.dlt
        visualize = calibration_properties.visualize
        energies = calibration_properties.energies
    except:
        raise Exception('Wrong calibration_properties object')
        
    try:
        calibration_properties.title
    except:
        calibration_properties.title = None
    
    #filter data - smooth and delete a background
    if filter_properties:
        xsample,sample = filter_spectr(xmin,xmax,xsample,sample,filter_properties)#5,smooth_wavelet='hanning',background_options='BACK1_ORDER2')
    
    #find peaks
    try:
        xpeaks,ypeaks = search_peaks1(xsample,sample,search_properties)#sigma=sigma,threshold=threshold) 
    except Search_peak_error as e:
        raise Calibration_Error( e)
    xpeaks = np.array(xpeaks)
    print(xpeaks)
    
    #delete too close peaks, select the highest of close peaks and sort them
    indx = np.concatenate(([True],np.abs(np.diff(xpeaks))>dlt)) 
    #print indx
    msv = []
    k = -1
    for i,j in enumerate(indx):
        if j > 0: 
            if k > -1:
                s = k+np.argmax(ypeaks[k:i]) if i-k>1 else k
                msv.append(s)
            k = i
    if k == len(xpeaks)-1:
        msv.append(k)
    else:
        msv.append(k+np.argmax(ypeaks[k:]))
    
    indx = np.array(msv)  
#    print 'x,y:'
#    for i,j in enumerate(xpeaks):
#        print i,j,ypeaks[i]
#    print 'indx:',indx     
#    print xpeaks[indx]
#    print ypeaks[indx]     
#    print type(xpeaks), 'indx:',indx, 'dlt: ',dlt,'diff: ',np.diff(xpeaks)
    xpeaks,ypeaks= xpeaks[indx],ypeaks[indx]
    
    #check if it's enough peaks detected
    if len(xpeaks) < 4:
        print( 'Error occured! Peaks found:',xpeaks,ypeaks)
        show_spectrs(xsample0,sample0,xsample,xmin,xmax,xpeaks,ypeaks,sample,energies,solution=None,title='Error report')
        raise Calibration_Error("Some peaks weren't indentificated in calibration procedure, you should change the parameters.")
        
    k = 15
    if len(xpeaks) > k:
        indx = ypeaks.argsort()
        ypeaks = ypeaks[indx][-k:]
        xpeaks = xpeaks[indx][-k:]
        
    ypeaks = ypeaks[xpeaks.argsort()] #sort lists of peaks
    xpeaks = sorted(xpeaks)
    xpeaks,ypeaks = np.array(xpeaks),np.array(ypeaks)
#    print 'afterproc',xpeaks

    #selecting valid peaks
    spectr_peak_dists = np.diff(np.array(energies,dtype = np.float64) ) / (energies[-1]-energies[0])
    spectr_length = (xpeaks[-1] - xpeaks[-2])/spectr_peak_dists[-1] 
    spectr_peak_dists1 = spectr_peak_dists*spectr_length # distances from right edge to points, proportion between them is the same as for energy calibration spectr
    spectr_peak_dists1 = xpeaks[-1] - spectr_peak_dists1[::-1].cumsum()
    x,y = [],[]
    x.append(xpeaks[-1])
    y.append(ypeaks[-1])
        
    def find_closest(msv,m0,dlt): #find a point from msv array, which is the closest to m0 inside dlt diapason
        msv=abs(msv-m0)
        i = msv.argmin()
        if msv[i] < dlt:
            return i
        else:
            return False
                
    for i in range(len(spectr_peak_dists1)):
        l =spectr_peak_dists1[i]
#        print 'l calculated peak:',l
        peak_ind = find_closest(xpeaks,l,dlt)
#        print 'pass:',peak_ind
#        if ypeaks[peak_ind] < 0.25*y[-1]:
#            continue
        if not peak_ind or xpeaks[peak_ind] == x[-1]: 
            ind = (xsample > l).argmax()
            x.append(xsample[ind])
            y.append(sample[ind])
            continue
        x.append(xpeaks[peak_ind])  
        y.append(ypeaks[peak_ind]) 
        spectr_length += l - x[-1] 
        spectr_peak_dists1 = spectr_peak_dists*spectr_length
        spectr_peak_dists1 = xpeaks[-1] - spectr_peak_dists1[::-1].cumsum()
               
    x.reverse()
    y.reverse()
    xpeaks,ypeaks= np.array(x,dtype=np.float64),np.array(y,dtype=np.float64) 
#    print 'subresult',xpeaks
    if len(xpeaks) < len(energies):
        print( 'Not enough peaks')
        print( 'Peaks found: ',xpeaks,ypeaks)

        raise Calibration_Error('Not enough valid peaks' )
         
    #correct xpeaks by calculating a weighted average of x using +-2*sigma window   
    if weighted_average_sigma:
        k = weighted_average_sigma
        #print 'start wa',k
        xpeaks1 = []
        for x in xpeaks:
            #print xsample, x in xsample	
            ind = (xsample >= x).argmax()#xsample.index(x)  
            hist_sum = sample[ind-k:ind+k+1]
            ind1 = hist_sum > 0
            hist_sum = hist_sum.sum()
            sample1 = sample[ind-k:ind+k+1]
            xsample1 = xsample[ind-k:ind+k+1]
            sample1 = sample1[ind1]
            xsample1 = xsample1[ind1]
            xpeaks1.append( (sample1*xsample1).sum()/hist_sum )
        xpeaks = np.array(xpeaks1) 
        
    #fitting by line
    def residuals(coef,y,x):
        return y - coef[0]*np.ones(len(x)) - coef[1]*x
    p0 = (-20,16) #init coefficients    
    
    xpeaks,ypeaks= np.array(xpeaks,dtype=np.float64),np.array(ypeaks,dtype=np.float64)
    solution = leastsq(residuals,p0, args = (energies,xpeaks) )[0]
    
    #output
    if visualize:
        show_spectrs(xsample0,sample0,xsample,xmin,xmax,xpeaks,ypeaks,sample,energies,solution=solution,title=calibration_properties.title)
    #print 'pr',xpeaks
    return xpeaks,solution    



def make_report(xpeaks,solution,ind,energies,filename=None):
    #energies = np.array([6040,6143,6264,6899,7137,7922,8699,9265]) 
    report = '\n%d    A = %2.5f ; B = %2.1f\n'%(ind,solution[1],solution[0])
    report += '%5s  %5s      %4s   %5s \n'%('Eexp','Ecal','Differ','Channel')
    S = 0
    for i,en in enumerate(energies):
        Ecalc = solution[0]+solution[1]*xpeaks[i]
        report += '%4.1f  %4.1f    %-5.1f    %-4.1f \n'%(en,Ecalc,Ecalc-en,xpeaks[i]) 
        S += (Ecalc-en)**2
    
    if filename:
        f = open(filename,'a')
        f.write(report)
        f.close()
    report += 'S/n = %3.1f \n' % (S**0.5/len(energies))
    return report, (S**0.5/len(energies))     
        


def calibrate_spectr(start,stop,xmin,xmax,hist,calibration_properties,filter_properties,search_properties,output_filename=None):
    good_results = 0  
    bad_results = []
    xpeaks_list = []
    coef_list = []
    for i in range(start,stop+1):
        print( 'strip: ',i)
        try:
            xpeaks,coef=calibrate_area(hist.index,hist[i],xmin,xmax,calibration_properties,filter_properties,search_properties)
            report, control = make_report(xpeaks,coef,i,energies=calibration_properties.energies,filename = output_filename)
            print( report)
            if control <= 10:
                good_results += 1
                xpeaks_list.append(xpeaks)
                coef_list.append(coef)
            else:
                bad_results.append(i)
        except Calibration_Error as e:
            bad_results.append(i)
            print( 'strip %1d was not calibrated correctly: %s' %(i,e))
            print( 'Error occured',e)
        except KeyError as e:
            print( 'Error occured: in hist array \n',e)
            bad_results.append(i)
            
    print( 'Good results:', good_results,'/',abs(stop - start)+1)
    print( 'Bad results:', bad_results)       
    return(xpeaks_list,coef_list)
        


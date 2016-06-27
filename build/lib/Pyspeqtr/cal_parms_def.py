def set_default():
    """
        Set default Calibration parameters in a way that makes it easy to redefine them.
        var is a variable name or list of variable names to set/reset. 
        # dlt dicritises the spectra to find peaks in each section. e.c
        # option to determine the conformity of the calculated position of the peaks 
        #(resp. calibration proportion of energy), and found the real peaks (method spectrum_tools.search_peaks)
        #, if the real peaks are close to settlement
        #параметр для определения области соответствия расчетного положения пиков (соотв. пропорции калибровочных энергий) 
        #и реальных найденных пиков (метод spectrum_tools.search_peaks), если реальные пики не находятся вблизи расчетных 
        #in dlt, they are replaced by estimates, otherwise selects the peak closest to the calculated. see calibration.calibrate_area.
        # в области dlt, то они заменяются на расчетные, в противном случае выбирается пик наиболее близкий к расчетному. 
    """   
    import numpy as np 
    #Define Parameter Classes
    class record: pass

    #energy scale check-calibration
    calibration_properties = record()
    calibration_properties.visualize=True
    calibration_properties.sigma = 10      #Assumed width of each peak. Count the channels! 
    calibration_properties.dlt = 20
    calibration_properties.energies = [466, 523.1,569.5, 650.3,665.4, 759.7, 848.8]
#[6040,6143,6264,6899.2,7137,7922,8699,9261]   #Other Options [7137,7922,8699,9261]#[6040,6143,6264,8699,9261]#   
    search_properties = record()                            #for smooth spectrum
    search_properties.widths=np.arange(1,10)
    search_properties.wavelet= 'ricker'
    search_properties.min_length=1.
    search_properties.min_snr=0.7
    search_properties.noise_perc= 0.1
  
    filter_properties = record()
    filter_properties.window_smooth=7
    filter_properties.smooth_wavelet='blackman'
    filter_properties.background_options='BACK1_ORDER8,BACK1_INCLUDE_COMPTON' 
    filter_properties.threshold = 4
 
    return(calibration_properties,filter_properties,search_properties)


def set_default_Yb():
    import numpy as np 
    #Define Parameter Classes 
    class record: pass
    #energy scale check-calibration
    calibration_properties = record()
    calibration_properties.visualize=True
    calibration_properties.sigma=6
    calibration_properties.dlt=4
    calibration_properties.energies=[6040, 6143, 6264, 6899.2, 7137, 7922, 8699, 9261]

    search_properties = record()               #for smooth spectrum
    search_properties.widths=np.arange(1,10)
    search_properties.wavelet='ricker'
    search_properties.min_length=1.0
    search_properties.min_snr=0.7
    search_properties.noise_perc=0.1

    filter_properties = record()
    filter_properties.window_smooth=5
    filter_properties.smooth_wavelet='blackman'
    filter_properties.background_options='BACK1_ORDER8,BACK1_INCLUDE_COMPTON'
    filter_properties.threshold = 6
    return(calibration_properties,filter_properties,search_properties)


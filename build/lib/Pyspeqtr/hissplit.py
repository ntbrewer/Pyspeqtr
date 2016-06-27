def his_split(his, gate, cycle, pars,t_bin=1):
        """
          Documentation strings go here.

        """
        import copy

        T0 = {'name' : 'T0', 'value' : cycle[0], 'vary' : False}
        T1 = {'name' : 'T1', 'value' : cycle[1], 'vary' : False}
        T2 = {'name' : 'T2', 'value' : cycle[-1], 'vary' : False}
        P1 = {'name' : 'P1', 'value' : 1000.0}
        t1 = {'name' : 't1', 'value' : pars['t1'].value, 'vary' : False}
        P2 = {'name' : 'P2', 'value' : 100.0}
        t2 = {'name' : 't2', 'value' : pars['t2'].value, 'vary' : False}    
        parameters = [T0, T1, T2, P1, t1, P2, t2]
        time_range=(cycle[0],cycle[-1])
        df = DecayFitter()
        if len(gate)==2:
            dE=1
        else:
            dE=gate[2]
        his1 = e.gy(his, gate_x=(gate[0],gate[1]), gate_y=time_range, bin_size=t_bin, plot=False)
        his2 = e.gy(his, gate_x=(gate[0],gate[1]), gate_y=time_range, bin_size=t_bin, plot=False)
        
        his1.histogram.weights*=0
        his2.histogram.weights*=0

        for i in range(gate[0],gate[1],dE):
            xgate = e.gx(his, gate_x=(i,i+dE), gate_y=time_range, bin_size=t_bin, plot=False)
            ygate = e.gy(his,time_range,gate_x=(i,i+dE-1), plot=False)
#            bckg = e.gx(his, gate_x=(i+dE,i+2*dE), gate_y=time_range, bin_size=t_bin, plot=False)    
            
            dyg = e._standard_errors_array(xgate.histogram.weights)
#            dyb = e._standard_errors_array(bckg.histogram.weights)

            gate_histo = copy.deepcopy(xgate)
#            gate_histo.weights = xgate.histogram.weights - bckg.histogram.weights
#        gate_histo.errors = numpy.sqrt(dyg**2 + dyb**2)
#            gate_histo.title = '{}: gx {} bg {} bin {}'.\
#                    format(his, gate[0], gate[1], t_bin)
#            plot_data = Plot(gate_histo, 'errorbar', True)

            t, n, p = df.fit(gate_histo.histogram.x_axis,
                gate_histo.histogram.weights, gate_histo.histogram.errors,
                'decay_only2', parameters)

            p1=p.get('P1').value
            p2=p.get('P2').value
            print(p1,p2)


            his1.histogram.weights[i:i+dE]+=ygate.histogram.weights*p1/(p1+p2)
            his2.histogram.weights[i:i+dE]+=ygate.histogram.weights*p2/(p1+p2)

        e.plotter.plot1d(his1, [cycle[0], cycle[-1]], None)
        e.plotter.plot1d(his2, [cycle[0], cycle[-1]], None)

        Experiment.plots.append(his1)
        Experiment.plots.append(his2)


def decay_only2(self, params, data_x):
        """Simplified bateman, decay part only, the second activity in the chain"""
        T0 = params['T0'].value
        T1 = params['T1'].value
        T2 = params['T2'].value
        N1 = params['P1'].value
        t1 = params['t1'].value
        N2 = params['P2'].value
        t2 = params['t2'].value
        y = []
        for t in data_x:
            if  t < T0:
                y.append(0)
            elif T0 <= t < T2:
                ts = t - T1
                y.append(N1 / (t1 - t2) * 
                        (numpy.exp(-ts / t1) - numpy.exp(-ts / t2)) +
                        N2 / t2 * numpy.exp(-ts/ t2))
            else:
                y.append(0)
        return numpy.array(y)      
  
def bateman(N1,N2,t1,t2):
        import numpy 
        y = []
        for t in range(0,130):
            y.append(N1 / (t1 - t2) * 
                (numpy.exp(-t / t1) - numpy.exp(-t / t2)) +
                N2 / t2 * numpy.exp(-t/ t2))
        return numpy.array(y)

def batemanA(N1,N2,t1,t2):
        import numpy 
        y = []
        for t in range(0,130):
            y.append(N1 / (t1 - t2) * 
                (numpy.exp(-t / t1) - numpy.exp(-t / t2)) +
                N2 / t2 * numpy.exp(-t/ t2) + N1 / t1 * numpy.exp(-t/ t1) )
        return numpy.array(y)

def bateman1(N1,N2,t1,t2):
        import numpy 
        y = []
        for t in range(0,130):
            y.append(N1 / (t1 - t2) * 
                (numpy.exp(-t / t1) - numpy.exp(-t / t2)))
        return numpy.array(y)

def bateman2(N1,N2,t1,t2):
        import numpy 
        y = []
        for t in range(0,130):
            y.append(N2 / t2 * numpy.exp(-t/ t2))
        return numpy.array(y)

def bateman3(N1,N2,t1,t2):
        import numpy 
        y = []
        for t in range(0,130):
            y.append(N1 / (t1) * (numpy.exp(-t / t1)) )
        return numpy.array(y)

def fit_decay(self, his, gate, cycle, 
                        t_bin=1, time_range=None,
                        model='grow_decay',
                        pars=None,
                        clear=True, mode="default"):
        """Fits decay time profile (grow-in/decay cycle):
        * his: is E-time histogram id
        * gate:  should be given in format:
            ((x0, x1), (bg0, bg1))
        * cycle: is list of beam start, beam stop, cycle end, e.g.
        (0, 100, 300)
        * t_bin: is a binning parameter (optional)
        * time_range: is a gate in time in (t0, t1) format (optional)
        * model: is a model used for fit (see decay_fitter)
                 (default is 'grow_decay')
        * pars is a list of dictionaries (one dict per each parameter)
        (optional, use if model is different than the default one, see
        decay_fitter for details)
        * if mode=="data" pass data to first variable (his) and background
              to be subtracted to second variable (gate).
            
        """
        if pars is None:
            T0 = {'name' : 'T0', 'value' : cycle[0], 'vary' : False}
            T1 = {'name' : 'T1', 'value' : cycle[1], 'vary' : False}
            T2 = {'name' : 'T2', 'value' : cycle[2], 'vary' : False}
            P1 = {'name' : 'P1', 'value' : 100.0}
            t1 = {'name' : 't1', 'value' : 100.0}
            parameters = [T0, T1, T2, P1, t1]
            if model == 'grow_decay2' or model == 'decay_only2':
                P2 = {'name' : 'P2', 'value' : 1000.0}
                t2 = {'name' : 't2', 'value' : 1000.0}
                parameters.append(P2)
                parameters.append(t2)
            
        else:
            parameters = pars

        df = DecayFitter()
        if mode == "default":
            xgate = e.gx(his, gate_x=gate[0], gate_y=time_range, bin_size=t_bin,
                                plot=False)
            bckg = e.gx(his, gate_x=gate[1], gate_y=time_range, bin_size=t_bin,
                              plot=False)
        elif mode =="data":
            xgate = his
            bckg = gate

        dyg = e._standard_errors_array(xgate.histogram.weights)
        dyb = e._standard_errors_array(bckg.histogram.weights)

        gate_histo = histogram.Histogram()
        gate_histo.x_axis = xgate.histogram.x_axis
        gate_histo.weights = xgate.histogram.weights - bckg.histogram.weights
        gate_histo.errors = numpy.sqrt(dyg**2 + dyb**2)
        if mode == "default":
            gate_histo.title = '{}: gx {} bg {} bin {}'.\
                    format(his, gate[0], gate[1], t_bin)
        elif mode == "data":
            gate_histo.title = 'data_mode'
        plot_data = Plot(gate_histo, 'errorbar', True)

        t, n, parameters = df.fit(gate_histo.x_axis, gate_histo.weights,
                                  gate_histo.errors, model, parameters)

        fit_histo = histogram.Histogram()
        fit_histo.x_axis = t
        fit_histo.weights = n
        fit_histo.title = e._replace_latex_chars('Fit: {}'.format(model))
        plot_fit = Plot(fit_histo, 'function', True)

        if clear:
            e.clear()

        e.plotter.plot1d(plot_fit, [cycle[0], cycle[2]], None)
        e.plotter.plot1d(plot_data, [cycle[0], cycle[2]], None)

        Experiment.plots.append(plot_fit)
        Experiment.plots.append(plot_data)

        return parameters
"""
        fit_histo = copy.deepcopy(xgate)
        fit_histo.histogram.x_axis = t
        fit_histo.histogram.weights = n
        fit_histo.histogram.title = e._replace_latex_chars('Fit: {}'.format(model))
        plot_fit = Plot(fit_histo, 'function', True)

        if clear:
            e.clear()

        e.plotter.plot1d(plot_fit, time_range, None)
        e.plotter.plot1d(plot_data, time_range, None)

        Experiment.plots.append(plot_fit)
        Experiment.plots.append(plot_data)

        return parameters
""" 

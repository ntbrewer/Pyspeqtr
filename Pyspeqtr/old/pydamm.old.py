#!/usr/bin/env python3
"""K. Miernik 2012
k.a.miernik@gmail.com
N.T Brewer 2014
nbrewer2@utk.edu
Distributed under GNU General Public Licence v3

This module is inteded to be loaded in an interactive interpreter session.
The ipython is strongly recommended. The pydamm is a python replacement for
DAMM programm.

"""

import math
import numpy
import sys
import matplotlib
#matplotlib.use('TKagg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import re

from collections import deque
from Pyspectr import hisfile as hisfile
from Pyspectr import histogram as histogram
from Pyspectr import plotter as plotter
from Pyspectr.exceptions import GeneralError as GeneralError
from Pyspectr.decay_fitter import DecayFitter as DecayFitter
from Pyspectr.peak_fitter import PeakFitter as PeakFitter
from Pyspectr.cal_parms_def import set_default as set_default
#from Pyspectr import QtPlotter as QtPlotter
class Plot:
    """
    The Plot class holds a set of data and parameters that are needed to
    display the data. 

    The bin_size attribute defines how much the histogram should be binned
    The norm attribute defines the normalization parameter used

    These parameters are altering the display only, the histogram 
    always keeps the original data.

    If the binned or normalized histogram data are needed for direct
    access, use functions avaible in the Histogram class.

    The mode defines the way the data are presented,
    'histogram' is displayed with steps-mid
    'function' with continuus line
    'errorbar' with yerrorbars
    'map' for 2D histograms

    """

    def __init__(self, histogram, mode, active):
        self.histogram = histogram
        self.mode = mode
        self.active = active

        self._bin_size = 1
        self._norm = 1


    @property
    def bin_size(self):
        return self._bin_size


    @bin_size.setter
    def bin_size(self, bs):
        if self.histogram.dim != 1:
            raise GeneralError('Currently only 1D histograms can be binned')

        if isinstance(bs, int):
            # You can only bin further the histogram, there is no
            # way back (load original data again)
            if bs > self._bin_size:
                self._bin_size = bs
                self.histogram = self.histogram.rebin1d(bs)
            elif bs <= self._bin_size:
                pass
            else:
                raise GeneralError('Attempt to set bin size to {}'.\
                        format(bs))
        else:
            raise GeneralError('Attempt to set bin size to {}'.\
                    format(bs))


    @property
    def norm(self):
        return self._norm


    @norm.setter
    def norm(self, n):
        if self.histogram.dim != 1:
            raise GeneralError('Currently only 1D histograms can be normalized')
        self.histogram = self.histogram.normalize1d(n, self.bin_size)



    def __str__(self):
        """Basic information Informative string"""
        string = 'Plot: {} bin {} norm {:.2e}'.\
                    format(self.histogram.title.strip(), self.bin_size,
                           self.norm)
        return string


    def __repr__(self):
        """More verbose string

        """
        string = 'Plot: "{}" bin {} norm {:.2e} active {} mode {}'.\
                    format(self.histogram.title.strip(), self.bin_size,
                           self.norm, self.active, self.mode)
        return string



class Experiment:
    """Main class for data visualization and analysis
       Now uses online keyword for toggling between plotting 
       with matplotlib and pyqtgraph. 
       usage: Experiment('his_name',online=True/False, size=#). 
    """
    # Deque lengths
    FIFO_1D = 50
    FIFO_2D = 5

    # These variables and registers are class wide, so
    # if more than one Experiment is open (one file = one Experiment)
    # they share the registers and auto-scaling is still working

    # Keeps current and past 1D plots
    # The right-most item is the current one
    plots = deque(maxlen=FIFO_1D)

    # Current and past 2D plots
    maps = deque(maxlen=FIFO_2D)

    # 1D plot ranges
    xlim = None
    ylim = None

    # 2D plot ranges
    xlim2d = None
    ylim2d = None
    logz = False

    # 1 for 1D, 2 for 2D
    _mode = 1

    #parameters for calibrating via wavelet method
    #import cal_parms_def as cpd
    cal_parms, fltr_parms, srch_parms =set_default()
    def __init__(self, file_name, online=False,size=11):
        """Initialize, open data file (his) and open plot window
        (size parameter decides on plot dimensions)

        """
        if online==False:
           print('')
        elif online!=True:
           online=False
        else:
           print('online mode')

        self.file_name = file_name
        # The current (active) file
        self.hisfile = None
        self.load(file_name)

        # Peaks for fitting
        self.peaks = []
        #Set online mode
        self.online = online
        # plotter front-end
        self.plotter = plotter.Plotter(size,online)
        #self.QtPlotter = QtPlotter.QtPlotter(size)


    def load(self, file_name):
        """Load his file (also tar gzipped files)"""
        self.hisfile = hisfile.HisFile(file_name)


    @property
    def mode(self):
        """ 1D or 2D plotting mode"""
        return Experiment._mode


    @mode.setter
    def mode(self, mode):
        """Deactivate all plots that have different mode (dimension)"""
        if mode not in [1, 2]:
            raise GeneralError('Only 1D and 2D plotting modes are possible')

        if mode == 2:
            self.plotter.ylin()

        Experiment._mode = mode
        for p in self.plots:
            if p.histogram.dim != mode:
                p.active = False


    def _replace_latex_chars(self, text):
        """Clear text from characters that are not accepted by latex"""
        replace_chars = [['_', '-'],
                            ['$', '\$'],
                            ['%', '\%'],
                            ['~', ' '],
                            ['"', "''"],
                            ['\\', ' ']]
        replaced_text = text
        for r_ch in replace_chars:
            replaced_text = replaced_text.replace(r_ch[0], r_ch[1])
        return replaced_text


    def show_registers(self):
        """Print the available registers"""
        i = -1
        print('1D histograms')
        print('{: <3} {: ^40} {: ^5} {: ^8} {: ^8}'.\
                format('i', 'Title', 'Bin', 'Norm', 'Active'))
        print('-' * 79)

        for p in reversed(Experiment.plots):
            print('{: >3} {: <40} {: >5} {: >5.2e} {: >5}'.\
                    format(i, p.histogram.title[:40], p.bin_size,
                           p.norm, p.active))
            i -= 1
        print()

        i = -1
        print('2D histograms')
        print('{: <3} {: ^40} {: ^5} {: ^8} {: ^8}'.\
                format('i', 'Title', 'Bin', 'Norm', 'Active'))
        print('-' * 79)

        for p in reversed(Experiment.maps):
            print('{: >3} {: <40} {: >5} {: >5.2e} {: >5}'.\
                    format(i, p.histogram.title[:40], p.bin_size,
                           p.norm, p.active))
            i -= 1
        print()


    def _expand_norm(self, norm, num_of_args):
        """Return normalization array of lenght equal to 
        num_of_args, expand integers to whole array, check 
        if list is of proper lenght

        """
        normalization = []
        if isinstance(norm, str):
            if norm.lower() == 'area':
                for i in range(num_of_args):
                    normalization.append('area')
            else:
                print("Normalization must be a float, ",
                    "list of floats or a 'area' string")
                return None
        elif isinstance(norm, float) or isinstance(norm, int):
            for i in range(num_of_args):
                normalization.append(norm)
        elif isinstance(norm, list):
            if len(norm) == num_of_args:
                normalization = norm
            else:
                print('List of normalization factors must be of the same' +
                      ' length as the list of histograms')
                return None
        else:
            print("Normalization must be a float, ",
                  "list of floats or a 'area' string")
            print(norm, ' was given')
            return None
        return normalization


    def _expand_bin_sizes(self, bin_size, num_of_args):
        """See _expand_norm"""
        bin_sizes = []
        if isinstance(bin_size, int):
            for i in range(num_of_args):
                bin_sizes.append(bin_size)
        elif isinstance(bin_size, list):
            if len(bin_size) == num_of_args:
                bin_sizes = bin_size
            else:
                print('List of bin sizes must be of the same' +
                      ' length as the list of histograms')
                return None
        else:
            print("Bin size must be an int or a list of ints")
            return None
        return bin_sizes


    def _expand_d_args(self, args):
        """Expand list of args to a list of histograms ids or Plot
        instances"""
        his_list = []
        for his in args:
            if isinstance(his, int):
                his_list.append(his)
            elif isinstance(his, str):
                try:
                    his_range = his.split('-')
                    his_range = [x for x in range(int(his_range[0]),
                                                int(his_range[1]) + 1)]
                    his_list += his_range
                except (ValueError, IndexError):
                    break
            elif isinstance(his, Plot):
                his_list.append(his)
            else:
                break
        else:
            return his_list
        print("Histogram list must be given in a 'x-y' format,",
                "where x and y are integers",
                "(note also quotation marks), e.g. '100-115'")
        return None



    def d(self, *args, norm=1, bin_size=1, clear=True):
        """
        Plot 1D histogram. 
        * args: is a list of histograms that may be given as:
              - positive integer: is interpreted as the histogram id
                                  from a currently open file
              - negative integer: is interpreted as the registry number
                                  (see (show_registers())
              - Plot object:       see Plot class
              - string:  in 'x-y' format where x and y are integers 
                        (note also mandatory quatation marks)
                        is interpreted as a range of histograms ids

        * norm: may be given as a single float or int value or an 'area' string,
                also a list of length matching the *args list may be used
                with any combination of the above accepted values
        * bin_size: must be an integer, a list of ints is 
                    also accepted (see norm)
        * clear: is True by default, which means that previous plot is 
                 cleared if False is given, the previous plots are not cleared.

        Example:
        e.d(100, plot1, '105-106', -3, bin_size=[1, 2, 1, 1, 10], clear=False)

        """
        plots = []

        his_list = self._expand_d_args(args)

        normalization = self._expand_norm(norm, len(his_list))
        if normalization is None:
            return None

        bin_sizes = self._expand_bin_sizes(bin_size, len(his_list))
        if bin_sizes is None:
            return None

        # Clear the plotting area (of clear is False, the currently
        # active plots are not deactivated, so they got replotted at
        # the end of this function)
        #if self.online == False:
        self.plotter.clear()

        # Switch mode to 1D
        self.mode = 1
        # Deactivate current plots if clear flag is used
        if clear:
            for p in Experiment.plots:
                p.active = False

        # Prepare data for plotting
        for i_plot, his in enumerate(his_list):
            if isinstance(his, int):
                # load histograms from the file
                if his > 0:
                    data = self.hisfile.load_histogram(his)
                    if data[0] != 1:
                        print('{} is not a 1D histogram'.format(his))
                        return None
                    title = self.hisfile.histograms[his]['title'].strip()
                    title = '{}:{}'.format(his, 
                                           self._replace_latex_chars(title))
                    histo = histogram.Histogram()
                    histo.title = title
                    histo.x_axis = data[1]
                    histo.weights = data[3]
                    histo.errors = self._standard_errors_array(data[3])
                    plot = Plot(histo, 'histogram', True)
                    plot.bin_size = bin_sizes[i_plot]
                    plot.norm = normalization[i_plot]
                    plots.append(plot)
                    Experiment.plots.append(plot)
                else:
                    # plot histograms from registry
                    # Numbered by negative numbers (-1 being the latest)
                    # Call show_registers for a list of available plots
                    try:
                        plot = Experiment.plots[his]
                        Experiment.plots[his].active = True
                        Experiment.plots[his].bin_size = bin_sizes[i_plot]
                        Experiment.plots[his].norm = normalization[i_plot]
                    except IndexError:
                        print('There is no plot in the registry under the',
                              'number', his, 'use show_registry() to see',
                              'available plots')
                        return None
                    plots.append(plot)
            elif isinstance(his, Plot):
                # If instance of Plot class is given, mark it active and add
                # to the deque (if not already there)
                # and to the array to be returned at the end
                his.active = True
                his.bin_size = bin_sizes[i_plot]
                his.norm = normalization[i_plot]
                plots.append(his)
                if his not in Experiment.plots:
                    Experiment.plots.append(his)

        # Count the number of active plots
        active_plots = 0
        for plot in Experiment.plots:
            if plot.active:
                active_plots += 1

        # Here the actual plotting happens
        i_plot = 0
        for plot in Experiment.plots:
            if plot.active:
                i_plot += 1
                # If ylim is not given explicitely, go through the
                # active plots to find the plot limits
                # This is run only for the last plot.
                # Note that this is neccesary as matplotlib is not
                # autoscaling Y axis when 
                # changing the X axis is being changed
                # If, in a future, the behaviour of matplotlib
                # changes, this part may dropped
                #if self.online==False:
                ylim = None
                if self.ylim is None and i_plot == active_plots:
                    ylim = self._auto_scale_y()
                else:
                    ylim = self.ylim
                #else:
                #    ylim=None
                # Note that ylim is autoscaled above if self.ylim is None
                # But we still keep self.ylim None, 
                # to indicate autoscaling
                self.plotter.plot1d(plot, Experiment.xlim, ylim)

        # Return plots that were added or activated
        return plots


    def _auto_scale_y(self):
        """Find the y limits taking into account all active plots """
        ylim = [None, None]
        for p in Experiment.plots:
            if p.active:
                histo = p.histogram
                #if p.bin_size > 1:
                    #histo = histo.rebin1d(p.bin_size)
                #if p.norm != 1:
                    #histo = histo.normalize1d(p.norm, p.bin_size)
                if Experiment.xlim is None:
                    ymin = min(histo.weights)
                    ymax = max(histo.weights)
                else:
                    i_xmin = Experiment.xlim[0] // p.bin_size - 1
                    if i_xmin < 0:
                        i_xmin = 0
                    i_xmax = Experiment.xlim[1] // p.bin_size + 1
                    try:
                        ymin = min(histo.weights[i_xmin:i_xmax])
                    except ValueError:
                        ymin = None
                    try:
                        ymax = max(histo.weights[i_xmin:i_xmax])
                    except ValueError:
                        ymax = None
                if ymin is not None:
                    if ylim[0] is not None:
                        if ymin < ylim[0]:
                            ylim[0] = ymin
                    else:
                        ylim[0] = ymin
                if ymax is not None:
                    if ylim[1] is not None:
                        if ymax > ylim[1]:
                            ylim[1] = ymax
                    else:
                        ylim[1] = ymax
        if ylim[0] is None or ylim[1] is None:
            return None
        else:
            return [ylim[0] - ylim[0] * 0.1, ylim[1] + ylim[1] * 0.1]


    def _auto_scale_x(self):
        """Find the x axis limits taking into account all active plots."""
        xlim = [None, None]
        for p in Experiment.plots:
            if p.active:
                histo = p.histogram
                if Experiment.xlim is None:
                    xmin = histo.x_axis[0]
                    xmax = histo.x_axis[-1]
                    if xlim[0] is not None:
                        if xmin < xlim[0]:
                            xlim[0] = xmin
                    else:
                        xlim[0] = xmin
                    if xlim[1] is not None:
                        if xmax > xlim[1]:
                            xlim[1] = xmax
                    else:
                        xlim[1] = xmax

        if xlim[0] is None or xlim[1] is None:
            return None
        else:
            return xlim


    def dl(self, x0=None, x1=None):
        """Change x range of 1D histogram"""
        if self.mode != 1:
            return None

        if x0 is None or x1 is None:
            Experiment.xlim = None
            self.plotter.xlim(self._auto_scale_x())
        else:
            Experiment.xlim = (x0, x1)
            self.plotter.xlim(Experiment.xlim)

        if self.ylim is None:
            self.plotter.ylim(self._auto_scale_y())


    def dmm(self, y0=None, y1=None):
        """Change yrange of 1D histogram """
        if self.mode != 1:
            return None

        if y0 is None or y1 is None:
            self.ylim = None
        else:
            self.ylim = (y0, y1)

        if self.ylim is None:
            self.plotter.ylim(self._auto_scale_y())
        else:
            self.plotter.ylim(self.ylim)


    def log(self):
        """Change y scale to log or z scale to log"""
        if self.mode == 1:
            self.plotter.ylog()
            self.dmm()
        elif self.mode == 2:
            Experiment.logz = True
            self.dd(-1, xc=Experiment.xlim2d, yc=Experiment.ylim2d)


    def lin(self):
        """Change y scale to linear or z scale to linear"""
        if self.mode == 1:
            self.plotter.ylin()
        if self.mode == 2:
            Experiment.logz = False
            self.dd(-1, xc=Experiment.xlim2d, yc=Experiment.ylim2d)


    def list(self, his_id=None):
        """List all histograms in the active data file
           or details on a selected histogram. Now accepts '1d','1D','2d',
           or '2D' as input for listing all 1 or 2 dimensional histograms. 
           Now also implements re to query the list (case insensitive)"""
        if his_id is None:
            for key in sorted(self.hisfile.histograms.keys()):
                print('{: <6} {}'.format(key, 
                                    self.hisfile.histograms[key]['title']))
        elif isinstance(his_id, str):
            if his_id is '1d' or his_id is'1D':
               for key in sorted(self.hisfile.histograms.keys()):
                   if self.hisfile.histograms[key]['dimension']==1:
                       print('{: <6} {}'.format(key, 
                                    self.hisfile.histograms[key]['title']))
            elif his_id is '2d' or his_id is'2D':
               for key in sorted(self.hisfile.histograms.keys()):
                   if self.hisfile.histograms[key]['dimension']==2:
                       print('{: <6} {}'.format(key, 
                                    self.hisfile.histograms[key]['title']))
            else:
               for key in sorted(self.hisfile.histograms.keys()):
                   if re.search(his_id,self.hisfile.histograms[key]['title'],re.I)!=None:
                      print('{: <6} {}'.format(key, 
                                    self.hisfile.histograms[key]['title']))
        else:
            try:
                dim = self.hisfile.histograms[his_id]['dimension']
                xmin = []
                xmax = []
                for i in range(dim):
                    xmin.append(self.hisfile.histograms[his_id]['minc'][0])
                    xmax.append(self.hisfile.histograms[his_id]['maxc'][0])
                print('{: <10} : {}'.format('ID', his_id))
                print('{: <10} : {}'.format('Title', 
                                    self.hisfile.histograms[his_id]['title']))
                print('{: <10} : {}'.format('Dimensions', dim))
                print('{: <10} : ({}, {})'.format('X range', xmin[0], xmax[0]))
                if dim > 1:
                    print('{: <10} : ({}, {})'.format('Y range', 
                                                      xmin[1], xmax[1]))
            except KeyError:
                print('Histogram id = {} not found'.format(his_id))

    def _standard_errors_array(self, data):
        """ Calculate standard error array (\sigma_i = \sqrt{n_i}),
           with a twist: if n_i = 0, the uncertainity is 1 (not 0)

        """
        errors = numpy.zeros(data.shape)
        for index, d in numpy.ndenumerate(data):
            if d == 0:
                errors[index] = 1
            else:
                errors[index] = math.sqrt(abs(d))
        return errors


    def _add_errors(self, error1, error2):
        """Add two error arrays
        \sigma = \sqrt{\sigma_1^2 + \sigma_2^2}

        """
        if error1.shape != error2.shape:
            raise GeneralError('Shape of array mismatches')
        errors = numpy.zeros(error1.shape)
        for index, d in numpy.ndenumerate(error1):
            errors[index] = math.sqrt(error1[index]**2 + error2[index]**2)
        return errors


    def gx(self, his, gate_x, gate_y=None, bg_gate=None, norm=1,
           bin_size=1, clear=True, plot=True):
        """Make projection on Y axis of 2D histogram with gate
        set on X (gate_x) and possibly on Y (gate_y)

        his: is a histogram id in a file
        gate_x: is range of bins in (x0, x1) format, this selects the
                range of X columns to be projected on Y axis
        gate_y: is a range of bins in (y0, y1) format (optional), this
                truncates the range of the projection along the Y axis
        bg_gate: is a range of bins in (x0, x1) format (optional), this
                selects the background gate that is subtracted from the
                selected gate_x
        norm: normalization factor (see d())
        bin_size: binning factor (see d())
        clear: True by default, clears previous plots
        plot: True by default, if False no plotting is taking place, 
              only the plot object is being returned
        
        """
        if gate_x is None or len(gate_x) != 2:
            print('Please select gate on X in a (min, max) format')
            return None
        if gate_y is not None and len(gate_y) != 2:
            print('Please select gate on Y in a (min, max) format')
            return None

        # If clear flag used, clear the plotting area
        if clear and plot:
            self.plotter.clear()

        # Switch mode to 1D
        self.mode = 1
        # Deactivate all plots if clear flag is used
        if clear and plot:
            for p in Experiment.plots:
                p.active = False

        data = self.hisfile.load_histogram(his)
        if data[0] != 2:
            print('{} is not a 2D histogram'.format(his))
            return None

        # x for x_axis data
        # y for y_axis data
        # w for weights
        # g for gate (result)
        # bg for background gate
        x = data[1]
        y = data[2]
        w = data[3]
        if gate_y is None:
            gate_y = [0, len(y)-2]
        y = y[gate_y[0]:gate_y[1]+1]
        g = w[gate_x[0]:gate_x[1]+1, gate_y[0]:gate_y[1]+1].sum(axis=0)
        dg = self._standard_errors_array(g)
        if bg_gate is not None:
            if (bg_gate[1] - bg_gate[0]) != (gate_x[1] - gate_x[0]):
                print('#Warning: background and gate of different widths')
            bg = w[bg_gate[0]:bg_gate[1]+1, gate_y[0]:gate_y[1]+1].sum(axis=0)
            g = g - bg
            # Note that since the gate is adding bins, the formula
            # used for standard error is no longer valid
            # This approximation should be good enough though
            dbg = self._standard_errors_array(bg)
            dg = self._add_errors(dg, dbg)

        title = '{}:{} gx({},{})'.format(his, self.hisfile.\
                                         histograms[his]['title'].strip(),
                                         gate_x[0], gate_x[1])
        if bg_gate is not None:
            title += ' bg ({}, {})'.format(bg_gate[0], bg_gate[1])
        title = self._replace_latex_chars(title)

        histo = histogram.Histogram()
        histo.title = title
        histo.x_axis = y
        histo.weights = g
        histo.errors = dg
        gate_plot = Plot(histo, 'histogram', True)
        gate_plot.bin_size = bin_size
        gate_plot.norm = norm

        if plot:
            Experiment.plots.append(gate_plot)
            ylim = None
            if self.ylim is None:
                ylim = self._auto_scale_y()
            else:
                ylim = self.ylim
            self.plotter.plot1d(gate_plot, Experiment.xlim, ylim)

        return gate_plot


    def gy(self, his, gate_y, gate_x=None, bg_gate=None, norm=1,
           bin_size=1, clear=True, plot=True):
        """Make projection on X axis of 2D histogram with gate
        set on Y (gate_y) and possibly on X (gate_x)
        
        see gx for more details
        """
        if gate_y is None or len(gate_y) != 2:
            print('Please select gate on Y in a (min, max) format')
            return None
        if gate_x is not None and len(gate_x) != 2:
            print('Please select gate on X in a (min, max) format')
            return None

        # If clear flag used, clear the plotting area
        if clear and plot:
            self.plotter.clear()

        # Switch mode to 1D
        self.mode = 1
        # Deactivate all plots if clear flag is used
        if clear and plot:
            for p in Experiment.plots:
                p.active = False

        data = self.hisfile.load_histogram(his)
        if data[0] != 2:
            print('{} is not a 2D histogram'.format(his))
            return None

        # x for x_axis data
        # y for y_axis data
        # w for weights  
        # g for gate (result)
        # bg for background gate
        x = data[1]
        y = data[2]
        w = data[3]
        if gate_x is None:
            gate_x = [0, len(x)-2]
        x = x[gate_x[0]:gate_x[1]+1]
        g = w[gate_x[0]:gate_x[1]+1, gate_y[0]:gate_y[1]+1].sum(axis=1)
        dg = self._standard_errors_array(g)
        if bg_gate is not None:
            if (bg_gate[1] - bg_gate[0]) != (gate_y[1] - gate_y[0]):
                print('#Warning: background and gate of different widths')

            bg = w[gate_x[0]:gate_x[1]+1, bg_gate[0]:bg_gate[1]+1].sum(axis=1)
            g = g - bg
            # Note that since the gate is adding bins, the formula
            # used for standard error is no longer valid
            # This approximation should be good enough though
            dbg = self._standard_errors_array(bg)
            dg = self._add_errors(dg, dbg)

        title = '{}:{} gy({},{})'.format(his, self.hisfile.\
                                         histograms[his]['title'].strip(),
                                         gate_y[0], gate_y[1])
        if bg_gate is not None:
            title += ' bg ({}, {})'.format(bg_gate[0], bg_gate[1])
        title = self._replace_latex_chars(title)

        histo = histogram.Histogram()
        histo.title = title
        histo.x_axis = x
        histo.weights = g
        histo.errors = dg
        gate_plot = Plot(histo, 'histogram', True)
        gate_plot.bin_size = bin_size
        gate_plot.norm = norm

        Experiment.plots.append(gate_plot)
        if plot:
            ylim = None
            if self.ylim is None:
                ylim = self._auto_scale_y()
            else:
                ylim = self.ylim
            self.plotter.plot1d(gate_plot, Experiment.xlim, ylim)

        return gate_plot


    def mark(self, x_mark):
        """Put vertical line on plot to mark the peak (or guide the eye)"""
        plt.axvline(x_mark, ls='--', c='black')


    def annotate(self, x, text, shiftx=0, shifty=0):
        """ Add arrow at x, with annotation text"""
        if self.mode != 1:
            print('Annotation works only for 1D histograms')
            return None

        length = 0.07 * (plt.ylim()[1] - plt.ylim()[0])
        y = self.plots[-1].histogram.weights[x // self.plots[-1].bin_size]
        plt.annotate(text, xy=(x, y),
                    xytext=(x + shiftx, y + length + shifty),
                    rotation=90.,
                    xycoords='data',
                    fontsize=9,
                    verticalalignment='bottom',
                    horizontalalignment='center',
                    arrowprops=dict(width=1, facecolor='black', headwidth=5,
                                    shrink=0.1))


    def load_gates(self, filename):
        """Load gamma gates from text file, the format is:
        # Comment line
        Name    x0  x1  bg0 bg1
        Example:
        110     111 113 115 117

        """
        gatefile = open(filename, 'r')
        lineN = 0
        gates = {}
        for line in gatefile:
            lineN += 1
            line = line.strip()
            if line.startswith('#'):
                continue
            items = line.split()
            if len(items) < 5:
                print('Warning: line {} bad data'.format(lineN))
                continue
            gates[int(items[0])] = ((int(items[1]), int(items[2])),
                                   (int(items[3]), int(items[4])))
        return gates


    def pk(self, *args, **kwargs):
        """Add peaks for gaussian fitting procedure. The args
        give a list of peak energy (approx.), the kwargs may include
        additional parameters (e.g. min or max, etc) used by peak_fitter"""
        for e in args:
            if isinstance(e, int) or isinstance(e, float):
                p = {'E' : e}
                p.update(kwargs)
                self.peaks.append(p)


    def pzot(self):
        """Clear all peaks """
        self.peaks=[]


    def dd(self, his, xc=None, yc=None, logz=None):
        """Plot 2D histogram,

        his may be a positive integer (loads histogram from the data file)
        negative integer (2D plots registry) or Plot instance (must be a 2D
        plot)

        xc is x range, yc is y range, that may be applied immediately, 
        see also xc() and yc() functions
        
        """
        self.mode = 2

        for p in Experiment.maps:
            p.active = False

        plot = None
        self.plotter.clear()

        if isinstance(his, int):
            if his > 0:
                data = self.hisfile.load_histogram(his)
                if data[0] != 2:
                    print('{} is not a 2D histogram'.format(his))
                    return None

                title = self.hisfile.histograms[his]['title'].strip()
                title = '{}:{}'.format(his, 
                                       self._replace_latex_chars(title))
                histo = histogram.Histogram(dim=2)
                histo.title = title
                histo.x_axis = data[1]
                histo.y_axis = data[2]
                histo.weights = data[3]
                plot = Plot(histo, 'map', True)
                Experiment.maps.append(plot)
            else:
                # plot histogram from the registry
                # Numbered by negative numbers (-1 being the latest)
                # Call show_registers for a list of available plots
                try:
                    plot = Experiment.maps[his]
                    Experiment.maps[his].active = True
                except IndexError:
                    print('There is no 2D plot in the registry under the',
                            'number', his, 'use show_registry() to see',
                            'available plots')
                    return None
        elif isinstance(his, Plot):
            # If instance of Plot class is given, mark it active and add
            # to the deque (if not already there)
            # and to the array to be returned at the end
            if his.histogram.dim != 2:
                print('This {} is not a 2D histogram'.format(his))
                return None
            his.active = True
            plot = his
            if his not in Experiment.maps:
                Experiment.maps.append(his)

        if xc is not None:
            Experiment.xlim2d = xc
        if yc is not None:
            Experiment.ylim2d = yc

        if logz is None:
            use_log = Experiment.logz
        else:
            use_log = logz
        if plot is not None:
            self.plotter.plot2d(plot, Experiment.xlim2d, 
                                Experiment.ylim2d, use_log)

        return [plot]


    def xc(self, x0=None, x1=None):
        """Change xrange of a 2D histograms"""
        if self.mode == 2:
            if x0 is None or x1 is None:
                Experiment.xlim2d = None
                xlim = None
                for p in Experiment.maps:
                    if p.active:
                        histo = p.histogram
                        xlim = (histo.x_axis[0], histo.x_axis[-1])
                        break
            else:
                Experiment.xlim2d = (x0, x1)
                xlim = (x0, x1)
            self.dd(-1, xc=xlim, yc=Experiment.ylim2d)


    def yc(self, y0=None, y1=None):
        """Change yrange of a 2D histogram"""
        if self.mode == 2:
            if y0 is None or y1 is None:
                Experiment.ylim2d = None
                ylim = None
                for p in Experiment.maps:
                    if p.active:
                        histo = p.histogram
                        ylim = (histo.y_axis[0], histo.y_axis[-1])
                        break
            else:
                Experiment.ylim2d = (y0, y1)
                ylim = (y0, y1)
            self.dd(-1, xc=Experiment.xlim2d, yc=ylim)



    def clear(self):
        self.plotter.clear()


    def color_map(self, cmap=None, clist=False):
        if self.mode == 2:
            self.plotter.color_map(cmap)
            self.dd(-1, xc=Experiment.xlim2d, yc=Experiment.ylim2d)
        if clist:
            maps=[m for m in cm.datad if not m.endswith("_r")]
            print(maps)
               
            

    def fit_peaks(self, his=None, rx=None, clear=True, model='linear'):
        """ 
        Fit gaussian peaks to 1D plot. If his is not given the
        current plot is used. If rx is not given, the current range is used
        Returns list of lists:
            [E, x0, dx, A, dA, s, Area]
        where E is name of the peak, x0, A and s are fitted parameters
        and d'something' is its uncertainity. Area is total calculated area.

        """
        import copy
        cdc = copy.deepcopy

        if rx is None:
            rx = Experiment.xlim
        if len(rx) != 2:
            print('Please use x range in format rx=(min, max), where',
                  'min and max are integers.')
            return None

        # Deactivate all the plots
        for p in Experiment.plots:
            if p.active:
                p.active = False

        peaks = []
        for p in self.peaks:
            if rx[0] <= p.get('E') <= rx[1]:
                peaks.append(p)

        PF = PeakFitter(peaks, model, '')

        if his is not None:
            if isinstance(his, int):
                if his > 0:
#                    data = self.hisfile.load_histogram(his)
                    if self.hisfile.load_histogram(his)[0] != 1:
                        print('{} is not a 1D histogram'.format(his))
                        return None                    
                    self.d(his)
                    data = cdc(Experiment.plots[-1])
                    x_axis = data.histogram.x_axis
                    weights = data.histogram.weights
                    title = data.histogram.title
#                    x_axis = data[1]
#                    weights = data[3]
#                    title = self.hisfile.histograms[his]['title'].strip()
#                    title = '{}:{}'.format(his,
#                                           self._replace_latex_chars(title))
                else:
                    try:
                        data = cdc(Experiment.plots[his])
                        x_axis = data.histogram.x_axis
                        weights = data.histogram.weights
                        title = data.histogram.title
                    except IndexError:
                        data=0
                        x_axis=0
                        weights=0
                        title='Err'
                        print('There is no plot in the registry under the',
                              'number', his, 'use show_registry() to see',
                              'available plots')
                        return None
            elif hasattr(his,'histogram'):
#                print('ok')
                data = cdc(his)
                x_axis = data.histogram.x_axis
                weights = data.histogram.weights
                title = data.histogram.title

        else:
            data = cdc(Experiment.plots[-1])
            x_axis = data.histogram.x_axis
            weights = data.histogram.weights
            title = data.histogram.title
            
        dweights = self._standard_errors_array(weights)

        if clear:
            self.clear()

#        histo_data = histogram.Histogram()
#        histo_data.x_axis = x_axis
#        histo_data.weights = weights
#        histo_data.errors = dweights
#        histo_data.title = title
#        plot_data = Plot(histo_data, 'histogram', True)
        # The histogram data is plotted here so the fit function
        # may be overlaid on in. However, the plot_data is appended 
        # to the registry after the fit functions so it is on top of the
        # registry.
#        self.plotter.plot1d(plot_data, xlim=rx)
        self.plotter.plot1d(data, xlim=rx)

        fit_result = PF.fit(x_axis[rx[0]:rx[1]], weights[rx[0]:rx[1]],
                            dweights[rx[0]:rx[1]])
#        if hasattr(his,'histogram'):
#            fit_result = PF.fit(histo_data.x_axis, histo_data.weights, histo_data.errors)
#        if hasattr(data,'histogram'):
#            fit_result = PF.fit(histo_data.x_axis, histo_data.weights, histo_data.errors)

#        histo_baseline = histogram.Histogram()
#        histo_baseline.x_axis = x_axis
        histo_baseline = cdc(data)
        histo_baseline.histogram.x_axis = x_axis[rx[0]:rx[1]]
        histo_baseline.histogram.weights = fit_result['baseline']
        histo_baseline.histogram.title = 'Baseline'
#        plot_baseline = Plot(histo_baseline, 'function', True)
        self.plotter.plot1d(histo_baseline)

#        histo_peaks = histogram.Histogram()
        histo_peaks = cdc(data)
        histo_peaks.histogram.x_axis = fit_result['x_axis']
        histo_peaks.histogram.weights = fit_result['fit']
        histo_peaks.histogram.title = 'Fit'
#        plot_peaks = Plot(histo_peaks, 'function', True)
        self.plotter.plot1d(histo_peaks)
        # Append all the plots to the registry, but
        # keep original data at the end, so the next fit_peaks()
        # call will use then again as default
        Experiment.plots.append(histo_baseline)
        Experiment.plots.append(histo_peaks)
        Experiment.plots.append(data)

        # Plot the last one with the auto_scale if needed
        if Experiment.ylim is None:
            ylim = self._auto_scale_y()
        else:
            ylim = Experiment.ylim

#        self.plotter.plot1d(histo_peaks)



        print('#{:^8} {:^8} {:^8} {:^8} {:^8} {:^8} {:^8}'
                .format('Peak', 'x0', 'dx', 'A', 'dA', 's', 'Area'))
        peak_data = []
        for i, peak in enumerate(peaks):
            if peak.get('ignore') == 'True':
                continue
            x0 = PF.params['x{}'.format(i)].value
            dx = PF.params['x{}'.format(i)].stderr
            A = PF.params['A{}'.format(i)].value
            dA = PF.params['A{}'.format(i)].stderr
            s = PF.params['s{}'.format(i)].value
            Area = PF.find_area(x_axis, i)
            print('{:>8} {:>8.2f} {:>8.2f} {:>8.1f} {:>8.1f} {:>8.3f} {:>8.1f}'
                    .format(peaks[i].get('E'), x0, dx, A, dA, s, Area))
            peak_data.append([peaks[i].get('E'), x0, dx, A, dA, s, Area])
        return peak_data


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
            if model == 'grow_decay2' or model == 'decay_only2' or model == 'decay_onlyAll'  or model == 'decay_onlyAll2':
                P2 = {'name' : 'P2', 'value' : 1000.0}
                t2 = {'name' : 't2', 'value' : 1000.0}
                parameters.append(P2)
                parameters.append(t2)
            if model == 'decay_onlyAll2':
                P3 = {'name' : 'P3', 'value' : 1000.0}
                parameters.append(P3)
            if model == 'grow_decay_offset':
                TOFF = {'name' : 'TOFF', 'value' : 0.0}
                parameters.append(TOFF)
            if model == 'decay_only_linbg':
                m = {'name' : 'm', 'value' : 1.0}
                b = {'name' : 'b', 'value' : 10.0}
                parameters.append(m)
                parameters.append(b)
            
        else:
            parameters = pars

        df = DecayFitter()
        if mode == "default":
            xgate = self.gx(his, gate_x=gate[0], gate_y=time_range, bin_size=t_bin,
                                plot=False)
            bckg = self.gx(his, gate_x=gate[1], gate_y=time_range, bin_size=t_bin,
                              plot=False)
        elif mode =="data":
            xgate = his
            bckg = gate

        dyg = self._standard_errors_array(xgate.histogram.weights)
        dyb = self._standard_errors_array(bckg.histogram.weights)

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

        t, n, parameters, rcs, pars = df.fit(gate_histo.x_axis, gate_histo.weights,
                                  gate_histo.errors, model, parameters)
  
        fit_histo = histogram.Histogram()
        fit_histo.x_axis = t
        fit_histo.weights = n
        fit_histo.title = self._replace_latex_chars('Fit: {}'.format(model))
        plot_fit = Plot(fit_histo, 'function', True)

        if clear:
            self.clear()

        self.plotter.plot1d(plot_fit, [cycle[0], cycle[2]], None)
        self.plotter.plot1d(plot_data, [cycle[0], cycle[2]], None)

        Experiment.plots.append(plot_fit)
        Experiment.plots.append(plot_data)

        return (parameters, rcs, pars)

    def gamma_gamma_spectra(self, gg_id, gate, bin_size=1):
        """ 
        Plots gamma-gamma gate broken into 4 subplots (0-600, 600-1200,
        1200-2000, 2000-4000. 
        gg_id is a 2D histogram id
        gate is in form ((x1, y1), (x2, y2)) where i=1 is gate on line, i=2
        is gate on background

        This special plot is not loaded into the registry in a 4 panel
        form, but as a 'standard' plot object
        """
        self.clear()
        plot = self.gy(gg_id, gate[0], bg_gate=gate[1], 
                       bin_size=bin_size, plot=False )
        ranges = (0, 600, 1200, 2000, 4000)
        self.plotter.plot1d_4panel(plot, ranges)

    def st(self,numx,numy,args):
        """
        This command is similar but different to the "st" or stack
        command in damm. Create an array of subplots numx by numy.
        Then disply a multiplot array of list *args). This should work for 
        1 and 2d histograms but ought to be cleaned up.
        """

        self.clear()
        for i in range(numx):
            for j in range(numy):
                if len(args)!=numx*numy: 
                   print('range mismatch')
                   break
                n = i+numx*(j-1) if numy!=1 else j+numy*(i-1) 
                ax = plt.subplot(numx,numy,n)
                ax.set_xlim([args[n].histogram.weights.nonzero()[0][0]-2,args[n].histogram.weights.nonzero()[0][-1]+10])
                ax.set_ylim([0,args[n].histogram.weights.max()*1.66])
                if args[n].histogram.dim==1:                 
                   ax.plot(args[n].histogram.weights,ls='steps-mid',label=args[n].histogram.title)
                   ax.legend()
                else:
                   self.plotter.plot2d(args[n])
        plt.tight_layout()


    def rebin(self,hisd):
        """
        Rebin the last 1 or 2d histogram as specified by hisd. Can be used after a "zoom" or "pan" in the canvas. TBD automated rebinning for faster displaying.
        """
        tup = (plt.xlim(),plt.ylim())
        if hisd==1:
           temp=self.plots[-1]
           tup=self.fence(temp,tup,hisd)
           self.dl(tup[0][0],tup[0][1])
        else:
           temp=self.maps[-1]
           tup=self.fence(temp,tup,hisd)
           self.dd(-1,xc=tup[0],yc=tup[1])


    def fence(self,args,tup,hisd):
        """
        Confine plot regions to the bounds of the data in the histogram.
        """
        x0 = tup[0][0]
        x1 = tup[0][1]
        y0 = tup[1][0]
        y1 = tup[1][1]
        if x0 <= args.histogram.x_axis[0]:
           x0 = args.histogram.x_axis[0]
        if x1 >= args.histogram.x_axis[-1]:
           x1 = args.histogram.x_axis[-1]
        if y0 <= 0:
           y0 = 0
        if hisd==1:
           if y1 >= args.histogram.weights.max():
              y1 = args.histogram.weights.max()*1.33
        else:
           if y1 >= args.histogram.y_axis[-1]:
              y1 = args.histogram.y_axis[-1]
        return ((x0,x1),(y0,y1))
   

    def scatterplot(self,plotobj):
        """
           A routine for replotting sparse data as a scatter plot. 
           plotobj should be a maps[register_number] not a dd(histogram_number) instance.
        """
        plt.subplots(1)
        arr=plotobj.histogram.weights.nonzero()
        y=arr[1]
        x=arr[0]
        cc = plotobj.histogram.weights[x,y]

        plt.scatter(x,y,s=50, c=cc, label=plotobj.histogram.title)
        plt.legend()
        plt.show()


    def gui_cal(self,numpeaks,hisnum, ouf='gui_cal.txt',plotrange=None, plotoptions=None):
        """
          A routine for graphically selecting peaks for 
          fitting peaks to calibration data. 
              numpeaks = the number of peaks, 
              ouf = is a file for output
              hisnum = histogram number or range of histograms 
              for fitting, range specified as an integer, 'start-stop',
              (start,stop), or (start,stop,increment).
          This routine takes the fitting range from the current 
          display range as a default. 
        """
        from pylab import ginput
        if plotrange==None:
           x0=0
           x1=16383
        elif len(plotrange)==2:
           try:
              [x0,x1]=plotrange
           except IndexError:
               print("Error in fitting range specification")
   
        if isinstance(hisnum,int):
           st=hisnum
           sp=hisnum
           icr=1
        elif isinstance(hisnum,str):
           st=eval(hisnum.split('-')[0])
           sp=eval(hisnum.split('-')[-1])
           icr=1
        elif len(hisnum)==2:
           st=hisnum[0]
           sp=hisnum[1]
           icr=1
        else:
           try:
               st=hisnum[0]
               sp=hisnum[1]
               icr=hisnum[2]
           except IndexError:
               print("Error in histogram (or range) specification")
        
        #file=open(ouf, 'a')
        #file.write('Position        Error        Area  \n')
        #file.close()

        for i in range(st,sp+icr,icr):
            linein=0
            lin=0
            self.dl(x0,x1)
            while lin<1:
                self.d(i)
                if i == st:
                    print('click '+str(numpeaks)+' peaks for calibration')
                    L=ginput(numpeaks)
#                    maxpk=self.plots[-1].histogram.weights.argmax()
                else:
                    same=input("take same peaks again? enter to accept")
                    if same == "":
                       print("ok")
                    else:
                        print('click '+str(numpeaks)+' peaks for calibration')
                        L=ginput(numpeaks)              
                self.pzot()
                p=[]
#                dp=[]
                for l in range(0,numpeaks):
                    p.append(L[l][0])
#                    if i != numpeaks:
#                        dp.append(L[l+1][0]-L[l][0])
                    self.pk(p[l])
                p=self.fit_peaks(i,(min(p)-100,max(p)+100))
                lin=eval(input("Press 0 to retry, 1 to continue, 2 to take peaks as clicked:"))
                if lin==1:
                    for k in p:
                       file = open(ouf, 'a')   
                       file.write(repr(i)+'   '+repr(k[1]) +'   '+repr(k[2]) +'   '+ repr(k[3])+ '\n')
                       file.close()
                if lin==2:
                    for m in L:
                       file = open(ouf, 'a')   
                       file.write(repr(i)+'   '+ repr(m[0]) +'   '+repr(1) +'   '+ repr(2)+ '\n')
                       file.close()

    def FitCalFile(self, callist,model='linear', inf='gui_cal.txt',ouf='gui_cal_fit.txt'):
        """
            A routine for fitting the file generated from gui_cal or an input file(inf=).
            An input file should contain at least 2 columns, the first is a histogram number or id,
            the second is a channel number corresponding to an energy in callist,
            and the third column may be used to include errors in the fit.
            
            callist is expected to be a 1 or 2 dimesional array of energies and optionally, errors.

            fitting model default is 'linear' but keywords 'linear', 'quadratic', and 'log' are supported.
            
            Example:
            callist=((E1,Err1),(E2,Err2),(E3,Err3)) or callist=(E1,E2,E3)
            Experiment.FitCalFile(callist,model='log',inf='myfile.txt')

            The routine returns a file (default ouf='gui_cal_fit.txt') with 
         """

        from lmfit import minimize, Parameters, Parameter, report_errors, report_ci, conf_interval
        import numpy as np
        import pylab

        inf=open(inf,'r',newline='')
        inl=np.loadtxt(inf)
        lcl=len(callist)
        tl=[]
        #validate data
        for k in np.arange(0,len(inl),lcl):
            j=0
            for i in range(0,lcl):
                if inl[k][0]!=inl[k+i][0]:
                    j=j+1
            if j>0:
                break
            tl.append(inl[k][0])
        
        k=0
        if len(tl)!=len(set(tl)):
            k=k+1

        if j>0:
            print('Error: mismatch in lengths of callist and input file' + repr(inf.name)+'.')
            return None
        if k>0:
            print('Error: Non-unique identifier found in first column in' + repr(inf.name)+'.')
            return None
            
        print('ok')
        #Add errors if not given.
        for i in range(0,len(inl)):
            if len(inl[i])<2:
                np.append(inl[i],0)
        
        for i in range(0,lcl):
            if not hasattr(callist[i],'count'):
                callist[i]=(callist[i],0)

        # define objective function: returns the array to be minimized
        if model=='linear':    
            if len(callist)<2:
               print('Model is incompatible with data given.')
               return None    
            def fcn4fit(params, x, data,weight):
                """ model line, subtract data, weight with error."""
                slope = params['slope'].value
                offset = params['offset'].value
                model = slope*x + offset
                resids = (model - data)
                return np.sqrt(resids**2/weight)
            params = Parameters()
            params.add('slope',  value= 5.0)
            params.add('offset', value= 0.1)

        elif model=='quadratic':
            if len(callist)<3:
               print('Model is incompatible with data given.')
               return None    
            def fcn4fit(params, x, data,weight):
                """ model line, subtract data, weight with error."""
                curve = params['curve'].value
                slope = params['slope'].value
                offset = params['offset'].value
                model = slope*x + offset
                resids = (model - data)
                return np.sqrt(resids**2/weight)
            params = Parameters()
            params.add('curve',  value= 0.2)
            params.add('slope',  value= 5.0)
            params.add('offset', value= 0.1)

        elif model=='log':
            if len(callist)<2:
               print('Model is incompatible with data given.')
               return None
            def fcn4fit(params, x, data,weight):
                """ model line, subtract data, weight with error."""
                amp = params['amp'].value
                const = params['const'].value
                model = const*x+amp
                resids = (model - np.log(data))
                return np.sqrt(resids**2/weight)
            params = Parameters()
            params.add('amp',   value= 5)
            params.add('const', value= .0001)
        
        else:
            print('Model not valid')
            return None
        id=[]
        x1=[]
        y1=[]
        xerr1=[]
        yerr1=[]
        for i in np.arange(0,len(inl),lcl):
            id.append(inl[i][0])
            for j in range(0,lcl):
                if j==0:
                   x1.append([inl[i][1]])
                   xerr1.append([inl[i][2]])
                else:
                   t=x1[-1]
                   t.append(inl[i+j][1])
                   t=xerr1[-1]
                   t.append(inl[i+j][2])
        
        for i in callist:
            y1.append(i[0])
            yerr1.append(i[1])
            
        weight1=[]
        for i in range(0,len(x1)):
            weight1.append([])
            t=weight1[-1]
            for j in range(0,len(x1[i])):
                if xerr1[i][j]==0:
                    xerr1[i][j]=0.0001
                if yerr1[j]==0:
                    yerr1[j]=0.0001
                t.append(((x1[i][j]/xerr1[i][j])**2+(y1[j]/yerr1[j])**2))
    
        x1=np.array(x1)
        y1=np.array(y1)
        weight1=np.array(weight1)
        pylab.close()
        for i in range(0,len(x1)):
            # do fit, here with leastsq model
            resfit = minimize(fcn4fit, params, args=(x1[i], y1, weight1[i]))
            file=open(ouf,'a')

            # calculate final result
            if model=='linear':
                if i==0:
                    file.write('id , m , b , (for y= m*x+b) \n')
                residfit= resfit.residual
                finalfit = resfit.params['slope'].value*x1[i] + resfit.params['offset'].value
                file.write(repr(id[i]) + ' , ' +repr(resfit.params['slope'].value) +' , '
                    +repr(resfit.params['offset'].value) + '\n')
            if model=='quadratic':
                if i==0:
                    file.write('id , a , b , c (for y= a*x**2+b*x+c) \n')
                residfit= resfit.residual
                finalfit = resfit.params['curve'].value*x1[i]**2+resfit.params['slope'].value*x1[i]+resfit.params['offset'].value
                file.write(repr(id[i]) + ' , ' +repr(resfit.params['curve'].value) +' , '
                           +repr(resfit.params['slope'].value) +' , '
                           +repr(resfit.params['offset'].value) + '\n')
            if model=='log':
                if i==0:
                    file.write('id , A , t (for y=A*exp(t*x)) \n')
                residfit= resfit.residual
                finalfit = np.exp(resfit.params['amp'].value)*np.exp(x1[i] * resfit.params['const'].value)
                file.write(repr(id[i]) + ' , ' +repr(resfit.params['amp'].value) +' , '
                           +repr(resfit.params['const'].value) + '\n')
            
            file.close

            # write error report
            report_errors(params)

            # try to plot results
            if resfit.success:
                pylab.errorbar(x1[i],y1,xerr=xerr1[i],yerr=yerr1)
                pylab.plot(x1[i], y1, '+')
                pylab.plot(x1[i], finalfit)
                pylab.show()
                print('Res(x)= ' + repr(residfit))
            else:
                print('fit failed')
        


    def WaveletCal(self,his,rx=None):
        """ 
            Method for using Root "wavelet" methods for automatic Calibrations.
            See cal_parms_def for the default parameters you will likely find 
            cal_parms[0].energies will need to be changed to suit your needs or set your own 
            defaults using the SaveParms method. cal_parms[1] are peak search parameters and cal_parms[2] are the filter
            filter parameters. 
        """
        import Pyspectr.calibration as cal
        import numpy as np
        from time import time
        
        if rx==None:
            try:
                xmin, xmax = self.xlim
            except:
                xmin, xmax = (10,8000)
        else:
            xmin,xmax = rx

        self.d(his)
        hist=self.plots[-1]
        xdat=hist.histogram.x_axis
        wdat=hist.histogram.weights

        xpeaks, coef = cal.calibrate_area(xdat,wdat,xmin,xmax,self.cal_parms, self.fltr_parms, self.srch_parms)
        report, control = cal.make_report(xpeaks,coef,1,self.cal_parms.energies,filename = None) 

        print(report,control)
        return(coef)

    def SaveParms(self, funcname='set_default_new',force='n'):

        """ 
            Once functional parameters have been found. Use this method to append 
            them to the cal_parms_def.py in the source directory. It will make a 
            local copy unless directed otherwise. Once build and install has been run
            It can then be used at any time via via the commmand 
            " from Pyspectr.cal_parms_def import funcname " funcname = " 
            set_default is protected unless force='y'.
            then run e.cal_parms, e.fltr_parms, e.srch_parms = set_default_Cd()
        """
        ouf=open('cal_parms_def.py','a')

        ouf.write("def " +funcname+ "():"+"\n" +
                  "    import numpy as np \n    #Define Parameter Classes \n" +
                  "    class record: pass\n")

        ouf.write("    #energy scale check-calibration\n    calibration_properties = record()\n"+
                  "    calibration_properties.visualize="+repr(self.cal_parms.visualize)+"\n"
                  "    calibration_properties.sigma="+repr(self.cal_parms.sigma)+"\n"
                  "    calibration_properties.dlt="+repr(self.cal_parms.dlt)+"\n"
                  "    calibration_properties.energies="+repr(self.cal_parms.energies)+"\n\n"
                  "    search_properties = record()               #for smooth spectrum\n"
                  "    search_properties.widths=np.arange(1,10)\n"
                  "    search_properties.wavelet="+ repr(self.srch_parms.wavelet)+"\n"
                  "    search_properties.min_length="+ repr(self.srch_parms.min_length)+"\n"
                  "    search_properties.min_snr="+ repr(self.srch_parms.min_snr)+"\n"
                  "    search_properties.noise_perc="+ repr(self.srch_parms.noise_perc)+"\n\n"
                  "    filter_properties = record()\n"
                  "    filter_properties.window_smooth="+ repr(self.fltr_parms.window_smooth)+"\n"
                  "    filter_properties.smooth_wavelet="+ repr(self.fltr_parms.smooth_wavelet)+"\n"
                  "    filter_properties.background_options="+ repr(self.fltr_parms.background_options)+"\n"
                  "    filter_properties.threshold = "+ repr(self.fltr_parms.threshold)+"\n")
        ouf.write("    return(calibration_properties,filter_properties,search_properties)")



if __name__ == "__main__":
    pass

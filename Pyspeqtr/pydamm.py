#!/usr/bin/env python3
"""
N.T Brewer 2016
brewer.nathant@gmail.com
Distributed under GNU General Public Licence v3

This module is inteded to be loaded in an interactive interpreter session.
The ipython is strongly recommended. The pydamm is a python replacement for
DAMM programm.

"""

import math
import numpy
import sys
import re

from collections import deque

from Pyspeqtr import hisfile as hisfile
from Pyspeqtr import histogram as histogram
from Pyspeqtr.exceptions import GeneralError as GeneralError
from Pyspeqtr.decay_fitter import DecayFitter as DecayFitter
from Pyspeqtr.peak_fitter import PeakFitter as PeakFitter


"""
from PyQt4.Qt import QImage, QPainter, QBuffer, QIODevice, QByteArray
from pyqtgraph import GraphicsLayoutWidget
from pyqtgraph.opengl import GLViewWidget, GLMeshItem
from pyqtgraph.Qt import QtGui
%gui qt
"""

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


    def __init__(self, file_name, size=11,gui="Online"):
        """Initialize, open data file (his) and open plot window
        (size parameter decides on plot dimensions)

        """
        if gui.lower() == "online" or gui.lower() == "offline":
            self.gui=gui.lower()
        else:
            print("using default: online. options are \"offline\" or \"online\".")
            self.gui="online"

        if self.gui == "online":
            from pyqtgraph import exit
            from IPython import get_ipython
            from Pyspeqtr import plotter_Qt as plotter
            ip3 = get_ipython()
            ip3.magic("gui 'qt'")

        if self.gui == "offline":
            import matplotlib.pyplot as plt
            import matplotlib.cm as cm
            from Pyspeqtr import plotter_mpl as plotter
            
        self.file_name = file_name
        # The current (active) file
        self.hisfile = None
        self.load(file_name)

        # Peaks for fitting
        self.peaks = []
        
        # plotter front-end
        self.plotter = plotter.Plotter(size)
        

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

#        if mode == 2:
#            self.plotter.ylin()

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
###            self.plotter.plot1d(gate_plot, Experiment.xlim, ylim)

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
        #if clear and plot:
        #    self.plotter.clear()

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
###            self.plotter.plot1d(gate_plot, Experiment.xlim, ylim)

        return gate_plot

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

    def d(self, *args, norm=1, bin_size=1, clear=True):
    #def scatterplot(self,plotobj):
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
                """ylim = None
                if self.ylim is None and i_plot == active_plots:
                    ylim = self._auto_scale_y()
                else:
                    ylim = self.ylim
                """
                # Note that ylim is autoscaled above if self.ylim is None
                # But we still keep self.ylim None, 
                # to indicate autoscaling
                self.plotter.plot1d(plot, self.xlim, self.ylim)

        # Return plots that were added or activated
        return plots
    def _auto_scale_x(self):
        """Find the x limits taking into account all active plots """
        xr=self.plotter._widget.plotItem.viewRange()[0]
        return((xr[0],xr[1]))
    def _auto_scale_y(self):
        """Find the y limits taking into account all active plots """
        yr=self.plotter._widget.plotItem.viewRange()[1]
        return((yr[0],yr[1]))

    def dl(self, x0=None, x1=None):
        """Change x range of 1D histogram"""
        if self.mode != 1:
            return None

        if x0 is None and x1 is None:
            Experiment.xlim = None
#            self.plotter.xlim(self._auto_scale_x())
            self.plotter.xlim(None)
        else:
            Experiment.xlim = [x0, x1]
            self.plotter.xlim(Experiment.xlim)
        """
        elif x0 is None:
#            xr=self._auto_scale_x()
#            Experiment.xlim = (xr[0],x1)
            self.plotter.xlim((None,x1))
            Experiment.xlim=self.plotter._widget.plotItem.viewRange()[0]
        elif x1 is None:
#            xr=self._auto_scale_x()
#            Experiment.xlim = (x0,xr[1])
            self.plotter.xlim((x0,None))
            Experiment.xlim=self.plotter._widget.plotItem.viewRange()[0]
        """


#        if Experiment.ylim is None:
#            self.plotter._widget.(self._auto_scale_y())

    def dmm(self, y0=None, y1=None):
        """Change y range of 1D histogram"""
        if self.mode != 1:
            return None
      
        if y0 is None and y1 is None:
            Experiment.ylim = None
            self.plotter.ylim(None)
        else:
            Experiment.ylim = [y0, y1]
            self.plotter.ylim(Experiment.ylim)
        """
        if y0 is None and y1 is None:
            Experiment.ylim = None
            self.plotter.ylim(self._auto_scale_y())
        elif y0 is None:
            yr=self._auto_scale_y()
            Experiment.ylim = (yr[0],y1)
            self.plotter.ylim(Experiment.ylim)
        elif y1 is None:
            yr=self._auto_scale_y()
            Experiment.ylim = (y0,yr[1])
            self.plotter.ylim(Experiment.ylim)
        else:
            Experiment.ylim = (y0, y1)
            self.plotter.ylim(Experiment.ylim)

#        if Experiment.xlim is None:
#            self.plotter.xlim(self._auto_scale_x())
        """

    def log(self):
        if self.mode == 1:
            self.plotter.ylog()
     #   elif self.mode == 2:
     #       Experiment.logz = True
     #       self.dd(-1, xc=Experiment.xlim2d, yc=Experiment.ylim2d)
    def lin(self):
        if self.mode == 1:
            self.plotter.ylin()
     #   if self.mode == 2:
     #       Experiment.logz = False
     #       self.dd(-1, xc=Experiment.xlim2d, yc=Experiment.ylim2d)

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
        #self.plotter.clear()

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

"""




    #def scatterplot(self,plotobj):

    def xc(self, x0=None, x1=None):

    def yc(self, y0=None, y1=None):

    def clear(self):

    def fit_peaks(self, his=None, rx=None, clear=True):

    def fit_decay(self, his, gate, cycle, 
                        t_bin=1, time_range=None,
                        model='grow_decay',
                        pars=None,
                        clear=True):
"""     
if __name__ == "__main__":
    pass

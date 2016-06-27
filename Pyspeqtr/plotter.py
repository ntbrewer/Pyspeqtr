#!/usr/bin/env python3
"""K. Miernik 2012
k.a.miernik@gmail.com
Distributed under GNU General Public Licence v3

This module provides simple front-end to matplotlib

"""

import math
import numpy as np
import copy
import pyqtgraph as pg
from PyQt4.Qt import QImage, QPainter, QBuffer, QIODevice, QByteArray
from pyqtgraph import GraphicsLayoutWidget, PlotWidget
from pyqtgraph.opengl import GLViewWidget, GLMeshItem
from pyqtgraph.Qt import QtGui, QtCore


class Plotter:
    """ This class communicates with the matplotlib library
    and plot the data

    """

    def __init__(self, size):
        """Initialize the plot window, size defines the shape and size
        of the figure
        0 - None, 
        1 - 8x6, 
        11 (default) - 12x8,
        2 - 2 figs 8x8,
        12 - 2 figs 12x8

        """
#        print('ok')
        class CustomViewBox(pg.ViewBox):
           def __init__(self, *args, **kwds):
               pg.ViewBox.__init__(self, *args, **kwds)
               self.setMouseMode(self.RectMode)
        
            ## reimplement right-click to zoom out
           def mouseClickEvent(self, ev):
               if ev.button() == QtCore.Qt.RightButton and QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                   self.autoRange()
               else:
                   pg.ViewBox.mouseClickEvent(self, ev)
            
           def mouseDragEvent(self, ev):
               mod = QtGui.QApplication.keyboardModifiers()
               if mod == QtCore.Qt.ControlModifier:
                   self.setMouseMode(self.PanMode)
               else:
                   self.setMouseMode(self.RectMode)
               
               if ev.button() == QtCore.Qt.RightButton:
                   pg.ViewBox.mouseDragEvent(self, ev)
               else:
                   pg.ViewBox.mouseDragEvent(self, ev)

        self.vb = CustomViewBox()
# not necessary?    def clear(self):

    def ylog(self):
        """
           Change y-scale to log
        """
        vr = self._widget.plotItem.viewRange()
        xr = vr[0]
        yr = vr[1]
#        print(xr,yr)
        if yr[0]<=0:
            yr[0]=1
        yr[0]=np.log10(yr[0])
        yr[1]=np.log10(yr[1])
#        print(xr,yr)
        self._widget.plotItem.disableAutoRange()
        self._widget.plotItem.setLogMode(x=False, y=True)
        self._widget.plotItem.setRange(xRange=xr, yRange=yr,padding=None)

    def ylin(self):
        """
           Change y-scale to linear
        """
        vr = self._widget.plotItem.viewRange()
        xr = vr[0]
        yr = vr[1]
#        print(xr,yr)
        yr[0]=10**yr[0]
        yr[1]=10**yr[1]
#        print(xr,yr)
        self._widget.plotItem.disableAutoRange()
        self._widget.plotItem.setLogMode(x=False, y=False)
        self._widget.plotItem.setRange(xRange=xr, yRange=yr,padding=None)

    def plot1d(self, plot, xlim=None, ylim=None):
        """ Plot 1D histogram
            The mode defines the way the data are presented,
            'histogram' is displayed with steps
            'function' with continuus line
            'errorbar' with yerrorbars

            The norm (normalization factor) and bin_size are given
            for the display purposes only. The histogram is not altered.

        """

        histo = plot.histogram

        self._widget = PlotWidget(viewBox=self.vb)
        plt = self._widget.plotItem
        plt.setTitle(histo.title)
        plt.plot(histo.x_axis, histo.weights)
        self._widget.show()
#        if plot.mode == 'histogram':
#            current_plot.plot()
#        win.show()
#        app.processEvents()

    def _repr_png_(self):

        self._widget.hide()

        QtGui.QApplication.processEvents()
        
        try:
            self.image = QImage(self._widget.viewRect().size().toSize(),
                                QImage.Format_RGB32)
        except AttributeError:
            self._widget.updateGL()
            self.image = self._widget.grabFrameBuffer()
            
        painter = QPainter(self.image)
        self._widget.render(painter)

        byte_array = QByteArray()
        buffer = QBuffer(byte_array)
        buffer.open(QIODevice.ReadWrite)
        self.image.save(buffer, 'PNG')
        buffer.close()        

        return bytes(byte_array)

    def xlim(self, x_range):
        """
        Set x range of plot preserving y limits.
        """

        if x_range==None:
            yar=self._widget.plotItem.getViewBox().autoRangeEnabled()[1]
            if yar==False:
                yr=self._widget.viewRange()[1]
            else:
                yr=None
#            self._widget.plotItem.autoRange()
            self._widget.plotItem.enableAutoRange()
            self._widget.plotItem.setAutoVisible(x=True,y=True)
            if yr!=None:
                self._widget.plotItem.setRange(yRange=yr,padding=None)
            return None
        else:
            if x_range[0]==None:
               x_range[0]=0
            if x_range[1]==None:
               self._widget.plotItem.setAutoVisible(x=True,y=False)
               xr=self._widget.viewRange()[0]
               x_range[1]=xr[1]            
            self._widget.plotItem.setXRange(x_range[0],x_range[1])
            self._widget.plotItem.setAutoVisible(x=False,y=True)

    def ylim(self, y_range):
        """
        Set y range of plot preserving x limits.
        """
        if y_range==None:
            xar=self._widget.plotItem.getViewBox().autoRangeEnabled()[0]
            if xar==False:
                xr=self._widget.viewRange()[0]
            else:
                xr=None
            self._widget.plotItem.enableAutoRange()
            self._widget.plotItem.setAutoVisible(x=True,y=True)
            if xr!=None:
                self._widget.plotItem.setRange(xRange=xr,padding=None)

            return None
        else:
            if y_range[0]==None:
               y_range[0]=0
            if y_range[1]==None:
               self._widget.plotItem.setAutoVisible(y=True,x=False)
               yr=self._widget.viewRange()[1]
               y_range[1]=yr[1]            
            self._widget.plotItem.setYRange(y_range[0],y_range[1])
            self._widget.plotItem.setAutoVisible(y=False,x=True)


#    def plot1d_4panel(self, plot, ranges):
        """
        Special 1D histogram plot. The plot is broken into 4 panels (stacked verically)
        the ranges variable should be given in a (x0, x1, x2, x3, x4) format, where
        xi defines the ranges of the subplots (x0-x1, x1-x2, x2-x3, x3-x4)

        """
        
    def plot2d(self, plot, xc=None, yc=None, logz=False):
        """Plot 2D histogram 
        xc is x range, yc is y range 
   	"""     


        if plot.histogram.dim != 2:
            raise GeneralError('plot2d function needs a 2D histogram!')
        x = plot.histogram.x_axis
        y = plot.histogram.y_axis
        w = plot.histogram.weights

        if xc is not None:
            x = x[xc[0]:xc[1]]
            w = w[xc[0]:xc[1],:]

        if yc is not None:
            y = y[yc[0]:yc[1]]
            w = w[:, yc[0]:yc[1]]

        title = plot.histogram.title
        # If logaritmic scale is used, mask values <= 0
        if logz:
            w = numpy.ma.masked_where(w <= 0, numpy.log10(w))
            title += ' (log10)'

        self._widget = pg.ImageView(view=pg.PlotItem(title=title))
        gv = self._widget.getView()
        gv.invertY(False)
        self._widget.setImage(w,pos=[x[0]-0.5,y[0]-0.5])
        self._widget.show()



#!/usr/bin/env python
#########################################################################################
#
# Visualizer for MRI volumes
#
#
# ---------------------------------------------------------------------------------------
# Copyright (c) 2015 Polytechnique Montreal <www.neuro.polymtl.ca>
# Authors: Benjamin De Leener
# Created: 2015-01-30
# Modified: 2015-02-02
#
# About the license: see the file LICENSE.TXT
#########################################################################################
import sys
from msct_parser import Parser
from msct_image import Image

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


class TrioPlot:
    """
    This class manages mouse events on the three image subplots.
    """
    def __init__(self, viewer_parent, ax_axial, fig_axial, ax_frontal, fig_frontal, ax_sagittal, fig_sagittal, volume):
        self.viewer = viewer_parent
        self.ax_axial = ax_axial
        self.fig_axial = fig_axial
        self.ax_frontal = ax_frontal
        self.fig_frontal = fig_frontal
        self.ax_sagittal = ax_sagittal
        self.fig_sagittal = fig_sagittal
        self.volume = volume
        self.press = 0, 0

        self.cidpress_axial = None
        self.cidrelease_axial = None
        self.cidmotion_axial = None
        self.cidpress_frontal = None
        self.cidrelease_frontal = None
        self.cidmotion_frontal = None
        self.cidpress_sagittal = None
        self.cidrelease_sagittal = None
        self.cidmotion_sagittal = None

    def connect(self):
        """
        connect to all the events we need
        """
        self.cidpress_axial = self.fig_axial.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease_axial = self.fig_axial.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion_axial = self.fig_axial.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

        self.cidpress_frontal = self.fig_frontal.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease_frontal = self.fig_frontal.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion_frontal = self.fig_frontal.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

        self.cidpress_sagittal = self.fig_sagittal.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease_sagittal = self.fig_sagittal.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion_sagittal = self.fig_sagittal.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        self.press = event.xdata, event.ydata
        return

    def move(self, event):
        if event.xdata is None or event.ydata is None or self.press[0] is None or self.press[1] is None:
            return

        if event.xdata and abs(event.xdata - self.press[0]) < 1 and abs(event.ydata - self.press[0]) < 1:
            self.press = event.xdata, event.ydata
            return

        if event.inaxes == self.fig_axial.axes:
            self.fig_frontal.set_data(self.volume.data[event.xdata, :, :].T)
            self.fig_sagittal.set_data(self.volume.data[:, event.ydata, :].T)

            self.ax_frontal.draw_artist(self.fig_frontal)
            self.fig_frontal.figure.canvas.blit(self.ax_frontal.bbox)

            self.ax_sagittal.draw_artist(self.fig_sagittal)
            self.fig_sagittal.figure.canvas.blit(self.ax_sagittal.bbox)

            self.press = event.xdata, event.ydata
            return

        elif event.inaxes == self.fig_frontal.axes:
            self.fig_axial.set_data(self.volume.data[:, :, event.ydata].T)
            self.fig_sagittal.set_data(self.volume.data[:, event.xdata, :].T)

            self.ax_axial.draw_artist(self.fig_axial)
            self.fig_axial.figure.canvas.blit(self.ax_axial.bbox)

            self.ax_sagittal.draw_artist(self.fig_sagittal)
            self.fig_sagittal.figure.canvas.blit(self.ax_sagittal.bbox)

            self.press = event.xdata, event.ydata
            return

        elif event.inaxes == self.fig_sagittal.axes:
            self.fig_axial.set_data(self.volume.data[:, :, event.ydata].T)
            self.fig_frontal.set_data(self.volume.data[event.xdata, :, :].T)

            self.ax_axial.draw_artist(self.fig_axial)
            self.fig_axial.figure.canvas.blit(self.ax_axial.bbox)

            self.ax_frontal.draw_artist(self.fig_frontal)
            self.fig_frontal.figure.canvas.blit(self.ax_frontal.bbox)

            self.press = event.xdata, event.ydata
            return

        else:
            self.press = event.xdata, event.ydata
            return

    def on_motion(self, event):
        if event.button == 1:
            return self.move(event)
        else:
            return

    def on_release(self, event):
        if event.button == 1:
            return self.move(event)
        else:
            return

    def update_min_max_contrast(self, min_value, max_value):
        self.fig_axial.set_clim(min_value, max_value)
        self.fig_frontal.set_clim(min_value, max_value)
        self.fig_sagittal.set_clim(min_value, max_value)

        self.ax_axial.draw_artist(self.fig_axial)
        self.fig_axial.figure.canvas.blit(self.ax_axial.bbox)
        self.ax_frontal.draw_artist(self.fig_frontal)
        self.fig_frontal.figure.canvas.blit(self.ax_frontal.bbox)
        self.ax_sagittal.draw_artist(self.fig_sagittal)
        self.fig_sagittal.figure.canvas.blit(self.ax_sagittal.bbox)

    def disconnect(self):
        """
        disconnect all the stored connection ids
        """
        self.fig_axial.figure.canvas.mpl_disconnect(self.cidpress_axial)
        self.fig_axial.figure.canvas.mpl_disconnect(self.cidrelease_axial)
        self.fig_axial.figure.canvas.mpl_disconnect(self.cidmotion_axial)

        self.fig_frontal.figure.canvas.mpl_disconnect(self.cidpress_frontal)
        self.fig_frontal.figure.canvas.mpl_disconnect(self.cidrelease_frontal)
        self.fig_frontal.figure.canvas.mpl_disconnect(self.cidmotion_frontal)

        self.fig_sagittal.figure.canvas.mpl_disconnect(self.cidpress_sagittal)
        self.fig_sagittal.figure.canvas.mpl_disconnect(self.cidrelease_sagittal)
        self.fig_sagittal.figure.canvas.mpl_disconnect(self.cidmotion_sagittal)


class VolViewer(object):
    """
    This class is a visualizer for volumes (3D images).
    """
    def __init__(self, image_input):
        if isinstance(image_input, Image):
            from numpy import percentile
            self.image = image_input
            self.image.change_orientation('LAS')
            self.min_contrast = percentile(self.image.data[:], 1)
            self.max_contrast = percentile(self.image.data[:], 99)
            self.min_contrast_init = self.min_contrast  # alternative: min(self.image.data[:])
            self.max_contrast_init = self.max_contrast  # alternative: max(self.image.data[:])
            self.trio = None
            self.fig = None
        else:
            print "Error, the image is actually not an image"

    def update_min_contrast(self, val):
        self.min_contrast = val
        self.trio.update_min_max_contrast(self.min_contrast, self.max_contrast)

    def update_max_contrast(self, val):
        self.max_contrast = val
        self.trio.update_min_max_contrast(self.min_contrast, self.max_contrast)

    def show(self):
        self.fig = plt.figure(facecolor='black')
        self.fig.subplots_adjust(bottom=0.1, left=0.1)

        im_size = self.image.data.shape

        ax_axial = self.fig.add_subplot(221)
        im_plot_axial = ax_axial.imshow(self.image.data[:, :, int(im_size[2]/2)].T, vmin=self.min_contrast,
                                        vmax=self.max_contrast)
        im_plot_axial.set_cmap('gray')
        im_plot_axial.set_interpolation('nearest')

        ax_frontal = self.fig.add_subplot(222)
        im_plot_frontal = ax_frontal.imshow(self.image.data[int(im_size[0]/2), :, :].T, vmin=self.min_contrast,
                                            vmax=self.max_contrast)
        im_plot_frontal.set_cmap('gray')
        im_plot_frontal.set_interpolation('nearest')

        ax_sagittal = self.fig.add_subplot(223)
        im_plot_sagittal = ax_sagittal.imshow(self.image.data[:, int(im_size[1]/2), :].T, vmin=self.min_contrast,
                                              vmax=self.max_contrast)
        im_plot_sagittal.set_cmap('gray')
        im_plot_sagittal.set_interpolation('nearest')

        self.trio = TrioPlot(self, ax_axial, im_plot_axial, ax_frontal, im_plot_frontal, ax_sagittal, im_plot_sagittal,
                             self.image)
        self.trio.connect()

        from matplotlib.widgets import Slider
        position_slider_min = self.fig.add_axes([0.55, 0.1, 0.35, 0.03])
        slider_min = Slider(position_slider_min, 'Min', self.min_contrast_init, self.max_contrast_init,
                            valinit=self.min_contrast)
        slider_min.on_changed(self.update_min_contrast)
        slider_min.label.set_color('white')
        slider_min.valtext.set_color('white')
        position_slider_max = self.fig.add_axes([0.55, 0.05, 0.35, 0.03])
        slider_max = Slider(position_slider_max, 'Max', self.min_contrast_init, self.max_contrast_init,
                            valinit=self.max_contrast, slidermin=slider_min)
        slider_max.on_changed(self.update_max_contrast)
        slider_max.label.set_color('white')
        slider_max.valtext.set_color('white')
        slider_min.__setattr__('slidermax', slider_max)

        plt.show()

# ======================================================================================================================
# Start program
# ======================================================================================================================
if __name__ == "__main__":
    parser = Parser(__file__)
    parser.usage.set_description('Volume Viewer')
    parser.add_option("-i", "file", "file", True)
    arguments = parser.parse(sys.argv[1:])

    image = Image(arguments["-i"])
    viewer = VolViewer(image)
    viewer.show()

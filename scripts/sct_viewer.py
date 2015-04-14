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
        self.axial_xlim_origin = self.ax_axial.get_xlim()
        self.axial_ylim_origin = self.ax_axial.get_ylim()
        self.fig_axial = fig_axial
        self.ax_frontal = ax_frontal
        self.frontal_xlim_origin = self.ax_frontal.get_xlim()
        self.frontal_ylim_origin = self.ax_frontal.get_ylim()
        self.fig_frontal = fig_frontal
        self.ax_sagittal = ax_sagittal
        self.sagittal_xlim_origin = self.ax_sagittal.get_xlim()
        self.sagittal_ylim_origin = self.ax_sagittal.get_ylim()
        self.fig_sagittal = fig_sagittal
        self.volume = volume
        self.press = 0, 0
        self.position = 0, 0

        self.cidpress_axial = None
        self.cidrelease_axial = None
        self.cidmotion_axial = None
        self.cidpress_frontal = None
        self.cidrelease_frontal = None
        self.cidmotion_frontal = None
        self.cidpress_sagittal = None
        self.cidrelease_sagittal = None
        self.cidmotion_sagittal = None

    def draw_axial(self):
        self.ax_axial.draw_artist(self.fig_axial)
        self.fig_axial.figure.canvas.blit(self.ax_axial.bbox)

    def draw_frontal(self):
        self.ax_frontal.draw_artist(self.fig_frontal)
        self.fig_frontal.figure.canvas.blit(self.ax_frontal.bbox)

    def draw_sagittal(self):
        self.ax_sagittal.draw_artist(self.fig_sagittal)
        self.fig_sagittal.figure.canvas.blit(self.ax_sagittal.bbox)

    def connect(self):
        """
        connect to all the events we need
        """
        self.cidpress_axial = self.fig_axial.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease_axial = self.fig_axial.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion_axial = self.fig_axial.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.cidscroll_axial = self.fig_axial.figure.canvas.mpl_connect('scroll_event', self.on_scroll)

        self.cidpress_frontal = self.fig_frontal.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease_frontal = self.fig_frontal.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion_frontal = self.fig_frontal.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.cidscroll_frontal = self.fig_frontal.figure.canvas.mpl_connect('scroll_event', self.on_scroll)

        self.cidpress_sagittal = self.fig_sagittal.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease_sagittal = self.fig_sagittal.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion_sagittal = self.fig_sagittal.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.cidscroll_sagittal = self.fig_sagittal.figure.canvas.mpl_connect('scroll_event', self.on_scroll)

    def on_press(self, event):
        self.press = event.xdata, event.ydata
        return

    def move(self, event):
        if event.inaxes == self.fig_axial.axes:
            self.fig_frontal.set_data(self.volume.data[event.xdata, :, :].T)
            self.fig_sagittal.set_data(self.volume.data[:, event.ydata, :].T)

            self.draw_frontal()
            self.draw_sagittal()

            self.press = event.xdata, event.ydata
            return

        elif event.inaxes == self.fig_frontal.axes:
            self.fig_axial.set_data(self.volume.data[:, :, event.ydata].T)
            self.fig_sagittal.set_data(self.volume.data[:, event.xdata, :].T)

            self.draw_axial()
            self.draw_sagittal()

            self.press = event.xdata, event.ydata
            return

        elif event.inaxes == self.fig_sagittal.axes:
            self.fig_axial.set_data(self.volume.data[:, :, event.ydata].T)
            self.fig_frontal.set_data(self.volume.data[event.xdata, :, :].T)

            self.draw_axial()
            self.draw_frontal()

            self.press = event.xdata, event.ydata
            return

        else:
            self.press = event.xdata, event.ydata
            return

    def on_scroll(self, event):
        if event.inaxes not in [self.fig_axial.axes, self.fig_frontal.axes, self.fig_sagittal.axes]:
            return

        base_scale = 0.95
        # get the current x and y limits
        cur_xlim = event.inaxes.get_xlim()
        cur_ylim = event.inaxes.get_ylim()
        if event.button == 'up':
            # deal with zoom in
            scale_factor = 1 / base_scale
        elif event.button == 'down':
            # deal with zoom out
            scale_factor = base_scale
        else:
            # deal with something that should never happen
            scale_factor = 1

        new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
        new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor
        relative_x = (cur_xlim[1] - event.xdata) / (cur_xlim[1] - cur_xlim[0])
        relative_y = (cur_ylim[1] - event.ydata) / (cur_ylim[1] - cur_ylim[0])
        event.inaxes.set_xlim([event.xdata - new_width * (1 - relative_x), event.xdata + new_width * relative_x])
        event.inaxes.set_ylim([event.ydata - new_height * (1 - relative_y), event.ydata + new_height * relative_y])

        if event.inaxes == self.fig_axial.axes:
            self.ax_axial.draw_artist(self.ax_axial.patch)
            self.draw_axial()
        elif event.inaxes == self.fig_frontal.axes:
            self.ax_frontal.draw_artist(self.ax_frontal.patch)
            self.draw_frontal()
        elif event.inaxes == self.fig_sagittal.axes:
            self.ax_sagittal.draw_artist(self.ax_sagittal.patch)
            self.draw_sagittal()

        self.position = event.xdata, event.ydata

        return

    def on_motion(self, event):
        """
        This function controls the movement of the mouse over panels.
        If the left button is pressed, the figures will display the slices according to the position of the mouse.
        If the right button is pressed, the mouse movements control the zoom on the corresponding panel.
        :return
        """
        if event.xdata is None or event.ydata is None or self.press[0] is None or self.press[1] is None:
            return

        if event.xdata and abs(event.xdata - self.press[0]) < 1 and abs(event.ydata - self.press[0]) < 1:
            self.press = event.xdata, event.ydata
            return

        if event.button == 1:  # left button
            return self.move(event)
        elif event.button == 3:  # right button
            return self.zoom_pan(event)
        else:
            return

    def on_release(self, event):
        if event.xdata is None or event.ydata is None or self.press[0] is None or self.press[1] is None:
            return

        if event.xdata and abs(event.xdata - self.press[0]) < 1 and abs(event.ydata - self.press[0]) < 1:
            self.press = event.xdata, event.ydata
            return

        if event.button == 1:
            return self.move(event)
        elif event.button == 3:  # right button
            return self.zoom_pan(event)
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
        self.fig_axial.figure.canvas.mpl_disconnect(self.cidscroll_axial)

        self.fig_frontal.figure.canvas.mpl_disconnect(self.cidpress_frontal)
        self.fig_frontal.figure.canvas.mpl_disconnect(self.cidrelease_frontal)
        self.fig_frontal.figure.canvas.mpl_disconnect(self.cidmotion_frontal)
        self.fig_frontal.figure.canvas.mpl_disconnect(self.cidscroll_frontal)

        self.fig_sagittal.figure.canvas.mpl_disconnect(self.cidpress_sagittal)
        self.fig_sagittal.figure.canvas.mpl_disconnect(self.cidrelease_sagittal)
        self.fig_sagittal.figure.canvas.mpl_disconnect(self.cidmotion_sagittal)
        self.fig_sagittal.figure.canvas.mpl_disconnect(self.cidscroll_sagittal)


class VolViewer(object):
    """
    This class is a visualizer for volumes (3D images).
    """
    def __init__(self):
        """
        Constructor
        :param image_input:
        :return:
        """
        self.list_image = []
        self.trio = None
        self.fig = None

    def add_image(self, image_input):
        if isinstance(image_input, Image):
            self.list_image.append(ImageViewer(image_input))
        else:
            print "Error, the image is actually not an image"

    def update_min_contrast(self, val):
        """
        Test
        :param val:
        :return
        """
        self.min_contrast = val  # type: int
        self.trio.update_min_max_contrast(self.min_contrast, self.max_contrast)

    def update_max_contrast(self, val):
        self.max_contrast = val
        self.trio.update_min_max_contrast(self.min_contrast, self.max_contrast)

    def show(self):
        self.fig = plt.figure(facecolor='black')
        self.fig.subplots_adjust(bottom=0.1, left=0.1)

        ax_axial = self.fig.add_subplot(221)
        ax_axial.patch.set_facecolor('black')
        for image_itr in self.list_image:
            im_plot_axial = ax_axial.imshow(image_itr.data[:, :, int(image_itr.data.shape[2]/2)].T,
                                            vmin=image_itr.min_contrast, vmax=image_itr.max_contrast)
            im_plot_axial.set_cmap('gray')
            im_plot_axial.set_interpolation('nearest')

        ax_frontal = self.fig.add_subplot(222)
        ax_frontal.patch.set_facecolor('black')
        for image_itr in self.list_image:
            im_plot_frontal = ax_frontal.imshow(image_itr.data[int(image_itr.data.shape[0]/2), :, :].T,
                                                vmin=image_itr.min_contrast, vmax=image_itr.max_contrast)
            im_plot_frontal.set_cmap('gray')
            im_plot_frontal.set_interpolation('nearest')

        ax_sagittal = self.fig.add_subplot(223)
        ax_sagittal.patch.set_facecolor('black')
        for image_itr in self.list_image:
            im_plot_sagittal = ax_sagittal.imshow(image_itr.data[:, int(image_itr.data.shape[1]/2), :].T,
                                                  vmin=image_itr.min_contrast, vmax=image_itr.max_contrast)
            im_plot_sagittal.set_cmap('gray')
            im_plot_sagittal.set_interpolation('nearest')

        # TODO: change TrioPlot to consider multiple image
        # TODO: make alpha possible
        # TODO: make masking possible when data is zero (zero is transparent)
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


class ImageViewer(Image):
    def __init__(self, image):
        super(ImageViewer, self).__init__(image)
        from numpy import percentile
        self.min_contrast = percentile(self.data[:], 1)
        self.max_contrast = percentile(self.data[:], 99)
        self.min_contrast_init = self.min_contrast  # alternative: min(self.image.data[:])
        self.max_contrast_init = self.max_contrast  # alternative: max(self.image.data[:])
        self.change_orientation('LAS')

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

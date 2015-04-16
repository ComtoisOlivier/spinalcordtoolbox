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
    def __init__(self, viewer_parent):
        self.viewer = viewer_parent
        self.press = 0, 0
        self.position = 0, 0

    def on_press(self, event):
        self.press = event.xdata, event.ydata
        return

    def move(self, event):
        if event.inaxes == self.viewer.list_image[0].axial_plot.axes:
            for image_itr in self.viewer.list_image:
                image_itr.frontal_plot.set_data(image_itr.data[event.xdata, :, :].T)
                image_itr.sagittal_plot.set_data(image_itr.data[:, event.ydata, :].T)

                image_itr.draw_frontal(self.viewer.ax_frontal)
                image_itr.draw_sagittal(self.viewer.ax_sagittal)

            self.press = event.xdata, event.ydata
            return

        elif event.inaxes == self.viewer.list_image[0].frontal_plot.axes:
            for image_itr in self.viewer.list_image:
                image_itr.axial_plot.set_data(image_itr.data[:, :, event.ydata].T)
                image_itr.sagittal_plot.set_data(image_itr.data[:, event.xdata, :].T)

                image_itr.draw_axial(self.viewer.ax_axial)
                image_itr.draw_sagittal(self.viewer.ax_sagittal)

            self.press = event.xdata, event.ydata
            return

        elif event.inaxes == self.viewer.list_image[0].sagittal_plot.axes:
            for image_itr in self.viewer.list_image:
                image_itr.axial_plot.set_data(image_itr.data[:, :, event.ydata].T)
                image_itr.frontal_plot.set_data(image_itr.data[event.xdata, :, :].T)

                image_itr.draw_axial(self.viewer.ax_axial)
                image_itr.draw_frontal(self.viewer.ax_frontal)

            self.press = event.xdata, event.ydata
            return

        else:
            self.press = event.xdata, event.ydata
            return

    def on_scroll(self, event):
        if event.inaxes not in [self.viewer.list_image[0].axial_plot.axes, self.viewer.list_image[0].frontal_plot.axes,
                                self.viewer.list_image[0].sagittal_plot.axes]:
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

        if event.inaxes == self.viewer.list_image[0].axial_plot.axes:
            self.viewer.ax_axial.draw_artist(self.viewer.ax_axial.patch)
            for image_itr in self.viewer.list_image:
                image_itr.draw_axial(self.viewer.ax_axial)
        elif event.inaxes == self.viewer.list_image[0].frontal_plot.axes:
            self.viewer.ax_frontal.draw_artist(self.viewer.ax_frontal.patch)
            for image_itr in self.viewer.list_image:
                image_itr.draw_frontal(self.viewer.ax_frontal)
        elif event.inaxes == self.viewer.list_image[0].sagittal_plot.axes:
            self.viewer.ax_sagittal.draw_artist(self.viewer.ax_sagittal.patch)
            for image_itr in self.viewer.list_image:
                image_itr.draw_sagittal(self.viewer.ax_sagittal)

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
        else:
            return

    def connect(self):
        for image_itr in self.viewer.list_image:
            image_itr.connect(self)

    def disconnect(self):
        for image_itr in self.viewer.list_image:
            image_itr.disconnect(self)


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
        self.list_image[0].min_contrast = val  # type: int
        self.list_image[0].update_min_max_contrast(self.ax_axial, self.ax_frontal, self.ax_sagittal)

    def update_max_contrast(self, val):
        self.list_image[0].max_contrast = val
        self.list_image[0].update_min_max_contrast(self.ax_axial, self.ax_frontal, self.ax_sagittal)

    def show(self):
        self.fig = plt.figure(facecolor='black')
        self.fig.subplots_adjust(bottom=0.1, left=0.1)

        # compute extent based on first image
        shape_first_image = self.list_image[0].data.shape
        axial_extent = 0, 0, shape_first_image[0], shape_first_image[1]
        frontal_extent = 0, 0, shape_first_image[1], shape_first_image[2]
        sagittal_extent = 0, 0, shape_first_image[0], shape_first_image[2]

        self.ax_axial = self.fig.add_subplot(221)
        self.ax_axial.patch.set_facecolor('black')
        self.ax_axial.hold(True)
        for image_itr in self.list_image:
            image_itr.axial_plot = self.ax_axial.imshow(image_itr.data[:, :, int(image_itr.data.shape[2] / 2)].T,
                                                        interpolation='nearest', cmap='gray',
                                                        vmin=image_itr.min_contrast, vmax=image_itr.max_contrast)

        self.ax_frontal = self.fig.add_subplot(222)
        self.ax_frontal.patch.set_facecolor('black')
        self.ax_frontal.hold(True)
        for image_itr in self.list_image:
            image_itr.frontal_plot = self.ax_frontal.imshow(image_itr.data[int(image_itr.data.shape[0]/2), :, :].T,
                                                            interpolation='nearest', cmap='gray',
                                                            vmin=image_itr.min_contrast, vmax=image_itr.max_contrast)

        self.ax_sagittal = self.fig.add_subplot(223)
        self.ax_sagittal.patch.set_facecolor('black')
        self.ax_sagittal.hold(True)
        for image_itr in self.list_image:
            image_itr.sagittal_plot = self.ax_sagittal.imshow(image_itr.data[:, int(image_itr.data.shape[1]/2), :].T,
                                                              interpolation='nearest', cmap='gray',
                                                              vmin=image_itr.min_contrast, vmax=image_itr.max_contrast)

        # TODO: change TrioPlot to consider multiple image
        # TODO: make alpha possible
        # TODO: make masking possible when data is zero (zero is transparent)
        self.trio = TrioPlot(self)
        self.trio.connect()

        from matplotlib.widgets import Slider
        position_slider_min = self.fig.add_axes([0.55, 0.1, 0.35, 0.03])
        slider_min = Slider(position_slider_min, 'Min', self.list_image[0].min_contrast_init,
                            self.list_image[0].max_contrast_init, valinit=self.list_image[0].min_contrast)
        slider_min.on_changed(self.update_min_contrast)
        slider_min.label.set_color('white')
        slider_min.valtext.set_color('white')
        position_slider_max = self.fig.add_axes([0.55, 0.05, 0.35, 0.03])
        slider_max = Slider(position_slider_max, 'Max', self.list_image[0].min_contrast_init,
                            self.list_image[0].max_contrast_init, valinit=self.list_image[0].max_contrast,
                            slidermin=slider_min)
        slider_max.on_changed(self.update_max_contrast)
        slider_max.label.set_color('white')
        slider_max.valtext.set_color('white')
        slider_min.__setattr__('slidermax', slider_max)

        plt.show()


class ImageViewer(Image):
    def __init__(self, image_input):
        super(ImageViewer, self).__init__(image_input)
        from numpy import percentile
        self.min_contrast = percentile(self.data[:], 1)
        self.max_contrast = percentile(self.data[:], 99)
        # check if min != max
        if self.min_contrast == self.max_contrast:
            from numpy import min, max
            self.min_contrast = min(self.data[:])
            self.max_contrast = max(self.data[:])

        self.min_contrast_init = self.min_contrast
        self.max_contrast_init = self.max_contrast
        self.change_orientation('LAS')

        self.axial_plot = None
        self.frontal_plot = None
        self.sagittal_plot = None

    def draw_axial(self, ax_axial):
        ax_axial.draw_artist(self.axial_plot)
        self.axial_plot.figure.canvas.blit(ax_axial.bbox)

    def draw_frontal(self, ax_frontal):
        ax_frontal.draw_artist(self.frontal_plot)
        self.frontal_plot.figure.canvas.blit(ax_frontal.bbox)

    def draw_sagittal(self, ax_sagittal):
        ax_sagittal.draw_artist(self.sagittal_plot)
        self.sagittal_plot.figure.canvas.blit(ax_sagittal.bbox)

    def update_min_max_contrast(self, ax_axial, ax_frontal, ax_sagittal):
        self.axial_plot.set_clim(self.min_contrast, self.max_contrast)
        self.frontal_plot.set_clim(self.min_contrast, self.max_contrast)
        self.sagittal_plot.set_clim(self.min_contrast, self.max_contrast)

        self.draw_axial(ax_axial)
        self.draw_frontal(ax_frontal)
        self.draw_sagittal(ax_sagittal)

    def connect(self, plot):
        """
        connect to all the events we need
        """
        self.cidpress_axial = self.axial_plot.figure.canvas.mpl_connect('button_press_event', plot.on_press)
        self.cidrelease_axial = self.axial_plot.figure.canvas.mpl_connect('button_release_event', plot.on_release)
        self.cidmotion_axial = self.axial_plot.figure.canvas.mpl_connect('motion_notify_event', plot.on_motion)
        self.cidscroll_axial = self.axial_plot.figure.canvas.mpl_connect('scroll_event', plot.on_scroll)

        self.cidpress_frontal = self.frontal_plot.figure.canvas.mpl_connect('button_press_event', plot.on_press)
        self.cidrelease_frontal = self.frontal_plot.figure.canvas.mpl_connect('button_release_event', plot.on_release)
        self.cidmotion_frontal = self.frontal_plot.figure.canvas.mpl_connect('motion_notify_event', plot.on_motion)
        self.cidscroll_frontal = self.frontal_plot.figure.canvas.mpl_connect('scroll_event', plot.on_scroll)

        self.cidpress_sagittal = self.sagittal_plot.figure.canvas.mpl_connect('button_press_event', plot.on_press)
        self.cidrelease_sagittal = self.sagittal_plot.figure.canvas.mpl_connect('button_release_event', plot.on_release)
        self.cidmotion_sagittal = self.sagittal_plot.figure.canvas.mpl_connect('motion_notify_event', plot.on_motion)
        self.cidscroll_sagittal = self.sagittal_plot.figure.canvas.mpl_connect('scroll_event', plot.on_scroll)

    def disconnect(self):
        """
        disconnect all the stored connection ids
        """
        self.axial_plot.figure.canvas.mpl_disconnect(self.cidpress_axial)
        self.axial_plot.figure.canvas.mpl_disconnect(self.cidrelease_axial)
        self.axial_plot.figure.canvas.mpl_disconnect(self.cidmotion_axial)
        self.axial_plot.figure.canvas.mpl_disconnect(self.cidscroll_axial)

        self.frontal_plot.figure.canvas.mpl_disconnect(self.cidpress_frontal)
        self.frontal_plot.figure.canvas.mpl_disconnect(self.cidrelease_frontal)
        self.frontal_plot.figure.canvas.mpl_disconnect(self.cidmotion_frontal)
        self.frontal_plot.figure.canvas.mpl_disconnect(self.cidscroll_frontal)

        self.sagittal_plot.figure.canvas.mpl_disconnect(self.cidpress_sagittal)
        self.sagittal_plot.figure.canvas.mpl_disconnect(self.cidrelease_sagittal)
        self.sagittal_plot.figure.canvas.mpl_disconnect(self.cidmotion_sagittal)
        self.sagittal_plot.figure.canvas.mpl_disconnect(self.cidscroll_sagittal)

# ======================================================================================================================
# Start program
# ======================================================================================================================
if __name__ == "__main__":
    parser = Parser(__file__)
    parser.usage.set_description('Volume Viewer')
    parser.add_option("-i", [[','], 'file'], "file", True)
    arguments = parser.parse(sys.argv[1:])

    viewer = VolViewer()
    for image_str in arguments["-i"]:
        image = Image(image_str)
        viewer.add_image(image)
    viewer.show()

#!/usr/bin/env python
########################################################################################################################
#
# Asman et al. groupwise multi-atlas segmentation method implementation
#
# ----------------------------------------------------------------------------------------------------------------------
# Copyright (c) 2014 Polytechnique Montreal <www.neuro.polymtl.ca>
# Author: Augustin Roux
# Modified: 2014-11-20
#
# About the license: see the file LICENSE.TXT
########################################################################################################################

# TODO change 'target' by 'input'

from scipy.misc import toimage
from msct_pca import PCA
import numpy as np
from math import sqrt
from math import exp
from msct_image import Image
from msct_parser import *
import matplotlib.pyplot as plt
import sct_utils as sct


class Param:
    def __init__(self):
        self.debug = 1
        self.path_dictionnary = '/home/django/aroux/Desktop/data_asman/Dictionnary/'
        self.patient_id = ['09', '24', '30', '31', '32', '25', '10', '08', '11', '16', '17', '18']
        #self.patient_id = ['09', '24', '30', '31', '32']
        self.include_GM = 0
        self.split_data = 1  # this flag enables to duplicate the image in the right-left direction in order to have more dataset for the PCA
        self.verbose = 0


########################################################################################################################
######-------------------------------------------------- MAIN ----------------------------------------------------######
########################################################################################################################

def main():
    v = param.verbose

    if param.debug:
        print '\n*** WARNING: DEBUG MODE ON ***\n'
        fname_input = "/home/django/aroux/Desktop/data_asman/Dictionnary/errsm_34.nii.gz"
        fname_input = "/home/django/aroux/Desktop/data_asman/Dictionnary/errsm_34_seg_in.nii.gz"
    else:
        parser = Parser(__file__)
        parser.usage.set_description('Project all the imput image slices on a PCA generated from set of t2star images')
        parser.add_option("-input", "file", "t2star you want to project", True, "t2star.nii.gz")

        # Getting the arguments
        arguments = parser.parse(sys.argv[1:])
        fname_input = arguments["-input"]

    # construct target image
    target_image = Image(fname_input, split=param.split_data)

    # build the appearance model
    appearance_model = AppearanceModel(target_image)

    appearance_model.pca.show(param.split_data)

    appearance_model.plot_omega()

    appearance_model.show_projected_target()


########################################################################################################################
######------------------------------------------------- Classes --------------------------------------------------######
########################################################################################################################

# ----------------------------------------------------------------------------------------------------------------------
# APPEARANCE MODEL -----------------------------------------------------------------------------------------------------
class AppearanceModel:
    def __init__(self, target_image=None, param=None):
        if param is None:
            self.param = Param()
        else:
            self.param = param
            # Load all the images's slices from param.path_dictionnary
        self.list_atlas_seg = self.load_dictionnary(self.param.split_data)
        # Construct a dataset composed of all the slices
        dataset = self.construct_dataset()
        sct.printv("The shape of the dataset used for the PCA is {}".format(dataset.shape), verbose=self.param.verbose)
        # Instantiate a PCA object given the dataset just build
        self.pca = PCA(dataset)
        # Get the target image
        self.target = target_image
        # coord_projected_target is a list of all the coord of the target's projected slices
        self.coord_projected_target = self.pca.project(target_image) if target_image is not None else None
        self.beta = self.compute_beta()

    # ------------------------------------------------------------------------------------------------------------------
    # Load the dictionary:
    # each slice of each patient will be load separately with its corresponding GM segmentation
    # they will be stored as tuples in list_atlas_seg
    def load_dictionnary(self, split_data):
        # init
        list_atlas_seg = []
        # loop across all the volume
        for id in param.patient_id:
            #atlas = Image(self.param.path_dictionnary + 'errsm_' + id + '.nii.gz')
            atlas = Image(param.path_dictionnary + 'errsm_' + id + '_seg_in.nii.gz')

            if split_data:
                if self.param.include_GM:
                    seg = Image(self.param.path_dictionnary + 'errsm_' + id + '_GMr.nii.gz')
                    index_s = 0
                    for slice in atlas.data:
                        left_slice, right_slice = split(slice)
                        seg_slice = seg.data[index_s]
                        left_slice_seg, right_slice_seg = split(seg_slice)
                        list_atlas_seg.append((left_slice, left_slice_seg))
                        list_atlas_seg.append((right_slice, right_slice_seg))
                        index_s += 1
                else:
                    index_s = 0
                    for slice in atlas.data:
                        left_slice, right_slice = split(slice)
                        list_atlas_seg.append((left_slice, None))
                        list_atlas_seg.append((right_slice, None))
                        index_s += 1
            else:
                if param.include_GM:
                    seg = Image(param.path_dictionnary + 'errsm_' + id + '_GMr.nii.gz')
                    index_s = 0
                    for slice in atlas.data:
                        seg_slice = seg.data[index_s]
                        list_atlas_seg.append((slice, seg_slice))
                        index_s += 1
                else:
                    for slice in atlas.data:
                        list_atlas_seg.append((slice, None))
        return list_atlas_seg

    # ------------------------------------------------------------------------------------------------------------------
    # in order to build the PCA from all the J atlases, we must construct a matrix of J columns and N rows,
    # with N the dimension of flattened images
    def construct_dataset(self):
        list_atlas_seg = self.list_atlas_seg
        dataset = []
        for atlas_slice in list_atlas_seg:
            dataset.append(atlas_slice[0].flatten())
        return np.asarray(dataset).T

    # ------------------------------------------------------------------------------------------------------------------
    def show_projected_target(self):
        # Retrieving projected image from the mean image & its coordinates
        import copy

        index = 0
        fig1 = plt.figure()
        fig2 = plt.figure()
        # loop across all the projected slices coord
        for coord in self.coord_projected_target:
            img_reducted = copy.copy(self.pca.mean_image)
            # loop across coord and build projected image
            for i in range(0, coord.shape[0]):
                img_reducted += int(coord[i][0]) * self.pca.W.T[i].reshape(self.pca.N, 1)

            if self.param.split_data:
                n = int(sqrt(self.pca.N * 2))
            else:
                n = int(sqrt(self.pca.N))

            # plot mean image
            # if self.param.split_data:
            #     imgplot = plt.imshow(self.pca.mean_image.reshape(n / 2, n))
            # else:
            #     imgplot = plt.imshow(self.pca.mean_image.reshape(n, n))
            # imgplot.set_interpolation('nearest')
            # imgplot.set_cmap('gray')
            # plt.title('Mean Image')
            # plt.show()
            #
            # Plot original image
            orig_ax = fig1.add_subplot(10, 3, index)
            orig_ax.set_title('original slice {} '.format(index))
            if self.param.split_data:
                imgplot = orig_ax.imshow(self.target.data[index, :, :].reshape(n / 2, n))
            else:
                imgplot = orig_ax.imshow(self.target.data[index].reshape(n, n))
            imgplot.set_interpolation('nearest')
            imgplot.set_cmap('gray')
            # plt.title('Original Image')
            # plt.show()

            index += 1
            # Plot projected image image
            proj_ax = fig2.add_subplot(10, 3, index)
            proj_ax.set_title('slice {} projected'.format(index))
            if self.param.split_data:
                imgplot = proj_ax.imshow(img_reducted.reshape(n / 2, n))
                #imgplot = plt.imshow(img_reducted.reshape(n / 2, n))
            else:
                # imgplot = plt.imshow(img_reducted.reshape(n, n))
                imgplot = proj_ax.imshow(img_reducted.reshape(n, n))
            imgplot.set_interpolation('nearest')
            imgplot.set_cmap('gray')
            # plt.title('Projected Image')
            # plt.show()
        plt.show()

    # ------------------------------------------------------------------------------------------------------------------
    # plot the pca and the target projection if target is provided
    def plot_omega(self):
        self.pca.plot_omega(target_coord=self.coord_projected_target) if self.coord_projected_target is not None \
            else self.pca.plot_omega()


# ----------------------------------------------------------------------------------------------------------------------
# RIGID REGISTRATION ---------------------------------------------------------------------------------------------------
class RigidRegistration:
    def __init__(self, appearance_model):
        self.beta = compute_beta(appearance_model)
        self.mu = appearance_model.pca.omega.dot(self.beta)
        self.sigma = self.compute_sigma(appearance_model)

    def compute_sigma(self, appearance_model):
        sigma = []
        j = 0
        for w_v in appearance_model.pca.omega.T:
            sig = 0
            for w_j in w_v:
                sig += self.beta[j]*(w_j - self.mu[j])
            sigma.append(sig)
        return sigma



# ----------------------------------------------------------------------------------------------------------------------
# beta is the model similarity between all the individual images and our input image
# beta = (1/Z)exp(-theta*square_norm(omega-omega_j))
# Z is the partition function that enforces the constraint tha sum(beta)=1
def compute_beta(appearance_model):
    beta = []
    theta = 1
    if appearance_model.coord_projected_target:
        # in omega matrix, each column correspond to the pojection of one of the original data image,
        # the transpose operator .T enable the loop to iterate over all the images coord
        for omega_j in appearance_model.pca.omega.T:
            square_norm = np.linalg.norm((omega_j - appearance_model.coord_projected_target), 2)
            beta.append(exp(-theta*square_norm))
    else:
        raise Exception("No projected input in the appearance model")
    Z = sum(beta)
    beta = np.asarray((1/Z)*beta)
    return beta

########################################################################################################################
######------------------------------------------------ FUNCTIONS -------------------------------------------------######
########################################################################################################################


# ----------------------------------------------------------------------------------------------------------------------
# Split a slice in two slices, used to deal with actual loss of data
def split(slice):
    left_slice = []
    right_slice = []
    column_length = slice.shape[1]
    i = 0
    for column in slice:
        if i < column_length / 2:
            left_slice.append(column)
        else:
            right_slice.insert(0, column)
        i += 1
    left_slice = np.asarray(left_slice)
    right_slice = np.asarray(right_slice)
    assert (left_slice.shape == right_slice.shape), \
        str(left_slice.shape) + '==' + str(right_slice.shape) + \
        'You should check that the first dim of your image (or slice) is an odd number'
    return left_slice, right_slice


# ----------------------------------------------------------------------------------------------------------------------
def show(coord_projected_img, pca, target):
    import matplotlib.pyplot as plt

    # Retrieving projected image from the mean image & its coordinates
    import copy

    img_reducted = copy.copy(pca.mean_image)
    for i in range(0, coord_projected_img.shape[0]):
        img_reducted += int(coord_projected_img[i][0]) * pca.W.T[i].reshape(pca.N, 1)

    if param.split_data:
        n = int(sqrt(pca.N * 2))
    else:
        n = int(sqrt(pca.N))
    if param.split_data:
        imgplot = plt.imshow(pca.mean_image.reshape(n, n / 2))
    else:
        imgplot = plt.imshow(pca.mean_image.reshape(n, n))
    imgplot.set_interpolation('nearest')
    imgplot.set_cmap('gray')
    plt.title('Mean Image')
    plt.show()
    if param.split_data:
        imgplot = plt.imshow(target.reshape(n, n / 2))
    else:
        imgplot = plt.imshow(target.reshape(n, n))
    imgplot.set_interpolation('nearest')
    #imgplot.set_cmap('gray')
    plt.title('Original Image')
    plt.show()
    if param.split_data:
        imgplot = plt.imshow(img_reducted.reshape(n, n / 2))
    else:
        imgplot = plt.imshow(img_reducted.reshape(n, n))
    imgplot.set_interpolation('nearest')
    #imgplot.set_cmap('gray')
    plt.title('Projected Image')
    plt.show()


# ----------------------------------------------------------------------------------------------------------------------
# This little loop save projection through several pcas with different k i.e. different number of modes
def save(dataset, list_atlas_seg):
    import scipy
    import copy

    betas = [0.6, 0.7, 0.75, 0.8, 0.82, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95]
    target = list_atlas_seg[8][0].flatten()
    for beta in betas:
        pca = PCA(dataset, beta)
        coord_projected_img = pca.project(target)
        img_reducted = copy.copy(pca.mean_image)
        n = int(sqrt(pca.N * 2))
        for i in range(0, coord_projected_img.shape[0]):
            img_reducted += int(coord_projected_img[i][0]) * pca.W.T[i].reshape(pca.N, 1)
        scipy.misc.imsave("/home/django/aroux/Desktop/pca_modesInfluence/" + str(pca.kept) + "modes.jpeg",
                          img_reducted.reshape(n, n / 2))


if __name__ == "__main__":
    param = Param()
    main()


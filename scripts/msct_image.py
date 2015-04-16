#!/usr/bin/env python
#########################################################################################
#
# Image class implementation
#
#
# ---------------------------------------------------------------------------------------
# Copyright (c) 2015 Polytechnique Montreal <www.neuro.polymtl.ca>
# Authors: Augustin Roux, Benjamin De Leener
# Modified: 2015-02-20
#
# About the license: see the file LICENSE.TXT
#########################################################################################


class Image(object):
    """

    """
    def __init__(self, param=None, hdr=None, orientation=None, absolute_path="", verbose=1):
        from numpy import zeros, ndarray, generic
        from sct_utils import extract_fname

        # initialization of all parameters
        self.data = None
        self.hdr = None
        self.orientation = None
        self.absolute_path = ""
        self.path = ""
        self.file_name = ""
        self.ext = ""
        self.dim = None

        # load an image from file
        if type(param) is str:
            self.load_from_path(param, verbose)
        # copy constructor
        elif type(param) is Image:
            self.copy(param)
        # create an empty image (full of zero) of dimension [dim]. dim must be [x,y,z] or (x,y,z). No header.
        elif type(param) is list:
            self.data = zeros(param)
            self.dim = param
            self.hdr = hdr
            self.orientation = orientation
            self.absolute_path = absolute_path
            self.path, self.file_name, self.ext = extract_fname(absolute_path)
        # create a copy of im_ref
        elif isinstance(param, (ndarray, generic)):
            self.data = param
            self.dim = self.data.shape
            self.hdr = hdr
            self.orientation = orientation
            self.absolute_path = absolute_path
            self.path, self.file_name, self.ext = extract_fname(absolute_path)
        else:
            raise TypeError(' Image constructor takes at least one argument.')

    def __deepcopy__(self, memo):
        from copy import deepcopy
        return type(self)(deepcopy(self.data, memo), deepcopy(self.hdr, memo), deepcopy(self.orientation, memo), deepcopy(self.absolute_path, memo))

    def copy(self, image_to_copy=None):
        from copy import deepcopy
        from sct_utils import extract_fname
        if image_to_copy is not None:
            self.data = deepcopy(image_to_copy.data)
            self.dim = deepcopy(image_to_copy.dim)
            self.hdr = deepcopy(image_to_copy.hdr)
            self.orientation = deepcopy(image_to_copy.orientation)
            self.absolute_path = deepcopy(image_to_copy.absolute_path)
            self.path, self.file_name, self.ext = extract_fname(self.absolute_path)
        else:
            return deepcopy(self)

    def load_from_path(self, path, verbose):
        """
        This function load an image from an absolute path using nibabel library
        :param path: path of the file from which the image will be loaded
        :return:
        """
        from nibabel import load, spatialimages
        from sct_utils import check_file_exist, printv, extract_fname
        from sct_orientation import get_orientation

        check_file_exist(path, verbose=verbose)
        try:
            im_file = load(path)
            self.orientation = get_orientation(path)
            self.data = im_file.get_data()
            self.hdr = im_file.get_header()
            self.absolute_path = path
            self.path, self.file_name, self.ext = extract_fname(path)
        except spatialimages.ImageFileError:
            printv('Error: make sure ' + path + ' is an image.')

    def setFileName(self, filename):
        from sct_utils import extract_fname
        self.absolute_path = filename
        self.path, self.file_name, self.ext = extract_fname(filename)

    def change_type(self, voxel_datatype=''):
        """
        Change the voxel type of the image
        :param type:    if not set, the image is saved in standard type
                        if 'minimize', image space is minimize
                        if 'minimize_int', image space is minimize and values are approximated to integers
                        (2, 'uint8', np.uint8, "NIFTI_TYPE_UINT8"),
                        (4, 'int16', np.int16, "NIFTI_TYPE_INT16"),
                        (8, 'int32', np.int32, "NIFTI_TYPE_INT32"),
                        (16, 'float32', np.float32, "NIFTI_TYPE_FLOAT32"),
                        (32, 'complex64', np.complex64, "NIFTI_TYPE_COMPLEX64"),
                        (64, 'float64', np.float64, "NIFTI_TYPE_FLOAT64"),
                        (256, 'int8', np.int8, "NIFTI_TYPE_INT8"),
                        (512, 'uint16', np.uint16, "NIFTI_TYPE_UINT16"),
                        (768, 'uint32', np.uint32, "NIFTI_TYPE_UINT32"),
                        (1024,'int64', np.int64, "NIFTI_TYPE_INT64"),
                        (1280, 'uint64', np.uint64, "NIFTI_TYPE_UINT64"),
                        (1536, 'float128', _float128t, "NIFTI_TYPE_FLOAT128"),
                        (1792, 'complex128', np.complex128, "NIFTI_TYPE_COMPLEX128"),
                        (2048, 'complex256', _complex256t, "NIFTI_TYPE_COMPLEX256"),
        :return:
        """
        if voxel_datatype == '':
            voxel_datatype = self.hdr.get_data_dtype()

        if voxel_datatype == 'minimize' or voxel_datatype == 'minimize_int':
            from numpy import nanmax, nanmin
            # compute max value in the image and choose the best pixel type to represent all the pixels within smallest
            # memory space
            # warning: does not take intensity resolution into account, neither complex voxels
            max_vox = nanmax(self.data)
            min_vox = nanmin(self.data)

            # check if voxel values are real or integer
            is_integer = True
            if voxel_datatype == 'minimize':
                for vox in self.data.flatten():
                    if int(vox) != vox:
                        is_integer = False
                        break

            if is_integer:
                if min_vox >= 0:  # unsigned
                    from numpy import iinfo, uint8, uint16, uint32, uint64
                    if max_vox <= iinfo(uint8).max:
                        voxel_datatype = 'uint8'
                    elif max_vox <= iinfo(uint16):
                        voxel_datatype = 'uint16'
                    elif max_vox <= iinfo(uint32).max:
                        voxel_datatype = 'uint32'
                    elif max_vox <= iinfo(uint64).max:
                        voxel_datatype = 'uint64'
                    else:
                        raise ValueError("Maximum value of the image is to big to be represented.")
                else:
                    from numpy import iinfo, int8, int16, int32, int64
                    if max_vox <= iinfo(int8).max and min_vox >= iinfo(int8).min:
                        voxel_datatype = 'int8'
                    elif max_vox <= iinfo(int16).max and min_vox >= iinfo(int16).min:
                        voxel_datatype = 'int16'
                    elif max_vox <= iinfo(int32).max and min_vox >= iinfo(int32).min:
                        voxel_datatype = 'int32'
                    elif max_vox <= iinfo(int64).max and min_vox >= iinfo(int64).min:
                        voxel_datatype = 'int64'
                    else:
                        raise ValueError("Maximum value of the image is to big to be represented.")
            else:
                from numpy import finfo, float32, float64
                # if max_vox <= np.finfo(np.float16).max and min_vox >= np.finfo(np.float16).min:
                #    type = 'np.float16' # not supported by nibabel
                if max_vox <= finfo(float32).max and min_vox >= finfo(float32).min:
                    voxel_datatype = 'float32'
                elif max_vox <= finfo(float64).max and min_vox >= finfo(float64).min:
                    voxel_datatype = 'float64'

        # print "The image has been set to "+type+" (previously "+str(self.hdr.get_data_dtype())+")"
        # change type of data in both numpy array and nifti header
        type_build = eval(voxel_datatype)
        self.data = type_build(self.data)
        self.hdr.set_data_dtype(voxel_datatype)

    def save(self, datatype=''):
        """
        Write an image in a nifti file
        :param datatype:    if not set, the image is saved in standard type
                        if 'minimize', image space is minimize
                        (2, 'uint8', np.uint8, "NIFTI_TYPE_UINT8"),
                        (4, 'int16', np.int16, "NIFTI_TYPE_INT16"),
                        (8, 'int32', np.int32, "NIFTI_TYPE_INT32"),
                        (16, 'float32', np.float32, "NIFTI_TYPE_FLOAT32"),
                        (32, 'complex64', np.complex64, "NIFTI_TYPE_COMPLEX64"),
                        (64, 'float64', np.float64, "NIFTI_TYPE_FLOAT64"),
                        (256, 'int8', np.int8, "NIFTI_TYPE_INT8"),
                        (512, 'uint16', np.uint16, "NIFTI_TYPE_UINT16"),
                        (768, 'uint32', np.uint32, "NIFTI_TYPE_UINT32"),
                        (1024,'int64', np.int64, "NIFTI_TYPE_INT64"),
                        (1280, 'uint64', np.uint64, "NIFTI_TYPE_UINT64"),
                        (1536, 'float128', _float128t, "NIFTI_TYPE_FLOAT128"),
                        (1792, 'complex128', np.complex128, "NIFTI_TYPE_COMPLEX128"),
                        (2048, 'complex256', _complex256t, "NIFTI_TYPE_COMPLEX256"),
        """
        from nibabel import Nifti1Image, save

        if datatype != '':
            self.change_type(datatype)

        self.hdr.set_data_shape(self.data.shape)
        img = Nifti1Image(self.data, None, self.hdr)
        print 'saving ' + self.path + self.file_name + self.ext + '\n'
        save(img, self.path + self.file_name + self.ext)

    def flatten(self):
        """
        Flatten the array in a single dimension vector, its shape will be (d, 1) compared to the flatten built in method
         which would have returned (d,)
        :return: image data flattened as a numpy array
        """
        # return self.data.flatten().reshape(self.data.flatten().shape[0], 1)
        return self.data.flatten()

    #
    def slices(self):
        """
        :return: a list of the image slices flattened
        """
        return [slc.flatten() for slc in self.data]

    def getNonZeroCoordinates(self, sorting=None, reverse_coord=False):
        """
        This function return all the non-zero coordinates that the image contains.
        Coordinate list can also be sorted by x, y, z, or the value with the parameter sorting='x', sorting='y', sorting='z' or sorting='value'
        If reverse_coord is True, coordinate are sorted from larger to smaller.
        """
        from msct_types import Coordinate

        x, y, z = (self.data > 0).nonzero()
        list_coordinates = [Coordinate([x[i], y[i], z[i], self.data[x[i], y[i], z[i]]]) for i in range(0, len(x))]

        if sorting is not None:
            if reverse_coord not in [True, False]:
                raise ValueError('reverse_coord parameter must be a boolean')

            if sorting == 'x':
                list_coordinates = sorted(list_coordinates, key=lambda obj: obj.x, reverse=reverse_coord)
            elif sorting == 'y':
                list_coordinates = sorted(list_coordinates, key=lambda obj: obj.y, reverse=reverse_coord)
            elif sorting == 'z':
                list_coordinates = sorted(list_coordinates, key=lambda obj: obj.z, reverse=reverse_coord)
            elif sorting == 'value':
                list_coordinates = sorted(list_coordinates, key=lambda obj: obj.value, reverse=reverse_coord)
            else:
                raise ValueError("sorting parameter must be either 'x', 'y', 'z' or 'value'")

        return list_coordinates

    # crop the image in order to keep only voxels in the mask, therefore the mask's slices must be squares or
    # rectangles of the same size
    # This method is called in sct_crop_over_mask script
    def crop_from_square_mask(self, mask):
        from numpy import asarray

        array = self.data
        data_mask = mask.data
        print 'ORIGINAL SHAPE: ', array.shape, '   ==   ', data_mask.shape
        array = asarray(array)
        data_mask = asarray(data_mask)
        new_data = []
        buffer_temp = []
        buffer_mask = []
        s = 0
        r = 0
        ok = 0
        for slice_image in data_mask:
            # print 'SLICE ', s, slice
            for row in slice_image:
                if sum(row) > 0:
                    buffer_mask.append(row)
                    buffer_temp.append(array[s][r])
                    # print 'OK1', ok
                    ok += 1
                r += 1
            new_slice_mask = asarray(buffer_mask).T
            new_slice = asarray(buffer_temp).T
            r = 0
            buffer_temp = []
            for row in new_slice_mask:
                if sum(row) != 0:
                    buffer_temp.append(new_slice[r])
                r += 1
            # print buffer_temp
            new_slice = asarray(buffer).T
            r = 0
            buffer_mask = []
            buffer_temp = []
            new_data.append(new_slice)
            s += 1
        new_data = asarray(new_data)
        # print data_mask
        print 'SHAPE ', new_data.shape
        self.data = new_data

    def invert(self):
        self.data = self.data.max() - self.data
        return self

    def change_orientation(self, orientation='RPI'):
        """
        This function changes the orientation of the data by swapping the image axis.
        Warning: the nifti image header is not change!!!
        :param orientation: string of three character representing the new orientation (ex: AIL, default: RPI)
        :return:
        """
        opposite_character = {'L': 'R', 'R': 'L', 'A': 'P', 'P': 'A', 'I': 'S', 'S': 'I'}

        if self.orientation is None:
            from sct_orientation import get_orientation
            self.orientation = get_orientation(self.file_name)

        # change the orientation of the image
        perm = [0, 1, 2]
        inversion = [1, 1, 1]
        for i, character in enumerate(self.orientation):
            try:
                perm[i] = orientation.index(character)
            except ValueError:
                perm[i] = orientation.index(opposite_character[character])
                inversion[i] = -1

        # axes inversion
        self.data = self.data[::inversion[0], ::inversion[1], ::inversion[2]]

        # axes manipulations
        from numpy import swapaxes
        if perm == [1, 0, 2]:
            self.data = swapaxes(self.data, 0, 1)
        elif perm == [2, 1, 0]:
            self.data = swapaxes(self.data, 0, 2)
        elif perm == [0, 2, 1]:
            self.data = swapaxes(self.data, 1, 2)
        elif perm == [2, 1, 0]:
            self.data = swapaxes(self.data, 0, 2)
        elif perm == [2, 0, 1]:
            self.data = swapaxes(self.data, 0, 2)  # transform [2, 0, 1] to [1, 0, 2]
            self.data = swapaxes(self.data, 0, 1)  # transform [1, 0, 2] to [0, 1, 2]
        elif perm == [1, 2, 0]:
            self.data = swapaxes(self.data, 0, 2)  # transform [1, 2, 0] to [0, 2, 1]
            self.data = swapaxes(self.data, 1, 2)  # transform [0, 2, 1] to [0, 1, 2]
        elif perm == [0, 1, 2]:
            # do nothing
            pass
        else:
            print 'Error: wrong orientation'

        self.orientation = orientation

    def show(self):
        from matplotlib.pyplot import imshow, show
        imgplot = imshow(self.data)
        imgplot.set_cmap('gray')
        imgplot.set_interpolation('nearest')
        show()

# ======================================================================================================================
# Start program
# ======================================================================================================================
if __name__ == "__main__":
    from msct_parser import Parser
    import sys

    parser = Parser(__file__)
    parser.usage.set_description('Image')
    parser.add_option("-i", "file", "file", True)
    arguments = parser.parse(sys.argv[1:])

    image = Image(arguments["-i"])
    image.change_type('minimize')
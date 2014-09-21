#!/usr/bin/env python
#########################################################################################
#
# Compute magnetization transfer ratio (MTR).
#
# ---------------------------------------------------------------------------------------
# Copyright (c) 2014 Polytechnique Montreal <www.neuro.polymtl.ca>
# Authors: Julien Cohen-Adad
# Modified: 2014-09-21
#
# About the license: see the file LICENSE.TXT
#########################################################################################


import sys
import os
import getopt
import commands
import sct_utils as sct
import time

# DEFAULT PARAMETERS
class param:
    ## The constructor
    def __init__(self):
        self.debug = 1
        self.register = 1
        self.verbose = 1
        self.file_out = 'mtr'
        self.remove_tmp_files = 1


# main
#=======================================================================================================================
def main():

    # Initialization
    fname_mt0 = ''
    fname_mt1 = ''
    file_out = param.file_out
    register = param.register
    fsloutput = 'export FSLOUTPUTTYPE=NIFTI; ' # for faster processing, all outputs are in NIFTI
    remove_tmp_files = param.remove_tmp_files
    verbose = param.verbose

    # get path of the toolbox
    status, path_sct = commands.getstatusoutput('echo $SCT_DIR')

    # Parameters for debug mode
    if param.debug:
        print '\n*** WARNING: DEBUG MODE ON ***\n'
        fname_mt0 = path_sct+'/testing/data/errsm_23/mt/mt0.nii.gz'
        fname_mt1 = path_sct+'/testing/data/errsm_23/mt/mt1.nii.gz'
        register = 1
        verbose = 1

    # Check input parameters
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hi:j:r:v:x')
    except getopt.GetoptError:
        usage()
    for opt, arg in opts:
        if opt == '-h':
            usage()
        elif opt in ('-i'):
            param.fname_mt0 = arg
        elif opt in ('-j'):
            param.fname_mt1 = arg
        elif opt in ('-r'):
            param.remove_tmp_files = int(arg)
        elif opt in ('-v'):
            param.verbose = int(arg)
        elif opt in ('-x'):
            param.register = int(arg)

    # display usage if a mandatory argument is not provided
    if fname_mt0 == '' or fname_mt1 == '':
        sct.printv('ERROR: All mandatory arguments are not provided. See usage.', 1, 'error')
        usage()

    # display input parameters
    sct.printv('\nInput parameters:', verbose)
    sct.printv('  mt0 ...................'+fname_mt0, verbose)
    sct.printv('  mt1 ...................'+fname_mt1, verbose)
    sct.printv('  register ..............'+str(register), verbose)

    # check existence of input files
    sct.check_file_exist(fname_mt0, verbose)
    sct.check_file_exist(fname_mt1, verbose)

    # Extract path/file/extension
    path_mt0, file_mt0, ext_mt0 = sct.extract_fname(fname_mt0)
    path_out, file_out, ext_out = '', file_out, ext_mt0

    # create temporary folder
    path_tmp = sct.slash_at_the_end('tmp.'+time.strftime("%y%m%d%H%M%S"), 1)
    sct.run('mkdir '+path_tmp, verbose)

    # Copying input data to tmp folder and convert to nii
    sct.printv('\nCopying input data to tmp folder and convert to nii...', verbose)
    sct.run('sct_c3d '+fname_mt0+' -o '+path_tmp+'mt0.nii')
    sct.run('sct_c3d '+fname_mt1+' -o '+path_tmp+'mt1.nii')

    # go to tmp folder
    os.chdir(path_tmp)

    # register MTC0 on MTC1
    if register:
        sct.printv('\nRegister mt0 on mt1...', verbose)
        sct.run("sct_register_multimodal -i mt0.nii -d mt1.nii", verbose)
        file_mt0_tmp = "mt0_reg.nii"
    else:
        file_mt0_tmp = "mt0.nii"

    # compute MTR
    sct.printv('\nCompute MTR...', verbose)
    sct.run(fsloutput+'fslmaths -dt double '+file_mt0_tmp+' -sub mt1.nii -mul 100 -div '+file_mt0+' -thr 0 -uthr 100 mtr.nii', verbose)

    # come back to parent folder
    os.chdir('..')

    # Generate output files
    sct.printv('\nGenerate output files...', verbose)
    sct.generate_output_file(path_tmp+'mtr.nii', path_out, file_out, ext_out)
    if register:
        sct.generate_output_file(path_tmp+'mt0_reg.nii', path_out, file_mt0+'_reg', ext_out)

    # Remove temporary files
    if remove_tmp_files == 1:
        print('\nRemove temporary files...')
        sct.run('rm -rf '+path_tmp)

    # to view results
    print '\nDone! To view results, type:'
    print 'fslview '+fname_mt0+' '+fname_mt1+' '+file_out+' &'
    print


# Print usage
# ==========================================================================================
def usage():
    print """
"""+os.path.basename(__file__)+"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Part of the Spinal Cord Toolbox <https://sourceforge.net/projects/spinalcordtoolbox>

DESCRIPTION
  Compute magnetization transfer ratio (MTR). Output is given in percentage.

USAGE
  """+os.path.basename(__file__)+""" -i <mt0> -j <mt1>

MANDATORY ARGUMENTS
  -i <mt0>         image without MT pulse
  -j <mt1>         image with MT pulse

OPTIONAL ARGUMENTS
  -x {0,1}         register mt0 on mt1. Default="""+str(param.register)+"""
  -r {0,1}         remove temporary files. Default="""+str(param.remove_tmp_files)+"""
  -v {0,1}         verbose. Default="""+str(param.verbose)+"""
  -h               help. Show this message

EXAMPLE
  """+os.path.basename(__file__)+""" -i mt0.nii.gz -j mt1.nii.gz -x 1\n"""

    # exit program
    sys.exit(2)


#=======================================================================================================================
# Start program
#=======================================================================================================================
if __name__ == "__main__":
    # initialize parameters
    param = param()
    # call main function
    main()
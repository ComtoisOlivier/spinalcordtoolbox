#!/usr/bin/env python
#########################################################################################
#
# Validation of WM atlas
#
# This script generates a synthetic volume from an atlas of white matter
# It then estimates the np.mean value within each tract using different methods.
# The estimation is afterwards substracted to the real value in the tracts to find the absolute error
# For each method, this process is repeated a certain number of iterations(bootstrap_iterations)
# The np.mean of all absolute deviations and the respective np.std is calculated for each method
#
# ---------------------------------------------------------------------------------------
# Copyright (c) 2014 Polytechnique Montreal <www.neuro.polymtl.ca>
# Author: Charles Naaman, Julien Cohen-Adad
# Modified: 2014-11-25
#
# About the license: see the file LICENSE.TXT
#########################################################################################

# TODO : Make a function for which the estimation is made from a particular method(remove code duplication in bootstrap loop and parameters initialization)
# TODO : Make a function to print results depending on which methods are used
# TODO : Add usage 
# TODO : Add class initialization of default parameters (not working for some reason)


# Import common Python libraries
import os
import sys
import time
import commands
import shutil
import numpy as np
status, path_sct = commands.getstatusoutput('echo $SCT_DIR')
# append path that contains scripts, to be able to load modules
sys.path.append(path_sct + '/scripts')
import sct_utils as sct
from generate_phantom import phantom_generation, get_tracts, save_3D_nparray_niftii

def main():
    # Parameters

    # param for validation
    bootstrap_iter = 200
    folder_atlas = '../WM_atlas_generation/WMtracts_outputs/final_results/'  # path to atlas. add / at the end
    user_tract = 'charles_tract_,julien_tract_,tanguy_tract_,simon_tract_'
    std_noise_list = [0.00001, 1, 5, 10, 20, 50]  # standard deviation of the noise added to the generated phantom
    range_tract_list = [0, 5, 10, 20, 50]
    fixed_range = 10
    fixed_noise = 10
    results_folder = 'results/'  # add / at the end

    # create output folder
    create_folder(results_folder)

    # loop across noise levels
    range_tract = fixed_range
    for std_noise in std_noise_list:
        results_file = 'results_noise'+str(std_noise)+'_range'+str(range_tract)+'.txt'
        validate_atlas(folder_atlas, bootstrap_iter, std_noise, range_tract, results_folder+results_file, user_tract)

    # loop across tract ranges
    std_noise = fixed_noise
    for range_tract in range_tract_list:
        results_file = 'results_noise'+str(std_noise)+'_range'+str(range_tract)+'.txt'
        validate_atlas(folder_atlas, bootstrap_iter, std_noise, range_tract, results_folder+results_file, user_tract)


def validate_atlas(folder_atlas, bootstrap_iterations, std_noise, range_tract, results_file, user_tract):
    # Parameters
    #bootstrap_iterations = 2  # number of bootstrap iterations. Default=200
    #folder_atlas = '../WM_atlas_generation/WMtracts_outputs/final_results/'  # add / at the end
    folder_cropped_atlas = "cropped_atlas/"
    crop = 1  # crop atlas, default=1
    zcrop_ind = '10,110,210,310,410'
    generated_phantom = "WM_phantom.nii.gz"
    generated_phantom_noise = "WM_phantom_noise.nii.gz"
    tracts_sum_img = "tracts_sum.nii.gz"
    true_value = 40
    # np.std_noise = 10
    #range_tract = 10
    # spec_tracts = np.arange(30)
    spec_tracts = 2, 17
    metrics_estimation_results = "metric_label.txt"
    dorsal_column_labels = '0,1,15,16'
    #results_file = "atlas_validation.txt"
    partial_vol_corr = 0

    # Parameters for the manual estimation
    # These parameters are associated with the manually created masks
    dorsal_column_mask_index = 30
    man_mask_index = 2,17,dorsal_column_mask_index
    #mask_prefix = 'mask_tract_'
    mask_prefix = ['charles_tract_', 'julien_tract_', 'tanguy_tract_', 'simon_tract_']
    mask_folder = ['manual_masks/charles/', 'manual_masks/julien/', 'manual_masks/tanguy/', 'manual_masks/simon/']
    mask_ext = '.nii.gz'

    start_time = time.time() # save start time for duration

    # # check input parameters
    # try:
    #     opts, args = getopt.getopt(sys.argv[1:], 'ha:cpt:n:s:f:g:b:l:m:d:r:z:') # define flag
    # except getopt.GetoptError as err: # check if the arguments are defined
    #     print str(err) # error
    #    # usage(label_title, label_name, label_num,fname_tracts) # display usage
    # for opt, arg in opts: # explore flags
    #     elif opt in '-a': # path to atlas
    #         folder_atlas = arg
    #     elif opt in '-c': # crop atlas
    #         crop = 1
    #     if opt == '-h': # help option
    #         usage(label_title, label_name, label_num,fname_tracts) # display usage
    #     elif opt in '-t': # Mean true value in the generated phantom
    #         true_value = float(arg)
    #     elif opt in '-n': # standard deviation of the noise added to the generated phantom
    #         np.std_noise = float(arg)
    #     elif opt in '-s': # range of the random values given to the tracts
    #         range_tract = float(arg)
    #     elif opt in '-p': # Correction of the influence of partial volume on the estimations
    #         partial_vol_corr = 1
    #     elif opt in '-f': # folder where the WM atlas is located
    #         folder_cropped_atlas = os.path.abspath(arg)
    #     elif opt in '-g': # name of the generated phantom
    #         generated_phantom = str(arg)
    #     elif opt in '-m': # name of the manual_mask
    #         mask_prefix = str(arg)
    #     elif opt in '-d': # name of the different manual_masks
    #         mask_prefix = str(arg)
    #         mask_prefix = mask_prefix.split(',')
    #     elif opt in '-b': # number of bootstrapping iterations
    #         bootstrap_iterations = int(arg)
    #     elif opt in '-l': # labels for which results are displayed
    #         spec_tracts = arg
    #         spec_tracts = spec_tracts.split(',')
    #     elif opt in '-r': # name of the text file where all results are saved
    #         results_file = str(arg)
    #     elif opt in '-z': # z index for which slices are extract in the atlas to create the cropped atlas
    #         zcrop_ind = arg
    #     else: # verify that all entries are correct
    #         print('\nERROR: Option {} unknown. Exit program.\n'.format(opt))
    #         sys.exit(2) # exit program
    
    # Crop the atlas
    if (crop == 1):
        sct.run('./crop_atlas.py -f '+folder_atlas+' -o '+folder_cropped_atlas+' -z '+str(zcrop_ind))
        # Copy the info_label.txt file in the cropped atlas' folder
        # This file needs to be there in order for the sct_extract_metric code to work
        sct.run('cp ../WM_atlas_generation/info_label.txt '+folder_cropped_atlas)
    else:
        folder_cropped_atlas = folder_atlas

    # Extract the tracts from the atlas' folder
    tracts = get_tracts(folder_cropped_atlas)
    numtracts = len(tracts)

    # Get ponderation of each tract for dorsal colum average
    # ponderation of each tract of the dorsal column
    pond_dc = np.zeros(numtracts)
    # sum of each 
    pond_sum = 0
    for i in dorsal_column_labels.split(','):
        i = int(i)
        # Sum tracts values which are higher than 0 in the tracts
        pond_dc[i]= sum(tracts[i, 0][tracts[i, 0]>0])
        pond_sum = pond_sum + pond_dc[i]
    # Normalize the sum of ponderations to 1 
    pond_dc = pond_dc / pond_sum

    # Initialize values of estimations
    [X_map, X_map_dc, D_map, D_map_dc] = init_values(numtracts, bootstrap_iterations)
    [X_wa, X_wa_dc, D_wa, D_wa_dc] = init_values(numtracts, bootstrap_iterations)
    [X_bin, X_bin_dc, D_bin, D_bin_dc] = init_values(numtracts, bootstrap_iterations)
    [X_man, X_man_dc, D_man, D_man_dc] = init_values(numtracts, bootstrap_iterations)
    [X_wath, X_wath_dc, D_wath, D_wath_dc] = init_values(numtracts, bootstrap_iterations)

    man_users_number = len(mask_prefix)
    X_man = np.zeros([man_users_number,numtracts+1])
    X_man_dc = np.zeros([man_users_number,bootstrap_iterations])
    D_man = np.zeros([man_users_number,numtracts,bootstrap_iterations])
    D_man_dc = np.zeros([man_users_number,1,bootstrap_iterations])

    # loop across bootstrap
    for i in range(0,bootstrap_iterations):
        print 'iteration :  ' + str(i+1) + '/' + str(bootstrap_iterations)   
        X_map = np.zeros([numtracts])
        X_wa = np.zeros([numtracts])
        X_bin = np.zeros([numtracts])
        X_wath = np.zeros([numtracts])
        X_man = np.zeros([man_users_number, numtracts+1])
        
        # Generate phantom
        [WM_phantom, WM_phantom_noise, values_synthetic_data, tracts_sum] = phantom_generation(tracts, std_noise, range_tract, true_value)
        # Save generated phantoms as nifti image (.nii.gz)
        save_3D_nparray_niftii(WM_phantom, generated_phantom)
        save_3D_nparray_niftii(WM_phantom_noise, generated_phantom_noise)
        
        if partial_vol_corr == 1:
            save_3D_nparray_niftii(tracts_sum, tracts_sum_img)
            # Creation of an inverse phantom
            # Substract 1 from the sum of tracts
            sct.run('fslmaths ' + tracts_sum_img + ' -sub 1 temp.nii.gz')
            # Multiply this image by -true_value
            sct.run('fslmaths temp.nii.gz -mul -' + str(true_value) + ' temp.nii.gz')
            # Add this image to the phantom with added noise so that the effect of partial volume is lowered
            sct.run('fslmaths ' + generated_phantom_noise + ' -add temp.nii.gz ' + generated_phantom_noise)
        
        # Get the np.mean of all values in dorsal column in the generated phantom
        dc_val_avg = 0
        for j in dorsal_column_labels.split(','):
            j = int(j)
            dc_val_avg = dc_val_avg + values_synthetic_data[j] * pond_dc[j]  
        dc_val_avg = float(dc_val_avg)
    
        # Perform maximum likelihood estimation in all tracts
        sct.run('sct_extract_metric -i ' + generated_phantom_noise + ' -f ' + folder_cropped_atlas+ ' -m ml ')
        X_map = read_results(metrics_estimation_results)
        # Get the absolute deviation with the real value in the phantom
        D_map[:,i] = abs(X_map.ravel() - values_synthetic_data)
        # Get metrics estimation in dorsal column
        sct.run('sct_extract_metric -i ' + generated_phantom_noise + ' -f ' + folder_cropped_atlas+ ' -m ml -l '+dorsal_column_labels +' -a' )
        X_map_dc = read_results(metrics_estimation_results)
        # Get the difference between this np.mean and the metrics estimated
        D_map_dc[0, i] = abs(X_map_dc - dc_val_avg)
        
        # Perform binary estimation
        sct.run('sct_extract_metric -i ' + generated_phantom_noise + ' -f ' + folder_cropped_atlas + ' -m bin')
        X_bin = read_results(metrics_estimation_results)
        D_bin[:,i] = abs(X_bin.ravel() - values_synthetic_data)
        # Get results in dorsal column
        sct.run('sct_extract_metric -i ' + generated_phantom_noise + ' -f ' + folder_cropped_atlas+ ' -m bin -l '+dorsal_column_labels +' -a' )
        X_bin_dc = read_results(metrics_estimation_results)
        D_bin_dc[0, i] = abs(X_bin_dc - dc_val_avg)
        
        # Perform weighted average estimation
        sct.run('sct_extract_metric -i ' + generated_phantom_noise + ' -f ' + folder_cropped_atlas + ' -m wa')
        X_wa = read_results(metrics_estimation_results)
        # Get the absolute deviation with the real value in the phantom
        D_wa[:,i] = abs(X_wa.ravel() - values_synthetic_data)
        # Get results in dorsal column
        sct.run('sct_extract_metric -i ' + generated_phantom_noise + ' -f ' + folder_cropped_atlas+ ' -m wa -l '+dorsal_column_labels +' -a' )
        X_wa_dc = read_results(metrics_estimation_results)
        D_wa_dc[0, i] = abs(X_wa_dc - dc_val_avg)
        
        # Perform thresholded weighted average estimation
        sct.run('sct_extract_metric -i ' + generated_phantom_noise + ' -f ' + folder_cropped_atlas + ' -m wath')
        X_wath = read_results(metrics_estimation_results)
        # Get the absolute deviation with the real value in the phantom
        D_wath[:,i] = abs(X_wath.ravel() - values_synthetic_data)
        # Get results in dorsal column
        sct.run('sct_extract_metric -i ' + generated_phantom_noise + ' -f ' + folder_cropped_atlas+ ' -m wath -l '+dorsal_column_labels +' -a' )
        X_wath_dc = read_results(metrics_estimation_results)
        D_wath_dc[0, i] = abs(X_wath_dc - dc_val_avg)

        # Manual estimation
        # Quantify image within mask
        # Perform a weighted average for all nonzero indices from the mask
        for l in range(0, man_users_number):
            for k in man_mask_index:
                if k < 10:
                    header_mask = mask_folder[l] + mask_prefix[l] + '0' + str(k) + mask_ext
                else:
                    header_mask = mask_folder[l] + mask_prefix[l] + str(k) + mask_ext                        
                status, output = sct.run('sct_average_data_within_mask -i ' + generated_phantom_noise + ' -m ' + header_mask + ' -v 0')
                X_man[l,k] = float(output)
                if k != dorsal_column_mask_index:
                    D_man[l,k,i] = abs(X_man[l,k] - values_synthetic_data[k])
                else:
                    D_man_dc[l, 0, i] = abs(X_man[l,k] - dc_val_avg)     

    # Calculate elapsed time
    elapsed_time = int(round(time.time() - start_time))

    # Extract time in minutes and seconds
    sec = elapsed_time % 60
    mte = (elapsed_time - sec) / 60

    # Open text file where results are printed
    results_text = open(results_file, 'w+')
    # Display time in seconds and minutes for more than one minute
    if mte != 0:
        print >>results_text, '\nFinished! Elapsed time: ' + str(mte) + 'min' + str(sec) + 's'

    # Display time in seconds for less than one minute
    else:
        print >>results_text, '\nFinished! Elapsed time: ' + str(sec) + 's'

    # Display parameters used in results_file
    print >>results_text, '\tsigma noise: ' + str(std_noise) + '% \trange tracts: (-' + str(range_tract) + '%:+' + str(range_tract) + '%)\ttrue_value: ' + str(true_value) + \
    '\n\tnumber of iterations: ' + str(bootstrap_iterations) + '\t\tpartial volume correction: ' + str(partial_vol_corr)
    # Print the np.mean absolute deviation and its associated standard deviation for each method

    print >>results_text,'=================================================================================================='
    print >>results_text,'                        Mean absolute error for each method                                   '
    print >>results_text,'=================================================================================================='
    print >>results_text,'\tLabel \t    ml\t\t    wa\t\t    wath \t  bin \t       '+  mask_prefix[0] + '  \t' + mask_prefix[1] + '\t' + mask_prefix[2] + '\t' + mask_prefix[3]

    for i in spec_tracts:
        print >>results_text,'\t' + str(i) +'\t' + str(round(np.mean(D_map[i,:]),3)) + '+/-' +  str(round(np.std(D_map[i,:]),3)) + \
         '\t' + str(round(np.mean(D_wa[i,:]),3)) + '+/-' +  str(round(np.std(D_wa[i,:]),3)) + \
         '\t' + str(round(np.mean(D_wath[i,:]),3)) + '+/-' +  str(round(np.std(D_wath[i,:]),3)) + \
         '\t' + str(round(np.mean(D_bin[i,:]),3)) + '+/-' +  str(round(np.std(D_bin[i,:]),3)) + \
         '\t' + str(round(np.mean(D_man[0, i,:]),3)) + '+/-' +  str(round(np.std(D_man[0, i,:]),3)) + \
         '\t' + str(round(np.mean(D_man[1, i,:]),3)) + '+/-' +  str(round(np.std(D_man[1, i,:]),3)) + \
         '\t' + str(round(np.mean(D_man[2, i,:]),3)) + '+/-' +  str(round(np.std(D_man[2, i,:]),3)) + \
         '\t' + str(round(np.mean(D_man[3, i,:]),3)) + '+/-' +  str(round(np.std(D_man[3, i,:]),3)) 		 

    print >>results_text,'  dorsal_column' +'\t' + str(round(np.mean(D_map_dc[0, :]),3)) + '+/-' +  str(round(np.std(D_map_dc[0,:]),3)) + \
     '\t' + str(round(np.mean(D_wa_dc[0, :]),3)) + '+/-' +  str(round(np.std(D_wa_dc[0, :]),3)) +\
     '\t' + str(round(np.mean(D_wath_dc[0, :]),3)) + '+/-' +  str(round(np.std(D_wath_dc[0, :]),3)) +\
     '\t' + str(round(np.mean(D_bin_dc[0, :]),3)) + '+/-' +  str(round(np.std(D_bin_dc[0, :]),3)) +\
     '\t' + str(round(np.mean(D_man_dc[0, 0, :]),3)) + '+/-' +  str(round(np.std(D_man_dc[0, 0, :]),3)) +\
     '\t' + str(round(np.mean(D_man_dc[1, 0, :]),3)) + '+/-' +  str(round(np.std(D_man_dc[1, 0, :]),3)) +\
     '\t' + str(round(np.mean(D_man_dc[2, 0, :]),3)) + '+/-' +  str(round(np.std(D_man_dc[2, 0, :]),3)) +\
     '\t' + str(round(np.mean(D_man_dc[3, 0, :]),3)) + '+/-' +  str(round(np.std(D_man_dc[3, 0, :]),3)) 	 

    results_text.close()
    status, output = sct.run('cat ' + results_file)
    print output


def create_folder(folder):
    """create folder-- delete if already exists"""
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.mkdir(folder)


def read_results(fname_metrics):
    # Read file
    f = open(fname_metrics)

    # Extract all lines in the results file from sct_extract_metric
    # Do not extract lines which start with #
    lines = [lines for lines in f.readlines() if lines.strip() if not lines.startswith("#")]

    # read each line
    metrics_results = []
    for i in range(0, len(lines)):
        line = lines[i].split(',')
        # Get np.mean value of metric from column 2 of the text file
        metrics_results.append(line[2][:-1].replace(" ", ""))
    # Transform the column into a numpy array
    metrics_results = np.array(metrics_results)
    # Change the type of the values in the numpy array to float
    metrics_results = metrics_results.astype(np.float)
    return metrics_results


def init_values(number_of_tracts, number_of_iterations):
    # Mean estimation in tract
    X_ = np.zeros([number_of_tracts])
    # Mean estimation in dorsal column
    X_dc = np.zeros([1, 1])
    # Mean absolute error between np.mean estimation and true value in a tract
    D_ = np.zeros([number_of_tracts,number_of_iterations])
    # Mean absolute error between np.mean estimation and true value in dorsal column
    D_dc = np.zeros([1, number_of_iterations])
    return [X_, X_dc, D_, D_dc]


if __name__ == "__main__":
    main()
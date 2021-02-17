#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 16:43:33 2021

@author: labuser
"""
import sys
import os
from math import floor
from numpy import savez_compressed
import numpy as np
import _pickle as cPickle
import pandas as pd

def _check_line(line,
                contig_idx,
                position_idx,
                reference_kmer_idx,
                model_kmer_idx,
                samples_idx):
    """
    Do a regex search against all defined regexes and
    return the key and match result of the first matching regex

    """
    line_split = line.rstrip().split('\t')
    # check model and reference are the same
    if line_split[reference_kmer_idx] != line_split[model_kmer_idx]:
        return None
    # check the model is not NNNNN
    if line_split[model_kmer_idx] == 'NNNNN':
        return None
    # check if there is any A in the model
    if 'A' not in line_split[model_kmer_idx]:
        return None
    
    return line_split


def _parse_kmers(checked_line,
                 contig_idx,
                 position_idx,
                 reference_kmer_idx,
                 model_kmer_idx,
                 samples_idx,
                 file_object):
    """
    """
    kmer_lines = {}
    position_1 = int(checked_line[position_idx])
    positions = [position_1, 
                 position_1+1,
                 position_1+2,
                 position_1+3,
                 position_1+4]

    while len(kmer_lines.keys()) < 5:
        
        samples = [float(i) for i in checked_line[samples_idx].split(',')]
        
        if int(checked_line[position_idx]) in kmer_lines:
            kmer_lines[int(checked_line[position_idx])] = [checked_line[model_kmer_idx], 
                                                          kmer_lines[int(checked_line[position_idx])][1] +\
                                                          samples,
                                                          checked_line[contig_idx]]
        else:
            kmer_lines[int(checked_line[position_idx])] = [checked_line[model_kmer_idx], 
                                                           samples,
                                                           checked_line[contig_idx]]
        line = file_object.readline()
        if line == '':
            break
        checked_line = _check_line(line, 
                                   contig_idx,
                                   position_idx,
                                   reference_kmer_idx,
                                   model_kmer_idx,
                                   samples_idx)
        
    if list(kmer_lines.keys()) == positions:
        return kmer_lines
    else:
        return False


def _smooth_kmer(parsed_kmer, model_kmer_dict, lenght_event):
    '''
    Smooth the signals to fix the lenght
    '''
    kmer_name = sorted(parsed_kmer.items())[0][1][0]+\
                sorted(parsed_kmer.items())[1][1][0][-1]+\
                sorted(parsed_kmer.items())[2][1][0][-1]+\
                sorted(parsed_kmer.items())[3][1][0][-1]+\
                sorted(parsed_kmer.items())[4][1][0][-1]
    
    id_kmer = list(parsed_kmer.values())[0][-1]+'_'+str(sorted(parsed_kmer)[0])+'_'+kmer_name
    
    signal_smoothed = []
    
    for pos in sorted(parsed_kmer):
        event = parsed_kmer[pos][1]
        event_smoothed = smooth_event(event, lenght_event) # smooth the event
        signal_smoothed += event_smoothed  # add the event to the signal
    
    # create an expected signal according to the kmer
    expected_smoothed = make_expected(model_kmer_dict, 
                                      kmer_name,
                                      lenght_event)
    # claculate distance between expected and actual signal
    distance_vector = distance_calculator(expected_smoothed,
                                          signal_smoothed)
    
    return signal_smoothed, distance_vector, id_kmer
    

def make_expected(model_kmer_dict, kmer, event_lenght):
    '''
    '''
    expected_signal = []
    for i in range(5):
        expected_signal += [model_kmer_dict[kmer[i:i+5]]]*event_lenght
    return expected_signal


def distance_calculator(signal_expected, event_smoothed):
    '''
    '''
    vector_distance = list(np.round(abs(np.array(signal_expected) - \
                                        np.array(event_smoothed)), 3))    
    return vector_distance

def smooth_event(raw_signal, lenght_events):
    '''
    smmoth the signal 
    '''
    raw_signal_events = []
    
    if len(raw_signal) < lenght_events:
        event = top_median(raw_signal, lenght_events)
        raw_signal_events = [round(i, 3) for i in event]
        
    else:
        division = floor(len(raw_signal)/lenght_events)
        new_event = []
        for i in range(0, len(raw_signal), division):
            new_event.append(np.median(raw_signal[i:i+division]))
            if len(new_event) == lenght_events:
                break
        if len(new_event) < lenght_events:
            new_event = top_median(new_event, lenght_events)
        raw_signal_events = [round(i, 3) for i in new_event]
    return raw_signal_events

def top_median(array, lenght):
    '''
    This function top an array until some specific lenght
    '''
    extra_measure = [np.median(array)]*(lenght-len(array))
    array += extra_measure
    return array


def _combine_vectors(smooth_signal,
                     smooth_distance,
                     ):
    '''
    combine signals and distance vectors
    '''
    events = np.array(smooth_signal).reshape(len(smooth_signal), 1)
    distances = np.array(smooth_distance).reshape(len(smooth_distance), 1)

    combined = np.concatenate((events,
                               distances), 
                               axis=1)
    return combined


def parse_nanopolish(filepath, model_kmer_dict, lenght_event, directory_out):
    """
    Parse nanopolish

    Parameters
    ----------
    filepath : str
        Filepath 

    Returns
    -------
    data :

    """
    # open the file and read through it line by line
    grab_kmer = None
    line_counter = 0
    with open(filepath, 'r') as file_object:
        line = file_object.readline()
        # check the header
        header = line.rstrip().split('\t')
        if header[-1] != 'samples':
            print('nanopolish samples are not found, please run nanopolish with flag --samples')
            sys.exit()
        # get the column
        try:
            contig_idx = header.index('contig')
            position_idx = header.index('position')
            reference_kmer_idx = header.index('reference_kmer')
            model_kmer_idx = header.index('model_kmer')
            samples_idx = header.index('samples')
        except:
            print('Some nanopolish columns are not found')
    
        while line != '':  # The EOF char is an empty string
            line = file_object.readline()
            if line == '':
                break
            # check line is fordward and does not have NNNNN in the model
            checked_line = _check_line(line, 
                                       contig_idx,
                                       position_idx,
                                       reference_kmer_idx,
                                       model_kmer_idx,
                                       samples_idx)
            if checked_line:
                # If the kmer 
                if checked_line[model_kmer_idx][-1] == 'A':
                    parsed_kmer = _parse_kmers(checked_line,
                                               contig_idx,
                                               position_idx,
                                               reference_kmer_idx,
                                               model_kmer_idx,
                                               samples_idx,
                                               file_object)
                    
                    if parsed_kmer:
                        smooth_signal, smooth_distance, ID = _smooth_kmer(parsed_kmer,
                                                                          model_kmer_dict,
                                                                          lenght_event)
                        
                        combined_signals = _combine_vectors(smooth_signal,
                                                            smooth_distance,
                                                            )
                        
                        with open(directory_out+'/'+'signals.P', "ab") as sig_out:
                            cPickle.dump(combined_signals, sig_out)
                            
                        with open(directory_out+'/'+'IDs.P', "ab") as id_out:
                            cPickle.dump(ID, id_out)
                            
                        
    return True
                                                  
                        
                    
                    
                    
                    
                    
        
if __name__ == '__main__':
    
    
    nanopolish = sys.argv[1]
    filepath = '/media/labuser/Data/nanopore/m6A_classifier/data/yeast/nanopolish_reads/head_200_nanopolish.txt'
                 
    #filepath =  '/home/pablo/lib/MILONGAS/m6Aclassifier/nanopolish_test.txt'
    
    model_kmer_path = '/media/labuser/Data/nanopore/xpore/xpore/diffmod/model_kmer.csv'
    model_kmer = pd.read_csv(model_kmer_path,
                             sep=',')
    
    # create a dictionary with each kmer and its current value
    model_kmer_dict = dict(zip(model_kmer['model_kmer'], model_kmer['model_mean']))
    lenght_event = 20
    
    directory_out = '/media/labuser/Data/nanopore/m6A_classifier/data/yeast/nanopolish_reads/test_MILOGNAS_preprocess'
    # create directory if it does not exits
    if os.path.exists(directory_out) is False:
        os.makedirs(directory_out)
        
    data = parse_nanopolish(filepath, 
                            model_kmer_dict, 
                            lenght_event, 
                            directory_out)
    
    
    ### load IDs from file
    with open(directory_out+'/'+'signals.P', 'rb') as signal_in:
        with open(directory_out+'/'+'IDs.P', 'rb') as id_in:
            while True:
                try:
                    event_id = cPickle.load(id_in)
                    event_signal = cPickle.load(signal_in)
                    print(event_id, event_signal)
                except:
                    print('All signals have been processed')
                    break
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
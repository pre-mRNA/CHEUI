#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 10:16:37 2021

@author: labuser
"""

import argparse

parser = argparse.ArgumentParser(prog='predict_MILONGAS v0.1', description=
                                 """ 
                                 This script takes a ID and signal files and generate predictions every motif with and A
                                 
                                 """, usage='python predict_MILONGAS.py -s <path_to_signals_file> -i <path_to_ID_file> '\
                                            '-m <path_to_DL_model> -o <file_out> \nversion: %(prog)s')

OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')


REQUIRED.add_argument("-s", "--signals_input",
                      help="path to the signal file",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-i", "--ids_input",
                      help="path to the IDs file",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-m", "--DL_model",
                      help="path to DL model",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-o", "--file_out",
                      help="Path to the output file",
                      metavar='\b',
                      required=True)

OPTIONAL.add_argument('-v', '--version', 
                        action='version', 
                        version='%(prog)s')                  

parser._action_groups.append(OPTIONAL)

ARGS = parser.parse_args()

# required arg
signals_input = ARGS.signals_input
ids_input = ARGS.ids_input
DL_model = ARGS.DL_model
file_out = ARGS.file_out


from tensorflow.keras import Input
from tensorflow.keras.models import Model
import _pickle as cPickle
from DL_models import build_Jasper
import pandas as pd
import numpy as np


# load the trainned model 
inputs = Input(shape=(100, 2))
output = build_Jasper(inputs,Deep=True)
model = Model(inputs=inputs, outputs=output)
    
# test the NN trainned with only 1 A
model.load_weights(DL_model)

### load the stored data 
counter = 0
IDs = []
signals = []
with open(signals_input, 'rb') as signal_in:
    with open(ids_input, 'rb') as id_in:
        while True:
            try:
                counter +=1
                IDs.append(cPickle.load(id_in))
                signals.append(cPickle.load(signal_in))
                # to avoid loading everything predict every 10k singnals
                if counter%10000 == 0:
                    print(counter, 'signals predicted')
                    predictions = model.predict(np.array(signals))
                    predictions_df = pd.DataFrame.from_dict({'KMER': IDs,
                                                            'Prediction': predictions.reshape(len(predictions)).tolist()}
                                                            )
                    
                    predictions_df.to_csv(file_out,
                                          mode='a',
                                          header=False,
                                          sep='\t', 
                                          index=False)
                    IDs = []
                    signals = []
                    predictions_df = pd.DataFrame({'KMER': [], 'Prediction': []})
            except:
                if IDs:
                    predictions = model.predict(np.array(signals))
                    predictions_df = pd.DataFrame.from_dict({'KMER': IDs,
                                                            'Prediction': predictions.reshape(len(predictions)).tolist()}
                                                            )
                    
                    predictions_df.to_csv(file_out,
                                          mode='a',
                                          header=False,
                                          sep='\t', 
                                          index=False)
                    IDs = []
                    signals = []
                    predictions_df = pd.DataFrame({'KMER': [], 'Prediction': []})
                print('All signals have been processed')
                break




















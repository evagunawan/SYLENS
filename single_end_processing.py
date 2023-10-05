#!/usr/bin/env python3

import logging
import random

from write_output_file import write_output_file

def process_single_end_sampling(argsRead1, argsRead2, argsSubsample, argsOutput, argsCompress, fastqDictionary1, argsSeed):
    
    output_object_2 = None
    input_file_name_2 = None

    #If no subsampling is desired
    if argsSubsample != None:

        logging.info('Downsampling beginning...')
        
        Random_IDs = random.sample(list(fastqDictionary1), argsSubsample)

        output_object_1 = [fastqDictionary1[info] for info in Random_IDs]

        input_file_name_1 = f'{argsSeed}_downsampled_{argsRead1}'

        logging.info('Writing to file...')

    #If subsampling is desired
    if argsSubsample == None:
        
        logging.info('No downsampling of single-end file.')
        
        All_IDs = list(fastqDictionary1)

        output_object_1 = [fastqDictionary1[info] for info in All_IDs]

        input_file_name_1 = f'non_downsampled_{argsRead1}'

        logging.info('Writing to file...')

    #Writes output file
    write_output_file(argsRead1, argsRead2, argsCompress, input_file_name_1, input_file_name_2, output_object_1, output_object_2, argsOutput)

    logging.info(f'Finished processing of {argsRead1}. Output is stored in {input_file_name_1}')
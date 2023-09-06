#!/usr/bin/env python3

import logging
import random

import read_fastq_file
from write_output_file import write_output_file

def process_single_end_sampling(argsRead1, argsRead2, argsSubsample, argsOutput, argsCompress):

    dictionary_Read1 = read_fastq_file.fastqDictionary1
    
    output_object_2 = None
    input_file_name_2 = None

    if argsSubsample != None:

        logging.info('Downsampling beginning...')
        
        Random_IDs = random.sample(list(dictionary_Read1), argsSubsample)

        output_object_1 = [dictionary_Read1[info] for info in Random_IDs]

        input_file_name_1 = f'downsampled_{argsRead1}'

        logging.info('Writing to file...')

    if argsSubsample == None:
        
        logging.info('No downsampling of single-end file.')
        
        All_IDs = list(dictionary_Read1)

        output_object_1 = [dictionary_Read1[info] for info in All_IDs]

        input_file_name_1 = f'non_downsampled_{argsRead1}'

        logging.info('Writing to file...')

    write_output_file(argsRead1, argsRead2, argsCompress, input_file_name_1, input_file_name_2, output_object_1, output_object_2, argsOutput)

    logging.info(f'Finished processing of {argsRead1}. Output is stored in {input_file_name_1}')
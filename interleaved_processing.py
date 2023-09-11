#!/usr/bin/env python3

import logging
import re 
import random

from write_output_file import write_interleaved_file

def process_interleaved_sampling(argsRead1, argsSubsample, argsOutput, argsCompress, fastqDictionary1, argsSeed, formatExpression):

    logging.info('Starting interleaved processing.')

    read1_list = []
    read2_list = []

    #Creating read 1 and read 2 lists to reference
    for element in list(fastqDictionary1):

        if re.search(formatExpression, element) and element.endswith('1'):

            read1_list.append(element)

        if re.search(formatExpression, element) and element.endswith('2'):

            read2_list.append(element)

    #Subsampling detected
    if argsSubsample != None:
        
        logging.debug('In interleaved with subsampling')

        #Creating a zipped list to randomly subsample both ids from 
        zipped_ids_list = list(zip(read1_list, read2_list))
        random_ids = random.sample(zipped_ids_list, argsSubsample)
        r1_ids, r2_ids = zip(*random_ids)
        
        #Storing subsampled IDs and info from dict into output objects 
        output_1 = [fastqDictionary1[info] for info in r1_ids]
        output_2 = [fastqDictionary1[info] for info in r2_ids]  

        #Creating output file name
        input_file_name = f'{argsSeed}_downsampled_{argsRead1}'

        logging.info('Writing to file...')

    #No subsampling
    if argsSubsample == None:

        logging.debug('In interleaved without subsampling')

        output_1 = [fastqDictionary1[info] for info in read1_list]
        output_2 = [fastqDictionary1[info] for info in read2_list]
        
        input_file_name = f'non_downsampled_{argsRead1}'

    write_interleaved_file(argsRead1, argsOutput, argsCompress, output_1, output_2, input_file_name)

    logging.info(f'Finished processing of {argsRead1}. Output is stored in {input_file_name}')
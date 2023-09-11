#!/usr/bin/env python3
import random
import logging
import re

from write_output_file import write_output_file

def process_paired_end_sampling(argsRead1, argsRead2, argsSubsample, argsOutput, argsCompress, fastqDictionary1, fastqDictionary2, formatExpression, argsSeed):

    logging.info('No errors detected. Analyzing paired end files.')

    read1_list = []
    read2_list = []

    #Creating read 1 and read 2 lists to reference
    for element in list(fastqDictionary1):

        if re.search(formatExpression, element) and element.endswith('1'):

            read1_list.append(element)

    for element in list(fastqDictionary2):

        if re.search(formatExpression, element) and element.endswith('2'):

            read2_list.append(element)

    #Subsampling detected
    if argsSubsample != None:

        #Creating a zipped list to randomly subsample both ids from 
        zipped_ids_list = list(zip(read1_list, read2_list))
        random_ids = random.sample(zipped_ids_list, argsSubsample)
        r1_ids, r2_ids = zip(*random_ids)

        #Storing subsampled IDs and info from dict into output objects 
        output_1 = [fastqDictionary1[info] for info in r1_ids]
        output_2 = [fastqDictionary2[info] for info in r2_ids]

        #Creating output file name
        input_file_name_1 = f'{argsSeed}_downsampled_{argsRead1}'
        input_file_name_2 = f'{argsSeed}_downsampled_{argsRead2}'

        input_file_name_1 = f'{argsSeed}_downsampled_{argsRead1}'

    #No subsampling    
    if argsSubsample == None:

        output_1 = [fastqDictionary1[info] for info in read1_list]
        output_2 = [fastqDictionary2[info] for info in read2_list]

        input_file_name_1 = f'non_downsampled_{argsRead1}'
        input_file_name_2 = f'non_downsampled_{argsRead2}'

    write_output_file(argsRead1, argsRead2, argsCompress, input_file_name_1, input_file_name_2, output_1, output_2, argsOutput)

    logging.info(f'Finished processing of {argsRead1} and {argsRead2}. Output is stored in {input_file_name_1} and {input_file_name_2}')    
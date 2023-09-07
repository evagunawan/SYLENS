#!/usr/bin/env python3
import random
import logging

import read_fastq_file
from write_output_file import write_output_file

def process_paired_end_sampling(argsRead1, argsRead2, argsSubsample, argsOutput, argsCompress):

    logging.info('No errors detected. Analyzing paired end files.')

    dictionary_Read1 = read_fastq_file.fastqDictionary1
    dictionary_Read2 = read_fastq_file.fastqDictionary2

    read1_list = []
    read2_list = []

    for element in list(dictionary_Read1):

        if element.endswith('1'):

            read1_list.append(element)

    for element in list(dictionary_Read2):

        if element.endswith('2'):

            read2_list.append(element)

    if argsSubsample != None:

        zipped_ids_list = list(zip(read1_list, read2_list))
        random_ids = random.sample(zipped_ids_list, argsSubsample)
        r1_ids, r2_ids = zip(*random_ids)

        output_1 = [dictionary_Read1[info] for info in r1_ids]
        output_2 = [dictionary_Read2[info] for info in r2_ids]

        input_file_name_1 = f'downsampled_{argsRead1}'
        input_file_name_2 = f'downsampled_{argsRead2}'
    
    if argsSubsample == None:

        output_1 = [dictionary_Read1[info] for info in read1_list]
        output_2 = [dictionary_Read2[info] for info in read2_list]

        input_file_name_1 = f'non_downsampled_{argsRead1}'
        input_file_name_2 = f'non_downsampled_{argsRead2}'

    write_output_file(argsRead1, argsRead2, argsCompress, input_file_name_1, input_file_name_2, output_1, output_2, argsOutput)
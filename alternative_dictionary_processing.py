#! /usr/bin/env python3

import logging
import re
import sys
import gzip

from Bio import SeqIO


def process_alternative_dictionary(argsRead1, argsRead2, argsFiletype, format_dictionary, format_dictionary_2):

    def check_dict_format(format_dictionary, format_dictionary_2, first_ID):

        for pattern in format_dictionary:

            if re.search(pattern, first_ID):

                formatExpression = pattern

                format = format_dictionary[pattern]

                logging.debug(f'Regular expression found fastq ID format is {formatExpression}')

                logging.info(f'Fastq Formatting type is: {format}')

        for pattern in format_dictionary_2:

            if re.search(pattern, first_ID):

                formatExpression = pattern

                format = format_dictionary_2[pattern]

                logging.debug(f'Regular expression found fastq ID format is {formatExpression}')

                logging.info(f'Fastq Formatting type is: {format}')

        return format, formatExpression
 
    def create_original_dict(argsRead2, argsFiletype):

        if argsRead2.endswith('.gz'):

            with gzip.open(argsRead2, 'rt') as infile: 

                fastqDictionary2 = SeqIO.to_dict(SeqIO.parse(infile, argsFiletype)) 

        else:
        
            with open(argsRead2, 'rt') as infile:        

                fastqDictionary2 = SeqIO.to_dict(SeqIO.parse(infile, argsFiletype))

        return fastqDictionary2

    logging.debug('Entered into process_alternative_dictionary')

    if argsRead2 == None:

        logging.debug('No read 2')

        fastqDictionary1 = SeqIO.to_dict(SeqIO.parse(argsRead1, argsFiletype), key_function = lambda rec : rec.description)
        fastqDictionary2 = {None:None, None:None}

        first_ID_1 = list(fastqDictionary1) [0]
        last_ID_1 = list(fastqDictionary1) [-1]

        first_ID_2 = list(fastqDictionary2) [0]
        last_ID_2 = list(fastqDictionary2) [-1]

        logging.debug(f'Formatting new dictionary with new first read 1: {first_ID_1} and new last read 1: {last_ID_1}')

        for pattern in format_dictionary:

            if re.search(pattern, first_ID_1):

                formatExpression = pattern

                format = format_dictionary[pattern]

                logging.debug(f'Regular expression found fastq ID format is {formatExpression}')

                logging.info(f'Fastq Formatting type is: {format}')

        for pattern in format_dictionary_2:

            if re.search(pattern, first_ID_1):

                logging.critical('File appears to be a read 2 file instead of read 1 or interleaved file. Program terminating...')

                sys.exit(1)                    

    if argsRead2 != None:

        logging.debug('Read 2 detected, creating new dictionary for read 1')

        fastqDictionary1 = SeqIO.to_dict(SeqIO.parse(argsRead1, argsFiletype), key_function = lambda rec : rec.description)
        
        first_ID_1 = list(fastqDictionary1) [0]
        last_ID_1 = list(fastqDictionary1) [-1]

        logging.debug(f'Formatting new dictionary with new first read 1: {first_ID_1} and new last read 1: {last_ID_1}')

        format1, formatExpression1 = check_dict_format(format_dictionary, format_dictionary_2, first_ID_1)

        logging.debug('Read 2 detected, creating traditional dictionary and checking if format is not found.')
        
        fastqDictionary2 = create_original_dict(argsRead2, argsFiletype)

        first_ID_2 = list(fastqDictionary2) [0]
        last_ID_2 = list(fastqDictionary2) [-1]

        format2, formatExpression2 = check_dict_format(format_dictionary, format_dictionary_2, first_ID_2)

        if format1 == format2:

            pass

        else:

            fastqDictionary2 = SeqIO.to_dict(SeqIO.parse(argsRead2, argsFiletype), key_function = lambda rec : rec.description)

            first_ID_2 = list(fastqDictionary2) [0]
            last_ID_2 = list(fastqDictionary2) [-1]

            format2, formatExpression2 =  check_dict_format(format_dictionary, format_dictionary_2, first_ID_2)

            if format1 != format2:

                logging.critical('Files appear to be different formats. Program terminating...')

                sys.exit(1)
        
        logging.debug(f'Formatting new dictionary with new first read 2: {first_ID_2} and new last read 2: {last_ID_2}')
    
    return fastqDictionary1, fastqDictionary2, first_ID_1, last_ID_1, first_ID_2, last_ID_2, format, formatExpression
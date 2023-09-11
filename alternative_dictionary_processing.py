#! /usr/bin/env python3
import logging
import re

from Bio import SeqIO


def process_alternative_dictionary(argsRead1, argsRead2, argsFiletype, format_dictionary):

    if argsRead2 == None:

        fastqDictionary1 = SeqIO.to_dict(SeqIO.parse(argsRead1, argsFiletype), key_function = lambda rec : rec.description)

        first_ID_1 = list(fastqDictionary1) [0]
        last_ID_1 = list(fastqDictionary1) [-1]

        logging.debug(f'Formatting new dictionary with new first read 1: {first_ID_1} and new last read 1: {last_ID_1}')

        for pattern in format_dictionary:

            if re.search(pattern, first_ID_1):
                    
                formatExpression = pattern

                format = format_dictionary[pattern]

                logging.debug(f'Regular expression found fastq ID format is {formatExpression}')

                logging.info(f'Fastq Formatting type is: {format}')

                return fastqDictionary1, first_ID_1, last_ID_1, format, formatExpression

    if argsRead2 != None:

        fastqDictionary1 = SeqIO.to_dict(SeqIO.parse(argsRead1, argsFiletype), key_function = lambda rec : rec.description)

        fastqDictionary2 = SeqIO.to_dict(SeqIO.parse(argsRead2, argsFiletype), key_function = lambda rec : rec.description)
                    
        first_ID_1 = list(fastqDictionary1) [0]
        last_ID_1 = list(fastqDictionary1) [-1]

        first_ID_2 = list(fastqDictionary2) [0]
        last_ID_2 = list(fastqDictionary2) [-1]

        logging.debug(f'Formatting new dictionary with new first read 1: {first_ID_1} and new last read 1: {last_ID_1}')
        logging.debug(f'Formatting new dictionary with new first read 2: {first_ID_2} and new last read 2: {last_ID_2}')

        for pattern in format_dictionary:

            if re.search(pattern, first_ID_1):

                formatExpression = pattern

                format = format_dictionary[pattern]

                logging.debug(f'Regular expression found fastq ID format is {formatExpression}')

                logging.info(f'Fastq Formatting type is: {format}')

                return fastqDictionary1, fastqDictionary2, first_ID_1, last_ID_1, first_ID_2, last_ID_2, format, formatExpression
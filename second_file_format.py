#! /user/bin/env python3

import logging
import re
import sys

from Bio import SeqIO

from alternative_dictionary_processing import process_alternative_dictionary

def determine_second_file_format(argsRead2, argsFiletype, first_ID_2, format_dictinoary_1, format_dictionary_2, ID_1_format):

    completed = False

    logging.debug('Determining if second file and first file have same format')
    
    for pattern in format_dictionary_2:

        try:

            if re.search(pattern, first_ID_2):

                format = format_dictionary_2[pattern]

                completed = True

                break
    
        except TypeError:

            pass
    
    if completed == False:

        fastqDictionary2 = SeqIO.to_dict(SeqIO.parse(argsRead2, argsFiletype), key_function = lambda rec : rec.description)

        first_ID_2 = list(fastqDictionary2) [0]
    
    for pattern in format_dictinoary_1:

        if re.search(pattern, first_ID_2):

            format = format_dictinoary_1[pattern]

            break
    
    for pattern in format_dictionary_2:

        if re.search(pattern, first_ID_2):

            format = format_dictionary_2[pattern]

            break
    
    if ID_1_format == format:

        logging.debug('Both file formats are the same.')

    if ID_1_format != format:

        logging.critical(f'File formats are different. Program terminating...')

        sys.exit(1)
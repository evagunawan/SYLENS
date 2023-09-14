#! /user/bin/env python3

import logging
import re
import sys

from Bio import SeqIO

def determine_second_file_format(argsRead2, argsFiletype, first_ID_2, format_dictionary_1, format_dictionary_2, ID_1_format):

    format = None
    completed = False

    logging.debug('Determining if second file and first file have same format')
    
    for pattern in format_dictionary_2:

        if re.search(pattern, first_ID_2):

            format = format_dictionary_2[pattern]

            completed = True

            break

        for pattern in format_dictionary_1:

            if re.search(pattern, first_ID_2):

                format = format_dictionary_1[pattern]

                completed = True

                break

 
        if completed != True:

            fastqDictionary2 = SeqIO.to_dict(SeqIO.parse(argsRead2, argsFiletype), key_function = lambda rec : rec.description)

            first_ID_2 = list(fastqDictionary2) [0]
    
            for pattern in format_dictionary_1:

                if re.search(pattern, first_ID_2):

                    format = format_dictionary_1[pattern]

                    break
    
            for pattern in format_dictionary_2:

                if re.search(pattern, first_ID_2):

                    format = format_dictionary_2[pattern]

                    break
    
    if ID_1_format == format:

        logging.debug('Both file formats are the same.')

    if ID_1_format != format:

        logging.critical(f'File formats are different. First file format is {ID_1_format} but second file format is {format}. Program terminating...')

        sys.exit(1)
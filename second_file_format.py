#! /user/bin/env python3

import logging
import re
import sys

from Bio import SeqIO

#Writing function to determine file formats
def determine_second_file_format(argsRead2, argsFiletype, first_ID_2, format_dictionary_1, format_dictionary_2, ID_1_format):

    format = None
    format_expression_2 = None
    completed = False

    logging.debug('Determining if second file and first file have same format')
    
    #For every pattern in the read 2 dictionary, looks for a regular expression match to first_ID_2's first read
    for pattern in format_dictionary_2:

        if re.search(pattern, first_ID_2):

            format_expression_2 = pattern

            format = format_dictionary_2[pattern]

            completed = True

            break

    #If no format was found in format dictionary 2, looks in format dictionary 1
    if format == None:

        for pattern in format_dictionary_1:

            if re.search(pattern, first_ID_2):

                format_expression_2 = pattern

                format = format_dictionary_1[pattern]

                completed = True

                break

    #If a pattern is not found still, will process an alternative dictionary with key as description
    if completed != True:

        fastqDictionary2 = SeqIO.to_dict(SeqIO.parse(argsRead2, argsFiletype), key_function = lambda rec : rec.description)

        first_ID_2 = list(fastqDictionary2) [0]
    
        #Tries to determine which pattern first_ID_2 matches to in format dictionary 1
        for pattern in format_dictionary_1:

            if re.search(pattern, first_ID_2):

                format_expression_2 = pattern

                format = format_dictionary_1[pattern]

                break
    
        #Tries to determine which pattern first_ID_2 matches to in format dictionary 2        
        for pattern in format_dictionary_2:

            if re.search(pattern, first_ID_2):

                format_expression_2 = pattern
                
                format = format_dictionary_2[pattern]

                break
    
    if ID_1_format == format:

        logging.debug('Both file formats are the same.')

    if ID_1_format != format:

        logging.critical(f'File formats are different. First file format is {ID_1_format} but second file format is {format}. Program terminating...')

        sys.exit(1)
    
    return format_expression_2
#! /usr/bin/env python3

import logging
import re
import sys
import gzip

from Bio import SeqIO

#Gets the key from the dictionary and returns it
def get_key(format1, format_dictionary_2):

    for key, value in format_dictionary_2.items():

        if format1 == value:

            formatExpression2 = key

            return formatExpression2

def process_alternative_dictionary(argsRead1, argsRead2, argsFiletype, format_dictionary, format_dictionary_2):

    #Checks the dictionary format of the fastq file
    def check_dict_format(format_dictionary, format_dictionary_2, first_ID):

        #Looking through all the patterns in the dictionary 1
        for pattern in format_dictionary:

            #If pattern is found in the first ID, sets expression and returns the proper format and format expression
            if re.search(pattern, first_ID):

                formatExpression = pattern

                format = format_dictionary[pattern]

                logging.debug(f'Regular expression found fastq ID format is {formatExpression}')

                logging.info(f'Fastq Formatting type is: {format}')

                return format, formatExpression

        #Looking through all the patterns in dictionary 2
        for pattern in format_dictionary_2:

            if re.search(pattern, first_ID):

                formatExpression = pattern

                format = format_dictionary_2[pattern]

                logging.debug(f'Regular expression found fastq ID format is {formatExpression}')

                logging.info(f'Fastq Formatting type is: {format}')

                return format, formatExpression
 
    #Used if a dictionary needs to be created for incorrect read 1 and read 2 pairings i.e. casava read 1 paired with NCBI read 2
    def create_original_dict(argsRead2, argsFiletype):

        if argsRead2.endswith('.gz'):

            with gzip.open(argsRead2, 'rt') as infile: 

                fastqDictionary2 = SeqIO.to_dict(SeqIO.parse(infile, argsFiletype)) 

        else:
        
            with open(argsRead2, 'rt') as infile:        

                fastqDictionary2 = SeqIO.to_dict(SeqIO.parse(infile, argsFiletype))

        return fastqDictionary2

    logging.debug('Entered into process_alternative_dictionary')

    format2 = None
    formatExpression2 = None

    #If no read 2 is present
    if argsRead2 == None:
        
        logging.debug('No read 2')

        #Creates dictionary with description as Key and second dictionary is blank
        fastqDictionary1 = SeqIO.to_dict(SeqIO.parse(argsRead1, argsFiletype), key_function = lambda rec : rec.description)
        fastqDictionary2 = {None:None, None:None}

        #Gets first and last reads of each dictionary
        first_ID_1 = list(fastqDictionary1) [0]
        last_ID_1 = list(fastqDictionary1) [-1]

        first_ID_2 = list(fastqDictionary2) [0]
        last_ID_2 = list(fastqDictionary2) [-1]

        logging.debug(f'Formatting new dictionary with new first read 1: {first_ID_1} and new last read 1: {last_ID_1}')

        #Checks to see if format pattern of file matches a read 2 format, if so terminates 
        for pattern in format_dictionary_2:

            if re.search(pattern, first_ID_1):

                logging.critical('Required positional file appears to be a reverse read file. Positional file should be a forward read or interleaved file. Program terminating...')

                sys.exit(1)   

        #If no format pattern is found in dictionary 2, looks for pattern in dictionary 1
        for pattern in format_dictionary:

            if re.search(pattern, first_ID_1):

                formatExpression1 = pattern

                format1 = format_dictionary[pattern]

                logging.debug(f'Regular expression found fastq ID format is {formatExpression1}')

                logging.info(f'Fastq Formatting type is: {format1}')

        #Gets second format expression in case there is an interleaved file that was entered
        formatExpression2 = get_key(format1, format_dictionary_2)                 

    #If two files are uploaded
    if argsRead2 != None:

        logging.debug('Read 2 detected, creating new dictionary for read 1')

        #Creates a dictionary with description as the key
        fastqDictionary1 = SeqIO.to_dict(SeqIO.parse(argsRead1, argsFiletype), key_function = lambda rec : rec.description)
        
        #Gets first and last values of dictionary
        first_ID_1 = list(fastqDictionary1) [0]
        last_ID_1 = list(fastqDictionary1) [-1]

        logging.debug(f'Formatting new dictionary with new first read 1: {first_ID_1} and new last read 1: {last_ID_1}')

        #Checks to grab format and format expression from the first file
        format1, formatExpression1 = check_dict_format(format_dictionary, format_dictionary_2, first_ID_1)

        logging.debug('Read 2 detected, creating traditional dictionary and checking if format is not found.')

        #Exception to try and create an original dictionary with read 2. Will only pass if wrong filetype is present
        try:

            fastqDictionary2 = create_original_dict(argsRead2, argsFiletype)

            first_ID_2 = list(fastqDictionary2) [0]
            last_ID_2 = list(fastqDictionary2) [-1]

            #Gets format and format expression of second file
            format2, formatExpression2 = check_dict_format(format_dictionary, format_dictionary_2, first_ID_2)

        #If try fails, will create an alternative dictionary using description as key
        except:

            fastqDictionary2 = SeqIO.to_dict(SeqIO.parse(argsRead2, argsFiletype), key_function = lambda rec : rec.description)

            first_ID_2 = list(fastqDictionary2) [0]
            last_ID_2 = list(fastqDictionary2) [-1]

            format2, formatExpression2 =  check_dict_format(format_dictionary, format_dictionary_2, first_ID_2)

        #IF both formats are the same, the file is passed to the rest of the program
        if format1 == format2:

            pass

        #If the second file is a reverse file, terminates
        if formatExpression1.endswith('(2)'):

            logging.critical('Required positional file appears to be a reverse read file. Positional file should be a forward read or interleaved file. Program terminating...')
    
            sys.exit(1)

        #If file formates are different, program terminates
        if format1 != format2:

            logging.critical('Files appear to be different formats. Program terminating...')

            sys.exit(1)
        
        logging.debug(f'Formatting new dictionary with new first read 2: {first_ID_2} and new last read 2: {last_ID_2}')
    
        #If both format expressions are the same there are either two interleaved files or single end files. 
        if formatExpression1 == formatExpression2:

            logging.critical('Files appear to be the same file type I.E. two interleaved files or single end forward files. Program terminating...')

            sys.exit(1)

    return fastqDictionary1, fastqDictionary2, first_ID_1, last_ID_1, first_ID_2, last_ID_2, format1, formatExpression1, formatExpression2
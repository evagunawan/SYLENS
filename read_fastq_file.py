#!/usr/bin/env python3

#TODO Finish creating output class to output into main script? Think about if I truly need this class

import gzip
import logging
import re 

from Bio import SeqIO
from process_paired_end_mistakes import process_mistakes

#Function to create a dictionary, can be used when needed just by importing it
def create_dictionary(input_file, filetype):

    if input_file.endswith('.gz'):

        with gzip.open(input_file, 'rt') as infile: 

            fastqDictionary = SeqIO.to_dict(SeqIO.parse(infile, filetype)) 

    else:

        with open(input_file, 'rt') as infile:

            fastqDictionary = SeqIO.to_dict(SeqIO.parse(infile, filetype))
    
    return fastqDictionary

'''
Class to make a custom FastqFile object. This gives me the ability to make my own 
data blueprint. This object can then be applied throughout the main script
'''

class FastqFile:

    def __init__(self, path1, path2, filetype):

        self.path1 = path1
        self.path2 = path2
        self.filetype = filetype   
        
    def reading_fastq_file(self):

        global first_ID_1, last_ID_1
        global first_ID_2, last_ID_2
        global fastqDictionary1
        global fastqDictionary2

        logging.debug('Entered into reading_fastq')

        if self.path2 == None:
            
            fastqDictionary1 = create_dictionary(self.path1, self.filetype)
            
            first_ID_1 = list(fastqDictionary1) [0]
            last_ID_1 = list(fastqDictionary1) [-1]

            first_ID_2 = None
            last_ID_2 = None
            
            logging.debug(f'{first_ID_1} {last_ID_1}')
        
        else:
            
            fastqDictionary1 = create_dictionary(self.path1, self.filetype)
            fastqDictionary2 = create_dictionary(self.path2, self.filetype)
            
            first_ID_1 = list(fastqDictionary1) [0]
            last_ID_1 = list(fastqDictionary1) [-1]

            first_ID_2 = list(fastqDictionary2) [0]
            last_ID_2 = list(fastqDictionary2) [-1]

            logging.debug(f'Reads are : {first_ID_1} {last_ID_1} {first_ID_2} {last_ID_2}')
            
        print(determine_fastq_ID_formatting(first_ID_1))

        print(determine_paired_single_interleaved(first_ID_1, last_ID_1, first_ID_2, last_ID_2))

        fastqReadOutputObject = FastqReadOutput(determined_filetype, format)

        fastqReadOutputObject.output()

#Determining format of fastq file to properly figure out if R1 and/or R2
def determine_fastq_ID_formatting(ID):

    global format
    global formatExpression
    
    #Taking the 3 known formats and creating regular expressions with them. Stored in a dictinoary 
    # Illumina and Casava have same format, for now.
    format_dictionary = {
        '(.+)(\d) [1]' :'IlluminaAndCasava',
        '(^SRR)(\w+).(\d+).([1])': 'NCBI'
        }
    
    for pattern in format_dictionary:

        if re.search(pattern, ID):
            
            formatExpression = pattern

            format = format_dictionary[pattern]

            logging.debug(f'Regular expression found fastq ID format is {formatExpression}')

    logging.info(f'Fastq Formatting type is: {format}')

    return format, formatExpression

#Figures out if input file is single, paired, or interleaved
def determine_paired_single_interleaved(firstR1, lastR1, firstR2, lastR2):

    global determined_filetype

    if firstR2 == None and lastR2 == None:
        
        single_file = True
        paired_file = False

    if firstR2 != None and lastR2 != None:

        single_file = False
        paired_file = True

    if single_file == True:
    
        if firstR1.endswith('1') and lastR1.endswith('2'):

            determined_filetype = 'Interleaved'

            logging.info('File is an interleaved file.')

        if firstR1.endswith('1') and lastR1.endswith('1'):
                
            determined_filetype = 'Single-end'

            logging.info('File is a single-end file.')

    if paired_file == True:

        process_mistakes(firstR1, lastR1, firstR2, lastR2)

        if firstR1.endswith('1') and lastR1.endswith('1'):
            
            logging.debug('Forward file is correct.')

            if firstR2.endswith('2') and lastR2.endswith('2'):

                logging.debug('Reverse file is correct.')
 
                determined_filetype = 'Paired-end'

    return determined_filetype

#Creating class object to return final information 
class FastqReadOutput:

    def __init__(self, determined_filetype, determined_formatting):

        self.determined_filetype = determined_filetype
        self.determined_formatting = determined_formatting
    
    def output(self):
        
        logging.info('The filetype is: ' + self.determined_filetype)

        logging.info('The format is: ' + self.determined_formatting)

        return self.determined_filetype, self.determined_formatting
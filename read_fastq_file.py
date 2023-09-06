#!/usr/bin/env python3

import gzip
import logging
import re 

from Bio import SeqIO

logging.debug('Moving to read_fastq_file.py')


def reading_fastq_file(input_file, input_format):
    
    #Gets the path and filetype
    path = input_file
    filetype = input_format

    #This checks if compressed
    if path.endswith('.gz'):
        gzipped = True
        logging.debug('Ended with .gz')
        
    else:
        gzipped = False
        logging.debug('Not .gz')
        
    #This will read the data from the file into memory instead of creating a new file to read from
    if gzipped == True:
        logging.debug('Opening file.gz')
        with gzip.open(path, 'rt') as infile: 
            fastqDictionary = SeqIO.to_dict(SeqIO.parse(infile, filetype))   
    else:
        logging.debug('Opening file')
        with open(path, 'rt') as infile:
            fastqDictionary = SeqIO.to_dict(SeqIO.parse(infile, filetype))

    first_ID = list(fastqDictionary) [0]
    last_ID = list(fastqDictionary) [-1]
    
    return first_ID and last_ID
    

import sys
def verifying_file_type(read1, read2):

    logging.debug('Entering second definition')
    if read2 == None:
        if re.search(r"(^\w+.)(\w+)(.)(1$)", first_ID) and re.search(r"(^\w+.)(\w+)(.)(2$)", last_ID): 
            file_format = 'Interleaved'
        else:
            file_format = 'Single_end'
    if read2 != None:
        if re.search(r"(^\w+.)(\w+)(.)(1$)", read1.first_ID) and re.search(r"(^\w+.)(\w+)(.)(1$)", read1.last_ID):
            if re.search(r"(^\w+.)(\w+)(.)(1$)", read2.first_ID) and re.search(r"(^\w+.)(\w+)(.)(1$)", read2.last_ID):
                file_format = 'Paired_end'
    print(file_format)
    
    
    sys.exit()
    
    
    #Taking the 3 known formats and creating regular expressions with them. Illumina and Casava have same format, right now.
    seqIDformat = {
        'IlluminaAndCasava' : "r'^@(.+) (\d)'",
        'NCBI' : "r'^@(\w+).(\d).([12]) '"
    }

    for format in seqIDformat:
        if re.search(format, first_ID):
            format_type = seqIDformat
            print(format_type)
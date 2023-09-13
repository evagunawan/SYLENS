#!/usr/bin/env python3
import gzip
import logging
import re 
import sys

from Bio import SeqIO

from alternative_dictionary_processing import process_alternative_dictionary
from process_paired_end_mistakes import process_mistakes
from single_end_processing import process_single_end_sampling
from interleaved_processing import process_interleaved_sampling
from paired_end_processing import process_paired_end_sampling
from second_file_format import determine_second_file_format

#Function to create a dictionary, can be used when needed just by importing it

'''
Class to make a custom FastqFile object. This gives me the ability to make my own 
data blueprint. This object can then be applied throughout the main script
'''

class FastqFileData:

    def __init__(self, argsRead1, argsRead2, argsSubsample, argsOutput, argsCompress, argsFiletype, argsSeed):

        self.argsRead1 = argsRead1
        self.argsRead2 = argsRead2
        self.argsSubsample = argsSubsample
        self.argsOutput = argsOutput
        self.argsCompress = argsCompress 
        self.argsFiletype = argsFiletype   
        self.argsSeed = argsSeed
  
    def reading_fastq_file(self):

        def create_dictionary(input_file, argsFiletype):

            if input_file.endswith('.gz'):
            
                with gzip.open(input_file, 'rt') as infile: 
                
                    fastqDictionary = SeqIO.to_dict(SeqIO.parse(infile, argsFiletype)) 
    
            else:
     
                with open(input_file, 'rt') as infile:        

                    fastqDictionary = SeqIO.to_dict(SeqIO.parse(infile, argsFiletype))
    
            return fastqDictionary

        if self.argsRead2 == None:

            logging.info('Beginning to process file.')

            #Made exception. First, try creating a dictionary. With single illumina/casava files, dictionary will be created but will not have read number 
            try:

                self.fastqDictionary1 = create_dictionary(self.argsRead1, self.argsFiletype)
                self.fastqDictionary2 = {None:None, None:None}

            #In cases of interleaved illumina/casava, duplicate keys are created since read # is found in description, will instead create blank dictionaries 
            except ValueError:

                self.fastqDictionary1 = {None:None, None:None}
                self.fastqDictionary2 = {None:None, None:None}

        else:

            logging.info('Beginning to process files.')

            #Made exception. First try making a dict. If illumina/casava both dict will have same keys. Read number is found in description
            try:

                self.fastqDictionary1 = create_dictionary(self.argsRead1, self.argsFiletype)
                self.fastqDictionary2 = create_dictionary(self.argsRead2, self.argsFiletype)

            #In cases where incorrect files were uploaded, i.e. two casava read 2s, this creates blank dict to send to alternative dict process where problem can be identified.
            except ValueError:

                self.fastqDictionary1 = {None:None, None:None}
                self.fastqDictionary2 = {None:None, None:None}

        return self.fastqDictionary1, self.fastqDictionary2

    #Determining format of fastq file to properly figure out if R1 and/or R2
    def determine_fastq_ID_formatting(self):

        #Taking the top 3 known formats and creating regular expressions to identify IDs. Stored REs in a dictinoary 
        # Illumina and Casava have same format, for now.
        format_dictionary = {
            '(^SRR)(\w+)[.+](\d+)[.+](1)': 'NCBI',
            '(.+)(\d) (1)' :'IlluminaAndCasava'
            }

        format_dictionary_2 = {
            '(^SRR)(\w+)[.+](\d+)[.+](2)': 'NCBI',
            '(.+)(\d) (2)' :'IlluminaAndCasava'
            }

        self.formatExpression = None
        self.format = None

        completed = False

        self.first_ID_1 = list(self.fastqDictionary1) [0]
        self.last_ID_1 = list(self.fastqDictionary1) [-1]

        self.first_ID_2 = list(self.fastqDictionary2) [0]
        self.last_ID_2 = list(self.fastqDictionary2) [-1]        

        #Exception to first try finding a pattern in first_ID
        try:

            for pattern in format_dictionary:

                if re.search(pattern, self.first_ID_1):

                    self.formatExpression = pattern

                    self.format = format_dictionary[pattern]

                    logging.debug(f'Regular expression found fastq ID format is {self.formatExpression}')                

                    completed = True

                    #Keeps for loop from trying again if pattern is found
                    break

                else:

                    #Added an error to catch read 2 file instead of read 1 or interleaved file
                    for pattern in format_dictionary_2:

                        if re.search(pattern, self.first_ID_1):

                            logging.critical('File appears to be a read 2 file instead of read 1 or interleaved file. Program terminating...')

                            sys.exit(1)  

            #If a dictionary was created for a file but a pattern is not found because the read ID was found in the description, feeds into alternative dictionary
            if completed == False:

                self.fastqDictionary1, self.fastqDictionary2, self.first_ID_1, self.last_ID_1, self.first_ID_2, self.last_ID_2, self.format, self.formatExpression = process_alternative_dictionary(self.argsRead1, self.argsRead2, self.argsFiletype, format_dictionary, format_dictionary_2)

            if self.first_ID_2 != None:

                determine_second_file_format(self.argsRead2, self.argsFiletype, self.first_ID_1, format_dictionary, format_dictionary_2, self.format)

            return self.first_ID_1, self.first_ID_2, self.last_ID_1, self.last_ID_2, self.format, self.formatExpression

        #Used if dictionary is made from an interleaved illumina/casava file has no matches which creates a TypeError 
        except TypeError:

            logging.debug('TypeError path')

            self.fastqDictionary1, self.fastqDictionary2, self.first_ID_1, self.last_ID_1, self.first_ID_2, self.last_ID_2, self.format, self.formatExpression = process_alternative_dictionary(self.argsRead1, self.argsRead2, self.argsFiletype, format_dictionary, format_dictionary_2)
 
            if self.first_ID_2 != None:

                determine_second_file_format(self.argsRead2, self.argsFiletype, self.first_ID_1, format_dictionary, format_dictionary_2, self.format)


            return self.first_ID_1, self.first_ID_2, self.last_ID_1, self.last_ID_2, self.format, self.formatExpression

    #Figures out if input file is single, paired, or interleaved
    def determine_paired_single_interleaved(self):

        if self.format == 'IlluminaAndCasava':

            self.formatExpression2 = '(.+)(\d) (2)'

        if self.format == 'NCBI':

            self.formatExpression2 = '(^SRR)(\w+)[.+](\d+)[.+](2)'

        if self.argsRead2 == None:

            self.single_file = True
            self.paired_file = False

        if self.argsRead2 != None:

            self.single_file = False
            self.paired_file = True

        if self.single_file == True:

            if re.search(self.formatExpression, self.first_ID_1) and re.search(self.formatExpression2, self.last_ID_1):

                self.determined_filetype = 'Interleaved'

                logging.info('File is an interleaved file.')

            if re.search(self.formatExpression, self.first_ID_1) and re.search(self.formatExpression, self.last_ID_1):

                self.determined_filetype = 'Single-end'

                logging.info('File is a single-end file.')

            if re.search(self.formatExpression2, self.first_ID_1):
                
                logging.critical('Read 2 was supplied sintead of read 1 or interleaved file. Program terminating...')
                sys.exit(1)

        if self.paired_file == True:

            process_mistakes(self.first_ID_1, self.last_ID_1, self.first_ID_2, self.last_ID_2, self.formatExpression, self.formatExpression2)

            if re.search(self.formatExpression, self.first_ID_1) and re.search(self.formatExpression, self.last_ID_1):

                logging.debug('Forward file is correct.')

                if re.search(self.formatExpression2, self.first_ID_2) and re.search(self.formatExpression2, self.last_ID_2):

                    logging.debug('Reverse file is correct.')
    
                    self.determined_filetype = 'Paired-end'

        return self.determined_filetype

    #Writing output of determined filetype
    def processing_filetype(self):

        if self.determined_filetype == 'Single-end':

            process_single_end_sampling(self.argsRead1, self.argsRead2, self.argsSubsample, self.argsOutput, self.argsCompress, self.fastqDictionary1, self.argsSeed)

        elif self.determined_filetype == 'Interleaved':

            process_interleaved_sampling(self.argsRead1, self.argsSubsample, self.argsOutput, self.argsCompress, self.fastqDictionary1, self.argsSeed, self.formatExpression, self.formatExpression2)

        elif self.determined_filetype == 'Paired-end':

            process_paired_end_sampling(self.argsRead1, self.argsRead2, self.argsSubsample, self.argsOutput, self.argsCompress, self.fastqDictionary1, self.fastqDictionary2, self.formatExpression, self.formatExpression2, self.argsSeed)
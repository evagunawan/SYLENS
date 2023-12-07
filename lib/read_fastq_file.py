#!/usr/bin/env python3
import gzip
import logging
import re 
import sys
import tempfile as temp
import shutil

from Bio import SeqIO

from lib.alternative_dictionary_processing import process_alternative_dictionary, get_key
from lib.single_end_processing import process_single_end_sampling
from lib.interleaved_processing import process_interleaved_sampling
from lib.paired_end_processing import process_paired_end_sampling
from lib.second_file_format import determine_second_file_format

'''
Class to make a custom FastqFile object. This gives me the ability to make my own 
data blueprint. This object can then be applied throughout the main script
'''

class FastqFileData:

    #Creating init class for all the arguments
    def __init__(self, Read1Path, Read2Path, SubsampleLevel, OutputFormat, CompressOutput, Filetype, Seed, SubsamplePercentage):

        self.Read1Path = Read1Path
        self.Read2Path = Read2Path
        self.SubsampleLevel = SubsampleLevel
        self.OutputFormat = OutputFormat
        self.CompressOutput = CompressOutput 
        self.Filetype = Filetype   
        self.Seed = Seed
        self.SubsamplePercentage = SubsamplePercentage

    #Setting up function to process fastq file
    def reading_fastq_file(self):

        #Creating function to call to create dictionary
        def create_dictionary(input_file, Filetype):

            if input_file.endswith('.gz'):
            
                with gzip.open(input_file, 'rt') as infile: 
                
                    fastqDictionary = SeqIO.to_dict(SeqIO.parse(infile, Filetype)) 
    
            else:
     
                with open(input_file, 'rt') as infile:        

                    fastqDictionary = SeqIO.to_dict(SeqIO.parse(infile, Filetype))
    
            return fastqDictionary
        
        #If only one file was entered
        if self.Read2Path == None:

            logging.info('Beginning to process one input file.')

            #Made exception. First, try creating a dictionary. With single illumina/casava files, dictionary will be created but will not have read number 
            try:

                logging.debug('Beginning Try statement')
                self.fastqDictionary1 = create_dictionary(self.Read1Path, self.Filetype)
                self.fastqDictionary2 = {None:None, None:None}

            #In cases of interleaved illumina/casava, duplicate keys are created since read # is found in description, will instead create blank dictionaries 
            except ValueError:

                logging.debug('Beginning ValueError statement')
                self.fastqDictionary1 = {None:None, None:None}
                self.fastqDictionary2 = {None:None, None:None}

        #If two files are entered
        else:

            logging.info('Beginning to process two input files.')

            #Made exception. First try making a dict. If illumina/casava both dict will have same keys, but read number is found in description
            try:

                logging.debug('Beginning Try statement')               
                self.fastqDictionary1 = create_dictionary(self.Read1Path, self.Filetype)
                self.fastqDictionary2 = create_dictionary(self.Read2Path, self.Filetype)

            #In cases where incorrect files were uploaded, i.e. two casava read 2s, this creates blank dict to send to alternative dict process where problem can be identified.
            except ValueError:

                logging.debug('Beginning ValueError statement')
                self.fastqDictionary1 = {None:None, None:None}
                self.fastqDictionary2 = {None:None, None:None}

        logging.debug('Returning fastqDictionary1 and fastqDictionary2')

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

        #Grabbing first and last IDs from the input file dictionaries
        self.first_ID_1 = list(self.fastqDictionary1) [0]
        self.last_ID_1 = list(self.fastqDictionary1) [-1]

        self.first_ID_2 = list(self.fastqDictionary2) [0]
        self.last_ID_2 = list(self.fastqDictionary2) [-1]        

        logging.debug(f'Acquired first and last IDs: {self.first_ID_1}, {self.last_ID_1}, {self.first_ID_2}, {self.last_ID_2}')

        #Exception to first try finding a pattern in first_ID. This should find a pattern in any NCBI fomatted file
        try:

            logging.debug('Starting try statment in determine fastq ID formatting')

            for pattern in format_dictionary:

                #If re finds pattern in dictionary with read 1 IDs saves format and expression
                if re.search(pattern, self.first_ID_1):

                    self.formatExpression = pattern

                    self.format = format_dictionary[pattern]

                    logging.debug(f'Regular expression found fastq ID format is {self.formatExpression}')                

                    completed = True

                    self.formatExpression2 = get_key(self.format, format_dictionary_2)
                    
                    #Used to stop processing of two different file formats are found, i.e. for NCBI with rev Illum files
                    if not re.search(self.formatExpression2, self.first_ID_2):
                        
                        logging.critical('There is an issue with file formats. i.e. two files with two different file formats like NCBI paired with Illumina, two interleaved files, or two forward read files. Program terminating...')

                        sys.exit(1)

                    logging.debug('Returning self.formatExpression and self.formatExpression2')

                    return self.formatExpression, self.formatExpression2

                #If a pattern cannot be found in forward dictionary
                else:

                    #Added an error to catch read 2 file instead of read 1 or interleaved file
                    for pattern in format_dictionary_2:

                        if re.search(pattern, self.first_ID_1):

                            logging.critical('Required positional file appears to be a reverse read file. Positional file should be a forward read or interleaved file. Program terminating...')

                            sys.exit(1)  

            #If a dictionary was created for a file but a pattern is not found because the read ID was found in the description, feeds into alternative dictionary
            if completed == False:

                logging.debug('No format expression was found. Now creating a new dictionary to access file input.')

                self.fastqDictionary1, self.fastqDictionary2, self.first_ID_1, self.last_ID_1, self.first_ID_2, self.last_ID_2, self.format, self.formatExpression, self.formatExpression2 = process_alternative_dictionary(self.Read1Path, self.Read2Path, self.Filetype, format_dictionary, format_dictionary_2)

            #If a second file is uploaded and a first ID for that file is found makes sure files have same formatting
            if self.first_ID_2 != None:

                logging.debug('Determine if second file matches first file formatting')

                self.formatExpression2 = determine_second_file_format(self.Read2Path, self.Filetype, self.first_ID_2, format_dictionary, format_dictionary_2, self.format)

                logging.debug('Returning self.first_ID_1, self.first_ID_2, self.last_ID_1, self.last_ID_2, self.format, self.formatExpression, self.formatExpression2')
            
            return self.first_ID_1, self.first_ID_2, self.last_ID_1, self.last_ID_2, self.format, self.formatExpression, self.formatExpression2

        #Creates a dictionary with keys that are the descriptions. Occurs if dictionary is made from an interleaved illumina/casava file and has no matches which creates a TypeError 
        except TypeError:

            if self.format != None:

                pass
            
            else:

                logging.debug('TypeError path for dictionary production. Occurs when duplicate keys are made in dictionary from interleaved illumina/casava files.')

                self.fastqDictionary1, self.fastqDictionary2, self.first_ID_1, self.last_ID_1, self.first_ID_2, self.last_ID_2, self.format, self.formatExpression, self.formatExpression2 = process_alternative_dictionary(self.Read1Path, self.Read2Path, self.Filetype, format_dictionary, format_dictionary_2)
 
                logging.debug('Looking to see if second uploaded file matches the first file format.')

            #Checking file format to ensure both files are same file format
            if self.first_ID_2 != None:

                self.formatExpression2 = determine_second_file_format(self.Read2Path, self.Filetype, self.first_ID_1, format_dictionary, format_dictionary_2, self.format)

            logging.debug('Return self.first_ID_1, self.first_ID_2, self.last_ID_1, self.last_ID_2, self.format, self.formatExpression, self.formatExpression2')

            return self.first_ID_1, self.first_ID_2, self.last_ID_1, self.last_ID_2, self.format, self.formatExpression, self.formatExpression2, self.fastqDictionary1, self.fastqDictionary2

    #Figures out if input file is single, paired, or interleaved
    def determine_paired_single_interleaved(self):

        logging.debug('Setting up second file expression')
        
        if self.Read2Path == None:

            self.single_file = True
            self.paired_file = False

        if self.Read2Path != None:

            self.single_file = False
            self.paired_file = True

        #If only one file was entered, figured out if it is interleaved or single end based on RE 
        if self.single_file == True:
            
            if self.formatExpression2 != None:

                if re.search(self.formatExpression, self.first_ID_1) and re.search(self.formatExpression2, self.last_ID_1):

                    self.determined_filetype = 'Interleaved'

                    logging.info('File is an interleaved file.')

            if re.search(self.formatExpression, self.first_ID_1) and re.search(self.formatExpression, self.last_ID_1):

                self.determined_filetype = 'Single-end'

                logging.info('File is a single-end file.')

        #If two files are entered, will initially process mistakes, i.e. 2 rev files. If passes mistakes, determines if paired end
        if self.paired_file == True:

            if re.search(self.formatExpression, self.first_ID_1) and re.search(self.formatExpression, self.last_ID_1):

                logging.debug('Forward file is correct.')

                if re.search(self.formatExpression2, self.first_ID_2) and re.search(self.formatExpression2, self.last_ID_2):

                    logging.debug('Reverse file is correct.')
    
                    self.determined_filetype = 'Paired-end'

        logging.debug('Returned self.determined_filetype')

        return self.determined_filetype

    #Writing output of determined filetype
    def processing_filetype(self):

        #Calculating percentage, if desired
        def percentage(subsample, dictionary):

            logging.debug('Determined number to sample if percentage == True')

            if self.SubsamplePercentage == True:

                length = len(dictionary)

                self.SubsampleLevel = int(length) * int(subsample) / int(100)
 
                self.SubsampleLevel = int(self.SubsampleLevel)

                logging.debug(f'Returning self.SubsampleLevel as {self.SubsampleLevel}')
                
                return self.SubsampleLevel
            
            else:

                return subsample
        
        self.SubsampleLevel = percentage(self.SubsampleLevel, self.fastqDictionary1)
        
        if self.determined_filetype == 'Single-end':

            process_single_end_sampling(self.Read1Path, self.Read2Path, self.SubsampleLevel, self.OutputFormat, self.CompressOutput, self.fastqDictionary1, self.Seed)

        elif self.determined_filetype == 'Interleaved':

            process_interleaved_sampling(self.Read1Path, self.SubsampleLevel, self.OutputFormat, self.CompressOutput, self.fastqDictionary1, self.Seed, self.formatExpression, self.formatExpression2)

        elif self.determined_filetype == 'Paired-end':

            process_paired_end_sampling(self.Read1Path, self.Read2Path, self.SubsampleLevel, self.OutputFormat, self.CompressOutput, self.fastqDictionary1, self.fastqDictionary2, self.formatExpression, self.formatExpression2, self.Seed)
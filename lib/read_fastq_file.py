import gzip
import logging
import re 
import sys
import uuid
import shutil
import os

from Bio import SeqIO

from lib.alternative_dictionary_processing import process_alternative_dictionary, get_key
from lib.single_end_processing import process_single_end_sampling
from lib.interleaved_processing import process_interleaved_sampling
from lib.paired_end_processing import process_paired_end_sampling

from lib.fastq_format_check import determine_fastq_ID_formatting

'''
Class to make a custom FastqFile object. This gives me the ability to make my own 
data blueprint. This object can then be applied throughout the main script
'''

class FastqFileData:

    Read1Temp = None
    Read2Temp = None

    Read1Index = None
    Read2Index = None

    Read1Format = None
    Read2Format = None

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

        ####
        #### If files are compressed create temp files and decompress gzip into temp files
        ####

        if self.Read1Path and self.Read1Path.endswith('.gz'):
            
            logging.debug(f"Decompressing Read1 {self.Read1Path}")
            
            self.Read1Temp = uuid.uuid4().hex
            with open(self.Read1Temp, "wb") as tmp:
                shutil.copyfileobj(gzip.open(self.Read1Path), tmp)
        
        if self.Read2Path and self.Read2Path.endswith('.gz'):

            logging.debug(f"Decompressing Read2 {self.Read2Path}")

            self.Read2Temp = uuid.uuid4().hex
            with open(self.Read2Temp, "wb") as tmp:
                shutil.copyfileobj(gzip.open(self.Read2Path), tmp)




        ####
        #### Creating fastq index
        ####

        logging.info(f"Creating Read Indexes")
        
        # Read 1
        logging.debug('Beginning to index read 1.')
        if self.Read1Temp:
            logging.debug('Read 1 was compressed, indexing temp file.')
            self.Read1Index = SeqIO.index(self.Read1Temp, Filetype)
        
        elif self.Read1Path:
            logging.debug('Indexing read 1.')
            self.Read1Index = SeqIO.index(self.Read1Path, Filetype)
        
        else:
            logging.error('Missing read 1.')
            sys.exit(1)

        # Read 2
        logging.debug('Beginning to index read 2 (if we have it).')
        if self.Read2Temp:
            logging.debug('Read 2 was compressed, indexing temp file.')
            self.Read2Index = SeqIO.index(self.Read2Temp, Filetype)
        
        elif self.Read2Path:
            logging.debug('Indexing read 2.')
            self.Read2Index = SeqIO.index(self.Read2Path, Filetype)
        
        else:
            logging.debug('Read 2 not found, continuing.')



        
        ####
        #### Get Fastq formats
        ####

        if self.Read1Index:
            logging.debug('Getting read 1 fastq format.')
            self.Read1Format = determine_fastq_ID_formatting(self.Read1Index)

        if self.Read2Index:
            logging.debug('Getting read 2 fastq format.')
            self.Read2Format = determine_fastq_ID_formatting(self.Read2Index)
    
    
    
    
    ####
    #### Clean up function to delete temp files after completion
    ####
    def cleanUP(self):
        logging.info("Cleaning Temp Files.")

        if self.Read1Temp:
            os.remove(self.Read1Temp)

        if self.Read2Temp:
            os.remove(self.Read2Temp)


    

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
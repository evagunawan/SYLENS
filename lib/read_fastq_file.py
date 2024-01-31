import gzip
import logging
import sys
import uuid
import shutil
import os

from lib.fastq_indexing import index

from lib.fastq_format_check import determine_fastq_ID_formatting, get_interleaved_ids

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

    Interleaved_Read1_IDs = None
    Interleaved_Read2_IDs = None

    R1_Total_Reads = None
    R2_Total_Reads = None

    ####
    #### Creating init class for all the arguments
    ####
    def __init__(self, Read1Path, Read2Path, Filetype):

        self.Read1Path = Read1Path
        self.Read2Path = Read2Path
        self.Filetype = Filetype   


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
            self.Read1Index = index(self.Read1Temp,Filetype)
        
        elif self.Read1Path:
            logging.debug('Indexing read 1.')
            self.Read1Index = index(self.Read1Path,Filetype)
        
        else:
            logging.error('Missing read 1.')
            sys.exit(1)

        # Read 2
        logging.debug('Beginning to index read 2 (if we have it).')
        if self.Read2Temp:
            logging.debug('Read 2 was compressed, indexing temp file.')
            self.Read2Index = index(self.Read2Temp,Filetype)
        
        elif self.Read2Path:
            logging.debug('Indexing read 2.')
            self.Read2Index = index(self.Read2Path,Filetype)
        
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
        #### Get list of IDs if interleaved
        ####
        if self.Read1Path and not self.Read2Path:
            self.Interleaved_Read1_IDs, self.Interleaved_Read2_IDs = get_interleaved_ids(self.Read1Index)

        ####
        #### Determine Number of Reads
        ####
        logging.debug('Calculating read length.')
        if self.Interleaved_Read1_IDs:
            self.R1_Total_Reads = len(self.Interleaved_Read1_IDs)
        
        else:
            self.R1_Total_Reads = len(self.Read1Index.keys())

        if self.Interleaved_Read2_IDs:
            self.R2_Total_Reads = len(self.Interleaved_Read2_IDs)
        
        if self.Read2Index:
            self.R2_Total_Reads = len(self.Read2Index.keys())


    ####
    #### Clean up function to delete temp files after completion
    ####
    def cleanUP(self):
        logging.info("Cleaning Temp Files.")

        if self.Read1Temp:
            os.remove(self.Read1Temp)

        if self.Read2Temp:
            os.remove(self.Read2Temp)
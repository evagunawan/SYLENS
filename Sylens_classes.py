#!/usr/bin/env python3
import logging
import gzip
import argparse
import re
import sys

from Bio import SeqIO
# from DualInputs import dual_dictionary
# from SingularInput import single_dictionary 

def parse_FormatPair():
    
    #Taking the 3 known formats and creating regular expressions with them. Illumina and Casava have same format, right now.
        global seqIDformat

        seqIDformat = {

        'IlluminaAndCasava' : "r'^@(.+) (\d)'",
        'NCBI' : "r'^@(\w+).(\d).([12]) '"
        }

        # print([(information, seqIDformat[information]) for information in seqIDformat])


#Think of class as a noun and definition as a verb
class FastqFile:
    
    def __init__(self, read1):
          self.read1 = read1
    def printing(self):
          return self.read1
sys.exit()
    # def determine_SeqID_format(self):

    #         for ID in self:
                
    #             count = 0
    #             if re.search(parse_FormatPair[count], ID):
                    
    #                 format_of_ID = parse_FormatPair()
    #                 print(f'Format is: {format_of_ID}')

    # def gathering_file_info(self, path, filetype):
        
    #     #This is getting the path and filetype
    #     self.path = path
    #     self.filetype = filetype

    #     #This checks if compressed
    #     if self.path.endswith('.gz'):
    #         self.path = True
    #         logging.debug('Ended with .gz')
        
    #     else:
    #         self.path = False
    #         logging.debug('Not .gz')
        
    #     #This will read the data from the file into memory instead of creating a new file
    #     if self.path == True:
    #         with gzip.open(self.path, 'rt') as infile:
    #             self.readsDictionary = SeqIO.to_dict(SeqIO.parse(infile, self.filetype))
    #             logging.debug('created dictionary for .gz')
   
    #     else:
    #         with open(self.path, 'rt') as infile:
    #             self.readsDictionary = SeqIO.to_dict(SeqIO.parse(infile, self.filetype))
    #             logging.debug('created dictionary for not .gz')
        
    #     #Looks at the first and last reads of the file
    #     #scan through readIDs to look for mix of read1 and read2
    #     global self_readIDs
    #     self_readIDs = self.readsDictionary.keys() 

    # determine_SeqID_format(self_readIDs)
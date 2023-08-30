#!/usr/bin/env python3

'''
TODO Think about how to do the other single end and interleaved etc. 
TODO Far future: Maybe add sam/bam file
TODO Not care about switching r1 and r2
IN_PROCESS Learn about exceptions
DONE? Logging library https://docs.python.org/3/library/logging.html#logging-levels
'''

#SYLENS: Sampling Yielding LEss Noticeable Samples

###Importing libraries### 
import argparse
import re
import random
import sys 
import time
import gzip
import shutil
import os
import logging

from Bio import SeqIO

# Creating class that displays help if no information is entered in parser
class DownsamplerParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help()
        sys.stderr.write(f'\nERROR: {message}\n')
        sys.exit(1)

# Downsampler program description 
parser = DownsamplerParser(prog = 'Subsampler for FASTQ file(s)',
    description = 'Enter in FASTQ file(s) and down sample based on a user supplied integer. If no user input is added the entire file is sampled and output as chosen filetype.',
    epilog = "Additional information and a README.md can be found at https://github.com/evagunawan/SYLENS"
    )
parser.add_argument('Read1', 
    help = 'Add R1/Single-end/Interleaved FASTQ file here. I.E. filename.fastq'
    )
parser.add_argument('Read2', 
    nargs = '?',
    help = 'Add additional FASTQ file here if paired-end. I.E. R2.fastq'
    )
parser.add_argument('-s','--subsample', 
    type = int,
    help = 'Enter an integer which will be the total number of down sampling of FASTQ files occuring. Leave blank if no subsampling is desired and file conversion is needed. I.E. --subsample 10000.'
    )
parser.add_argument('-f', '--filetype',
    choices = ['fastq', 'fastq-solexa'],
    default = 'fastq',
    help = "Add what type of fastq file is being input. fastq/fastq-sanger-> uses ASCII offset of 33 whereas fastq-solexa-> uses ASCII of 64. I.E -f fastq-solexa"
    )
parser.add_argument('--seed', 
    type = int,
    default = int(time.time()),
    help = 'Enter the seed number if you would like to reproduce previous results. I.E. --seed 1691696502'
    )
parser.add_argument('-o', '--output',
    choices = ['fastq', 'fastq-solexa'],
    default = 'fastq',
    help = "Add what type of fastq file is desired at output. I.E --output fastq-solexa"
    )
parser.add_argument('-c', '--compress',
    default = False,
    action = 'store_true',
    help = "Compress fastq file into fastq.gz file on output. Common typing mistakes are included. I.E. -c no"
    )

#Runs the parser and allows you to call arguments downstream
args = parser.parse_args()

#Format for logging debug-critical information
logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

'''
Set seed value to time since epoch that way each time program is run, 
seed value is different. If user enters seed value, program will 
reproduce the same results 
'''
random.seed(args.seed)

print(args.Read1, args.Read2, args.subsample, args.filetype, args.seed, args.output, args.compress)





# if args.Read1 == True:
#     FastqFile(args.Read1)


# sys.exit()

# # if args.Read2:
# #     read2 = ReadingFile(args.Read2, args.filetype)
# #     print('True')

# # if not args.Read2 and not args.interleaved:
# #     read1.singleEnd = True
# #     print('TRue')


# # if args.Read2:
# #     read2 = ReadingFile(args.Read2, args.filetype, args.subsample)

# #Sends read 1 through reading file argument
# read1 = ReadingFile(args.Read1, args.filetype)

# if args.Read2 != None:
#     read2 = ReadingFile(args.Read2, args.filetype)

# if args.Read2 == None:
#     read1.analysis_of_single_file()

# ### Defining a few empty lists ###
# subsample_max_IDs = []
# subsample_max_IDs_2 = []
# interleaved_R1_IDs = []
# interleaved_R2_IDs = []
# r1_IDS = []
# r2_IDS = []
# r1_split_ids = []
# r2_split_ids = []
# paired_end_1 = []
# paired_end_2 = []

# ############################################### Creating a bunch of custom functions for single/paired/interleaved interpretation, zipping, etc. ###############################################

# # function to unzip file 
# def extract_gzip_file(input_file): 
#     logging.info('File is compressed. Opening...')
#     gzip.open(input_file)


# '''
# If subsampling is occurring, the generated dictionary is randomly sampled. The randomly sampled IDs are saved in a list. Using that list, the IDs are mapped back to
# their respective entires in the called dictionary. SeqIO writes the output file with the filetype indicated by the user. IF no downsampling occurs, the fastq is turned
# into a sequence_record object and written to the output file.
# '''
# #function to downsample single files
# def single_end_file(dictionary_name, subsample_max):
#     if args.subsample != None:
#         logging.info("Downsampling occuring.")
#         Random_IDs = random.sample(list(dictionary_name), subsample_max)
#         single_file_output = [dictionary_name[info] for info in Random_IDs]
#         logging.info('Writing to file...')    
        
#         if args.Read1[0].endswith('.gz'):
#             with gzip.open(f"down_sampled_single_end_{args.Read1[0]}", 'w') as IDs_Seq_File:
#                 SeqIO.write(single_file_output, IDs_Seq_File, args.output)
#         else:
#             if args.compress == 'no' or args.compress == 'NO' or args.compress == 'No' or args.compress == 'nO':
#                 with open(f"down_sampled_{args.Read1[0]}", 'w') as IDs_Seq_File:
#                     SeqIO.write(single_file_output, IDs_Seq_File, args.output)
#             if args.compress == 'yes' or args.compress == 'YES' or args.compress == 'Yes' or args.compress == 'yES': 
#                 with gzip.open(f"down_sampled_single_end_{args.Read1[0]}", 'w') as IDs_Seq_File:
#                     SeqIO.write(single_file_output, IDs_Seq_File, args.output)
        
#         logging.info(f"Done writing to: down_sampled_single_end_{args.Read1[0]}")
#         logging.debug(f'SUCCESS: finished down_sampled_single_end_{args.Read1[0]}')   
    
#     if args.subsample == None:
#         logging.info("No downsampling occurring.")
#         All_IDs = list(dictionary_name)
#         All_output = [dictionary_name[info] for info in All_IDs]
#         logging.info('Writing to file...')
        
#         if args.Read1[0].endswith('.gz'):
#             with gzip.open(f"non_down_sampled_single_end_{args.Read1[0]}", 'wb') as IDs_Seq_File:
#                 SeqIO.write(All_output, IDs_Seq_File, args.output)
        
#         else:
#             if args.compress == 'no' or args.compress == 'NO' or args.compress == 'No' or args.compress == 'nO':
#                 with open(f"non_down_sampled_single_end_{args.Read1[0]}", 'w') as IDs_Seq_File:
#                     SeqIO.write(All_output, IDs_Seq_File, args.output)
#             if args.compress == 'yes' or args.compress == 'YES' or args.compress == 'Yes' or args.compress == 'yES': 
#                 with gzip.open(f"non_down_sampled_single_end_{args.Read1[0]}", 'wb') as IDs_Seq_File:
#                     SeqIO.write(All_output, IDs_Seq_File, args.output)
        
#         logging.info(f"Done writing to: non_down_sampled_single_end_{args.Read1[0]}")
#         logging.debug(f'SUCCESS: finished non_downsample_single_end_{args.Read1[0]}')


# '''
# Regular expressions store R1/R2 IDs into separate lists. Lists zipped to randomly sample tuples. Rezip randomly sampled 
# IDs into separate lists. Search their respective dictionary. Pull input that matches R1/R2's randomly sampled IDS. Store as
# output. Separate interleaved fastq with R1 first and R2 second using counts for R1 and R2 to index the output list. Write indexed value
# to output file. Count = 1 so R1 record is written first. Added 1 total length since the count < int is exclusive.
# '''
# #function to downsample interleaved files
# def interleaved_file(dictionary_name, subsample_max):
#     if subsample_max != None:
#         logging.info("File is interleaved. Processing interleaved file for downsampling.")
#         for element in list(dictionary_name):
#             if re.search(r"(^\w+.)(\w+)(.)(2$)", element):
#                 interleaved_R2_IDs.append(element)
#             if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
#                 interleaved_R1_IDs.append(element)
#         zipped_ids_list = list(zip(interleaved_R1_IDs, interleaved_R2_IDs))
#         random_ids = random.sample(zipped_ids_list, subsample_max)
#         r1_ids, r2_ids = zip(*random_ids)
#         output_1 = [dictionary_name[info] for info in r1_ids]
#         output_2 = [dictionary_name[info] for info in r2_ids]  
#         logging.info('Writing to file...')
#         count = 1
#         R1_count = 0
#         R2_count = 0
        
#         if args.Read1[0].endswith('.gz'):
#             with gzip.open(f"down_sampled_interleaved_{args.Read1[0]}", 'w') as IDs_Seq_File:
#                 while count < (1 + 2*len(output_1)):
#                     if count % 2 != 0 :
#                         count += 1
#                         SeqIO.write(output_1[R1_count], IDs_Seq_File, args.output)
#                         R1_count += 1
#                     else:
#                         if count % 2 == 0:
#                             count+= 1
#                             SeqIO.write(output_2[R2_count], IDs_Seq_File, args.output)
#                             R2_count += 1
        
#         else:
#             if args.compress == 'no' or args.compress == 'NO' or args.compress == 'No' or args.compress == 'nO':
#                 with open(f"down_sampled_interleaved_{args.Read1[0]}", 'w') as IDs_Seq_File:
#                     while count < (1 + 2*len(output_1)):
#                         if count % 2 != 0 :
#                             count += 1
#                             SeqIO.write(output_1[R1_count], IDs_Seq_File, args.output)
#                             R1_count += 1
#                         else:
#                             if count % 2 == 0:
#                                 count+= 1
#                                 SeqIO.write(output_2[R2_count], IDs_Seq_File, args.output)
#                                 R2_count += 1
            
#             if args.compress == 'yes' or args.compress == 'YES' or args.compress == 'Yes' or args.compress == 'yES': 
#                 with gzip.open(f"down_sampled_interleaved_{args.Read1[0]}", 'w') as IDs_Seq_File:
#                     while count < (1 + 2*len(output_1)):
#                         if count % 2 != 0 :
#                             count += 1
#                             SeqIO.write(output_1[R1_count], IDs_Seq_File, args.output)
#                             R1_count += 1
#                         else:
#                             if count % 2 == 0:
#                                 count+= 1
#                                 SeqIO.write(output_2[R2_count], IDs_Seq_File, args.output)
#                                 R2_count += 1
        
#         logging.info(f'Done writing to: down_sampled_interleaved_{args.Read1[0]}')
#         logging.debug(f'SUCCESS: finished down_sample_interleaved_{args.Read1[0]}')
    
    
#     if subsample_max == None:
#         logging.info("File is interleaved. No downsampling occurring. Processing interleaved file...")
#         for element in list(dictionary_name):
#             if re.search(r"(^\w+.)(\w+)(.)(2$)", element):
#                 interleaved_R2_IDs.append(element)
#             if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
#                 interleaved_R1_IDs.append(element)
#         output_1 = [dictionary_name[info] for info in interleaved_R1_IDs]
#         output_2 = [dictionary_name[info] for info in interleaved_R2_IDs]  
#         logging.info('Writing to file...')
#         count = 1
#         R1_count = 0
#         R2_count = 0
        
#         if args.Read1[0].endswith('.gz'):
#             with gzip.open(f"non_down_sampled_interleaved_{args.Read1[0]}", 'w') as IDs_Seq_File:
#                 while count < (1 + 2*len(output_1)):
#                     if count % 2 != 0 :
#                         count += 1
#                         SeqIO.write(output_1[R1_count], IDs_Seq_File, args.output)
#                         R1_count += 1
#                     else:
#                         if count % 2 == 0:
#                             count+= 1
#                             SeqIO.write(output_2[R2_count], IDs_Seq_File, args.output)
#                             R2_count += 1
        
#         else:
#             if args.compress == 'no' or args.compress == 'NO' or args.compress == 'No' or args.compress == 'nO':
#                 with open(f"non_down_sampled_interleaved_{args.Read1[0]}", 'w') as IDs_Seq_File:
#                     while count < (1 + 2*len(output_1)):
#                         if count % 2 != 0 :
#                             count += 1
#                             SeqIO.write(output_1[R1_count], IDs_Seq_File, args.output)
#                             R1_count += 1
#                         else:
#                             if count % 2 == 0:
#                                 count+= 1
#                                 SeqIO.write(output_2[R2_count], IDs_Seq_File, args.output)
#                                 R2_count += 1
            
#             if args.compress == 'yes' or args.compress == 'YES' or args.compress == 'Yes' or args.compress == 'yES': 
#                 with open(f"non_down_sampled_interleaved_{args.Read1[0]}", 'w') as IDs_Seq_File:
#                     while count < (1 + 2*len(output_1)):
#                         if count % 2 != 0 :
#                             count += 1
#                             SeqIO.write(output_1[R1_count], IDs_Seq_File, args.output)
#                             R1_count += 1
#                         else:
#                             if count % 2 == 0:
#                                 count+= 1
#                                 SeqIO.write(output_2[R2_count], IDs_Seq_File, args.output)
#                                 R2_count += 1
        
#         logging.info("Done writing to: non_down_sampled_interleaved")
#         logging.debug('SUCCESS: finished no_downsample_interleaved_file')


# '''
# Creates two lists and searches for every .1 and .2 element in the given dictionary. Then creates two separate files. A list of all of the IDs
# are zipped together and randomly sampled as tuples. Rezipped apart into two lists. Saves the output and writes it to the filetype indicated by
# the user in two distinct separate files. 
# '''
# #function to downsample paired end files

# '''
# #function to downsample switched paired end file
# # def switched_paired_end_files(dictionary_name_R2, dictionary_name_R1, subsample_max):
# #     if subsample_max != None:
# #         logging.info("Reads 1 and 2 supplied. Downsampling occurring. Processing Read 1 and Read 2 files...")
# #         for element in list(dictionary_name_R1):
# #             if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
# #                 paired_end_1.append(element)
# #         for element in list(dictionary_name_R2): 
# #             if re.search(r"(^\w+.)(\w+)(.)(2$)", element):
# #                 paired_end_2.append(element)
# #         zipped_ids_list = list(zip(paired_end_1, paired_end_2))
# #         random_ids = random.sample(zipped_ids_list, subsample_max)
# #         paired_end_1_ids, paired_end_2_ids = zip(*random_ids)
# #         output_2 = [dictionary_name_R2[info] for info in paired_end_2_ids]
# #         output_1 = [dictionary_name_R1[info] for info in paired_end_1_ids]  
# #         logging.info('Writing to file...')
# #         if args.Read1[0].endswith('.gz'):
# #             with gzip.open(f"down_sampled_{args.Read1[0]}", 'w') as IDs_Seq_File_2:
# #                 SeqIO.write(output_2, IDs_Seq_File_2, args.output)
        

# #         if args.Read2.endswith('.gz'):
# #             with gzip.open(f'down_sample_{args.Read2}', 'w') as IDs_Seq_File_1:
# #                 SeqIO.write(output_1, IDs_Seq_File_1, args.output)
        

# #         if args.Read1[0].endswith('.fastq'):
# #             if args.compress == 'no' or args.compress == 'NO' or args.compress == 'No' or args.compress == 'nO':
# #                 with open(f"down_sampled_{args.Read1[0]}", 'w') as IDs_Seq_File_2:
# #                     SeqIO.write(output_2, IDs_Seq_File_2, args.output)
# #             if args.compress == 'yes' or args.compress == 'YES' or args.compress == 'Yes' or args.compress == 'yES': 
# #                 with gzip.open(f"down_sampled_{args.Read1[0]}", 'w') as IDs_Seq_File_2:
# #                     SeqIO.write(output_2, gzip(IDs_Seq_File_2), args.output)
        

# #         if args.Read2.endswith('.fastq'):
# #             if args.compress == 'no' or args.compress == 'NO' or args.compress == 'No' or args.compress == 'nO':
# #                 with open(f"down_sampled_{args.Read2}", 'w') as IDs_Seq_File_1:
# #                     SeqIO.write(output_1, IDs_Seq_File_1, args.output)
# #             if args.compress == 'yes' or args.compress == 'YES' or args.compress == 'Yes' or args.compress == 'yES': 
# #                 with gzip.open(f"down_sampled_{args.Read2}", 'w') as IDs_Seq_File_1:
# #                     SeqIO.write(output_1, gzip(IDs_Seq_File_1), args.output)
        
# #         logging.info(f"Done writing to: down_sampled_{args.Read1[0]} and down_sampled_{args.Read2}") 
# #         logging.debug('SUCCESS: finished switched paired end')
    

# #     if subsample_max == None:
# #         logging.info("Reads 2 and 1 supplied. No downsampling occurring. Processing Read 1 and Read 2 files...")
# #         for element in list(dictionary_name_R1):
# #             if re.search(r"(^\w+.)(\w+)(.)(2$)", element):
# #                 paired_end_2.append(element)
# #         for element in list(dictionary_name_R2): 
# #             if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
# #                 paired_end_1.append(element)
# #         output_2 = [dictionary_name_R1[info] for info in paired_end_2]
# #         output_1 = [dictionary_name_R2[info] for info in paired_end_1]  
# #         logging.info('Writing to file...')
        
        
# #         if args.Read2.endswith('.fastq'):
# #             if args.compress == 'no' or args.compress == 'NO' or args.compress == 'No' or args.compress == 'nO':
# #                 with open(f"non_down_sampled__{args.Read2}", 'w') as IDs_Seq_File_1:
# #                     SeqIO.write(output_1, IDs_Seq_File_1, args.output)
# #             if args.compress == 'yes' or args.compress == 'YES' or args.compress == 'Yes' or args.compress == 'yES': 
# #                 with gzip.open(f"non_down_sampled_{args.Read2}", 'w') as IDs_Seq_File_1:
# #                     SeqIO.write(output_1, IDs_Seq_File_1, args.output)
        
        
# #         if args.Read1[0].endswith('.fastq'):
# #             if args.compress == 'no' or args.compress == 'NO' or args.compress == 'No' or args.compress == 'nO':
# #                 with open(f"non_down_sampled_{args.Read1[0]}", 'w') as IDs_Seq_File_2:
# #                     SeqIO.write(output_2, IDs_Seq_File_2, args.output)
# #             if args.compress == 'yes' or args.compress == 'YES' or args.compress == 'Yes' or args.compress == 'yES': 
# #                 with gzip.open(f"non_down_sampled_{args.Read1[0]}", 'w') as IDs_Seq_File_2:
# #                     SeqIO.write(output_2, IDs_Seq_File_2, args.output)


# #         if args.Read1[0].endswith('.gz'):
# #             with gzip.open(f"non_down_sampled_{args.Read1[0]}", 'w') as IDs_Seq_File_2:
# #                 SeqIO.write(output_2, IDs_Seq_File_2, args.output)
        
        
# #         if args.Read2.endswith('.gz'):
# #             with gzip.open(f'non_down_sampled_{args.Read2}', 'w') as IDs_Seq_File_1:
# #                 SeqIO.write(output_1, IDs_Seq_File_1, args.output)
        
        
# #         logging.info(f"Done writing to: non_down_sampled_{args.Read1[0]} and non_down_sampled_{args.Read2}")
# #         logging.debug(f'SUCCESS: finished switched paired-end')


# #function to zip files
# def gzip_output_file(original_file, new_output_file):
#     with open(original_file, 'rb') as input_file:
#         with gzip.open(new_output_file, 'wb') as output_file:
#             shutil.copyfileobj(input_file, output_file)
#             logging.info(f'Renaming file: {new_output_file}')
#             if os.path.exists(original_file):
#                 os.remove(original_file)
#             logging.debug(f'SUCCESS: finished gzip {new_output_file}')
# '''




# ''' 
# These statements check to make sure user input has been pulled and allows to user to know their input was taken
# They are not necessary, just somewhat helpful 
# ''' 
# if args.Read2 == None:
#     logging.info(f'The FASTQ filename is: {args.Read1[0]}')
# else:
#     logging.info(f'The FASTQ filenames are: {args.Read1[0]} and  {args.Read2}')

# logging.info(f'The input file type is: {args.filetype}')

# logging.info(f"The run's seed number is: {args.seed}")

# logging.info(f'The output file type is: {args.output}')

# logging.info(f'Compressing on output: {args.compress}')


# '''
# Ensures that the subsample integer is pulled in as a usable int not str variable. Also allows for no subsampling and just final conversion.
# Makes sure subsample size is not larger than total file size. If so terminates program 
# '''
# if args.subsample != None:
#     subsample_max = int(args.subsample)
#     logging.info(f'The amount to subsample is: {args.subsample}')
#     logging.info('Checking to see if subsample size is larger than file')
#     if args.Read1[0].endswith('.gz'):
#         with gzip.open(args.Read1[0], 'rt') as input_file:
#             Dict_for_subsample_max = SeqIO.to_dict(SeqIO.parse(input_file, args.filetype))
#         logging.debug('SUCCESS: made dictionary for subsampled zipped file')
#     else:
#         Dict_for_subsample_max = SeqIO.to_dict(SeqIO.parse(args.Read1[0], args.filetype))
#         logging.debug('SUCCESS: made dictionary for subsampled unzipped file')
#     subsample_max_IDs = []
#     subsample_max_IDs_2 = []
#     for element in list(Dict_for_subsample_max):
#         if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
#             subsample_max_IDs.append(element)
#         if re.search(r"(^\w+.)(\w+)(.)(2$)", element):
#             subsample_max_IDs_2.append(element)
#     logging.debug('SUCCESS: made .1 list for subsample max')
#     logging.debug('SUCCESS: made .2 list for subsample max')
#     total_length = len(subsample_max_IDs)
#     if total_length == 0:
#         total_length = len(subsample_max_IDs_2) 
#         logging.warning('File 1 and 2 may be switched. Recalculating length of file to ensure subsample size is not larger than file size...')
#     if total_length < subsample_max:
#         logging.debug('FAILURE: Total length is less than subsample max')
#         logging.critical(f'Subsample size of {subsample_max} is larger than total number of samples in file. Please provide a subsample size less than {total_length}. Terminating...')
#         sys.exit(1)


# if args.subsample == None:
#     if args.Read1[0].endswith('.gz'):
#         with gzip.open(args.Read1[0], 'rt') as input_file:
#             Dict_for_subsample_max = SeqIO.to_dict(SeqIO.parse(input_file, args.filetype))
#         logging.debug('SUCCESS: made dictionary for not subsampled zipped file')
#     else:
#         Dict_for_subsample_max = SeqIO.to_dict(SeqIO.parse(args.Read1[0], args.filetype))
#         logging.debug('SUCCESS: made dictionary for not subsampled file')
#     for element in list(Dict_for_subsample_max):
#         if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
#             subsample_max_IDs.append(element)
#         else:
#             if re.search(r"(^\w+.)(\w+)(.)(2$)", element):    
#                 subsample_max_IDs_2.append(element)
#         logging.debug('SUCCESS: made .1 list for subsampled')
#         logging.debug('SUCCESS: made .2 list for subsampled')
#     subsample_max = len(subsample_max_IDs)
#     if subsample_max == 0:
#         subsample_max = len(subsample_max_IDs_2)
#     logging.info(f'Total amount of reads in the file: {subsample_max}')
    

# ############################################### Only Read 1 supplied ###############################################


# # This sees that no read 2 was entered. Then converts given fastq file into dictionary. Grabs the first and last records.
# if args.Read2 == None:
#     logging.info("No read 2 entered, checking if single-end or interleaved file...")
#     if args.Read1[0].endswith('.gz'):
#         R1_dict = SeqIO.to_dict(SeqIO.parse(extract_gzip_file(args.Read1[0]), args.filetype))
#         logging.debug('SUCCESS: made dictionary for subsampled zipped single end file')
#     else:
#         R1_dict =  SeqIO.to_dict(SeqIO.parse(args.Read1[0], args.filetype))
#         logging.debug('SUCCESS: made dictionary for subsampled unzipped single end file')
#     last_record = list(R1_dict) [-1]
#     first_record = list(R1_dict) [0]


# # If regular expression finds that both the first and last record end with '.2', it will terminate the program
#     if re.search(r"(^\w+.)(\w+)(.)(2$)", first_record):
#         logging.debug('FAILURE: Found .2 read in first record indicating .2 file was provided instead of .1 file.')
#         logging.critical("This appears to be a read 2 file instead of a read 1 file. Program terminating...")
#         sys.exit(1)    


# # Regular expression looks if an interleaved file was supplied. If first entry ends with ".1" and the last record ends with a ".2", interleaved file. If a subsample size was added, the
# # downsample_interleaved_file function is run. If a subsample size is not indicated, the no_downsample_interleaved_file is run. 
#     elif re.search(r"(^\w+.)(\w+)(.)(2$)", last_record):
#         if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record):
#             interleaved_file(R1_dict, subsample_max)
            
 
# # Uses a regular expression to look at the first and last entry of the fastq. If first and last entry end with ".1", single end file. Moves into whether or not downsampling is occurring
#     elif re.search(r"(^\w+.)(\w+)(.)(1$)", last_record):
#         if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record):
#             logging.info("File is single-end.")
#             single_end_file(R1_dict, subsample_max)


# ############################################### Both Read 1 and 2 supplied ###############################################


# # Two files supplied. Looks if a read or both are .gz, decompresses them, and assigns each file to a dictionary. Then grabs the first and last IDs of both files. 
# if args.Read2 != None:
#     if args.Read1[0].endswith('.gz'):
#         with gzip.open(args.Read1[0, 'rt']) as input_1:
#             R1_dict =  SeqIO.to_dict(SeqIO.parse(input_1, args.filetype))
#         logging.debug('SUCCESS: made dictionary for subsampled zipped R1 file')
#     else:
#         R1_dict =  SeqIO.to_dict(SeqIO.parse(args.Read1[0], args.filetype)) 
#         logging.debug('SUCCESS: made dictionary for subsampled unzipped R1 file')
#     if args.Read2.endswith('.gz'):  
#         with gzip.open(args.Read2, 'rt') as input_2:
#             R2_dict =  SeqIO.to_dict(SeqIO.parse(input_2, args.filetype))
#         logging.debug('SUCCESS: made dictionary for subsampled zipped R2 file')
#     else:
#         R2_dict =  SeqIO.to_dict(SeqIO.parse(args.Read2, args.filetype))
#         logging.debug('SUCCESS: made dictionary for subsampled unzipped R1 file')
#     last_record_1 = list(R1_dict) [-1]
#     first_record_1 = list(R1_dict) [0]
#     last_record_2 = list(R2_dict) [-1]
#     first_record_2 = list(R2_dict) [0]


# #Errors if wrong files are supplied and deletes files that were made preparing for decompression, if applicable.
#     if re.search(r"(^\w+.)(\w+)(.)(2$)", last_record_2):
#         if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record_2):
#             if re.search(r"(^\w+.)(\w+)(.)(2$)", last_record_1):
#                 if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record_1):
#                     logging.critical('Both supplied files appear to be interleaved. Program terminating...')
#                     sys.exit(1)
#     if re.search(r"(^\w+.)(\w+)(.)(1$)", last_record_1):
#         if re.search(r"(^\w+.)(\w+)(.)(1$)", last_record_2):
#             logging.critical("Both supplied files appear to be Read 1 files. Program terminating...")
#             sys.exit(1)
#     if re.search(r"(^\w+.)(\w+)(.)(2$)", first_record_1):
#         if re.search(r"(^\w+.)(\w+)(.)(2$)", first_record_2):
#             logging.critical("Both supplied files appear to be Read 2 files. Program terminating...")
#             sys.exit(1)
#     if re.search(r"(^\w+.)(\w+)(.)(2$)", last_record_1):
#         if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record_1):
#             logging.critical("First input file is interleaved. Program terminating...")
#             sys.exit(1)
#     if re.search(r"(^\w+.)(\w+)(.)(2$)", last_record_2):
#         if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record_2):
#             logging.critical("Second input file is interleaved. Program terminating...")
#             sys.exit(1)
    
    

# # Switched input files R2 first instead of R1
#     if re.search(r"(^\w+.)(\w+)(.)(1$)", last_record_2):
#         if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record_2):
#             if re.search(r"(^\w+.)(\w+)(.)(2$)", last_record_1):
#                 if re.search(r"(^\w+.)(\w+)(.)(2$)", first_record_1):
#                     logging.warning("Read 1 and Read 2 appear to be flipped. Switching inputs...")
#                     switched_paired_end_files(R1_dict, R2_dict, subsample_max)
                                        
    
# # As long as none of the exceptions above match, it moves on to downsample or not to downsample paired end 
#     if re.search(r"(^\w+.)(\w+)(.)(1$)", last_record_1):
#         if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record_1):
#             if re.search(r"(^\w+.)(\w+)(.)(2$)", last_record_2):
#                 if re.search(r"(^\w+.)(\w+)(.)(2$)", first_record_2):
#                     paired_end_files(R1_dict, R2_dict, subsample_max)
                    

# '''
# Program has finished without issue. 
# '''
# logging.info('Sylens has finished processing. Closing...  ')
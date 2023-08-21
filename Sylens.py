#!/usr/bin/env python3
'''
TODO Far future: Maybe add sam/bam file
DONE Come up with name :D :D :D 
DONE Add ability to compress as a flag
DONE Create repo in git for this project
IN_PROCESS Learn about exceptions
DONE? Logging library https://docs.python.org/3/library/logging.html#logging-levels
'''

#SYLENS: Sampling Yielding LEss Noticeable Samples

###Importing libraries### 
import argparse
from Bio import SeqIO
import re
import random
import sys 
import time
import gzip
import shutil
import os
import logging

'''
Creating a bunch of custom functions for downsampling and for not downsampling, zipping, etc. Just makes them easier to reference later
'''

# function to unzip file 
def unzip(input_file): 
    logging.info('File is compressed. Opening...')
    with gzip.open(input_file, 'rb') as fileobj:
        with open('out_file.fastq', 'wb') as output:
            shutil.copyfileobj(fileobj, output)
            
def unzip_2(input_file_2):
    logging.info('Second file is compressed. Opening...')
    with gzip.open(input_file_2, 'rb') as fileobj:
        with open('out_file_2.fastq', 'wb') as output:
            shutil.copyfileobj(fileobj, output)



'''A generated dictionary is randomly sampled. The randomly sampled IDs are saved in a list. Using that list, the IDs are mapped back to
their respective entires in the called dictionary. SeqIO writes the output file with the filetype indicated by the user. 
'''
#function to downsample single files
def downsample_single_end_file(dictionary_name, subsample_max):
    logging.info("File is single-end.")
    IDs_Seq_File = open(f"down_sampled_{args.Read1[0]}", 'w')
    Random_IDs = random.sample(list(dictionary_name), subsample_max)
    single_file_output = [dictionary_name[info] for info in Random_IDs]
    logging.info('Writing to file...')
    SeqIO.write(single_file_output, IDs_Seq_File, args.output)
    logging.info(f"Done writing to: down_sampled_{args.Read1[0]}")


'''
Regular expressions store R1/R2 IDs into separate lists. Lists zipped to randomly sample tuples. Rezip randomly sampled 
IDs into separate lists. Search their respective dictionary. Pull input that matches R1/R2's randomly sampled IDS. Store as
output. Separate interleaved fastq with R1 first and R2 second using counts for R1 and R2 to index the output list. Write indexed value
to output file. Count = 1 so R1 record is written first. Added 1 total length since the count < int is exclusive.
'''
#function to downsample interleaved files
def downsample_interleaved_file(dictionary_name, subsample_max):
    logging.info("File is interleaved. Processing interleaved file...")
    IDs_Seq_File = open(f"down_sampled_{args.Read1[0]}", 'w')
    for element in list(dictionary_name):
        if re.search(r"(^\w+.)(\w+)(.)(2$)", element):
            interleaved_R2_IDs.append(element)
        if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
            interleaved_R1_IDs.append(element)
    zipped_ids_list = list(zip(interleaved_R1_IDs, interleaved_R2_IDs))
    random_ids = random.sample(zipped_ids_list, subsample_max)
    r1_ids, r2_ids = zip(*random_ids)
    output_1 = [dictionary_name[info] for info in r1_ids]
    output_2 = [dictionary_name[info] for info in r2_ids]  
    logging.info('Writing to file...')
    count = 1
    R1_count = 0
    R2_count = 0
    while count < (1 + 2*len(output_1)):
        if count % 2 != 0 :
            count += 1
            SeqIO.write(output_1[R1_count], IDs_Seq_File, args.output)
            R1_count += 1
        else:
            if count % 2 == 0:
                count+= 1
                SeqIO.write(output_2[R2_count], IDs_Seq_File, args.output)
                R2_count += 1
    logging.info(f"Done writing to: down_sampled_{args.Read1[0]}")


'''
Creates two lists and searches for every .1 and .2 element in the given dictionary. Then creates two separate files. A list of all of the IDs
are zipped together and randomly sampled as tuples. Rezipped apart into two lists. Saves the output and writes it to the filetype indicated by
the user in two distinct separate files. 
'''
#function to downsample paired end files
def downsample_paired_end_files(dictionary_name_1, dictionary_name_2, subsample_max):
    logging.info("Read 1 and Read 2 supplied. Downsampling starting...")
    for element in list(dictionary_name_1):
        if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
            paired_end_1.append(element)
    for element in list(dictionary_name_2): 
        if re.search(r"(^\w+.)(\w+)(.)(2$)", element):
            paired_end_2.append(element)
    IDs_Seq_File = open(f"down_sampled_{args.Read1[0]}", 'w')
    IDs_Seq_File_2 = open(f"down_sampled_{args.Read2}", 'w')
    zipped_ids_list = list(zip(paired_end_1, paired_end_2))
    random_ids = random.sample(zipped_ids_list, subsample_max)
    paired_end_1_ids, paired_end_2_ids = zip(*random_ids)
    output_1 = [dictionary_name_1[info] for info in paired_end_1_ids]
    output_2 = [dictionary_name_2[info] for info in paired_end_2_ids]  
    logging.info('Writing to file...')
    SeqIO.write(output_1, IDs_Seq_File, args.output)
    SeqIO.write(output_2, IDs_Seq_File_2, args.output)
    logging.info(f"Done writing to: down_sampled_{args.Read1[0]} and down_sampled_{args.Read2}") 

#function to downsample switched paired end file
def switched_downsample_paired_end_files(dictionary_name_1, dictionary_name_2, subsample_max):
    logging.info("Reads 1 and 2 supplied. Downsampling occurring. Processing Read 1 and Read 2 files...")
    for element in list(dictionary_name_2):
        if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
            paired_end_1.append(element)
    for element in list(dictionary_name_1): 
        if re.search(r"(^\w+.)(\w+)(.)(2$)", element):
            paired_end_2.append(element)
    IDs_Seq_File_2 = open(f"down_sampled_{args.Read1[0]}", 'w')
    IDs_Seq_File = open(f"down_sampled_{args.Read2}", 'w')
    zipped_ids_list = list(zip(paired_end_1, paired_end_2))
    random_ids = random.sample(zipped_ids_list, subsample_max)
    paired_end_1_ids, paired_end_2_ids = zip(*random_ids)
    output_2 = [dictionary_name_1[info] for info in paired_end_2_ids]
    output_1 = [dictionary_name_2[info] for info in paired_end_1_ids]  
    logging.info('Writing to file...')
    SeqIO.write(output_1, IDs_Seq_File, args.output)
    SeqIO.write(output_2, IDs_Seq_File_2, args.output)
    logging.info(f"Done writing to: down_sampled_{args.Read1[0]} and down_sampled_{args.Read2}") 


'''
A generated dictionary is stored in list form, the IDS are saved. The IDs are then mapped back to the dicionary and written out as the user 
indicated. The All_IDs/All_output step is used to change the entry to a seqrecord object to make it usedable with SeqIO.
''' 
#function to not downsample single end file
def no_downsample_single_end_file(dictionary_name):
    logging.info("File is single-end. No downsampling occurring")
    IDs_Seq_File = open(f"non_down_sampled_{args.Read1[0]}", 'w')
    All_IDs = list(dictionary_name)
    All_output = [dictionary_name[info] for info in All_IDs]
    logging.info('Writing to file...')
    SeqIO.write(All_output, IDs_Seq_File, args.output)
    logging.info(f'Done writing to: non_down_sampled_{args.Read1[0]}')


'''
For all of the elements in the dictionary, if it ends with .2, R2 list IDs are saved, if .1, R1 IDs are saved. The lists are then used to 
create seqrecord objects by mapping the information from the indicated dictionary back to the corresponding IDs. Used count to index the
lists and write every other function to the file
'''
#function to not downsample interleaved file
def no_downsample_interleaved_file(dictionary_name):
    logging.info("File is interleaved. No downsampling occurring. Processing interleaved file...")
    IDs_Seq_File = open(f"non_down_sampled_{args.Read1[0]}", 'w')
    for element in list(dictionary_name):
        if re.search(r"(^\w+.)(\w+)(.)(2$)", element):
            interleaved_R2_IDs.append(element)
        if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
            interleaved_R1_IDs.append(element)
    output_1 = [dictionary_name[info] for info in interleaved_R1_IDs]
    output_2 = [dictionary_name[info] for info in interleaved_R2_IDs]  
    logging.info('Writing to file...')
    count = 1
    R1_count = 0
    R2_count = 0
    while count < (1 + 2*len(output_1)):
        if count % 2 != 0 :
            count += 1
            SeqIO.write(output_1[R1_count], IDs_Seq_File, args.output)
            R1_count += 1
        else:
            if count % 2 == 0:
                count+= 1
                SeqIO.write(output_2[R2_count], IDs_Seq_File, args.output)
                R2_count += 1
    logging.info(f"Done writing to: non_down_sampled_{args.Read1[0]}")


'''
Creates two lists and searches for every .1 and .2 element in the dictionary. Two lists are mapped back to their corresponding dictionaries 
to produce seqrecord objects. Saves the output and writes it to the filetype indicated by the user in two distinct separate files.
'''
#function to not downsample paired end files
def no_downsample_paired_end_files(dictionary_name_1, dictionary_name_2):
    logging.info("Read 1 and Read 2 supplied. No downsampling occurring. Processing Read 1 and Read 2 files...")
    for element in list(dictionary_name_1):
        if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
            paired_end_1.append(element)
    for element in list(dictionary_name_2): 
        if re.search(r"(^\w+.)(\w+)(.)(2$)", element):
            paired_end_2.append(element)
    IDs_Seq_File = open(f"non_down_sampled_{args.Read1[0]}", 'w')
    IDs_Seq_File_2 = open(f"non_down_sampled_{args.Read2}", 'w')
    output_1 = [dictionary_name_1[info] for info in paired_end_1]
    output_2 = [dictionary_name_2[info] for info in paired_end_2]  
    logging.info('Writing to file...')
    SeqIO.write(output_1, IDs_Seq_File, args.output)
    SeqIO.write(output_2, IDs_Seq_File_2, args.output)
    logging.info(f"Done writing to: non_down_sampled_{args.Read1[0]} and non_down_sampled_{args.Read2}")


#function to not downsample switched paired end files
def switched_no_downsample_paired_end_files(dictionary_name_1, dictionary_name_2):
    logging.info("Reads 1 and 2 supplied. No downsampling occurring. Processing Read 1 and Read 2 files...")
    for element in list(dictionary_name_1):
        if re.search(r"(^\w+.)(\w+)(.)(2$)", element):
            paired_end_2.append(element)
    for element in list(dictionary_name_2): 
        if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
            paired_end_1.append(element)
    IDs_Seq_File_2 = open(f"non_down_sampled_{args.Read1[0]}", 'w')
    IDs_Seq_File = open(f"non_down_sampled_{args.Read2}", 'w')
    output_2 = [dictionary_name_1[info] for info in paired_end_2]
    output_1 = [dictionary_name_2[info] for info in paired_end_1]  
    logging.info('Writing to file...')
    SeqIO.write(output_1, IDs_Seq_File, args.output)
    SeqIO.write(output_2, IDs_Seq_File_2, args.output)
    logging.info(f"Done writing to: non_down_sampled_{args.Read2} and non_down_sampled_{args.Read1[0]}")


#function to zip files
def gzip_output_file(original_file, new_output_file):
    with open(original_file, 'rb') as input_file:
        with gzip.open(new_output_file, 'wb') as output_file:
            shutil.copyfileobj(input_file, output_file)

#function to gzip a file if -c == yes ###
def compress(original_file, new_output_file):
    if args.compress == 'yes': 
        if args.Read1[0].endswith('.fastq'):
            gzip_output_file(original_file, new_output_file)
            os.remove(original_file)      

#function to remove generated files
def remove_nonsense_files():
    if os.path.exists('out_file.fastq'):
        os.remove('out_file.fastq')
    if os.path.exists('out_file_2.fastq'):
        os.remove('out_file_2.fastq')


### Defining a few empty lists ###
interleaved_R1_IDs = []
interleaved_R2_IDs = []
r1_IDS = []
r2_IDS = []
r1_split_ids = []
r2_split_ids = []
paired_end_1 = []
paired_end_2 = []

### Creating the parser to take user inputs and adding a description of the input ###
# Downsampler program description 
parser = argparse.ArgumentParser( 
    prog = 'Subsampler for FASTQ file(s)',
    description= 'Enter in FASTQ file(s) and down sample based on a user supplied integer. If no user input is added the entire file is sampled and output as chosen filetype.',
    epilog = 'This is supposed to show the epilog at the bottom of help, but I do not know how to write an epilog, haha')


# Creating a function which takes in a filename.fastq or R1.fastq. The number of arguments needed here is 1
parser.add_argument(
    'Read1', 
    nargs = 1,
    help = 'Add FASTQ file here. I.E. filename.fastq or R1.fastq'
    )


# Creating optional argument which can be R2.fastq for situations where separate paired-end files are. The number of arguments here can be 0 or more
parser.add_argument(
    'Read2', 
    nargs = '?',
    help = 'Add second FASTQ file here if paired-end. I.E. R2.fastq'
    )


# Creating input which allows the user to input an integer for # of samples to downsample to 
parser.add_argument('-s','--subsample', 
    type = int,
    help = 'Enter an integer which will be the total number of down sampling of FASTQ files occuring. I.E. -s 10000. Leave blank if no subsampling needs to occur and only file converion needs to occur.'
    )


# Creating argument to denote what type of fastq file is being looked at
parser.add_argument('-f', '--filetype',
    choices = ['fastq', 'fastq-solexa'],
    default = 'fastq',
    help = "Add what type of fastq file is being input. fastq/fastq-sanger-> uses ASCII offset of 33 whereas fastq-solexa-> uses ASCII of 64. I.E fastq-solexa"
    )


# Creating part of the argument which allows the user to input an int for subsampling 
parser.add_argument('--seed', 
    type = int,
    default = int(time.time()),
    help = 'Enter the seed number if you would like to reproduce previous results. I.E. --seed 1691696502'
    )


# Creating argument to denote what type of fastq file is wanted during output
parser.add_argument('-o', '--output',
    choices = ['fastq', 'fastq-solexa'],
    default = 'fastq',
    help = "Add what type of fastq file is desired at output. I.E fastq-solexa"
    )

# Creating argument to denote what type of fastq file is wanted during output
parser.add_argument('-c', '--compress',
    choices = ['yes', 'no'],
    default = 'no',
    help = "Compress fastq file into fastq.gz file on output. I.E. --c no"
    )

args = parser.parse_args()


logging.basicConfig(level = logging.INFO, format = '%(levelname)s - %(message)s')


### Set seed value to time since epoch that way each time program is run, ###
### seed value is different. If user enters seed value, program will ###
### reproduce the same results ###
random.seed(args.seed)

### These statements check to make sure user input has been pulled and allows to user to know their input was taken###
### They are not necessary, just somewhat helpful ### 
if args.Read2 == None:
    logging.info('The FASTQ filename is: %s', ''.join(args.Read1))
else:
    logging.info('The FASTQ filenames are: %s', (args.Read1[0], args.Read2))

logging.info('The input file type is: %s', args.filetype)

logging.info("The run's seed number is: %s", args.seed)

logging.info('The output file type is: %s', args.output)

logging.info('Compressing on output: %s', args.compress)

#Ensures that the subsample integer is pulled in as a usable int not str variable or to allow for no subsampling.
if args.subsample != None:
    subsample_max = int(args.subsample)
    logging.info('The amount to subsample is: %s', args.subsample)
    logging.debug('Checking to see if subsample size is larger than file')
    if args.Read1[0].endswith('.gz'):
        unzip(args.Read1[0])
        Dict_for_subsample_max = SeqIO.to_dict(SeqIO.parse('out_file.fastq', args.filetype))
    else:
        Dict_for_subsample_max = SeqIO.to_dict(SeqIO.parse(args.Read1[0], args.filetype))
    subsample_max_IDs = []
    for element in list(Dict_for_subsample_max):
        if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
            subsample_max_IDs.append(element)
    total_length = len(subsample_max_IDs)
    if total_length < subsample_max:
        logging.critical('Subsample size is larger than total number of samples in file. Terminating...')
        sys.exit(1)
    
if args.subsample == None:
    if args.Read1[0].endswith('.gz'):
        unzip(args.Read1[0])
        Dict_for_subsample_max = SeqIO.to_dict(SeqIO.parse('out_file.fastq', args.filetype))
    else:
        Dict_for_subsample_max = SeqIO.to_dict(SeqIO.parse(args.Read1[0], args.filetype))
    subsample_max_IDs = []
    subsample_max_IDs_2 = []
    for element in list(Dict_for_subsample_max):
        if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
            subsample_max_IDs.append(element)
        else:
            if re.search(r"(^\w+.)(\w+)(.)(2$)", element):    
                subsample_max_IDs_2.append(element)
    subsample_max = len(subsample_max_IDs)
    if subsample_max == 0:
        subsample_max = len(subsample_max_IDs_2)
    logging.info('The total amount of sampling occurring: %s', subsample_max)
    




############################################### Only Read 1 supplied ###############################################


# This sees that no read 2 was entered. Then converts given fastq file into dictionary. Grabs the first and last records.
if args.Read2 == None:
    logging.info("No read 2 entered, checking if single-end or interleaved file...")
    if args.Read1[0].endswith('.gz'):
        unzip(args.Read1[0])
        R1_dict = SeqIO.to_dict(SeqIO.parse('out_file.fastq', args.filetype))
    else:
        R1_dict =  SeqIO.to_dict(SeqIO.parse(args.Read1[0], args.filetype))
    last_record = list(R1_dict) [-1]
    first_record = list(R1_dict) [0]


# If regular expression finds that both the first and last record end with '.2', it will terminate the program
    if re.search(r"(^\w+.)(\w+)(.)(2$)", first_record):
        logging.critical("This appears to be a read 2 file instead of a read 1 file. Program terminating...")
        sys.exit(1)    


# Regular expression looks if an interleaved file was supplied. If first entry ends with ".1" and the last record ends with a ".2", interleaved file. If a subsample size was added, the
# downsample_interleaved_file function is run. If a subsample size is not indicated, the no_downsample_interleaved_file is run. 
    elif re.search(r"(^\w+.)(\w+)(.)(2$)", last_record):
        if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record):
            if args.subsample != None:
                downsample_interleaved_file(R1_dict, subsample_max)
                if args.Read1[0].endswith('.gz'):
                    gzip_output_file(f'down_sampled_{args.Read1[0]}', f'down_sampled_{args.Read1[0]}')
                compress(f'down_sampled_{args.Read1[0]}', f'down_sampled_{args.Read1[0]}')
            if args.subsample == None:
                no_downsample_interleaved_file(R1_dict)
                if args.Read1[0].endswith('.gz'):
                    gzip_output_file(f'non_down_sampled_{args.Read1[0]}', f'non_down_sampled_{args.Read1[0]}.gz')
                compress(f'non_down_sampled_{args.Read1[0]}', f'non_down_sampled_{args.Read1[0]}.gz')
            
 
# Uses a regular expression to look at the first and last entry of the fastq. If first and last entry end with ".1", single end file. Moves into whether or not downsampling is occurring
    elif re.search(r"(^\w+.)(\w+)(.)(1$)", last_record):
        logging.info("Making sure file is single-end...")
        if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record):
            if args.subsample != None:
                downsample_single_end_file(R1_dict, subsample_max)
                if args.Read1[0].endswith('.gz'):
                    gzip_output_file(f'down_sampled_{args.Read1[0]}', f'down_sampled_{args.Read1[0]}.gz')
                    os.remove(f'down_sampled_{args.Read1[0]}')
                    os.remove('out_file.fastq')
                compress(f'down_sampled_{args.Read1[0]}', f'down_sampled_{args.Read1[0]}.gz')
            if args.subsample == None:
                no_downsample_single_end_file(R1_dict)
                if args.Read1[0].endswith('.gz'):
                    gzip_output_file(f'non_down_sampled_{args.Read1[0]}', f'non_down_sampled_{args.Read1[0]}.gz')
                    os.remove(f'non_down_sampled_{args.Read1[0]}')
                    os.remove('out_file.fastq')
                compress(f'non_down_sampled_{args.Read1[0]}', f'non_down_sampled_{args.Read1[0]}.gz')






############################################### Both Read 1 and 2 supplied ###############################################


# Two files supplied. Looks if a read or both are .gz, decompresses them, and assigns each file to a dictionary. Then grabs the first and last IDs of both files. 
if args.Read2 != None:
    if args.Read1[0].endswith('.gz'):
        unzip(args.Read1[0])
        R1_dict =  SeqIO.to_dict(SeqIO.parse('out_file.fastq', args.filetype))
    else:
        R1_dict =  SeqIO.to_dict(SeqIO.parse(args.Read1[0], args.filetype)) 
    if args.Read2.endswith('.gz'):  
        unzip_2(args.Read2)
        R2_dict =  SeqIO.to_dict(SeqIO.parse('out_file_2.fastq', args.filetype))
    else:
        R2_dict =  SeqIO.to_dict(SeqIO.parse(args.Read2, args.filetype))
    last_record_1 = list(R1_dict) [-1]
    first_record_1 = list(R1_dict) [0]
    last_record_2 = list(R2_dict) [-1]
    first_record_2 = list(R2_dict) [0]


#Errors if wrong files are supplied and deletes files that were made preparing for decompression, if applicable.
    if re.search(r"(^\w+.)(\w+)(.)(2$)", last_record_1):
        if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record_1):
            if re.search(r"(^\w+.)(\w+)(.)(2$)", last_record_2):
                if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record_2):
                    logging.critical("Both supplied files appear to be interleaved. Program terminating...")
                    remove_nonsense_files()
                    sys.exit(1)
    if re.search(r"(^\w+.)(\w+)(.)(1$)", last_record_1):
        if re.search(r"(^\w+.)(\w+)(.)(1$)", last_record_2):
            logging.critical("Both supplied files appear to be Read 1 files. Program terminating...")
            remove_nonsense_files()
            sys.exit(1)
    if re.search(r"(^\w+.)(\w+)(.)(2$)", last_record_1):
        if re.search(r"(^\w+.)(\w+)(.)(2$)", first_record_1):
            if re.search(r"(^\w+.)(\w+)(.)(2$)", last_record_2):
                if re.search(r"(^\w+.)(\w+)(.)(2$)", first_record_2):
                    logging.critical("Both supplied files appear to be Read 2 files. Program terminating...")
                    remove_nonsense_files()
                    sys.exit(1)
    if re.search(r"(^\w+.)(\w+)(.)(2$)", last_record_1):
        if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record_1):
            logging.critical("First input file is interleaved. Program terminating...")
            remove_nonsense_files()
            sys.exit(1)
    if re.search(r"(^\w+.)(\w+)(.)(2$)", last_record_2):
        if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record_2):
            logging.critical("Second input file is interleaved. Program terminating...")
            remove_nonsense_files()
            sys.exit(1)
    
    

# Switched input files R2 first instead of R1
    if re.search(r"(^\w+.)(\w+)(.)(1$)", last_record_2):
        if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record_2):
            if re.search(r"(^\w+.)(\w+)(.)(2$)", last_record_1):
                if re.search(r"(^\w+.)(\w+)(.)(2$)", first_record_1):
                    logging.warning("Read 1 and Read 2 appear to be flipped. Switching inputs...")
                    if args.subsample != None:
                        switched_downsample_paired_end_files(R1_dict, R2_dict, subsample_max)
                        if args.Read1[0].endswith('.gz'):
                            gzip_output_file(f'down_sampled_{args.Read2}', f'down_sampled_{args.Read2}.gz')
                            os.remove(f'down_sampled_{args.Read2}')
                            os.remove(f'out_file.fastq')
                        if args.Read2.endswith('.gz'):
                            gzip_output_file(f'down_sampled_{args.Read1[0]}', f'down_sampled_{args.Read1[0]}.gz')
                            os.remove(f'down_sampled_{args.Read1[0]}')
                            os.remove('out_file_2.fastq')
                        compress(f'down_sampled_{args.Read1[0]}', f'down_sampled_{args.Read1[0]}.gz')
                        compress(f'down_sampled_{args.Read2}', f'down_sampled_{args.Read2}.gz')
                    if args.subsample == None:
                        switched_no_downsample_paired_end_files(R1_dict, R2_dict)
                        if args.Read1[0].endswith('.gz'):
                            gzip_output_file(f'non_down_sampled_{args.Read2}', f'non_down_sampled_{args.Read2}.gz')
                            os.remove(f'non_down_sampled_{args.Read2}')
                            os.remove('out_file.fastq')
                        if args.Read2.endswith('.gz'):
                            gzip_output_file(f'non_down_sampled_{args.Read1[0]}', f'non_down_sampled_{args.Read1[0]}.gz')
                            os.remove(f'non_down_sampled_{args.Read1[0]}')
                            os.remove('out_file_2.fastq')
                        compress(f'non_down_sampled_{args.Read1[0]}', f'non_down_sampled_{args.Read1[0]}.gz')
                        compress(f'non_down_sampled_{args.Read2}', f'non_down_sampled_{args.Read2}.gz')

    
# As long as none of the exceptions above match, it moves on to downsample or not to downsample paired end 
    if re.search(r"(^\w+.)(\w+)(.)(1$)", last_record_1):
        if re.search(r"(^\w+.)(\w+)(.)(1$)", first_record_1):
            if re.search(r"(^\w+.)(\w+)(.)(2$)", last_record_2):
                if re.search(r"(^\w+.)(\w+)(.)(2$)", first_record_2):
                    if args.subsample != None:
                        downsample_paired_end_files(R1_dict, R2_dict, subsample_max)
                        if args.Read1[0].endswith('.gz'):
                            gzip_output_file(f'down_sampled_{args.Read1[0]}', f'down_sampled_{args.Read1[0]}.gz')
                            os.remove(f'down_sampled_{args.Read1[0]}')
                            os.remove('out_file.fastq')
                        if args.Read2.endswith('.gz'):
                            gzip_output_file(f'down_sampled_{args.Read2}', f'down_sampled_{args.Read2}.gz')
                            os.remove(f'down_sampled_{args.Read2}')
                            os.remove('out_file_2.fastq')
                        compress(f'down_sampled_{args.Read1[0]}', f'down_sampled_{args.Read1[0]}.gz')
                        compress(f'down_sampled_{args.Read2}', f'down_sampled_{args.Read2}.gz')
                    if args.subsample == None:
                        no_downsample_paired_end_files(R1_dict,R2_dict)
                        if args.Read1[0].endswith('.gz'):
                            gzip_output_file(f'non_down_sampled_{args.Read1[0]}', f'non_down_sampled_{args.Read1[0]}.gz')
                            os.remove(f'non_down_sampled_{args.Read1[0]}')
                            os.remove('out_file.fastq')
                        if args.Read2.endswith('.gz'):
                            gzip_output_file(f'non_down_sampled_{args.Read2}', f'non_down_sampled_{args.Read2}.gz')
                            os.remove(f'non_down_sampled_{args.Read2}')
                            os.remove('out_file_2.fastq')
                        compress(f'non_down_sampled_{args.Read1[0]}', f'non_down_sampled_{args.Read1[0]}.gz')
                        compress(f'non_down_sampled_{args.Read2}', f'non_down_sampled_{args.Read2}.gz')


### Everything is doooooonnnnnnnnnnneeeeeeeeeeeeeeeeeeeeeeeeeee ###
logging.info('Yay, all finished! :)')
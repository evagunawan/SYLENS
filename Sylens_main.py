#!/usr/bin/env python3

'''
TODO Flesh out readme
TODO GNU GPLv3* license research
TODO reading up containerize Sylens 
TODO Put together request for new laptop
'''

#SYLENS: Sampling Yielding LEss Noticeable Samples

###Importing libraries### 
import argparse
import random
import sys 
import time
import logging

from read_fastq_file import FastqFileData


# Creating class that displays help if no information is entered in parser and changes how the error is displayed when something is entered wrong
class DownsamplerParser(argparse.ArgumentParser):

    def error(self, message):

        self.print_help()
        sys.stderr.write(f'\nERROR DETECTED: {message}\n')

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
    help = "Compress fastq file into fastq.gz file on output. I.E. -c"
    )
parser.add_argument('-p', '--percentage',
    default = False,
    action = 'store_true',
    help = "With the -p flag, subsampling is done as a percentage of reads instead of an indicated number of reads. Percentage of reads should be an integer between 1-100. I.E. -p -s 15"
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

logging.info(f"This run's seed number is {args.seed}")

logging.debug('Starting processing of file(s)')

#Created fastq file object using the parameters specified in secondary script 
fastq_data_object = FastqFileData(args.Read1, args.Read2, args.subsample, args.output, args.compress, args.filetype, args.seed, args.percentage)

#Analyzing fastq_information_object
logging.debug('Starting reading_fastq_file from main script')
fastq_data_object.reading_fastq_file()

logging.debug('Starting determine_Fastq_ID_formatting from main script')
fastq_data_object.determine_fastq_ID_formatting()

logging.debug('Starting determine_paired_single_interleaved from main script')
fastq_data_object.determine_paired_single_interleaved()

logging.debug('Starting processing_filetype from main script')
fastq_data_object.processing_filetype()

logging.info('Sylens has finished processing files. Closing.')

sys.exit()
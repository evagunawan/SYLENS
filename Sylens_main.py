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
import random
import sys 
import time
import logging

import read_fastq_file

from read_fastq_file import FastqFile
from single_end_processing import process_single_end_sampling
from interleaved_processing import process_interleaved_sampling
from paired_end_processing import process_paired_end_sampling

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
    help = "Compress fastq file into fastq.gz file on output. I.E. -c"
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
fastq_information_object = FastqFile(args.Read1, args.Read2, args.filetype)

fastq_information_object.reading_fastq_file()

#Writing output of determined filetype
if read_fastq_file.determined_filetype == 'Single-end':

    process_single_end_sampling(args.Read1, args.Read2, args.subsample, args.output, args.compress)

elif read_fastq_file.determined_filetype == 'Interleaved':

    process_interleaved_sampling(args.Read1, args.subsample, args.output, args.compress)

elif read_fastq_file.determined_filetype == 'Paired-end':

    process_paired_end_sampling(args.Read1, args.Read2, args.subsample, args.output, args.compress)

logging.info('Sylens has finished processing files. Closing')

sys.exit()
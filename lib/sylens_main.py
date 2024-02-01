#SYLENS: Sampling Yielding LEss Noticeable Samples

###Importing libraries### 
import argparse
import sys 
import time
import logging
import os

from lib.subsampling import subsample_single, subsample_paired
from lib.read_fastq_file import FastqFileData
from lib.write_output_file import write_single, write_interleaved_paired_end



# Creating class that displays help if no information is entered in parser and changes how the error is displayed when something is entered wrong
class DownsamplerParser(argparse.ArgumentParser):

    def error(self, message):

        self.print_help()
        sys.stderr.write(f'\nERROR DETECTED: {message}\n')

        sys.exit(1)

def main():
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
    parser.add_argument('-o', '--outputFormat',
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
        type = int,
        help = "Subsampling is done as a percentage of reads instead of an indicated number of reads. Percentage of reads should be an integer between 1-100. I.E. -p 15"
        )
    parser.add_argument('--output_type',
        choices = ['separate', 'joined', 'interleaved'],
        default = 'separate',
        help = "While writing the output file for interleaved and paired end files, output can be written with all R1 reads first and all R2 reads after (joined), written in different files (separate), or written with an alternating R1 first and R2 second pattern (interleaved). By default the output will be separate."
        )

    #Runs the parser and allows you to call arguments downstream
    args = parser.parse_args()

    #Format for logging debug-critical information
    logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

    '''
    Set seed value to time since e`poch that way each time program is run, 
    seed value is different. If user enters seed value, program will 
    reproduce the same results 
    '''
    
    logging.info(f"This run's seed number is {args.seed}")

    logging.debug('Starting processing of file(s)')

    #Created fastq file object using the parameters specified in secondary script 
    fastq_data_object = FastqFileData(args.Read1, args.Read2, args.filetype)

    # If no subsampling is desired
    if args.subsample == None and args.percentage == None:
        args.percentage = 100

    # Grabbing just the file name if a file path is included.
    args.Read1 = os.path.basename(args.Read1)
    if args.Read2:
        args.Read2 = os.path.basename(args.Read2)
    
    #Creating file output name
    if args.percentage == 100:

        file_naming_convention_1 = f'non_downsampled_{args.Read1}'    

        if args.Read2 != None:
            file_naming_convention_2 = f'non_downsampled_{args.Read2}'  
        else:
            file_naming_convention_2 = f'non_downsampled_Read2_{args.Read1}' 

    else:

        file_naming_convention_1 = f'{args.seed}_downsampled_{args.Read1}'

        if args.Read2 != None:
            file_naming_convention_2 = f'{args.seed}_downsampled_{args.Read2}'
        else:
            file_naming_convention_2 = f'{args.seed}_downsampled_Read2_{args.Read1}'

    #If percentage given determine subsample level
    if args.percentage:
        subsample_level = int((fastq_data_object.R1_Total_Reads)*(args.percentage) / (100))

    else:
        subsample_level = args.subsample

    logging.debug(f'The amount to subsample is: {subsample_level}.')
    logging.debug(f'The total amount of reads: {fastq_data_object.R1_Total_Reads}.')

    if fastq_data_object.Interleaved_Read1_IDs:
        Read1_IDs, Read2_IDs = subsample_paired(
            fastq_data_object.Interleaved_Read1_IDs,
            fastq_data_object.Interleaved_Read2_IDs,
            subsample_level,
            args.seed
        )
        write_interleaved_paired_end(fastq_data_object, Read1_IDs, Read2_IDs, args.outputFormat, args.compress, file_naming_convention_1, file_naming_convention_2, args.output_type)

    elif fastq_data_object.Read2Index:
        Read1_IDs, Read2_IDs = subsample_paired(
            list(fastq_data_object.Read1Index.keys()),
            list(fastq_data_object.Read2Index.keys()),
            subsample_level,
            args.seed
        )
        write_interleaved_paired_end(fastq_data_object, Read1_IDs, Read2_IDs, args.outputFormat, args.compress, file_naming_convention_1, file_naming_convention_2, args.output_type)
    
    else:
        Read_IDs = subsample_single(
            list(fastq_data_object.Read1Index.keys()),
            subsample_level,
            args.seed
        )
        write_single(fastq_data_object, Read_IDs, args.outputFormat, args.compress, file_naming_convention_1)
    fastq_data_object.cleanUP()
    sys.exit(0)
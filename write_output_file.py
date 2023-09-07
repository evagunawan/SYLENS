#!/usr/bin/env python3
import gzip
import logging

from Bio import SeqIO

def no_compress(file_name, object_output, argsOutput):

    with open(f'{file_name}', 'wt') as IDs_Seq_File:

        SeqIO.write(object_output, IDs_Seq_File, argsOutput)  

def compress_gz(file_name, output_object, argsOutput,):

    with gzip.open(f'{file_name}', 'wt') as IDs_Seq_File:

        SeqIO.write(output_object, IDs_Seq_File, argsOutput)

def compress_interleaved(argsOutput, output_1, output_2, input_file_name):

    with gzip.open(f"{input_file_name}", 'wt') as IDs_Seq_File:

        logging.debug('Writing to gzipped file.')

        #Writing every other entry to output file
        while count < (1 + 2*len(output_1)):

            if count % 2 != 0 :

                count += 1
                SeqIO.write(output_1[R1_count], IDs_Seq_File, argsOutput)
                R1_count += 1

            else:

                if count % 2 == 0:

                    count+= 1
                    SeqIO.write(output_2[R2_count], IDs_Seq_File, argsOutput)
                    R2_count += 1

def write_interleaved_file(argsRead1, argsOutput, argsCompress, output_1, output_2, input_file_name):

    logging.debug('In writing interleaved file')

    count = 1
    R1_count = 0
    R2_count = 0

    #Compression desired
    if argsRead1.endswith('.gz'):

        compress_interleaved(argsOutput, output_1, output_2, input_file_name)

    if argsCompress == True and not argsRead1.endswith('.gz'):
        
        input_file_name = input_file_name + '.gz'

        compress_interleaved(argsOutput, output_1, output_2, input_file_name)

    #No compression desired
    if argsCompress != True and not argsRead1.endswith('.gz'):

        with open(f'{input_file_name}', 'wt') as IDs_Seq_File:

            logging.debug('Writing to file.')

            while count < (1 + 2*len(output_1)):
                    
                if count % 2 != 0 :
                    
                    count += 1
                    SeqIO.write(output_1[R1_count], IDs_Seq_File, argsOutput)
                    R1_count += 1
                    
                else:
                    
                       if count % 2 == 0:

                        count+= 1
                        SeqIO.write(output_2[R2_count], IDs_Seq_File, argsOutput)
                        R2_count += 1

def write_output_file(argsRead1, argsRead2, argsCompress, file_name_1, file_name_2, output_1, output_2, argsOutput):

    if argsRead2 == None:

        logging.debug('Entered into argsread2 == None in write_output_file')

        if argsRead1.endswith('.gz'):

            compress_gz(file_name_1, output_1, argsOutput)

        elif argsCompress == True and not argsRead1.endswith('.gz'):

            file_name_1 = file_name_1 + '.gz'

            compress_gz(file_name_1, output_1, argsOutput)

        elif argsCompress != True:

            no_compress(file_name_1, output_1, argsOutput)               

    if argsRead2 != None:

        logging.debug('Entered into argsread2 != None in write_output_file')

        if argsRead1.endswith('.gz'):

            compress_gz(file_name_1, output_1, argsOutput)

        if argsRead2.endswith('.gz'):
            
            compress_gz(file_name_2, output_2, argsOutput)

        if argsCompress == True and not argsRead1.endswith('.gz'):

            file_name_1 = file_name_1 + '.gz'
            compress_gz(file_name_1, output_1, argsOutput)

        if argsCompress == True and not argsRead2.endswith('.gz'):

            file_name_2 = file_name_2 + '.gz'
            compress_gz(file_name_2, output_2, argsOutput)

        if argsCompress != True and not argsRead1.endswith('.gz'):

            no_compress(file_name_1, output_1, argsOutput)

        if argsCompress != True and not argsRead2.endswith('.gz'):

            no_compress(file_name_2, output_2, argsOutput)     

    logging.debug('Finished writing to file.')
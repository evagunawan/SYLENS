#!/usr/bin/env python3
import gzip
import logging

from Bio import SeqIO

#If no compression is desired
def no_compress(file_name, object_output, argsOutput):

    with open(f'{file_name}', 'wt') as IDs_Seq_File:

        SeqIO.write(object_output, IDs_Seq_File, argsOutput)  

#If compression is desired
def compress_gz(file_name, output_object, argsOutput,):

    with gzip.open(f'{file_name}', 'wt') as IDs_Seq_File:

        SeqIO.write(output_object, IDs_Seq_File, argsOutput)

#If file is an interleaved file that needs to be compressed 
def compress_interleaved(argsOutput, output_1, output_2, input_file_name, count, R1_count, R2_count):

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

#If no compression is desired on an interleaved file
def no_compress_interleaved(argsOutput, output_1, output_2, input_file_name, count, R1_count, R2_count):
    logging.debug('Writing to file.')

    with open(f'{input_file_name}', 'wt') as IDs_Seq_File:

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

#How to write the interleaved file
def write_interleaved_file(argsRead1, argsOutput, argsCompress, output_1, output_2, input_file_name):

    logging.debug('In writing interleaved file')

    count = 1
    R1_count = 0
    R2_count = 0

    #Compression desired of interleaved file
    if argsRead1.endswith('.gz'):

        compress_interleaved(argsOutput, output_1, output_2, input_file_name, count, R1_count, R2_count)

    #No compression desired of interleaved file
    if argsCompress == True and not argsRead1.endswith('.gz'):
        
        input_file_name = input_file_name + '.gz'

        compress_interleaved(argsOutput, output_1, output_2, input_file_name, count, R1_count, R2_count)

    #No compression desired
    if argsCompress != True and not argsRead1.endswith('.gz'):

        no_compress_interleaved(argsOutput, output_1, output_2, input_file_name, count, R1_count, R2_count)

#Writing output file 
def write_output_file(argsRead1, argsRead2, argsCompress, file_name_1, file_name_2, output_1, output_2, argsOutput):

    #If there is no read 2
    if argsRead2 == None:

        logging.debug('Entered into argsread2 == None in write_output_file')

        #Compression inherent since OG file ends with .gz
        if argsRead1.endswith('.gz'):

            compress_gz(file_name_1, output_1, argsOutput)

        #Compression desired and OG file does not end with .gz
        elif argsCompress == True and not argsRead1.endswith('.gz'):

            file_name_1 = file_name_1 + '.gz'

            compress_gz(file_name_1, output_1, argsOutput)

        #No compression desired
        elif argsCompress != True:

            no_compress(file_name_1, output_1, argsOutput)               

    #If read 2 is present
    if argsRead2 != None:

        logging.debug('Entered into argsread2 != None in write_output_file')

        #Compression is inherent since OG file ends with .gz
        if argsRead1.endswith('.gz'):

            compress_gz(file_name_1, output_1, argsOutput)

        if argsRead2.endswith('.gz'):
            
            compress_gz(file_name_2, output_2, argsOutput)

        #Compression desired and OG files does not end with .gz
        if argsCompress == True and not argsRead1.endswith('.gz'):

            file_name_1 = file_name_1 + '.gz'
            compress_gz(file_name_1, output_1, argsOutput)

        if argsCompress == True and not argsRead2.endswith('.gz'):

            file_name_2 = file_name_2 + '.gz'
            compress_gz(file_name_2, output_2, argsOutput)

        #No compression desired
        if argsCompress != True and not argsRead1.endswith('.gz'):

            no_compress(file_name_1, output_1, argsOutput)

        if argsCompress != True and not argsRead2.endswith('.gz'):

            no_compress(file_name_2, output_2, argsOutput)     

    logging.debug('Finished writing to file.')
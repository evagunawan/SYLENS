#!/usr/bin/env python3
import gzip
import logging
import sys
import re

from Bio import SeqIO

def write_reads(fastq_object, Read1_IDs, Read2_IDs, outputFormat, Compression, outputType, file_name_1, file_name_2):

    def get_seq_from_IDs(input_fastq, id_list):
        
        logging.debug(f'Getting sequences from {input_fastq}')

        seq_records = []
        
        for each_ID in id_list:
            
            each_ID = each_ID.split()[0]
            
            for record in SeqIO.parse(input_fastq, 'fastq'):
                
                if re.search(each_ID, record.id):
                    seq_records.append(record)

        logging.debug(f'Finished getting sequences from {input_fastq}')

        return seq_records

    def write_compressed(file_name, sequences, outputFormat):

        with gzip.open(f'{file_name}', 'wt') as IDs_Seq_File: 

            SeqIO.write(sequences, IDs_Seq_File, outputFormat)
            logging.info(f"Writing to {file_name}")

    def write_not_compressed(file_name, sequences, outputFormat):

        with open(f'{file_name}', 'wt') as IDs_Seq_File:

            SeqIO.write(sequences, IDs_Seq_File, outputFormat)
            logging.info(f"Writing to {file_name}")

    if Read2_IDs == None:
        
        logging.debug("Processing single-end file")
        
        sequences_to_write = get_seq_from_IDs(fastq_object.Read1Temp, Read1_IDs)
        
        if fastq_object.Read1Path.endswith(".gz") or Compression == True:
            
            write_compressed(file_name_1, sequences_to_write, outputFormat)

        else:
            
            write_not_compressed(file_name_1, sequences_to_write, outputFormat)

    else:
        
        logging.debug("Processing interleaved/paired end file(s)")
        
        sequences_to_write_1 = get_seq_from_IDs(fastq_object.Read1Temp, Read1_IDs)
        sequences_to_write_2 = get_seq_from_IDs(fastq_object.Read2Temp, Read2_IDs)
        
        if outputType == "separate":
            
            if fastq_object.Read1Path.endswith(".gz") or Compression == True:
                write_compressed(file_name_1, sequences_to_write_1, outputFormat)
        
            else:
                write_not_compressed(file_name_1, sequences_to_write_1, outputFormat)

            if fastq_object.Read2Path.endswith(".gz") or Compression == True:
                write_compressed(file_name_2, sequences_to_write_2, outputFormat)
            
            else:
                write_not_compressed(file_name_2, sequences_to_write_2, outputFormat)
            
        if outputType == "interleaved":
            
            R1_count = 0
            R2_count = 0
            count = 1

            if fastq_object.Read1Path.endswith(".gz") or Compression == True:
                
                file_name = file_name_1.split(".")[0]
                file_name = file_name + "_interleaved.fastq.gz"

                with gzip.open(f'{file_name}', 'wt') as output:
                    
                    while count < (1 + 2*len(Read1_IDs)):
                    
                        if count % 2 != 0 :
                            count += 1
                            SeqIO.write(sequences_to_write_1[R1_count], output, outputFormat)
                            R1_count += 1
                    
                        else:
                    
                            if count % 2 == 0:
                                count+= 1
                                SeqIO.write(sequences_to_write_2[R2_count], output, outputFormat)
                                R2_count += 1
            else:

                file_name = file_name_1.split(".")[0]
                file_name = file_name + "_interleaved.fastq"

                with open(f'{file_name}', 'wt') as output:

                        while count < (1 + 2*len(Read1_IDs)):

                            if count % 2 != 0 :
                                count += 1
                                SeqIO.write(sequences_to_write_1[R1_count], output, outputFormat)
                                R1_count += 1

                            else:

                                if count % 2 == 0:
                                    count+= 1
                                    SeqIO.write(sequences_to_write_2[R2_count], output, outputFormat)
                                    R2_count += 1


        # Write R1 first and R2 second
        if outputType == "joined":    
            
            if fastq_object.Read1Path.endswith(".gz") or Compression == True:

                file_name = file_name_1.split(".")[0]
                file_name = file_name + "_joined.fastq.gz"

                with gzip.open(f'{file_name}', 'wt') as IDs_Seq_File: 

                    SeqIO.write(sequences_to_write_1, IDs_Seq_File, outputFormat)    
                    SeqIO.write(sequences_to_write_2, IDs_Seq_File, outputFormat)

            else:
                
                file_name = file_name_1.split(".")[0]
                file_name = file_name + "_joined.fastq"

                with open(f'{file_name}', 'wt') as IDs_Seq_File: 

                    SeqIO.write(sequences_to_write_1, IDs_Seq_File, outputFormat)    
                    SeqIO.write(sequences_to_write_2, IDs_Seq_File, outputFormat)
            
            logging.info(f"Writing to {file_name}")

    logging.info("Finished processing, Sylens exiting.")























# #If compression is desired
# def compress_gz(file_name, output_object, argsOutput):
#     with gzip.open(f'{file_name}', 'wt') as IDs_Seq_File:
#         SeqIO.write(output_object, IDs_Seq_File, argsOutput)

# #If no compression is desired
# def no_compress(file_name, object_output, argsOutput):
#     with open(f'{file_name}', 'wt') as IDs_Seq_File:
#         SeqIO.write(object_output, IDs_Seq_File, argsOutput)  

# #If file is an interleaved file that needs to be compressed 
# def write_interleaved(outputFormat, output_1, output_2, input_file_name, count, R1_count, R2_count, interleaved_output):

#     # Writing every other entry to output file
#     if interleaved_output == "every_other":
#         while count < (1 + 2*len(output_1)):
#             if count % 2 != 0 :
#                 count += 1
#                 SeqIO.write(output_1[R1_count], IDs_Seq_File, outputFormat)
#                 R1_count += 1
#             else:
#                 if count % 2 == 0:
#                     count+= 1
#                     SeqIO.write(output_2[R2_count], IDs_Seq_File, outputFormat)
#                     R2_count += 1
    
#     # Write R1 first and R2 second
#     if interleaved_output == "joined":    
#         SeqIO.write(output_1, IDs_Seq_File, outputFormat)    
#         SeqIO.write(output_2, IDs_Seq_File, outputFormat)

#     # Write in separate files
#     if interleaved_output == "separate":    
#         IDs_Seq_File_2 = "f{IDs_Seq_File_1}_2"
#         SeqIO.write(output_1, IDs_Seq_File, outputFormat)    
#         SeqIO.write(output_2, IDs_Seq_File_2, outputFormat)

# #How to write the interleaved file
# def write_interleaved_file(argsRead1, argsOutput, argsCompress, output_1, output_2, input_file_name):

#     logging.debug('In writing interleaved file')

#     count = 1
#     R1_count = 0
#     R2_count = 0

#     #Compression desired of interleaved file
#     if argsRead1.endswith('.gz'):

#         compress_interleaved(argsOutput, output_1, output_2, input_file_name, count, R1_count, R2_count)

#     #No compression desired of interleaved file
#     if argsCompress == True and not argsRead1.endswith('.gz'):
        
#         input_file_name = input_file_name + '.gz'

#         compress_interleaved(argsOutput, output_1, output_2, input_file_name, count, R1_count, R2_count)

#     #No compression desired
#     if argsCompress != True and not argsRead1.endswith('.gz'):

#         no_compress_interleaved(argsOutput, output_1, output_2, input_file_name, count, R1_count, R2_count)

# #Writing output file 
# def write_output_file(argsRead1, argsRead2, argsCompress, file_name_1, file_name_2, output_1, output_2, argsOutput):

#     #If there is no read 2
#     if argsRead2 == None:

#         logging.debug('Entered into argsread2 == None in write_output_file')

#         #Compression inherent since OG file ends with .gz
#         if argsRead1.endswith('.gz'):

#             compress_gz(file_name_1, output_1, argsOutput)

#         #Compression desired and OG file does not end with .gz
#         elif argsCompress == True and not argsRead1.endswith('.gz'):

#             file_name_1 = file_name_1 + '.gz'

#             compress_gz(file_name_1, output_1, argsOutput)

#         #No compression desired
#         elif argsCompress != True:

#             no_compress(file_name_1, output_1, argsOutput)               

#     #If read 2 is present
#     if argsRead2 != None:

#         logging.debug('Entered into argsread2 != None in write_output_file')

#         #Compression is inherent since OG file ends with .gz
#         if argsRead1.endswith('.gz'):

#             compress_gz(file_name_1, output_1, argsOutput)

#         if argsRead2.endswith('.gz'):
            
#             compress_gz(file_name_2, output_2, argsOutput)

#         #Compression desired and OG files does not end with .gz
#         if argsCompress == True and not argsRead1.endswith('.gz'):

#             file_name_1 = file_name_1 + '.gz'
#             compress_gz(file_name_1, output_1, argsOutput)

#         if argsCompress == True and not argsRead2.endswith('.gz'):

#             file_name_2 = file_name_2 + '.gz'
#             compress_gz(file_name_2, output_2, argsOutput)

#         #No compression desired
#         if argsCompress != True and not argsRead1.endswith('.gz'):

#             no_compress(file_name_1, output_1, argsOutput)

#         if argsCompress != True and not argsRead2.endswith('.gz'):

#             no_compress(file_name_2, output_2, argsOutput)     

#     logging.debug('Finished writing to file.')
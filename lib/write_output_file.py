#!/usr/bin/env python3
import gzip
import logging
import os

from Bio import SeqIO

def process_and_write_reads(fastq_object, compression, read_id, file_name, outputFormat, index):

    # If compression is desired or if original file was compressed
        if compression or fastq_object.Read1Path.endswith(".gz"):

            if fastq_object.Read1Path.endswith(".gz"):
                pass
            else:
                file_name += ".gz"

            # If the file path does not exist, create and write to the file
            if not os.path.isfile(file_name):
                with gzip.open(f'{file_name}', 'wt') as IDs_Seq_File:                 
                    SeqIO.write(index[read_id], IDs_Seq_File, outputFormat)

            # If the file path does exist, append to the file
            else:
                with gzip.open(f'{file_name}', 'at') as IDs_Seq_File:
                    IDs_Seq_File.write(str(SeqIO.write(index[read_id], IDs_Seq_File, outputFormat)))

        # If no compression is desired
        else:

            # If the file path does not exist, create and write to the file
            if not os.path.isfile(file_name):
                with open(f'{file_name}', 'wt') as IDs_Seq_File:
                    SeqIO.write(index[read_id], IDs_Seq_File, outputFormat)

            # If the file path does exist, append to the file
            else:
                with open(f'{file_name}', 'at') as IDs_Seq_File:
                    SeqIO.write(index[read_id], IDs_Seq_File, outputFormat)

def write_single(fastq_object, Read1_IDs, outputFormat, compression, file_name_1):

    for each_read_id in Read1_IDs:

        process_and_write_reads(fastq_object, compression, each_read_id, file_name_1, outputFormat, fastq_object.Read1Index)

def process_joined(fastq_object, Read1_IDs, Read2_IDs, outputFormat, compression, file_name_1):

    file_name_1 = "joined_" + file_name_1

    for each_read_id in Read1_IDs:
        process_and_write_reads(fastq_object, compression, each_read_id, file_name_1, outputFormat, fastq_object.Read1Index)

    for each_read_id in Read2_IDs:

        # If Read 2 index doesn't exist (Interleaved)
        if not fastq_object.Read2Index:
            process_and_write_reads(fastq_object, compression, each_read_id, file_name_1, outputFormat, fastq_object.Read1Index)

        # If Read 2 Index exists (Paired)
        else:
            process_and_write_reads(fastq_object, compression, each_read_id, file_name_1, outputFormat, fastq_object.Read2Index)

def process_separate(fastq_object, Read1_IDs, Read2_IDs, outputFormat, compression, file_name_1, file_name_2):
    
    for each_read_id in Read1_IDs:
        process_and_write_reads(fastq_object, compression, each_read_id, file_name_1, outputFormat, fastq_object.Read1Index)

    for each_read_id in Read2_IDs:

        # If Read 2 index doesn't exist (Interleaved)
        if not fastq_object.Read2Index:
            process_and_write_reads(fastq_object, compression, each_read_id, file_name_2, outputFormat, fastq_object.Read1Index)

        # If Read 2 Index exists (Paired)
        else:
            process_and_write_reads(fastq_object, compression, each_read_id, file_name_2, outputFormat, fastq_object.Read2Index)

def process_interleaved(fastq_object, Read1_IDs, Read2_IDs, outputFormat, compression, file_name_1):

    file_name_1 = "interleaved_" + file_name_1

    count = 1
    read1_count = 0
    read2_count = 0

    while count < (1 + 2*len(Read1_IDs)):

        if count % 2 != 0 :

            read1 = Read1_IDs[read1_count]

            # If compression is desired or if original file was compressed
            if compression or fastq_object.Read1Path.endswith(".gz"):

                if not fastq_object.Read1Path.endswith(".gz"):
                    file_name += ".gz"

                # If the file path does not exist, create and write to the file
                if not os.path.isfile(file_name_1):
                    with gzip.open(f'{file_name_1}', 'wt') as IDs_Seq_File: 
                        SeqIO.write(fastq_object.Read1Index[read1], IDs_Seq_File, outputFormat)
                
                # If the file path does exist, append to the file
                else:
                    with gzip.open(f'{file_name_1}', 'at') as IDs_Seq_File:
                        IDs_Seq_File.write(str(SeqIO.write(fastq_object.Read1Index[read1], IDs_Seq_File, outputFormat)))
            
            # If no compression is desired
            else:

                # If the file path does not exist, create and write to the file
                if not os.path.isfile(file_name_1):
                    with open(f'{file_name_1}', 'wt') as IDs_Seq_File:
                        SeqIO.write(fastq_object.Read1Index[read1], IDs_Seq_File, outputFormat)
                
                # If the file path does  exist, append to the file                
                else:
                    with open(f'{file_name_1}', 'at') as IDs_Seq_File:
                        SeqIO.write(fastq_object.Read1Index[read1], IDs_Seq_File, outputFormat)
            count += 1
            read1_count += 1

        else:

            if count % 2 == 0:
                read2 = Read2_IDs[read2_count]

                # If Read 2 Index does not exist (Interleaved)
                if not fastq_object.Read2Index:

                    # If compression is desired or if original file is compressed
                    if compression or fastq_object.Read1Path.endswith(".gz"):

                        # File will already exists, so appending to file
                        with gzip.open(f'{file_name_1}', 'at') as IDs_Seq_File:                 
                            SeqIO.write(fastq_object.Read1Index[read2], IDs_Seq_File, outputFormat)

                    # If no compression is desired
                    else:
                        with open(f'{file_name_1}', 'at') as IDs_Seq_File:
                            SeqIO.write(fastq_object.Read1Index[read2], IDs_Seq_File, outputFormat)

                # If read 2 index does exist (paired end files)
                else:

                    # If compression is desired or if original file is compressed
                    if compression or fastq_object.Read1Path.endswith(".gz"):

                        # File will already exists, so appending to file
                        with gzip.open(f'{file_name_1}', 'at') as IDs_Seq_File:                 
                            SeqIO.write(fastq_object.Read2Index[read2], IDs_Seq_File, outputFormat)

                    # If no compression is desired
                    else:
                        with open(f'{file_name_1}', 'at') as IDs_Seq_File:
                            SeqIO.write(fastq_object.Read2Index[read2], IDs_Seq_File, outputFormat)

                read2_count += 1
                count += 1

def write_interleaved_paired_end(fastq_object, Read1_IDs, Read2_IDs, outputFormat, compression, file_name_1, file_name_2, outputType):

    if outputType == 'separate':
        process_separate(fastq_object, Read1_IDs, Read2_IDs, outputFormat, compression, file_name_1, file_name_2)

    if outputType == 'joined':
        process_joined(fastq_object, Read1_IDs, Read2_IDs, outputFormat, compression, file_name_1)

    if outputType == 'interleaved':
        process_interleaved(fastq_object, Read1_IDs, Read2_IDs, outputFormat, compression, file_name_1)

    logging.info("Finished processing, Sylens exiting.")
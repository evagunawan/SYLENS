#!/usr/bin/env python3
import gzip
import logging

from Bio import SeqIO
import lib.fastq_format_check

def write_reads(fastq_object, Read1_IDs, Read2_IDs, outputFormat, compression, outputType, file_name_1, file_name_2, input_type):

    def get_seq_from_IDs(input_fastq, id_list, input_type):

        logging.debug(f'Getting sequences from {input_fastq}')
        seq_records = {}

        # Using comprehension to extract just the prefix of the id for every id in the id list
        # if lib.fastq_format_check.determine_fastq_ID_formatting == "IlluminaAndCasava":
        #     id_set = set(each_ID.split() for each_ID in id_list)
        # else:
        #     
        id_set = set(each_ID.split()[0] for each_ID in id_list)

# This works better for smaller datasets instead of larger datasets 
# Reframe how to think about this: want to keep one read in memory at a time!
# Write them immediately after gathering them one by one instead of storing all in memory 
# SeqIO.Index function, read1index dictionary like 
# Just use index and read IDs and output file style 
        for record in SeqIO.parse(input_fastq, input_type):
            if record.id in id_set:
                seq_records[record.id] = record

        logging.debug(f'Finished getting sequences from {input_fastq}')
        return list(seq_records.values())

    def write_output(file_name, sequences, outputFormat, compression):

        logging.debug(f"Writing output to {file_name}")

        if compression or fastq_object.Read1Path.endswith(".gz"):
            if fastq_object.Read1Path.endswith(".gz"):
                pass
            else:
                file_name += ".gz"
            with gzip.open(f'{file_name}', 'wt') as IDs_Seq_File: 
                SeqIO.write(sequences, IDs_Seq_File, outputFormat)

        else:
            with open(f'{file_name}', 'wt') as IDs_Seq_File:
                SeqIO.write(sequences, IDs_Seq_File, outputFormat)

    def write_interleaved(sequences_1, sequences_2,output, outputFormat):

        R1_count = 0
        R2_count = 0
        count = 1
        
        while count < (1 + 2*len(Read1_IDs)):
            if count % 2 != 0 :
                count += 1
                SeqIO.write(sequences_1[R1_count], output, outputFormat)
                R1_count += 1

            else:
                if count % 2 == 0:
                    count+= 1
                    SeqIO.write(sequences_2[R2_count], output, outputFormat)
                    R2_count += 1

    if Read2_IDs == None:

        logging.debug("Processing single-end file")
        
        if fastq_object.Read1Temp:
            sequences_to_write = get_seq_from_IDs(fastq_object.Read1Temp, Read1_IDs, input_type)
            write_output(file_name_1, sequences_to_write, outputFormat, compression)
        else:
            sequences_to_write = get_seq_from_IDs(fastq_object.Read1Path, Read1_IDs, input_type)
            write_output(file_name_1, sequences_to_write, outputFormat, compression)

        logging.info(f"Finished writing to output: {file_name_1}")

    else:
        
        logging.debug("Processing interleaved/paired end file(s)")

        # For interleaved processing        
        if fastq_object.Interleaved_Read1_IDs != None and fastq_object.Interleaved_Read2_IDs != None:
            if fastq_object.Read1Temp:
                sequences_to_write_1 = get_seq_from_IDs(fastq_object.Read1Temp, Read1_IDs, input_type)
                sequences_to_write_2 = get_seq_from_IDs(fastq_object.Read1Temp, Read2_IDs, input_type)
            else:
                sequences_to_write_1 = get_seq_from_IDs(fastq_object.Read1Path, Read1_IDs, input_type)
                sequences_to_write_2 = get_seq_from_IDs(fastq_object.Read1Path, Read2_IDs, input_type)
            if outputType == "separate":
                file_name_2 = file_name_1.split(".")[0] + "_Read_2.fastq"
                file_name_1 = file_name_1.split(".")[0] + "_Read_1.fastq"
                if compression or fastq_object.Read1Path.endswith('.gz'):
                    file_name_1 += ".gz"
                    file_name_2 += ".gz"
                    

        # For paired end processing
        if fastq_object.Interleaved_Read1_IDs == None and fastq_object.Interleaved_Read2_IDs == None:
            if fastq_object.Read1Temp:
                sequences_to_write_1 = get_seq_from_IDs(fastq_object.Read1Temp, Read1_IDs, input_type)
            else:
                sequences_to_write_1 = get_seq_from_IDs(fastq_object.Read1Path, Read1_IDs, input_type)

            if fastq_object.Read2Temp:
                sequences_to_write_2 = get_seq_from_IDs(fastq_object.Read2Temp, Read2_IDs, input_type)
            else:
                sequences_to_write_2 = get_seq_from_IDs(fastq_object.Read2Path, Read2_IDs, input_type)
        
        if outputType == "separate":

            write_output(file_name_1, sequences_to_write_1, outputFormat, compression)
            logging.info(f"Finished writing to output: {file_name_1}")
            
            write_output(file_name_2, sequences_to_write_2, outputFormat, compression)           
            logging.info(f"Finished writing to output: {file_name_2}")
            
        if outputType == "interleaved":

            file_name = file_name_1.split(".")[0]
            file_name = file_name + "_interleaved.fastq"

            if fastq_object.Read1Path.endswith(".gz") or compression:
                if fastq_object.Read1Path.endswith(".gz"):
                    file_name += ".gz"

                with gzip.open(f'{file_name}', 'wt') as output:
                    write_interleaved(sequences_to_write_1, sequences_to_write_2, output, outputFormat)

            else:
                with open(f'{file_name}', 'wt') as output:
                    write_interleaved(sequences_to_write_1, sequences_to_write_2, output, outputFormat)

            logging.info(f"Finished writing to output: {file_name}")

        if outputType == "joined":    
            if fastq_object.Read1Path.endswith(".gz") or compression:
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

            logging.info(f"Finished writing to output: {file_name}")

    
    for item in Read1_IDs:
        print(fastq_object.Read1Index[item].id)
        print(fastq_object.Read1Index[item].seq)
        print(fastq_object.Read1Index[item].format('fastq'))
        
    logging.info("Finished processing, Sylens exiting.")
#!/usr/bin/env python3

import re
import sys
import gzip
import logging
from Bio import SeqIO


def pairedEndProcess(dictionary_name_1, dictionary_name_2, subsample_max):
    print('paired end')

sys.exit()
    # if subsample_max != None: 
    #     logging.info("Read 1 and Read 2 supplied. Downsampling starting...")
    #     for element in list(dictionary_name_1):
    #         if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
    #             paired_end_1.append(element)
    #     for element in list(dictionary_name_2): 
    #         if re.search(r"(^\w+.)(\w+)(.)(2$)", element):
    #             paired_end_2.append(element)
    #     zipped_ids_list = list(zip(paired_end_1, paired_end_2))
    #     random_ids = random.sample(zipped_ids_list, subsample_max)
    #     paired_end_1_ids, paired_end_2_ids = zip(*random_ids)
    #     output_1 = [dictionary_name_1[info] for info in paired_end_1_ids]
    #     output_2 = [dictionary_name_2[info] for info in paired_end_2_ids]  
    #     logging.info('Writing to file...')
        
        
    #     if args.Read1[0].endswith('.gz'):
    #         with gzip.open(f"down_sampled_{args.Read1[0]}", 'w') as IDs_Seq_File_1:
    #             SeqIO.write(output_1, IDs_Seq_File_1, args.output)
        

    #     if args.Read2.endswith('.gz'):
    #         with gzip.open(f'down_sample_{args.Read2}', 'w') as IDs_Seq_File_2:
    #             SeqIO.write(output_2, IDs_Seq_File_2, args.output)
        

    #     if args.Read1[0].endswith('.fastq'):
    #         if args.compress == 'no' or args.compress == 'NO' or args.compress == 'No' or args.compress == 'nO':
    #             with open(f"down_sampled_{args.Read1[0]}", 'w') as IDs_Seq_File_1:
    #                 SeqIO.write(output_1, IDs_Seq_File_1, args.output)
    #         if args.compress == 'yes' or args.compress == 'YES' or args.compress == 'Yes' or args.compress == 'yES': 
    #             with gzip.open(f"down_sampled_{args.Read1[0]}", 'w') as IDs_Seq_File_1:
    #                 SeqIO.write(output_1, IDs_Seq_File_1, args.output)
        

    #     if args.Read2.endswith('.fastq'):
    #         if args.compress == 'no' or args.compress == 'NO' or args.compress == 'No' or args.compress == 'nO':
    #             with open(f"down_sampled_{args.Read2}", 'w') as IDs_Seq_File_2:
    #                 SeqIO.write(output_2, IDs_Seq_File_2, args.output)
    #         if args.compress == 'yes' or args.compress == 'YES' or args.compress == 'Yes' or args.compress == 'yES': 
    #             with gzip.open(f"down_sampled_{args.Read2}", 'w') as IDs_Seq_File_2:
    #                 SeqIO.write(output_2, IDs_Seq_File_2, args.output)
        
    #     logging.info(f"Done writing to: down_sampled_{args.Read1[0]} and down_sampled_{args.Read2}") 
    #     logging.debug('SUCCESS: finished downsample_paired_end_files')
    
    
    # if subsample_max == None:
    #     logging.info("Read 1 and Read 2 supplied. No downsampling occurring. Processing Read 1 and Read 2 files...")
    #     for element in list(dictionary_name_1):
    #         if re.search(r"(^\w+.)(\w+)(.)(1$)", element):
    #             paired_end_1.append(element)
    #     for element in list(dictionary_name_2): 
    #         if re.search(r"(^\w+.)(\w+)(.)(2$)", element):
    #             paired_end_2.append(element)
    #     output_1 = [dictionary_name_1[info] for info in paired_end_1]
    #     output_2 = [dictionary_name_2[info] for info in paired_end_2]  
    #     logging.info('Writing to file...')

        
    #     if args.Read1[0].endswith('.gz'):
    #         with gzip.open(f"down_sampled_{args.Read1[0]}", 'w') as IDs_Seq_File_1:
    #             SeqIO.write(output_1, IDs_Seq_File_1, args.output)
        

    #     if args.Read2.endswith('.gz'):
    #         with gzip.open(f'down_sample_{args.Read2}', 'w') as IDs_Seq_File_2:
    #             SeqIO.write(output_2, IDs_Seq_File_2, args.output)
        

    #     if args.Read1[0].endswith('.fastq'):
    #         if args.compress == 'no' or args.compress == 'NO' or args.compress == 'No' or args.compress == 'nO':
    #             with open(f"down_sampled_{args.Read1[0]}", 'w') as IDs_Seq_File_1:
    #                 SeqIO.write(output_1, IDs_Seq_File_1, args.output)
    #         if args.compress == 'yes' or args.compress == 'YES' or args.compress == 'Yes' or args.compress == 'yES': 
    #             with gzip.open(f"down_sampled_{args.Read1[0]}", 'w') as IDs_Seq_File_1:
    #                 SeqIO.write(output_1, IDs_Seq_File_1, args.output)
        

    #     if args.Read2.endswith('.fastq'):
    #         if args.compress == 'no' or args.compress == 'NO' or args.compress == 'No' or args.compress == 'nO':
    #             with open(f"down_sampled_{args.Read2}", 'w') as IDs_Seq_File_2:
    #                 SeqIO.write(output_2, IDs_Seq_File_2, args.output)
    #         if args.compress == 'yes' or args.compress == 'YES' or args.compress == 'Yes' or args.compress == 'yES': 
    #             with gzip.open(f"down_sampled_{args.Read2}", 'w') as IDs_Seq_File_2:
    #                 SeqIO.write(output_2, IDs_Seq_File_2, args.output)


    #     logging.info("Done writing to: non_down_sampled_R1 and non_down_sampled_R2")
    #     logging.debug('SUCCESS: finished no_downsample_paired_end_files')

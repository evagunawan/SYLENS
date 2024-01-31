import re
import logging

#Taking the top 3 known formats and creating regular expressions to identify IDs. Stored REs in a dictinoary 
# Illumina and Casava have same format, for now.

formatList = (
    ("NCBI", '(^SRR)(\w+)[.+](\d+)[.+](1)'),
    ('IlluminaAndCasava', '(.+)(\d) (1)'),
    ("NCBI", '(^SRR)(\w+)[.+](\d+)[.+](2)'),
    ('IlluminaAndCasava', '(.+)(\d) (2)')
)

#Determining format of fastq file to properly figure out if R1 and/or R2
def determine_fastq_ID_formatting(ReadIndex):

    #then get first and last ids
    first_ID = list(ReadIndex.keys())[0]
    last_ID = list(ReadIndex.keys())[-1]

    first_ID_Format = None
    last_ID_Format = None

    #loop over format expressions and save format type if matches
    for format, formatExpression in formatList:

        if re.search(formatExpression, first_ID):
            first_ID_Format = format
        
        if re.search(formatExpression, last_ID):
            last_ID_Format = format
    
    #check we have a matching format and return it
    if first_ID_Format is not None and first_ID_Format == last_ID_Format:
        return first_ID_Format
    
    elif first_ID_Format != last_ID_Format:
        logging.error("The ID of the first read in the fastq file does not match the same format as the ID of the last read in the fastq file.")
    
    else:
        logging.error("The Format of the fastq file did not match a known type. Please contact the developers if a new format type needs to be added.")




#Determining if the fastq file is interleaved
#   Note interleaved could mean R1 followed by R2 or all R1 then all R2.
def get_interleaved_ids(ReadIndex):

    #get first and second and last ids
    first_ID = list(ReadIndex.keys())[0]
    second_ID = list(ReadIndex.keys())[1]
    last_ID = list(ReadIndex.keys())[-1]
    
    first_ID_formatExpression = None
    second_ID_formatExpression = None
    last_ID_formatExpression = None

    detected_format = None

    #loop over formatExpressions and save if one matches
    for format, formatExpression in formatList:

        if re.search(formatExpression, first_ID):
            first_ID_formatExpression = formatExpression
            detected_format = format
        
        if re.search(formatExpression, second_ID):
            second_ID_formatExpression = formatExpression
        
        if re.search(formatExpression, last_ID):
            last_ID_formatExpression = formatExpression

    #if all format expressions are the same then the file is not interleaved
    if first_ID_formatExpression == second_ID_formatExpression and second_ID_formatExpression == last_ID_formatExpression:
        return None, None

    logging.info('Getting read IDs from interleaved file.')
    
    #get list of all IDs that are Read 1 and Read 2
    Read1_IDs = []
    Read2_IDs = []

    for id in ReadIndex.keys():
        for format, formatExpression in formatList:
            if format == detected_format:
                if re.search(formatExpression, id):
                    if "1" in formatExpression:
                        Read1_IDs.append(id)
                    if "2" in formatExpression:
                        Read2_IDs.append(id)

    if len(Read1_IDs) == len(Read2_IDs):
        return Read1_IDs,Read2_IDs
    
    else:
        logging.warn("There were a differing number of R1 to R2 reads in the interleaved fastq.")

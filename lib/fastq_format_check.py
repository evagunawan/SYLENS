import re
import chardet
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

    #These are messy but its the only way to get the pure fastq ids
    #   we have to use the index to get the raw read, then get the first line of the read and convert to str

    #first we get char encoding
    encoding = chardet.detect(ReadIndex.get_raw(list(ReadIndex.keys())[0]).split(b'\n')[0])

    #then get raw first and last ids and decode them
    first_ID = ReadIndex.get_raw(list(ReadIndex.keys())[0]).split(b'\n')[0].decode(encoding['encoding'])
    last_ID = ReadIndex.get_raw(list(ReadIndex.keys())[-1]).split(b'\n')[0].decode(encoding['encoding'])
    
    first_ID_Format = None
    last_ID_Format = None

    for format, formatExpression in formatList:

        if re.search(formatExpression, first_ID):
            first_ID_Format = format
        
        if re.search(formatExpression, last_ID):
            last_ID_Format = format
    
    if first_ID_Format is not None and first_ID_Format == last_ID_Format:
        return first_ID_Format
    
    elif first_ID_Format != last_ID_Format:
        raise Exception("The ID of the first read in the fastq file does not match the same format as the ID of the last read in the fastq file.")
    
    else:
        raise Exception("The Format of the fastq file did not match a known type. Please contact the developers if a new format type needs to be added.")
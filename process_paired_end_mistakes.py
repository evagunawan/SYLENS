#! usr/bin/env/ python3
import logging
import re
import sys

def process_mistakes(firstR1, lastR1, firstR2, lastR2, formatExpression, formatExpression2):

    logging.info('Checking for common errors.')

    #Errors if wrong files are supplied and terminates process
    if re.search(formatExpression2, firstR1) and re.search(formatExpression2, lastR1):

        if re.search(formatExpression, firstR2) and re.search(formatExpression, lastR2):

            logging.critical('Supplied files appear to be switched with Read 2 first and Read 1 second. Please flip files and restart Sylens. Program terminating...')

            sys.exit(1)

    if re.search(formatExpression, firstR1) and re.search(formatExpression2, lastR1):

        if re.search(formatExpression, firstR2) and re.search(formatExpression2,lastR2):

            logging.critical('Both supplied files appear to be interleaved. Program terminating...')

            sys.exit(1)

    if re.search(formatExpression,firstR1) and re.search(formatExpression, lastR1):

        if re.search(formatExpression, firstR2) and re.search(formatExpression, lastR2):

            logging.critical("Both supplied files appear to be Read 1 files. Program terminating...")

            sys.exit(1)
    
    if re.search(formatExpression2, firstR1) and re.search(formatExpression2, lastR1):

        if re.search(formatExpression2, firstR2) and re.search(formatExpression2, lastR2):

            logging.critical("Both supplied files appear to be Read 2 files. Program terminating...")

            sys.exit(1)

    if re.search(formatExpression, firstR1) and re.search(formatExpression2, lastR1):

        logging.critical("First input file is interleaved. Program terminating...")

        sys.exit(1)    

    if re.search(formatExpression, firstR2) and re.search(formatExpression2, lastR2):

        logging.critical("Second input file is interleaved. Program terminating...")

        sys.exit(1)
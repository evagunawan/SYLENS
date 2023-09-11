#! usr/bin/env/ python3
import logging
import sys

def process_mistakes(firstR1, lastR1, firstR2, lastR2):

    logging.info('Checking for common errors.')

    #Errors if wrong files are supplied

    if firstR1.endswith('2') and lastR1.endswith('2'):

        if firstR2.endswith('1') and lastR2.endswith('1'):

            logging.critical('Supplied files appear to be switched with Read 2 first and Read 1 second. Please flip files and restart Sylens. Program terminating...')
            sys.exit(1)

    if firstR1.endswith('1') and lastR1.endswith('2'):

        if firstR2.endswith('1') and lastR2.endswith('2'):

            logging.critical('Both supplied files appear to be interleaved. Program terminating...')
            sys.exit(1)

    if firstR1.endswith('1') and lastR1.endswith('1'):

        if firstR2.endswith('1') and lastR2.endswith('1'):

            logging.critical("Both supplied files appear to be Read 1 files. Program terminating...")
            sys.exit(1)
    
    if firstR1.endswith('2') and lastR1.endswith('2'):

        if firstR2.endswith('2') and lastR2.endswith('2'):

            logging.critical("Both supplied files appear to be Read 2 files. Program terminating...")
            sys.exit(1)

    if firstR1.endswith('1') and lastR1.endswith('2'):

        logging.critical("First input file is interleaved. Program terminating...")
        sys.exit(1)    

    if firstR2.endswith('1') and lastR2.endswith('2'):

        logging.critical("Second input file is interleaved. Program terminating...")
        sys.exit(1)
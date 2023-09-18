# SYLENS
![Sylens Logo](Assets/Sylens_logo.png)

**S**ampling **Y**ielding **LE**ss **N**oticeable **S**amples

## **SUMMARY**

Sylens is a Python program designed to intake user inputs through [**argparse**](https://docs.python.org/3/library/argparse.html) and alter the output with [**Bio.SeqIO**](https://biopython.org/wiki/SeqIO). This program not only allows for file output conversions based on user preferences, but can randomly down sample fastq files that have millions of reads for future analyses. 

-------------------------------------------------------------------------------------------------------------------------------------------------------------
## **TABLE OF CONTENTS:**

Summary

Usage

Program Outline

Output

Authors

-------------------------------------------------------------------------------------------------------------------------------------------------------------

## **USAGE**

Sylens requiers **Python 3.8.10** or greater to use. Decompressed and compressed paired-end, single-end, and interleaved fastq files are all file types that can be analyzed with Sylens. 

To begin Sylens with a paired-end file use:
```
Sylens_main.py FILE1.fastq FILE2.fastq
```

To begin with a single-end or interleaved file, use:
```
Sylens_main.py INTERLEAVED.fastq
```

File input type by default is fastq. However, if the input file format is not fastq, use the flag `-f` or `--filetype` with the input file's correct formatting. Currently Sylens supports fastq-sanger, also known as fastq format (ASCII 33), and fastq-solexa (ASCII 64).
```
Sylens_main.py FILE1.fastq -f fastq-solexa
```

Subsampling with Sylens is done through the `-s` or `--subsample` flag with the desired integer to down sample to.
```
Sylens_main.py FILE1.fastq -s 1000
```

By default Sylens will subsample to the exact integer indicated after the -s flag. If a percentage is preferred, use the flag -p or --percentage with the subsample integer to take a percentage of the samples.
```
Sylens_main.py FILE1.fastq -p -s 10
```

Compressing or not compressing a file on output is done by using the `-c` flag. If a .gz file is input, the output will be .gz. By default, no compression occurs on output.
```
Sylens_main.py FILE1.fastq -c
```

By default, files output by Sylens are in fastq format. Changing output file formats is done by adding the `-o` or `--output` flag with the output file type you would like to convert to.
```
Sylens_main.py FILE1.fastq -o fastq-solexa
```

For reproducibility, Sylens provides a seed number. To denote a seed generated from a previous run, use the `--seed` flag with the seed number.
```
Sylens_main.py FILE1.fastq --seed 1691696502
```

If any additional explanations are needed, use the `-h` or `--help` flag.
```
Sylens_main.py FILE1.fastq --help
```

Multiple flags can be utilized in one line of code, if desired.
```
Sylens_main.py FILE1.fastq FILE2.fastq -s 1000 -c --seed 1691696502 -f fastq-solexa -o fastq
```

-------------------------------------------------------------------------------------------------------------------------------------------------------------
## **PROGRAM OUTLINE**
![Sylens Program Map](Assets/Program_Map_Sylens.PNG)

### **Legend**
![Sylens Legend](Assets/legend_Sylens.PNG)
-------------------------------------------------------------------------------------------------------------------------------------------------------------

## **OUTPUT**
Output files by default will be fastq files. If the output filetype indicated is differet than the input format, Bio.SeqIO will write it to the desired output.

## **AUTHORS**
[Eva Gunawan](https://github.com/evagunawan), Bioinformatics Fellow through APHL

[Kelsey Florek](https://github.com/k-florek), WSLH Senior Genomics and Data Scientist

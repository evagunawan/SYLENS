# SYLENS
![text-1692889436721](https://github.com/evagunawan/SYLENS/assets/124393795/4ae4bf26-1471-40b5-b53d-02358534bbb1)

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

Sylens requiers **Python 3.8.10** or greater to use.

This program takes in 'fastq-sanger' or 'fastq-solexa' files. To begin the pipeline with a paired-end file use:
```
Sylens.py FILE1.fastq FILE2.fastq
```

To begin the pipeline with a single-end or interleaved file, use:
```
Sylens.py INTERLEAVED.fastq
```

Subsampling with Sylens is done through the `-s` or `--subsample` flag with the integer you want to down sample to.
```
Sylens.py FILE1.fastq -s 1000
```

Compressing or not compressing a file on output is done by using the `-c` or `--compress`flag with yes or no attached. If a .gz file is input, the output will be .gz.
```
Sylens.py FILE1.fastq -c yes
```

By default, files output by Sylens are in fastq format. Changing file formats is done by adding the `-o` or `--output` flag with the output file type you would like to convert to.
```
Sylens.py FILE1.fastq -o fastq-solexa
```

File input type by default is fastq. If the input file format is not fastq, use the flag `-f` or `--filetype` with the input file's correct formatting.
```
Sylens.py FILE1.fastq -f fastq-solexa
```

For reproducibility, Sylens provides a seed number. To denote a seed generated from a previous run use the `--seed` flag with the seed number.
```
Sylens.py FILE1.fastq -s 1691696502
```

If any additional explanations are needed, use the `-h` or `--help` flag.
```
Sylens.py FILE1.fastq --help
```

Multiple flags can be utilized in one line of code, if desired.
```
Sylens.py FILE1.fastq FILE2.fastq -s 1000 -c yes --seed 1691696502 -f fastq-solexa -o fastq-solexa
```

-------------------------------------------------------------------------------------------------------------------------------------------------------------
## **PROGRAM OUTLINE**

![Program_Map](https://github.com/evagunawan/SYLENS/assets/124393795/4f3995f0-3d2f-404f-9d03-9cb80d6defff)

### **Legend**

![legend_Sylens](https://github.com/evagunawan/SYLENS/assets/124393795/d8fb7491-d199-4303-b820-64095f70ac51)


-------------------------------------------------------------------------------------------------------------------------------------------------------------

## **OUTPUT**
Output files by default will be fastq files. If the output filetype indicated is differet than the input format, Bio.SeqIO will write it to the desired output. 

## **AUTHORS**
[Eva Gunawan](https://github.com/evagunawan), Bioinformatics Fellow through APHL

[Kelsey Florek](https://github.com/k-florek), WSLH Senior Genomics and Data Scientist

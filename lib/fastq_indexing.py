from io import StringIO
from typing import (IO,Iterator,Union,AnyStr)
from os import PathLike

from Bio.File import _open_for_random_access
from Bio.File import _IndexedSeqFileProxy
from Bio.File import _IndexedSeqFileDict
from Bio.Seq import Seq

from Bio.SeqIO.QualityIO import FastqPhredIterator, FastqGeneralIterator
from Bio.SeqRecord import SeqRecord


_IOSource = Union[IO[AnyStr], PathLike, str, bytes]
_TextIOSource = _IOSource[str]



SANGER_SCORE_OFFSET = 33
SOLEXA_SCORE_OFFSET = 64



####
#### Custom Changes to Biopython SeqIO Iterators to keep full ID
####

class CustomFastqPhredIterator(FastqPhredIterator):

    def iterate(self, handle: IO[str]) -> Iterator[SeqRecord]:
        """Parse the file and generate SeqRecord objects."""
        assert SANGER_SCORE_OFFSET == ord("!")
        # Originally, I used a list expression for each record:
        #
        # qualities = [ord(letter)-SANGER_SCORE_OFFSET for letter in quality_string]
        #
        # Precomputing is faster, perhaps partly by avoiding the subtractions.
        q_mapping = {
            chr(letter): letter - SANGER_SCORE_OFFSET
            for letter in range(SANGER_SCORE_OFFSET, 94 + SANGER_SCORE_OFFSET)
        }

        for title_line, seq_string, quality_string in FastqGeneralIterator(handle):
            descr = title_line
            id = descr
            name = id
            record = SeqRecord(Seq(seq_string), id=id, name=name, description=descr)
            try:
                qualities = [q_mapping[letter2] for letter2 in quality_string]
            except KeyError:
                raise ValueError("Invalid character in quality string") from None
            # For speed, will now use a dirty trick to speed up assigning the
            # qualities. We do this to bypass the length check imposed by the
            # per-letter-annotations restricted dict (as this has already been
            # checked by FastqGeneralIterator). This is equivalent to:
            # record.letter_annotations["phred_quality"] = qualities
            dict.__setitem__(record._per_letter_annotations, "phred_quality", qualities)
            yield record


def CustomFastqSolexaIterator(source: _TextIOSource, alphabet: None = None,) -> Iterator[SeqRecord]:
    
    if alphabet is not None:
        raise ValueError("The alphabet argument is no longer supported")

    q_mapping = {
        chr(letter): letter - SOLEXA_SCORE_OFFSET
        for letter in range(SOLEXA_SCORE_OFFSET - 5, 63 + SOLEXA_SCORE_OFFSET)
    }

    for title_line, seq_string, quality_string in FastqGeneralIterator(source):
        descr = title_line
        id = descr
        name = id
        record = SeqRecord(Seq(seq_string), id=id, name=name, description=descr)
        try:
            qualities = [q_mapping[letter2] for letter2 in quality_string]
        # DO NOT convert these into PHRED qualities automatically!
        except KeyError:
            raise ValueError("Invalid character in quality string") from None
        # Dirty trick to speed up this line:
        # record.letter_annotations["solexa_quality"] = qualities
        dict.__setitem__(record._per_letter_annotations, "solexa_quality", qualities)
        yield record


def CustomFastqIlluminaIterator(source: _TextIOSource, alphabet: None = None,) -> Iterator[SeqRecord]:
    
    if alphabet is not None:
        raise ValueError("The alphabet argument is no longer supported")

    q_mapping = {
        chr(letter): letter - SOLEXA_SCORE_OFFSET
        for letter in range(SOLEXA_SCORE_OFFSET, 63 + SOLEXA_SCORE_OFFSET)
    }

    for title_line, seq_string, quality_string in FastqGeneralIterator(source):
        descr = title_line
        id = descr
        name = id
        record = SeqRecord(Seq(seq_string), id=id, name=name, description=descr)
        try:
            qualities = [q_mapping[letter2] for letter2 in quality_string]
        except KeyError:
            raise ValueError("Invalid character in quality string") from None
        # Dirty trick to speed up this line:
        # record.letter_annotations["phred_quality"] = qualities
        dict.__setitem__(record._per_letter_annotations, "phred_quality", qualities)
        yield record




####
#### Dictionary of Custom Functions/Classes
####
CustomFormatToIterator = {
    "fastq": CustomFastqPhredIterator,
    "fastq-sanger": CustomFastqPhredIterator,
    "fastq-solexa": CustomFastqSolexaIterator,
    "fastq-illumina": CustomFastqIlluminaIterator,
}




####
#### Custom Changes Seq file randomized access to keep full ID
####
class SeqFileRandomAccess(_IndexedSeqFileProxy):
    def __init__(self, filename, format):
        self._handle = _open_for_random_access(filename)
        self._format = format
        # __getitem__ call:
        self._iterator = CustomFormatToIterator[format]

    def get(self, offset):
        """Return SeqRecord."""
        # Should be overridden for binary file formats etc:
        return next(self._iterator(StringIO(self.get_raw(offset).decode())))




class FastqRandomAccess(SeqFileRandomAccess):
    """Random access to a FASTQ file (any supported variant).

    With FASTQ the records all start with a "@" line, but so can quality lines.
    Note this will cope with line-wrapped FASTQ files.
    """

    def __iter__(self):
        """Iterate over the sequence records in the file."""
        handle = self._handle
        handle.seek(0)
        id = None
        start_offset = handle.tell()
        line = handle.readline()
        if not line:
            # Empty file!
            return
        if line[0:1] != b"@":
            raise ValueError(f"Problem with FASTQ @ line:\n{line!r}")
        while line:
            # assert line[0]=="@"
            # This record seems OK (so far)
            id = line[1:].rstrip().split(b'\n', 1)[0]
            # Find the seq line(s)
            seq_len = 0
            length = len(line)
            while line:
                line = handle.readline()
                length += len(line)
                if line.startswith(b"+"):
                    break
                seq_len += len(line.strip())
            if not line:
                raise ValueError("Premature end of file in seq section")
            # assert line[0]=="+"
            # Find the qual line(s)
            qual_len = 0
            while line:
                if seq_len == qual_len:
                    if seq_len == 0:
                        # Special case, quality line should be just "\n"
                        line = handle.readline()
                        if line.strip():
                            raise ValueError(
                                f"Expected blank quality line, not {line!r}"
                            )
                        length += len(line)  # Need to include the blank ling
                    # Should be end of record...
                    end_offset = handle.tell()
                    line = handle.readline()
                    if line and line[0:1] != b"@":
                        raise ValueError(f"Problem with line {line!r}")
                    break
                else:
                    line = handle.readline()
                    qual_len += len(line.strip())
                    length += len(line)
            if seq_len != qual_len:
                raise ValueError("Problem with quality section")
            yield id.decode(), start_offset, length
            start_offset = end_offset

    def get_raw(self, offset):
        """Return the raw record from the file as a bytes string."""
        handle = self._handle
        handle.seek(offset)
        line = handle.readline()
        data = line
        if line[0:1] != b"@":
            raise ValueError(f"Problem with FASTQ @ line:\n{line!r}")
        # Find the seq line(s)
        seq_len = 0
        while line:
            line = handle.readline()
            data += line
            if line.startswith(b"+"):
                break
            seq_len += len(line.strip())
        if not line:
            raise ValueError("Premature end of file in seq section")
        assert line[0:1] == b"+"
        # Find the qual line(s)
        qual_len = 0
        while line:
            if seq_len == qual_len:
                if seq_len == 0:
                    # Special case, quality line should be just "\n"
                    line = handle.readline()
                    if line.strip():
                        raise ValueError(f"Expected blank quality line, not {line!r}")
                    data += line
                # Should be end of record...
                line = handle.readline()
                if line and line[0:1] != b"@":
                    raise ValueError(f"Problem with line {line!r}")
                break
            else:
                line = handle.readline()
                data += line
                qual_len += len(line.strip())
        if seq_len != qual_len:
            raise ValueError("Problem with quality section")
        return data



def index(filename, format, key_function=None):

    return _IndexedSeqFileDict(FastqRandomAccess(filename,format), key_function, "index(%r, %r)" % (filename, format), "SeqRecord")
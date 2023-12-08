import re

from io import BytesIO
from io import StringIO

from Bio.File import _open_for_random_access
from Bio.File import _IndexedSeqFileProxy
from Bio.File import _IndexedSeqFileDict

class SeqFileRandomAccess(_IndexedSeqFileProxy):
    def __init__(self, filename, format):
        self._handle = _open_for_random_access(filename)
        self._format = format
        # __getitem__ call:
        self._iterator = FastqRandomAccess

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

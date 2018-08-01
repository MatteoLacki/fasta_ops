from    Bio import SeqIO
import  re
from    Py.filter_proteins import filter_proteins

filepath = "data/human_reviewed.fasta"

def iter_proteins(filepath):
    """Iterate over proteins."""
    with open(filepath, "r") as f:
        for r in SeqIO.parse(f, "fasta"):
            yield r

it      = iter_proteins(filepath)
r       = next(it)
seq     = str(r.seq)
pattern = re.compile("R")

list((m.start(), m.end())for m in pattern.finditer(seq))

def iter_findings(seq, pattern, left_offset=10, right_offset=10):
    for m in pattern.finditer(seq):
        s, e = m.start(), m.end()
        yield (s + 1,
               seq[max(s-left_offset, 0):s],
               seq[s:min(s+10, len(seq))])

list(iter_findings(seq, pattern))



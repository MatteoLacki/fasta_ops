import re
from multiprocessing import Pool
from itertools      import product

from fasta_ops.dump_2_csv import write2csv
from fasta_ops.fasta_ops  import iter_proteins


def iter_findings(seq, pattern, left_offset=10, right_offset=10):
    """Iterate over findings.

    Args:
        seq (Bio.SeqRecord.SeqRecord): A biopython record sequence with protein info.
        patter (_sre.SRE_Pattern):     A compiled re pattern to search for.
        left_offset (int):             How much AA to include before the found sequence.
        right_offset (int):            How much AA to include together with the sequence, starting from the first AA in the sequence. 
    Yields:
        tuple: the found sequence, its position in the fasta, left and right sequences.
    """
    for m in pattern.finditer(seq):
        s, e = m.start(), m.end()
        yield (seq[s:e],
               s + 1,
               seq[max(s-left_offset, 0):s],
               seq[s:min(s+10, len(seq))],
               seq)

def iter_findings_across_sequences(filepath,
                                   searched_sequences,
                                   left_offset  = 10,
                                   right_offset = 10):
    p = re.compile("|".join(searched_sequences))
    for prot in iter_proteins(filepath):
        for f in iter_findings(seq          = str(prot.seq),
                               pattern      = p,
                               left_offset  = left_offset,
                               right_offset = right_offset):
            yield f


def iter_solutions(proteins, peptides, verbose = False, iter_ring = 1000):
    K = len(proteins) * len(peptides)
    i = 0
    for prot, pep in product(proteins, peptides):
        for f in iter_findings(str(prot.seq), re.compile(pep)):
            yield f
        if verbose:
            div, mod = divmod(i, iter_ring)
            if mod == 0:
                www = (i, K, 100*(i/K))
                print("Got {} out of {} = {:.2f}%.".format(*www))
        i += 1

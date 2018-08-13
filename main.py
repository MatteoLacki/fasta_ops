import os
from   os.path import join
import re
from   time.time import time

from fasta_ops.fasta_ops import iter_proteins
from fasta_ops.misc import iter_findings,\
                           iter_findings_across_sequences,\
                           iter_solutions
from fasta_ops.dump_2_csv import write2csv
from fasta_ops.array_ops  import percentile_pairs_of_N_integers

# filepath = "data/mus_musculus.fasta"
# searched_sequences = [line.rstrip('\n') for line in open('data/peptides.txt')]
# peptides = [s for s in searched_sequences if len(s) > 0]
# proteins = list(iter_proteins(filepath))
# problems = iter_solutions(proteins, peptides, verbose=True)
# solution = list(problems)
# write2csv("res/mus_musculus.csv", 
#           ['seq',
#            'index',
#            'Ten_AAs_before',
#            'AAs_after_and_including_seq',
#            'sequence_we_investigate'],
          # solution.__iter__())

filepath  = "data/mus_musculus.fasta"
seq_files = [p for p in os.listdir('data') if '.seq' in p]
searched  = [[line.rstrip('\n') for line in open(join('data', f)) if line[0] != '>'] for f in seq_files]
peptides  = [[s for s in ss if len(s) > 0] for ss in searched]
proteins  = list(iter_proteins(filepath))

def foo(searched_against, searched_for,
        left_offset=10, right_offset=10, verbose = False):
    searched_for = "|".join(searched_for)
    RE = re.compile(searched_for)
    N  = len(searched_against)
    for k, p in enumerate(searched_against):
        seq = str(p.seq)
        for f in RE.finditer(seq):
            s = f.start()
            e = f.end()
            S = max(s-left_offset, 0)
            E = min(s+10, len(seq))
            yield f.group(), s+1, seq[S:s], seq[s:E], p.id
        if verbose:
            div, mod = divmod(k, 100)
            if mod == 0:
                print(f"{k}/{N}")

def foo_split(searched_against, searched_for, K):
    N = len(searched_for)
    for s, e in percentile_pairs_of_N_integers(N, K):
        yield searched_against, searched_for[s:e]

def get_res(proteins, peptides, K=20, verbose=True):
    res = []
    for pep_no, p in enumerate(peptides):
        p_res = []
        for i, (sa, sf) in enumerate(foo_split(proteins, p, K)):
            for r in foo(sa, sf, verbose=verbose):
                p_res.append(r)
            if verbose:
                print(f"{i}/{K} for peptide {pep_no}/{len(peptides)}")
        res.append(p_res)
    return res

res2 = get_res(proteins, peptides, K=20, verbose=True)

colnames = ['seq',
            'index',
            'ten_AAs_before',
            'AAs_after_and_including_seq',
            'accession']

for f, r in zip(seq_files, res):
    f = join('res', f.replace('.seq', '.csv'))
    write2csv(f, colnames, r)

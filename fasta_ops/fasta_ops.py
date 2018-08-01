from   Bio import SeqIO



def iter_proteins(filepath):
    """Iterate over proteins."""
    with open(filepath, "r") as f:
        for r in SeqIO.parse(f, "fasta"):
            yield r

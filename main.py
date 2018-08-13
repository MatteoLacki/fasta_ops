from fasta_ops.misc import iter_findings,\
                          iter_findings_across_sequences,\
                          iter_solutions,\
                          write2csv

filepath = "data/mus_musculus.fasta"
searched_sequences = [line.rstrip('\n') for line in open('data/peptides.txt')]
peptides = [s for s in searched_sequences if len(s) > 0]
proteins = list(iter_proteins(filepath))



problems = iter_solutions(proteins, peptides, verbose=True)
solution = list(problems)
write2csv("res/mus_musculus.csv", 
          ['seq',
           'index',
           'Ten_AAs_before',
           'AAs_after_and_including_seq',
           'sequence_we_investigate'],
          solution.__iter__())



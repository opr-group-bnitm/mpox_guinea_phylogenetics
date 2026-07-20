#!/usr/bin/env python

"""
Prepare the alignment to be used in BEAST.
Remove the outgroup from the multiple sequence alignment and add the annotations from a metadata file.
"""

import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main():
    outdir = '../results/IIa/beast_in'
    fname_alignment_in = '../results/IIa/squirrel/sequences_all.aln.fasta'
    fname_annotations = '../data/IIa/meta_all.tsv'
    outgroup = 'NC_063383.1'
    os.makedirs(outdir, exist_ok=True)
    fname_alignment_out = os.path.join(outdir, 'cladeIIa.annotated.aln.fasta')
    meta = pd.read_csv(fname_annotations, sep='\t').set_index('Accession')

    print(meta)

    with open(fname_alignment_out, 'w') as f_out:
        for seq in SeqIO.parse(fname_alignment_in, 'fasta'):
            if seq.id != outgroup:
                row = meta.loc[seq.id]
                seq_id = seq.id.replace('_', '')
                new_id = f'{seq_id}_{row.Date}'
                annotated_seq = SeqRecord(seq.seq, id=new_id, name=new_id, description='')
                SeqIO.write(annotated_seq, f_out, 'fasta')


if  __name__ == '__main__':
    main()

#!/bin/python
path = "/disk/yt/kaks/fa/blast/"
species = ["GRCh38", "GRCm39"]

for s in species:
    blast_file = path + s + ".txt"
    gff_output = path + s + '_20aa.gff'

    with open(blast_file, 'r') as blast_res, open(gff_output, 'w') as out:
        for line in blast_res:
            line_data = line.strip().split('\t')

            pident = line_data[0]
            chromosome = line_data[1]
            feature_id = line_data[2]
            start = line_data[3]
            end = line_data[4]

            if '100' in pident:
                out.write(f"{chromosome}\t{feature_id}\t{start}\t{end}\n")
from Bio import GenBank

import glob
import os
import sys

i = 0
out_path = input('Print output path ')
out_file = open(out_path, 'w')
for file in glob.glob('*gbk'):
    if file == 'ncbi_bact_chr_repr.gbk':
        continue
    i = i + 1
    print(i)
    with open(file) as handle:
        for record in GenBank.parse(handle):
        #record = GenBank.read(handle)
            genome_id = ''.join(record.version)
            print(genome_id)
            coord1 = int(record.structured_comment['antiSMASH-Data']['Orig. start'])
            coord2 = int(record.structured_comment['antiSMASH-Data']['Orig. end'])
            #print(coord1)
            #print(coord2)
            for feature in record.features:
                if feature.key == 'protocluster':
                    intern_coord = feature.location
                    c1,c2=str(intern_coord).split('..')
                    #print(c1,c2)
                    real_coord1 = coord1 + int(c1)
                    real_coord2 = coord1 + int(c2)
                    #print(feature.qualifiers)
                    for q in feature.qualifiers:
                        if q.key == '/product=':
                            product = q.value
                            #print(product)
                            out_file.write(f'{genome_id}\t{real_coord1}\t{real_coord2}\t{product}\n')
    #print(x)
out_file.close()

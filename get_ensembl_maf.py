#!/usr/bin env python3

'''
Timothy Isonio
July 2016
'''

import requests, sys, csv

def chunks(l, n):
    """Yield successive n-sized chunks from l.
    http://stackoverflow.com/a/312464/5287504"""
    for i in range(0, len(l), n):
        yield l[i:i+n]

MAX_RSIDS = 100
rs_ids = []
mafs = []

# get list of rs_ids from input file
with open(sys.argv[1], 'r') as tsv_in:
    tsv_contents = csv.reader(tsv_in, delimiter='\t')
    i = 0
    for row in tsv_contents:
        rs_ids.append(row[0])
        i += 1
        if i > MAX_RSIDS:
            break

# break rs_ids into chunks under 1000 per Ensembl guidelines
rs_ids_chunks = list(chunks(rs_ids, 10))
num_chunks = len(rs_ids_chunks)


# Ensembl parameters
server = "http://grch37.rest.ensembl.org"
ext = "/variation/homo_sapiens"
headers={ "Content-Type" : "application/json", "Accept" : "application/json"}


# query Ensembl for each chunk of 1000
for rs_ids_chunk in rs_ids_chunks:   
    r = requests.post(server+ext, headers=headers, data='{ "ids" : ["' + '", "'.join(rs_ids_chunk) + '"]}')
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()

    # keep track of MAFs
    for i in range(len(rs_ids_chunk)):
        mafs.append(decoded[rs_ids_chunk[i]]['MAF'])
    

# output to new tsv
with open("out_" + sys.argv[1], 'w') as out_tsv:
    with open(sys.argv[1], 'r') as in_tsv:
        line_index = 0
        for line in in_tsv:
            print(line.strip() + '\t' + mafs[line_index], file=out_tsv)
            line_index += 1







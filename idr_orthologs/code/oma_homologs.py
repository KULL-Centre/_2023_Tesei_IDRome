#!/usr/bin/env python3
# (C) 2023 by Kristoffer E. Joahnsson <kristoffer.johansson@bio.ku.dk>

import os
import sys
import csv
import gzip
import pickle
import time

__version__ = 1.00

def is_aa_one_nat(sequence, additional=""):
    for aa in sequence.upper():
        if not (aa in "ACDEFGHIKLMNPQRSTVWY" or aa in additional.upper()):
            return(False)
    return(True)

def read_fasta(filename, comment_char="#;", extra_aa=""):
    """Flexible FASTA file reader without dependencies"""
    seq_list = []
    file_id = filename.split("/")[-1].split(".")[0]
    reading_fasta = False
    with open(filename,"r") as file_handle:
        for line in file_handle:
            words = line.split()
            if len(words) == 0:
                continue
            if words[0][0] in comment_char:
                continue
            if words[0][0] == ">":
                if reading_fasta:
                    seq_list.append((seq_id,seq))
                if len(words[0]) == 1 and len(words) > 1:
                    # this is the case where a space is mistakenly put between '>' and the identifier
                    seq_id = words[1]
                else:
                    seq_id = words[0][1:]
                seq = ""
                reading_fasta = True
            elif reading_fasta:
                if len(words) > 1:
                    print("WARNING: Found FASTA line with more white space separated fields:\n%s" % (line), file=sys.stderr)
                seq = seq + words[0]
            else:
                if len(words) == 1:
                    seq = words[0]
                    if not is_aa_one_nat(seq, extra_aa):
                        print("ERROR: Non-FASTA single-column file should have a protein sequence in first column:", file=sys.stderr)
                        print(line, file=sys.stderr)
                        return(None)
                    seq_id = file_id+"%05d" % (len(seq_list))
                    seq_list.append((seq_id,seq))
                elif len(words) > 1:
                    seq = words[1]
                    if not is_aa_one_nat(seq, extra_aa):
                        print("ERROR: Non-FASTA multi-column file should have a protein sequence in second column:", file=sys.stderr)
                        print(line, file=sys.stderr)
                        return(None)
                    seq_list.append((words[0],seq))
        if reading_fasta:
            seq_list.append((seq_id,seq))
    return(seq_list)


# full sequences of proteins containing IDR's
input_filename = 'idr_proteins.csv'
assert os.path.exists(input_filename)

# wget https://omabrowser.org/All/oma-uniprot.txt.gz
oma_map_filename = 'oma-uniprot.txt.gz'

# wget https://omabrowser.org/All/oma-protein-annotations.txt.gz
oma_protein_filename = 'oma-protein-annotations.txt.gz'

# wget https://omabrowser.org/All/oma-seqs.fa.gz
# gunzip oma-seqs.fa.gz
oma_seq_filename = 'oma-seqs.fa'

# wget https://omabrowser.org/All/oma-pairs.txt.gz
oma_pairs_filename = 'oma-pairs.txt.gz'

outdir = 'human_homologs/'
assert os.path.exists(outdir)

# time total execution
tot_t0 = time.time()

# make a dict of all uniprot entries to map
unip_seq_dict = {}
with open(input_filename, 'r') as csvfile:
    proteins = csv.reader(csvfile)
    
    # first row is the header
    header = next(proteins)
    assert header[0]=='uniprot' and header[1]=='seq'

    for row in proteins:
        assert is_aa_one_nat(row[1])
        unip_seq_dict[row[0]] = row[1]
print("Got %d uniprot sequences from %s" % (len(unip_seq_dict),input_filename))
print("")


print("Reading OMA protein annotations from %s" % (oma_protein_filename))
t0 = time.time()
oma_main_isoform = {}
with gzip.open(oma_protein_filename, 'rt') as protfile:
    for line in protfile:
        # Format: OMA ID<tab>Original ID<tab>Chromosome/Scaffold<tab>Location<tab>MainIsoform<tab>Description
        fields = line.split('\t')
        if fields[0][0] == '#':
            continue
        oma = fields[0]
        # store all human proteins
        if oma[0:5] == 'HUMAN':
            assert not oma in oma_main_isoform.keys()
            if fields[4] == 'self':
                main_isoform = oma
            else:
                main_isoform = fields[4]
            oma_main_isoform[oma] = main_isoform

print("Processed in %.2f seconds" % (time.time()-t0))

with open('oma_isoforms.pickle', 'wb') as pickle_file:
    pickle.dump(oma_main_isoform, pickle_file)
# with open('oma_isoforms.pickle', 'rb') as pickle_file:
#     oma_main_isoform = pickle.load(pickle_file)

        
print("Reading OMA mapping of uniprot entries from %s" % (oma_map_filename))
t0 = time.time()
unip_oma_main_dict = {}
with gzip.open(oma_map_filename, 'rt') as mapfile:
    for line in mapfile:
        # Format: OMA ID<tab>UniProt ID
        # this should split on both tab and whitespace and strip whitespaces
        #   this is most effective then line.split('\t') when assuming a single word between tabs
        fields = line.split()
        if fields[0][0] == '#':
            continue
        oma = fields[0]
        unip = fields[1]
        if unip in unip_seq_dict.keys():
            assert oma[0:5]=='HUMAN'
            # map to main isoform - only these heve annotations of pairs
            oma_main = oma_main_isoform[oma]
            if unip in unip_oma_main_dict.keys():
                if not oma_main in unip_oma_main_dict[unip]:
                    unip_oma_main_dict[unip].append(oma_main)
            else:
                unip_oma_main_dict[unip] = [oma_main]
                                
print("Processed in %.2f seconds" % (time.time()-t0))

with open('unip_oma_dict.pickle', 'wb') as pickle_file:
    pickle.dump(unip_oma_main_dict, pickle_file)
# with open('unip_oma_dict.pickle', 'rb') as pickle_file:
#     unip_oma_main_dict = pickle.load(pickle_file)

print("Found OMA map for %d of %d (%.2f%%) uniprot entries" % (len(unip_oma_main_dict),len(unip_seq_dict),len(unip_oma_main_dict)*1.0/len(unip_seq_dict)*100))
print("")


print("Reading OMA homologs (pairs) from %s" % (oma_pairs_filename))
t0 = time.time()
oma_pairs = {}
with gzip.open(oma_pairs_filename, 'rt') as mapfile:
    for line in mapfile:
        # Format: Protein 1<tab>Protein 2<tab>Orthology type<tab>OMA group (if any)
        # Every pair is listed only once, and in no particular order.
        fields = line.split()
        if fields[0][0] == '#':
            continue
        if fields[0][0:5]=='HUMAN':
            if fields[0] in oma_pairs.keys():
                oma_pairs[fields[0]].append(fields[1])
            else:
                oma_pairs[fields[0]] = [fields[1]]
        elif fields[1][0:5]=='HUMAN':
            if fields[1] in oma_pairs.keys():
                oma_pairs[fields[1]].append(fields[0])
            else:
                oma_pairs[fields[1]] = [fields[0]]

print("Processed in %.2f seconds" % (time.time()-t0))

with open('oma_pairs.pickle', 'wb') as pickle_file:
    pickle.dump(oma_pairs, pickle_file)
# with open('oma_pairs.pickle', 'rb') as pickle_file:
#     oma_pairs = pickle.load(pickle_file)

n_pairs = [len(homol_list) for homol_list in oma_pairs.values()]
mean_homol = sum([n for n in n_pairs])*1.0/len(oma_pairs)
print("Found on average %.2f homologs for %d proteins" % (mean_homol, len(oma_pairs)))
print("")


print("Reading OMA sequences from %s" % (oma_seq_filename))
# It is by far the most effective to read in all sequence and extract them from a dict
# Alternatives tested is to only store needed sequences, from oma_pairs, or using list.index()
t0 = time.time()
seq_list = read_fasta(oma_seq_filename)
size = sys.getsizeof(seq_list) + sum([sys.getsizeof(seq_id)+sys.getsizeof(seq) for (seq_id,seq) in seq_list])
print("Cached all %d OMA sequences in %.2f GB memory" % (len(seq_list), size/10**9))

oma_seq_dict = dict(seq_list)

print("Processed in %.2f seconds" % (time.time()-t0))
print("")

with open('oma_seq_dict.pickle', 'wb') as pickle_file:
    pickle.dump(oma_seq_dict, pickle_file)
# with open('oma_seq_dict.pickle', 'rb') as pickle_file:
#     oma_seq_dict = pickle.load(pickle_file)

size = sys.getsizeof(oma_seq_dict) + sum([sys.getsizeof(seq) for seq in oma_seq_dict.values()])
print("Cached all %d OMA sequences in %.2f GB memory" % (len(oma_seq_dict), size/10**9))


# assemble alignments for mapped proteins and dump files
counter = 0
output_dict = {}
for (unip,oma_list) in unip_oma_main_dict.items():
    t0 = time.time()
    unip_seq = unip_seq_dict[unip]

    # all OMA entries are mapped to main isoforms so all sould have homologs (unless there are none)
    oma_match = None
    oma_with_homologs = [oma for oma in oma_list if oma in oma_pairs.keys()]
    oma_seq_list = []
    if len(oma_with_homologs) == 0:
        print("WARNING: No homologs for %s mapped to %s - skipping" % (unip,str(oma_list)))
        continue
    elif len(oma_with_homologs) == 1:
        oma_match = oma_with_homologs[0]
    else:
        for oma in oma_with_homologs:
            oma_seq = oma_seq_dict[oma]
            # expensive lookup - store for later
            oma_seq_list.append(oma_seq)
            if oma_seq == unip_seq:
                if oma_match is None:
                    oma_match = oma
                else:
                    oma_match_homologs = oma_pairs[oma_match]
                    oma_homologs = oma_pairs[oma]
                    if len(oma_homologs) > len(oma_match_homologs):
                        print("WARNING: Sequence of %s match both %s (%d homologs) and %s (%d homologs), using %s" %
                              (unip,oma_match,len(oma_match_homologs),oma,len(oma_homologs),oma))
                        oma_match = oma
                    else:
                        print("WARNING: Sequence of %s match both %s (%d homologs) and %s (%d homologs), using the first" %
                              (unip,oma_match,len(oma_match_homologs),oma,len(oma_homologs)))
        if oma_match is None:
            assert len(oma_with_homologs) == len(oma_seq_list)
            oma_homologs_list = []
            scores = []
            print("More sets of homologs found for %s but no sequence match" % (unip))
            for i in range(len(oma_with_homologs)):
                oma = oma_with_homologs[i]
                oma_seq = oma_seq_list[i]
                oma_homologs = oma_pairs[oma]
                oma_homologs_list.append(oma_homologs)
                matching_len = len(oma_seq) == len(unip_seq)
                matching_res = 0
                if matching_len:
                    matching_res = sum([aa1==aa2 for (aa1,aa2) in zip(unip_seq,oma_seq)])
                # most homologs wins but if equal, a matching sequences length wins. If more OMS's have both equal, best sequence match wins
                scores.append(len(oma_homologs)*100 + matching_len + matching_res)
                print("    %s: score=%d,  n_homo=%d,  match_len=%d,  match_Res=%d" % (oma,scores[-1],len(oma_homologs),matching_len,matching_res))
            i = scores.index(max(scores))
            oma_match = oma_with_homologs[i]
            print("    select %d %s with score %d" % (i,oma_match,scores[i]))
    oma = oma_match
    oma_seq = oma_seq_dict[oma_match]

    # write fasta file
    homologs = oma_pairs[oma]
    assert len(homologs) > 0
    output_dict[unip] = len(homologs)
    counter += 1
    outfile_filename = "%s/%s_%s.fasta" % (outdir,unip,oma)    
    with open(outfile_filename, 'wt') as outfile:
        outfile.write(">%s_%s\n" % (unip,oma))
        # use original uniprot sequence for query
        outfile.write("%s\n" % (unip_seq))
        for homo_oma in oma_pairs[oma]:
            outfile.write(">%s\n" % (homo_oma))
            outfile.write("%s\n" % (oma_seq_dict[homo_oma]))

total_homologs = sum([val for val in output_dict.values()])
print("")
print("Done writing files for %d of %d (%.2f%%) proteins containing %d homologs. Total time %d seconds" %
      (counter, len(unip_seq_dict), counter*1.0/len(unip_seq_dict)*100, total_homologs,time.time()-tot_t0))

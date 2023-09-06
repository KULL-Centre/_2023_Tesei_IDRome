#!/usr/bin/env python3
# (C) 2023 by Kristoffer E. Joahnsson <kristoffer.johansson@bio.ku.dk>

import sys
import os
import csv
import glob
from math import ceil

__version__ = 1.3

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


if __name__ == "__main__":

    min_homolog_length = 30
    max_homolog_length_factor = 3
    max_homolog_identity = 0.90
    dist_fac_redun = 5.0
    max_total_subst_per_res = 20

    verbose = 2
    # 0: silent
    # 1: report things that should not happen
    # 2: report result per IDR
    # 3: max verbose for debugging
    
    out_filename = 'idr_orthologs.csv'
    alignments_dir = 'human_alignments/'
    
    # make a dict of all uniprot entries to map
    input_filename = 'idr.csv'
    
    idr_list = []
    with open(input_filename, 'rt') as csvfile:
        idr = csv.reader(csvfile, delimiter=';')
        
        # first row is the header
        header = next(idr)
        assert header[0]=='uniprot' and header[11]=='name'

        for row in idr:
            assert is_aa_one_nat(row[6])
            assert row[11] == "%s_%s_%s" % (row[0],int(row[3]),int(row[4]))
            idr_list.append((row[0],int(row[1]),int(row[3]),int(row[4]),row[6]))
    print("Got %d IDR's from %s" % (len(idr_list),input_filename))
    print("")
    
    outfile = open(out_filename,'wt')
    outfile.write("# Script version %.2f\n" % (__version__))
    outfile.write("# Filters:\n")
    outfile.write("#   Sequences with non-natural amino acids removed\n")
    outfile.write("#   min_homolog_length: %d\n" % (min_homolog_length))
    outfile.write("#   max_homolog_length_factor: %.2f\n" % (max_homolog_length_factor))
    outfile.write("#   max_homolog_identity: %.2f\n" % (max_homolog_identity))
    outfile.write("#   dist_fac_redun: %.2f\n" % (dist_fac_redun))
    outfile.write("#   max_total_subst_per_res: %d\n" % (max_total_subst_per_res))
    outfile.write("# Residue numbering is 1-based first and last residue\n")
    outfile.write("#\n")
    outfile.write("%s,%s,%s,%s,%s\n" % ('idr','idr_full_len','idr_ortholog','ortholog_full_len','ortholog_seq'))

    homologs_dict = {}
    tot_homologs_dict = {}
    discarded_is_aa = 0
    discarded_length = 0
    discarded_len_factor = 0
    discarded_similarity = 0    
    discarded_similarity_incl = 0
    discarded_redun = 0
    discarded_tot_subst = 0
    for (unip,nres_unip,resi_first,resi_last,idr_seq) in idr_list:
        idr_name = "%s_%d_%d" % (unip,resi_first,resi_last)
        
        # Read alignment
        search_str = "%s/%s_*_aln.fasta" % (alignments_dir,unip)
        file_list = glob.glob(search_str)
        if len(file_list) == 0:
            if verbose > 0:
                print("WARNING: No homolog file found for %s: %s" % (idr_name,search_str))
            continue
        if len(file_list) > 1:
            if verbose > 0:
                print("WARNING: More homolog files found for %s - using first: %s" % (idr_name,str(file_list)))
        alignment_filename = file_list[0]
        homol_list = read_fasta(alignment_filename)
        if len(homol_list) < 2:
            if verbose > 0:
                print("WARNING: Ortholog file %s only contains %d sequences (first should be query)" % (alignment_filename,len(homol_list)))
            continue

        # Process query
        query = homol_list[0][1]
        n_align_len = len(query)
        assert len(query.replace('-','')) == nres_unip
        
        # Make a map between unaligned residue indices and indices in alignment
        query_map = [i for i in range(len(query)) if query[i] != '-']
        assert len(query_map) == nres_unip
        if verbose > 2:
            print("=== Query for %s: %s with %d residues, alignment length %d" % (idr_name,homol_list[0][0], nres_unip, n_align_len))
            print("  Got %d homologs from %s" % (len(homol_list)-1,alignment_filename))
        
        # Residue numbering is 1-based first and last (counting style), alignment numbering is zero-based begin and end (iterator style)
        align_begin = query_map[resi_first-1]
        align_end = query_map[resi_last-1]+1
        if verbose > 2:
            print("  Original res %d to %d in alignment %d to %d (1-based first and last residue)" % (resi_first,resi_last,align_begin+1,align_end+1))
    
        query_idr = query[align_begin:align_end]
        query_idr_nodash = query_idr.replace('-','')
        if verbose > 2:
            print("  Query:                           %s" % (query_idr_nodash))
        if not query_idr_nodash == idr_seq:
            if verbose > 0:
                print("ERROR: Bad IDR sequence")
                print("  input: %s" % (idr_seq))
                print("  query: %s" % (query_idr_nodash))
            continue

        # Look at all homologs 
        homologs = []
        dist_dict = {}
        tot_subst_dist = 0
        tot_homologs_dict[idr_name] = len(homol_list)-1
        for (seq_id,seq) in homol_list[1:]:
            # In aligned fasta, all sequences should have the same length
            assert len(seq) == n_align_len
            
            # slice homologs to IDR
            homolog_idr = seq[align_begin:align_end]
            # homolog_idr_nodash = ''.join( homolog_idr[i] for i in range(len(homolog_idr)) if homolog_idr[i] != '-')
            homolog_idr_nodash = homolog_idr.replace('-','')
            # homolog_idr_resi_first = len(''.join(seq[:align_begin]).replace('-','')) + 1
            homolog_idr_resi_first = len(seq[:align_begin].replace('-','')) + 1
            homolog_idr_resi_last = homolog_idr_resi_first + len(homolog_idr_nodash) - 1
            homolog_idr_name = "%s_%d_%d" % (seq_id,homolog_idr_resi_first,homolog_idr_resi_last)
            # check if anything aligns
            if len(homolog_idr_nodash) < 1:
                if verbose > 2:
                    print("  Homolog %s of %s has zero length - skipping" % (homolog_idr_name,idr_name))
                continue
            
            # count identical amino acids in aligned IDR's
            identities = sum([homolog_idr[i]==query_idr[i] for i in range(len(query_idr)) if query_idr[i] != '-'])
            # percentage is out of residues in query
            identity = identities*1.0/min(len(query_idr_nodash),len(homolog_idr_nodash))
        
            # Filters
            if not is_aa_one_nat(homolog_idr_nodash):
                discarded_is_aa += 1
                if verbose > 2:
                    print("  Skip homolog %s for IDR %s with non-natural amino acids" % (homolog_idr_name,idr_name))
                    print("    %s" % (homolog_idr_nodash))
                continue    
            if len(homolog_idr_nodash) < min_homolog_length:
                discarded_length += 1
                if verbose > 2:
                    print("  Skip homolog %s for IDR %s with %d < %d residues" % (homolog_idr_name,idr_name,len(homolog_idr_nodash),min_homolog_length))
                continue
            elif len(homolog_idr_nodash) > len(query_idr_nodash)*max_homolog_length_factor or len(homolog_idr_nodash) < ceil(len(query_idr_nodash)/max_homolog_length_factor):
                discarded_len_factor += 1
                if verbose > 2:
                    print("  Skip homolog %s for IDR %s with length %d very different from IDR length %d" %
                          (homolog_idr_name,idr_name,len(homolog_idr_nodash),len(query_idr_nodash)))  
                continue
            elif identity > max_homolog_identity:
                if verbose > 2:
                    print("  Skip homolog %s for IDR %s with identity %.2f > %.2f over %d residues" %
                          (homolog_idr_name,idr_name,identity,max_homolog_identity,min([len(homolog_idr_nodash),len(query_idr_nodash)])))  
                discarded_similarity += 1
                continue

            # output homolog
            homol_full_len = len(seq.replace('-',''))
            homologs.append((homolog_idr_name,homolog_idr_nodash,homolog_idr,identity,identities,homol_full_len))

        if len(homologs) < 1:
            if verbose > 1:
                print("IDR %17s rejected all %d homologs" % (idr_name,len(homol_list)-1))
            continue

        # order homologs according to identity
        def get_identity(tup):
            return tup[3]        
        homologs.sort(reverse=True, key=get_identity)
        if verbose > 2:
            print("  Ordered %d homologs and added closest homolog %s with %.3f identity" % (len(homologs),homologs[0][0],homologs[0][3]))
            
        homologs2 = []
        homol_considered = 0
        for (homolog_idr_name,homolog_idr_nodash,homolog_idr,identity,identities,homol_full_len) in homologs:
            # if everything checks, count identical amino acids to any homolog included so far
            homol_considered += 1
            exclude_homol = False
            for incl_homol_tupl in homologs2:
                incl_homol_idr = incl_homol_tupl[2]
                identities_incl = sum([homolog_idr[i]==incl_homol_idr[i] for i in range(len(homolog_idr)) if homolog_idr[i] != '-'])
                # percentage is out of residues in query
                identity_incl = identities_incl*1.0/min(len(incl_homol_tupl[1]),len(homolog_idr_nodash))
                dist_incl = 1-identity_incl
                dist_query = (1-identity)
                if dist_incl * dist_fac_redun <= dist_query:
                    exclude_homol = True
                    discarded_redun += 1
                    if verbose > 2:
                        print("  Skip redundant homolog %s for IDR %s with dist %.2f*%.2f=%.2f <= %.2f to previously included homolog %s" %
                              (homolog_idr_name,idr_name,dist_incl,dist_fac_redun,dist_incl*dist_fac_redun,dist_query,incl_homol_tupl[0]))
                    break
                
                if identity_incl > max_homolog_identity:
                    exclude_homol = True
                    discarded_similarity_incl += 1
                    if verbose > 2:
                        print("  Skip homolog %s for IDR %s with identity %.2f > %.2f over %d residues to previously included homolog %s" %
                              (homolog_idr_name,idr_name,identity_incl,max_homolog_identity,min([len(homolog_idr_nodash),len(query_idr_nodash)]),incl_homol_tupl[0]))
                    break
                
            if exclude_homol:
                continue
            
            tot_subst_dist += len(query_idr_nodash) - identities
            homologs2.append((homolog_idr_name,homolog_idr_nodash,homolog_idr,identity))
            outfile.write("%s,%d,%s,%d,%s\n" % (idr_name,len(query.replace('-','')),homolog_idr_name,homol_full_len,homolog_idr_nodash))

            # test total number of substitutions among IDR homologs
            if tot_subst_dist > max_total_subst_per_res*len(query_idr_nodash):
                if verbose > 2:
                    print("  IDR %s length %d now have %d total substitutions (%.2f per site) among %d homologs (of %d) - skipping remaining" %
                          (idr_name,len(query_idr_nodash),tot_subst_dist,tot_subst_dist/len(query_idr_nodash),len(homologs2),len(homol_list)-1))
                discarded_tot_subst += len(homologs) - homol_considered
                break

        homologs_dict[idr_name] = homologs2
        if verbose > 1:
            avg_length = sum([len(seq) for (name,seq,seq_aln,identity) in homologs2])*1.0/len(homologs2)
            print("IDR %17s length %3d accepted %3d of %4d (%6.2f%%) homologs of average length %6.2f" %
                  (idr_name,len(query_idr_nodash),len(homologs2),len(homol_list)-1,len(homologs2)/(len(homol_list)-1)*100,avg_length))

    outfile.close()
    accepted_homol = sum([len(homol) for homol in homologs_dict.values()])
    total_homol = sum([n_homol for n_homol in tot_homologs_dict.values()])
    reject_homol = total_homol - accepted_homol
    reject_homol_sum = discarded_is_aa + discarded_length + discarded_len_factor + discarded_similarity + discarded_similarity_incl + discarded_redun + discarded_tot_subst
    print("")
    print("Done writing %d homologs for %d of %d (%.2f%%) IDR's" %
          (accepted_homol,len(homologs_dict),len(idr_list),len(homologs_dict)*1.0/len(idr_list)*100))
    print("Filtered out %d of %d (%.2f%%) orthologs" % (reject_homol,total_homol,reject_homol*1.0/total_homol*100))
    print("  Natural amino acid filter: %5d (%5.2f%%)" % (discarded_is_aa,discarded_is_aa*1.0/total_homol*100))
    print("  Absolute length filter:    %5d (%5.2f%%)" % (discarded_length,discarded_length*1.0/total_homol*100))
    print("  Length factor filter:      %5d (%5.2f%%)" % (discarded_len_factor,discarded_len_factor*1.0/total_homol*100))
    print("  Similarity filter:         %5d (%5.2f%%)" % (discarded_similarity,discarded_similarity*1.0/total_homol*100))
    print("  Similarity to incl filter: %5d (%5.2f%%)" % (discarded_similarity_incl,discarded_similarity_incl*1.0/total_homol*100))
    print("  Redundancy filter:         %5d (%5.2f%%)" % (discarded_redun,discarded_redun*1.0/total_homol*100))
    print("  Enough total subst:        %5d (%5.2f%%)" % (discarded_tot_subst,discarded_tot_subst*1.0/total_homol*100))
    print("  Total:                     %5d (%5.2f%%)" % (reject_homol_sum,reject_homol_sum*1.0/total_homol*100))


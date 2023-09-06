options(width=160, stringsAsFactors=F)

settings = list(version=6.0)

# Make a data frame with one segment per row - this is probably as fast as can get
segregate = function(rome, pome, get_seq=TRUE) {
    # For each protein, get a (string) vector containing indices of the last segment, order to match pome
    agg = aggregate(rome$category, by=list(rome$uniprot), function(v){ r=rle(v); paste0(cumsum(r$lengths), collapse=";") })
    agg = agg[match(pome$uniprot,agg$Group.1), ]

    # a list of same length as pome with vectors of last and first indices of segments
    last_resi_per_protein = lapply(strsplit(agg$x, ";"), as.numeric)
    first_resi_per_protein = lapply(last_resi_per_protein, function(v){ c(1, v[0:(length(v)-1)]+1) })

    # rome indices of proteins and segments
    ri_first_res_per_protein = c(1, cumsum(pome$nres)[1:(nrow(pome)-1)]+1)
    ri_last_res_per_seg = unlist(lapply(seq_along(last_resi_per_protein), function(i){ last_resi_per_protein[[i]] + ri_first_res_per_protein[i] - 1 }))
    nseg = length(ri_last_res_per_seg)
    ri_first_res_per_seg = c(1, ri_last_res_per_seg[1:(nseg-1)] + 1)

    first_resi_per_seg = unlist(first_resi_per_protein)
    last_resi_per_seg = unlist(last_resi_per_protein)

    segments = data.frame(uniprot = rome[ri_last_res_per_seg,"uniprot"],
                         nres_unip = rome[ri_last_res_per_seg,"nres"],
		         nres_seg = last_resi_per_seg - first_resi_per_seg + 1,
                         first = first_resi_per_seg,
                         last = last_resi_per_seg,
		         category = rome[ri_last_res_per_seg,"category"], row.names=NULL)
			 
    if (get_seq) {
        segments$sequence = substring(paste0(rome$aa,collapse=""), ri_first_res_per_seg, ri_last_res_per_seg)
    }
    
    return(segments)
}


# Make a data frame with segments information per residue
segregate_res = function(rome, pome) {
    seg = segregate(rome, pome, get_seq=FALSE)
    
    prev_cat = c("", seg$cat[1:(nrow(seg)-1)])
    prev_cat[which(seg$first==1)] = "N"
    
    next_cat = c(seg$cat[2:nrow(seg)], "")
    next_cat[which(seg$last==seg$nres_unip)] = "C"
    
    seg_res = data.frame(uniprot = rome$uniprot,
                         category = rome$category,
		         segl = rep(seg$nres_seg, seg$nres_seg),
			 prev_cat = rep(prev_cat, seg$nres_seg),
			 next_cat = rep(next_cat, seg$nres_seg))
    return(seg_res)
}


plot_seg = function(seg, xlim=NA, ylim=NA, ...) {
    max_nres = max(seg$nres_seg)
    h_d = hist(seg[which(seg$category=="D"),"nres_seg"], breaks=seq(0,max_nres+1), plot=FALSE)
    h_f = hist(seg[which(seg$category=="F"),"nres_seg"], breaks=seq(0,max_nres+1), plot=FALSE)
    h_g = hist(seg[which(seg$category=="g"),"nres_seg"], breaks=seq(0,max_nres+1), plot=FALSE)
    if (any(is.na(xlim))) { xlim=c(0,max_nres) }
    if (any(is.na(ylim))) { ylim=c(0,max(c(h_d$counts, h_f$counts))) }
    # if (any(is.na(ylim))) { ylim=c(0,max(c(h_d$counts, h_f$counts, h_g$counts))) }
    plot(0,0,col=0, xlim=xlim, ylim=ylim, xlab="Number of residues", ylab="Counts", ...)
    lines(h_g$mids, h_g$counts, lwd=1, col=8)
    lines(h_d$mids, h_d$counts, lwd=2, col=2)
    lines(h_f$mids, h_f$counts, lwd=2, col=3)
    abline(v=c(settings$min_idr_len,settings$max_idr_len), lty=2, col=1)
    legend("topright", c("Disordered","Folded","Gap","Filter"), lty=c(1,1,1,2), lwd=c(2,2,1,1), col=c(2,3,8,1), bg="white")
}


####################################################################
## Proteome
####################################################################
# Load human proteome from uniprot version 2021_4
pome = read.table(gzfile("UP000005640_9606.seq.gz"))
colnames(pome) = c("entry","sequence")
split_list = strsplit(pome$entry, "|", fixed=T)
pome$uniprot = sapply(split_list, "[[", 2)
pome$db = sapply(split_list, "[[", 1)
pome$gene = sapply(split_list, "[[", 3)
pome$nres = nchar(pome$sequence)

## Alphafold pLDDT scores from Akdel et al
# Window smoothed 1-pLDDT from Akdel et al 2021. Smooth should be average of +- 20 residues, 41 res in total
akdel = read.table(gzfile("SDataset1_human_pLDDT_windowed.20.tdt.gz"))
colnames(akdel) = c("uniprot","pLDDT")

i_akdel = match(pome$uniprot,akdel$uniprot)
pome$akdel = akdel[i_akdel,"pLDDT"]
print(sprintf("Mapped %d of %d (%.2f%%) proteins from Akdel pLDDT data to uniprot proteome using uniprot id",
              length(na.omit(i_akdel)), nrow(akdel), length(na.omit(i_akdel))/nrow(akdel)*100))

# i_akdel_unmatched = setdiff(seq(nrow(akdel)),i_akdel)
# print(sprintf("The following %d proteins in the Akdel data are not matched to uniprot proteome", length(i_akdel_unmatched)))
# print(akdel[i_akdel_unmatched,"uniprot"])

# # some proteins have pLDDT value "NA", expand this to a "nan" per residue ("nan,nan,...") to match the format of assigned proteins
# i_na = which(is.na(pome$akdel))
# pome[i_na,"akdel"] = sapply(pome[i_na,"nres"], function(n) {paste(rep("nan",n),collapse=",")})

# # I'm not sure why these are included in the reviewed proteome since they are not
# i_tr = which(pome$db=="tr")

# # There are only pLDDT assignments for sequences >15 aa
# i_short = which(pome$nres <= 15)

# # The UGA codon is made to encode selenocysteine by the presence of a selenocysteine insertion sequence (SECIS) in the mRNA.
# i_sec = unname(which(sapply(pome$sequence, function(s){ any(strsplit(s,"")[[1]]=="U") })))

# # some sequences have 'X'
# i_x = unname(which(sapply(pome$sequence, function(s){ any(strsplit(s,"")[[1]]=="X") })))

# print(sprintf("The following uniprot are from TrEMBL, have selenocysteine (U) or 'X' or <=15 residues (total %d) but have pLDDT assigned (total %d)",
#               length(unique(c(i_short,i_sec,i_tr,i_x))), length(i_na)))
# print( pome[setdiff(c(i_short,i_sec,i_tr,i_x),i_na), "uniprot"] )

# i_unknown_na = setdiff(i_na,c(i_short,i_sec,i_tr,i_x))
# print(sprintf("The following %d uniprot are not matched in Akdel pLDDT data for unknown reasons",
#               length(i_unknown_na)))
# print( pome[i_unknown_na, "uniprot"] )
# print("These may represent mismatch between proteome version but some I don't understand")

# Check that all proteins have same length in uniprot and pLDDT data
akdel_lengths = sapply(strsplit(pome$akdel,","), length)

i_mismatch = which(pome$nres != akdel_lengths)
if (length(i_mismatch) > 0) {
    # print("pLDDT assignments for the following uniprots are omitted due to protein length mismatch between uniprot and Akdel pLDDT data")
    # print(pome[i_mismatch,"uniprot"])
    # print("Lengths in uniprot")
    # print( pome[i_mismatch,"nres"] )
    # print("Lengths of Akdel pLDDT data")
    # print( akdel_lengths[i_mismatch] )
    pome[i_mismatch,"akdel"] = sapply(pome[i_mismatch,"nres"], function(n) {paste(rep("nan",n),collapse=",")})
}

# check that all sequences and pLDDT lists have the same lengths
stopifnot( all( pome$nres == sapply(strsplit(pome$akdel,","), length) ))

has_akdel = sapply(strsplit(pome$akdel,","), "[", 1) != 'nan'
print(sprintf("Mapped %d of %d (%.2f%%) proteins from Akdel to %d of %d (%.2f%%) proteins in uniprot proteome",
              sum(has_akdel), nrow(akdel), sum(has_akdel)/nrow(akdel)*100, sum(has_akdel), nrow(pome), sum(has_akdel)/nrow(pome)*100))

rm(akdel)


####################################################################
## Residue-ome: Data frame with one row per residue
####################################################################
print("Build per residue data frame")
rome = data.frame(uniprot = rep(pome$uniprot, pome$nres),
                  aa = unlist(strsplit(pome$sequence,"")),
		  resi = unlist(sapply(pome$nres,seq)),
		  nres = rep(pome$nres, pome$nres),
		  pLDDT_w20 = 1 - as.numeric( unlist(strsplit(pome$akdel,",")) ))

rome$id = paste0(rome$uniprot, "-", rome$aa, rome$resi)
rome_id_noaa = paste0(rome$uniprot, "-", rome$resi)


####################################################################
## Assign SPOT-Disorder1
####################################################################
# read data from https://github.com/IPritisanac/AF2.IDR/tree/main/FIG1
# here I don't have sequences so I cannot check this
spot1 = read.table(gzfile("UP000005640_9606_SPOTD.out.txt.gz"))
# residue indices seem to 1-based (no index zero observed) and first-last (many IDRs are expected to cover the CT but
#   only few ranges goes beyond the number of residues)
colnames(spot1) = c("name","first","last")
# parse uniprot entry
splitlist = strsplit(spot1$name, "|", fixed=T)
spot1$uniprot = sapply(splitlist, "[", 2)
spot1$nres = pome[match(spot1$uniprot,pome$uniprot),"nres"]

res_vec = function(v){
    resi = seq(v[2],v[3])
    paste0(v[1],"-",seq(v[2],v[3]))
}
spot1_res = unlist(apply(spot1[,c("uniprot","first","last")], MARGIN=1, res_vec))
rome$spot1 = ! is.na(match(rome_id_noaa,spot1_res))

# these are much longer than nres and many-many more have last==nres
i_mismatch = which(spot1$last > spot1$nres)
unip_mismatch = unique(spot1[i_mismatch,"uniprot"])
rome[which(rome$uniprot %in% unip_mismatch),"spot1"] = NA

n_match = nrow(spot1) - length(i_mismatch)
print(sprintf("Matched %d of %d (%.2f%%) SPOT-Disorder1 IDR's to %d (%.2f%%) proteome residues with %d (%.2f%%) disordered residues",
              n_match, nrow(spot1), n_match/nrow(spot1)*100, sum(! is.na(rome$spot1)), sum(! is.na(rome$spot1))/nrow(rome)*100,
	      sum(rome$spot1, na.rm=T), sum(rome$spot1, na.rm=T)/sum(!is.na(rome$spot1))*100))
rm(spot1, spot1_res)


####################################################################
## Category based disorder annotation - dual1
####################################################################
# everything with pLDDT==NA ends up being folded?
print("Dual threshold disorder categories from Alphafold")

settings$disorder_threshold = 0.7
settings$folded_threshold = 0.8

rome$category = "g"
rome[which(rome$pLDDT_w20 <= settings$disorder_threshold),"category"] = "D"
rome[which(rome$pLDDT_w20 >= settings$folded_threshold),"category"] = "F"

quartz(width=8, height=12)
par(mfcol=c(3,1), bg="white")
seg = segregate(rome, pome) # per segment (not residue) needed for plotting
plot_seg(seg, xlim=c(1,600), main="Before filtering")

seg = segregate_res(rome, pome)

# print("Make cache"); save(pome,rome,seg,settings,file="idr_cache_dual.rda")
# print("Read cache"); load("idr_cache_dual.rda")
# print("Done")

# Remove small folded or disordered regions
settings$min_segl = 10
i1 = which(seg$segl < settings$min_segl & seg$cat %in% c("F","D"))
print(sprintf("Remove %d folded or disordered segments shorter then %d residues",length(i1), settings$min_segl))
rome[i1,"category"] = "g"

seg = segregate(rome, pome) # per segment (not residue) needed for plotting
plot_seg(seg, xlim=c(1,600), main="Short removed")

seg = segregate_res(rome, pome)

i1_recalc = which(seg$segl < settings$min_segl & seg$cat %in% c("F","D"))
print(sprintf("Recalc %d",length(i1_recalc)))

# Gap regions between disordered annotated as disordered
i2 = which(seg$cat == "g" & seg$prev_cat %in% c("N","D") & seg$next_cat %in% c("C","D"))
print(sprintf("Remove %d gap segments between disordered regions",length(i2)))
rome[i2,"category"] = "D"

# Gap regions between folded annotated as folded - also covers proteins that only has gap category residues
i3 = which(seg$cat == "g" & seg$prev_cat %in% c("N","F") & seg$next_cat %in% c("C","F"))
print(sprintf("Remove %d gap segments between folded regions",length(i3)))
rome[i3,"category"] = "F"

# Gap regions between folded and disordered are kept as gap regions

seg = segregate_res(rome, pome)
i2_recalc = which(seg$cat == "g" & seg$prev_cat %in% c("N","D") & seg$next_cat %in% c("C","D"))
i3_recalc = which(seg$cat == "g" & seg$prev_cat %in% c("N","F") & seg$next_cat %in% c("C","F"))
print(sprintf("Recalc %d and %d",length(i2_recalc),length(i3_recalc)))


####################################################################
## Data frame of disordered regions
####################################################################

seg = segregate(rome, pome)

plot_seg(seg, xlim=c(1,600), main="Intermediate gaps removed")
quartz.save("filtering.png", type="png")

idr = seg[which(seg$category=="D"),]

# Filter IDR on length
settings$min_idr_len = 30
settings$max_idr_len = 1500

# use IDR within thresholds as they are
idr_filter = idr[which(idr$nres_seg >= settings$min_idr_len & idr$nres_seg <= settings$max_idr_len),]


### Assign a IDR segment name to each residue
stopifnot("id" %in% colnames(rome))
id_list = apply(idr_filter[,c("uniprot","first","last","sequence")], MARGIN=1, function(l){ paste0(l[1],"-",strsplit(l[4],"")[[1]],seq(l[2],l[3])) })
id_list = unname(unlist(id_list))

idr_filter$name = paste0(idr_filter$uniprot,"_",idr_filter$first,"_",idr_filter$last)
idr_list = apply(idr_filter[,c("name","nres_seg")], MARGIN=1, function(l){ rep(l[1],l[2]) })
idr_list = unname(unlist(idr_list))
stopifnot(length(id_list) == length(idr_list))
rome[match(id_list,rome$id),"idr"] = idr_list

# add per-residue assignments
i = which(rome$resi==1)
protein_offset = i-1
names(protein_offset) = rome[i,"uniprot"]

line_field_collapsed = function(i, field, collapse=":") {
    i1 = protein_offset[idr_filter[i,"uniprot"]] + idr_filter[i,"first"]
    i2 = protein_offset[idr_filter[i,"uniprot"]] + idr_filter[i,"last"]
    paste0(rome[i1:i2,field], collapse=collapse)
}
idr_filter$pLDDT_w20 = sapply(seq(nrow(idr_filter)), line_field_collapsed, field="pLDDT_w20")

# Load QCDpred predictions
# write.table(idr_filter[,c("name","sequence")], file="idr.seq", quote=F, row.names=F, col.names=F)
# ../../tools/QCDpred.py idr.seq > idr_qcdpred.txt
qcdpred = read.table(gzfile("idr_qcdpred.txt.gz"))
colnames(qcdpred) = c("name","seq","score","aa","resi")

# assume they are ordered!
agg = aggregate(qcdpred$score, by=list(qcdpred$name), paste0, collapse=":")

# consider adding 8 NA's to start and end!
idr_filter$QCDpred = agg[match(idr_filter$name,agg[,1]),2]

# check that all IDRs have something assigned
stopifnot(sum(is.na(idr_filter$QCDpred)) == 0)


####################################################################
## Store 
####################################################################

# Store all data as r-data object
save(idr_filter, seg, rome, pome, settings, file="idr.rda")

# # Store un-filtered IDR as csv
# write.table(idr, "idr_prefilter.csv", row.names=F, quote=F, sep=";")

# Store filtered IDR as csv
write.table(idr_filter, "idr.csv", row.names=F, quote=F, sep=";")

# Dump full the full sequences of proteins containing IDR
idr_proteins = data.frame(uniprot=unique(idr_filter$uniprot))
idr_proteins$seq = pome[match(idr_proteins$uniprot,pome$uniprot),"sequence"]
write.csv(idr_proteins, file="idr_proteins.csv", quote=F, row.names=F)


####################################################################
## Alternative IDRome definition
####################################################################
# save original categories
rome$cat_plddt = rome$category
rome$category = ifelse(rome$spot1, "D", "F")

quartz(width=8, height=12)
par(mfcol=c(3,1), bg="white")
seg = segregate(rome, pome) # per segment (not residue) needed for plotting
plot_seg(seg, xlim=c(1,600), main="ALT: Before filtering")

# seg = segregate_res(rome, pome)

seg = segregate(rome, pome)

idr_alt = seg[which(seg$category=="D"),]

# use IDR within thresholds as they are
idr_alt_filter = idr_alt[which(idr_alt$nres_seg >= settings$min_idr_len & idr_alt$nres_seg <= settings$max_idr_len),]

idr_alt_filter$name = paste0(idr_alt_filter$uniprot,"_",idr_alt_filter$first,"_",idr_alt_filter$last)

# add pLDDT strings 
first_aa_per_seg = substring(idr_alt_filter$sequence,1,1)
ri_first_per_seg = match(paste0(idr_alt_filter$uniprot,"-",first_aa_per_seg,idr_alt_filter$first),rome$id)
stopifnot(all(! is.na(ri_first_per_seg)))

nc = nchar(idr_alt_filter$sequence)
last_aa_per_seg = substring(idr_alt_filter$sequence,nc,nc)
ri_last_per_seg = match(paste0(idr_alt_filter$uniprot,"-",last_aa_per_seg,idr_alt_filter$last),rome$id)
stopifnot(all(! is.na(ri_last_per_seg)))

idr_alt_filter$pLDDT_w20 = apply(cbind(ri_first_per_seg,ri_last_per_seg), MARGIN=1, function(v){ paste0(rome[seq(v[1],v[2]),"pLDDT_w20"],collapse=":")})
# idr_alt_filter$plddt_cat = apply(cbind(ri_first_per_seg,ri_last_per_seg), MARGIN=1, function(v){ paste0(rome[seq(v[1],v[2]),"cat_plddt"],collapse="")})

# dump
write.table(idr_alt_filter, "idr_spot1_simple.csv", row.names=F, quote=F, sep=";")

# # Store un-filtered IDR as csv
# write.table(idr_alt, "idr_spot1_simple_prefilter.csv", row.names=F, quote=F, sep=";")

print(sprintf("Alternative IDRome has %d (vs %d) IDR's with %d of %d (%.2f%%) residues disordered",
              nrow(idr_alt_filter), nrow(idr_filter), sum(idr_alt_filter$nres_seg), sum(!is.na(rome$category)),
	      sum(idr_alt_filter$nres_seg)/sum(!is.na(rome$category))*100 ))
print("IDR length distribution of pLDDT and SPOT-Disorder1 IDRome:")
print(summary(idr_filter$nres_seg))
print(summary(idr_alt_filter$nres_seg))

unip_new = unique(idr_alt_filter[which(! idr_alt_filter$uniprot %in% unique(idr_filter$uniprot)),"uniprot"])
print(sprintf("Alternative IDRome has %d additional uniprot entries", length(unip_new)))


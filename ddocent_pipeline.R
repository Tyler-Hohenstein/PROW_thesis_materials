library(HardyWeinberg)
library(dplyr)
library(genetics)
library(hierfstat)
library(ggplot2)
library(pcadapt)
library(qvalue)
library(Biostrings)
library(raster)
library(sp)
library(psych)
library(vegan)
library(egg)
library(geosphere)
library(Rsamtools)
library(genbankr)
library(rtracklayer)



##### functions #####

# observed heterozygosity
obsH = function(dat) {
  obs.count = length(which(dat == 1))
  total.counts = table(is.na(dat))
  comp.dat = total.counts[[1]]
  return(obs.count/comp.dat)
}

# expected heterozygosity
expH = function(dat) {
  snpv = as.vector(table((dat[!is.na(dat)])))
  q = maf(snpv)
  p = 1 - q
  two_pq = 2*p*q
  
  return(two_pq)
}

# F inbreeding coefficient
Fic = function(dat) {
  ans = (expH(dat) - obsH(dat)) / expH(dat) 
  return(ans)
}

# Fic for all SNPs, output to a new data frame
F.inbred.all = function(inp, out) {
  gen.names = colnames(inp)
  for(i in 1:length(gen.names)) {
    #print(i)
    F_i = Fic(inp[,i])
    out$F[i] = F_i
  }
  return(out)
}



##### Global Variables #####

`%ni%` = Negate(`%in%`)



##### set wd and inputs #####

#home_dir = readline(prompt="Enter home directory: ")
#input = readline((prompt="Enter 012 input (eg. '45_ref'): "))

home_dir = "C:/Users/Tyler" #edit this to your home directory
input = "ddocent_60"

wd = paste(home_dir, "/Desktop/Graduate School/VCU/Thesis/Scripts/R_env", sep="") #edit this to your working directory
setwd(wd)

SNPs_fh = paste(input, "_raw.012", sep="")
pos_fh = paste(input, "_raw.012.pos", sep="")
indv_fh = paste(input, "_raw.012.indv", sep="")
depth_fh = paste(input, "_raw.gdepth", sep="")
geo_fh = "PROW_lat_long.csv" 
SNPs_out_fh = paste("PROW_", input, "_filtered.012", sep="")
pos_out_fh = paste("PROW_", input, "_filtered.012.pos", sep="")
indv_out_fh = paste("PROW_", input, "_filtered.012.indv", sep="")
SNP_names_out_fh = paste("PROW_", input, "_filtered.012.SNPs_names", sep="")
struct_out_fh = paste("PROW_", input, "_structure.txt", sep="")
ped_out_fh = paste(input, "_filtered.ped", sep="")
map_out_fh = paste(input, "_filtered.map", sep="")



##### read in and clean up files #####

print("reading in and cleaning up files...")

SNPs = read.table(SNPs_fh, header=FALSE, sep="\t")
SNPs = SNPs[,2:length(colnames(SNPs))]

depth = read.table(depth_fh, header=TRUE, sep="\t")

pos = read.table(pos_fh, header=FALSE, sep="\t")

indv = read.table(indv_fh, header=FALSE, sep="\t")

lat_long = read.table('PROW_lat_long.txt', header=TRUE, sep="\t")

raw_consensus_fasta = readDNAStringSet("reference.fasta")
raw_contig_names = names(raw_consensus_fasta)
raw_contig_seqs = paste(raw_consensus_fasta)
  

# set row and column names

rownames(SNPs) = indv$V1

SNP_names = c()
for(i in 1:length(colnames(SNPs))) {
  cur_SNP = paste("SNP_", i, sep="")
  #print(cur_SNP)
  SNP_names[i] = cur_SNP
}

colnames(SNPs) = SNP_names
rownames(pos) = SNP_names
colnames(pos) = c("contig", "position") # subpopulation labels corresponding to sample site have been added to the .indv file in Excel
colnames(indv) = c("indv", "subpop")
rownames(indv) = rownames(SNPs)
rownames(depth) = SNP_names

rownames(lat_long) = lat_long$Site_Name

raw_SNPs = SNPs
raw_depth = depth
Unfiltered_indv = indv
Unfiltered_pos = pos



##### Set all -1 values to NAs #####

SNPs[SNPs == -1] = NA
Unfiltered_SNPs = SNPs
depth[depth == -1] = NA
Unfiltered_depth = depth



##### FILTERING #####

#filter for depth of 4 or greater for homozygous SNP calls

#print("filtering for read depth of 4 or greater for all homozygous genotype calls...")

depth_slim = depth[,3:length(depth[1,])]

#for(i in 1:length(rownames(SNPs))) {
#  #print(i)
#  for(j in 1:length(colnames(SNPs))) {
#    #print(paste(i, j, sep=" "))
#    if(!is.na(depth_slim[j,i])) {
#      if(depth_slim[j,i] < 4) {
#        if(SNPs[i,j] != 1) {
#          SNPs[i,j] = NA
#        }
#      }
#    }
#  }
#}



# remove contigs that don't align to PROW reference genome

print("removing SNPs on unplaced contigs...")

unplaced_reads = readDNAStringSet("ddocent_mywa_contigs_unplaced.fa")
unplaced_contig_names = names(unplaced_reads)

pos = pos[pos$contig %ni% unplaced_contig_names,]
SNPs = SNPs[, rownames(pos)]
SNPs_alignment_to_keep = colnames(SNPs)
depth = depth[SNPs_alignment_to_keep,]
depth_slim = depth_slim[SNPs_alignment_to_keep,]



# remove any sample with >90% missing data

print("filtering out all samples with >90% missing data...")

sample_stats = data.frame(matrix(nrow=length(rownames(SNPs)), ncol=0))
rownames(sample_stats) = indv$indv

sample_stats$percent_missing = NA
for(i in 1:length(rownames(SNPs))) {
  #print(SNPs[,i])
  n = 0
  for(j in SNPs[i,]) {
    #print(j)
    if(is.na(j)) {
      n = n + 1
    }
  }
  perc_miss = n/length(colnames(SNPs))
  #print(perc_miss)
  sample_stats$percent_missing[i] = perc_miss
}

samples_pm_to_toss = filter(sample_stats, percent_missing > 0.9)
sample_stats = filter(sample_stats, percent_missing <= 0.9)
samples_pm_to_keep = rownames(sample_stats)
SNPs = SNPs[samples_pm_to_keep,]
indv = indv[samples_pm_to_keep,]
samples_pm_rm = rownames(samples_pm_to_toss)
depth = subset(depth, select=colnames(depth) %ni% samples_pm_rm)
depth_slim = subset(depth_slim, select=colnames(depth_slim) %ni% samples_pm_rm)



# remove any SNP with >50% missing data

print("filtering out all SNPs with >50% missing data...")

SNP_stats = data.frame(matrix(nrow=length(colnames(SNPs)), ncol=0))
rownames(SNP_stats) = colnames(SNPs)
SNP_stats$contig = pos$contig
SNP_stats$position = pos$position
SNP_stats$percent_missing = NA

for(i in 1:length(colnames(SNPs))) {
  #print(SNPs[,i])
  n = 0
  for(j in SNPs[,i]) {
    #print(j)
    if(is.na(j)) {
      n = n + 1
    }
  }
  perc_miss = n/length(rownames(SNPs))
  #print(perc_miss)
  SNP_stats$percent_missing[i] = perc_miss
}

SNPs_pm_to_toss = filter(SNP_stats, percent_missing > 0.5)
SNP_stats = filter(SNP_stats, percent_missing <= 0.5)
SNPs_pm_to_keep = rownames(SNP_stats)
SNPs = subset(SNPs, select=SNPs_pm_to_keep)
depth = depth[SNPs_pm_to_keep,]
depth_slim = depth_slim[SNPs_pm_to_keep,]
pos = pos[SNPs_pm_to_keep,]



# remove all SNPs with MAF < 0.01

print("filtering for SNPs with MAF < 0.01...")

SNP_stats$MAF = NA

for(i in 1:length(colnames(SNPs))) {
  #print(i)
  SNP = SNPs[,i]
  snpv = as.vector(table((SNP[!is.na(SNP)])))
  x = maf(snpv)
  SNP_stats$MAF[i] = x   
}

SNPs_maf_to_toss = filter(SNP_stats, MAF < 0.01)
SNP_stats = filter(SNP_stats, MAF >= 0.01)
SNPs_maf_to_keep = rownames(SNP_stats)
SNPs = subset(SNPs, select=SNPs_maf_to_keep)
depth = depth[SNPs_maf_to_keep,]
depth_slim = depth_slim[SNPs_maf_to_keep,]
pos = pos[SNPs_maf_to_keep,]



# remove any SNP for which F < -0.5

print("filtering out all SNPs for which F < -0.5...")

SNP_stats = F.inbred.all(SNPs, SNP_stats)

SNPs_F_to_toss = filter(SNP_stats, F < -0.5)
SNP_stats = filter(SNP_stats, F >= -0.5)
SNPs_F_to_keep = rownames(SNP_stats)
SNPs = subset(SNPs, select=SNPs_F_to_keep)
depth = depth[SNPs_F_to_keep,]
depth_slim = depth_slim[SNPs_F_to_keep,]
pos = pos[SNPs_F_to_keep,]



##### CALCULATIONS #####

# calculate Fst

print("calculating Fst and Fis for all SNPs...")

subpop = indv$subpop
SNPs = cbind(subpop, SNPs)

SNP_stats$Fst = NA
SNP_stats$Fis = NA
for(i in 2:length(SNPs[1,])) {
  #print(i)
  SNP_stats$Fst[i-1] = wc(ndat=SNPs[,c(1,i)])$FST
  SNP_stats$Fis[i-1] = wc(ndat=SNPs[,c(1,i)])$FIS
}

SNPs = SNPs[,2:length(colnames(SNPs))]

ggplot(SNP_stats, aes(F)) + geom_area(stat="bin")
ggplot(SNP_stats, aes(Fis)) + geom_area(stat="bin")
ggplot(SNP_stats, aes(Fst)) + geom_area(stat="bin")

print(cat("Fst of SNPs for which F > 0.5:", 
          round(SNP_stats[which(SNP_stats$F > 0.5),"Fst"], digits=4), sep="\t"))
print(cat("Fst of SNPs for which F < 0.5:", 
          round(SNP_stats[which(SNP_stats$F <= 0.5),"Fst"], digits=4), sep="\t"))

Fst_for_F_over_0.5 = data.frame("Fst"=SNP_stats[which(SNP_stats$F > 0.5),"Fst"])
Fst_for_F_under_0.5 = data.frame("Fst"=SNP_stats[which(SNP_stats$F <= 0.5),"Fst"])

ggplot(Fst_for_F_over_0.5, aes(Fst)) + geom_area(stat="bin")
ggplot(Fst_for_F_under_0.5, aes(Fst)) + geom_area(stat="bin")

boxplot(SNP_stats[which(SNP_stats$F <= 0.5),"Fst"], 
        SNP_stats[which(SNP_stats$F > 0.5), "Fst"], 
        names=c("F < 0.5", "F > 0.5"), ylab="Fst")
boxplot(SNP_stats[which(SNP_stats$F <= 0.5),"Fis"], 
        SNP_stats[which(SNP_stats$F > 0.5), "Fis"], 
        names=c("F < 0.5", "F > 0.5"), ylab="Fis")



# calculate average depth

print("calculating average read depth for all SNPs...")

SNP_stats$avg_depth = NA
for(i in 2:(length(colnames(SNPs))+1)) {
  #print(paste("i:", i, sep=" "))
  depth_list = c()
  for(j in 1:length(rownames(SNPs))) {
    #print(paste("j:", j, sep=" "))
    #print(paste("depth:" , depth[i-1,j+2], sep=" "))
    depth_list = append(depth_list, depth[i-1,j+2])
  }
  #print(depth_list)
  #print(mean(depth_list, na.rm=TRUE))
  SNP_stats$avg_depth[i-1] = mean(depth_list, na.rm=TRUE)
}

ggplot(SNP_stats, aes(avg_depth)) + geom_area(stat="bin")+
  xlab("Avg. SNP Read Depth")

print(cat("Avg. Read Depth of SNPs for which F > 0.5:", 
          round(SNP_stats[which(SNP_stats$F > 0.5), "avg_depth"], digits=4), sep="\t"))
print(cat("Avg. Read Depth of SNPs for which F < 0.5:", 
          round(SNP_stats[which(SNP_stats$F <= 0.5),"Fst"], digits=4), sep="\t"))

Avg_depth_for_F_over_0.5 = data.frame("avg_depth"=SNP_stats[which(SNP_stats$F > 0.5), "avg_depth"])
Avg_depth_for_F_under_0.5 = data.frame("avg_depth"=SNP_stats[which(SNP_stats$F <= 0.5), "avg_depth"])

ggplot(Avg_depth_for_F_over_0.5, aes(avg_depth)) + geom_area(stat="bin")+
  xlab("SNP Avg. Depth")
ggplot(Avg_depth_for_F_under_0.5, aes(avg_depth)) + geom_area(stat="bin") +
  xlab("SNP Avg. Depth")

boxplot(SNP_stats[which(SNP_stats$F <= 0.5),"avg_depth"], 
        SNP_stats[which(SNP_stats$F > 0.5), "avg_depth"], 
        names=c("F < 0.5", "F > 0.5"), 
        ylab="Average Coverage Depth per SNP")



# remove any SNP for which F > 0.5

print("filtering out all SNPs for which F > 0.5...")

SNPs_F_to_toss = filter(SNP_stats, F > 0.5)
SNP_stats = filter(SNP_stats, F <= 0.5)
SNPs_F_to_keep = rownames(SNP_stats)
SNPs = subset(SNPs, select=SNPs_F_to_keep)
depth = depth[SNPs_F_to_keep,]
depth_slim = depth_slim[SNPs_F_to_keep,]
pos = pos[SNPs_F_to_keep,]



# calculate pairwise linkage disequilibrium

#print("calculating pairwise linkage disequilibrium...")

#SNPs_genotypes = makeGenotypes(SNPs, convert=colnames(SNPs), method=as.genotype.allele.count)

#SNP_1 = c()
#SNP_2 = c()
#SNPs_old = c()
#linkage = c()
#rownames_LD = c()
#n = c()

#for(i in rownames(SNP_stats)) {
#  SNPs_old = c(SNPs_old, i)
#  for(j in rownames(SNP_stats)) {
#    if(i != j) {
#      if(SNP_stats[i, "contig"] == SNP_stats[j, "contig"]) {
#        if(j %ni% SNPs_old) {
#          #print(paste(i, j, sep="-"))
#          SNP_1 = c(SNP_1, i)
#          SNP_2 = c(SNP_2, j)
#          x = LD(SNPs_genotypes[,i], SNPs_genotypes[,j])
#          linkage = c(linkage, x$`R^2`)
#          rownames_LD = c(rownames_LD, paste(i, j, sep="-"))
#        }
#      }
#    }
#  }
#}

#LD_table = data.frame(SNP_1, SNP_2, linkage)
#rownames(LD_table) = rownames_LD

#print("Filtering SNPs with linkage R^2 > 0.45, retaining SNPs out of pairs with the least missing data...")

#SNPs_LD_to_toss = c()
#for(i in 1:length(rownames(LD_table))) {
#  if(LD_table$linkage[i] > 0.45) {
#    SNPs_pair = c(LD_table$SNP_1[i], LD_table$SNP_2[i])
#    perc_miss_pair = c(SNP_stats[SNPs_pair[1], "percent_missing"], SNP_stats[SNPs_pair[2], "percent_missing"])
#    print(paste(SNPs_pair[1], perc_miss_pair[1], SNPs_pair[2], perc_miss_pair[2], sep = " "))
#    if(perc_miss_pair[1] < perc_miss_pair[2]) {
#      SNPs_LD_to_toss = c(SNPs_LD_to_toss, SNPs_pair[1])
#    }
#    if(perc_miss_pair[2] < perc_miss_pair[1]) {
#      SNPs_LD_to_toss = c(SNPs_LD_to_toss, SNPs_pair[2])
#    }
#  }
#}

#SNPs_LD_to_toss
#SNPs = subset(SNPs, select=names(SNPs) %ni% SNPs_LD_to_toss)
#SNPs_LD_to_keep = colnames(SNPs)
#depth = depth[SNPs_LD_to_keep,]
#depth_slim = depth_slim[SNPs_LD_to_keep,]
#SNP_stats = SNP_stats[SNPs_LD_to_keep,]
#pos = pos[SNPs_LD_to_keep,]

#sample_stats.1 = sample_stats
#SNPs.1 = SNPs
#indv.1 = indv
#depth.1 = depth
#depth_slim.1 = depth_slim



##### filter for one SNP per contig #####

print("Filtering for one SNP per contig...")

contigs = unique(SNP_stats$contig)
SNPs_old = c()
SNPs_rad_to_toss = c()

for(n in contigs) {
  for(i in rownames(SNP_stats)[which(SNP_stats$contig == n)]) {
    SNPs_old = c(SNPs_old, i)
    for(j in rownames(SNP_stats)[which(SNP_stats$contig == n)]) {
      if(i != j) {
        if(j %ni% SNPs_old) {
          if(SNP_stats[i, "percent_missing"] > SNP_stats[j, "percent_missing"]) {
            SNPs_rad_to_toss = c(SNPs_rad_to_toss, i)
          }
          if(SNP_stats[i, "percent_missing"] < SNP_stats[j, "percent_missing"]) {
            SNPs_rad_to_toss = c(SNPs_rad_to_toss, j)
          }
          if(SNP_stats[i, "percent_missing"] == SNP_stats[j, "percent_missing"]) {
            coinflip = sample(c(0,1), size=1, prob=c(0.5,0.5))
            if(coinflip == 0) {
              SNPs_rad_to_toss = c(SNPs_rad_to_toss, i)
            }
            if(coinflip == 1) {
              SNPs_rad_to_toss = c(SNPs_rad_to_toss, j)
            }
          }
        }
      }
    }
  }
}

#SNPs_rad_to_toss
SNPs = subset(SNPs, select=names(SNPs) %ni% SNPs_rad_to_toss)
SNPs_rad_to_keep = colnames(SNPs)
depth = depth[SNPs_rad_to_keep,]
depth_slim = depth_slim[SNPs_rad_to_keep,]
SNP_stats = SNP_stats[SNPs_rad_to_keep,]
pos = pos[SNPs_rad_to_keep,]

sample_stats.1 = sample_stats
SNPs.1 = SNPs
indv.1 = indv
depth.1 = depth
depth_slim.1 = depth_slim
pos.1 = pos



##### generate multi-locus PCA plots #####

print("generating multi-locus PCA plots...")

u_count = c()
for(j in 1:length(colnames(SNPs))) {
  u_j = sum(SNPs[,j], na.rm=TRUE) / sum(!is.na(SNPs[,j]))
  u_count = c(u_count, u_j)
}

SNPs_c_s = SNPs
for(j in 1:length(colnames(SNPs))) {
  u_j = u_count[j]
  p_j = u_j / 2
  for(i in 1:length(rownames(SNPs))) {
    c_ij = SNPs_c_s[i,j]
    SNPs_c_s[i,j] = (c_ij - u_j) / (sqrt(p_j * (1 - p_j)))
  }
}

SNPs_c_s[is.na(SNPs_c_s)] = 0

pca = prcomp(x=SNPs_c_s, center=FALSE)
eig_sdev = pca$sdev
eig_var = eig_sdev^2

# percent variance explained
PVE = (eig_var/sum(eig_var))*100
print(cat("Percent Variance Explained:", 
          round(PVE, digits=4), sep="\t"))

# Plot 1st & 2nd PCs
prin_comps = data.frame(pca$x)
prin_comps = cbind(indv$subpop, prin_comps)
names(prin_comps)[1] = "subpop"

for(i in 1:length(prin_comps$subpop)) {
  x = prin_comps[i,1]
  if(x=="PopN"|x=="PopM"|x=="PopL"|x=="PopZ"|x=="PopK"|x=="PopJ"|x=="PopO") {
    prin_comps[i,1] = "West"
  }
  else {
    prin_comps[i,1] = "East"
  }
}

prin_comps = cbind(indv$subpop, prin_comps)
names(prin_comps)[1] = "state"

for(i in 1:length(prin_comps$subpop)) {
  x = prin_comps[i,1]
  if(x=="PopO") {
    prin_comps[i,1] = "Wisconsin"
  }
  if(x=="PopJ") {
    prin_comps[i,1] = "Ohio"
  }
  if(x=="PopN") {
    prin_comps[i,1] = "Arkansas"
  }
  if(x=="PopM" | x=="PopL" | x=="PopK") {
    prin_comps[i,1] = "Louisiana"
  }
  if(x=="PopZ") {
    prin_comps[i,1] = "Unknown"
  }
  if(x=="PopH") {
    prin_comps[i,1] = "South Carolina"
  }
  if(x=="PopG") {
    prin_comps[i,1] = "North Carolina"
  }
  if(x=="PopA" | x=="PopB" | x=="PopE" | x=="PopF") {
    prin_comps[i,1] = "Virginia"
  }
}

prin_comps = cbind(indv$subpop, prin_comps)
names(prin_comps)[1] = "site"

for(i in 1:length(prin_comps$site)) {
  x = prin_comps[i,1]
  if(x=="PopA") {
    prin_comps[i,1] = "VA, Deep Bottom"
  }
  if(x=="PopB") {
    prin_comps[i,1] = "VA, Dragon Run"
  }
  if(x=="PopE") {
    prin_comps[i,1] = "VA, Fort AP Hill"
  }
  if(x=="PopF") {
    prin_comps[i,1] = "VA, Great Dismal Swamp"
  }
  if(x=="PopG") {
    prin_comps[i,1] = "NC, Holt Lake"
  }
  if(x=="PopH") {
    prin_comps[i,1] = "SC, Francis Beidler Forest"
  }
  if(x=="PopJ") {
    prin_comps[i,1] = "OH, Hoover Reservoir"
  }
  if(x=="PopK") {
    prin_comps[i,1] = "LA, Palmetto Island"
  }
  if(x=="PopL") {
    prin_comps[i,1] = "LA, Bluebonnet Swamp"
  }
  if(x=="PopM") {
    prin_comps[i,1] = "LA, Barataria Preserve"
  }
  if(x=="PopN") {
    prin_comps[i,1] = "AR, White River"
  }
  if(x=="PopO") {
    prin_comps[i,1] = "WI, Sugar River"
  }
}

prin_comps$subpop = as.factor(prin_comps$subpop)
prin_comps$state = as.factor(prin_comps$state)
prin_comps$state = factor(prin_comps$state, 
                          levels=c("Virginia", "North Carolina",
                                   "South Carolina", "Louisiana", 
                                   "Unknown", "Arkansas", "Ohio",
                                   "Wisconsin"))
prin_comps$site = as.factor(prin_comps$site)
prin_comps$site = factor(prin_comps$site, 
                           levels=c("VA, Deep Bottom", "VA, Dragon Run",
                                    "VA, Fort AP Hill", "VA, Great Dismal Swamp",
                                    "NC, Holt Lake", "SC, Francis Beidler Forest",
                                    "LA, Palmetto Island", "LA, Bluebonnet Swamp",
                                    "LA, Barataria Preserve", "AR, White River",
                                    "OH, Hoover Reservoir", "WI, Sugar River"))

#plots
ggplot(prin_comps, 
       aes(x=PC1, y=PC2,
           color = as.factor(subpop))) + 
  xlab("PC1 (0.896)") +
  ylab("PC2 (0.853%)") + 
  geom_point(alpha=0.8, size=2) +
  geom_point(shape=1, alpha=0.3, color="black", size=2) +
  labs(color="Subregion") + 
  geom_hline(yintercept=0 , color="black", lty="dashed") +
  geom_vline(xintercept=0 , color="black", lty="dashed") 
ggplot(prin_comps, 
       aes(x=PC1, y=PC2,
           color = as.factor(state))) + 
  xlab("PC1 (0.896%)") +
  ylab("PC2 (0.853%)") + 
  geom_point(alpha=0.8, size=3.3) + 
  geom_point(shape=1, alpha=0.3, color="black", size=3.3) +
  labs(color="State") + 
  scale_color_manual(values = c("#d53e4f", "#fc8d59", 
                                "#fee08b", "#ffffbf", 
                                "#e6f598", "#99d594",
                                "#3288bd")) +
  geom_hline(yintercept=0 , color="black", lty="dashed") +
  geom_vline(xintercept=0 , color="black", lty="dashed") 
ggplot(prin_comps, 
       aes(x=PC1, y=PC2,
           color=as.factor(site))) + 
  xlab("PC1 (0.896%)") +
  ylab("PC2 (0.853%)") + 
  geom_point(size=3.3) + 
  geom_point(shape=1, alpha=0.3, color="black", size=3.3) + 
  labs(color="Site") +
  geom_hline(yintercept=0 , color="black", lty="dashed") +
  geom_vline(xintercept=0 , color="black", lty="dashed")



##### Filter for individuals > 6 stdevs away from the mean for the first 10 PCs #####

print("filtering out individuals > 6 stdevs away from the mean for the first 10 PCs...")

for(i in 1:4) {
  pc.means = apply(X=pca$x[,1:10], MARGIN=2, FUN=mean)
  pc.sdevs = apply(X=pca$x[,1:10], MARGIN=2, FUN=sd)

  pc.thresh.up = pc.sdevs*6
  pc.thresh.down = pc.sdevs*-6

  prin_comps_f = data.frame(pca$x)

  samples_to_toss_PC = c()
  for(j in colnames(prin_comps_f[,1:10])) {
    for(i in rownames(prin_comps_f)) {
      #print(prin_comps_f[i,j])
      if(prin_comps_f[i,j] > pc.thresh.up[j]) {
        #print(prin_comps_f[i,j])
        samples_to_toss_PC = c(samples_to_toss_PC, i)
      }
      if(prin_comps_f[i,j] < pc.thresh.down[j]) {
        #print(prin_comps_f[i,j])
        samples_to_toss_PC = c(samples_to_toss_PC, i)
      }
    }
  }

  samples_to_toss_PC
  samples_to_keep_PC = (indv %>% filter(indv %ni% samples_to_toss_PC))$indv
  sample_stats = sample_stats %>% filter(rownames(sample_stats) %in% samples_to_keep_PC)
  SNPs = SNPs[samples_to_keep_PC,]
  indv = indv[samples_to_keep_PC,]
  depth = subset(depth, select=colnames(depth) %ni% samples_to_toss_PC)
  depth_slim = subset(depth_slim, select=colnames(depth_slim) %ni% samples_to_toss_PC)

  u_count = c()
  for(j in 1:length(colnames(SNPs))) {
    u_j = sum(SNPs[,j], na.rm=TRUE) / sum(!is.na(SNPs[,j]))
    u_count = c(u_count, u_j)
  }

  SNPs_c_s = SNPs
  for(j in 1:length(colnames(SNPs))) {
    u_j = u_count[j]
    p_j = u_j / 2
    for(i in 1:length(rownames(SNPs))) {
      c_ij = SNPs_c_s[i,j]
      SNPs_c_s[i,j] = (c_ij - u_j) / (sqrt(p_j * (1 - p_j)))
    }
  }

  SNPs_c_s[is.na(SNPs_c_s)] = 0
  pca = prcomp(x=SNPs_c_s, center=FALSE)
}
 


##### generate more multi-locus PCA plots #####

print("generating more multi-locus PCA plots...")

u_count = c()
for(j in 1:length(colnames(SNPs))) {
  u_j = sum(SNPs[,j], na.rm=TRUE) / sum(!is.na(SNPs[,j]))
  u_count = c(u_count, u_j)
}

SNPs_c_s.2 = SNPs
for(j in 1:length(colnames(SNPs))) {
  u_j = u_count[j]
  p_j = u_j / 2
  for(i in 1:length(rownames(SNPs))) {
    c_ij = SNPs_c_s.2[i,j]
    SNPs_c_s.2[i,j] = (c_ij - u_j) / (sqrt(p_j * (1 - p_j)))
  }
}

SNPs_c_s.2[is.na(SNPs_c_s.2)] = 0

pca.2 = prcomp(x=SNPs_c_s.2, center=FALSE)
eig_sdev.2 = pca.2$sdev
eig_var.2 = eig_sdev.2^2

# percent variance explained
PVE.2 = (eig_var.2/sum(eig_var.2))*100
print(cat("Percent Variance Explained:", 
          round(PVE.2, digits=4), sep="\t"))

# Plot 1st & 2nd PCs
prin_comps.2 = data.frame(pca.2$x)
prin_comps.2 = cbind(indv$subpop, prin_comps.2)
names(prin_comps.2)[1] = "subpop"

for(i in 1:length(prin_comps.2$subpop)) {
  x = prin_comps.2[i,1]
  if(x=="PopN"|x=="PopM"|x=="PopL"|x=="PopK"|x=="PopJ"|x=="PopO") {
    prin_comps.2[i,1] = "West"
  }
  else {
    prin_comps.2[i,1] = "East"
  }
}

prin_comps.2 = cbind(indv$subpop, prin_comps.2)
names(prin_comps.2)[1] = "state"

for(i in 1:length(prin_comps.2$subpop)) {
  x = prin_comps.2[i,1]
  if(x=="PopO") {
    prin_comps.2[i,1] = "Wisconsin"
  }
  if(x=="PopJ") {
    prin_comps.2[i,1] = "Ohio"
  }
  if(x=="PopN") {
    prin_comps.2[i,1] = "Arkansas"
  }
  if(x=="PopM" | x=="PopL" | x=="PopK") {
    prin_comps.2[i,1] = "Louisiana"
  }
  if(x=="PopH") {
    prin_comps.2[i,1] = "South Carolina"
  }
  if(x=="PopG") {
    prin_comps.2[i,1] = "North Carolina"
  }
  if(x=="PopA" | x=="PopB" | x=="PopE" | x=="PopF") {
    prin_comps.2[i,1] = "Virginia"
  }
}

prin_comps.2 = cbind(indv$subpop, prin_comps.2)
names(prin_comps.2)[1] = "site"

for(i in 1:length(prin_comps.2$site)) {
  x = prin_comps.2[i,1]
  if(x=="PopA") {
    prin_comps.2[i,1] = "VA, Deep Bottom"
  }
  if(x=="PopB") {
    prin_comps.2[i,1] = "VA, Dragon Run"
  }
  if(x=="PopE") {
    prin_comps.2[i,1] = "VA, Fort AP Hill"
  }
  if(x=="PopF") {
    prin_comps.2[i,1] = "VA, Great Dismal Swamp"
  }
  if(x=="PopG") {
    prin_comps.2[i,1] = "NC, Holt Lake"
  }
  if(x=="PopH") {
    prin_comps.2[i,1] = "SC, Francis Beidler Forest"
  }
  if(x=="PopJ") {
    prin_comps.2[i,1] = "OH, Hoover Reservoir"
  }
  if(x=="PopK") {
    prin_comps.2[i,1] = "LA, Palmetto Island"
  }
  if(x=="PopL") {
    prin_comps.2[i,1] = "LA, Bluebonnet Swamp"
  }
  if(x=="PopM") {
    prin_comps.2[i,1] = "LA, Barataria Preserve"
  }
  if(x=="PopN") {
    prin_comps.2[i,1] = "AR, White River"
  }
  if(x=="PopO") {
    prin_comps.2[i,1] = "WI, Sugar River"
  }
}

prin_comps.2$subpop = as.factor(prin_comps.2$subpop)
prin_comps.2$state = as.factor(prin_comps.2$state)
prin_comps.2$state = factor(prin_comps.2$state, 
                            levels=c("Virginia", "North Carolina",
                                     "South Carolina", "Louisiana", 
                                     "Arkansas", "Ohio", "Wisconsin"))
prin_comps.2$site = as.factor(prin_comps.2$site)
prin_comps.2$site = factor(prin_comps.2$site, 
                            levels=c("VA, Deep Bottom", "VA, Dragon Run",
                                     "VA, Fort AP Hill", "VA, Great Dismal Swamp",
                                     "NC, Holt Lake", "SC, Francis Beidler Forest",
                                     "LA, Palmetto Island", "LA, Bluebonnet Swamp",
                                     "LA, Barataria Preserve", "AR, White River",
                                     "OH, Hoover Reservoir", "WI, Sugar River"))

pca2.2 = ggplot(prin_comps.2, 
       aes(x=PC1, y=PC2,
           color = as.factor(subpop))) + 
  xlab("PC1 (0.869%)") +
  ylab("PC2 (0.818%)") + 
  geom_point(alpha=0.8, size=3.3) + 
  geom_point(shape=1, alpha=0.3, color="black", size=3.3) +
  labs(color = "Subregion") +
  geom_hline(yintercept=0 , color="black", lty="dashed") +
  geom_vline(xintercept=0 , color="black", lty="dashed") 
pca2.1 = ggplot(prin_comps.2, 
       aes(x=PC1, y=PC2,
           color = as.factor(state))) + 
  xlab("PC1 (0.869%)") +
  ylab("PC2 (0.818%)") + 
  geom_point(alpha=0.8, size=3.3) + 
  geom_point(shape=1, alpha=0.3, color="black", size=3.3) +
  labs(color = "State") + 
  scale_color_manual(values = c("#d53e4f", "#fc8d59", 
                                "#fee08b", "#ffffbf", 
                                "#e6f598", "#99d594",
                                "#3288bd")) +
  geom_hline(yintercept=0 , color="black", lty="dashed") +
  geom_vline(xintercept=0 , color="black", lty="dashed")
ggplot(prin_comps.2, 
       aes(x=PC1, y=PC2,
           color = as.factor(site))) + 
  xlab("PC1 (0.869%)") +
  ylab("PC2 (0.818%)") + 
  geom_point(size=2) + 
  #geom_point(shape=1, alpha=0.2, color="black") + 
  labs(color = "Site") +
  geom_hline(yintercept=0 , color="black", lty="dashed") +
  geom_vline(xintercept=0 , color="black", lty="dashed")

ggarrange(pca2.1, pca2.2, ncol=1, 
          labels=c("A", "B"))



##### preliminary outlier SNP visualization #####

hist(pca.2$rotation[,1])
PC1_weights = pca.2$rotation[,1]
PC1_weight_bounds = quantile(pca.2$rotation[,1], probs=c(0.025,0.975))
PC1_weight_bounds
SNPs_PC1_down = PC1_weights[which(PC1_weights < PC1_weight_bounds[1])] 
SNPs_PC1_up = PC1_weights[which(PC1_weights > PC1_weight_bounds[2])]
SNPs_PC1_down
SNPs_PC1_up



##### calculate and visualizing missing data rate by subregion #####

print("calculating and visualizing missing data rate by subregion...")

sample_stats$subpop = prin_comps.2$subpop

pm_pop_1 = data.frame("percent_missing"=sample_stats$percent_missing[which(sample_stats$subpop=="West")])
pm_pop_2 = data.frame("percent_missing"=sample_stats$percent_missing[which(sample_stats$subpop=="East")])

ggplot(pm_pop_1, aes(percent_missing)) + geom_area(stat="bin") + 
  xlab("Percent Missing")
ggplot(pm_pop_2, aes(percent_missing)) + geom_area(stat="bin") + 
  xlab("Percent Missing")

print(paste("Pop 1 percent missing: ", mean(pm_pop_1$percent_missing), sep=""))
print(paste("Pop 2 percent missing: ", mean(pm_pop_2$percent_missing), sep=""))



##### hierarchical F stats #####

levels = data.frame(prin_comps.2$subpop, indv$subpop)
colnames(levels) = c("subregion", "site")
levels$subregion = as.factor(levels$subregion)
levels$site = as.factor(levels$site)
varcomp.out = varcomp.glob(levels=levels, loci=SNPs, diploid=TRUE)

F.varcomp = varcomp.out$F
loc.varcomp = data.frame(varcomp.out$loc)

Fct = loc.varcomp[1] / (loc.varcomp[1] + loc.varcomp[2] + loc.varcomp[3] + loc.varcomp[4])
Fst = (loc.varcomp[1] + loc.varcomp[2]) / (loc.varcomp[1] + loc.varcomp[2] +
                                             loc.varcomp[3] + loc.varcomp[4])
Fsc = loc.varcomp[2] / (loc.varcomp[2] + loc.varcomp[3] + loc.varcomp[4])
F.Eck = loc.varcomp[1] / (loc.varcomp[1] + loc.varcomp[2])

SNP_F_stats = data.frame(SNP_stats[,1:2], Fct, Fst, Fsc, F.Eck)
colnames(SNP_F_stats) = c("contig", "position", "Fct", "Fst", "Fsc", "F.Eck")

write.table(F.varcomp, "hier.F.stats.txt", col.names=TRUE, row.names=TRUE, sep="\t")



##### create structure table #####

print("creating structure table...")

SNPs[is.na(SNPs)] = -1

struct = data.frame(rownames(SNPs), prin_comps.2$site)

struct_SNPs = c()
for(i in colnames(SNPs)) {
  cur_SNP = i
  SNP_1 = paste(i, ".1", sep="")
  SNP_2 = paste(i, ".2", sep="")
  struct_SNPs = c(struct_SNPs, SNP_1, SNP_2)
}

for(i in struct_SNPs) {
  struct[[i]] = NA
}

rownames(struct) = rownames(SNPs)

for(i in colnames(SNPs)) {
  cur_SNP = i
  SNP_1 = paste(i, ".1", sep="")
  SNP_2 = paste(i, ".2", sep="")
  for(j in rownames(SNPs)) {
    #print(SNPs[j,i])
    if(SNPs[j,i] == -1) {
      struct[j, SNP_1] = -9
      struct[j, SNP_2] = -9
    }
    if(SNPs[j,i] == 0) {
      struct[j, SNP_1] = 0
      struct[j, SNP_2] = 0
    }
    if(SNPs[j,i] == 1) {
      struct[j, SNP_1] = 0
      struct[j, SNP_2] = 1
    }
    if(SNPs[j,i] == 2) {
      struct[j, SNP_1] = 1
      struct[j, SNP_2] = 1
    }
  }
}

names(struct)[1] = "individual"
names(struct)[2] = "site"



##### write .012 and structure output files #####

print("writing .012 and structure output files...")

struct_out = struct
indv_out = data.frame(indv, prin_comps.2$subpop)
pos_out = data.frame(SNP_stats$contig, SNP_stats$position)
colnames(pos_out) = NULL
colnames(indv_out) = NULL
SNPs_out = SNPs
rownames(SNPs_out) = NULL
colnames(SNPs_out) = NULL
colnames(struct_out) = NULL
rownames(struct_out) = NULL

# indv
write.table(indv_out, indv_out_fh, quote=FALSE, sep="\t", row.names=FALSE)

# pos
write.table(pos_out, pos_out_fh, quote=FALSE, sep="\t", row.names=FALSE)

# SNP_names
write.table(colnames(SNPs), SNP_names_out_fh, quote=FALSE, sep="\t", row.names=FALSE)

# SNPs
write.table(SNPs_out, SNPs_out_fh, quote=FALSE, sep="\t", row.names=FALSE)

#contigs
write.table(pos$contig, "dDocent_SNP_contig_ids.txt", quote=FALSE, sep="\t", row.names=FALSE)

# structure
write.table(struct_out, struct_out_fh, quote=FALSE, sep=" ", row.names=FALSE)



##### create .ped table #####

print("creating .ped table...")

fam_id = seq(1, length(rownames(SNPs)))
indv_id = rownames(SNPs)
pat_id = rep(0, length(rownames(SNPs)))
mat_id = pat_id
sex = pat_id
phenotype = pat_id

ped = struct[,3:ncol(struct)]
ped[ped == 1] = 2
ped[ped == 0] = 1
ped[ped == -9] = 0 

ped = cbind(fam_id, indv_id, pat_id, mat_id, sex, phenotype, ped)

rownames(ped) = rownames(SNPs)



##### create .map table #####

print("creating .map table...")

chrom_map = c() 
for(i in SNP_stats$contig) {
  chrom_map = c(chrom_map, 1)
}

SNP_id = colnames(SNPs)
Pos_morgans = rep(0, length(colnames(SNPs)))

BP_coord = c()
n = 1
for(i in 1:length(SNP_stats$contig)) {
  BP_coord = c(BP_coord, n)
  n = n + 10000
}

map = data.frame(chrom_map, SNP_id, Pos_morgans, BP_coord)



##### write .ped and .map output files #####

print("writing .ped and .map output files")

ped_out = ped
map_out = map

colnames(ped_out) = NULL
rownames(ped_out) = NULL
colnames(map_out) = NULL

# ped
write.table(ped_out, ped_out_fh, quote=FALSE, sep="\t", row.names=FALSE)

# map
write.table(map_out, map_out_fh, quote=FALSE, sep="\t", row.names=FALSE)

# convert to .bed in PLINK



##### pcadapt #####

print("identifying outlier SNPs through pcadapt...")

pcadapt.filename = read.pcadapt("ddocent_60_filtered.bed", type = "bed")
outlier.object = pcadapt(input=pcadapt.filename, K=100)

plot(outlier.object, option = "screeplot")

sitelist.int = as.numeric(prin_comps.2$site)
sitelist.names = prin_comps.2$site

statelist.int = as.numeric(prin_comps.2$state)
statelist.names = prin_comps.2$state

poplist.int = as.numeric(prin_comps.2$subpop)
poplist.names = prin_comps.2$subpop

plot(outlier.object, option = "scores", pop = sitelist.names)
plot(outlier.object, option = "scores", pop = statelist.names)
plot(outlier.object, option = "scores", pop = poplist.names)

plot(outlier.object, option = "scores", i = 3, j = 4, pop = sitelist.names)
plot(outlier.object, option = "scores", i = 3, j = 4, pop = statelist.names)
plot(outlier.object, option = "scores", i = 3, j = 4, pop = poplist.names)

plot(outlier.object, option = "scores", i = 4, j = 5, pop = sitelist.names)
plot(outlier.object, option = "scores", i = 4, j = 5, pop = statelist.names)
plot(outlier.object, option = "scores", i = 4, j = 5, pop = poplist.names)

pcadapt_scores = data.frame(indv$indv, prin_comps.2$site, prin_comps.2$state, 
                            prin_comps.2$subpop, outlier.object$scores[,1:5])
colnames(pcadapt_scores) = c("individual", "site", "state", "subregion",
                             "PC1", "PC2", "PC3", "PC4", "PC5")
ggplot(pcadapt_scores, 
       aes(x=PC1, y=PC3,
           color = as.factor(site),
           label=individual)) +
  geom_text(hjust=0, vjust=0) + 
  xlab("PC1") +
  ylab("PC2") + 
  geom_point(alpha=0.8, size=2) + 
  geom_point(shape = 1, alpha=0.3, color = "black", size=2) +
  labs(color = "Subregion") +
  geom_hline(yintercept= 0 , color = "black", lty = "dashed") +
  geom_vline(xintercept= 0 , color = "black", lty = "dashed") 

p1 = ggplot(pcadapt_scores, 
            aes(x=PC1, y=PC2,
                color = as.factor(state))) + 
  xlab("PC1") +
  ylab("PC2") + 
  ggtitle("Projection onto PC1 and PC2") +  
  scale_color_manual(values = c("#d53e4f", "#fc8d59", 
                                "#fee08b", "#ffffbf", 
                                "#e6f598", "#99d594",
                                "#3288bd")) +
  geom_point(alpha=0.8, size=2) + 
  geom_point(shape = 1, alpha=0.3, color = "black", size=2) +
  labs(color = "Subregion") +
  geom_hline(yintercept= 0 , color = "black", lty = "dashed") +
  geom_vline(xintercept= 0 , color = "black", lty = "dashed") +
  theme(legend.position="none")
p2 = ggplot(pcadapt_scores, 
            aes(x=PC3, y=PC2,
                color = as.factor(state))) + 
  xlab("PC3") +
  ylab("PC2") +
  ggtitle("Projection onto PC2 and PC3") +  
  scale_color_manual(values = c("#d53e4f", "#fc8d59", 
                                "#fee08b", "#ffffbf", 
                                "#e6f598", "#99d594",
                                "#3288bd")) +
  geom_point(alpha=0.8, size=2) + 
  geom_point(shape = 1, alpha=0.3, color = "black", size=2) +
  labs(color = "Subregion") +
  geom_hline(yintercept= 0 , color = "black", lty = "dashed") +
  geom_vline(xintercept= 0 , color = "black", lty = "dashed")
p3 = ggplot(pcadapt_scores, 
            aes(x=PC3, y=PC4,
                color = as.factor(state))) + 
  xlab("PC3") +
  ylab("PC4") +
  ggtitle("Projection onto PC3 and PC4") +  
  scale_color_manual(values = c("#d53e4f", "#fc8d59", 
                                "#fee08b", "#ffffbf", 
                                "#e6f598", "#99d594",
                                "#3288bd")) +
  geom_point(alpha=0.8, size=2) + 
  geom_point(shape = 1, alpha=0.3, color = "black", size=2) +
  labs(color = "Subregion") +
  geom_hline(yintercept= 0 , color = "black", lty = "dashed") +
  geom_vline(xintercept= 0 , color = "black", lty = "dashed") +
  theme(legend.position="none")
p4 = ggplot(pcadapt_scores, 
            aes(x=PC5, y=PC4,
                color = as.factor(state))) + 
  xlab("PC5") +
  ylab("PC4") +
  ggtitle("Projection onto PC4 and PC5") +  
  scale_color_manual(values = c("#d53e4f", "#fc8d59", 
                                "#fee08b", "#ffffbf", 
                                "#e6f598", "#99d594",
                                "#3288bd")) +
  geom_point(alpha=0.8, size=2) + 
  geom_point(shape = 1, alpha=0.3, color = "black", size=2) +
  labs(color = "Subregion") +
  geom_hline(yintercept= 0 , color = "black", lty = "dashed") +
  geom_vline(xintercept= 0 , color = "black", lty = "dashed") +
  theme(legend.position="none")
ggarrange(p1, p2, p3, p4, ncol=2, 
          labels=c("A", "B", "C", "D"))


outlier.object = pcadapt(input=pcadapt.filename, K = 5)
summary(outlier.object)

pvals = data.frame(outlier.object$pvalues)
colnames(pvals) = c("pvalues")

p1 = plot(outlier.object, option="manhattan")
p2 = plot(outlier.object, option = "qqplot")
p3 = ggplot(pvals, aes(pvalues)) + 
  geom_area(stat="bin") + 
  ggtitle("Histogram") +
  xlab("p-values")
#hist(outlier.object$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
p4 = plot(outlier.object, option = "stat.distribution")
ggarrange(p1, p2, p3, p4, ncol=2, 
          labels=c("A", "B", "C", "D"))


alpha = 0.05

#benjamini-hochberg procedure
padj_BH = p.adjust(outlier.object$pvalues,method="BH")
outliers_BH = which(padj_BH < alpha)
length(outliers_BH)

#bonferroni correction
padj_BC = p.adjust(outlier.object$pvalues,method="bonferroni")
outliers_BC = which(padj_BC < alpha)
length(outliers_BC)

#associations between PCs and outliers
snp_pc_BH = data.frame(get.pc(outlier.object, outliers_BH))
snp_pc_BC = data.frame(get.pc(outlier.object, outliers_BC))

pos$outlier_BH = NA
pos$outlier_BC = NA
for(i in 1:length(rownames(pos))) {
  #print(i)
  if(i %in% outliers_BC) {
    pos$outlier_BC[i] = 1
  }
  else {
    pos$outlier_BC[i] = 0
  }
  if(i %in% outliers_BH) {
    pos$outlier_BH[i] = 1
  }
  else {
    pos$outlier_BH[i] = 0
  }
}

pos.outlier.BH = pos[which(pos$outlier_BH == 1), 1:2] 
pos.outlier.BC = pos[which(pos$outlier_BC == 1), 1:2]

intersect(pos.outlier.BH$contig, pos.outlier.BC$contig)

# investigate weird outlier individuals on PC1
PopJ_004_SNPs = colnames(SNPs[which(SNPs["PopJ_004",] != -1)])
PopJ_017_SNPs = colnames(SNPs[which(SNPs["PopJ_017",] != -1)])

PopJ_008_SNPs = colnames(SNPs[which(SNPs["PopJ_008",] != -1)])
PopA_022_SNPs = colnames(SNPs[which(SNPs["PopA_022",] != -1)])

bonferroni_pc1 = rownames(pos[which(pos$outlier_BC == 1),])
BH_pc1 = rownames(pos[which(pos$outlier_BH == 1),])

# outlier individuals
PopJ_004_BC = intersect(PopJ_004_SNPs, bonferroni_pc1)
PopJ_004_BH = intersect(PopJ_004_SNPs, BH_pc1)
PopJ_017_BC = intersect(PopJ_017_SNPs, bonferroni_pc1)
PopJ_017_BH = intersect(PopJ_017_SNPs, BH_pc1)

# compared to non-outlier individuals
PopJ_008_BC = intersect(PopJ_008_SNPs, bonferroni_pc1)
PopJ_008_BH = intersect(PopJ_008_SNPs, BH_pc1)
PopA_022_BC = intersect(PopA_022_SNPs, bonferroni_pc1)
PopA_022_BH = intersect(PopA_022_SNPs, BH_pc1)

# compare missing data between J004 and J017
missing_j004 = colnames(SNPs["PopJ_004", which(SNPs["PopJ_004",] == -1)])
missing_j017 = colnames(SNPs["PopJ_017", which(SNPs["PopJ_017",] == -1)])
missing_h007 = colnames(SNPs["PopH_007", which(SNPs["PopH_007",] == -1)])

missing_j_overlap = intersect(missing_j004, missing_j017)
missing_test = intersect(missing_j017, missing_h007)

SNPs_j004_j017 = rbind(SNPs["PopJ_004",], SNPs["PopJ_017",])

#top p-value SNPs
SNP_P_BH = data.frame(colnames(SNPs), pos$contig, padj_BH)
SNP_P_BC = data.frame(colnames(SNPs), pos$contig, padj_BC)
colnames(SNP_P_BH) = c("SNP", "contig", "P") 
colnames(SNP_P_BC) = c("SNP", "contig", "P")

SNP_P_BH = SNP_P_BH[with(SNP_P_BH,order(P)),]
SNP_P_BC = SNP_P_BC[with(SNP_P_BC,order(P)),]

top_10_BH = SNP_P_BH[1:10,]
top_10_BC = SNP_P_BC[1:10,]



##### RDA #####

print("conducting RDA to identify outlier SNPs associated with environment...")

# pull environmental data from worldclim
wld_cli = getData("worldclim", var="bio", res=2.5)
points = data.frame(lat_long[,5], lat_long[,4])
rownames(points) = rownames(lat_long)
colnames(points) = c("x", "y")
clim_data = extract(wld_cli, points)
rownames(clim_data) = rownames(lat_long)

env_pred = data.frame()
geo_condition = data.frame()
for(i in prin_comps.2$site) {
  env_pred = rbind(env_pred, clim_data[i,])
  geo_condition = rbind(geo_condition, points[i,])
}
colnames(env_pred) = colnames(clim_data)
rownames(env_pred) = rownames(SNPs_c_s.2)
colnames(geo_condition) = c("long", "lat")
rownames(geo_condition) = rownames(SNPs_c_s.2)

# first rda
PROW.rda = rda(SNPs_c_s.2, env_pred, geo_condition, scale.=TRUE) 
PROW.rda

RsquareAdj(PROW.rda)
summary(eigenvals(PROW.rda, model = "constrained"))
screeplot(PROW.rda, main="")
signif.full = anova.cca(PROW.rda, parallel=getOption("mc.cores"))
signif.full
signif.axis = anova.cca(PROW.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis
vif.cca(PROW.rda)

sub = prin_comps.2$state
#sub = prin_comps.2$site

bg = c("#d53e4f", "#fc8d59", 
       "#fee08b", "#ffffbf", 
       "#e6f598", "#99d594",
       "#3288bd")
#bg = c("#a6cee3", "#1f78b4",
      #"#b2df8a", "#33a02c",
      #"#fb9a99", "#e31a1c",
      #"#fdbf6f", "#ff7f00",
      #"#cab2d6", "#6a3d9a",
      #"#ffff99", "#b15928")

plot(PROW.rda, type="n", scaling=1)
points(PROW.rda, display="species", pch=3, cex=0.05, col="gray", scaling=1)  
points(PROW.rda, display="sites", pch=21, cex=1.2, col="gray32", scaling=1, bg=bg[sub]) 
text(PROW.rda, scaling=1, display="bp", col="black", cex=1.5)

# PCA of environmental predictors
pc.env = prcomp(x=env_pred, center=TRUE, scale.=TRUE)

env.prin_comps = data.frame(pc.env$x)

env_sdev = pc.env$sdev
env_var = env_sdev^2

# env percent variance explained
env_PVE = (env_var/sum(env_var))*100
sum(env_PVE[1:2])
sum(env_PVE[1:3])

env.pc.rotation = data.frame(pc.env$rotation[,1:3])

par(mfrow=c(3,1))
barplot(pc.env$rotation[,1], 
        main="Environmental Variable Loadings: PC1")
barplot(pc.env$rotation[,2], 
        main="Environmental Variable Loadings: PC2")
barplot(pc.env$rotation[,3], 
        main="Environmental Variable Loadings: PC3")
par(mfrow=c(1,1))

env_pcs = env.prin_comps[,1:3]
colnames(env_pcs) = c("ENV_PC1", "ENV_PC2", "ENV_PC3")

# second rda
PROW.pc.rda = rda(SNPs_c_s.2, env_pcs, geo_condition, scale=T) 
PROW.pc.rda

RsquareAdj(PROW.pc.rda)
summary(eigenvals(PROW.pc.rda, model = "constrained"))
screeplot(PROW.pc.rda, main="")
signif.full.pc = anova.cca(PROW.pc.rda, parallel=getOption("mc.cores"))
signif.full.pc
signif.axis.pc = anova.cca(PROW.pc.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis.pc
vif.cca(PROW.pc.rda)

plot(PROW.pc.rda, type="n", scaling=1)
points(PROW.pc.rda, display="species", pch=3, cex=0.05, col="gray", scaling=1)  
points(PROW.pc.rda, display="sites", pch=21, cex=1.2, col="gray32", scaling=1, bg=bg[sub]) 
text(PROW.pc.rda, scaling=1, display="bp", col="black", cex=1)

plot(PROW.pc.rda, type="n", scaling=1, choices=c(1,3))
points(PROW.pc.rda, display="species", pch=3, cex=0.05, col="gray", scaling=1, choices=c(1,3))  
points(PROW.pc.rda, display="sites", pch=21, cex=1.2, col="gray32", scaling=1, bg=bg[sub], choices=c(1,3)) 
text(PROW.pc.rda, scaling=1, display="bp", col="black", cex=1, choices=c(1,3))

plot(PROW.pc.rda, type="n", scaling=1, choices=c(2,3))
points(PROW.pc.rda, display="species", pch=3, cex=0.05, col="gray", scaling=1, choices=c(2,3))  
points(PROW.pc.rda, display="sites", pch=21, cex=1.2, col="gray32", scaling=1, bg=bg[sub], choices=c(2,3)) 
text(PROW.pc.rda, scaling=1, display="bp", col="black", cex=1, choices=c(2,3))



##### pairwise fst #####

print("calculating pairwise Fst and plotting against pairwise geographic distance...")

pairwise_SNP_table = data.frame(prin_comps.2$site, SNPs)
names(pairwise_SNP_table)[1] = "site"

pairwise_fst = pairwise.WCfst(pairwise_SNP_table, diploid=TRUE)

site_1 = c()
site_2 = c()
site_old = c()
  
for(i in rownames(lat_long)) {
  site_old = c(site_old, i)
  for(j in rownames(lat_long)) {
    if(i != j) {
      if(j %ni% site_old) {
        #print(paste(i, j, sep="-"))
        site_1 = c(site_1, i)
        site_2 = c(site_2, j)
      }
    }
  }
}

site_stats = data.frame(site_1, site_2)

fst_pair = c()
for(i in 1:length(site_stats$site_1)) {
  #print(paste(site_stats$site_1[i], site_stats$site_2[i]))
  x = pairwise_fst[site_stats$site_1[i], site_stats$site_2[i]]
  #print(x)
  fst_pair = c(fst_pair, x)
}

distance = c()
for(i in 1:length(site_stats$site_1)) {
  #print(paste(site_stats$site_1[i], site_stats$site_2[i]))
  geo_dist = distVincentyEllipsoid(points[site_stats$site_1[i],], points[site_stats$site_2[i],]) / 1000
  #print(geo_dist)
  distance = c(distance, geo_dist)
}

cc = as.matrix(points)
dist_mat = distm(cc, fun=distVincentyEllipsoid)

site_stats = data.frame(site_stats, fst_pair, distance)
names(site_stats)[3] = "fst"
site_stats$`fst / 1-fst` = site_stats$fst / (1-site_stats$fst)
site_stats$`ln fst` = log(site_stats$fst)
site_stats$log10_dist = log10(site_stats$distance)


ggplot(site_stats, 
       aes(x=distance, y=`fst / 1-fst`)) + 
  xlab("Distance (km)") +
  ylab("Fst / 1-Fst") + 
  geom_point(alpha=0.5) + 
  geom_point(shape=1) +
  geom_smooth(method=lm, se=FALSE)

ggplot(site_stats, 
       aes(x=log10_dist, y=`fst / 1-fst`)) + 
  xlab("log10(Geographical Distance (km))") +
  ylab("Fst / 1-Fst") + 
  geom_point(alpha=0.5, size=3.3) + 
  geom_point(shape=1, size=3.3) +
  geom_smooth(method=lm, se=TRUE)

pairwise_gen_dist = pairwise_fst / (1-pairwise_fst)

mantel(log10(dist_mat), pairwise_gen_dist, method="pearson", permutations=1000)

write.table(pairwise_fst, "pairwise_fst_table.txt", col.names=TRUE, row.names=TRUE, sep="\t")



##### third RDA #####

print("conducting RDA for breeding season environmental predictors...")

env.rast.list = list.files(path=paste(wd, "/RDA_predictors", sep=""), pattern=".tif$", all.files=TRUE, full.names=TRUE)
env.rasters = stack(env.rast.list)
breeding.env = data.frame(extract(env.rasters, points))
rownames(breeding.env) = rownames(lat_long)
colnames(breeding.env) = c('prec_4', 'prec_5', 'prec_6', 'prec_7',
                           'srad_4', 'srad_5', 'srad_6', 'srad_7',
                           'tmax_4', 'tmax_5', 'tmax_6', 'tmax_7',
                           'tmin_4', 'tmin_5', 'tmin_6', 'tmin_7',
                           'vapr_4', 'vapr_5', 'vapr_6', 'vapr_7')

env_pred_2 = data.frame()
for(i in prin_comps.2$site) {
  env_pred_2 = rbind(env_pred_2, breeding.env[i,])
}
rownames(env_pred_2) = rownames(SNPs_c_s.2)

# PCA of environmental predictors
pc.env.2 = prcomp(x=env_pred_2, center=TRUE, scale.=TRUE)

env.prin_comps.2 = data.frame(pc.env.2$x)

env_sdev.2 = pc.env.2$sdev
env_var.2 = env_sdev.2^2

# env percent variance explained
env_PVE.2 = (env_var.2/sum(env_var.2))*100
sum(env_PVE.2[1:2])
sum(env_PVE.2[1:3])

env.pc.rotation.2 = data.frame(pc.env.2$rotation[,1:3])

par(mfrow=c(3,1))
barplot(pc.env.2$rotation[,1], 
        main="Environmental Variable Loadings: PC1")
barplot(pc.env.2$rotation[,2], 
        main="Environmental Variable Loadings: PC2")
barplot(pc.env.2$rotation[,3], 
        main="Environmental Variable Loadings: PC3")
par(mfrow=c(1,1))

env_pcs.2 = env.prin_comps.2[,1:3]
colnames(env_pcs.2) = c("ENV_PC1", "ENV_PC2", "ENV_PC3")

PROW.pc.rda.2 = rda(SNPs_c_s.2, env_pcs.2, geo_condition, scale=T) 
PROW.pc.rda.2

RsquareAdj(PROW.pc.rda.2)
summary(eigenvals(PROW.pc.rda.2, model = "constrained"))
screeplot(PROW.pc.rda.2, main="")
signif.full.pc.2 = anova.cca(PROW.pc.rda.2, parallel=getOption("mc.cores"))
signif.full.pc.2
signif.axis.pc.2 = anova.cca(PROW.pc.rda.2, by="axis", parallel=getOption("mc.cores"))
signif.axis.pc.2
vif.cca(PROW.pc.rda.2)

plot(PROW.pc.rda.2, type="n", scaling=1)
points(PROW.pc.rda.2, display="species", pch=3, cex=0.05, col="gray", scaling=1)  
points(PROW.pc.rda.2, display="sites", pch=21, cex=1.5, col="gray32", scaling=1, bg=bg[sub]) 
text(PROW.pc.rda.2, scaling=1, display="bp", col="black", cex=1)

plot(PROW.pc.rda.2, type="n", scaling=1, choices=c(1,3))
points(PROW.pc.rda.2, display="species", pch=3, cex=0.05, col="gray", scaling=1, choices=c(1,3))  
points(PROW.pc.rda.2, display="sites", pch=21, cex=1.5, col="gray32", scaling=1, bg=bg[sub], choices=c(1,3)) 
text(PROW.pc.rda.2, scaling=1, display="bp", col="black", cex=1, choices=c(1,3))

plot(PROW.pc.rda.2, type="n", scaling=1, choices=c(2,3))
points(PROW.pc.rda.2, display="species", pch=3, cex=0.05, col="gray", scaling=1, choices=c(2,3))  
points(PROW.pc.rda.2, display="sites", pch=21, cex=1.5, col="gray32", scaling=1, bg=bg[sub], choices=c(2,3)) 
text(PROW.pc.rda.2, scaling=1, display="bp", col="black", cex=1, choices=c(2,3))

# identify outlier snps
load.rda = scores(PROW.pc.rda.2, choices=c(1:3), display="species") 

hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3")

outliers = function(x, z) {
  lims = mean(x) + c(-1, 1) * z * sd(x)    
  x[x < lims[1] | x > lims[2]]
}

cand1 = outliers(load.rda[,1],3)
cand2 = outliers(load.rda[,2],3)
cand3 = outliers(load.rda[,3],3)
ncand = length(cand1) + length(cand2) + length(cand3)
ncand

cand1 = cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 = cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 = cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3)<- c("axis","snp","loading")

cand = rbind(cand1, cand2, cand3)
cand$snp = as.character(cand$snp)

foo = matrix(nrow=(ncand), ncol=3)
colnames(foo) = colnames(env.prin_comps[,1:3])

for (i in 1:length(cand$snp)) {
  #print(paste("i", i, sep=" "))
  nam = cand[i,2]
  #print(paste("nam", nam, sep=" "))
  snp.gen = SNPs_c_s.2[,nam]
  #print(paste("snp.gen", snp.gen, sep=" "))
  foo[i,] = apply(env.prin_comps[,1:3], 2, function(x) cor(x, snp.gen))
  #print(foo[i,])
}

cand = cbind.data.frame(cand,foo)  
head(cand)

colnames(cand) = c("RDA axis", "SNP", "Loading", colnames(env_pcs.2))

cand_overlap = cand[which(cand$SNP %in% rownames(pos.outlier.BC)), "SNP"]
cand_overlap
cand_overlap_loadings = cand[which(cand$SNP %in% rownames(pos.outlier.BC)), "Loading"]
cand_overlap_loadings

pos$outlier_rda = NA
pos$outlier_overlap = NA
for(i in rownames(pos)) {
  #print(i)
  if(i %in% cand$SNP) {
    pos[i, 'outlier_rda'] = 1
  }
  else {
    pos[i, 'outlier_rda'] = 0
  }
  if(i %in% cand_overlap) {
    pos[i, 'outlier_overlap'] = 1
  }
  else {
    pos[i, 'outlier_overlap'] = 0
  }
}

pos.outlier.rda = pos[which(pos$outlier_rda == 1), 1:2] 
pos.outlier.overlap = pos[which(pos$outlier_overlap == 1), 1:2]

cand$abs_loading = abs(cand$Loading)
cand_sorted = cand[with(cand,order(-abs_loading)),]
rda_top_10 = pos[cand_sorted[1:10, "SNP"], "contig"]



##### write outlier contig files #####

print("writing outlier contig files...")

write.table(pos.outlier.BH$contig, "pcadapt.BH.outlier.contigs.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(pos.outlier.BC$contig, "pcadapt.BC.outlier.contigs.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(pos.outlier.rda$contig, "rda.outlier.contigs.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(pos.outlier.overlap$contig, "overlap.outlier.contigs.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(top_10_BH$contig, "pcadapt.BH.top10.contigs.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(top_10_BC$contig, "pcadapt.BC.top10.contigs.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(rda_top_10, "rda.top10.contigs.txt", quote=FALSE, sep="\t", row.names=FALSE)



##### read in .bam alignment file #####

print("reading in .bam alignment file...")

what = c("rname", "strand", "pos", "seq", "qual", "maps")
param = ScanBamParam(what=scanBamWhat(), flag=scanBamFlag(isUnmappedQuery=FALSE))
bam = scanBam("ddocent_mywa_contigs.bam", param=param)

bam.pos.info = data.frame(bam[[1]]$qname, bam[[1]]$rname, bam[[1]]$pos)
colnames(bam.pos.info) = c("contig", "chromosome", "contig_pos")

bam.pos.info = bam.pos.info[which(bam.pos.info$contig %in% pos$contig),]
bam.SNPs = c()
for(i in bam.pos.info$contig) {
  #print(i)
  bam.SNPs = c(bam.SNPs, rownames(pos)[pos$contig == i])
}
rownames(bam.pos.info) = bam.SNPs
bam.SNP.pos = c()
for(i in rownames(bam.pos.info)) {
  #print(i)
  bam.SNP.pos = c(bam.SNP.pos, pos[i, "position"])
}
bam.pos.info$SNP_pos = bam.SNP.pos
bam.pos.info$chrom_pos = bam.pos.info$contig_pos + 
  (bam.pos.info$SNP_pos - 1)

pos.outlier.rda = data.frame(pos.outlier.rda, 
                             bam.pos.info[intersect(rownames(bam.pos.info), rownames(pos.outlier.rda)), "chromosome"],
                             bam.pos.info[intersect(rownames(bam.pos.info), rownames(pos.outlier.rda)), "chrom_pos"])
colnames(pos.outlier.rda) = c("contig", "pos", "chromosome", "chrom_pos")
pos.outlier.BH = data.frame(pos.outlier.BH, 
                             bam.pos.info[intersect(rownames(bam.pos.info), rownames(pos.outlier.BH)), "chromosome"],
                             bam.pos.info[intersect(rownames(bam.pos.info), rownames(pos.outlier.BH)), "chrom_pos"])
colnames(pos.outlier.BH) = c("contig", "pos", "chromosome", "chrom_pos")
pos.outlier.BC = data.frame(pos.outlier.BC, 
                             bam.pos.info[intersect(rownames(bam.pos.info), rownames(pos.outlier.BC)), "chromosome"],
                             bam.pos.info[intersect(rownames(bam.pos.info), rownames(pos.outlier.BC)), "chrom_pos"])
colnames(pos.outlier.BC) = c("contig", "pos", "chromosome", "chrom_pos")
pos.outlier.overlap = data.frame(rownames(pos.outlier.overlap), pos.outlier.overlap$contig, 
                             SNP_P_BC[which(SNP_P_BC$SNP %in% rownames(pos.outlier.overlap)), "P"],
                             cand[which(cand$SNP %in% rownames(pos.outlier.overlap)), c("RDA axis", "Loading")],
                             bam.pos.info[intersect(rownames(bam.pos.info), rownames(pos.outlier.overlap)), "chromosome"],
                             bam.pos.info[intersect(rownames(bam.pos.info), rownames(pos.outlier.overlap)), "chrom_pos"])
colnames(pos.outlier.overlap) = c("SNP", "contig", "P", "RDA axis", "loading", "chromosome", "chrom_pos")
top_10_BC = data.frame(top_10_BC, bam.pos.info[intersect(rownames(bam.pos.info), top_10_BC$SNP), "chromosome"],
                       bam.pos.info[intersect(rownames(bam.pos.info), top_10_BC$SNP), "chrom_pos"])
colnames(top_10_BC) = c("SNP", "contig", "P", "chromosome", "chrom_pos")
top_10_BH = data.frame(top_10_BH, bam.pos.info[intersect(rownames(bam.pos.info), top_10_BH$SNP), "chromosome"],
                       bam.pos.info[intersect(rownames(bam.pos.info), top_10_BH$SNP), "chrom_pos"])
colnames(top_10_BH) = c("SNP", "contig", "P", "chromosome", "chrom_pos")
pos.rda.top_10 = data.frame(cand_sorted[1:10, "SNP"], 
                            bam.pos.info[intersect(rownames(bam.pos.info), cand_sorted[1:10, "SNP"]), "contig"], 
                            cand_sorted[1:10, c("RDA axis", "Loading")],
                            bam.pos.info[intersect(rownames(bam.pos.info), cand_sorted[1:10, "SNP"]), "chromosome"],
                            bam.pos.info[intersect(rownames(bam.pos.info), cand_sorted[1:10, "SNP"]), "chrom_pos"])
colnames(pos.rda.top_10) = c("SNP", "contig", "RDA axis", "loading", "chromosome", "chrom_pos")



##### read in and parse genbank file #####

print("reading in and parsing genbank file...")

gb = readGFF("mywagenomev2.1.all.noseq.gff")

genes = gb@listData[["Name"]][which(gb@listData[["type"]]=="gene")]
gene_chrom = gb@listData[["seqid"]][which(gb@listData[["type"]]=="gene")]
gene_chrom = as.character(gene_chrom)
gene_start = gb@listData[["start"]][which(gb@listData[["type"]]=="gene")]
gene_end = gb@listData[["end"]][which(gb@listData[["type"]]=="gene")]

gene = data.frame(genes, gene_chrom, gene_start, gene_end)

proteins = gb@listData[["Name"]][which(gb@listData[["type"]]=="protein_match")]
prot_chrom = gb@listData[["seqid"]][which(gb@listData[["type"]]=="protein_match")]
prot_chrom = as.character(prot_chrom)
prot_start = gb@listData[["start"]][which(gb@listData[["type"]]=="protein_match")]
prot_end = gb@listData[["end"]][which(gb@listData[["type"]]=="protein_match")]

prot = data.frame(proteins, prot_chrom, prot_start, prot_end)

chrom_conversion_A = c('chr1',
                       'chr1a',
                       'chr2',
                       'chr3',
                       'chr4',
                       'chr5',
                       'chr4a',
                       'chr6',
                       'chr7',
                       'chr8',
                       'chr9',
                       'chr10',
                       'chr11',
                       'chr12',
                       'chr13',
                       'chr14',
                       'chr15',
                       'chr17',
                       'chr18',
                       'chr19',
                       'chr20',
                       'chr21',
                       'chr22',
                       'chr23',
                       'chr24',
                       'chr25',
                       'chr26',
                       'chr27',
                       'chr28',
                       'chr29',
                       'chrz',
                       'mito')
chrom_conversion_B = c('CM027507.1',
                       'CM027536.1',
                       'CM027508.1',
                       'CM027509.1',
                       'CM027510.1',
                       'CM027511.1',
                       'CM027537.1',
                       'CM027512.1',
                       'CM027513.1',
                       'CM027514.1',
                       'CM027515.1',
                       'CM027516.1',
                       'CM027517.1',
                       'CM027518.1',
                       'CM027519.1',
                       'CM027520.1',
                       'CM027521.1',
                       'CM027522.1',
                       'CM027523.1',
                       'CM027524.1',
                       'CM027525.1',
                       'CM027526.1',
                       'CM027527.1',
                       'CM027528.1',
                       'CM027529.1',
                       'CM027530.1',
                       'CM027531.1',
                       'CM027532.1',
                       'CM027533.1',
                       'CM027534.1',
                       'CM027535.1',
                       'CM027538.1')

chrom_conversion = data.frame(chrom_conversion_A, chrom_conversion_B)
colnames(chrom_conversion) = c("A", "B")

for(i in 1:length(chrom_conversion$A)) {
  A = chrom_conversion[i, 1]
  B = chrom_conversion[i, 2]
  gene$gene_chrom[gene$gene_chrom == A] = B
  prot$prot_chrom[prot$prot_chrom == A] = B
}



##### annotate outlier SNPs #####

print("annotating outlier SNPs...")

gene_match_SNP = c()
gene_match_dist = c()
gene_match_gene = c()

for(i in 1:length(pos.rda.top_10$SNP)) {
  SNP = pos.rda.top_10$SNP[i]
  chrom = pos.rda.top_10$chromosome[i]
  position = pos.rda.top_10$chrom_pos[i]
  print(chrom)
  #print(position)
  for(j in rownames(gene)[which(gene$gene_chrom == chrom)]) {
    g = gene[j, "genes"]
    #print(g)
    start = gene[j, "gene_start"]
    end = gene[j, "gene_end"]
    #print(start)
    #print(end)
    if(position >= start & position <= end) {
    gene_match_SNP = c(gene_match_SNP, SNP)
    gene_match_dist = c(gene_match_dist, 0)
    gene_match_gene = c(gene_match_gene, g)
    }
    else if(position >= (start - 30000) & position <= (end + 30000)) {
    gene_match_SNP = c(gene_match_SNP, SNP)
    gene_match_gene = c(gene_match_gene, g)
    gene_dist_1 = abs(start - position)
    gene_dist_2 = abs(end - position)
      if(gene_dist_1 < gene_dist_2) {
        gene_match_dist = c(gene_match_dist, gene_dist_1)
      }
      if(gene_dist_2 < gene_dist_1) {
        gene_match_dist = c(gene_match_dist, gene_dist_2)
      }
    }
  }
}

gene_match_rda = data.frame(gene_match_SNP, gene_match_dist, gene_match_gene)



gene_match_SNP = c()
gene_match_dist = c()
gene_match_gene = c()

for(i in 1:length(top_10_BH$SNP)) {
  SNP = top_10_BH$SNP[i]
  chrom = top_10_BH$chromosome[i]
  position = top_10_BH$chrom_pos[i]
  print(chrom)
  #print(position)
  for(j in rownames(gene)[which(gene$gene_chrom == chrom)]) {
    g = gene[j, "genes"]
    #print(g)
    start = gene[j, "gene_start"]
    end = gene[j, "gene_end"]
    #print(start)
    #print(end)
    if(position >= start & position <= end) {
      gene_match_SNP = c(gene_match_SNP, SNP)
      gene_match_dist = c(gene_match_dist, 0)
      gene_match_gene = c(gene_match_gene, g)
    }
    else if(position >= (start - 30000) & position <= (end + 30000)) {
      gene_match_SNP = c(gene_match_SNP, SNP)
      gene_match_gene = c(gene_match_gene, g)
      gene_dist_1 = abs(start - position)
      gene_dist_2 = abs(end - position)
      if(gene_dist_1 < gene_dist_2) {
        gene_match_dist = c(gene_match_dist, gene_dist_1)
      }
      if(gene_dist_2 < gene_dist_1) {
        gene_match_dist = c(gene_match_dist, gene_dist_2)
      }
    }
  }
}

gene_match_pcadapt = data.frame(gene_match_SNP, gene_match_dist, gene_match_gene)



gene_match_SNP = c()
gene_match_dist = c()
gene_match_gene = c()

for(i in 1:length(rownames(pos.outlier.overlap))) {
  #SNP = rownames(pos.outlier.overlap)[i]
  SNP = pos.outlier.overlap$SNP[i]
  chrom = pos.outlier.overlap$chromosome[i]
  position = pos.outlier.overlap$chrom_pos[i]
  print(chrom)
  #print(position)
  for(j in rownames(gene)[which(gene$gene_chrom == chrom)]) {
    g = gene[j, "genes"]
    #print(g)
    start = gene[j, "gene_start"]
    end = gene[j, "gene_end"]
    #print(start)
    #print(end)
    if(position >= start & position <= end) {
      gene_match_SNP = c(gene_match_SNP, SNP)
      gene_match_dist = c(gene_match_dist, 0)
      gene_match_gene = c(gene_match_gene, g)
    }
    else if(position >= (start - 30000) & position <= (end + 30000)) {
      gene_match_SNP = c(gene_match_SNP, SNP)
      gene_match_gene = c(gene_match_gene, g)
      gene_dist_1 = abs(start - position)
      gene_dist_2 = abs(end - position)
      if(gene_dist_1 < gene_dist_2) {
        gene_match_dist = c(gene_match_dist, gene_dist_1)
      }
      if(gene_dist_2 < gene_dist_1) {
        gene_match_dist = c(gene_match_dist, gene_dist_2)
      }
    }
  }
}

gene_match_overlap = data.frame(gene_match_SNP, gene_match_dist, gene_match_gene)



rownames(gene) = gene$genes

chr = c()
start = c()
end = c()

for(i in gene_match_rda$gene_match_gene) {
  chr = c(chr, gene[i, "gene_chrom"])
  start = c(start, gene[i, "gene_start"])
  end = c(end, gene[i, "gene_end"])
}

gene_match_rda = data.frame(gene_match_rda, chr, start, end)

chr = c()
start = c()
end = c()

for(i in gene_match_pcadapt$gene_match_gene) {
  chr = c(chr, gene[i, "gene_chrom"])
  start = c(start, gene[i, "gene_start"])
  end = c(end, gene[i, "gene_end"])
}

gene_match_pcadapt = data.frame(gene_match_pcadapt, chr, start, end)

chr = c()
start = c()
end = c()

for(i in gene_match_overlap$gene_match_gene) {
  chr = c(chr, gene[i, "gene_chrom"])
  start = c(start, gene[i, "gene_start"])
  end = c(end, gene[i, "gene_end"])
}

gene_match_overlap = data.frame(gene_match_overlap, chr, start, end)



##### write gene match output tables #####

print("writing gene match output tables...")

colnames(pos.outlier.overlap) = c('SNP', 'contig', 'P', 'pRDA axis', 'loading', 'chromosome', 'position')
colnames(top_10_BC) = c("SNP", "contig", "P", "chromosome", "position")
colnames(pos.rda.top_10) = c("SNP", "contig", "pRDA axis", "loading", "chromosome", "position")
colnames(gene_match_rda) = c("SNP", "gene_distance", "gene", "chromosome", "start", "end")
colnames(gene_match_pcadapt) = c("SNP", "gene_distance", "gene", "chromosome", "start", "end")
colnames(gene_match_overlap) = c("SNP", "gene_distance", "gene", "chromosome", "start", "end")

write.table(top_10_BC, "pcadapt.top10.pos.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(pos.rda.top_10, "rda.top10.pos.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(pos.outlier.overlap, "outlier.overlap.pos.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(gene_match_rda, "gene_match_rda.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(gene_match_pcadapt, "gene_match_pcadapt.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(gene_match_overlap, "gene_match_overlap.txt", quote=FALSE, sep="\t", row.names=FALSE)



##### write outlier SNP lists #####

print("writing outlier SNP lists...")

rda_contig = c()
rda_chrom = c()
rda_pos = c()
for(i in cand_sorted$SNP) {
  bam_row = bam.pos.info[i,]
  rda_contig = c(rda_contig, bam_row$contig)
  rda_chrom = c(rda_chrom, as.character(bam_row$chromosome))
  rda_pos = c(rda_pos, bam_row$chrom_pos)
}

rda_list_out = data.frame(cand_sorted$SNP, rda_contig, rda_chrom, rda_pos, 
                          cand_sorted$`RDA axis`, cand_sorted$Loading,
                          cand_sorted$ENV_PC1, cand_sorted$ENV_PC2,
                          cand_sorted$ENV_PC3)

colnames(rda_list_out) = c("SNP", "contig", "chromosome", "chromosome position", 
                           "pRDA axis", "loading", "ENV_PC1", "ENV_PC2", 
                           "ENV_PC3")

BH_SNP_list = rownames(pos[which(pos$outlier_BH == 1),])

pcadapt_contig = c()
pcadapt_chrom = c()
pcadapt_pos = c()
BC_list = c()
for(i in BH_SNP_list) {
  bam_row = bam.pos.info[i,]
  pcadapt_contig = c(pcadapt_contig, bam_row$contig)
  pcadapt_chrom = c(pcadapt_chrom, as.character(bam_row$chromosome))
  pcadapt_pos = c(pcadapt_pos, bam_row$chrom_pos)
  if(pos[i, "outlier_BC"] == 1) {
    BC_list = c(BC_list, 1)
  }
  else {
    BC_list = c(BC_list, 0)
  }
}

pcadapt_list_out = data.frame(BH_SNP_list, pcadapt_contig, pcadapt_chrom,
                              pcadapt_pos, BC_list)
colnames(pcadapt_list_out) = c("SNP", "contig", "chromosome", "chromosome position",
                               "Bonferroni correction")

write.table(rda_list_out[,1:4], "rda_list_out.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(pcadapt_list_out, "pcadapt_list_out.txt", quote=FALSE, sep="\t", row.names=FALSE)


##### plot structure run ln prob. of data #####

print("plotting structure est. ln prob. of data for each run...")

structure_probs = read.csv("structure_ln_prob_data.csv")
colnames(structure_probs) = c("Run", "K", "Est. Ln Prob. of Data")

ggplot(structure_probs, 
       aes(x=K, y=`Est. Ln Prob. of Data`)) + 
  xlab("K") +
  ylab("Est. Ln Prob. of Data") + 
  geom_point(alpha=0.5, size=4) + 
  geom_point(shape=1, size=4) +
  scale_x_continuous(breaks=seq(1, 13, by=1))



save.image(file="ddocent_aligned.RData")

#Genomics Project Roadmap and Pipeline

#aggregate data (identify homologs?) !!!!!!find out groups to root tree (I think I did this already)?

#clustering?

#perform multiple sequence alignment
#load required library
library(bio3d)

#read in unaligned fasta file of data (set your filepath as I have here)
dat_raw <- read.fasta("/Users/nicholashuron/Google Drive/Temple/BIOL_5403/project_data/HSF1_mRNA_all.fasta",
                      rm.dup = T)

#align the data with MUSCLE and save
##the function hates directories that have spaces in them for some reason
#you need to download muscle and know the filepath for it (insert in exefile)
#outfile is where you want the alignment to be saved on your hd
dat_aln <- seqaln(aln=dat_raw, exefile = "/Applications/Phylogenetics/MUSCLE-3.8.31/muscle3.8.31_i86darwin64",
       outfile = "/Users/nicholashuron/Desktop/tester/HSF1_mRNA_alignment.fasta",
       protein=FALSE,
       seqgroup = TRUE,
       verbose = TRUE
       )


#alternative with bioconductor package
#install biconductor
source("https://bioconductor.org/biocLite.R")
biocValid()
biocLite()
#install and load multiple sequence alignment package
biocLite("msa")
#after the first install, you only ever need to call msa library command to load the package
library("msa")

#get to work on reading in data and aligning
#filepath for your raw data from genbank
dat_raw <- Biostrings::readDNAStringSet(filepath = "/Users/nicholashuron/Google Drive/Temple/BIOL_5403/project_data/HSF1_mRNA_all.fasta",
                                        use.names = T)

#do the alignment from the data your wrote in with MUSCLE
dat_aln <- msaMuscle(inputSeqs = dat_raw,
          order = "input",
          verbose = TRUE
          )

#save the alignment as a FASTA file to the directory listed
writeXStringSet(x = as(dat_aln, "DNAStringSet"), 
                filepath = "/Users/nicholashuron/Google Drive/Temple/BIOL_5403/project_data/HSF1_mRNA_all_aln.fasta")

#optional convert to nexus format: some phylo programs want it in this format
library(seqinr)
library(ape)
#read in fasta saved alignment (aln); unnecessary if you have it saved as the object dat_aln above
dat_aln_fa <- seqinr::read.fasta(file = "/Users/nicholashuron/Google Drive/Temple/BIOL_5403/project_data/HSF1_mRNA_all_aln.fasta",
                                 seqtype = "DNA")

#save the same thing as a nexus file now!
write.nexus.data(x = dat_aln_fa,
                 file = "/Users/nicholashuron/Google Drive/Temple/BIOL_5403/project_data/HSF1_mRNA_all_aln.nex",
                 format = "dna")

#model selection for phylo tree

#phylo tree construction

#assess model?

#detect selective regimes with PAML

#gene tree - species tree reconciliation?
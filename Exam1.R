# Jonathan Bahamon 
# Febuary 20th, 2024
# Bioinformatics with Professor Johnson

# Packages that will be used for the Exam ####
# comment out any package install lines so that you don't accidentally 
# re-install them, which is very time consuming
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("Biostrings")

# install.packages("UniprotR")
# install.packages("protti")
# install.packages("r3dmol")
# install.packages("msa")
# install.packages("seqinr")
# install.packages("biomaRt")

# BiocManager::install("GenomicAlignments")

# Check to see if all the packages are correctly installed in R ####

library(msa)
library(UniprotR)
library(protti)
library(r3dmol)
library(GenomicAlignments)
library(Biostrings)
library(seqinr)
library(biomaRt)

# Importing and aligning the DNA sequences ####

# mySequenceFile1 <- system.file("Sequence_1",
#                                "Data/sequences.fasta",
#                                format = "fasta")
mySequences1 <- readDNAStringSet("sequences.fasta") # I didn't see any "Data" folder
mySequences1

# I see a possible mutation in Homo_sapien_4. It goes from C to A. 

print(mySequences1,
      show = "complete")
nchar(mySequences1)

detail(mySequences1)

msahum <- msa(mySequences1)
nchar(msahum)
print(msahum,
      show = "complete")

# Letter frequency function was not working for me 
NiceLetters <- alphabetFrequency(msahum)
print(NiceLetters)

# I wanted to see the middle ones specifically to compare them to each other as 
# well, which is why I separated them this way. 
mySequences1[6]
# There appears to be a point mutation in sequence 6 
mySequences1[7]
mySequences1[8]
mySequences1[9]
mySequences1[10]
# Sequence 10 has 2 point mutations from A to T. 
mySequences1[11]
mySequences1[12]
mySequences1[13]
mySequences1[14]
mySequences1[15]

# I'm not sure what the purpose of these next few lines of code is
# but they do run correctly
# Concentrate the sequences into a single string
entire_sequence <- paste(mySequences1, collapse = "NA")

# Print the entire DNA sequence without breaks
cat(entire_sequence)

#Use ClustalW for more information (msapackage)
# duplicate code. msa already assigned to variable 'msahum'
myClustalWAlignment <- msa(mySequences1,
                           "ClustalW")
myClustalWAlignment

print(myClustalWAlignment,
      show = "complete")
#Homo sapien 6 has a deletion on shown above consensus 
# if a comment line is very long, carry it over to the next line

# Number 2: ####
## The DNA sequences are typically very similar, expect for the point mutations 
# and deletion. Their width is all relatively similar as well. DNA sequences 4,6, 
# and 10 have point mutations. 
# but which of the three is the most different? Using the 'dist' function that 
# we used in the homework would have told you this.
# although I see that you figured it out in a circular way in the next question

# Number 3: ####
# I used blastn to find the gene that this may be. I exported it from the data 
# file itself that I have saved in the data folder(it should all be in Github).
# It deals with the hemoglobin beta chain in Homo sapiens (hbb gene). The accession 
# number that was an exact match LC1212775.1 (Query Cover 100%, E.value = 0.0, 
# Percent Identity = 100%). I input the different DNA sequences into blast and 
# compared it to the accession code LC1212775.1. Homo_sapien_6 had the lowest percent identity (99.84%). 

# Question 4: Translating the sequence to a protein ####
# Begin with separating the desired sequence from the entire fasta file (using Sequence 6 from mySequences[6]. 

class(mySequences1)

HOMO6 <- mySequences1$Homo_sapiens_6
HOMO6

# Use the package "Biostrings" to convert the isolated DNA sequence into an amino acid sequence.
# The last 2 bases were ignored with this code, but it should not affect the overall translation.

AAseq <- Biostrings::translate(HOMO6)
AAseq

# Writing the amino acid sequence to a FASTA file
# Begin by writing a header for the fasta file 

header <- "Translated.Amino.Acid.Sequence"

# Write the amino acid to a blank file 
# add a file type (fasta) to the file name
# this file didn't exist in GitHub, did you run it?
write.fasta(Biostrings::translate(HOMO6),
            "Homo_sapiens_6",
            "AAseq.fasta",
            open = "w",
            nbchar = 60,
            as.string = FALSE)
 
# What does the protein match to: ####
# A0A0J9YWK4 and it has catalytic activity in the HBB gene 
# Component: Chromosome 11 with the identifier "UP000005640"
# Although I didnt generate the 3D image, the image will be in the output folder 
# of the Exam 1 folder. It will also be in the Github repository. 

# Cross-referencing this with uniProt indicate beta Thalassemia (genetic blood 
# disorder in the HBB gene). This results in many complications such as: anemia, 
# bone problems, heart problems, and even growth and developmental issues. 
 
 









# Try2 ####
#DNA <- DNAString("AATCTACTCCCAGGAGCAGGGAGGGCAGGAGCCAGGGCTGGGCATGAAAGTCAGGGCAGAGCCATCTATTGCTTACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCACCTGACTCCTGTGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGGTGAGTCTATGGGACCCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAGTTCATGTCATAGGAAGGGG", nchar = NA)
#print(DNA)

#AA <- translate(DNA)

#writeXStringSet(XStringSet(AA),
#                file = "Translated_sequence.fasta",
#                format = "fasta",
#                names = header)
# The following code I want to use to perform a pairwise sequence alignment, but it wasn't working ####

# aligned_sequences <- pairwiseAlignment(mySequences2 , substitutionMatrix = "BLOSUM50", gapOpening = 2, gapExtension = 4) 
# print(aligned_sequences)

# I wanted to use the Blocks substitution matrix to assign a score to each possible substitution of a nucleotide for another in an alignment. Although the one I used for this exam is for protein sequence alignments optimized for alignments of closely related sequences. Based on DNAstringset and AAstringset, there seems to be some that were similar, but not all of them. This could be why it wasn't working. 
# The "gap0pening" is for a parameter to specift the penalty score assigned to opening a gap in the alignment. I changed the numeric value of it, but nothing occurred. When I took it out completely, the code also didn't work. There doesnt seem to be any gaps, but I will leave that part here for future alignments that may have it. 
# "gapExtension" is just to specify the penalty score for extending an existing gap. (DELETE)
                








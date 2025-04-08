# If you have a windows PC this code should be run in a https://posit.cloud/ R instance.
# There is a problem with the functions msaPrettyPrint() running on a windows machines.

# Install packages, if needed, remove leading # ----

#install.packages("BiocManager")
#BiocManager::install("msa")
#BiocManager::install("Biostrings")


## Load required libraries ----
library(Biostrings)
library(msa)

## Load fasta file ----
myFastaFile <- "https://raw.githubusercontent.com/JakeSaunders/BIO2910-Bioinformatics/refs/heads/main/cytochrome_c_oxidase_subunit_I.fasta"
mySequences <- readDNAStringSet(myFastaFile)

#DNA or RNA alignments ----
## Transforming sequences using biostrings package---- 

mySequences

# DNA string set is a special class of objects defined in the Biostrings package 
# that holds information about DNA, RNA or AA sequences. 
class(mySequences)

# this class holds the sequences and a lot of descriptive information about them 
# you can see a summary of this information by looking at the structure of the object
str(mySequences)

# the names() function will return the names of all the sequences in these objects
names(mySequences)

# you can use the $ to acces a specific sequence
mySequences$`C_001700.1:6216-7760 Felis catus mitochondrion, complete genome`

# yuou can also use [] and a range of numbers to call specific sequences
mySequences[6:7]
mySequences[c(2,5)]

#other information like the length of the seuqences can be accessed by the @
mySequences@ranges
# here is one way you can access the length of the sequences
mySequences@ranges@width

# You can also process the sequences in other ways
# you can calculate the complementry DNA
mySequences
reverseComplement(mySequences)

# you can translate from DNA to protein
translate(mySequences)

## Renameing our sequences to something more managable ----
names(mySequences)

short.names <- c("Felis catus (domestic cat)","Halichoerus grypus (grey seal)", 
  "Rattus rattus (common rat)", "Cavia porcellus (domestic guinea pig)", 
  "Tursiops truncatus (bottlenose dolphin)", "Balaenoptera musculus (blue whale)",
  "Loxodonta africana (African elephant)", "Oryctolagus cuniculus (European rabbit)",
  "Ornithorhynchus anatinus (platypus)" )

names(mySequences) <- paste("COXI", short.names)

names(mySequences)

## Align sequences ----

?msa

# There are different methods for doing MSAs, the most common ones (which are all
# supported by the msa package) are "ClustalW", "ClustalOmega", and "Muscle"
# verbose = True displays the processes being undertaken for the 

myFirstAlignment <- msa(mySequences, method ="Muscle", verbose = T)

# ClustalW: A widely used, progressive alignment algorithm, considered reliable but can be slower for large datasets compared to newer options like MUSCLE. 
# ClustalOmega: An updated version of ClustalW with a more efficient algorithm, often faster for large alignments, and uses a different method to calculate the guide tree. 
# MUSCLE: An algorithm known for its speed and accuracy, often considered the preferred choice for large-scale alignments due to its ability to quickly process large numbers of sequences while maintaining good alignment quality. 

### Save alignment as PDF ----

?msaPrettyPrint

msaPrettyPrint(myFirstAlignment,
               file="COXI-mammals.dvi", output="dvi",
               showNames = "left",
               askForOverwrite=FALSE, verbose=FALSE,
               shadingMode = "similar",consensusColors = "BlueRed",
               showLogo="none",)

msaPrettyPrint(myFirstAlignment,
               file="COXI-mammals.pdf", output="pdf",
               showNames="left", 
               askForOverwrite=FALSE, verbose=FALSE,
               shadingMode = "similar",consensusColors = "BlueRed",
               # Logo is one way to look at a summary of the similar sequences
               showLogo="top")

# Protein alignments ----

## Load fasta file ----
myProteins <- "https://raw.githubusercontent.com/JakeSaunders/BIO2910-Bioinformatics/refs/heads/main/myoglobin_hemoglobin.fasta"
myAASequences <- readAAStringSet(myProteins)

myAASequences
names(myAASequences)

## Align sequences ----
myProteinAlignment <- msa(myAASequences, method ="Muscle", verbose = T)



BrowseSeqs(myXStringSet = myAASequences)

### Save alignment as PDF ----

?msaPrettyPrint

msaPrettyPrint(myProteinAlignment,
               file="GlobinProteins.pdf", output="pdf",
               showNames = "left",
               askForOverwrite=FALSE, verbose=T,
               shadingMode = "similar",
               #consensusColors = "BlueRed",
               showLogo="top",
               logoColors = "hydropathy")

# logoColors=c("chemical", "rasmol", "hydropathy", "structure", "standard area", "accessible area"),

# for information on msaPrettyPrint opitions 
# https://bioconductor.org/packages/release/bioc/html/msa.html


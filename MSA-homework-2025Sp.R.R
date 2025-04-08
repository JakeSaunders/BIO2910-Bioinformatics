## Load required libraries ----
library(Biostrings)
library(msa)

## Load fasta file ----
myFastaFile <- "https://raw.githubusercontent.com/JakeSaunders/BIO2910-Bioinformatics/refs/heads/main/mouse_mGluR_Tas1R.fasta"
myAASequences <- readAAStringSet(myFastaFile)

# view sequences and full names ----
myAASequences
names(myAASequences)

# change names to shorter names

# save full names for reference
long.names <- names(myAASequences)

# assign variable of short names
short.names <- c("XP_006512612.1 mGluR1", "NP_001153825.1 mGluR2","XP_036020611.1 mGluR3",
                 "XP_011244777.1 mGluR4","XP_036008452.1 mGluR5","XP_017169707.1 mGluR6","XP_017176787.1 mGluR7",         
                 "NP_001348054.1 mGluR8","XP_006538537.1 Tas1R1","NP_114079.1 Tas1R2","NP_114078.1 Tas1R3")

# asssign short.names as names of biostring object
names(myAASequences) <- short.names

# drop taste receptors and mGluR1 and mGluR5 to make a better MSA
myAASequences <- myAASequences[c(2:4,6:8)]

## Align sequences ----
myProteinAlignment <- msa(myAASequences, method ="Muscle", verbose = T)

# view alignment 
myProteinAlignment

print(myProteinAlignment, show="complete",nameWidth=22)

# Printing MSA's -----

# note that the shading modes for this MSA has been sent to color code amino acids for hydropathy
# see key at bottom of MSA, and remember that hydrophobic non polar amino acids will be the ones
# that pass through the lipid membrane of cells. 

msaPrettyPrint(
  x = myProteinAlignment, file = "MetabotropicGlutamateReceptor.pdf",
  output="pdf", askForOverwrite=FALSE, #showNames="left",
  shadingMode = "functional",
  shadingModeArg = "hydropathy",
  showLogo = "top",
  logoColors= "hydropathy",
  showLogoScale="left",
  showLegend = TRUE,
)


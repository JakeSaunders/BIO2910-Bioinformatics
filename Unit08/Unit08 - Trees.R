# BIO 2910 · Unit 08 · Building and drawing phylogenetic trees in R
# Revised November 2026.
#
# What this script does
#   Example 1  Coronavirus RdRp genes: align, build a tree, ROOT it, BOOTSTRAP it, annotate by host
#   Example 2  Mammal COX1 genes: tree with PhyloPic silhouettes and habitat annotations
#   Example 3  Globin proteins: a gene tree with gene duplications labeled
#
# Run it one line at a time (Ctrl + Enter) so you can see what each step does.

# Install packages (run once, then leave commented out) ----
# install.packages("BiocManager")
# install.packages(c("tidyverse", "ape", "rphylopic"))
# BiocManager::install(c("Biostrings", "msa", "ggtree"))

library(tidyverse)
library(Biostrings)
library(msa)
library(ape)
library(ggtree)
library(rphylopic)

# All course data files live in the data/ folder of the course repo.
BASE <- "https://raw.githubusercontent.com/JakeSaunders/BIO2910-Bioinformatics/main/data/"


# =============================================================================
# EXAMPLE 1 · Coronaviruses ----
# =============================================================================

## 1.1 Load the sequences ----
# Partial RdRp (RNA-dependent RNA polymerase) genes from coronaviruses that infect
# different hosts. RdRp copies the viral genome, so every coronavirus has one.
cov.url <- paste0(BASE, "coronaviruses_RdRp_gene.fasta")
cov <- readDNAStringSet(cov.url)
cov
names(cov)

## 1.3 Multiple sequence alignment (Unit 07) ----
cov.aln <- msa(cov, method = "ClustalOmega")
cov.aln

## 1.4 Distance matrix ----
# Convert the alignment to ape's DNAbin format, then measure distances.
# model = "K80" corrects for sites that changed more than once.
# pairwise.deletion = TRUE ignores gaps pair by pair, because these are partial genes
# that cover different stretches of RdRp.
cov.dna <- msaConvert(cov.aln, type = "ape::DNAbin")
cov.dist <- dist.dna(cov.dna, model = "K80", pairwise.deletion = TRUE)
round(as.matrix(cov.dist), 2)

## 1.5 Build a neighbor-joining tree ----
cov.tree <- nj(cov.dist)
cov.tree          # an UNROOTED tree with 9 tips
is.rooted(cov.tree)

## 1.6 Root the tree with an outgroup ----
# The two Miniopterus bat viruses are ALPHAcoronaviruses. Every other virus here is a
# BETAcoronavirus, so the alphas fall outside the group we care about: an outgroup.
outgroup <- str_subset(cov.tree$tip.label, "Mi-CoV-1|Miniopterus")
outgroup
cov.tree <- root(cov.tree, outgroup = outgroup, resolve.root = TRUE)
is.rooted(cov.tree)

## 1.7 Bootstrap: how well supported is each clade? ----
# boot.phylo() resamples alignment columns 100 times, rebuilds the tree each time
# (using exactly the same steps as above), and counts how often each clade reappears.
make.cov.tree <- function(x) {
  root(nj(dist.dna(x, model = "K80", pairwise.deletion = TRUE)),
       outgroup = outgroup, resolve.root = TRUE)
}
set.seed(2910)    # makes the random resampling the same for everyone
cov.boot <- boot.phylo(cov.tree, cov.dna, FUN = make.cov.tree, B = 100, rooted = TRUE)
cov.boot          # one number per internal node: how many of 100 trees had that clade
cov.tree$node.label <- cov.boot

## 1.8 Draw the tree with ggtree ----
# ggtree is built on ggplot2. ggtree() is like ggplot(): it sets up the plot.
# Manual: https://yulab-smu.top/treedata-book/
ggtree(cov.tree)

# geom_tiplab() adds the tip labels
ggtree(cov.tree) +
  geom_tiplab()

# size makes the labels smaller; hexpand() adds room on the right so long labels fit
ggtree(cov.tree) +
  geom_tiplab(size = 3) +
  hexpand(0.4)

# color changes the label color; align = TRUE lines the labels up with dotted leaders
ggtree(cov.tree) +
  geom_tiplab(size = 3, color = "blue", align = TRUE) +
  hexpand(0.4)

# geom_nodelab() prints the bootstrap values at the nodes;
# geom_treescale() adds a scale bar (substitutions per site)
ggtree(cov.tree) +
  geom_tiplab(size = 3, align = TRUE) +
  geom_nodelab(size = 3, color = "darkred", hjust = 1.3, vjust = -0.5) +
  geom_treescale() +
  hexpand(0.4)

## 1.9 Same tree, different shapes ----
# The relationships are stored in the tree object. layout only changes the drawing.
ggtree(cov.tree, layout = "roundrect") + geom_tiplab(size = 3) + hexpand(0.4)
ggtree(cov.tree, layout = "slanted")   + geom_tiplab(size = 3) + hexpand(0.4)
ggtree(cov.tree, layout = "ellipse")   + geom_tiplab(size = 3) + hexpand(0.4)
ggtree(cov.tree, layout = "circular")  + geom_tiplab(size = 3)
ggtree(cov.tree, layout = "fan", open.angle = 90) + geom_tiplab(size = 3)
ggtree(cov.tree, layout = "unrooted")  + geom_tiplab(size = 3)   # ignores the root

## 1.10 Annotate the tips with a data frame ----
# The FIRST column must match the tip labels exactly. Here we pull the host out of
# each label, so the table can never get out of order with the tree.
cov.info <- tibble(
  label = cov.tree$tip.label,
  host  = case_when(
    str_detect(label, "human")    ~ "Human",
    str_detect(label, "bat")      ~ "Bat",
    str_detect(label, "pangolin") ~ "Pangolin",
    str_detect(label, "hare")     ~ "Hare",
    TRUE                          ~ "Other"
  )
)
cov.info

# %<+% joins the data frame to the tree; then host can be used inside aes()
ggtree(cov.tree) %<+% cov.info +
  geom_tiplab(size = 3, align = TRUE) +
  geom_tippoint(aes(shape = host, color = host), size = 4) +
  hexpand(0.4)

## 1.11 Finished figure ----
ggtree(cov.tree) %<+% cov.info +
  geom_tiplab(size = 3, align = TRUE, offset = 0.02) +
  geom_tippoint(aes(shape = host, color = host), size = 5) +
  geom_nodelab(size = 3, color = "grey30", hjust = 1.3, vjust = -0.5) +
  geom_treescale() +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  hexpand(0.5) +
  labs(title = "Coronavirus RdRp genes: neighbor-joining tree, rooted on alphacoronaviruses",
       subtitle = "Node labels = bootstrap support (100 replicates)",
       caption = "Scale bar: substitutions per site",
       shape = "Host", color = "Host")

## 1.12 Save it ----
# ggsave() saves the last plot you made
ggsave(filename = paste0("Tree - Coronavirus - ", Sys.Date(), ".pdf"),
       width = 11, height = 8.5, units = "in")


# =============================================================================
# EXAMPLE 2 · Mammal COX1 genes with silhouettes ----
# =============================================================================

## 2.1 Load the sequences (from the Databases unit) ----
cox1 <- readDNAStringSet(paste0(BASE, "cytochrome_c_oxidase_subunit_I.fasta"))
names(cox1)

## 2.2 Short tip labels ----
# Headers look like "NC_001700.1:6216-7760 Felis catus mitochondrion, complete genome".
# Keep just the scientific name, e.g. "Felis catus".
names(cox1) <- word(names(cox1), 2, 3)
names(cox1)

## 2.3 A table of information about each species ----
# Matched by name, not by row order, so the order of rows doesn't matter.
cox1.info <- tribble(
  ~label,                      ~common.name,          ~habitat,                           ~region,
  "Felis catus",               "Domestic cat",        "Human-associated",                 "Worldwide",
  "Halichoerus grypus",        "Grey seal",           "Marine",                           "North Atlantic",
  "Rattus rattus",             "Black rat",           "Human-associated",                 "Worldwide",
  "Cavia porcellus",           "Guinea pig",          "Human-associated",                 "South America",
  "Tursiops truncatus",        "Bottlenose dolphin",  "Marine",                           "Warm oceans",
  "Balaenoptera musculus",     "Blue whale",          "Marine",                           "All oceans",
  "Loxodonta africana",        "African elephant",    "Savanna",                          "Africa",
  "Oryctolagus cuniculus",     "European rabbit",     "Grassland",                        "Europe",
  "Ornithorhynchus anatinus",  "Platypus",            "Freshwater",                       "Australia"
)
# Check: every tip has a row in the table (should print character(0))
setdiff(names(cox1), cox1.info$label)

## 2.4 Align, measure, build, root ----
cox1.aln  <- msa(cox1, method = "ClustalOmega")
cox1.dna  <- msaConvert(cox1.aln, type = "ape::DNAbin")
cox1.tree <- nj(dist.dna(cox1.dna, model = "K80"))
# The platypus is a monotreme; every other species here is a placental mammal.
cox1.tree <- root(cox1.tree, outgroup = "Ornithorhynchus anatinus", resolve.root = TRUE)

set.seed(2910)
cox1.tree$node.label <- boot.phylo(cox1.tree, cox1.dna, B = 100, rooted = TRUE,
  FUN = function(x) root(nj(dist.dna(x, model = "K80")),
                         outgroup = "Ornithorhynchus anatinus", resolve.root = TRUE))

ggtree(cox1.tree) +
  geom_tiplab(size = 3.5, fontface = "italic", align = TRUE) +
  geom_nodelab(size = 3, color = "darkred", hjust = 1.3, vjust = -0.5) +
  hexpand(0.4)

## 2.5 Add PhyloPic silhouettes ----
# get_uuid() asks phylopic.org for a silhouette ID for each species.
# https://rphylopic.palaeoverse.org/reference/get_uuid.html
# If a species has no silhouette, PhyloPic returns its closest match.
cox1.info$uuid <- map_chr(cox1.info$label, ~ get_uuid(name = .x, n = 1))

ggtree(cox1.tree) %<+% cox1.info +
  geom_tiplab(aes(label = common.name), size = 3.5, align = TRUE, offset = 0.02) +
  geom_phylopic(aes(uuid = uuid), height = 0.6) +
  hexpand(0.5) +
  labs(title = "COX1 genes of nine mammals")

# Silhouettes are a geom like any other: change color and transparency
ggtree(cox1.tree) %<+% cox1.info +
  geom_tiplab(aes(label = common.name), size = 3.5, align = TRUE, offset = 0.02) +
  geom_phylopic(aes(uuid = uuid, color = habitat, x = x + 0.15), height = 0.6, alpha = 0.7) +
  scale_color_brewer(palette = "Dark2") +
  hexpand(0.5) +
  labs(title = "COX1 genes of nine mammals", color = "Habitat")

## 2.6 Finished figure: two annotations at once ----
# color = habitat (silhouettes), fill = region (tip points). Using color and fill
# for different variables keeps their two legends separate.
ggtree(cox1.tree) %<+% cox1.info +
  geom_tiplab(aes(label = common.name), size = 3.5, align = TRUE, offset = 0.04) +
  geom_tippoint(aes(color = habitat), size = 5, alpha = 0) +                 # hidden; builds the Habitat legend
  geom_tippoint(aes(fill = region), shape = 21, size = 4, color = "black") +
  geom_phylopic(aes(uuid = uuid, color = habitat, x = x + 0.10),
                height = 0.6, alpha = 0.7, show.legend = FALSE) +
  geom_nodelab(size = 3, color = "grey30", hjust = 1.3, vjust = -0.5) +
  geom_treescale() +
  scale_color_brewer(palette = "Dark2", name = "Habitat (animal silhouette)") +
  scale_fill_brewer(palette = "Set3", name = "Region (tip point)") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 16)),
         fill  = guide_legend(override.aes = list(shape = 21, size = 4))) +
  hexpand(0.6) +
  labs(title = "COX1 genes of nine mammals",
       subtitle = "Neighbor-joining, rooted on platypus; node labels = bootstrap support",
       caption = "Scale bar: substitutions per site")

ggsave(filename = paste0("Tree - COX1 - ", Sys.Date(), ".pdf"),
       width = 11, height = 8.5, units = "in")

# Question to think about: are the two rodents (rat and guinea pig) in one clade?
# Is that clade well supported? What does that tell you about using one gene?


# =============================================================================
# EXAMPLE 3 · Globins: a gene tree with duplications ----
# =============================================================================

## 3.1 Load protein sequences ----
glob <- readAAStringSet(paste0(BASE, "myoglobin_hemoglobin.fasta"))
names(glob)

## 3.2 Short labels: gene + species ----
# "NP_000508.1 hemoglobin subunit alpha [Homo sapiens]"  ->  "HbA Homo sapiens"
gene <- case_when(
  str_detect(names(glob), "myoglobin") ~ "Mb",
  str_detect(names(glob), "alpha")     ~ "HbA",
  str_detect(names(glob), "beta")      ~ "HbB"
)
species <- str_extract(names(glob), "(?<=\\[).+(?=\\])")
names(glob) <- paste(gene, species)
names(glob)

## 3.3 Align, measure, build ----
glob.aln  <- msa(glob, method = "ClustalOmega")
glob.aa   <- msaConvert(glob.aln, type = "ape::AAbin")
glob.tree <- nj(dist.aa(glob.aa, pairwise.deletion = TRUE, scaled = TRUE))

## 3.4 Root on the myoglobins ----
# Myoglobin split off before hemoglobin alpha and beta split from each other,
# so the myoglobins are the outgroup for the hemoglobins.
mb <- str_subset(glob.tree$tip.label, "^Mb")
glob.tree <- root(glob.tree, outgroup = mb, resolve.root = TRUE)

set.seed(2910)
glob.tree$node.label <- boot.phylo(glob.tree, glob.aa, B = 100, rooted = TRUE,
  FUN = function(x) root(nj(dist.aa(x, pairwise.deletion = TRUE, scaled = TRUE)),
                         outgroup = mb, resolve.root = TRUE))

ggtree(glob.tree) +
  geom_tiplab(size = 3) +
  geom_nodelab(size = 3, color = "darkred", hjust = 1.3, vjust = -0.5) +
  hexpand(0.3)

## 3.5 Find the node for each gene family ----
# getMRCA() (from ape) returns the node number of the most recent common ancestor of a set of tips.
# Finding nodes this way still works if the tree changes; typing node numbers does not.
hba.node <- getMRCA(glob.tree, str_subset(glob.tree$tip.label, "^HbA"))
hbb.node <- getMRCA(glob.tree, str_subset(glob.tree$tip.label, "^HbB"))
mb.node  <- getMRCA(glob.tree, mb)
c(HbA = hba.node, HbB = hbb.node, Mb = mb.node)

# To see every node number on the tree:
ggtree(glob.tree) + geom_tiplab(size = 3) + geom_text(aes(label = node), hjust = -0.3, size = 3)

## 3.6 Highlight and label the clades ----
ggtree(glob.tree) +
  geom_hilight(node = mb.node,  fill = "darkseagreen", alpha = 0.4) +
  geom_hilight(node = hba.node, fill = "gold",         alpha = 0.4) +
  geom_hilight(node = hbb.node, fill = "dodgerblue",   alpha = 0.3) +
  geom_tiplab(size = 3) +
  geom_nodelab(size = 2.8, color = "grey30", hjust = 1.3, vjust = -0.5) +
  geom_cladelab(node = mb.node,  label = "Myoglobin",        offset = 0.25, textcolor = "darkgreen") +
  geom_cladelab(node = hba.node, label = "Hemoglobin\nalpha", offset = 0.25, textcolor = "darkgoldenrod") +
  geom_cladelab(node = hbb.node, label = "Hemoglobin\nbeta",  offset = 0.25, textcolor = "blue") +
  geom_treescale() +
  hexpand(0.5) +
  labs(title = "The evolution of globin proteins in vertebrates",
       subtitle = "Neighbor-joining on protein p-distance, rooted on myoglobin; node labels = bootstrap support",
       caption = "Scale bar: proportion of amino acid sites that differ")

ggsave(filename = paste0("Tree - Globin - ", Sys.Date(), ".pdf"),
       width = 8.5, height = 8.5, units = "in")

# Question to think about: within hemoglobin beta, does the order in which the animals
# branch off match the species tree (fish, then frog, then chicken, then mammals)?

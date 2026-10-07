# BIO 2910 - Unit 05 - NMDS of dead wood fungal communities ----
# Data: presence (1) / absence (0) of saproxylic fungi at 14 New Jersey parks, 2017-2018
# vegan does the statistics; ggplot2 does all of the plotting.

# Load packages ----
# If a package is missing, remove the # and run the install line once.
# install.packages(c("tidyverse", "vegan", "ggrepel", "usmap", "ggspatial"))

library(tidyverse)
library(vegan)     # vegdist(), metaMDS(), envfit(), scores()
library(ggrepel)   # geom_text_repel(): labels that don't overlap

# Load and explore data ----
# read.csv() can read straight from a URL; no download step needed
deadwood_class <- read.csv(
  "https://raw.githubusercontent.com/JakeSaunders/BIO2910-Bioinformatics/refs/heads/main/data/NJ.fungus.csv"
)

# View() opens a spreadsheet tab; you can also click the object in the Environment pane
View(deadwood_class)

dim(deadwood_class)
head(deadwood_class)
tail(deadwood_class)
str(deadwood_class)

# names() lists the column names, of which there are many
names(deadwood_class)

# QUESTIONS 1-4 ----

# Clean up the site information ----
# The locality names are too long to use as plot labels, so add a short name column.
# Adding it as a column (instead of a separate vector) keeps each name attached to its row.
# habitat has a stray space ("Hardwood "); trimws() removes it so labels look right.
deadwood_class <- deadwood_class %>%
  mutate(habitat = trimws(habitat),
         site = c("Wawayanda", "Stokes", "Weiss", "Meadowood", "Teetertown", "Schiff",
                  "Thompson", "Ocean County", "Forest REC", "Rancocas", "Franklin Parker",
                  "Wells Mills", "Belleplain", "Cattus Island"))

deadwood_class %>% select(site, locality, habitat)

# Keep only the species columns ----
# NMDS compares communities using a distance metric. Communities that are most
# similar in species composition should end up closest together on the plot.
# vegan needs a data frame of ONLY numbers: one row per site, one column per species.

# Option 1: square brackets, keeping columns 12 to 260
deadwood1 <- deadwood_class[ , 12:260]

# Option 2: select() and drop the descriptive columns with a minus sign
deadwood2 <- deadwood_class %>%
  select(-locality, -habitat, -decimalLatitude, -decimalLongitude, -geodeticDatum,
         -eventDate, -basisOfRecord, -countryCode, -stateProvince, -kingdom,
         -recordedBy, -site)

# Option 3: select() a range of consecutive columns with first:last
deadwood3 <- deadwood_class %>%
  select(Abortiporus.biennis:Xylobolus.subpileatus)

# identical() returns TRUE if two objects are exactly the same
identical(deadwood1, deadwood2)
identical(deadwood1, deadwood3)

# Distance matrix ----
# Jaccard distance = 1 - (species shared / total species in both sites)
# 0 = identical communities, 1 = no species in common
wood.dis <- vegdist(deadwood1, method = "jaccard")
wood.dis

# Run the NMDS ----
# metaMDS() rebuilds the distance matrix itself, then tries many random starting
# layouts and keeps the one with the lowest stress (the best fit).
# Stress: < 0.05 excellent, < 0.1 great, < 0.2 OK, > 0.3 poor.

# NMDS starts from random layouts. set.seed() makes the "random" numbers repeat,
# so everyone in class gets the same answer.
set.seed(2910)
wood.mds <- metaMDS(deadwood1, distance = "jaccard", k = 2, trymax = 100)

# Printing the object summarizes the analysis
wood.mds

# vegan made its own object class to hold the results and metadata
class(wood.mds)

# QUESTIONS 5 & 6 ----

# Stress plot (Shepard plot) ----
# Does the ordination keep the rank order of the original distances?
# If the fit is good, the blue points follow the red step line closely.

# vegan's built-in version, in base R graphics:
stressplot(wood.mds,
  main = "Stress Plot of Dead Wood Dataset")


# QUESTION 7 ----

# Get the NMDS coordinates into data frames ----
# vegan's own plots (ordiplot, ordihull) are base R graphics. 

ordiplot(wood.mds)
ordihull(wood.mds)

# To use ggplot2 we pull out the coordinates ("scores") and add them to a data frame.

# one row per site: NMDS1, NMDS2, plus the site information
site_scores <- as.data.frame(scores(wood.mds, display = "sites")) %>%
  bind_cols(deadwood_class %>% select(site, locality, habitat,
                                      decimalLatitude, decimalLongitude))
site_scores

# one row per species. 26 species were never found at these 14 sites (all 0s),
# so they have no position; filter() drops those NA rows.
species_scores <- as.data.frame(scores(wood.mds, display = "species")) %>%
  rownames_to_column("species") %>%
  filter(!is.na(NMDS1))
head(species_scores)

# Plot the NMDS ----

## Sites and species together (vegan's version is ordiplot(wood.mds))
ggplot() +
  geom_point(data = species_scores, aes(x = NMDS1, y = NMDS2),
             shape = 3, colour = "red") +
  geom_point(data = site_scores, aes(x = NMDS1, y = NMDS2)) +
  coord_equal() +
  theme_classic()

## Sites only, with different point sizes
site_scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 6) +
  coord_equal() +
  theme_classic()

site_scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 3) +
  coord_equal() +
  theme_classic()

site_scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 0.5) +
  coord_equal() +
  theme_classic()

# QUESTION 8 ----

## Colour the points by habitat
site_scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2, colour = habitat)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("forestgreen", "sienna")) +
  coord_equal() +
  theme_classic()

# Add hulls around each habitat ----
# A hull is the smallest polygon that encloses a group of points.
# chull() returns which rows are the corners, so slice() keeps just those rows,
# separately for each habitat because of group_by().
hulls <- site_scores %>%
  group_by(habitat) %>%
  slice(chull(NMDS1, NMDS2))
hulls

site_scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2, colour = habitat)) +
  geom_polygon(data = hulls, aes(fill = habitat), alpha = 0.2, linewidth = 1) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("forestgreen", "sienna")) +
  scale_fill_manual(values = c("forestgreen", "sienna")) +
  coord_equal() +
  theme_classic()

# Add site labels ----
# geom_text_repel() nudges labels so they don't overlap the points or each other
nmds_plot <- site_scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2, colour = habitat)) +
  geom_polygon(data = hulls, aes(fill = habitat), alpha = 0.2, linewidth = 1) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = site), colour = "black", size = 3) +
  scale_colour_manual(values = c("forestgreen", "sienna")) +
  scale_fill_manual(values = c("forestgreen", "sienna")) +
  labs(title = "Fungal communities differ by habitat type",
       subtitle = paste("NMDS, Jaccard distance, stress =", round(wood.mds$stress, 3)),
       colour = "Habitat", fill = "Habitat") +
  coord_equal() +
  theme_classic()

nmds_plot

# QUESTION 9 ----

# Fit environmental variables ----
# Do latitude or longitude line up with the community differences?
# envfit() tests each variable against the NMDS layout with permutations.
variables <- deadwood_class %>% select(decimalLatitude, decimalLongitude)
variables

set.seed(2910)
ef <- envfit(wood.mds, variables, permutations = 999)
ef

# QUESTION 10 ----

# Add the significant variables as arrows ----
# scores() gives each arrow's direction (already scaled by its correlation).
# The arrows are short, so stretch the longest one to about 80% of the plot's size.
arrows <- as.data.frame(scores(ef, display = "vectors")) %>%
  rownames_to_column("variable") %>%
  mutate(pval = ef$vectors$pvals,
         variable = str_remove(variable, "decimal")) %>%   # "Latitude" fits better
  filter(pval <= 0.1)

stretch <- 0.8 * max(abs(c(site_scores$NMDS1, site_scores$NMDS2))) /
  max(sqrt(arrows$NMDS1^2 + arrows$NMDS2^2))
arrows <- arrows %>%
  mutate(NMDS1 = NMDS1 * stretch,
         NMDS2 = NMDS2 * stretch)
arrows

# add layers to the plot we already saved
nmds_plot +
  geom_segment(data = arrows,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               colour = "blue", linewidth = 1,
               arrow = arrow(length = unit(0.3, "cm")),
               inherit.aes = FALSE) +
  # hjust/vjust tuck the label beside the arrow tip so it isn't cut off
  geom_text(data = arrows,
            aes(x = NMDS1, y = NMDS2, label = variable,
                hjust = ifelse(NMDS1 < 0, -0.2, 1.2),
                vjust = ifelse(NMDS2 < 0, 1.2, -0.5)),
            colour = "blue", inherit.aes = FALSE)

# QUESTION 11 ----
# Save the final plot for QUESTION 12 (ggsave saves the last plot shown)
ggsave("NMDS_plot.png", width = 7, height = 6)

# Map the sites ----
library(usmap)
library(ggspatial)

# plot_usmap() draws maps of the US as ggplot objects
plot_usmap(regions = "states")

plot_usmap(include = .northeast_region, labels = TRUE)

plot_usmap("counties", include = "NJ", labels = TRUE)

plot_usmap("counties", include = "NJ", labels = TRUE) +
  # coord_sf() lets the map use longitude and latitude
  coord_sf(crs = usmap_crs()) +
  # geom_spatial_ functions plot using longitude and latitude as x and y
  geom_spatial_point(data = deadwood_class,
                     aes(x = decimalLongitude, y = decimalLatitude))

plot_usmap("counties", include = "NJ", labels = FALSE) +
  coord_sf(crs = usmap_crs()) +
  geom_spatial_point(data = deadwood_class,
                     aes(x = decimalLongitude, y = decimalLatitude, colour = habitat),
                     size = 6, alpha = 0.5) +
  scale_colour_manual(values = c("forestgreen", "sienna")) +
  theme(legend.position = "top")

plot_usmap("counties", include = "NJ", labels = FALSE) +
  coord_sf(crs = usmap_crs()) +
  geom_spatial_point(data = deadwood_class,
                     aes(x = decimalLongitude, y = decimalLatitude, colour = habitat),
                     size = 6, alpha = 0.4) +
  geom_spatial_text_repel(data = deadwood_class,
                          aes(x = decimalLongitude, y = decimalLatitude, label = site),
                          size = 3) +
  scale_colour_manual(values = c("forestgreen", "sienna")) +
  labs(title = "Location of Deadwood Survey Sites",
       subtitle = "New Jersey, United States",
       colour = "Habitat Type") +
  theme(legend.position = "top")

# QUESTION 12: save the map too
ggsave("NMDS_map.png", width = 6, height = 8)


# Panel plot: NMDS and map side by side ----
# Run this AFTER the rest of the script; it reuses nmds_plot, arrows
# and deadwood_class. ggarrange() comes from ggpubr (used in Unit 04).
library(ggpubr)

# OPTIONAL: Combine plots
library(ggpubr)

# each one must be saved as an object with <-
# A: NMDS with hulls, labels and the latitude arrow
nmds_final <- nmds_plot +
  geom_segment(data = arrows,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               colour = "blue", linewidth = 1,
               arrow = arrow(length = unit(0.3, "cm")),
               inherit.aes = FALSE) +
  geom_text(data = arrows,
            aes(x = NMDS1, y = NMDS2, label = variable,
                hjust = ifelse(NMDS1 < 0, -0.2, 1.2),
                vjust = ifelse(NMDS2 < 0, 1.2, -0.5)),
            colour = "blue", inherit.aes = FALSE)

# B: map of the sites
site_map <- plot_usmap("counties", include = "NJ", labels = FALSE) +
  coord_sf(crs = usmap_crs()) +
  geom_spatial_point(data = deadwood_class,
                     aes(x = decimalLongitude, y = decimalLatitude, colour = habitat),
                     size = 6, alpha = 0.4) +
  geom_spatial_text_repel(data = deadwood_class,
                          aes(x = decimalLongitude, y = decimalLatitude, label = site),
                          size = 3) +
  scale_colour_manual(values = c("forestgreen", "sienna")) +
  labs(title = "Location of Deadwood Survey Sites",
       subtitle = "New Jersey, United States",
       colour = "Habitat Type") +
  theme(legend.position = "none") # panel A already has the legend

# Arrange: NMDS on the left (wider), map on the right
panel <- ggarrange( site_map, nmds_final,
                   ncol = 2, 
                   widths = c(1, 3), 
                   labels = c("A", "B"))
panel

ggsave("NMDS_panel.png", panel, width = 10, height = 6.5, bg = "white")



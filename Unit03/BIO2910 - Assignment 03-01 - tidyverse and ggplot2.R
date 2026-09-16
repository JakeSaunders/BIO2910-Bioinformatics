# Tidyverse and download data ----
# install the package tidyverse
install.packages("tidyverse")
install.packages("ggthemes")

# load tidyverse packages 
library(tidyverse)
library(ggthemes)

# Run the function below to download a large data set on finches 
download.file(
  url = "https://raw.githubusercontent.com/JakeSaunders/BIO2910-Bioinformatics/refs/heads/main/data/Bio2910-Finches.csv",
  destfile = "Bio2910-Finches.csv"
)

# ...And read the file into R with the following code
# In this table Weight is recorded in grams, all other measurements are cm.

finches <- read.csv("Bio2910-Finches.csv")

# What is the Tidyverse ----

# This is the pipe, you can think of it as meaning "...and then do "

# %>% 

# Clt + Shift + M is a hotkey that will write a %>% 

# Tidyverse was primarly developed by a biostatstiction named Hadley Wickham
# the goal was to make coding in R more linear/intuitive and most of the 
# books he has written on R-programing are avaiable for free (see links on his website)
  # https://hadley.nz/
  # https://github.com/hadley
  # https://en.wikipedia.org/wiki/Hadley_Wickham

# The tidyverse is built around a few principles about how data should 
# be entered in spreadsheets or data frames
  # https://vita.had.co.nz/papers/tidy-data.pdf

# Principle 1: Each variable is stored in a single column
# Principle 2: Each observation of a variable is stored in a different row
# Principle 3: Each cell stores a single value

# Using the Tidyverse to summerize data ----

finches %>%
  # group_by is a function that selects the variables (columns) you are interested in
  # think of these as possible independent variables 
  group_by( Species, Sex, Year) %>% 
  # summarise is a function that allows you to calculate summary numbers for dependent variables
  # and name the column these calculated values will fit in 
  # all unneed columns are dropped
  summarise( 
    # the format for nex variables is as follows:
    # new.column.name = function to calculate the summary value
    Weight.avg = mean(Weight), 
    Weight.sd = sd(Weight), 
    # n() is a function that just returns the replicate count
    N = n()
  )

# Notice that all the replicates in this data set are from the same species 
# Therefore We can simplify the data dropping the species from the list of 
# grouping variables like so:

finches %>%
  group_by( Sex, Year) %>% 
  summarise( 
    Weight.avg = mean(Weight), 
    Weight.sd = sd(Weight), 
    N = n()
  )

# when you get the data summarized how you want just add an assignment operator
# at the start of the pipeline and name a new object

df.weight <- finches %>%
  group_by( Sex, Year) %>% 
  summarise( 
    Weight.avg = mean(Weight), 
    Weight.sd = sd(Weight), 
    N = n() 
  )

# ggplot2 ----

# Watch this video explaining how to use ggplot2 
# https://www.youtube.com/watch?v=FdVy57oGJuc&t=3s

# cheat sheets
# https://rstudio.github.io/cheatsheets/
# ggplot cheat sheet:
# https://rstudio.github.io/cheatsheets/html/data-visualization.html


## bar plot ggplot ----
ggplot(data = df.weight, aes(x = Year, y = Weight.avg, fill = Sex)) +
  geom_col(position = "dodge") +
  xlab("Year") + 
  ylab("Weight (g)") +
  labs(title = "Finiches of all Sexes Gained Weight between 1977 and 1978") 


ggplot(data = df.weight, aes(x = Year, y = Weight.avg, fill = Sex)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes( ymin = Weight.avg - Weight.sd, ymax = Weight.avg + Weight.sd),
                position = "dodge") +
  xlab("Year") + 
  ylab("Weight (g)") +
  labs(title = "Finiches of all Sexes Gained Weight between 1977 and 1978") 


## scatter plot ggplot ----

ggplot(data = finches, aes(x = Weight, y = Wing)) +
  geom_point() 


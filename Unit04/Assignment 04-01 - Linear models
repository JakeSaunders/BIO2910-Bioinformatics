# BIO 2910 - Unit 04 - Linear models ----
# Data: Darwin's medium ground finch (Geospiza fortis), Daphne Major, 1977-1978
# Weight is in grams (g). All other measurements are in millimeters (mm).

# Load packages ----
# If a package is missing, remove the # and run the install line once.
# install.packages(c("tidyverse", "ggpubr", "GGally"))

library(tidyverse)
library(ggpubr)   # theme_classic2(), stat_cor(), stat_regline_equation()

# Load data ----
# read.csv() can read straight from a URL; no download step needed
finches <- read.csv(
  "https://raw.githubusercontent.com/JakeSaunders/BIO2910-Bioinformatics/main/data/Bio2910-Finches.csv"
)

# Year is a label, not a quantity, so treat it (and Sex) as categorical factors.
# Doing this once here means every plot and model below agrees.
finches <- finches %>%
  mutate(Year = factor(Year),
         Sex  = factor(Sex))

str(finches)

# Scatter plots and trend lines in ggplot ----

## Basic scatter plot
finches %>%
  ggplot(aes(x = Wing, y = Weight)) +
  geom_point(alpha = 0.5) +
  labs(title = "Heavier finches have longer wings",
       x = "Wing (mm)", y = "Weight (g)") +
  theme_classic2()

## Add a trend line: the default is a LOESS curve (local smoothing), not a straight line
finches %>%
  ggplot(aes(x = Wing, y = Weight)) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  labs(title = "Heavier finches have longer wings",
       x = "Wing (mm)", y = "Weight (g)") +
  theme_classic2()

## Add a straight trend line: method = "lm" fits a linear model
finches %>%
  ggplot(aes(x = Wing, y = Weight)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x) +
  labs(title = "Heavier finches have longer wings",
       x = "Wing (mm)", y = "Weight (g)") +
  theme_classic2()

## Add an Excel-style equation, R-squared and p-value
finches %>%
  ggplot(aes(x = Wing, y = Weight)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x) +
  stat_regline_equation(label.y = 21.5) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")),
           label.y = 20.8) +
  labs(title = "Heavier finches have longer wings",
       x = "Wing (mm)", y = "Weight (g)",
       # always tell the reader what annotations mean
       caption = "Shaded area is the 95% confidence interval") +
  theme_classic2()

## ggplot makes it easy to add another variable: color by Year
finches %>%
  ggplot(aes(x = Wing, y = Weight, color = Year)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x) +
  stat_regline_equation(aes(label = paste(after_stat(eq.label), after_stat(rr.label),
                                          sep = "~~~")),
                        label.y = c(21.5, 20.8)) +
  labs(title = "Finches were heavier in 1978 at every wing length",
       x = "Wing (mm)", y = "Weight (g)",
       caption = "Shaded area is the 95% confidence interval") +
  theme_classic2()

# Explore every pair of variables at once ----

morph <- c("Weight", "Wing", "Tarsus", "Beaklength", "Beakdepth", "Beakwidth")
pairs(finches[,3:10])

GGally::ggpairs(finches)

# Linear models with lm() ----
# Formula syntax is  response ~ predictor(s),  read as "Weight is predicted by Wing"

## Model 1: does wing length predict weight? ----
lm1 <- lm(Weight ~ Wing, data = finches)
lm1           # just the intercept and slope
summary(lm1)  # slope, p-values, R-squared
anova(lm1)    # the same model as an ANOVA table
class(lm1)    # a fitted model is its own object class: "lm"

### Check the model's assumptions: four diagnostic plots on one page
# 1. Residuals vs Fitted (linearity): points should scatter randomly around 0; a curve means the relationship isn't linear
# 2. Normal Q-Q (normality of residuals): points should follow the diagonal; S-shapes or strong tail deviation are a problem
# 3. Scale-Location (equal variance): red line should be roughly flat; a funnel or upward trend means variance changes with y
# 4. Residuals vs Leverage (influential points): watch for points beyond the dashed Cook's distance lines; R labels extreme points by row number, e.g. finches[72, ]
par(mfrow = c(2, 2))
plot(lm1)
par(mfrow = c(1, 1))



## Model 2: do wing AND tarsus length predict weight? ----
lm2 <- lm(Weight ~ Wing + Tarsus, data = finches)
summary(lm2)
anova(lm2)   # note: terms are tested in order, so Weight ~ Tarsus + Wing gives a different table

par(mfrow = c(2, 2))
plot(lm2)
par(mfrow = c(1, 1))

## Model 3: add Sex, a categorical predictor ----
# R turns a factor into 0/1 "dummy" variables. "female" is the reference level
# (first alphabetically), so "Sexmale" = difference between males and females.
# Note: about a third of these birds have Sex = "unknown", which R treats as its own group.
lm3 <- lm(Weight ~ Wing + Tarsus + Sex, data = finches)
summary(lm3)
anova(lm3)

## Model 4: add Year ----
lm4 <- lm(Weight ~ Wing + Tarsus + Year, data = finches)
summary(lm4)
anova(lm4)

par(mfrow = c(2, 2))
plot(lm4)
par(mfrow = c(1, 1))

# Which model is best? ----
# Let's get in the weeds here, this isn't something you are expected to do on your own for the exam

# anova(small, big) compares two models with an F-test.
# It asks: does adding predictors reduce the unexplained error (RSS) more than expected by chance?
# The models must be NESTED: the big model contains every predictor in the small one, plus more.
#   lm1 (Wing) is nested in lm2 (Wing + Tarsus), so they can be compared.
#   lm3 (+ Sex) and lm4 (+ Year) are NOT nested in each other, so anova() can't compare them.
# How to read the output:
#   RSS    = residual sum of squares, the error left over (smaller = better fit)
#   Df     = number of predictors added
#   Pr(>F) = p-value; if < 0.05, the bigger model fits significantly better, so keep the extra predictor(s)
anova(lm1, lm2)   # does adding Tarsus help? (p = 0.0002, yes)
anova(lm2, lm4)   # does adding Year help? (p = 0.018, yes)

# AIC (Akaike Information Criterion) scores each model's fit, with a penalty for every predictor added.
# This stops you from rewarding a model just for being more complicated.
# Lower AIC = better. Differences of less than about 2 are basically a tie.
# Unlike anova(), AIC can compare non-nested models, as long as they use the same data and response.
#   df = number of parameters estimated (the intercept + the predictors + the residual error)
AIC(lm1, lm2, lm3, lm4)
# lm4 (Wing + Tarsus + Year) has the lowest AIC, so it is the best of these four.
# lm3 scores worse than lm2: Sex added 2 parameters but didn't improve fit enough to justify them.

# A t-test is a linear model ----
# Same p-value both ways. Look at the Year1978 row of the lm summary.
t.test(Weight ~ Year, data = finches, var.equal = TRUE)
summary(lm(Weight ~ Year, data = finches))

# Showing several variables in one plot ----

## BAD: too many trend lines ----
# shape = Sex is in ggplot(aes()), so geom_smooth() inherits it and fits one line
# per Year x Sex combination: 2 x 3 = 6 lines instead of 2.
finches %>%
  ggplot(aes(x = Wing, y = Weight, color = Year, shape = Sex)) +
  geom_point(aes(size = Tarsus), alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  stat_regline_equation(label.x = 62, label.y = seq(23, 18, length.out = 6)) +
  stat_cor(label.x = 66.5, label.y = seq(23, 18, length.out = 6)) +
  labs(title = "Finch weight by wing length, tarsus, sex, and year",
       x = "Wing (mm)", y = "Weight (g)", size = "Tarsus (mm)") +
  theme_classic2()

## FIX: one trend line per year ----
# Aesthetics in ggplot(aes()) apply to every layer; aesthetics in geom_point(aes()) apply only to the points.
# Moving shape = Sex into geom_point() means geom_smooth() only splits by Year.
finches %>%
  ggplot(aes(x = Wing, y = Weight, color = Year)) +
  geom_point(aes(size = Tarsus, shape = Sex), alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  stat_regline_equation(label.x = 62, label.y = c(22, 21.3)) +
  stat_cor(label.x = 65.5, label.y = c(22, 21.3)) +
  labs(title = "Finch weight by wing length, tarsus, sex, and year",
       x = "Wing (mm)", y = "Weight (g)", size = "Tarsus (mm)") +
  theme_classic2()


## Another Fix: one trend line total ----
# override  aes arguments in geom_smooth
finches %>%
  ggplot(aes(x = Wing, y = Weight, color = Year)) +
  geom_point(aes(size = Tarsus, shape = Sex), alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black") +
  stat_regline_equation(label.x = 62, label.y = 22, color = "black") +
  stat_cor(label.x = 65.5, label.y = 22, color = "black") +
  labs(title = "Finch weight by wing length, tarsus, sex, and year",
       x = "Wing (mm)", y = "Weight (g)", size = "Tarsus (mm)") +
  theme_classic2()


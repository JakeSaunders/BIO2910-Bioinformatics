## Exploring the iris data set ----

# iris is a dataset R has for ppl to practice on
iris

#check iris help file
?iris

# assign iris to a var call dat
dat <- iris

# explore dat

# class function check what class an object is
class(dat)
# putting a ? before a function opens help
?class

# how big is dat, dim gives number of rows and columns
dim(dat)

# what kind of data is stored in dat, use structure func str()
str(dat)

# plot ever column against every other column
pairs(dat)


## Histogram ----

# look up hist help
?hist

#select Sepal.Length column
# 2 ways 
# first way, sub set using []
# var[ #row , #col  ]
dat[ , 1 ]

dat[ , 1:2 ]

dat[ 1:10 , 1:2 ]

# second way, sub set using $
# var$column.name

dat$Sepal.Length


# Make histogram of sepal length

hist(
  # x is an argument for the values to be plotted
  x = dat$Sepal.Length, 
  # number of 
  breaks = 10,
  # change bar color of bars and boarders
  col = "gold",border = "navyblue",
  # main titel
  main = "Sepal Length of Three Species of Iris",
  # x-axis label
  xlab = "Sepal Length (cm)"
)


## barplot ----

#colMeans calculates means for each column of a data frame
colMeans(dat[1:4])

?barplot

barplot(
  # make bar heights equal to column means
  height = colMeans(dat[1:4]),
  # use the c() function to provide labels for the bars
  names.arg = c("Sepal Length", "Sepal Width", "Petal Length", "Petal Width"),
  # use the c() function to pick colors for the bars
  col = c("blue","blue","red","red"),
  main = "This is my main title", 
  xlab = "This is my x-axis label", 
  ylab = "This is my y-axis label",
)
 
## boxplot ----
?boxplot

boxplot( 
  dat$Sepal.Length ~ dat$Species,
  col = "purple",
  main = "This is my main title", 
  sub = "", 
  xlab = "Species", 
  ylab = "Sepal.Length (cm)",
  notch = TRUE
  )


## x - y plot ----



# remember that you don't have to pick the same 
# things to graph that I did in class use the 
# pairs function to explore the data and pick
# two interesting variables to plot.

pairs(dat)

?plot

plot(x = dat$Petal.Width, y = dat$Petal.Length,
     col = dat$Species,
     pch = 19,
     main = "Petal Width and Length are Proportional", 
     sub = "Dot color represents different species", 
     xlab = "Petal Width (cm)", 
     ylab = "Petal Length (cm)",
     cex = 2)





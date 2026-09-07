# Before you do anything else in R we need to make sure 
# you have the files required for this lab and
# tell R which directory the files are in (called setting your working directory).
# We have written code which does this for you but you need to run it.
# To run it follow the directions on line 7, 12 and 14 to 16.

# <---------    Highlight from this line...
download.file(url = "https://github.com/JakeSaunders/mri.lab/raw/master/MRIlab.zip",
              destfile = "MRIlab.zip")
unzip(zipfile = "MRIlab.zip",exdir = "MRIlab")
setwd("MRIlab/")
# <---------    to this line.

# 1) Press and hold CTRL, 
# 2) Pressing and release ENTER once
# 3) Release CTRL

# That is how you run code in R, we will go over it again 
# because you will be doing it a lot. 
# Now you are able to start the exercise.

#### Part 1) R, the very basics ----

# This window or panel is the R Script Editor, this is where you write code to be run later
# You run code from this window by highlighting it and pressing CTRL+ENTER.
# Highlight the code on the line below and run it by highlighting & pressing CTRL+ENTER.

10+7

# Now look in the window below this one to find the console.
# The code you ran appeared in the console with the calculated answer.
# You can do simple and complex math in R using of +, -, * , / and ^ 
# to do addition, subtraction, multiplication, division and calculate exponents. 
# Write three to five equations below using any combination of these operators 
# and run the code by highlighting & pressing CTRL+ENTER.





#1) Copy your input and output from the terminal below into the answer sheet
# or word document. 

# When coding you often want to save a number or series of numbers as a descriptive variable.
# In R we assign numbers to a variable using the assignment operator <- like so:

x <- 17

# When you highlight and run this line of code you will notice that x appears in the 
# list of variables in the Environment panel in the upper right, along with a short description
# now run the code below:

x

# you can do math with a variable just like you could a number, as in the code below. 

17^x

# you can also save a string of characters as a variable, as below

txt <- "This is my text variable"

# and the entire string of characters will be returned if you run the variable alone 

txt

# R is really good at storing and processing lists of things. 
# To make a list (or as R calls it a Vector) of numbers use the first and last number in
# a list separated by a : like so

1:100

# You can also make a vector of non-consecutive numbers by 
# putting them inside the combine or concatenate function c(...) 
# The three dots is just short hand for a list of things separated by commas like so:

c(1,7,9,100,80,280,2,59)

# Combine is one of the simplest commands or functions in R, and just like every command or function 
# in R you can get a help with it by putting a ? in front of it in running it like so:

?c

# Look to the bottom right panel for the details or manual on how to use the c() function. 
# Read through the c( ) documentation and see if you have any questions. As you can see the
# information is very detailed but all the R documentation entries have the same format so
# the secret to figuring out how to do new things in R is often to read the documentation.


# As you assigned a single number of character string to a variable. 
# You can also assign a list to a variable:

my.numbers <- c(1,7,9,100,80,280,2,59)

my.numbers

my.strings <- c("You","Can","Even","Make","A","List","of","strings")

my.other.strings <- c("These","Are","My", "Other", "Strings")

my.strings

my.other.strings

# If a list is made up of all numeric elements you can do math with the whole list, as below:

my.numbers/10

# The result of which is to divide every number in the list by 10.

# You can even do math with two numeric lists, as below:


my.numbers/10:17



# The result of this code is to divide every number in the 
# my.numbers variable with the corresponding number in the
# list of numbers 10, 11, 12, 13, 14, 15, 16, 17.
# Check a couple of these by hand so you know I am correct.

my.numbers/my.numbers
# this code divides the list by the list resulting in a list of 1's

# You can also assign the output of code to a new variable

my.new.numbers <- my.numbers/10:17

my.new.numbers

# Make sure the value of the thing on the right is what you want
# the value of the thing on the left to become. The assignment 
# operator ( <- ) is an arrow because of this reason. The equals sign
# also works this way but it is also used for other things in R, so 
# it is best if you don't use the equals sign for both things. 

# In R there are also commands or functions that can process character data. 
# One of the simples functions is paste(...)
# which pastes together any two things passed to it.

paste("Two","Things")

# paste can also work on a list

paste(my.strings, "Adding This")

# and on two lists

paste(my.strings, my.other.strings)

# Some functions only work on one type of information (e.g. numbers or strings) 
# and will try to force or coerce one type of information into the other to make the
# command work. For example if you give paste numeric information it forces it into a string

paste(1,2)

paste( my.strings, my.numbers)


# Some functions will have options that you can pick in addition to the data you give them.
# These are called arguments and can be found in the help page for a function. 
# Write a line of code to look up the help (?) for the paste function to see the arguments
# and what they do.



# Then after reading the help document use paste(), my.strings, my.numbers and 
# one of the arguments to write a line of code to produce the following output:

"You.1-Can.7-Even.9-Make.100-A.80-List.280-of.2-strings.59"

# and paste this code in the worksheet under entry 2. 




# This best way to approach this problem is to 
# write code you think might work and run it to
# see what's the result. You can't break anything
# by doing this, the worse that will happen is 
# that you will see a red error message in the console
# if things don't work. If you can't figure out the 
# answer with the help document, try looking for the 
# answer using http://www.google.com or 
# http://r.789695.n4.nabble.com/ or 
# https://stackoverflow.com/ before 
# asking for help from the instructor.


#### Part 2) Working with data sets in R ----

# Packages are groups of R functions or commands people have written 
# to analyze or work with particular kinds of data. Packages are 
# loaded into R with the library function.

# The next line of code loads a package with sample data sets. 
library("datasets")

# You have been working with excel spreadsheets 
# to store data and R uses a similar format called 
# a data.frame or a matrix. A data.frame is just a 
# list of lists that are combined together into rows and columns. 

# You can view a data.frame from the datasets library by running the code below

DNase

# To get the size of the data frame you can use the function dim() to return
# a list with the number of rows and the number of columns of the data frame.

dim(DNase)

# Now write two lines of code.
# The first should assign the DNase to a new variable name. 
# The second should use the dim() function to check that 
# it is the same size as DNase. Paste both the lines 
# into answer #3 on the worksheet. 



# There are two ways to select different parts of a data frame.
# One way selects the numbers of the columns or rows using a 
# sequence of numbers like
# DNase[ numeric-list-of-rows, numeric-list-of-columns]
# The following code selects row 20 to row 30 of the DNase dataframe.

DNase[20:30, ]

# You can also subset by columns.
# The following code selects column 2 and 3

DNase[ , 2:3]

# You can also subset by both rows and columns like so:
DNase[40:45, 2:3]

# Before you move on, notice the commas in the above lines of code.
# Comma placement is very important in R, as you might have already realized
# from the exercise with the paste() function. A good first step when
# you are troubleshooting things in R is to double check your commas. 

# While comma placement is very important, spaces and even line returns (e.g. "enters") 
# are generally not important. For example run the following pieces of code and notice 
# that each gives you the same results. 

DNase[40:45, 2:3]

      DNase[              40:45      ,      2:3       ]

DNase[
    40:45,
    2:3
]

# You can use the white space to make your code more readable when you are writing it.


# Now write a line of code that selects the last 5 rows from the DNAase data.frame and 
# paste the output it produces into answer #4 of the worksheet.




# The elements of lists can have names, and often data.frames rows and columns are named.
# Look up the function names() and write a line of code using it to see the names of 
# the columns of the DNase object. 



# These names can make it You can also select a column by name using the $ like so...

DNase$conc

# this code returns all the numbers in the column "conc," now you write a line of code 
# that returns all the numbers in the "density" column



# This nomenclature makes it really easy to make graphs from data frames 
# with R. Much of the power and usefulness of R comes in its ability 
# to analyze data and produce graphs quickly. For example, 
# if I wanted a graph summarizing the data in DNase I could 
# do it in one line of code. 

plot(x = DNase$conc, y = DNase$density)

# Notice how I used the $ to select individual columns of
# the data set. Remembering that will make the last part 
# of this exercise easier. 

# From here I could make it a publication quality graph 
# with a title and more descriptive X and Y axes labels
# by looking at ?plot and filling in the other arguments below. 
# (assume concentration is in mM and density is lumens)
# If you want to know more about how the DNase data set used
# for this step was generated check using the help 
# feature (e.g. ?DNase ).

# For more information on scatterplots see: 
# https://www.statmethods.net/graphs/scatterplot.html

# Once you are done, paste your updated code into #5 on the worksheet
plot(
  DNase$density ~ DNase$conc,
  main = "",
  xlab = "",
  ylab = ""
)

# Then just above the graph in the bottom right panel 
# click on Export > Save as Image to save the graph 
# you generated as a JPEG file and paste 
# it under your code in #5 on the work sheet.


# Before we start the next section we what to clear our work space of the variables
# we have assigned. To see all the variables currently assigned use the ls() function
# or look in the "Environment" table in the upper right. You can even hover on these 
# variables to see details about them.

ls()

# To remove or delete the variables we use the rm() function. Write a line of code 
# looking up the manual or help for the rm() function and then write a second line of code
# to remove the variables currently assigned. Paste this line of code into #6 on the worksheet.


#### Part 3) Get some practice with R ----

# One of the reasons so many scientists use R is that it's free and
# people have written groups of small programs or functions to deal with 
# all different types of data or to accomplish many different tasks.
# These groups of programs are call packages. 

# Several of the most used packages are loaded when you install R and you can
# view a list them by looking under the "Packages" tab in the bottom right.
# Take a look.

# But not everybody needs every package so specialized packages sometimes have
# to be installed before they can be used. And that is accomplished with a function
# install.packages(...)

#There is even a package called swirl that teaches you how to use R (https://swirlstats.com). 
# install it by runing the code below:

install.packages("swirl")
install.packages("jpeg")

# After the package is installed, you then have to tell R that you are going to be using it.
# this is called loading a package and it is accomplished with the library function. 
# run the code below to load the swirl. 

library("swirl")
# Note: Run the next two lines, then follow the instructions in the console

# then we are going to install some specific exercises in swirl by running the following code

install_course("R Programming")
install_course("Exploratory_Data_Analysis")

# wait for the courses to downloads then...

# Running this following command should start instructions in the console (bottom left panel). 

swirl()

# Follow these red text instructions in the console until it asks which course you would 
# like to choose then select the number for "R Programming" 
# Then complete the following lessons 
    # 7: Matrices and Data Frames    
    # 8: Logic
    # 15: Base Graphics

# After completing those lessons, select "Exploratory Data Analysis" within swirl and 
# complete "5: Base Plotting System" 

# Here are some specific commands that you might want to remember while in swirl:
# | When you are at the R prompt (>):
# | -- Typing skip() allows you to skip the current question.
# | -- Typing play() lets you experiment with R on your own; swirl will ignore what you do...
# | -- UNTIL you type nxt() which will regain swirl's attention.
# | -- Typing bye() causes swirl to exit. Your progress will be saved.
# | -- Typing main() returns you to swirl's main menu.
# | -- Typing info() displays these options again.


#### Part 4) Plotting your own data ----

# Once you are done with these tutorials you should have a 
# good idea of how to make graphs of your own data in R. 

# Run the function below to download a data set from the experiments last week
download.file(
    url = "https://raw.githubusercontent.com/JakeSaunders/mri.lab/master/SensorimotorLab.csv",
    destfile = "SensorimotorLab.csv"
    )

# ...And read the file into R with the following code

lab <- read.csv("SensorimotorLab.csv")

# look at the data you assigned to the variable lab



# Now using what you know write code to make 3 graphs 
# Save those graphs using the export command in the "Plots" panel
# You might want to do some research to figure out how to add some color 
# to your graphs to make them more interesting and unique
    
    # A boxplot showing the relationship between Trial and Time




    # A boxplot showing the relationship between Trial and Errors




    # A scatterplot showing the relationship between Time and Errors







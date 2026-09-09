#### [Optional] If you would like to get some more practice with R ----

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

# Running the following line of code should start instructions in the console (bottom left panel). 
# From here on out type your responses directly in the console.
swirl()

# Follow these red text instructions in the console until it asks which course you would 
# like to choose then select the number for "R Programming" 
# You can do the lessons in order or skip around. If you skip around I recommend the following courses:
# 2: Workspace and Files
# 7: Matrices and Data Frames    
# 8: Logic
# 12: Looking at Data 
# 15: Base Graphics

# Here are some specific commands that you might want to remember while in swirl:
#  When you are at the R prompt (>):
#  -- Typing skip() allows you to skip the current question.
#  -- Typing play() lets you experiment with R on your own; swirl will ignore what you do...
#  -- UNTIL you type nxt() which will regain swirl's attention.
#  -- Typing bye() causes swirl to exit. Your progress will be saved.
#  -- Typing main() returns you to swirl's main menu.
#  -- Typing info() displays these options again.

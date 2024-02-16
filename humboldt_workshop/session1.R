##### 1.1 Introduction
# R is an interpreted programming language, meaning that input to a terminal is
# evaluated directly. This is different from other programming languages, such
# as C or Java, where code must be compiled before running. Eliminating the
# preliminary step of compilation allows for R to be run interactively, which
# has several advantages. For example, interactive mode allows us to manipulate
# and run calculations on our data in real time, displaying results
# (somewhat) instantly
#
# To demonstrate the concept of the interactivity in R, run the following code:
# This can either copy and pasted or typed directly into an interactive R
# terminal. Press enter/return to evaluate the line.

2 + 2

# The output should be something along the lines of
# [1] 4
# Here, the bracketed 1, [1] is a vector index. Our output vector is length 1

##### 1.2 Variables in R and Vectorization
# If you've used other programming languages, you may be familiar with "=" as
# the usual assignment operator, however, it conventional in R to use "<-"
myVar <- 3

# The numeric value 3 is now stored in "myVar" we can type the variable name
# into the terminal to see it's values. 
myVar

# There are other types of variables, besides numberic variables.
myCharacter <- "Hello World"
myNumericVector <- 1:100
myFactor <- factor("cat", "dog", "dog")
myMatrix <- matrix(1:9, nrow=3, ncol=3, byrow=TRUE)
myComplex <- 1 + 1i

# It is important to know that R is a "loosely typed" language, which means
# that the programmer does not neccesarily tell R what class the variable is
# going to be. Lower level languages require this. If you want to keep track
# of the class, or "type" of a variable, you can use
class(myVar)

# One of R's most powerful features is that it is vectorized. This means that
# R readily inerprets operations on vectors. For example
1:10 #a simple numeric vector

3 + 1:10 #adds three to each element of the vector

# An operation like this in a lower level programming language would usually
# require a for loop. This is not only more concise but also saves valuable
# programming time

## 1.2 Exercises
# 1.2.1 R has a built in dataset called iris. Determine the class of this 
#   object
# 1.2.2 Determine the class of every column of the aforementioned dataset

##### 1.3 Functions and Arguments
# Functions an important part of any programming language. The factor, matrix,
# and class functions have already been used so far. Functions can have input
# and output. In R, the input to a function is called an argument and it is 
# placed inside a pair of parenthesis next to the function name. For example:

print("Hello World!")

# Takes the character string "Hello World!" and prints it to the screen. What
# is the input, what is the output (if any)? 
#
# Some functions take multiple arguments, like the simple c function, the
# purpose of this function is to "c"oncatenate multiple variables into a vector\
c(1, 2, 7)

# Keep in mind the output of a function can be saved into a variable
X <- c(3, 1, 4, 1, 5)
X
Y <- c(3, 1, 4, 1, 5)
Y

# We can make our own functions as well. This will save us from copy and 
# pasting and make our lives vastly more efficient. Function are a powerful
# tool, they are hard not to use. Here is function that evaluate a simple
# second degree polynomial. The function is saved into a variable.
myFunc <- function(x, A, B, C) {
    return(A*x^2 + B*x + C)
}

# One line function can be made more succinct be removing the return statement
# this will come in handy later in life when we learn about apply statements
anotherFunc <- function(x, A, B, C) A*x^2 + B*x + C

## 1.3 Exercises
# 1.3.1 Quiz yourself on the input(s) and output(s) to the above functions
# 1.3.2 Use if statements to write a function that tells you congratulations
#   if its first argument is divisible by 7, unless it's second argument is
#   FALSE. The second argument should default to TRUE, however.

##### 1.4 Libraries in R
# You may have already done something like install.packages("something"). This
# is installing a library that has useful functions, data structures, and 
# dataset for you to use. Try installing the package colorspace with the 
# following command. You will be prompted to select a CRAN mirror, pick the 
# closest one to you to optimize download speed
install.packages("colorspace")

# Now that the package is installed, it must be loaded to be used
library(colorspace)

# Now you can access the function rainbow_hcl from the colorspace library
# you wouldn't have been able to before. Try using the following (it wouldn't)
# have been possible without first loading the library
rainbow_hcl(4)

## 1.4 Excerises
# 1.4.1 Discuss the following three packages: tuneR, seewave, and soundecology
#   why might these be useful? Do you have them installed?

##### * The following sections are adapted from 
##### *     Soundscape Ecology R Course by Amandine Gasc

##### 1.5 Read and play a wav file 

# Before we begin this section, make sure you have the following files in your
# workshop folder

# ACRRIS.wav and ACRRIS.mp3: song of the marsh warbler, 
#   Acrocephalus palustris (Deroussen and Jiguet,2006)
# ANTTRI.wav: song of the tree pipit, 
#   Anthus trivialis (Deroussen and Jiguet, 2006)

# Load a sound file from your computer with R: 
# First, open a file by indicating the working directory; in our
# case the directory where the recording is located. To do this, use the
# function setwd( ), this will depend your personal computer setup
setwd("C:/Users/jvanscha/Desktop/workshop/") #windows
setwd("/home/jack/costa_rica/workshop") #linux


#Next, use the function readWave( ) to open a .wav file or readMP3( ) to open
# an .mp3 file. Assign it to an object which will be easily manipulated later.
# You can name this object as you like, for example: ObjectWave
readWave("ACRRIS.wav")
ObjectWave <- readWave("ACRRIS.wav")
ObjectMP3 <- readMP3("ACRRIS.wav")

# Most of the time, you do not want to analyze the entire audio file. It is 
# possible to select a part of the file using arguments of the function
# readWave(). In this example we will select only the first 20 seconds of
# the file
ObjectWave <- readWave("ACRRIS.wav", from=0, to=20, units="seconds")

# Check the information associated with this object by simply typing it's name
# into the console
ObjectWave

## 1.5 Excerises
# 1.5.1 Read a wave file from the Costa Rica dataset on the flashdrive
# 1.5.2 Do the above but try playing around with different readWave parameters
#   such as wl (window length

###### 1.6 Manipulate a wav file
# The oscillogram, the spectrogram and the mean spectrum are the three sound
# objects usually manipulated. In detail:
# - Oscillogram: a 2-dimensional representation of time and amplitude 
# - Spectrogram, a 3-dimensional representation of time, frequency, and amp 
# - Mean Spectrum, a 2-dimensional representation of amplitude and frequency
#
# To obtain these three objects from a .wav file, we will use different R 
# functions that will create a table with the value of the object and a 
# graphical representation of the object. In this example we will use {seewave}
# to work with a  sound of  the neotropical sparrow Zonotrichia Capensis,
# recorded by Thierry Aubin. To load this sound, use the function data().
# This function loads data that is built into an R package. To create the 
# oscillogram, use the funciton oscillo(). 

data(tico)
Oscillogram <- oscillo(tico)
Oscillogram

# To create the spectrogram use the function spectro(). 
spectrogram1 <- spectro(tico)
spectrogram1

# Finally, use the function meansepc() to create the mean spectrum. 
mean_frequency_spectrum <- meanspec(tico)
mean_frequency_spectrum

#Note that the argument osc=TRUE in the function spectr() allows you to draw 
# both the spectrogram and the oscillogram. 
spectrogram2 <- spectro(tico, osc=TRUE)
spectrogram2

#In case of the error "margin too large" do not close the graphical window but
# increase the size of the graphical window and run the code again. If you have
# any other graphical problem, close the previous graphical window and run the
# code again

# It is also possible to splice together peices of a wave file using pastew 
# and cutw. 
acrris_paste <- pastew(
    cutw(acrris, from=0, to=4,output="Wave"), 
    cutw(acrris, from =10, to = 12,output="Wave"),
    output="Wave")
    
## 1.6 Exercises
# 1.6.1 Splice together two costa rice files from different years
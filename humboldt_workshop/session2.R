###### 2.1 Loading Test Data
# Before we jump into massive amounts of calculations, lets first load a
# practice sound file from the library soundecology

library(soundecology)
library(seewave)
library(tuneR)
data(tropicalsound)

# We can see that the data we have loaded is a 20 second wave files

tropicalsound

# Looking at the spectrogram and oscillogram tell us abaout the structure of
# the sound file. We can see some frequency bands that potentially represent
# three different species

spectro(tropicalsound, osc=TRUE)

# We can also look back at the mean spectrum
meanspec(tropicalsound)

##### 2.2 Running an index: Number of Peaks
# Now we can start analyzing and "taking measurements" on our sound files. We
# start with an index, Number of peaks, as proposed by Gasc et al 20XX. This 
# index calculates the number of mean spectrum peaks above a certain threshold
# Different indices have different inputs, some require wav file objects, while
# others may require spectrograms or mean spectrums. We use our tropicalsound
# test data, as before
ms <- meanspec(tropicalsound)

# Now we run the index function, fpeaks
fpeaks(ms, amp=c(1/90,1/90))

# run help to get an idea for the different parameters
?fpeaks

# We see that amp is a parameter that controls how steep our peaks need to be
# Try changing this parameter and running fpeaks again
fpeaks(ms, amp=c(1/60,1/60))

# Setting the plot=FALSE option, we are left with a matrix of the peak
# locations. 
fpeaks(ms, amp=c(1/60,1/60), plot=FALSE)

# If all we want is the number of peaks, we can get the number of rows of the
# output matrix of fpeaks. This might be useful if we need to programmatically
# run this index on a large range of file, and don't want to be cluttered
# with a bunch of excess plots
NP_index <- nrow(fpeaks(ms,amp=c(1/90,1/90), plot=FALSE))

## 2.2 Exercises
# 2.2.1 Try changing some of the parameters of the fpeaks function, and look at
#   the resulting plots. What might be the ecological significance of changing
#   these parameters.

###### 2.3 Trying other indices
# Let's calculate the ADI (Acoustic Diversity Index)
ADI <- acoustic_diversity(tropicalsound)
ADI

# If you recall from earlier, the tropicalsound file is a mono file. If stereo
# ADI would return values for both channels, here it just stores the mono file
# index value where the left channel output would be.
#
# Notice that this time we are useing the wav file as input as opposed to then
# meanspectrum
ADI_index <- acoustic_diversity(tropicalsound)$adi_left

# We then do something similar with another index, ACI
ACI_index <- acoustic_complexity(tropicalsound)$AciTotAll_left

# The H index has two components. It is the multiplicaiton of the shannon
# entropy of the meansecptrum and the amplitude envelope. We can calculate
# each of these seperately
Hf <- sh(ms)
Ht <- th(env(tropicalsound, plot=FALSE))
H_index <- Hf * Ht

# One more index, acoustic evenness. Similar to ADI and ACI
AEI_index <- acoustic_evenness(tropicalsound)$aei_left

#We can save all these indices into a vector
vector_indices <- c(NP_index,ADI_index,ACI_index,H_index)
vector_indices

## 2.3 Exercises
# 2.3.1 Take several non overlapping subsamples from the tropicalsound wav file
#   or costa rice data file of your chosing. Recalculate the indices and see 
#   how they compare
# 2.3.2 Load few files from the Costa Rica dataseta and calculate the ACI 
#   index for file. How do that data values compares?

##### 2.4 Apply Statements and loops in R
# While lower level languages tend to use looping statements for iterative
# procedures, R offers what are called apply statements. These statements take
# the place of for loops. In R, apply statements are actually reccomnded over
# for loops, because they are faster and follow the vectorized paradigm
#
# Apply statements are functions themselves, whose arguments are a vector
# and a function. The function is applied to elements of the vector in varying
# ways. For example lapply applies the function to each element and returns
# a list
lapply(1:10, function(x) x^2)

# with sapply, R does automatic coercion on output data types, but is otherwise
# similar to lapply in function
sapply(c("Hello", "World"), function(word) paste(word, "!", sep=""))

# mapply is multivariate apply and will vectorize the nth element of all of its
# input vectorize, and apply the function to the resulting vector
A <- c(1, 1, 1)
B <- c(1, 2, 2)
C <- c(1, 1, 3)
mapply(sum, A, B, C)

# Just so you know how, here is how you would do a for loop in R. Notice how it
# is not vectorized.
for (x in 1:10) {
    print(x^2)
}

## 2.4 Exercises
# 2.4.1 Use an sapply statement to count the number of integers between 1 and
#   1000 inclusive that are divisible by both 3 and 5
# 2.4.2 Learn to use tapply. Use tapply to find the mean sepal width of each
#   unique species in the iris dataset
# 2.4.3 (Challenge) Use an sapply function to transpose a matrix

###### 2.5 Coding an Index from scratch
# We will be coding the AcouOccupancy index. It is important to understand the
# inner working of these indices so we know why they might be ecologically 
# significant. AcouOccupancy is the proportion of points in time that have a
# spectral amplitude over a certain threshold in at least on frequency bin.
# Why might this be ecologically significant.
#
# We start by generating a spectrogram, letting the spectrogram function know
# the sampling rate of the recording, that we don't want a plot, and that we
# want to use the "A" decibal weighting (for consistency sake)
sp <- spectro(tropicalsound, f=22050, plot=FALSE, dB="A")

# notice the spectrogram object spec, has multiple components. The amplitude
# values are what we are interested in
names(sp)
dim(sp$amp)

# First, lets take a look at the frequency bins we are working with
sp$freq

# Maybe, we are only interested in frequency bins between 300 Hz and 8 kHz. 
# Lets grab the vector indices of these frqeuncy bins
inds <- which(sp$freq > .3 & sp$freq < 8)

# Now we reduce the spectrogram matrix to these frequency bins. This is
# similar to applying a filter
specM <- sp$amp[inds,]

# Now we count the number of points in time where at least spectrum is above 
# a certain threshold, here we use -50
count <- sapply(1:ncol(specM), function(i) any(specM[,i] > -30))

# Now all that's left to do is do the propotion of these counts throughout
# the whole file
acouoc <- sum(count)/ncol(specM)
acouoc

# Let's complete the process by putting this in a function for easy reuse
AcouOc <- function(wav,f=44100, dBtype="A", minkHz=0.3, maxkHz=8, th=-30) {
    sp <- spectro(wav, f=f, plot=FALSE, dB=dBtype)
    M <- sp$amp[which(sp$freq > minkHz & sp$freq < maxkHz),]
    count <- sapply(1:ncol(M), function(i) any(M[,i] > th))
    return(sum(count)/ncol(M))
}

# Now compare with our result from earlier
AcouOc(tropicalsound, f=22050)
acouoc

## 2.5 Exercises
# 2.5.1 Using sapply and your AcouOc function, calculate the AcouOc on every 
#   file the in the Costa Rica datset. Also calculate ACI and H. Use R's
#   built in tools such as the hist() function to compare the distributions
#   of index values for the different years

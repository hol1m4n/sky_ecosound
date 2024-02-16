
##### 3.1 Introduction to GGPLOT
# ggplot is a powerful R library for data visualization. The idea behind ggplot
# is to generate a plot based on 3 things: a dataset, visual representations of
# datapoints, and a coordinate system 

# load ggplot with 
library(ggplot2)

# we will also load the colorspace library for use in 3.7
library(colorspace)

## 3.1 Excercises
# 3.2.1 Find the documentation for ggplot
# 3.2.2 Take a look at the ggplot cheat sheet at
#   https://www.rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf

##### 3.2 A First Plot
# For ggplot we need a dataset. R has several built in datasets, one of which 
# is stored in the variable iris. This variable is already loaded into your
# R terminal. We can pass it to some function to better understand what is 
# going on in the dataset
class(iris)

# iris is a data frame, which is the usual class of data input for ggplot
head(iris)

# Here we can see variable names and the first six values for each variable
# lets try some plotting with the Sepal.Length and Sepal.Width variables
ggplot(iris, aes(x=Sepal.Length, y=Sepal.Width)) + geom_point()

# Let's breakdown the components:
#     ggplot(...) --  this is simply the function that returns a plot object
#     iris -- this argument corresponds to the dataset we are using
#     aes(...) -- creates an aesthetic object for plotting
#     x=Sepal.Length -- tells the aes object we want Sepal Length horizontal
#     y=Sepal.Width -- and we want Sepal Width verical
#     + -- we can alter the plot by adding outputs of other ggplot functions
#     geom_point() -- plots points at the given x-y coordinates
#
# We can make a slight alteration to jitter the data. This attempts to make
# overlapping datapoint more visible by adding slight noise to the data
ggplot(iris, aes(x=Sepal.Length, y=Sepal.Width)) + geom_jitter()

# We can take this a step further and add color to the points. We will color 
# each point by it's corresponding species
ggplot(iris, aes(x=Sepal.Length, y=Sepal.Width, color=Species)) + geom_jitter()

# That's a pretty cool first plot! Hopefully it can give us some insight 
# into the data

## 3.2 Exercises
# 3.3.1 Generated a colored plot similar to above, but try different variables
#   such as Petal.Length and Petal.Width. What is different?
# 3.3.2 Find a way to adjust the transparency of points on the plot? Why might
#   this be helpful for displaying data?
# 3.3.3 Use another variable other than species to color in the points
# 3.3.4 Use the facet_wrap function to make seperate plots based on species, 
#   while coloring data point color as in Ex 3.2.3

##### 3.3 Plotting Soundscape Index Data
# The methodology for plotting soundscape index data is quite similar to above.
# the only difference is that now we have more data points, more variables,
# and more options for visualization of the data
#
# Load the index data into variable "indexData" with the read.table function
# then remove any NA's (due to corrupted files)

indexData <- read.table("indexData.Rda")
indexData <- indexData[complete.cases(indexData),]

# Run functions to get a feel for the data

class(indexData)
dim(indexData)
head(indexData)

# We can see that indexData is a data frame with 8,450 observations of 19
# variables. There are columns for time, location, and several index values.
# This is the dataset from the Borneo 2014 study. Alpha index values were run
# on the first minute of every recording taken from Transect 2 during the
# Feb-Mar recording period
#
# Now we can get into some plotting with ggplot, how about some boxplots
# to summarize the data
ggplot(indexData, aes(factor(hourMin), ACI)) + geom_boxplot()

# Well that looks like something! However it's hard to see daily trends because
# the three outliars are overtaking the scaling of the entire plot. Lets try
# again, with another arguments to custumoize the scaling of the y axis. Keep
# in mind it is okay to have multi-line statements with operaters so long as
# the operators are kept at the end of lines as opposed to the beginnings
ggplot(indexData, aes(factor(hourMin), ACI)) + 
    geom_boxplot() +
    ylim(1400, 2000) 
    
# Those axis labels aren't that attractive. Let's change theme and a main title
ggplot(indexData, aes(factor(hourMin), ACI)) + 
    geom_boxplot() +
    ylim(1400, 2000) +
    xlab("Hour of Day") +
    ylab("ACI value") +
    ggtitle("Distribution of Alpha Index (ACI) at different times of day")

# One last change. lets fix the overlapping axis labels on horizontal axis and
# make the font size slightly larger
hours <- sort(unique(indexData$hourMin))
labs <- sapply(1:length(hours), function(i) 
    if((i%%4)==1){ return(hours[i]) }
    else{ return("") })
ggplot(indexData, aes(factor(hourMin), ACI)) + 
    geom_boxplot() +
    ylim(1400, 2000) +
    xlab("Hour of Day") +
    ylab("ACI value") +
    ggtitle("Distribution of ACI at different times of day") +
    scale_x_discrete(breaks=hours, labels=labs) +
    theme(text = element_text(size=16))

# You can also save your ggplot object into a variable, to use if later
myplot <- ggplot(indexData, aes(factor(hourMin), ACI)) + 
    geom_boxplot() +
    ylim(1400, 2000) +
    xlab("Hour of Day") +
    ylab("ACI value") +
    ggtitle("Distribution of ACI at different times of day") +
    scale_x_discrete(breaks=hours, labels=labs) +
    theme(text = element_text(size=16))

# Then view it by simply typing the variable name into terminal, as with 
# any other variable

myplot

## 3.3 Exercises
# 3.4.1 Try doing boxplots of a different histogram
# 3.4.2 Instead of boxplots, try with a violin plot of dotplot
    
##### 3.4 Saving plots
# It is possible to save plots to the hard disk instead of opening them 
# directly in your R terminal. This may be beneifical for several reasons.
# Maybe you are generating a large quantity of plots and don't have the time
# to screenshot and save every single one. Also, you could be running R on
# a remote server that doesn't have support for a remote window system. If so
# one can savea the plots to a file, grab it with FTP or the lik###### 2.4 Make your own indexe, and open
# locally.
#
# First check your working directory, so you know you will be saving your files
# to. This can be changed with the setwd() function
getwd()

# We start by itializing another R graphics device. Use the png function with
# its required filename argument 
png("first_plot.png")
myplot
dev.off()

# Any graphical objects that would normally be displayed onto your screen are now
# written to first_plot.png. The dev.off() command lets R know that you are done
# using the previously invoked graphics device. Check out the files
#
# That's nice. We can change the size if we want. The png function takes 
# arguments in pixels
png("second_plot.png", height=600, width=600)
myplot
dev.off()

# There are also several other formats. For example pdf. This function uses
# inches for it's width and height arguments
pdf("third_plot.pdf", height=6, width=6)
myplot
dev.off()

## 3.4 Exercises
# 3.5.1 Try saving a bmp, jpeh and tiff file as well.
# 3.5.2 Save the iris plot from earlier
# 3.5.3 (Challenge) Programmaticaly save boxplots a la 3.3 for every single 
#   index value. You shouldn't copy and paste any code
# 3.5.4 Sometimes it is useful to show multiple plots at once. Try displaying
#   multiple boxplots with the mulitplot or par function. Then save them.

##### 3.5 Temporal Trends
# Sometimes, we will use custom aesthetic that aren't passed directly to the 
# the main gglplot function. This demonstrates that ggplot isn't confined to
# simple scatterplots
indexDataSA1 <- indexData[indexData$transectName == "SA1",]

ggplot(indexDataSA1, aes(Bioac)) + 
    geom_rect(aes(xmin=day, 
        ymin=hourMin, 
        xmax=day+1, 
        ymax=hourMin+0.5, 
        color=Bioac, 
        fill=Bioac)) + 
    xlab("Day") + 
    ylab("Hour of Day") + 
    ggtitle("Bioac Index for Every T2 Recording")

# It's not centered properly, lets zoom in
ggplot(indexDataSA1, aes(Bioac)) + 
    geom_rect(aes(xmin=day, 
        ymin=hourMin, 
        xmax=day+1, 
        ymax=hourMin+0.5, 
        color=Bioac, 
        fill=Bioac)) + 
    xlab("Day") + 
    ylab("Hour of Day") + 
    ggtitle("Bioac Index for Every T2 Recording") +
    xlim(50, 66) +
    ylim(0, 24)
    
# One last thing, lets use a different color gradient
ggplot(indexDataSA1, aes(Bioac)) + 
    geom_rect(aes(xmin=day, 
        ymin=hourMin, 
        xmax=day+1, 
        ymax=hourMin+0.5, 
        color=Bioac, 
        fill=Bioac)) + 
    xlab("Day") + 
    ylab("Hour of Day") + 
    ggtitle("Bioac Index for Every T2 Recording") +
    xlim(50, 66) +
    ylim(0, 24) +
    scale_color_gradient(low="white", high="black") +
    scale_fill_gradient(low="white", high="black")

## 3.5 Excercises 
# 3.6.1 Try this again with another index
# 3.6.2 Notice that we only look at one sensor. If we loooked at all 14, then
#   we would have a problem with overlapping points. Can you come up with a way
#   to display all 13 sensors at once? Or maybe average them? Use the full
#   indexData dataset for this. You may have to somehow reduce it into another 
#   dataset

##### 3.6 Mean Spectrograms
# Here we will see how to make an Average Mean Spectrum, as seen in the borneo
# presentation. First we load the R object that has all the spectrograms we 
# to average. It is a matrix named "M". We will save an examine it

load("3D_spec_output_t2.Rdata")
t2 <- M
class(t2)
dim(t2)
head(rownames(t2))
head(colnames(t2))

# We can see that each row corresponds to a different time and each column
# corresponds to a different frequency bin. We need to convert the rownames
# into a more meaniningful format, and then extract meaningful daytime, 
# because that is what we are averaging over

t2date <- as.POSIXct(as.numeric(rownames(t2)), origin="1970-01-01")
t2hour <- as.POSIXlt(t2date)$hour + as.POSIXlt(t2date)$min/60
hours <- sort(unique(t2hour[(t2hour * 10) %% 5 == 0]))

# The image function in R will make a colorful plot given an input matrix
# of numeric values. We will be using the function, and first generating
# a placeholder matrix that we will fill up with averages
imageT2 <- matrix(NA, nrow=length(hours), ncol=256)

#now use a combination of for and sapply to do the averaging calculations
for(h in 1:length(hours)) {
    imageT2[h,] <- sapply(1:256, function(i)
        mean(t2[which(t2hour == hours[h]), i], na.rm=TRUE)
    )
}

par(cex=1.5)
image(hours, as.numeric(colnames(t2)), imageT2, 
    main="Mean Mean Spectrum T2", 
    xlab="Time of day", 
    ylab="Frequency", 
    col=heat_hcl(32))

#NOTE: This following page describes more color palettes to choose from
# http://www.inside-r.org/packages/cran/colorspace/docs/rainbow_hcl
# and this one describes default R color pallettes (my second choices)
# http://www.r-bloggers.com/color-palettes-in-r/

## 3.6 Exercises
# 3.7.1 Try changing the color pallette to one that you like. Try one from
#   both colorspace and one of the default R pallettes
# 3.7.2 You will see in workshop folder that the mean spectrums have been
#   provided for transect 3 as well. Repeat this activity with the transect 3
#   data
# 3.7.3 Instead of plotting the mean, plot another statistic such as standard
#   deviation
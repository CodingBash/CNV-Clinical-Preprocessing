install.packages("CORE")

# A table of DNA copy number gain events observed in 100 individual tumor cells
generateInputCORE <- function(){

}

# A table of chromosome boundary positions for DNA copy number analysis
generateInputBoundaries <- function(){

}

testInputCORE <- generateInputCORE()
testInputBoundaries <- generateInputBoundaries()

#Compute 3 cores and perform no randomization
#(meaningless for estimate of significance).
data(testInputCORE)
data(testInputBoundaries)
myCOREobj<-CORE(dataIn=testInputCORE,maxmark=3,nshuffle=0,
boundaries=testInputBoundaries,seedme=123)
## Not run:
#Extend this computation to a much larger number of randomizations,
#using 2 cores of a host computer.
newCOREobj<-CORE(dataIn=myCOREobj,keep=c("maxmark","seedme","boundaries"),
nshuffle=20,distrib="Rparallel",njobs=2)
#When using "Grid", make sure you have write premission to the current
#work space.
newCOREobj<-CORE(dataIn=myCOREobj,keep=c("maxmark","seedme","boundaries"),
nshuffle=20,distrib="Grid",njobs=2)
## End(Not run)



library(shapeR)
library(here)

#Example from package

## Not run: 
# This example has two sections: (1) Demonstration of how a shapeR object
# is analyzed and (2) How to create a shapeR object from an archive of
# image files.
#-----------------------------------------
# Section 1: Analyzing a shapeR object
data(shape)
#Standardize coefficients
shape = stdCoefs(shape,"pop","length_cm")
#Visualize Wavelet and Fourier coefficients
plotWavelet(shape,level=5,class.name= "pop",useStdcoef=TRUE)
plotFourier(shape,class.name= "pop",useStdcoef=TRUE)
#Examine the mean shapes
plotWaveletShape(shape, "pop",show.angle = TRUE,lwd=2,lty=1)
plotFourierShape(shape, "pop",show.angle = TRUE,lwd=2,lty=1)


#Canonical analysis
library(vegan)
cap.res = capscale(getStdWavelet(shape) ~ getMasterlist(shape)$pop)
anova(cap.res)
#Visualize the canonical scores
eig=eigenvals(cap.res,constrained=TRUE)
eig.ratio = eig/sum(eig)
cluster.plot(scores(cap.res)$sites[,1:2],getMasterlist(shape)$pop
,plotCI=TRUE
,xlab=paste("CAP1 (",round(eig.ratio[1]*100,1),"%)",sep="")
,ylab=paste("CAP2 (",round(eig.ratio[2]*100,1),"%)",sep="")
,main="Canonical clustering"
)
#Only analyze Icelandic and Norwegian samples
shape = setFilter(shape, getMasterlist(shape, useFilter = FALSE)$pop %in% c("NO","IC"))
#Classifier on standardized wavelet
lda.res.w = lda(getStdWavelet(shape),getMasterlist(shape)$pop,CV=TRUE)
ct.w = table(getMasterlist(shape)$pop,lda.res.w$class)
diag(prop.table(ct.w, 1))
# Total percent correct
sum(diag(prop.table(ct.w)))
cap.res = capscale(getStdWavelet(shape) ~ getMasterlist(shape)$pop)
anova(cap.res)
#Classifier on canoncial values
lda.res.w = lda(scores(cap.res)$sites,getMasterlist(shape)$pop,CV=TRUE)
ct.w = table(getMasterlist(shape)$pop,lda.res.w$class)
diag(prop.table(ct.w, 1))
# Total percent correct
sum(diag(prop.table(ct.w)))


#-----------------------------------------
# Section 2: Creating a shapeR object from image files
# The following example requires the user to download an archive of JPEG
# files from https://github.com/lisalibungan/shapeR/
# place the ShapeAnalysis directory inside the working directory.

shape = shapeR(project.path = here("ShapeAnalysis"),
            info.file = "FISH.csv")

shape = detect.outline(shape, threshold=0.2,write.outline.w.org = FALSE)
shape = generateShapeCoefficients(shape)
shape = enrich.master.list(shape)
## End(Not run)
## Not run: 



## End(Not run)

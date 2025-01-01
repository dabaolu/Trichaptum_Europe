#Based on script by Julia M.I. Barth

#Read the file into R
#imiss <- read.table("illumina1_hardfilter_exclude_DP3_GQ20_nomono_noindel_biallelic.imiss")
imiss <- read.table("illumina1_hardfilter_include_DP3_GQ20_nomono_noindel_biallelic.imiss")

#Add a third column containing the percentage of missing data per individual
imiss <- cbind(imiss, V3=((imiss$V2*100)/3779458)) #with cbind we combine "imiss" with a third column "V3=", 
#containing the value in the second column "imiss$V2", devided by the total amount of reads.

#Sort the individuals according to the amount of missing data
imiss <- imiss[order(imiss$V3),]

#Make a barplot of missing data
#pdf("illumina1_hardfilter_exclude_DP3_GQ20_nomono_noindel_biallelic.imiss.pdf", height=25, width=10) #save the plot in a PDF file
pdf("illumina1_hardfilter_include_DP3_GQ20_nomono_noindel_biallelic.imiss.pdf", height=25, width=10) #save the plot in a PDF file
par(mar=c(3, 5, 3, 1)) #enlargen figure margins (bottom, left, top, and right)
barplot(imiss$V3, main="Missingness per individual", horiz=T, xlim=c(0,100), names.arg=imiss$V1,
        cex.names=0.5, las=1, space=5) #main is the title, we plot the bars horizontally, the x-axis maximum is 100%,
#column names are in column 1 of imiss, the names should be in small font (cex), the orientation of lables (las), 
#and the space and width of the bars is reduced (space).
abline(v=70, col="red") #add a line at the cut-off
dev.off() #close the plot

#Inspect the plots by opening the PDF through Rstudio (File > Open file...) and return to the activity web page.
#We set a cut-off at 70% missing sites, meaning that we will remove individuals with more than 
#70% missing data.

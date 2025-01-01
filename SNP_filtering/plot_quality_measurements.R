#Read the illumina_1_FS.SOR.MQRS.RPRS.QD.MQ.DP.txt file into R, interpret "." as missing data
illumina1 <- read.table("illumina_1_FS.SOR.MQRS.RPRS.QD.MQ.DP.txt",na.strings=c("."))

#Add a header line to the data frame
colnames(illumina1) <- c("FS","SOR","MQRS","RPRS","QD","MQ","DP")

#Check if the column headers were added correctly by printing the first lines
head(illumina1)

#Plot histograms of all parameters of interest
#Before saving all plots in one PDF file, execute the plotting commands individually
#for each quality measurement to adjust the plotting area. Adjust xlim and ylim with meaningful values for plotting.

pdf("illumina1_gvcf_quality.pdf", height=25, width=10) #save them in one PDF file
par(mfrow=c(7,2)) #make a multi-paneled plotting window with 7 rows and 2 columns
hist(illumina1$FS, main="FS")
hist(illumina1$FS, xlim = c(0,100), ylim = c(0,8000), breaks=1000,  main="FS scaled")
hist(illumina1$SOR, main="SOR org.")
hist(illumina1$SOR, xlim = c(0,10), ylim = c(0,800000), breaks=1000,  main="SOR scaled")
hist(illumina1$MQ, main="MQ org.")
hist(illumina1$MQ, xlim = c(20,60), ylim = c(0,100000), breaks=1000,  main="LG5 MQ scaled")
hist(illumina1$MQRS, main="MQRankSum org.")
hist(illumina1$MQRS, xlim = c(-10,10), ylim = c(0,80000), breaks=1000,  main="MQRankSum scaled")
hist(illumina1$QD, main="QD org.")
hist(illumina1$QD, xlim = c(0,50), ylim = c(0,70000), breaks=1000,  main="QD scaled")
hist(illumina1$RPRS, main="ReadPosRankSum org.")
hist(illumina1$RPRS, xlim = c(-5,5), ylim = c(0,60000), breaks=1000,  main="ReadPosRankSum scaled")
hist(illumina1$DP, main="DP org.")
hist(illumina1$DP, xlim = c(0,5000), ylim = c(0,500000), breaks=1000,  main="DP scaled")
dev.off() #close the plot

#You may be interested in the mean depth (DP). 
mean(illumina1$DP) #Note, this is the mean depth across 75 samples

#To decide for a filtering theshold, you may be interested in the percentage of variants with an FS > 60, or an FS > 40
((sum((illumina1$FS) > 60))/nrow(illumina1))*100 #0.024% of variants have an FS > 60
((sum((illumina1$FS) > 40))/nrow(illumina1))*100 #0.068% of variants have an FS > 40

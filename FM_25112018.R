library(tidyr)
library(dplyr)
library(data.table)

file1 <- read.csv(file = "FM_PSMS_29092018.csv")

file1 <- list(file1)
raw.input <- rbindlist(file1)
raw.input<-unique(raw.input)

raw.input <- unique(raw.input)
write.csv(raw.input, file = "raw.input.fm.csv",row.names = FALSE)

Run = as.character(unique(raw.input$Spectrum.File))
length(Run)
Channel  = c("126", "127","128", "129", "130", "131" )

conditionf <- read.csv(file = "annotation_file_FM.csv")
a <- c("a","a","a","a","a","a")
condition = list(a,a,a,a,a,a,a,a,a,a)

for (j in 1:10) {
  condition1 <- as.character(conditionf[(6*j-5):(6*j),3])
  condition[[j]] = condition1
}
condition
annov <- matrix(rep("a", 534*5), ncol=5)

colnames(annov) <- c("Run", "Channel", "Mixture", "BioReplicate", "Condition")

############################################################################
for(i in 1:length(Run)){
  if (substr(Run[i], 23, 25)=="R10") {
    annov[(1+6*(i-1)):(6+6*(i-1)),1]  <- Run[i]  # set run information
    annov[(1+6*(i-1)):(6+6*(i-1)),2]  <- Channel # set channel information
    annov[(1+6*(i-1)):(6+6*(i-1)),3] <- substr(Run[i], 23, 25) #mixture information
    annov[(1+6*(i-1)):(6+6*(i-1)),4] <- paste(annov[(1+6*(i-1)):(6+6*(i-1)),3], annov[(1+6*(i-1)):(6+6*(i-1)),2], sep =".") #bioreplicate information
    annov[(1+6*(i-1)):(6+6*(i-1)),5] <- condition[[10]] #condition
  }
  else{
    annov[(1+6*(i-1)):(6+6*(i-1)),1]  <- Run[i]  # set run information
    annov[(1+6*(i-1)):(6+6*(i-1)),2]  <- Channel # set channel information
    annov[(1+6*(i-1)):(6+6*(i-1)),3] <- substr(Run[i], 23, 24) #mixture information
    annov[(1+6*(i-1)):(6+6*(i-1)),4] <- paste(annov[(1+6*(i-1)):(6+6*(i-1)),3], annov[(1+6*(i-1)):(6+6*(i-1)),2], sep =".") #bioreplicate information
    n <- as.numeric(substr(Run[i], 24, 24))
    print(n)
    
    annov[(1+6*(i-1)):(6+6*(i-1)),5] <- condition[[n]] #condition
    #annov[(1+6*(i-1)):(6+6*(i-1)),5] <- i
  }
}

################################################################################
colnames(annov) <- c("Run", "Channel", "Mixture", "BioReplicate", "Condition")
annov <- as.data.frame(annov)
write.csv(annov, file = "annov1.csv",row.names = FALSE)

library(MSstatsTMT)

processed.input <- PDtoMSstatsTMTFormat(input = raw.input, annotation = annov, fraction = TRUE,
                                        which.proteinid = "Master.Protein.Accessions", useNumProteinsColumn = TRUE,
                                        useUniquePeptide = TRUE, rmPSM_withMissing_withinRun = FALSE,
                                        rmPSM_withfewMea_withinRun = TRUE, rmProtein_with1Feature = FALSE,
                                        summaryforMultipleRows = sum)


processed.input1 <- unique(processed.input)


save(processed.input1, file = "processed.input.fm.rda1")

write.csv(processed.input1, file = "processed.input1.fm.csv",row.names = FALSE)

quant.msstats.fm <- protein.summarization(processed.input1, method = "msstats", normalization = TRUE,
                                          MBimpute = TRUE, maxQuantileforCensored = NULL)

write.csv(quant.msstats.fm, file = "Normalized_protein_abundance_fm.csv", row.names = FALSE)


quant.msstats.fm1 <- protein.summarization(processed.input1, method = "msstats", normalization = FALSE,
                                           MBimpute = TRUE, maxQuantileforCensored = NULL)

write.csv(quant.msstats.fm1, file = "Non_normalized_protein.abundance.fm.csv", row.names = FALSE)

###########################################################################
#change of data format for quantile normalization

b <- with(quant.msstats.fm,tapply(Abundance, list(Protein,Condition), mean))
b1 <- as.data.frame(b)

write.csv(b, "prequantile_normalised.csv")
b1 <- as.matrix(b)

###########################################################################
#Quantile normalization

library(preprocessCore)
c <- normalize.quantiles(b1)
c1 <- as.data.frame(c)
colnames(c1) <- colnames(b1)
rownames(c1) <- rownames(b1)
write.csv(c1, "Quantile_normalized_fm.csv")

#########################################################
#quantile Normalized value have been replaced in Normalized_protein_abundance_fm.csv

quant.msstats.fm2 <- read.csv("Quantile_normalized_protein_abundance_fm.csv")


file <- quant.msstats.fm2
k=1
i=1
while (TRUE) {
  # print(i)
  k = k+1
  print(i)
  cur_pro<- as.character(file$Protein[1+60*(i-1)])
  
  for (j in 1:60){
    a=FALSE
    nex_pro <- as.character(file$Protein[j+60*(i-1)])
    # print(j)
    # print(cur_pro)
    # print(nex_pro)
    if (is.na(nex_pro) == TRUE || nex_pro != cur_pro) {
      a= TRUE
      n1 = 1+60*(i-1)
      print(n1)
      n2=  j-1+60*(i-1)
      file <- file[-c(n1:n2), ]
      break
    }
  }
  if (a==FALSE){
    i=i+1
  }
  if (a==TRUE) {
    i=i
  }
  
  print(k)
  if (k==2500){
    break
  }
}

write.csv(file, file = "normed_proteins_selected ones_fm.csv",row.names = FALSE)

sample=file
sample=read.csv("normed_proteins_selected ones_fm.csv")
type = sample$Condition[9]
sample$Condition<- as.character(sample$Condition)
for (i in 1:length(sample$Run)){
  
  type = sample$Condition[i]
  print(type)
  if (substr(type,1,2)=="HC"){
    sample$Condition[i]= "HC"
  }
  else if (substr(type,1,4)=="NSFM"){
    sample$Condition[i]="SF"
  }
  else if (substr(type,1,3)=="CSA"){
    sample$Condition[i]="control"
  }
  else if (substr(type,1,2)=="SA"){
    sample$Condition[i]="SF"
  }
  else if (substr(type,1,2)=="SF"){
    sample$Condition[i]="SF"
  }
  else if (substr(type,1,3)=="CCB"){
    sample$Condition[i]="control"
  }
  else if (substr(type,1,2)=="CB"){
    sample$Condition[i]="SF"
  }
  else{}
}
write.csv(sample,"channel_name_changed_for_comparison.csv")



test.pd <- groupComparison.TMT(data = sample, contrast.matrix = "pairwise",
                               remove_norm_channel = TRUE, moderated = TRUE, adj.method = "BH")


write.csv(test.pd, file = "testing.result_fm_HC vs Malaria.csv", row.names = FALSE)
SignificantProteins <- test.pd[test.pd$adj.pvalue < 0.05 ,]

nrow(SignificantProteins)

####################################################################
## Profile plot
dataProcessPlots.TMT(data.psm = processed.input1,
                     data.summarization = quant.msstats.fm2,
                     type='ProfilePlot',
                     width = 21,
                     height = 7)

## QC plot
dataProcessPlots.TMT(data.psm = processed.input1,
                     data.summarization = quant.msstats.fm2,
                     type='QCPlot',
                     width = 21,
                     height = 7)



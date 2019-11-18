library(tidyr)
library(dplyr)
library(data.table)
file1 <- read.csv(file = "DF_PSMs_30092018.csv")


file1 <- list(file1)
dengue.raw.input <- rbindlist(file1)
dengue.raw.input<-unique(dengue.raw.input)


write.csv(dengue.raw.input, file = "dengue.raw.input.1.csv",row.names = FALSE)

Run = as.character(unique(dengue.raw.input$Spectrum.File))
length(Run)
Channel  = c("126", "127","128", "129", "130", "131" )

conditionf <- read.csv(file = "annotation_file_dengue.csv")
a <- c("a","a","a","a","a","a")
condition = list(a,a,a,a) #4 will be replaced by no. of reaction (in this case total 4 reactions),6 will be replaced y no. of channel. 
for (j in 1:4) {
  condition1 <- as.character(conditionf[(6*j-5):(6*j),3])
  condition[[j]] = condition1
}
condition
annov <- matrix(rep("a", 144*5), ncol=5)
colnames(annov) <- c("Run", "Channel","Condition", "Mixture" , "BioReplicate")
for(i in 1:length(Run)){
  
  annov[(1+6*(i-1)):(6+6*(i-1)),1]  <- Run[i]  # set run information
  annov[(1+6*(i-1)):(6+6*(i-1)),2]  <- Channel # set channel information
  annov[(1+6*(i-1)):(6+6*(i-1)),4] <- paste('R',substr(Run[i], 21,21),sep='') #mixture information ,21 is the postion of reaction number in spectrum.file column
  annov[(1+6*(i-1)):(6+6*(i-1)),5] <- paste(annov[(1+6*(i-1)):(6+6*(i-1)),4], annov[(1+6*(i-1)):(6+6*(i-1)),2], sep =".") #bioreplicate information
  n <- as.numeric(substr(Run[i], 21, 21))
  print(n)
  
  annov[(1+6*(i-1)):(6+6*(i-1)),3] <- condition[[n]] #condition
  #annov[(1+6*(i-1)):(6+6*(i-1)),5] <- i
}



colnames(annov) <- c("Run", "Channel","Condition" , "Mixture", "BioReplicate")
annov <- as.data.frame(annov)
write.csv(annov, file = "annov1.csv",row.names = FALSE)

library(MSstatsTMT)


processed.input <- PDtoMSstatsTMTFormat(input = dengue.raw.input, annotation = annov, fraction = TRUE,
                                        which.proteinid = "Master.Protein.Accessions", useNumProteinsColumn = TRUE,
                                        useUniquePeptide = TRUE, rmPSM_withMissing_withinRun = FALSE,
                                        rmPSM_withfewMea_withinRun = TRUE, rmProtein_with1Feature = FALSE,
                                        summaryforMultipleRows = sum)

processed.input1 <- unique(processed.input)

save(processed.input1, file = "processed.input.1.rda1")

write.csv(processed.input1, file = "processed.input1.df.csv",row.names = FALSE)

quant.msstats.df <- protein.summarization(processed.input1, method = "msstats", normalization = TRUE,
                                          MBimpute = TRUE, maxQuantileforCensored = NULL)

write.csv(quant.msstats.df, file = "Normalized_protein_abundance_df.csv", row.names = FALSE)


quant.msstats.df1 <- protein.summarization(processed.input1, method = "msstats", normalization = FALSE,
                                           MBimpute = TRUE, maxQuantileforCensored = NULL)

write.csv(quant.msstats.df1, file = "Non_Normalized_protein_abundance_df.csv", row.names = FALSE)


###########################################################################

b <- with(quant.msstats.df,tapply(Abundance, list(Protein,Condition), mean))
b1 <- as.data.frame(b)

write.csv(b, "prequantile_normalised_df.csv")
b1 <- as.matrix(b)

###########################################################################
library(preprocessCore)
c <- normalize.quantiles(b1)
c1 <- as.data.frame(c)
colnames(c1) <- colnames(b1)
rownames(c1) <- rownames(b1)
write.csv(c1, "Quantile_normalized_df.csv")

#########################################################

quant.msstats.df2 <- read.csv("Quantile_normalized_protein_abundance_df.csv")


file <- quant.msstats.df2

k=1
i=1
while (TRUE) {
  # print(i)
  k = k+1
  print(i)
  cur_pro<- as.character(file$Protein[1+24*(i-1)])
  
  for (j in 1:24){
    a=FALSE
    nex_pro <- as.character(file$Protein[j+24*(i-1)])
    # print(j)
    # print(cur_pro)
    # print(nex_pro)
    if (is.na(nex_pro) == TRUE || nex_pro != cur_pro) {
      a= TRUE
      n1 = 1+24*(i-1)
      print(n1)
      n2=  j-1+24*(i-1)
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
  if (k == 800){
    break
  }
}

write.csv(file, file = "Quantlie_normalized_proteins_selected ones_df.csv",row.names = FALSE)


sample=file
sample=read.csv("Quantlie_normalized_proteins_selected ones_df.csv")
type = sample$Condition[3]
sample$Condition<- as.character(sample$Condition)
for (i in 1:length(sample$Run)){
  
  type = sample$Condition[i]
  print(type)
  if (substr(type,1,2)=="HC"){
    sample$Condition[i]= "HC"
  }
  else if (substr(type,1,2)=="SD"){
    sample$Condition[i]="DF"
  }
  else if (substr(type,1,2)=="DF"){
    sample$Condition[i]="DF"
  }
  else{}
}
write.csv(sample,"channel_name_changed_for_comparison.csv")


test.pd <- groupComparison.TMT(data = sample, contrast.matrix = "pairwise",
                               remove_norm_channel = TRUE, moderated = TRUE, adj.method = "BH")


write.csv(test.pd, file = "testing.result_dengue.csv", row.names = FALSE)
SignificantProteins <- test.pd[test.pd$adj.pvalue < 0.05 ,]

nrow(SignificantProteins)

####################################################################
## Profile plot
dataProcessPlots.TMT(data.psm = processed.input1,
                     data.summarization = quant.msstats.df1,
                     type='ProfilePlot',
                     width = 21,
                     height = 7)

## QC plot
dataProcessPlots.TMT(data.psm = processed.input1,
                     data.summarization = quant.msstats.df1,
                     type='QCPlot',
                     width = 21,
                     height = 7)


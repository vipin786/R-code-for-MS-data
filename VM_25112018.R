library(tidyr)
library(dplyr)
library(data.table)

file1 <- read.csv(file = "VM_PSMs_29092018.csv")


file1 <- list(file1)
raw.input <- rbindlist(file1)
raw.input<-unique(raw.input)


vivax.raw.input <- unique(raw.input)
write.csv(vivax.raw.input, file = "raw.input.vivax.csv",row.names = FALSE)

Run = as.character(unique(raw.input$Spectrum.File))
length(Run)
Channel  = c("126","127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C","131")

conditionf <- read.csv(file = "annotation_file_vivax.csv")
a <- c("a","a","a","a","a","a","a","a","a","a")
condition = list(a,a,a,a)

for (j in 1:4) {
  condition1 <- as.character(conditionf[(10*j-9):(10*j),3])
  condition[[j]] = condition1
}
condition
annov <- matrix(rep("a", 250*5), ncol=5)

colnames(annov) <- c("Run", "Channel", "Mixture", "BioReplicate", "Condition")

for(i in 1:length(Run)){
  
  annov[(1+10*(i-1)):(10+10*(i-1)),1]  <- Run[i]  # set run information
  annov[(1+10*(i-1)):(10+10*(i-1)),2]  <- Channel # set channel information
  annov[(1+10*(i-1)):(10+10*(i-1)),3]  <- paste('R',substr(Run[i], 27,27),sep='') #mixture information
  annov[(1+10*(i-1)):(10+10*(i-1)),4]  <- paste(annov[(1+10*(i-1)):(10+10*(i-1)),3], annov[(1+10*(i-1)):(10+10*(i-1)),2], sep =".") #bioreplicate information
  n <- as.numeric(substr(Run[i], 27, 27))
  print(n)
  
  annov[(1+10*(i-1)):(10+10*(i-1)),5] <- condition[[n]] #condition
  #annov[(1+6*(i-1)):(6+6*(i-1)),5] <- i
}


colnames(annov) <- c("Run", "Channel", "Mixture", "BioReplicate", "Condition")
annov <- as.data.frame(annov)
write.csv(annov, file = "annov1.csv",row.names = FALSE)

library(MSstatsTMT)

processed.input <- PDtoMSstatsTMTFormat(input = vivax.raw.input, annotation = annov, fraction = TRUE,
                                        which.proteinid = "Master.Protein.Accessions", useNumProteinsColumn = TRUE,
                                        useUniquePeptide = TRUE, rmPSM_withMissing_withinRun = FALSE,
                                        rmPSM_withfewMea_withinRun = TRUE, rmProtein_with1Feature = FALSE,
                                        summaryforMultipleRows = sum)

processed.input1 <- unique(processed.input)

save(processed.input1, file = "processed.input1.vivax.rda1")

write.csv(processed.input1, file = "processed.input1.vivax1.csv",row.names = FALSE)

quant.msstats.vm <- protein.summarization(processed.input1, method = "msstats", normalization = TRUE,
                                          MBimpute = TRUE, maxQuantileforCensored = NULL)

write.csv(quant.msstats.vm, file = "Normlaized_protein_abundance_vm.csv", row.names = FALSE)

quant.msstats.vm1 <- protein.summarization(processed.input1, method = "msstats", normalization = FALSE,
                                           MBimpute = TRUE, maxQuantileforCensored = NULL)

write.csv(quant.msstats.vm1, file = "Non_Normlaized_protein_abundance_vm.csv", row.names = FALSE)



###########################################################################

b <- with(quant.msstats.vm,tapply(Abundance, list(Protein,Condition), mean))
b1 <- as.data.frame(b)

write.csv(b, "Prequantile_normalised_vm.csv")
b1 <- as.matrix(b)

###########################################################################
library(preprocessCore)
c <- normalize.quantiles(b1)
c1 <- as.data.frame(c)
colnames(c1) <- colnames(b1)
rownames(c1) <- rownames(b1)
write.csv(c1, "Quantile_normalized_vm.csv")

#########################################################

quant.msstats.vm2 <- read.csv("Quantile_normlaized_protein_abundance_VM.csv")

file <- quant.msstats.vm2
k=1
i=1
while (TRUE) {
  # print(i)
  k = k+1
  print(i)
  cur_pro <- as.character(file$Protein[1+40*(i-1)])
  
  for (j in 1:40){
    a=FALSE
    nex_pro <- as.character(file$Protein[j+40*(i-1)])
    # print(j)
    # print(cur_pro)
    # print(nex_pro)
    if (is.na(nex_pro) == TRUE || nex_pro != cur_pro) {
      a= TRUE
      n1 = 1+40*(i-1)
      print(n1)
      n2=  j-1+40*(i-1)
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
  if (k==800){
    break
  }
}

write.csv(file, file = "normed_proteins_selected ones_vm.csv",row.names = FALSE)



sample=file
sample=read.csv("normed_proteins_selected ones_vm.csv")
type = sample$Condition[3]
sample$Condition<- as.character(sample$Condition)
for (i in 1:length(sample$Run)){
  
  type = sample$Condition[i]
  print(type)
  if (substr(type,1,1)=="H"){
    sample$Condition[i]= "HC"
  }
  else if (substr(type,2,2)=="A"||substr(type,3,3)=="A"){
    sample$Condition[i]="SV"
  }
  else if (substr(type,1,2)=="SV"){
    sample$Condition[i]="SV"
  }
  else{}
}
sample
write.csv(sample,"channel_name_changed_for_comparison.csv")



test.pd <- groupComparison.TMT(data = sample, contrast.matrix = "pairwise",
                               remove_norm_channel = TRUE, moderated = TRUE, adj.method = "BH")


write.csv(test.pd, file = "testing.result_vivax_HC vs SV.csv", row.names = FALSE)
SignificantProteins <- 
  test.pd[test.pd$adj.pvalue < 0.05 ,]

nrow(SignificantProteins)

####################################################################
## Profile plot
dataProcessPlots.TMT(data.psm = processed.input1,
                     data.summarization = quant.msstats.vm2,
                     type='ProfilePlot',
                     width = 21,
                     height = 7)

## QC plot
dataProcessPlots.TMT(data.psm=processed.input1,
                     data.summarization=quant.msstats.vm2,
                     type='QCPlot',
                     width = 21,
                     height = 7)




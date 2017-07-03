library(ggplot2)


#Load 50 files
Y1 <- NULL
Y2 <- NULL
Y3 <- NULL
Y4 <- NULL
Y5 <- NULL
Y6 <- NULL
Y7 <- NULL
Y8 <- NULL
Y9 <- NULL
Y10 <- NULL
for (var in 1:10)
{
fileName <- paste("P:/SERP/Modelling/Evaluation_Modelling/HCV prison model/Results/out", var, ".csv", sep="")
tempFile <- read.csv(fileName, sep="\t",header=TRUE)

colnames(tempFile) <- c("TIME", "C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "Pmin0", "Pmin1", "Pmin2", "Pmin3", "Pmin4", "Pmin5", "Pmin6", "Pmin7", "Pmin8", "Pmin9", "Pmin10", "Pmin11", "Pmed0", "Pmed1", "Pmed2", "Pmed3", "Pmed4", "Pmed5", "Pmed6", "Pmed7", "Pmed8", "Pmed9", "Pmed10", "Pmed11", "Pmax0", "Pmax1", "Pmax2", "Pmax3", "Pmax4", "Pmax5", "Pmax6", "Pmax7", "Pmax8", "Pmax9", "Pmax10", "Pmax11", "NewCases", "PrisonHCV", "COMMUNITYHCV", "RELEASED", "ReleasedHCV", "ReleasedHCVCom", "ReleasedHCVPris","ClearHCV", "Death", "DeathHCV", "DeathHCVCom", "DeathHCVPris", "TotalPrisonPop", "E0", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18", "E19", "E20", "E21", "E22", "E23", "E24", "E25", "E26", "R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12")

Y1 <- rbind(Y1, tempFile[1:365,])
Y2 <- rbind(Y2, tempFile[366:730,])
Y3 <- rbind(Y3, tempFile[731:1095,])
Y4 <- rbind(Y4, tempFile[1096:1460,])
Y5 <- rbind(Y5, tempFile[1461:1825,])
Y6 <- rbind(Y6, tempFile[1826:2190,])
Y7 <- rbind(Y7, tempFile[2191:2555,])
Y8 <- rbind(Y8, tempFile[2556:2920,])
Y9 <- rbind(Y9, tempFile[2921:3285,])
Y10 <- rbind(Y10, tempFile[3286:3650,]) #YRS
}

#Create Total Columns for P1 P2 P3
Y1$PminTotal <- Y1$Pmin0+Y1$Pmin1+Y1$Pmin2+Y1$Pmin3+Y1$Pmin4+Y1$Pmin5+Y1$Pmin6+Y1$Pmin7+Y1$Pmin8+Y1$Pmin9+Y1$Pmin10+Y1$Pmin11
Y2$PminTotal <- Y2$Pmin0+Y2$Pmin1+Y2$Pmin2+Y2$Pmin3+Y2$Pmin4+Y2$Pmin5+Y2$Pmin6+Y2$Pmin7+Y2$Pmin8+Y2$Pmin9+Y2$Pmin10+Y2$Pmin11
Y3$PminTotal <- Y3$Pmin0+Y3$Pmin1+Y3$Pmin2+Y3$Pmin3+Y3$Pmin4+Y3$Pmin5+Y3$Pmin6+Y3$Pmin7+Y3$Pmin8+Y3$Pmin9+Y3$Pmin10+Y3$Pmin11
Y4$PminTotal <- Y4$Pmin0+Y4$Pmin1+Y4$Pmin2+Y4$Pmin3+Y4$Pmin4+Y4$Pmin5+Y4$Pmin6+Y4$Pmin7+Y4$Pmin8+Y4$Pmin9+Y4$Pmin10+Y4$Pmin11
Y5$PminTotal <- Y5$Pmin0+Y5$Pmin1+Y5$Pmin2+Y5$Pmin3+Y5$Pmin4+Y5$Pmin5+Y5$Pmin6+Y5$Pmin7+Y5$Pmin8+Y5$Pmin9+Y5$Pmin10+Y5$Pmin11
Y6$PminTotal <- Y6$Pmin0+Y6$Pmin1+Y6$Pmin2+Y6$Pmin3+Y6$Pmin4+Y6$Pmin5+Y6$Pmin6+Y6$Pmin7+Y6$Pmin8+Y6$Pmin9+Y6$Pmin10+Y6$Pmin11
Y7$PminTotal <- Y7$Pmin0+Y7$Pmin1+Y7$Pmin2+Y7$Pmin3+Y7$Pmin4+Y7$Pmin5+Y7$Pmin6+Y7$Pmin7+Y7$Pmin8+Y7$Pmin9+Y7$Pmin10+Y7$Pmin11
Y8$PminTotal <- Y8$Pmin0+Y8$Pmin1+Y8$Pmin2+Y8$Pmin3+Y8$Pmin4+Y8$Pmin5+Y8$Pmin6+Y8$Pmin7+Y8$Pmin8+Y8$Pmin9+Y8$Pmin10+Y8$Pmin11
Y9$PminTotal <- Y9$Pmin0+Y9$Pmin1+Y9$Pmin2+Y9$Pmin3+Y9$Pmin4+Y9$Pmin5+Y9$Pmin6+Y9$Pmin7+Y9$Pmin8+Y9$Pmin9+Y9$Pmin10+Y9$Pmin11
Y10$PminTotal <- Y10$Pmin0+Y10$Pmin1+Y10$Pmin2+Y10$Pmin3+Y10$Pmin4+Y10$Pmin5+Y10$Pmin6+Y10$Pmin7+Y10$Pmin8+Y10$Pmin9+Y10$Pmin10+Y10$Pmin11

Y1$PmedTotal <- Y1$Pmed0+Y1$Pmed1+Y1$Pmed2+Y1$Pmed3+Y1$Pmed4+Y1$Pmed5+Y1$Pmed6+Y1$Pmed7+Y1$Pmed8+Y1$Pmed9+Y1$Pmed10+Y1$Pmed11
Y2$PmedTotal <- Y2$Pmed0+Y2$Pmed1+Y2$Pmed2+Y2$Pmed3+Y2$Pmed4+Y2$Pmed5+Y2$Pmed6+Y2$Pmed7+Y2$Pmed8+Y2$Pmed9+Y2$Pmed10+Y2$Pmed11
Y3$PmedTotal <- Y3$Pmed0+Y3$Pmed1+Y3$Pmed2+Y3$Pmed3+Y3$Pmed4+Y3$Pmed5+Y3$Pmed6+Y3$Pmed7+Y3$Pmed8+Y3$Pmed9+Y3$Pmed10+Y3$Pmed11
Y4$PmedTotal <- Y4$Pmed0+Y4$Pmed1+Y4$Pmed2+Y4$Pmed3+Y4$Pmed4+Y4$Pmed5+Y4$Pmed6+Y4$Pmed7+Y4$Pmed8+Y4$Pmed9+Y4$Pmed10+Y4$Pmed11
Y5$PmedTotal <- Y5$Pmed0+Y5$Pmed1+Y5$Pmed2+Y5$Pmed3+Y5$Pmed4+Y5$Pmed5+Y5$Pmed6+Y5$Pmed7+Y5$Pmed8+Y5$Pmed9+Y5$Pmed10+Y5$Pmed11
Y6$PmedTotal <- Y6$Pmed0+Y6$Pmed1+Y6$Pmed2+Y6$Pmed3+Y6$Pmed4+Y6$Pmed5+Y6$Pmed6+Y6$Pmed7+Y6$Pmed8+Y6$Pmed9+Y6$Pmed10+Y6$Pmed11
Y7$PmedTotal <- Y7$Pmed0+Y7$Pmed1+Y7$Pmed2+Y7$Pmed3+Y7$Pmed4+Y7$Pmed5+Y7$Pmed6+Y7$Pmed7+Y7$Pmed8+Y7$Pmed9+Y7$Pmed10+Y7$Pmed11
Y8$PmedTotal <- Y8$Pmed0+Y8$Pmed1+Y8$Pmed2+Y8$Pmed3+Y8$Pmed4+Y8$Pmed5+Y8$Pmed6+Y8$Pmed7+Y8$Pmed8+Y8$Pmed9+Y8$Pmed10+Y8$Pmed11
Y9$PmedTotal <- Y9$Pmed0+Y9$Pmed1+Y9$Pmed2+Y9$Pmed3+Y9$Pmed4+Y9$Pmed5+Y9$Pmed6+Y9$Pmed7+Y9$Pmed8+Y9$Pmed9+Y9$Pmed10+Y9$Pmed11
Y10$PmedTotal <- Y10$Pmed0+Y10$Pmed1+Y10$Pmed2+Y10$Pmed3+Y10$Pmed4+Y10$Pmed5+Y10$Pmed6+Y10$Pmed7+Y10$Pmed8+Y10$Pmed9+Y10$Pmed10+Y10$Pmed11

Y1$PmaxTotal <- Y1$Pmax0+Y1$Pmax1+Y1$Pmax2+Y1$Pmax3+Y1$Pmax4+Y1$Pmax5+Y1$Pmax6+Y1$Pmax7+Y1$Pmax8+Y1$Pmax9+Y1$Pmax10+Y1$Pmax11
Y2$PmaxTotal <- Y2$Pmax0+Y2$Pmax1+Y2$Pmax2+Y2$Pmax3+Y2$Pmax4+Y2$Pmax5+Y2$Pmax6+Y2$Pmax7+Y2$Pmax8+Y2$Pmax9+Y2$Pmax10+Y2$Pmax11
Y3$PmaxTotal <- Y3$Pmax0+Y3$Pmax1+Y3$Pmax2+Y3$Pmax3+Y3$Pmax4+Y3$Pmax5+Y3$Pmax6+Y3$Pmax7+Y3$Pmax8+Y3$Pmax9+Y3$Pmax10+Y3$Pmax11
Y4$PmaxTotal <- Y4$Pmax0+Y4$Pmax1+Y4$Pmax2+Y4$Pmax3+Y4$Pmax4+Y4$Pmax5+Y4$Pmax6+Y4$Pmax7+Y4$Pmax8+Y4$Pmax9+Y4$Pmax10+Y4$Pmax11
Y5$PmaxTotal <- Y5$Pmax0+Y5$Pmax1+Y5$Pmax2+Y5$Pmax3+Y5$Pmax4+Y5$Pmax5+Y5$Pmax6+Y5$Pmax7+Y5$Pmax8+Y5$Pmax9+Y5$Pmax10+Y5$Pmax11
Y6$PmaxTotal <- Y6$Pmax0+Y6$Pmax1+Y6$Pmax2+Y6$Pmax3+Y6$Pmax4+Y6$Pmax5+Y6$Pmax6+Y6$Pmax7+Y6$Pmax8+Y6$Pmax9+Y6$Pmax10+Y6$Pmax11
Y7$PmaxTotal <- Y7$Pmax0+Y7$Pmax1+Y7$Pmax2+Y7$Pmax3+Y7$Pmax4+Y7$Pmax5+Y7$Pmax6+Y7$Pmax7+Y7$Pmax8+Y7$Pmax9+Y7$Pmax10+Y7$Pmax11
Y8$PmaxTotal <- Y8$Pmax0+Y8$Pmax1+Y8$Pmax2+Y8$Pmax3+Y8$Pmax4+Y8$Pmax5+Y8$Pmax6+Y8$Pmax7+Y8$Pmax8+Y8$Pmax9+Y8$Pmax10+Y8$Pmax11
Y9$PmaxTotal <- Y9$Pmax0+Y9$Pmax1+Y9$Pmax2+Y9$Pmax3+Y9$Pmax4+Y9$Pmax5+Y9$Pmax6+Y9$Pmax7+Y9$Pmax8+Y9$Pmax9+Y9$Pmax10+Y9$Pmax11
Y10$PmaxTotal <- Y10$Pmax0+Y10$Pmax1+Y10$Pmax2+Y10$Pmax3+Y10$Pmax4+Y10$Pmax5+Y10$Pmax6+Y10$Pmax7+Y10$Pmax8+Y10$Pmax9+Y10$Pmax10+Y10$Pmax11

#Summary
#Min
summary(Y1$PminTotal) #Y1
summary(Y2$PminTotal) #Y1
summary(Y3$PminTotal) #Y1
summary(Y4$PminTotal) #Y1
summary(Y5$PminTotal) #Y1
summary(Y6$PminTotal) #Y1
summary(Y7$PminTotal) #Y1
summary(Y8$PminTotal) #Y1
summary(Y9$PminTotal) #Y1
summary(Y10$PminTotal) #Y1

#Med
summary(Y1$PmedTotal) #Y1
summary(Y2$PmedTotal) #Y1
summary(Y3$PmedTotal) #Y1
summary(Y4$PmedTotal) #Y1
summary(Y5$PmedTotal) #Y1
summary(Y6$PmedTotal) #Y1
summary(Y7$PmedTotal) #Y1
summary(Y8$PmedTotal) #Y1
summary(Y9$PmedTotal) #Y1
summary(Y10$PmedTotal) #Y1

#Max
summary(Y1$PmaxTotal) #Y1
summary(Y2$PmaxTotal) #Y1
summary(Y3$PmaxTotal) #Y1
summary(Y4$PmaxTotal) #Y1
summary(Y5$PmaxTotal) #Y1
summary(Y6$PmaxTotal) #Y1
summary(Y7$PmaxTotal) #Y1
summary(Y8$PmaxTotal) #Y1
summary(Y9$PmaxTotal) #Y1
summary(Y10$PmaxTotal) #Y1

############Load 50 files from second file
CY1 <- NULL
CY2 <- NULL
CY3 <- NULL
CY4 <- NULL
CY5 <- NULL
CY6 <- NULL
CY7 <- NULL
CY8 <- NULL
CY9 <- NULL
CY10 <- NULL
for (var in 1:50)
{
fileName <- paste("~/Dropbox/HCV Model/HCV New Model/HCV Model/calibration/Final Calibration Sets/Movement Rates Sensitivity/Probability/Pure/out", var, ".csv", sep="")
tempFile <- read.csv(fileName, sep="\t",header=TRUE)

colnames(tempFile) <- c("TIME", "C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "Pmin0", "Pmin1", "Pmin2", "Pmin3", "Pmin4", "Pmin5", "Pmin6", "Pmin7", "Pmin8", "Pmin9", "Pmin10", "Pmin11", "Pmed0", "Pmed1", "Pmed2", "Pmed3", "Pmed4", "Pmed5", "Pmed6", "Pmed7", "Pmed8", "Pmed9", "Pmed10", "Pmed11", "Pmax0", "Pmax1", "Pmax2", "Pmax3", "Pmax4", "Pmax5", "Pmax6", "Pmax7", "Pmax8", "Pmax9", "Pmax10", "Pmax11", "NewCases", "PrisonHCV", "COMMUNITYHCV", "RELEASED", "ReleasedHCV", "ReleasedHCVCom", "ReleasedHCVPris","ClearHCV", "Death", "DeathHCV", "DeathHCVCom", "DeathHCVPris", "TotalPrisonPop", "E0", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18", "E19", "E20", "E21", "E22", "E23", "E24", "E25", "E26", "R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12")

CY1 <- rbind(CY1, tempFile[1:365,])
CY2 <- rbind(CY2, tempFile[366:730,])
CY3 <- rbind(CY3, tempFile[731:1095,])
CY4 <- rbind(CY4, tempFile[1096:1460,])
CY5 <- rbind(CY5, tempFile[1461:1825,])
CY6 <- rbind(CY6, tempFile[1826:2190,])
CY7 <- rbind(CY7, tempFile[2191:2555,])
CY8 <- rbind(CY8, tempFile[2556:2920,])
CY9 <- rbind(CY9, tempFile[2921:3285,])
CY10 <- rbind(CY10, tempFile[3286:3650,]) #YRS
}

#Create Total Columns for P1 P2 P3
CY1$PminTotal <- CY1$Pmin0+CY1$Pmin1+CY1$Pmin2+CY1$Pmin3+CY1$Pmin4+CY1$Pmin5+CY1$Pmin6+CY1$Pmin7+CY1$Pmin8+CY1$Pmin9+CY1$Pmin10+CY1$Pmin11
CY2$PminTotal <- CY2$Pmin0+CY2$Pmin1+CY2$Pmin2+CY2$Pmin3+CY2$Pmin4+CY2$Pmin5+CY2$Pmin6+CY2$Pmin7+CY2$Pmin8+CY2$Pmin9+CY2$Pmin10+CY2$Pmin11
CY3$PminTotal <- CY3$Pmin0+CY3$Pmin1+CY3$Pmin2+CY3$Pmin3+CY3$Pmin4+CY3$Pmin5+CY3$Pmin6+CY3$Pmin7+CY3$Pmin8+CY3$Pmin9+CY3$Pmin10+CY3$Pmin11
CY4$PminTotal <- CY4$Pmin0+CY4$Pmin1+CY4$Pmin2+CY4$Pmin3+CY4$Pmin4+CY4$Pmin5+CY4$Pmin6+CY4$Pmin7+CY4$Pmin8+CY4$Pmin9+CY4$Pmin10+CY4$Pmin11
CY5$PminTotal <- CY5$Pmin0+CY5$Pmin1+CY5$Pmin2+CY5$Pmin3+CY5$Pmin4+CY5$Pmin5+CY5$Pmin6+CY5$Pmin7+CY5$Pmin8+CY5$Pmin9+CY5$Pmin10+CY5$Pmin11
CY6$PminTotal <- CY6$Pmin0+CY6$Pmin1+CY6$Pmin2+CY6$Pmin3+CY6$Pmin4+CY6$Pmin5+CY6$Pmin6+CY6$Pmin7+CY6$Pmin8+CY6$Pmin9+CY6$Pmin10+CY6$Pmin11
CY7$PminTotal <- CY7$Pmin0+CY7$Pmin1+CY7$Pmin2+CY7$Pmin3+CY7$Pmin4+CY7$Pmin5+CY7$Pmin6+CY7$Pmin7+CY7$Pmin8+CY7$Pmin9+CY7$Pmin10+CY7$Pmin11
CY8$PminTotal <- CY8$Pmin0+CY8$Pmin1+CY8$Pmin2+CY8$Pmin3+CY8$Pmin4+CY8$Pmin5+CY8$Pmin6+CY8$Pmin7+CY8$Pmin8+CY8$Pmin9+CY8$Pmin10+CY8$Pmin11
CY9$PminTotal <- CY9$Pmin0+CY9$Pmin1+CY9$Pmin2+CY9$Pmin3+CY9$Pmin4+CY9$Pmin5+CY9$Pmin6+CY9$Pmin7+CY9$Pmin8+CY9$Pmin9+CY9$Pmin10+CY9$Pmin11
CY10$PminTotal <- CY10$Pmin0+CY10$Pmin1+CY10$Pmin2+CY10$Pmin3+CY10$Pmin4+CY10$Pmin5+CY10$Pmin6+CY10$Pmin7+CY10$Pmin8+CY10$Pmin9+CY10$Pmin10+CY10$Pmin11

CY1$PmedTotal <- CY1$Pmed0+CY1$Pmed1+CY1$Pmed2+CY1$Pmed3+CY1$Pmed4+CY1$Pmed5+CY1$Pmed6+CY1$Pmed7+CY1$Pmed8+CY1$Pmed9+CY1$Pmed10+CY1$Pmed11
CY2$PmedTotal <- CY2$Pmed0+CY2$Pmed1+CY2$Pmed2+CY2$Pmed3+CY2$Pmed4+CY2$Pmed5+CY2$Pmed6+CY2$Pmed7+CY2$Pmed8+CY2$Pmed9+CY2$Pmed10+CY2$Pmed11
CY3$PmedTotal <- CY3$Pmed0+CY3$Pmed1+CY3$Pmed2+CY3$Pmed3+CY3$Pmed4+CY3$Pmed5+CY3$Pmed6+CY3$Pmed7+CY3$Pmed8+CY3$Pmed9+CY3$Pmed10+CY3$Pmed11
CY4$PmedTotal <- CY4$Pmed0+CY4$Pmed1+CY4$Pmed2+CY4$Pmed3+CY4$Pmed4+CY4$Pmed5+CY4$Pmed6+CY4$Pmed7+CY4$Pmed8+CY4$Pmed9+CY4$Pmed10+CY4$Pmed11
CY5$PmedTotal <- CY5$Pmed0+CY5$Pmed1+CY5$Pmed2+CY5$Pmed3+CY5$Pmed4+CY5$Pmed5+CY5$Pmed6+CY5$Pmed7+CY5$Pmed8+CY5$Pmed9+CY5$Pmed10+CY5$Pmed11
CY6$PmedTotal <- CY6$Pmed0+CY6$Pmed1+CY6$Pmed2+CY6$Pmed3+CY6$Pmed4+CY6$Pmed5+CY6$Pmed6+CY6$Pmed7+CY6$Pmed8+CY6$Pmed9+CY6$Pmed10+CY6$Pmed11
CY7$PmedTotal <- CY7$Pmed0+CY7$Pmed1+CY7$Pmed2+CY7$Pmed3+CY7$Pmed4+CY7$Pmed5+CY7$Pmed6+CY7$Pmed7+CY7$Pmed8+CY7$Pmed9+CY7$Pmed10+CY7$Pmed11
CY8$PmedTotal <- CY8$Pmed0+CY8$Pmed1+CY8$Pmed2+CY8$Pmed3+CY8$Pmed4+CY8$Pmed5+CY8$Pmed6+CY8$Pmed7+CY8$Pmed8+CY8$Pmed9+CY8$Pmed10+CY8$Pmed11
CY9$PmedTotal <- CY9$Pmed0+CY9$Pmed1+CY9$Pmed2+CY9$Pmed3+CY9$Pmed4+CY9$Pmed5+CY9$Pmed6+CY9$Pmed7+CY9$Pmed8+CY9$Pmed9+CY9$Pmed10+CY9$Pmed11
CY10$PmedTotal <- CY10$Pmed0+CY10$Pmed1+CY10$Pmed2+CY10$Pmed3+CY10$Pmed4+CY10$Pmed5+CY10$Pmed6+CY10$Pmed7+CY10$Pmed8+CY10$Pmed9+CY10$Pmed10+CY10$Pmed11

CY1$PmaxTotal <- CY1$Pmax0+CY1$Pmax1+CY1$Pmax2+CY1$Pmax3+CY1$Pmax4+CY1$Pmax5+CY1$Pmax6+CY1$Pmax7+CY1$Pmax8+CY1$Pmax9+CY1$Pmax10+CY1$Pmax11
CY2$PmaxTotal <- CY2$Pmax0+CY2$Pmax1+CY2$Pmax2+CY2$Pmax3+CY2$Pmax4+CY2$Pmax5+CY2$Pmax6+CY2$Pmax7+CY2$Pmax8+CY2$Pmax9+CY2$Pmax10+CY2$Pmax11
CY3$PmaxTotal <- CY3$Pmax0+CY3$Pmax1+CY3$Pmax2+CY3$Pmax3+CY3$Pmax4+CY3$Pmax5+CY3$Pmax6+CY3$Pmax7+CY3$Pmax8+CY3$Pmax9+CY3$Pmax10+CY3$Pmax11
CY4$PmaxTotal <- CY4$Pmax0+CY4$Pmax1+CY4$Pmax2+CY4$Pmax3+CY4$Pmax4+CY4$Pmax5+CY4$Pmax6+CY4$Pmax7+CY4$Pmax8+CY4$Pmax9+CY4$Pmax10+CY4$Pmax11
CY5$PmaxTotal <- CY5$Pmax0+CY5$Pmax1+CY5$Pmax2+CY5$Pmax3+CY5$Pmax4+CY5$Pmax5+CY5$Pmax6+CY5$Pmax7+CY5$Pmax8+CY5$Pmax9+CY5$Pmax10+CY5$Pmax11
CY6$PmaxTotal <- CY6$Pmax0+CY6$Pmax1+CY6$Pmax2+CY6$Pmax3+CY6$Pmax4+CY6$Pmax5+CY6$Pmax6+CY6$Pmax7+CY6$Pmax8+CY6$Pmax9+CY6$Pmax10+CY6$Pmax11
CY7$PmaxTotal <- CY7$Pmax0+CY7$Pmax1+CY7$Pmax2+CY7$Pmax3+CY7$Pmax4+CY7$Pmax5+CY7$Pmax6+CY7$Pmax7+CY7$Pmax8+CY7$Pmax9+CY7$Pmax10+CY7$Pmax11
CY8$PmaxTotal <- CY8$Pmax0+CY8$Pmax1+CY8$Pmax2+CY8$Pmax3+CY8$Pmax4+CY8$Pmax5+CY8$Pmax6+CY8$Pmax7+CY8$Pmax8+CY8$Pmax9+CY8$Pmax10+CY8$Pmax11
CY9$PmaxTotal <- CY9$Pmax0+CY9$Pmax1+CY9$Pmax2+CY9$Pmax3+CY9$Pmax4+CY9$Pmax5+CY9$Pmax6+CY9$Pmax7+CY9$Pmax8+CY9$Pmax9+CY9$Pmax10+CY9$Pmax11
CY10$PmaxTotal <- CY10$Pmax0+CY10$Pmax1+CY10$Pmax2+CY10$Pmax3+CY10$Pmax4+CY10$Pmax5+CY10$Pmax6+CY10$Pmax7+CY10$Pmax8+CY10$Pmax9+CY10$Pmax10+CY10$Pmax11

#Summary
#Min
summary(CY1$PminTotal) #CY1
summary(CY2$PminTotal) #CY1
summary(CY3$PminTotal) #CY1
summary(CY4$PminTotal) #CY1
summary(CY5$PminTotal) #CY1
summary(CY6$PminTotal) #CY1
summary(CY7$PminTotal) #CY1
summary(CY8$PminTotal) #CY1
summary(CY9$PminTotal) #CY1
summary(CY10$PminTotal) #CY1

#Med
summary(CY1$PmedTotal) #CY1
summary(CY2$PmedTotal) #CY1
summary(CY3$PmedTotal) #CY1
summary(CY4$PmedTotal) #CY1
summary(CY5$PmedTotal) #CY1
summary(CY6$PmedTotal) #CY1
summary(CY7$PmedTotal) #CY1
summary(CY8$PmedTotal) #CY1
summary(CY9$PmedTotal) #CY1
summary(CY10$PmedTotal) #CY1

#Max
summary(CY1$PmaxTotal) #CY1
summary(CY2$PmaxTotal) #CY1
summary(CY3$PmaxTotal) #CY1
summary(CY4$PmaxTotal) #CY1
summary(CY5$PmaxTotal) #CY1
summary(CY6$PmaxTotal) #CY1
summary(CY7$PmaxTotal) #CY1
summary(CY8$PmaxTotal) #CY1
summary(CY9$PmaxTotal) #CY1
summary(CY10$PmaxTotal) #CY1



#Plots
#Population
# Define 2 vectors
DCSMin <- data.frame("Population"=c(3320, 3481, 3286, 3347, 3362, 3312, 3221, 3443, 3627))
DCSMin$floor <- 0
DCSMin$ceiling <- 0 
DCSMin$Year <- c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9") 
DCSMin$Security <- c("DCSMin")

Min <- data.frame("Population"=c(summary(Y1$PminTotal)["Mean"], summary(Y2$PminTotal)["Mean"], summary(Y3$PminTotal)["Mean"], summary(Y4$PminTotal)["Mean"], summary(Y5$PminTotal)["Mean"], summary(Y6$PminTotal)["Mean"], summary(Y7$PminTotal)["Mean"], summary(Y8$PminTotal)["Mean"], summary(Y9$PminTotal)["Mean"]))
Min$floor <- c(summary(Y1$PminTotal)["Mean"]-summary(Y1$PminTotal)["Min."], summary(Y2$PminTotal)["Mean"]-summary(Y2$PminTotal)["Min."], summary(Y3$PminTotal)["Mean"]-summary(Y3$PminTotal)["Min."], summary(Y4$PminTotal)["Mean"]-summary(Y4$PminTotal)["Min."], summary(Y5$PminTotal)["Mean"]-summary(Y5$PminTotal)["Min."], summary(Y6$PminTotal)["Mean"]-summary(Y6$PminTotal)["Min."], summary(Y7$PminTotal)["Mean"]-summary(Y7$PminTotal)["Min."], summary(Y8$PminTotal)["Mean"]-summary(Y8$PminTotal)["Min."], summary(Y9$PminTotal)["Mean"]-summary(Y9$PminTotal)["Min."])
Min$ceiling <- c(summary(Y1$PminTotal)["Max."]-summary(Y1$PminTotal)["Mean"], summary(Y2$PminTotal)["Max."]-summary(Y2$PminTotal)["Mean"], summary(Y3$PminTotal)["Max."]-summary(Y3$PminTotal)["Mean"], summary(Y4$PminTotal)["Max."]-summary(Y4$PminTotal)["Mean"], summary(Y5$PminTotal)["Max."]-summary(Y5$PminTotal)["Mean"], summary(Y6$PminTotal)["Max."]-summary(Y6$PminTotal)["Mean"], summary(Y7$PminTotal)["Max."]-summary(Y7$PminTotal)["Mean"], summary(Y8$PminTotal)["Max."]-summary(Y8$PminTotal)["Mean"], summary(Y9$PminTotal)["Max."]-summary(Y9$PminTotal)["Mean"])
Min$Year <- c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9")  
Min$Security <- c("Min (adjusted med efflux)")

#Min to compare
CMin <- data.frame("Population"=c(summary(CY1$PminTotal)["Mean"], summary(CY2$PminTotal)["Mean"], summary(CY3$PminTotal)["Mean"], summary(CY4$PminTotal)["Mean"], summary(CY5$PminTotal)["Mean"], summary(CY6$PminTotal)["Mean"], summary(CY7$PminTotal)["Mean"], summary(CY8$PminTotal)["Mean"], summary(CY9$PminTotal)["Mean"]))
CMin$floor <- c(summary(CY1$PminTotal)["Mean"]-summary(CY1$PminTotal)["Min."], summary(CY2$PminTotal)["Mean"]-summary(CY2$PminTotal)["Min."], summary(CY3$PminTotal)["Mean"]-summary(CY3$PminTotal)["Min."], summary(CY4$PminTotal)["Mean"]-summary(CY4$PminTotal)["Min."], summary(CY5$PminTotal)["Mean"]-summary(CY5$PminTotal)["Min."], summary(CY6$PminTotal)["Mean"]-summary(CY6$PminTotal)["Min."], summary(CY7$PminTotal)["Mean"]-summary(CY7$PminTotal)["Min."], summary(CY8$PminTotal)["Mean"]-summary(CY8$PminTotal)["Min."], summary(CY9$PminTotal)["Mean"]-summary(CY9$PminTotal)["Min."])
CMin$ceiling <- c(summary(CY1$PminTotal)["Max."]-summary(CY1$PminTotal)["Mean"], summary(CY2$PminTotal)["Max."]-summary(CY2$PminTotal)["Mean"], summary(CY3$PminTotal)["Max."]-summary(CY3$PminTotal)["Mean"], summary(CY4$PminTotal)["Max."]-summary(CY4$PminTotal)["Mean"], summary(CY5$PminTotal)["Max."]-summary(CY5$PminTotal)["Mean"], summary(CY6$PminTotal)["Max."]-summary(CY6$PminTotal)["Mean"], summary(CY7$PminTotal)["Max."]-summary(CY7$PminTotal)["Mean"], summary(CY8$PminTotal)["Max."]-summary(CY8$PminTotal)["Mean"], summary(CY9$PminTotal)["Max."]-summary(CY9$PminTotal)["Mean"])
CMin$Year <- c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9")  
CMin$Security <- c("Min (pure prison transition rates)")

DCSMed <- data.frame("Population"=c(2200, 2354, 2410, 2558, 2665, 2282, 1983, 1951, 2043))
DCSMed$floor <- 0
DCSMed$ceiling <- 0 
DCSMed$Year <- c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9") 
DCSMed$Security <- c("DCSMed")

Med <- data.frame("Population"=c(summary(Y1$PmedTotal)["Mean"], summary(Y2$PmedTotal)["Mean"], summary(Y3$PmedTotal)["Mean"], summary(Y4$PmedTotal)["Mean"], summary(Y5$PmedTotal)["Mean"], summary(Y6$PmedTotal)["Mean"], summary(Y7$PmedTotal)["Mean"], summary(Y8$PmedTotal)["Mean"], summary(Y9$PmedTotal)["Mean"]))
Med$floor <- c(summary(Y1$PminTotal)["Mean"]-summary(Y1$PminTotal)["Min."], summary(Y2$PmedTotal)["Mean"]-summary(Y2$PmedTotal)["Min."], summary(Y3$PmedTotal)["Mean"]-summary(Y3$PmedTotal)["Min."], summary(Y4$PmedTotal)["Mean"]-summary(Y4$PmedTotal)["Min."], summary(Y5$PmedTotal)["Mean"]-summary(Y5$PmedTotal)["Min."], summary(Y6$PmedTotal)["Mean"]-summary(Y6$PmedTotal)["Min."], summary(Y7$PmedTotal)["Mean"]-summary(Y7$PmedTotal)["Min."], summary(Y8$PmedTotal)["Mean"]-summary(Y8$PmedTotal)["Min."], summary(Y9$PmedTotal)["Mean"]-summary(Y9$PmedTotal)["Min."])
Med$ceiling <- c(summary(Y1$PminTotal)["Max."]-summary(Y1$PminTotal)["Mean"], summary(Y2$PmedTotal)["Max."]-summary(Y2$PmedTotal)["Mean"], summary(Y3$PmedTotal)["Max."]-summary(Y3$PmedTotal)["Mean"], summary(Y4$PmedTotal)["Max."]-summary(Y4$PmedTotal)["Mean"], summary(Y5$PmedTotal)["Max."]-summary(Y5$PmedTotal)["Mean"], summary(Y6$PmedTotal)["Max."]-summary(Y6$PmedTotal)["Mean"], summary(Y7$PmedTotal)["Max."]-summary(Y7$PmedTotal)["Mean"], summary(Y8$PmedTotal)["Max."]-summary(Y8$PmedTotal)["Mean"], summary(Y9$PmedTotal)["Max."]-summary(Y9$PmedTotal)["Mean"])
Med$Year <- c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9") 
Med$Security <- c("Med (adjusted med efflux)")

CMed <- data.frame("Population"=c(summary(CY1$PmedTotal)["Mean"], summary(CY2$PmedTotal)["Mean"], summary(CY3$PmedTotal)["Mean"], summary(CY4$PmedTotal)["Mean"], summary(CY5$PmedTotal)["Mean"], summary(CY6$PmedTotal)["Mean"], summary(CY7$PmedTotal)["Mean"], summary(CY8$PmedTotal)["Mean"], summary(CY9$PmedTotal)["Mean"]))
CMed$floor <- c(summary(CY1$PminTotal)["Mean"]-summary(CY1$PminTotal)["Min."], summary(CY2$PmedTotal)["Mean"]-summary(CY2$PmedTotal)["Min."], summary(CY3$PmedTotal)["Mean"]-summary(CY3$PmedTotal)["Min."], summary(CY4$PmedTotal)["Mean"]-summary(CY4$PmedTotal)["Min."], summary(CY5$PmedTotal)["Mean"]-summary(CY5$PmedTotal)["Min."], summary(CY6$PmedTotal)["Mean"]-summary(CY6$PmedTotal)["Min."], summary(CY7$PmedTotal)["Mean"]-summary(CY7$PmedTotal)["Min."], summary(CY8$PmedTotal)["Mean"]-summary(CY8$PmedTotal)["Min."], summary(CY9$PmedTotal)["Mean"]-summary(CY9$PmedTotal)["Min."])
CMed$ceiling <- c(summary(CY1$PminTotal)["Max."]-summary(CY1$PminTotal)["Mean"], summary(CY2$PminTotal)["Max."]-summary(CY2$PminTotal)["Mean"], summary(CY3$PmedTotal)["Max."]-summary(CY3$PmedTotal)["Mean"], summary(CY4$PmedTotal)["Max."]-summary(CY4$PmedTotal)["Mean"], summary(CY5$PmedTotal)["Max."]-summary(CY5$PmedTotal)["Mean"], summary(CY6$PmedTotal)["Max."]-summary(CY6$PmedTotal)["Mean"], summary(CY7$PmedTotal)["Max."]-summary(CY7$PmedTotal)["Mean"], summary(CY8$PmedTotal)["Max."]-summary(CY8$PmedTotal)["Mean"], summary(CY9$PmedTotal)["Max."]-summary(CY9$PmedTotal)["Mean"])
CMed$Year <- c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9") 
CMed$Security <- c("Med (pure prison transition rates)")

DCSMax <- data.frame("Population"=c(2880, 2974, 3313, 3627, 3480, 3634, 3627, 3793, 4179))
DCSMax$floor <- 0
DCSMax$ceiling <- 0 
DCSMax$Year <- c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9") 
DCSMax$Security <- c("DCSMax")



Max <- data.frame("Population"=c(summary(Y1$PmaxTotal)["Mean"], summary(Y2$PmaxTotal)["Mean"], summary(Y3$PmaxTotal)["Mean"], summary(Y4$PmaxTotal)["Mean"], summary(Y5$PmaxTotal)["Mean"], summary(Y6$PmaxTotal)["Mean"], summary(Y7$PmaxTotal)["Mean"], summary(Y8$PmaxTotal)["Mean"], summary(Y9$PmaxTotal)["Mean"]))
Max$floor <- c(summary(Y1$PmaxTotal)["Mean"]-summary(Y1$PmaxTotal)["Min."], summary(Y2$PmaxTotal)["Mean"]-summary(Y2$PmaxTotal)["Min."], summary(Y3$PmaxTotal)["Mean"]-summary(Y3$PmaxTotal)["Min."], summary(Y4$PmaxTotal)["Mean"]-summary(Y4$PmaxTotal)["Min."], summary(Y5$PmaxTotal)["Mean"]-summary(Y5$PmaxTotal)["Min."], summary(Y6$PmaxTotal)["Mean"]-summary(Y6$PmaxTotal)["Min."], summary(Y7$PmaxTotal)["Mean"]-summary(Y7$PmaxTotal)["Min."], summary(Y8$PmaxTotal)["Mean"]-summary(Y8$PmaxTotal)["Min."], summary(Y9$PmaxTotal)["Mean"]-summary(Y9$PmaxTotal)["Min."])
Max$ceiling <- c(summary(Y1$PmaxTotal)["Max."]-summary(Y1$PmaxTotal)["Mean"], summary(Y2$PmaxTotal)["Max."]-summary(Y2$PmaxTotal)["Mean"], summary(Y3$PmaxTotal)["Max."]-summary(Y3$PmaxTotal)["Mean"], summary(Y4$PmaxTotal)["Max."]-summary(Y4$PmaxTotal)["Mean"], summary(Y5$PmaxTotal)["Max."]-summary(Y5$PmaxTotal)["Mean"], summary(Y6$PmaxTotal)["Max."]-summary(Y6$PmaxTotal)["Mean"], summary(Y7$PmaxTotal)["Max."]-summary(Y7$PmaxTotal)["Mean"], summary(Y8$PmaxTotal)["Max."]-summary(Y8$PmaxTotal)["Mean"], summary(Y9$PmaxTotal)["Max."]-summary(Y9$PmaxTotal)["Mean"])
Max$Year <- c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9")  
Max$Security <- c("Max (adjusted med efflux)")

CMax <- data.frame("Population"=c(summary(CY1$PmaxTotal)["Mean"], summary(CY2$PmaxTotal)["Mean"], summary(CY3$PmaxTotal)["Mean"], summary(CY4$PmaxTotal)["Mean"], summary(CY5$PmaxTotal)["Mean"], summary(CY6$PmaxTotal)["Mean"], summary(CY7$PmaxTotal)["Mean"], summary(CY8$PmaxTotal)["Mean"], summary(CY9$PmaxTotal)["Mean"]))
CMax$floor <- c(summary(CY1$PmaxTotal)["Mean"]-summary(CY1$PmaxTotal)["Min."], summary(CY2$PmaxTotal)["Mean"]-summary(CY2$PmaxTotal)["Min."], summary(CY3$PmaxTotal)["Mean"]-summary(CY3$PmaxTotal)["Min."], summary(CY4$PmaxTotal)["Mean"]-summary(CY4$PmaxTotal)["Min."], summary(CY5$PmaxTotal)["Mean"]-summary(CY5$PmaxTotal)["Min."], summary(CY6$PmaxTotal)["Mean"]-summary(CY6$PmaxTotal)["Min."], summary(CY7$PmaxTotal)["Mean"]-summary(CY7$PmaxTotal)["Min."], summary(CY8$PmaxTotal)["Mean"]-summary(CY8$PmaxTotal)["Min."], summary(CY9$PmaxTotal)["Mean"]-summary(CY9$PmaxTotal)["Min."])
CMax$ceiling <- c(summary(CY1$PmaxTotal)["Max."]-summary(CY1$PmaxTotal)["Mean"], summary(CY2$PmaxTotal)["Max."]-summary(CY2$PmaxTotal)["Mean"], summary(CY3$PmaxTotal)["Max."]-summary(CY3$PmaxTotal)["Mean"], summary(CY4$PmaxTotal)["Max."]-summary(CY4$PmaxTotal)["Mean"], summary(CY5$PmaxTotal)["Max."]-summary(CY5$PmaxTotal)["Mean"], summary(CY6$PmaxTotal)["Max."]-summary(CY6$PmaxTotal)["Mean"], summary(CY7$PmaxTotal)["Max."]-summary(CY7$PmaxTotal)["Mean"], summary(CY8$PmaxTotal)["Max."]-summary(CY8$PmaxTotal)["Mean"], summary(CY9$PmaxTotal)["Max."]-summary(CY9$PmaxTotal)["Mean"])
CMax$Year <- c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9")  
CMax$Security <- c("Max (pure prison transition rates)")


pop <- rbind(DCSMin, Min)
pop <- rbind(pop, CMin)
pop <- rbind(pop, DCSMed)
pop <- rbind(pop, CMed)
pop <- rbind(pop, Med)
pop <- rbind(pop, DCSMax)
pop <- rbind(pop, CMax)
pop <- rbind(pop, Max)

ggplot(data=pop, aes(x=Year, y=Population, group=Security, col=Security, linetype=Security)) + geom_errorbar(aes(ymin=Population-floor, ymax=Population+ceiling), colour="black", width=.1) + geom_line(aes(linetype=Security, color=Security)) + geom_point(aes(color=Security)) + scale_linetype_manual(values=c("solid","solid","solid","twodash","dotted", "twodash","dotted","twodash", "dotted")) + scale_color_manual(values=c("#660000", "#006633", "#003366", "#993333", "#C0392B", "#009933", "#239B56", "#3399CC", "#2874A6")) + ylim(0, 5000) + ggtitle("NSW prison population") + theme(legend.position="right")
 
#Incidence and Prevalence
#Create total IDU columns
Y1$TotalIDU <- Y1$Pmin0+Y1$Pmin1+Y1$Pmin2+Y1$Pmin3+Y1$Pmin4+Y1$Pmin5+Y1$Pmed0+Y1$Pmed1+Y1$Pmed2+Y1$Pmed3+Y1$Pmed4+Y1$Pmed5+Y1$Pmax0+Y1$Pmax1+Y1$Pmax2+Y1$Pmax3+Y1$Pmax4+Y1$Pmax5
Y2$TotalIDU <- Y2$Pmin0+Y2$Pmin1+Y2$Pmin2+Y2$Pmin3+Y2$Pmin4+Y2$Pmin5+Y2$Pmed0+Y2$Pmed1+Y2$Pmed2+Y2$Pmed3+Y2$Pmed4+Y2$Pmed5+Y2$Pmax0+Y2$Pmax1+Y2$Pmax2+Y2$Pmax3+Y2$Pmax4+Y2$Pmax5
Y3$TotalIDU <- Y3$Pmin0+Y3$Pmin1+Y3$Pmin2+Y3$Pmin3+Y3$Pmin4+Y3$Pmin5+Y3$Pmed0+Y3$Pmed1+Y3$Pmed2+Y3$Pmed3+Y3$Pmed4+Y3$Pmed5+Y3$Pmax0+Y3$Pmax1+Y3$Pmax2+Y3$Pmax3+Y3$Pmax4+Y3$Pmax5
Y4$TotalIDU <- Y4$Pmin0+Y4$Pmin1+Y4$Pmin2+Y4$Pmin3+Y4$Pmin4+Y4$Pmin5+Y4$Pmed0+Y4$Pmed1+Y4$Pmed2+Y4$Pmed3+Y4$Pmed4+Y4$Pmed5+Y4$Pmax0+Y4$Pmax1+Y4$Pmax2+Y4$Pmax3+Y4$Pmax4+Y4$Pmax5
Y5$TotalIDU <- Y5$Pmin0+Y5$Pmin1+Y5$Pmin2+Y5$Pmin3+Y5$Pmin4+Y5$Pmin5+Y5$Pmed0+Y5$Pmed1+Y5$Pmed2+Y5$Pmed3+Y5$Pmed4+Y5$Pmed5+Y5$Pmax0+Y5$Pmax1+Y5$Pmax2+Y5$Pmax3+Y5$Pmax4+Y5$Pmax5
Y6$TotalIDU <- Y6$Pmin0+Y6$Pmin1+Y6$Pmin2+Y6$Pmin3+Y6$Pmin4+Y6$Pmin5+Y6$Pmed0+Y6$Pmed1+Y6$Pmed2+Y6$Pmed3+Y6$Pmed4+Y6$Pmed5+Y6$Pmax0+Y6$Pmax1+Y6$Pmax2+Y6$Pmax3+Y6$Pmax4+Y6$Pmax5
Y7$TotalIDU <- Y7$Pmin0+Y7$Pmin1+Y7$Pmin2+Y7$Pmin3+Y7$Pmin4+Y7$Pmin5+Y7$Pmed0+Y7$Pmed1+Y7$Pmed2+Y7$Pmed3+Y7$Pmed4+Y7$Pmed5+Y7$Pmax0+Y7$Pmax1+Y7$Pmax2+Y7$Pmax3+Y7$Pmax4+Y7$Pmax5
Y8$TotalIDU <- Y8$Pmin0+Y8$Pmin1+Y8$Pmin2+Y8$Pmin3+Y8$Pmin4+Y8$Pmin5+Y8$Pmed0+Y8$Pmed1+Y8$Pmed2+Y8$Pmed3+Y8$Pmed4+Y8$Pmed5+Y8$Pmax0+Y8$Pmax1+Y8$Pmax2+Y8$Pmax3+Y8$Pmax4+Y8$Pmax5
Y9$TotalIDU <- Y9$Pmin0+Y9$Pmin1+Y9$Pmin2+Y9$Pmin3+Y9$Pmin4+Y9$Pmin5+Y9$Pmed0+Y9$Pmed1+Y9$Pmed2+Y9$Pmed3+Y9$Pmed4+Y9$Pmed5+Y9$Pmax0+Y9$Pmax1+Y9$Pmax2+Y9$Pmax3+Y9$Pmax4+Y9$Pmax5
Y10$TotalIDU <- Y10$Pmin0+Y10$Pmin1+Y10$Pmin2+Y10$Pmin3+Y10$Pmin4+Y10$Pmin5+Y10$Pmed0+Y10$Pmed1+Y10$Pmed2+Y10$Pmed3+Y10$Pmed4+Y10$Pmed5+Y10$Pmax0+Y10$Pmax1+Y10$Pmax2+Y10$Pmax3+Y10$Pmax4+Y10$Pmax5

#Incidence
summary(Y1$NewCases) #Y1
summary(Y2$NewCases) #Y2
summary(Y3$NewCases) #Y3
summary(Y4$NewCases) #Y4
summary(Y5$NewCases) #Y5
summary(Y6$NewCases) #Y6
summary(Y7$NewCases) #Y7
summary(Y8$NewCases) #Y8
summary(Y9$NewCases) #Y9
summary(Y10$NewCases) #Y10

#Prevalence
summary(Y1$PrisonHCV) #Y1
summary(Y2$PrisonHCV) #Y1
summary(Y3$PrisonHCV) #Y1
summary(Y4$PrisonHCV) #Y1
summary(Y5$PrisonHCV) #Y1
summary(Y6$PrisonHCV) #Y1
summary(Y7$PrisonHCV) #Y1
summary(Y8$PrisonHCV) #Y1
summary(Y9$PrisonHCV) #Y1
summary(Y10$PrisonHCV) #Y1

#Y1
mean(Y1$NewCases)
mean(Y1$TotalPrisonPop)
mean(Y1$TotalIDU)
(mean(Y1$NewCases)/mean(Y1$TotalIDU))*100
(mean(Y1$PrisonHCV)/mean(Y1$TotalIDU))*100
#totalPop
(mean(Y1$NewCases)/mean(Y1$TotalPrisonPop))*100
(mean(Y1$PrisonHCV)/mean(Y1$TotalPrisonPop))*100

#TotalIDU incidence and prevalence
HCV <- data.frame("Incidence"=c(((mean(Y1$NewCases)/mean(Y1$TotalIDU))*100), ((mean(Y2$NewCases)/mean(Y2$TotalIDU))*100), ((mean(Y3$NewCases)/mean(Y3$TotalIDU))*100), ((mean(Y4$NewCases)/mean(Y4$TotalIDU))*100), ((mean(Y5$NewCases)/mean(Y5$TotalIDU))*100), ((mean(Y6$NewCases)/mean(Y6$TotalIDU))*100), ((mean(Y7$NewCases)/mean(Y7$TotalIDU))*100), ((mean(Y8$NewCases)/mean(Y8$TotalIDU))*100), ((mean(Y9$NewCases)/mean(Y9$TotalIDU))*100)))
HCV$IncidenceFloor <- c(((mean(Y1$NewCases)/mean(Y1$TotalIDU))*100)-((min(Y1$NewCases)/min(Y1$TotalIDU))*100), ((mean(Y2$NewCases)/mean(Y2$TotalIDU))*100)-((min(Y2$NewCases)/min(Y2$TotalIDU))*100), ((mean(Y3$NewCases)/mean(Y3$TotalIDU))*100)-((min(Y3$NewCases)/min(Y3$TotalIDU))*100), ((mean(Y4$NewCases)/mean(Y4$TotalIDU))*100)-((min(Y4$NewCases)/min(Y4$TotalIDU))*100), ((mean(Y5$NewCases)/mean(Y5$TotalIDU))*100)-((min(Y5$NewCases)/min(Y5$TotalIDU))*100), ((mean(Y6$NewCases)/mean(Y6$TotalIDU))*100)-((min(Y6$NewCases)/min(Y6$TotalIDU))*100), ((mean(Y7$NewCases)/mean(Y7$TotalIDU))*100)-((min(Y7$NewCases)/min(Y7$TotalIDU))*100), ((mean(Y8$NewCases)/mean(Y8$TotalIDU))*100)-((min(Y8$NewCases)/min(Y8$TotalIDU))*100), ((mean(Y9$NewCases)/mean(Y9$TotalIDU))*100)-((min(Y9$NewCases)/min(Y9$TotalIDU))*100))  
HCV$IncidenceCeiling <- c(((max(Y1$NewCases)/max(Y1$TotalIDU))*100)-((mean(Y1$NewCases)/mean(Y1$TotalIDU))*100), ((max(Y2$NewCases)/max(Y2$TotalIDU))*100)-((mean(Y2$NewCases)/mean(Y2$TotalIDU))*100), ((max(Y3$NewCases)/max(Y3$TotalIDU))*100)-((mean(Y3$NewCases)/mean(Y3$TotalIDU))*100), ((max(Y4$NewCases)/max(Y4$TotalIDU))*100)-((mean(Y4$NewCases)/mean(Y4$TotalIDU))*100), ((max(Y5$NewCases)/max(Y5$TotalIDU))*100)-((mean(Y5$NewCases)/mean(Y5$TotalIDU))*100), ((max(Y6$NewCases)/max(Y6$TotalIDU))*100)-((mean(Y6$NewCases)/mean(Y6$TotalIDU))*100), ((max(Y7$NewCases)/max(Y7$TotalIDU))*100)-((mean(Y7$NewCases)/mean(Y7$TotalIDU))*100), ((max(Y8$NewCases)/max(Y8$TotalIDU))*100)-((mean(Y8$NewCases)/mean(Y8$TotalIDU))*100), ((max(Y9$NewCases)/max(Y9$TotalIDU))*100)-((mean(Y9$NewCases)/mean(Y9$TotalIDU))*100))
HCV$Prevalence <- c(((mean(Y1$PrisonHCV)/mean(Y1$TotalIDU))*100), ((mean(Y2$PrisonHCV)/mean(Y2$TotalIDU))*100), ((mean(Y3$PrisonHCV)/mean(Y3$TotalIDU))*100), ((mean(Y4$PrisonHCV)/mean(Y4$TotalIDU))*100), ((mean(Y5$PrisonHCV)/mean(Y5$TotalIDU))*100), ((mean(Y6$PrisonHCV)/mean(Y6$TotalIDU))*100), ((mean(Y7$PrisonHCV)/mean(Y7$TotalIDU))*100), ((mean(Y8$PrisonHCV)/mean(Y8$TotalIDU))*100), ((mean(Y9$PrisonHCV)/mean(Y9$TotalIDU))*100))
HCV$PrevalenceFloor <- c(((mean(Y1$PrisonHCV)/mean(Y1$TotalIDU))*100)-((min(Y1$PrisonHCV)/min(Y1$TotalIDU))*100), ((mean(Y2$PrisonHCV)/mean(Y2$TotalIDU))*100)-((min(Y2$PrisonHCV)/min(Y2$TotalIDU))*100), ((mean(Y3$PrisonHCV)/mean(Y3$TotalIDU))*100)-((min(Y3$PrisonHCV)/min(Y3$TotalIDU))*100), ((mean(Y4$PrisonHCV)/mean(Y4$TotalIDU))*100)-((min(Y4$PrisonHCV)/min(Y4$TotalIDU))*100), ((mean(Y5$PrisonHCV)/mean(Y5$TotalIDU))*100)-((min(Y5$PrisonHCV)/min(Y5$TotalIDU))*100), ((mean(Y6$PrisonHCV)/mean(Y6$TotalIDU))*100)-((min(Y6$PrisonHCV)/min(Y6$TotalIDU))*100), ((mean(Y7$PrisonHCV)/mean(Y7$TotalIDU))*100)-((min(Y7$PrisonHCV)/min(Y7$TotalIDU))*100), ((mean(Y8$PrisonHCV)/mean(Y8$TotalIDU))*100)-((min(Y8$PrisonHCV)/min(Y8$TotalIDU))*100), ((mean(Y9$PrisonHCV)/mean(Y9$TotalIDU))*100)-((min(Y9$PrisonHCV)/min(Y9$TotalIDU))*100))
HCV$PrevalenceCeiling <- c(((mean(Y1$PrisonHCV)/mean(Y1$TotalIDU))*100)-((min(Y1$PrisonHCV)/min(Y1$TotalIDU))*100), ((mean(Y2$PrisonHCV)/mean(Y2$TotalIDU))*100)-((min(Y2$PrisonHCV)/min(Y2$TotalIDU))*100), ((mean(Y3$PrisonHCV)/mean(Y3$TotalIDU))*100)-((min(Y3$PrisonHCV)/min(Y3$TotalIDU))*100), ((mean(Y4$PrisonHCV)/mean(Y4$TotalIDU))*100)-((min(Y4$PrisonHCV)/min(Y4$TotalIDU))*100), ((mean(Y5$PrisonHCV)/mean(Y5$TotalIDU))*100)-((min(Y5$PrisonHCV)/min(Y5$TotalIDU))*100), ((mean(Y6$PrisonHCV)/mean(Y6$TotalIDU))*100)-((min(Y6$PrisonHCV)/min(Y6$TotalIDU))*100), ((mean(Y7$PrisonHCV)/mean(Y7$TotalIDU))*100)-((min(Y7$PrisonHCV)/min(Y7$TotalIDU))*100), ((mean(Y8$PrisonHCV)/mean(Y8$TotalIDU))*100)-((min(Y8$PrisonHCV)/min(Y8$TotalIDU))*100), ((mean(Y9$PrisonHCV)/mean(Y9$TotalIDU))*100)-((min(Y9$PrisonHCV)/min(Y9$TotalIDU))*100))
HCV$Year <- c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9")

#Incidence
ggplot(data=HCV, aes(x=Year, y=Incidence, group=1)) + geom_errorbar(aes(ymin=Incidence-IncidenceFloor, ymax=Incidence+IncidenceFloor), colour="black", width=.1) + geom_line() + geom_point() + ylim(0, 0.5) + ggtitle("HCV incidence in NSW prison") + theme(legend.position="right")

#Prevalence 
ggplot(data=HCV, aes(x=Year, y=Prevalence, group=1)) + geom_errorbar(aes(ymin=Prevalence-PrevalenceFloor, ymax=Prevalence+PrevalenceFloor), colour="black", width=.1) +  geom_line() + geom_point() + ylim(0, 10) + ggtitle("HCV Prevalence in NSW prison") + theme(legend.position="right")



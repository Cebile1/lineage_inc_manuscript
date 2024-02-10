library(dplyr)
require(epiR)

#input isolate data from Supplementary data

paper2 <- read.csv("../Downloads/merged_final_ZA.csv")

#Pen 
## Need to exclude vaccine periods so that only pre-PCV and late-PCV13 are compared. Additionally remove rows where there are NF values in Predicted_PEN_SIR
dat <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13" & Predicted_PEN_SIR!="NF"))

#dat$Vaccine_Period <- factor(dat$Vaccine_Period,levels(dat$Vaccine_Period)[c(2,1)])
## Need to select the correct levels in the correct chronilogical order ("PCV7","late-PCV13") and add in factor() (not sure why factors wasn't required previously, perhaps a change in versions/package updates)
dat$Vaccine_Period <- factor(dat$Vaccine_Period,levels(factor(dat$Vaccine_Period))[c(2,1)])
## Need to add in factor()
dat$Predicted_PEN_SIR <- factor(dat$Predicted_PEN_SIR,levels(factor(dat$Predicted_PEN_SIR))[c(2,1)])
glm_penNS <- glm(Predicted_PEN_SIR ~ Vaccine_Period, data = dat, family = binomial)
##Need to remove rows without SIR predictions ("NF")     
datNVT <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13" & Vaccine_Status=="NVT" & Predicted_PEN_SIR!="NF"))
#datNVT$Vaccine_Period <- factor(datNVT$Vaccine_Period,levels(datNVT$Vaccine_Period)[c(2,1)])
## Need to put the levels in the correct order and add in factor()
datNVT$Vaccine_Period <- factor(datNVT$Vaccine_Period,levels(factor(datNVT$Vaccine_Period))[c(2,1)])
#datNVT$Predicted_PEN_SIR <- factor(datNVT$Predicted_PEN_SIR,levels(datNVT$Predicted_PEN_SIR)[c(2,1)])
## Need to add in factor()
datNVT$Predicted_PEN_SIR <- factor(datNVT$Predicted_PEN_SIR,levels(factor(datNVT$Predicted_PEN_SIR))[c(2,1)])
glm_penNS_NVT <- glm(Predicted_PEN_SIR ~ Vaccine_Period, data = datNVT, family = binomial)

#repeat removal of rows with missing data for each antimicrobial and add in factor() throughout

#Ery
#dat <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13"))
dat <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13" & Predicted_ERY_SIR!="Flag"))
#dat$Vaccine_Period <- factor(dat$Vaccine_Period,levels(dat$Vaccine_Period)[c(2,1)])
dat$Vaccine_Period <- factor(dat$Vaccine_Period,levels(factor(dat$Vaccine_Period))[c(2,1)])
#dat$Predicted_ERY_SIR <- factor(dat$Predicted_ERY_SIR,levels(dat$Predicted_ERY_SIR)[c(2,1)])
dat$Predicted_ERY_SIR <- factor(dat$Predicted_ERY_SIR,levels(factor(dat$Predicted_ERY_SIR))[c(2,1)])
glm_EryNS <- glm(Predicted_ERY_SIR ~ Vaccine_Period, data = dat, family = binomial)

#datNVT <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13" & Vaccine_Status=="NVT")
datNVT <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13" & Vaccine_Status=="NVT"& Predicted_ERY_SIR!="Flag"))
#datNVT$Vaccine_Period <- factor(datNVT$Vaccine_Period,levels(datNVT$Vaccine_Period)[c(2,1)])
datNVT$Vaccine_Period <- factor(datNVT$Vaccine_Period,levels(factor(datNVT$Vaccine_Period))[c(2,1)])
#datNVT$Predicted_ERY_SIR <- factor(datNVT$Predicted_ERY_SIR,levels(datNVT$Predicted_ERY_SIR)[c(2,1)])
datNVT$Predicted_ERY_SIR <- factor(datNVT$Predicted_ERY_SIR,levels(factor(datNVT$Predicted_ERY_SIR))[c(2,1)])
glm_EryNS_NVT <- glm(Predicted_ERY_SIR ~ Vaccine_Period, data = datNVT, family = binomial)


#Tet
dat <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13"))
dat$Vaccine_Period <- factor(dat$Vaccine_Period,levels(factor(dat$Vaccine_Period))[c(2,1)])
dat$Predicted_TET_SIR <- factor(dat$Predicted_TET_SIR,levels(factor(dat$Predicted_TET_SIR))[c(2,1)])
glm_tetNS <- glm(Predicted_TET_SIR ~ Vaccine_Period, data = dat, family = binomial)


datNVT <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13" & Vaccine_Status=="NVT"))
datNVT$Vaccine_Period <- factor(datNVT$Vaccine_Period,levels(factor(datNVT$Vaccine_Period))[c(2,1)])
datNVT$Predicted_TET_SIR <- factor(datNVT$Predicted_TET_SIR,levels(factor(datNVT$Predicted_TET_SIR))[c(2,1)])
#glm_tetNS_NVT <- glm(Predicted_TET_SIR ~ Vaccine_Period + Country, data = datNVT, family = binomial)
## Country not required, factor() added in 
glm_tetNS_NVT <- glm(Predicted_TET_SIR ~ Vaccine_Period, data = datNVT, family = binomial)

#COT
#dat <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13"))
##Bionomial model expecting two levels e.g. S and R, here with intermediate values we need S and NS (R+I) so I created that column in the input file. I also filter out rows with missing data
dat <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13"  & Predicted_COT_SIR!="Flag"))
dat$Vaccine_Period <- factor(dat$Vaccine_Period,levels(factor(dat$Vaccine_Period))[c(2,1)])
dat$Predicted_COT_SIR <- factor(dat$Predicted_COT_SIR,levels(factor(dat$Predicted_COT_SIR))[c(2,1)])
glm_COT_NS <- glm(Predicted_COT_SIR ~ Vaccine_Period, data = dat, family = binomial)
#Country not required, added in factor, removed flag rows
datNVT <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13" & Vaccine_Status=="NVT" & Predicted_COT_SIR!="Flag"))
datNVT$Vaccine_Period <- factor(datNVT$Vaccine_Period,levels(factor(datNVT$Vaccine_Period))[c(2,1)])
datNVT$Predicted_COT_SIR <- factor(datNVT$Predicted_COT_SIR,levels(factor(datNVT$Predicted_COT_SIR))[c(2,1)])
glm_COT_NS_NVT <- glm(Predicted_COT_SIR ~ Vaccine_Period, data = datNVT, family = binomial)

#Chl
dat <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13" & Predicted_CHL_SIR!="Flag"))
dat$Vaccine_Period <- factor(dat$Vaccine_Period,levels(factor(dat$Vaccine_Period))[c(2,1)])
dat$Predicted_CHL_SIR <- factor(dat$Predicted_CHL_SIR,levels(factor(dat$Predicted_CHL_SIR))[c(2,1)])
glm_ChlNS <- glm(Predicted_CHL_SIR ~ Vaccine_Period, data = dat, family = binomial)

datNVT <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13" & Vaccine_Status=="NVT" & Predicted_CHL_SIR!="Flag"))
datNVT$Vaccine_Period <- factor(datNVT$Vaccine_Period,levels(factor(datNVT$Vaccine_Period))[c(2,1)])
datNVT$Predicted_CHL_SIR <- factor(datNVT$Predicted_CHL_SIR,levels(factor(datNVT$Predicted_CHL_SIR))[c(2,1)])
glm_ChlNS_NVT <- glm(Predicted_CHL_SIR ~ Vaccine_Period, data = datNVT, family = binomial)

#Cli
dat <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13" & Predicted_CLI_SIR!="Flag"))
dat$Vaccine_Period <- factor(dat$Vaccine_Period,levels(factor(dat$Vaccine_Period))[c(2,1)])
dat$Predicted_CLI_SIR <- factor(dat$Predicted_CLI_SIR,levels(factor(dat$Predicted_CLI_SIR))[c(2,1)])
glm_CliNS <- glm(Predicted_CLI_SIR ~ Vaccine_Period, data = dat, family = binomial)

datNVT <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13" & Vaccine_Status=="NVT" & Predicted_CLI_SIR!="Flag"))
datNVT$Vaccine_Period <- factor(datNVT$Vaccine_Period,levels(factor(datNVT$Vaccine_Period))[c(2,1)])
datNVT$Predicted_CLI_SIR <- factor(datNVT$Predicted_CLI_SIR,levels(factor(datNVT$Predicted_CLI_SIR))[c(2,1)])
glm_CliNS_NVT <- glm(Predicted_CLI_SIR ~ Vaccine_Period, data = datNVT, family = binomial)


#Number of classes
dat <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13"))
dat$Vaccine_Period <- factor(dat$Vaccine_Period,levels(factor(dat$Vaccine_Period))[c(2,1)])
dat$MDR <- factor(dat$MDR,levels(factor(dat$MDR))[c(1,2)])
glm_MDR <- glm(MDR ~ Vaccine_Period, data = dat, family = binomial)

datNVT <- droplevels(subset(paper2, Vaccine_Period!="PCV7" & Vaccine_Period!="early-PCV13" & Vaccine_Status=="NVT"))
datNVT$Vaccine_Period <- factor(datNVT$Vaccine_Period,levels(factor(datNVT$Vaccine_Period))[c(2,1)])
datNVT$MDR <- factor(datNVT$MDR,levels(factor(datNVT$MDR))[c(1,2)])
glm_MDR_NVT <- glm(MDR ~ Vaccine_Period, data = datNVT, family = binomial)

results <-rbind(
  c("PenNS_post",summary(glm_penNS)$coefficients['Vaccine_Periodlate-PCV13',]),
  c("ChlNS_post",summary(glm_ChlNS)$coefficients['Vaccine_Periodlate-PCV13',]),
  c("EryNS_post",summary(glm_EryNS)$coefficients['Vaccine_Periodlate-PCV13',]),
  c("CotNS_post",summary(glm_COT_NS)$coefficients['Vaccine_Periodlate-PCV13',]),
  c("TetNS_post",summary(glm_tetNS)$coefficients['Vaccine_Periodlate-PCV13',]),
  c("CliNS_post",summary(glm_CliNS)$coefficients['Vaccine_Periodlate-PCV13',]),
  c("MDR_post",summary(glm_MDR)$coefficients['Vaccine_Periodlate-PCV13',]),
  c("PenNS_NVTpost",summary(glm_penNS_NVT)$coefficients['Vaccine_Periodlate-PCV13',]),
  c("ChlNS_NVTpost",summary(glm_ChlNS_NVT)$coefficients['Vaccine_Periodlate-PCV13',]),
  c("EryNS_NVTpost",summary(glm_EryNS_NVT)$coefficients['Vaccine_Periodlate-PCV13',]),
  c("CotNS_NVTpost",summary(glm_COT_NS_NVT)$coefficients['Vaccine_Periodlate-PCV13',]),
  c("TetNS_NVTpost",summary(glm_tetNS_NVT)$coefficients['Vaccine_Periodlate-PCV13',]),
  c("CliNS_NVTpost",summary(glm_CliNS_NVT)$coefficients['Vaccine_Periodlate-PCV13',]),
  c("MDR_NVTpost",summary(glm_MDR_NVT)$coefficients['Vaccine_Periodlate-PCV13',])
)

colnames(results)[1] <- "Type"
write.csv(results, file ="Tab4_post_NVT_resistance.csv", row.names = FALSE)

#Fig2_glmIRR_adjcounts_perGPSC_model_select_divide_by_selection.R

require(tidyverse)
library(dplyr)
library(MASS)
library(pscl)
library(sandwich)
library(epiR)
library(readxl)
library(tidyr)

####Impute GPSCs based on proportion detected in sample

#Import datasets
#Main GERMS DB

germs <- read_excel("GERMS_MOTHER_54,199.xlsx", 
                    sheet = "<18y_all_cultures_17,121")

#Genotyping
geno <- read_excel("GERMS_GENO_3104.xlsx", 
                   sheet = "Sheet1")


#Population estimates
pop <- read_excel("mid Year pop.xlsx") %>%
filter(Age_category=="all") %>% 
  mutate(period = case_when(Year<=2008 ~ "Pre-PCV",
                            Year>=2009 & Year<=2010 ~ "PCV7",
                            Year>=2011 & Year<=2014 ~ "early-PCV13",
                            Year>=2015 ~ "late-PCV13")) %>% 
  rename("Year"="Year", "Population"="pop") %>% 
  dplyr::select(Year, period, Population)

#Summary of case numbers by age and Year
caseno_sum <- germs %>% 
  group_by(Age_category, Year) %>% 
  summarise(totalcases = n())

#Summary of case numbers per lineage, age group and Year
geno_no <- geno %>% 
  filter(!is.na(GPSC)) %>% 
  mutate(GPSC = factor(GPSC)) %>% 
  group_by(Age_category, Year) %>% 
  mutate(total = n()) %>% 
  ungroup() %>% 
  group_by(Age_category, Year, GPSC) %>% 
  summarise(caseno = n(),
            total = max(total)) %>% 
  complete(GPSC, fill = list(caseno = 0)) %>% 
  mutate(prop_geno = caseno/total) %>% 
  replace(is.na(.), 0)

#Merge total GERMS case numbers with lineage proportions
imp_no <- left_join(geno_no, caseno_sum, by=c("Age_category"="Age_category", "Year"="Year")) %>% 
  mutate(imp_cases = totalcases*prop_geno,
         period = case_when(Year<=2008 ~ "Pre-PCV",
                            Year>=2009 & Year<=2010 ~ "PCV7",
                            Year>=2011 & Year<=2014 ~ "early-PCV13",
                            Year>=2015 ~ "late-PCV13")) %>% 
  ungroup() %>% 
  group_by(Year, period, GPSC) %>% 
  summarise(Freq = sum(caseno),
            estimated.cases = round(sum(imp_cases))) %>% 
  rename("Year"="Year")

#period <- c("Post-PCV7","Post-PCV13")
period <- c("PCV7","early-PCV13","late-PCV13")

#######Calculate IRR per GPSC for different PCV periods############
GPSC_IRR_PCV7 <- matrix(data=NA,nrow=0,ncol=19)
GPSC_IRR_early_PCV13 <- matrix(data=NA,nrow=0,ncol=19)
GPSC_IRR_late_PCV13 <- matrix(data=NA,nrow=0,ncol=19)

for (post in period){

  one_pop <- subset(pop, period==post | period=="Pre-PCV") 
  
  #Get list of GPSCs
  GPSCs <- unique(imp_no$GPSC)
  #Calculate IRR pre-PCV vs post-PCV13 for each country
  for (cluster in GPSCs){
    #GPSC_dat <- droplevels(subset(dat, GPSC==cluster, select = c(Year)))
    tab <- subset(imp_no, GPSC==cluster & (period==post | period=="Pre-PCV"), select = c(Year, Freq, estimated.cases)) %>% 
      rename("GPSC_dat"="Year")
    #combine with genome counts per Year in NVT-GPSCs/VT-GPSCs
    pop_tab <- merge(one_pop, tab, by.y = "GPSC_dat", by.x = "Year", all.x=TRUE)
    #pop_tab <- merge(one_pop, tab, by.x = "Year", all.x=TRUE)
    pop_tab$period <- factor(pop_tab$period, levels = c("Pre-PCV", post))
    pop_tab$Freq[is.na(pop_tab$Freq)] <- 0
    pre_counts <- sum(subset(pop_tab, period=="Pre-PCV", select=c(Freq)))
    post_counts <- sum(subset(pop_tab, period==post, select=c(Freq)))
    pre_cases <- sum(subset(pop_tab, period=="Pre-PCV", select=c(estimated.cases)))
    post_cases <- sum(subset(pop_tab, period==post, select=c(estimated.cases)))
    pre_Years <- dim(subset(pop_tab, period=="Pre-PCV"))[1]
    post_Years <- dim(subset(pop_tab, period==post))[1]      
    pre_population_avg <- sum(subset(pop_tab, period=="Pre-PCV", select=c(Population)))/pre_Years
    post_population_avg <- sum(subset(pop_tab, period==post, select=c(Population)))/post_Years
    #round estimated counts
    pop_tab$estimated.cases <- round(pop_tab$estimated.cases)
    
    if  (sum(pop_tab$Freq)>=5){ #JK note: Why is this limiting to if there are >=5 cases in the full period?
      #Add one estimated case to pre or post if no genomes
      if (sum(subset(pop_tab, period=="Pre-PCV")['Freq'])==0 | sum(subset(pop_tab, period==post)['Freq'])==0){
        pop_tab$estimated.cases <- pop_tab$estimated.cases+1
        add <- 1
      } else {
        add <- 0
      }
      
      #calculate IRR using period averages
      IRR_by2 <- matrix(c(post_cases/post_Years,pre_cases/pre_Years,post_population_avg,pre_population_avg), nrow = 2, byrow = TRUE)
      rownames(IRR_by2) <- c("post_avg_population", "pre_avg_population"); colnames(IRR_by2) <- c("post_est_annual_cases", "pre_est_annual_cases")
      IRR_by2 <- round(IRR_by2)
      #if 0 in one period add 1 to both
      if (sum(IRR_by2[1,]==0,na.rm=T) > 0){
        IRR_by2[1,] <- IRR_by2[1,]+1
      }
      res <- epi.2by2(IRR_by2, method = "cross.sectional", conf.level = 0.95, units = 100,
                      outcome = "as.rows")  
      #IRR_calc <- res$res$IRR.strata.wald$est
      IRR_calc <- res$massoc.detail$OR.strata.wald$est
      #confi_lo_calc <- res$res$IRR.strata.wald$lower
      confi_lo_calc <-res$massoc.detail$OR.strata.wald$lower
      #confi_up_calc <- res$res$IRR.strata.wald$upper
      confi_up_calc <- res$massoc.detail$OR.strata.wald$upper
      #ps_calc <- res$res$chisq.strata$p.value
      ps_calc <- res$massoc.detail$chi2.strata.fisher$p.value.2s
      
      #runglm
      res = glm(estimated.cases ~ period + offset(log(Population)) , data=pop_tab, family = "poisson")
      #test fit
      GoFit <- 1 - pchisq(summary(res)$deviance, 
                          summary(res)$df.residual)
      #set zero inflation values to NA for models with no zeros
      ZI <- NA
      ZI_period <- NA
      
      #assess zero inflation total and period
      if (sum(pop_tab$estimated.cases ==0)>0){
        res.zip = zeroinfl(estimated.cases ~ period|period, data = pop_tab)
        ZI <- summary(res.zip)$coefficients$zero[1,4]
        ZI_period <- summary(res.zip)$coefficients$zero[2,4]
      }
      
      #If poisson does fit
      if (GoFit >0.05){
        model <- "poisson"
        converged <- res$converged
        pre_post <- res$coefficients %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
        IRR <- unname(pre_post[2])
        ps <-  unname(coef(summary(res))[,4][2])
        pre_postconfint <- try(confint(res) %>% exp)
        confi_lo <- pre_postconfint[2,1]
        confi_up <- pre_postconfint[2,2]
        #average incidence before vaccine (intercept)/100,000 population
        pre_inc <- pre_post[1]*100000
        #average incidence post vaccine/100,000 population derived from pre incidence and IRR
        post_inc <- (pre_post[1]*100000)*pre_post[2]
        
        #If poisson does NOT quite fit (minor violation)
      } else if (GoFit >=0.01 & GoFit <=0.05){
        #Calcultate robust standard errors
        model <- "poisson robust SE"
        converged <- res$converged
        pre_post <- res$coefficients %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
        IRR <- unname(pre_post[2])
        cov.res <- vcovHC(res, type="HC0")
        std.err <- sqrt(diag(cov.res))
        r.est <- cbind(Estimate= coef(res), "Robust SE" = std.err,
                       "Pr(>|z|)" = 2 * pnorm(abs(coef(res)/std.err), lower.tail=FALSE),
                       LL = coef(res) - 1.96 * std.err,
                       UL = coef(res) + 1.96 * std.err)
        ps <- r.est[2,3]
        #convert confidence intervals from log values
        confi_lo <- r.est[2,4] %>% exp
        confi_up <- r.est[2,5] %>% exp
        #average incidence before vaccine (intercept)/100,000 population
        pre_inc <- pre_post[1]*100000
        #average incidence post vaccine/100,000 population derived from pre incidence and IRR
        post_inc <- (pre_post[1]*100000)*pre_post[2]
        
        #If poisson does NOT fit
      } else {
        #fit negative bionomial (this gives theta.ml warning if zero inflation present, which is the next model to be fitted if glm.nb is <0.05)
        model <- "negative bionomial"
        res.nb = glm.nb(estimated.cases ~ period + offset(log(Population)) , data=pop_tab)
        #test fit
        GoFit <- 1 - pchisq(summary(res.nb)$deviance,
                            summary(res.nb)$df.residual)
        converged <- res.nb$converged
        pre_post <- res.nb$coefficients %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
        IRR <- unname(pre_post[2])
        ps <-  unname(coef(summary(res.nb))[,4][2])
        pre_postconfint <- try(confint(res.nb) %>% exp)
        confi_lo <- pre_postconfint[2,1]
        confi_up <- pre_postconfint[2,2]
        #average incidence before vaccine (intercept)/100,000 population
        pre_inc <- pre_post[1]*100000
        #average incidence post vaccine/100,000 population derived from pre incidence and IRR
        post_inc <- (pre_post[1]*100000)*pre_post[2]
        #assess zero inflation total and period
        if (sum(pop_tab$estimated.cases ==0)>0){
          res.zip = zeroinfl(estimated.cases ~ period|period, data = pop_tab)
          ZI <- summary(res.zip)$coefficients$zero[1,4]
          ZI_period <- summary(res.zip)$coefficients$zero[2,4]
        }
        #If negative bionomial does NOT fit
        if (GoFit <0.05 & !is.na(ZI) & ZI <0.05){
          #fit zero inflated negative binomial
          res.zinb = zeroinfl(estimated.cases ~ period|1 + offset(log(population)) , data=pop_tab, dist="negbin")
          #test fit using Log(theta)
          GoFit <- summary(res.zinb)$coefficients$count[3,4]
          #select model based on theta
          if(summary(res.zinb)$coefficients$count[3,4] <0.05){
            model <- "zero inflated negative binomial"
            converged <- res.zinb$converged
            ZI <- summary(res.zinb)$coefficients$zero[1,4]
            pre_post <- res.zinb$coefficients$count %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
            IRR <- unname(pre_post[2])
            ps <-  unname(coef(summary(res.zinb))$count[,4][2])
            pre_postconfint <- try(confint(res.zinb) %>% exp)
            confi_lo <- pre_postconfint[2,1]
            confi_up <- pre_postconfint[2,2]
            #average incidence before vaccine (intercept)/100,000 population
            pre_inc <- pre_post[1]*100000
            #average incidence post vaccine/100,000 population derived from pre incidence and IRR
            post_inc <- (pre_post[1]*100000)*pre_post[2]
          } else {
            model <- "zero inflated poisson"
            res.zip = zeroinfl(estimated.cases ~ period|1, data = pop_tab)
            ZI <- summary(res.zip)$coefficients$zero[1,4]
            converged <- res.zip$converged
            pre_post <- res.zip$coefficients$count %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
            IRR <- unname(pre_post[2])
            ps <-  unname(coef(summary(res.zip))$count[,4][2])
            pre_postconfint <- try(confint(res.zip) %>% exp)
            confi_lo <- pre_postconfint[2,1]
            confi_up <- pre_postconfint[2,2]
            #average incidence before vaccine (intercept)/100,000 population
            pre_inc <- pre_post[1]*100000
            #average incidence post vaccine/100,000 population derived from pre incidence and IRR
            post_inc <- (pre_post[1]*100000)*pre_post[2]
            
          }   
        }
      }
      if (GoFit <0.05 & model != "poisson robust SE"){
        #calculate IRR using period averages
        IRR_by2 <- matrix(c(post_cases/post_Years,pre_cases/pre_Years,post_population_avg,pre_population_avg), nrow = 2, byrow = TRUE)
        colnames(IRR_by2) <- c("post", "pre"); rownames(IRR_by2) <- c("est_annual_cases", "population")
        IRR_by2 <- round(IRR_by2)
        #if 0 in one period add 1 to both
        if (sum(IRR_by2[1,]==0,na.rm=T) > 0){
          IRR_by2[1,] <- IRR_by2[1,]+1
        }
        #calculate IRR using period averages
        model <- "none"
        #res <- epi.2by2(IRR_by2, method = "cross.sectional", conf.level = 0.95, units = 100,outcome = "as.rows")  
        IRR <- NA
        confi_lo <- NA
        confi_up <- NA
        ps <- NA
        GoFit <- NA
        converged <- NA
      }
      if (post == "PCV7"){
        GPSC_IRR_PCV7 <- rbind(GPSC_IRR_PCV7,c(cluster, model, GoFit, ZI, ZI_period, IRR, confi_lo, confi_up, ps, pre_counts, pre_cases, post_counts, post_cases,pre_inc, post_inc,IRR_calc,confi_lo_calc,confi_up_calc,ps_calc))
      } else if (post == "early-PCV13") {
        GPSC_IRR_early_PCV13 <- rbind(GPSC_IRR_early_PCV13,c(cluster, model, GoFit, ZI, ZI_period, IRR, confi_lo, confi_up, ps, pre_counts, pre_cases, post_counts, post_cases,pre_inc, post_inc,IRR_calc,confi_lo_calc,confi_up_calc,ps_calc))
      } else {
        GPSC_IRR_late_PCV13 <- rbind(GPSC_IRR_late_PCV13,c(cluster, model, GoFit, ZI, ZI_period, IRR, confi_lo, confi_up, ps, pre_counts, pre_cases, post_counts, post_cases,pre_inc, post_inc,IRR_calc,confi_lo_calc,confi_up_calc,ps_calc))
      }
    }
  }
}


#Changed to reflect columns in pop_Years_ZA.csv and change columns to adjust p-value because of removal of country column
#adjust model IRR p-value
GPSC_IRR_PCV7 <- cbind(GPSC_IRR_PCV7,p.adjust(GPSC_IRR_PCV7[,9], method = "BH", n=length(GPSC_IRR_PCV7[,ncol(GPSC_IRR_PCV7)])))
GPSC_IRR_early_PCV13 <- cbind(GPSC_IRR_early_PCV13,p.adjust(GPSC_IRR_early_PCV13[,9], method = "BH", n=length(GPSC_IRR_early_PCV13[,ncol(GPSC_IRR_early_PCV13)])))
GPSC_IRR_late_PCV13 <- cbind(GPSC_IRR_late_PCV13,p.adjust(GPSC_IRR_late_PCV13[,9], method = "BH", n=length(GPSC_IRR_late_PCV13[,ncol(GPSC_IRR_late_PCV13)])))
#adjust average IRR p-value
GPSC_IRR_PCV7 <- cbind(GPSC_IRR_PCV7,p.adjust(GPSC_IRR_PCV7[,19], method = "BH", n=length(GPSC_IRR_PCV7[,ncol(GPSC_IRR_PCV7)])))
GPSC_IRR_early_PCV13 <- cbind(GPSC_IRR_early_PCV13,p.adjust(GPSC_IRR_early_PCV13[,19], method = "BH", n=length(GPSC_IRR_early_PCV13[,ncol(GPSC_IRR_early_PCV13)])))
GPSC_IRR_late_PCV13 <- cbind(GPSC_IRR_late_PCV13,p.adjust(GPSC_IRR_late_PCV13[,19], method = "BH", n=length(GPSC_IRR_late_PCV13[,ncol(GPSC_IRR_late_PCV13)])))
#reorder columns
GPSC_IRR_PCV7 <- cbind(GPSC_IRR_PCV7[,1:9],GPSC_IRR_PCV7[,20],GPSC_IRR_PCV7[,10:19],GPSC_IRR_PCV7[,21])
GPSC_IRR_early_PCV13 <- cbind(GPSC_IRR_early_PCV13[,1:9],GPSC_IRR_early_PCV13[,20],GPSC_IRR_early_PCV13[,10:19],GPSC_IRR_early_PCV13[,21])
GPSC_IRR_late_PCV13 <- cbind(GPSC_IRR_late_PCV13[,1:9],GPSC_IRR_late_PCV13[,20],GPSC_IRR_late_PCV13[,10:19],GPSC_IRR_late_PCV13[,21])
#write files
colnames(GPSC_IRR_PCV7) <- c("GPSC","model","GOF","ZI", "ZI_period","IRR","lower","upper","p","adj.p","pre-genomes", "pre-estimated cases","post-genomes","post-estimated cases", "pre-avg-incidence","post-avg-incidence","IRR_average","lower","upper","p-value","adj.avg.p")
colnames(GPSC_IRR_early_PCV13) <- c("GPSC","model","GOF","ZI", "ZI_period","IRR","lower","upper","p","adj.p","pre-genomes", "pre-estimated cases","post-genomes","post-estimated cases", "pre-avg-incidence","post-avg-incidence","IRR_average","lower","upper","p-value","adj.avg.p")
colnames(GPSC_IRR_late_PCV13) <- c("GPSC","model","GOF","ZI", "ZI_period","IRR","lower","upper","p","adj.p","pre-genomes", "pre-estimated cases","post-genomes","post-estimated cases", "pre-avg-incidence","post-avg-incidence","IRR_average","lower","upper","p-value","adj.avg.p")
write.csv(GPSC_IRR_PCV7, file ="FigS2-4_PCV7_glmIRR_pseudo_adjfreq_model_select_divide_by_selection.csv", row.names = FALSE)
write.csv(GPSC_IRR_early_PCV13 , file ="FigS2-4_earlyPCV13_glmIRR_pseudo_adjfreq_model_select_divide_by_selection.csv", row.names = FALSE)
write.csv(GPSC_IRR_late_PCV13 , file ="FigS2-4_latePCV13_glmIRR_pseudo_adjfreq_model_select_divide_by_selection.csv", row.names = FALSE)


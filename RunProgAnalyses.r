library(survival)
library(dplyr)
library(data.table)

# InFile = '/home/ivm/Documents/Output/ProgressionPrelim/R10/PRSexp_20220927/InterveneIn'
args = commandArgs(trailingOnly=TRUE)
InFile = args[1]
OutFile = 'ProgAnaOut'

EndPt = c('COLORECTAL', 'PROSTATE', 'BREAST', 'CHD', 'HF', 'AF', 'CKD', 'AD', 'T2D', 'STR', 'CVD_HARD' )
SurvMin = 0.25

Input = fread(InFile)
Input[['BirthYear']] = as.numeric(format(as.Date(Input$DateOfBirth), "%Y"))
Input[['AgeStartFollow']] = as.numeric(as.difftime(as.Date(Input$FollowUpStart)-as.Date(Input$DateOfBirth), units = "days"))/365.25
Input[['AgeEndFollow']] = as.numeric(as.difftime(as.Date(Input$FollowUpEnd)-as.Date(Input$DateOfBirth), units = "days"))/365.25
Input = Input[Input$AgeEndFollow >= 0, ]
Input = Input[Input$AgeStartFollow >= 0, ]
Input = Input[Input$AgeEndFollow >= Input$AgeStartFollow, ]

Res = data.frame(matrix(ncol = 12, nrow = 0))
colnames(Res) = c("Disease", "Target", "Predictor", "Model", "MaxSurvTime", "Sample", "Beta", "SE", "Pval", "nEvent", "nTotSample", "AnaFlag")

for (item in EndPt) {
  Covar = c('BirthYear','Sex','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')
  CovarList = 'BirthYear+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'
  if (item == "BREAST" | item == "PROSTATE") {
    Covar = c('BirthYear','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')
    CovarList = 'BirthYear+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'
  }
  
  Input[[paste(item, "_SurvTime", sep = "")]] = Input$AgeEndFollow - Input[[paste(item, "_Age", sep = "")]]
  # ColIn = append(Covar, c(item, "AgeStartFollow", "LongScore", "Death", paste(item, "_PRS", sep = ""), paste(item, "_Age", sep = ""), paste(item, "_SurvTime", sep = ""), paste(item, "_Death", sep = ""), paste(item, "_DeathUI", sep = "")))
  ColIn = append(Covar, c(item, "AgeStartFollow", "LongScore", "Death", paste(item, "_PRS", sep = ""), paste(item, "_Age", sep = ""), paste(item, "_SurvTime", sep = ""), paste(item, "_Death", sep = "")))
  ItemIn = Input[, ..ColIn]
  ItemIn = ItemIn[ItemIn[[paste(item, "_Age", sep = "")]] >= 0, ]
  ItemIn = ItemIn[complete.cases(ItemIn), ]
  nSample = nrow(ItemIn)
  nEvent = sum(ItemIn[[item]])
  
  Predictor = paste("scale(", item, "_PRS)", sep = "")
  Flag = 0
  form = paste(item," ~ ", Predictor, " + ", CovarList, sep = "")
  m = withCallingHandlers(tryCatch({glm(as.formula(form), data = ItemIn, family = binomial, na.action=na.exclude)}, 
                                   error = function(e){Flag <<- 2; return(NULL)}),
                          warning=function(w) {Flag <<- 1; invokeRestart("muffleWarning")})
  if (!is.null(m)) {
    b = summary(m)$coefficient[Predictor, "Estimate"]
    se = summary(m)$coefficient[Predictor, "Std. Error"]
    p = summary(m)$coefficient[Predictor, "Pr(>|z|)"]
    Res[nrow(Res) + 1,] = c(item, "Incidence", Predictor, "Logit", 999, "All", b, se, p, nEvent, nSample, Flag)
    rm(m, b, se, p)
  }
  else {
    Res[nrow(Res) + 1,] = c(item, "Incidence", Predictor, "Logit", 999, "All", "-", "-", "-", nEvent, nSample, Flag)
  }

  Flag = 0
  form = paste("Surv(", item, "_Age, ", item, ") ~ ", Predictor, " + ", CovarList, sep = "")
  m = withCallingHandlers(tryCatch({coxph(as.formula(form), data = ItemIn, na.action=na.exclude)}, 
                                   error = function(e){Flag <<- 2; return(NULL)}),
                          warning=function(w) {Flag <<- 1; invokeRestart("muffleWarning")})
  if (!is.null(m)) {
    b = summary(m)$coefficient[Predictor,"coef"]
    se = summary(m)$coefficient[Predictor,"se(coef)"]
    p = summary(m)$coefficient[Predictor,"Pr(>|z|)"]
    Res[nrow(Res) + 1,] = c(item, "TimeToOnset", Predictor, "CoxPh", 999, "All", b, se, p, nEvent, nSample, Flag)
    rm(m, b, se, p)
  }
  else {
    Res[nrow(Res) + 1,] = c(item, "TimeToOnset", Predictor, "CoxPh", 999, "All", "-", "-", "-", nEvent, nSample, Flag)
  }
  
  ItemIn = ItemIn[ItemIn[[item]] == 1,]
  for (SurvMax in c(2, 5, 10, 999)) {
    for (Sample in c('All', 'DiagAfterJoin')) {
      data = ItemIn[(ItemIn[[paste(item, "_SurvTime", sep = "")]] > SurvMin) & (ItemIn[[paste(item, "_SurvTime", sep = "")]] < SurvMax), ]
      if (Sample == 'DiagAfterJoin') {
        data = data[data$AgeStartFollow > data[[paste(item, "_Age", sep = "")]], ]
      }
      data = data[complete.cases(data), ]
      nSample = nrow(data)
      for (Predictor in c(paste("scale(", item, "_PRS)", sep = ""), "scale(LongScore)")) {
        # for (Target in c("Death", paste(item, "_Death", sep = ""), paste(item, "_DeathUI", sep = ""))) {
        for (Target in c("Death", paste(item, "_Death", sep = ""))) {
          nEvent = sum(data[[Target]])
          for (AnaCovar in c(CovarList, paste(CovarList, "+", item, "_Age", sep = ""))) {
            model = "CoxPh"
            if (AnaCovar == CovarList) {
              model = "CoxPh-NoDiagAge"
            }
            Flag = 0
            form = paste("Surv(", item, "_SurvTime, ", Target,") ~ ", Predictor, " + ", AnaCovar, sep = "")
            m = withCallingHandlers(tryCatch({coxph(as.formula(form), data = data, na.action=na.exclude)}, 
                                             error = function(e){Flag <<- 2; return(NULL)}),
                                    warning=function(w) {Flag <<- 1; invokeRestart("muffleWarning")})
            if (!is.null(m)) {
              b = summary(m)$coefficient[Predictor,"coef"]
              se = summary(m)$coefficient[Predictor,"se(coef)"]
              p = summary(m)$coefficient[Predictor,"Pr(>|z|)"]
              Res[nrow(Res) + 1,] = c(item, Target, Predictor, model, SurvMax, Sample, b, se, p, nEvent, nSample, Flag)
              rm(m, b, se, p)
            }
            else {
              Res[nrow(Res) + 1,] = c(item, Target, Predictor, model, SurvMax, Sample, "-", "-", "-", nEvent, nSample, Flag)
            }
            
            model = "Logit"
            if (AnaCovar == CovarList) {
              model = "Logit-NoDiagAge"
            }
            Flag = 0
            form = paste(Target," ~ ", Predictor, " + ", AnaCovar, sep = "")
            m = withCallingHandlers(tryCatch({glm(as.formula(form), data = data, family = binomial, na.action=na.exclude)}, 
                                             error = function(e){Flag <<- 2; return(NULL)}),
                                    warning=function(w) {Flag <<- 1; invokeRestart("muffleWarning")})
            if (!is.null(m)) {
              b = summary(m)$coefficient[Predictor, "Estimate"]
              se = summary(m)$coefficient[Predictor, "Std. Error"]
              p = summary(m)$coefficient[Predictor, "Pr(>|z|)"]
              Res[nrow(Res) + 1,] = c(item, Target, Predictor, model, SurvMax, Sample, b, se, p, nEvent, nSample, Flag)
              rm(m, b, se, p)
            }
            else {
              Res[nrow(Res) + 1,] = c(item, Target, Predictor, model, SurvMax, Sample, "-", "-", "-", nEvent, nSample, Flag)
            }
          }
        }
      }
      nEvent = nSample
      Target = paste(item, "_Age", sep = "")
      Predictor = paste("scale(", item, "_PRS)", sep = "")
      Flag = 0
      form = paste(Target," ~ ", Predictor, " + ", CovarList, sep = "")
      m = withCallingHandlers(tryCatch({lm(as.formula(form), data = data, na.action=na.exclude)}, 
                                       error = function(e){Flag <<- 2; return(NULL)}),
                              warning=function(w) {Flag <<- 1; invokeRestart("muffleWarning")})
      if (!is.null(m)) {
        b = summary(m)$coefficient[Predictor, "Estimate"]
        se = summary(m)$coefficient[Predictor, "Std. Error"]
        p = summary(m)$coefficient[Predictor, "Pr(>|t|)"]
        Res[nrow(Res) + 1,] = c(item, Target, Predictor, "Linear", SurvMax, Sample, b, se, p, nEvent, nSample, Flag)
        rm(m, b, se, p)
      }
      else {
        Res[nrow(Res) + 1,] = c(item, Target, Predictor, "Linear", SurvMax, Sample, "-", "-", "-", nEvent, nSample, Flag)
      }
    }
  }
}
fwrite(Res, OutFile, quote = F, row.names = F, sep = "\t")
  
  
  

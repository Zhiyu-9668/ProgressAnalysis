library(survival)
library(dplyr)
library(data.table)

# InFile = '/home/ivm/Documents/Output/ProgressionPrelim/R10/PRSexp_20220927/InterveneIn'
args = commandArgs(trailingOnly=TRUE)
InFile = args[1]
OutFile = 'ProgAnaOut'

EndPt = c('COLORECTAL', 'PROSTATE', 'BREAST', 'CHD', 'HF', 'AF', 'CKD', 'AD', 'T2D', 'STR' )
SurvMin = 0.25

Input = fread(InFile)
Input[['BirthYear']] = as.numeric(format(as.Date(Input$DateOfBirth), "%Y"))
Input[['AgeStartFollow']] = as.numeric(as.difftime(as.Date(Input$FollowUpStart)-as.Date(Input$DateOfBirth), units = "days"))/365.25
Input[['AgeEndFollow']] = as.numeric(as.difftime(as.Date(Input$FollowUpEnd)-as.Date(Input$DateOfBirth), units = "days"))/365.25
Input = Input[Input$AgeEndFollow >= 0, ]
Input = Input[Input$AgeStartFollow >= 0, ]
Input = Input[Input$AgeEndFollow >= Input$AgeStartFollow, ]

Res = data.frame(matrix(ncol = 14, nrow = 0))
colnames(Res) = c("Disease", "Target", "Predictor", "Model", "MaxSurvTime", "Sample", "AgeGrp", "AgeGrpCutoff", "Beta", "SE", "Pval", "nEvent", "nTotSample", "AnaFlag")
Covar = c('BirthYear','Sex','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')

for (item in EndPt) {
  if(item %in% colnames(Input)) {
    CovarList = 'BirthYear+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'
    if (item == "BREAST" | item == "PROSTATE") {
      CovarList = 'BirthYear+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'
    }

    Input[[paste(item, "_Age", sep = "")]] = as.numeric(as.difftime(as.Date(Input[[paste(item, "_Date", sep = "")]])-as.Date(Input$DateOfBirth), units = "days"))/365.25
    Input[[paste(item, "_SurvTime", sep = "")]] = Input$AgeEndFollow - Input[[paste(item, "_Age", sep = "")]]
    # ColIn = append(Covar, c(item, "AgeStartFollow", "LongScore", "Death", paste(item, "_PRS", sep = ""), paste(item, "_Age", sep = ""), paste(item, "_SurvTime", sep = ""), paste(item, "_Death", sep = ""), paste(item, "_DeathUI", sep = "")))
    ColIn = append(Covar, c(item, "AgeStartFollow", "AgeEndFollow", "LongScore", "Death", paste(item, "_PRS", sep = ""), paste(item, "_Age", sep = ""), paste(item, "_SurvTime", sep = ""), paste(item, "_Death", sep = "")))
    EndPtAll = Input[, ..ColIn]
    EndPtAll = EndPtAll[EndPtAll[[paste(item, "_Age", sep = "")]] >= 0, ]
    EndPtAll = EndPtAll[complete.cases(EndPtAll), ]
    
    if (item == "BREAST") {
      EndPtAll = EndPtAll[EndPtAll$Sex=="female",]
    }
    if (item == "PROSTATE") {
      EndPtAll = EndPtAll[EndPtAll$Sex=="male",]
    }
    
    Cutoff = quantile(EndPtAll$AgeEndFollow, 0.5)
    AgeGrpCutoff = paste("Median_AgeEndFollow", Cutoff, sep = "=")
    for (AgeGrp in c("All", "Lower", "Upper")) {
      if (AgeGrp == "All") {
        ItemIn = EndPtAll
      }
      if (AgeGrp == "Lower") {
        ItemIn = EndPtAll[EndPtAll$AgeEndFollow < Cutoff,]
      }
      if (AgeGrp == "Upper") {
        ItemIn = EndPtAll[EndPtAll$AgeEndFollow >= Cutoff,]
      }
      nSample = nrow(ItemIn)
      nEvent = sum(ItemIn[[item]])

      Predictor = paste("scale(", item, "_PRS)", sep = "")
      Flag = 0
      form = paste(item," ~ ", Predictor, " + ", CovarList, sep = "")
      m = withCallingHandlers(tryCatch({glm(as.formula(form), data = ItemIn, family = binomial, na.action=na.exclude)}, 
                                       error = function(e){Flag <<- 2; return(NULL)}),
                              warning=function(w) {Flag <<- 1; invokeRestart("muffleWarning")})
      if (!is.null(m)) {
        b = try(summary(m)$coefficient[Predictor, "Estimate"], silent = TRUE)
        se = try(summary(m)$coefficient[Predictor, "Std. Error"], silent = TRUE)
        p = try(summary(m)$coefficient[Predictor, "Pr(>|z|)"], silent = TRUE)
        if (class(b)=="try-error" | class(se)=="try-error" | class(p)=="try-error") {
          Res[nrow(Res) + 1,] = c(item, "Incidence", Predictor, "Logit", 999, "All", AgeGrp, AgeGrpCutoff, "-", "-", "-", nEvent, nSample, Flag)
        } else {
          Res[nrow(Res) + 1,] = c(item, "Incidence", Predictor, "Logit", 999, "All", AgeGrp, AgeGrpCutoff, b, se, p, nEvent, nSample, Flag)
        }
        rm(m, b, se, p)
      } else {
        Res[nrow(Res) + 1,] = c(item, "Incidence", Predictor, "Logit", 999, "All", AgeGrp, AgeGrpCutoff, "-", "-", "-", nEvent, nSample, Flag)
      }

      Flag = 0
      form = paste("Surv(", item, "_Age, ", item, ") ~ ", Predictor, " + ", CovarList, sep = "")
      m = withCallingHandlers(tryCatch({coxph(as.formula(form), data = ItemIn, na.action=na.exclude)}, 
                                       error = function(e){Flag <<- 2; return(NULL)}),
                              warning=function(w) {Flag <<- 1; invokeRestart("muffleWarning")})
      if (!is.null(m)) {
        b = try(summary(m)$coefficient[Predictor,"coef"], silent = TRUE)
        se = try(summary(m)$coefficient[Predictor,"se(coef)"], silent = TRUE)
        p = try(summary(m)$coefficient[Predictor,"Pr(>|z|)"], silent = TRUE)
        if (class(b)=="try-error" | class(se)=="try-error" | class(p)=="try-error") {
          Res[nrow(Res) + 1,] = c(item, "TimeToOnset", Predictor, "CoxPh", 999, "All", AgeGrp, AgeGrpCutoff, "-", "-", "-", nEvent, nSample, Flag)
        } else {
          Res[nrow(Res) + 1,] = c(item, "TimeToOnset", Predictor, "CoxPh", 999, "All", AgeGrp, AgeGrpCutoff, b, se, p, nEvent, nSample, Flag)
        }
        rm(m, b, se, p)
      } else {
        Res[nrow(Res) + 1,] = c(item, "TimeToOnset", Predictor, "CoxPh", 999, "All", AgeGrp, AgeGrpCutoff, "-", "-", "-", nEvent, nSample, Flag)
      }
    }

    EndPtCase = EndPtAll[EndPtAll[[item]] == 1,]
    Cutoff = quantile(EndPtCase[[paste(item, "_Age", sep = "")]], 0.5)
    AgeGrpCutoff = paste("Median_AgeOnset", Cutoff, sep = "=")
    for (AgeGrp in c("All", "Lower", "Upper")) {
      if (AgeGrp == "All") {
        ItemIn = EndPtCase
      }
      if (AgeGrp == "Lower") {
        ItemIn = EndPtCase[EndPtCase[[paste(item, "_Age", sep = "")]] < Cutoff,]
      }
      if (AgeGrp == "Upper") {
        ItemIn = EndPtCase[EndPtCase[[paste(item, "_Age", sep = "")]] >= Cutoff,]
      }

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
                  b = try(summary(m)$coefficient[Predictor,"coef"], silent = TRUE)
                  se = try(summary(m)$coefficient[Predictor,"se(coef)"], silent = TRUE)
                  p = try(summary(m)$coefficient[Predictor,"Pr(>|z|)"], silent = TRUE)
                  if (class(b)=="try-error" | class(se)=="try-error" | class(p)=="try-error") {
                    Res[nrow(Res) + 1,] = c(item, Target, Predictor, model, SurvMax, Sample, AgeGrp, AgeGrpCutoff, "-", "-", "-", nEvent, nSample, Flag)
                  } else {
                    Res[nrow(Res) + 1,] = c(item, Target, Predictor, model, SurvMax, Sample, AgeGrp, AgeGrpCutoff, b, se, p, nEvent, nSample, Flag)
                  }
                  rm(m, b, se, p)
                } else {
                  Res[nrow(Res) + 1,] = c(item, Target, Predictor, model, SurvMax, Sample, AgeGrp, AgeGrpCutoff, "-", "-", "-", nEvent, nSample, Flag)
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
                  b = try(summary(m)$coefficient[Predictor, "Estimate"], silent = TRUE)
                  se = try(summary(m)$coefficient[Predictor, "Std. Error"], silent = TRUE)
                  p = try(summary(m)$coefficient[Predictor, "Pr(>|z|)"], silent = TRUE)
                  if (class(b)=="try-error" | class(se)=="try-error" | class(p)=="try-error") {
                    Res[nrow(Res) + 1,] = c(item, Target, Predictor, model, SurvMax, Sample, AgeGrp, AgeGrpCutoff, "-", "-", "-", nEvent, nSample, Flag)
                  } else {
                    Res[nrow(Res) + 1,] = c(item, Target, Predictor, model, SurvMax, Sample, AgeGrp, AgeGrpCutoff, b, se, p, nEvent, nSample, Flag)
                  }
                  rm(m, b, se, p)
                } else {
                  Res[nrow(Res) + 1,] = c(item, Target, Predictor, model, SurvMax, Sample, AgeGrp, AgeGrpCutoff, "-", "-", "-", nEvent, nSample, Flag)
                }
              }
            }
          }
        }
      }
    }

    ItemIn = EndPtCase
    data = ItemIn[(ItemIn[[paste(item, "_SurvTime", sep = "")]] > SurvMin), ]
    nSample = nrow(data)
    nEvent = nSample
    Target = paste(item, "_Age", sep = "")
    Predictor = paste("scale(", item, "_PRS)", sep = "")
    Flag = 0
    form = paste(Target," ~ ", Predictor, " + ", CovarList, sep = "")
    m = withCallingHandlers(tryCatch({lm(as.formula(form), data = data, na.action=na.exclude)}, 
                                     error = function(e){Flag <<- 2; return(NULL)}),
                            warning=function(w) {Flag <<- 1; invokeRestart("muffleWarning")})
    if (!is.null(m)) {
      b = try(summary(m)$coefficient[Predictor, "Estimate"], silent = TRUE)
      se = try(summary(m)$coefficient[Predictor, "Std. Error"], silent = TRUE)
      p = try(summary(m)$coefficient[Predictor, "Pr(>|t|)"], silent = TRUE)
      if (class(b)=="try-error" | class(se)=="try-error" | class(p)=="try-error") {
        Res[nrow(Res) + 1,] = c(item, Target, Predictor, "Linear", 999, "All", "All", AgeGrpCutoff, "-", "-", "-", nEvent, nSample, Flag)
      } else {
        Res[nrow(Res) + 1,] = c(item, Target, Predictor, "Linear", 999, "All", "All", AgeGrpCutoff, b, se, p, nEvent, nSample, Flag)
      }
      rm(m, b, se, p)
    } else {
      Res[nrow(Res) + 1,] = c(item, Target, Predictor, "Linear", 999, "All", "All", AgeGrpCutoff, "-", "-", "-", nEvent, nSample, Flag)
    }

    Target = paste(item, "_Death", sep = "")
    nEvent = sum(data[[Target]])
    InterAct = paste("scale(", item, "_PRS)*", item, "_Age", sep = "")
    Flag = 0
    form = paste(Target," ~ ", InterAct, " + ", CovarList, sep = "")
    m = withCallingHandlers(tryCatch({lm(as.formula(form), data = data, na.action=na.exclude)}, 
                                     error = function(e){Flag <<- 2; return(NULL)}),
                            warning=function(w) {Flag <<- 1; invokeRestart("muffleWarning")})
    if (!is.null(m)) {
      b = try(summary(m)$coefficient[paste("scale(", item, "_PRS):", item, "_Age", sep = ""), "Estimate"], silent = TRUE) 
      se = try(summary(m)$coefficient[paste("scale(", item, "_PRS):", item, "_Age", sep = ""), "Std. Error"], silent = TRUE)
      p = try(summary(m)$coefficient[paste("scale(", item, "_PRS):", item, "_Age", sep = ""), "Pr(>|t|)"], silent = TRUE)
      if (class(b)=="try-error" | class(se)=="try-error" | class(p)=="try-error") {
        Res[nrow(Res) + 1,] = c(item, Target, InterAct, "CoxPh", 999, "All", "All", AgeGrpCutoff, "-", "-", "-", nEvent, nSample, Flag)
      } else {
        Res[nrow(Res) + 1,] = c(item, Target, InterAct, "CoxPh", 999, "All", "All", AgeGrpCutoff, b, se, p, nEvent, nSample, Flag)
      }
      rm(m, b, se, p)
    } else {
      Res[nrow(Res) + 1,] = c(item, Target, InterAct, "CoxPh", 999, "All", "All", AgeGrpCutoff, "-", "-", "-", nEvent, nSample, Flag)
    }

    Target = "Death"
    nEvent = sum(data[[Target]])
    InterAct = paste("scale(LongScore)*", item, "_Age", sep = "")
    Flag = 0
    form = paste(Target," ~ ", InterAct, " + ", CovarList, sep = "")
    m = withCallingHandlers(tryCatch({lm(as.formula(form), data = data, na.action=na.exclude)}, 
                                     error = function(e){Flag <<- 2; return(NULL)}),
                            warning=function(w) {Flag <<- 1; invokeRestart("muffleWarning")})
    if (!is.null(m)) {
      b = try(summary(m)$coefficient[paste("scale(LongScore):", item, "_Age", sep = ""), "Estimate"], silent = TRUE)
      se = try(summary(m)$coefficient[paste("scale(LongScore):", item, "_Age", sep = ""), "Std. Error"], silent = TRUE)
      p = try(summary(m)$coefficient[paste("scale(LongScore):", item, "_Age", sep = ""), "Pr(>|t|)"], silent = TRUE)
      if (class(b)=="try-error" | class(se)=="try-error" | class(p)=="try-error") {
        Res[nrow(Res) + 1,] = c(item, Target, InterAct, "CoxPh", 999, "All", "All", AgeGrpCutoff, "-", "-", "-", nEvent, nSample, Flag)
      } else {
        Res[nrow(Res) + 1,] = c(item, Target, InterAct, "CoxPh", 999, "All", "All", AgeGrpCutoff, b, se, p, nEvent, nSample, Flag)
      }
      rm(m, b, se, p)
    } else {
      Res[nrow(Res) + 1,] = c(item, Target, InterAct, "CoxPh", 999, "All", "All", AgeGrpCutoff, "-", "-", "-", nEvent, nSample, Flag)
    }
  }
}
fwrite(Res, OutFile, quote = F, row.names = F, sep = "\t")

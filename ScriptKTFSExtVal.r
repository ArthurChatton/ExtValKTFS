############################################
## Script adapted by Arthur Chatton (Université de Montréal, QC) 
## from the internship of Emilie Pilote (Université de Montréal, QC)
## Contact: see email address at https://arthurchatton.netlify.app/
############################################


##########
# Packages
##########

library(survival)
library(rms)
library(riskRegression)
library(timeROC)
library(geepack)
library(readxl)

###########
# Functions
###########

#Computing Cox model performance
perfs <- function(data, Coef, time, evt, horizon, plot=FALSE){
  
  
  # Baseline survival at 8 years
  h0 <- 0.996233
  # Design matrix of predictors
  des_matr <- as.data.frame(model.matrix(~ I(CreatD > 190) + I(AgeR > 25) + I(Ngreffe > 2) 
                                         + Rejet + I(creat3/10) + I(sqrt(Creat12)) + SexeR + Prot12 
                                         + I(Prot12^2) + SexeR:Prot12 + SexeR:I(Prot12^2), 
                                         data = data))
  des_matr$`(Intercept)` <- NULL
  
  if(Coef=="ini") {
    coef <- c(-0.76811, -0.99039,  1.07866, 0.25468, -0.03844, 0.44031,  -0.86668,  0.55057, -0.02110, 0.50685, -0.07623)
  } else if (Coef =="optim") {
    coef <- c(-0.75072, -1.02316, 1.17295, 0.22288, 0.01881, 0.41551, -0.88001, 0.61121, 0.04077, 0.48605, -0.06115)
  } else {
    stop("Invalid coefficient set")
  }
  # Prognostic index (PI)
  data$PI <- as.vector(as.matrix(des_matr) %*% cbind(coef))
  
  # Estimated absolute risk at 8 years (1 - S(t), 1 - survival at time t)
  data$pred <-  as.vector(1 - h0**exp(data$PI))
  
  
  
  data <- data[!is.na(data$PI),]
  
  surv.obj <- Surv(data[,time], data[,evt])
  
  Uno_C <- concordance(surv.obj ~ PI, 
                       data, 
                       reverse = TRUE,
                       timewt = "n/G2")$concordance
  
  AUROC <-
    timeROC::timeROC(
      T = data[,time], 
      delta = data[,evt],
      marker = data$PI,
      cause = 1, 
      weighting = "marginal", 
      times = horizon-.01,
      iid = TRUE
    )$AUC[2]
  
  # Calibration
  obj <- summary(survfit(surv.obj ~ 1, 
                         data = data), 
                 times = horizon)
  
  pred <- data$pred
  
  mean.cal <- (1 - obj$surv) / mean(pred)
  
  #Weak calibration
  lp.val <- log(-log(1 - pred))
  
  f.val <- if(any(pred==1)) {
    pred2 <- pred + 1e-7 # Handle values close to 1 
    lp.val <- log(-log(1 - pred2)) 
    f.val <- coxph(surv.obj ~ lp.val)
  } else {
    f.val <- coxph(surv.obj ~ lp.val)
  }
  
  weak.cal <- f.val$coefficients[1]
  
  
  # Moderate calibration
  # Calibration plot --------
  # Basic model
  
  data$pred.cll <- log(-log(1 - data$pred))
  
  # Estimate actual risk - basic model
  vcal <- rms::cph(surv.obj ~ rcs(pred.cll, 3),
                   x = T,
                   y = T,
                   surv = T,
                   data = data
  ) 
  
  
  dat_cal <- cbind.data.frame(
    "obs" = 1 - rms::survest(vcal,
                             times = horizon,
                             newdata = data)$surv,
    
    "lower" = 1 - rms::survest(vcal,
                               times = horizon,
                               newdata = data)$upper,
    
    "upper" = 1 - rms::survest(vcal,
                               times = horizon,
                               newdata = data)$lower,
    "pred" = data$pred
    
  )
  
  
  # Flexible calibration curve 
  dat_cal <- dat_cal[order(dat_cal$pred), ]
  
  
  if(plot){
    plot(
      dat_cal$pred, 
      dat_cal$obs,
      type = "l", 
      lty = 1, 
      xlim = c(0, 1),
      ylim = c(0, 1), 
      lwd = 2,
      xlab = "Predicted risk from developed model",
      ylab = "Predicted risk from refitted model", bty = "n"
    )
    lines(dat_cal$pred, 
          dat_cal$lower, 
          type = "l", 
          lty = 2, 
          lwd = 2)
    lines(dat_cal$pred, 
          dat_cal$upper,
          type = "l", 
          lty = 2, 
          lwd = 2)
    abline(0, 1, lwd = 1, lty = 1, col = 'gray')
    legend("bottomright",
           c("Ideal calibration",
             "Estimated calibration",
             "95% confidence interval"),
           col = c('gray', 1, 1),
           lty = c(1, 1, 2),
           lwd = c(1, 2, 2),
           bty = "n",
           cex = 0.85)
  }
  
  
  # Numerical measures ---------------
  # Basic model
  absdiff_cph <- abs(dat_cal$pred - dat_cal$obs)
  
  moderate.cal <- c(
    "ICI" = mean(absdiff_cph, na.rm=T), 
    setNames(quantile(absdiff_cph, c(0.5, 0.9), na.rm=T), c("E50", "E90"))
  )
  
  brier <-
    riskRegression::Score(list("Validation" = data$pred),
                          formula = Surv(Tretour.annee, Eretour) ~ 1, 
                          data = data, 
                          conf.int = TRUE,
                          times = horizon-.01,
                          cens.model = "km", 
                          metrics = "brier",
                          summary = "ipa"
    )$Brier$score[2,c(3,7)]
  
  brier[[2]] <- 100*brier[[2]]
  
  
  return(c(Uno_C, AUROC, brier, mean.cal, weak.cal, moderate.cal) |> setNames(nm=c("Uno_C", 'tAUROC', 'BS', 'sBS', 'mean.cal', 'weak.cal', 'ICI', 'E50', 'E90')))
  
  
}

#Computing cause-specific Cox model performance; need a CSC model called fit_csh
perfsCC <- function(data, primary_evt, horizon, plot=FALSE){ #evt is 0,1,2, primary_evt is the one of interest (usually 1), evt_name is "cens", 'return', 'death' 
 
  
  # Calculate estimated risk for each patient (in validation data) by time horizon 
  pred <- predictRisk(
    object = fit_csh, 
    cause = primary_evt, 
    newdata = data, 
    times = horizon
  )
  
  
  # pseudo observations
  score_vdata <- Score(
    list("csh_validation" = fit_csh),
    formula = Hist(Tretour.annee, evt_num) ~ 1,
    cens.model = "km",
    data = data,
    conf.int = TRUE,
    times = horizon,
    metrics = c("auc", "brier"),
    summary = c("ipa"),
    cause = primary_evt,
    plots = "calibration"
  )

  pseudos <- data.frame(score_vdata$Calibration$plotframe)
  
  
  # Calibration (O/E) -------------------------------------------------------
  
  
  # First calculate Aalen-Johansen estimate (as 'observed')
  obj <- summary(survfit(Surv(Tretour.annee, evt_name) ~ 1, data = data), times = horizon)

  # Calculate O/E
  OE <- obj$pstate[, 3] / mean(pred) # Be careful here, value "3" depends of the ordering of the levels of evt_name
  
  
  # Calibration intercept/slope ---------------------------------------------
  
  
  # Add the cloglog risk estimates to dataset with pseudo-observations
  pseudos$risk <- pseudos$risk - 0.0000001 #because some risk are == 1
  
  pseudos$cll_pred <- log(-log(1 - pseudos$risk))  
  
  # Fit model for calibration intercept
  fit_cal_int <- geese(
    pseudovalue ~ offset(cll_pred), 
    data = pseudos,
    id = ID, 
    scale.fix = TRUE, 
    family = gaussian,
    mean.link = "cloglog",
    corstr = "independence", 
    jack = TRUE
  )
  
  # Fit model for calibration slope
  fit_cal_slope <- geese(
    pseudovalue ~ offset(cll_pred) + cll_pred, 
    data = pseudos,
    id = ID, 
    scale.fix = TRUE, 
    family = gaussian,
    mean.link = "cloglog",
    corstr = "independence", 
    jack = TRUE
  )
  
  # AUC, BS, and sBS in score_vdata---------------------------------------------
  
  # Moderate calibration -------------------------------------------------------
  
  calplot_pseudo <- plotCalibration(
    x = score_vdata,
    brier.in.legend = FALSE,
    auc.in.legend = FALSE,
    cens.method = "pseudo",
    bandwidth = 0.05, # leave as NULL for default choice of smoothing
    cex = 1,
    round = FALSE, # Important, keeps all unique risk estimates rather than rounding
    xlim = c(0, 0.6),
    ylim = c(0, 0.6),
    rug = TRUE,
    xlab = "Predictions",
    bty = "n",
    plot=FALSE
  )
  
  # We can extract predicted and observed, observed will depend on degree of smoothing (bandwidth)
  dat_pseudo <- calplot_pseudo$plotFrames$csh_validation
  
  # Calculate difference between predicted and observed (make sure to use all estimated risks, not just unique ones)
  diff_pseudo <- pred - dat_pseudo$Obs[match(pred, dat_pseudo$Pred)]
  
  # Collect all numerical summary measures
  numsum_pseudo <- c(
    "ICI" = mean(abs(diff_pseudo)),
    setNames(quantile(abs(diff_pseudo), c(0.5, 0.9)), c("E50", "E90")),
    "Emax" = max(abs(diff_pseudo)),
    "Root squared bias" = sqrt(mean(diff_pseudo^2))
  )
  
  # Calibration plot
  
  if(plot){
    pseudos <- pseudos[order(pseudos$risk), ]
    
    # Use linear loess (weighted local regression with polynomial degree = 1) smoothing
    smooth_pseudos <- predict(
      stats::loess(pseudovalue ~ risk, data = pseudos, degree = 1, span = 0.33), 
      se = TRUE
    )
    
    # First, prepare histogram of estimated risks for x-axis
    spike_bounds <- c(-0.075, 0)
    bin_breaks <- seq(0, 0.6, length.out = 100 + 1)
    freqs <- table(cut(pred, breaks = bin_breaks))
    bins <- bin_breaks[-1]
    freqs_valid <- freqs[freqs > 0]
    freqs_rescaled <- spike_bounds[1] + (spike_bounds[2] - spike_bounds[1]) * 
      (freqs_valid - min(freqs_valid)) / (max(freqs_valid) - min(freqs_valid))
    
    # Produce plot
    par(xaxs = "i", yaxs = "i", las = 1)
    plot(
      x = pseudos$risk, 
      y = pseudos$pseudovalue,
      xlim = c(0, 0.6), 
      ylim = c(spike_bounds[1], 0.6),
      yaxt = "n",
      frame.plot = FALSE,
      xlab = "Estimated risks",
      ylab = "Observed outcome proportions", 
      type = "n"
    )
    axis(2, seq(0, 0.6, by = 0.1), labels = seq(0, 0.6, by = 0.1))
    polygon(
      x = c(pseudos$risk, rev(pseudos$risk)),
      y = c(
        pmax(smooth_pseudos$fit - qt(0.975, smooth_pseudos$df) * smooth_pseudos$se, 0),
        rev(smooth_pseudos$fit + qt(0.975, smooth_pseudos$df) * smooth_pseudos$se)
      ),
      border = FALSE,
      col = "lightgray"
    )
    abline(a = 0, b = 1, col = "gray")
    lines(x = pseudos$risk, y = smooth_pseudos$fit, lwd = 2)
    segments(
      x0 = bins[freqs > 0], 
      y0 = spike_bounds[1], 
      x1 = bins[freqs > 0], 
      y1 = freqs_rescaled
    )
  }
  
  
  return(c(score_vdata$AUC$score$AUC, score_vdata$Brier$score$Brier[2], score_vdata$Brier$score$IPA[2]*100, OE, summary(fit_cal_int)$mean$estimate, summary(fit_cal_slope)$mean$estimate[2], numsum_pseudo) |> setNames(nm=c('tAUROC', 'BS', 'sBS', 'mean.cal', 'cal.int', "cal.slope", 'ICI', 'E50', 'E90', 'Emax', 'RSB')))
   
}

# Non-parametric bootstrap according to the KTFS weights
boot_perfs <- function(B, data, Coef=NULL, time="Tretour.annee", evt='Eretour', horizon=8, primary_evt=NULL, evt_name=NULL){
  
  
  if(is.null(Coef)){
    res <- data.frame(matrix(nrow=B, ncol=11))
    for(b in 1:B){
      print(b)
      dboot <- data[sample(1:nrow(data), nrow(data), replace=T),]
      
      res[b,] <- unlist(perfsCC(data=dboot, horizon=horizon, primary_evt=primary_evt, plot=F))
    }
  }else{
    res <- data.frame(matrix(nrow=B, ncol=9))
    for(b in 1:B){
      print(b)
      dboot <- data[sample(1:nrow(data), nrow(data), replace=T),]
      
      res[b,] <- unlist(perfs(data=dboot, Coef=Coef, time=time, evt=evt, horizon=horizon, plot=F))
    }
  }
  
  
  return(res)
  
}

###########
# Data
###########

## DIVAT ##

# Need data and KTFS R object from Foucher et al. (2010, Kidney Int)
data.auc <- read.csv("divat.csv")

Data.auc$evt_num <- ifelse(Data.auc$Eretour==1, 1, ifelse(Data.auc$Edeces==1, 2, 0))
Data.auc$evt_name <- factor(ifelse(Data.auc$Eretour==1, "return", ifelse(Data.auc$Edeces==1, "death", "cens")))

# Cause-specific Cox model
fit_csh <- CSC(
  formula = Hist(Tretour.annee, evt_num) ~ I(CreatD>190)  + I(AgeR>25) +  I(Ngreffe>2)  + Rejet + I(creat3) + I(sqrt(Creat12)) + SexeR + Prot12 +  I(Prot12^2)  + Prot12:SexeR + I(Prot12^2):SexeR, 
  data = Data.auc
)

### EKiTE  ###
db <- read.csv("ekite_completecase.csv")

db$time <- db$time/365.25
db <- db[db$time>=1,]


names(db) <- c("X.1", "X", "id_p", "id_tr", "start_yr", "AgeR", "SexeR", "h", "w", "Ngreffe", "cit", "HLA.A", "HLA.B", "HLA.DR", "center", "dgf", "tdeath", "treturn", "tmax", "t.are", "CreatD", "age_d", "gender_d", "hrt_bs", "dcsd", "tmax2", "Tretour.annee", "Eretour", "country", "Prot12", "prt_r_3", "prt_r_6", "Creat12", "creat3", "crt_r_6", "surv_at1yr", "bmi", "incomp", "first", "Ntrans", "age_r_cl", "eGFR", "cit_hrs", "crt_ld_cl", "Rejet", "ktfs")

db$lp <- predict(modele, newdata = db)

ekite <- db[!is.na(db$lp),]

ekite$evt_num <- ifelse(ekite$Eretour==1, 1, ifelse(!is.na(ekite$tdeath), 2, 0))
ekite$evt_name <- factor(ifelse(ekite$Eretour==1, "return", ifelse(!is.na(ekite$tdeath), "death", "cens")))


ekite2 <- data.frame("SexeR" = ekite$SexeR, "SexeD" = ekite$gender_d, "IMC" = ekite$bmi, "Ngreffe" = ekite$Ngreffe, "Ischemie" = ekite$cit_hrs, "IncompABDR" = ekite$incomp, "Tretour.annee" = ekite$Tretour.annee, "AgeR" = ekite$AgeR, "AgeD" = ekite$age_d, "DGF" = ekite$dgf, "Eretour" = ekite$Eretour, "CreatD" = ekite$CreatD, "creat3" = ekite$creat3, "creat6"=ekite$crt_r_6, "Creat12"=ekite$Creat12, "prot3"=ekite$prt_r_3, "prot6" = ekite$prt_r_6, "Prot12"=ekite$Prot12, "Rejet"=ekite$Rejet, "CL12"=ekite$eGFR, "country"=ekite$country, "evt_num"=ekite$evt_num, "evt_name"=ekite$evt_name)


### CHUM ###
chum <- as.data.frame(read_excel("chum.complete.xlsx"))
chum$Eretour <- as.integer(chum$Eretour)

chum$evt_num <- ifelse(chum$Eretour==1, 1, ifelse(chum$Edeces==1, 2, 0))
chum$evt_name <- factor(ifelse(chum$Eretour==1, "return", ifelse(chum$Edeces==1, "death", "cens")))

chum$country <- "Canada"

### Pooled external validation

pool_EV <- rbind(ekite2, chum[,-21])


####
# Performance (prefix and suffix refer to the dataset/subset and the KTFS weights, respectively)
####

pool_ini <- boot_perfs(2000, pool_EV, "ini")
write.csv(pool_ini, "pool_ini.csv")

pool_optim <- boot_perfs(2000, pool_EV, "optim")
write.csv(pool_optim, "pool_optim.csv")

FR_ini <- boot_perfs(2000, pool_EV[pool_EV$country=="France",], "ini")
write.csv(FR_ini, "FR_ini.csv")

FR_optim <- boot_perfs(2000, pool_EV[pool_EV$country=="France",], "optim")
write.csv(FR_optim, "FR_optim.csv")

NO_ini <- boot_perfs(2000, pool_EV[pool_EV$country=="Norway",], "ini")
write.csv(NO_ini, "NO_ini.csv")

NO_optim <- boot_perfs(2000, pool_EV[pool_EV$country=="Norway",], "optim")
write.csv(NO_optim, "NO_optim.csv")

BE_ini <- boot_perfs(2000, pool_EV[pool_EV$country=="Belgium",], "ini")
write.csv(BE_ini, "BE_ini.csv")

BE_optim <- boot_perfs(2000, pool_EV[pool_EV$country=="Belgium",], "optim")
write.csv(BE_optim, "BE_optim.csv")

QC_ini <- boot_perfs(2000, pool_EV[pool_EV$country=="Canada",], "ini")
write.csv(QC_ini, "QC_ini.csv")

QC_optim <- boot_perfs(2000, pool_EV[pool_EV$country=="Canada",], "optim")
write.csv(QC_optim, "QC_optim.csv")



pool_cc <- boot_perfs(2000, pool_EV, primary_evt=1)
write.csv(pool_cc, "pool_cc.csv")

FR_cc <- boot_perfs(2000, pool_EV[pool_EV$country=="France",], primary_evt=1)
write.csv(FR_cc, "FR_cc.csv")

NO_cc <- boot_perfs(2000, pool_EV[pool_EV$country=="Norway",], primary_evt=1)
write.csv(NO_cc, "NO_cc.csv")

BE_cc <- boot_perfs(2000, pool_EV[pool_EV$country=="Belgium",], primary_evt=1)
write.csv(BE_cc, "BE_cc.csv")

QC_cc <- boot_perfs(2000, pool_EV[pool_EV$country=="Canada",], primary_evt=1,)
write.csv(QC_cc, "QC_cc.csv")

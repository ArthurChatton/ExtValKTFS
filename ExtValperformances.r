##################
## Script by Arthur Chatton (Université de Montréal)
## Contact: See email address at https://arthurchatton.netlify.app/
##################


# R packages
library(survminer)
library(pmcalibration)
library(riskRegression)
library(geese)
library(dcurves)


## Performance of the Cause-specific Cox model
res <- read.csv("QC_cc.csv") 
res <- res[,-1] |> setNames(nm=c('tAUROC', 'BSnull', 'BS', 'sBSnull',  'sBS', 'mean.cal', 'cal.int', "cal.slope", 'ICI', 'E50', 'E90'))

res[8] <- 1+res[8]

colMeans(res[,c(1,6,8,9,5)])
apply(res[,c(1,6,8,9,5)], 2, function(x) quantile(x, probs=c(.025,.975)))


## Cox model (initial or discrimination optimized KTFS)
res <- read.csv("pool_ini.csv")

c(mean=mean(res[,5]), quantile(res[,5], probs=c(.025,.975))) |> round(1)


## Figures

# (need pool_EV from previous script)

# survival curve
survminer::ggsurvplot(
  survfit(Surv(Tretour.annee, Eretour)~country, data=pool_EV), 
  data=pool_EV,
  add.all=TRUE,
  palette=ggokabeito::palette_okabe_ito(c(9, 1, 7, 2, 3)),
  conf.int = TRUE,
  conf.int.alpha=0.25,
  size=1,
  risk.table = TRUE,
  risk.table.y.text.col = T, 
  risk.table.y.text = FALSE,
  legend.labs = c("Pool", "Belgium", "Canada", "France", "Norway"),    
  risk.table.height = 0.35, 
  xlim = c(0,8),
  ylim=c(0.7,1),
  censor=FALSE,
  xlab="Years since transplantation",
  ggtheme = theme_bw(),
  linetype = 'strata',
  risk.table.fontsize=3
)


# Evt repartition
by(pool_EV$evt_name, pool_EV$country, function(x) prop.table(table(x)))
prop.table(table((pool_EV$evt_name)))
prop.table(table((Data.auc$evt_name)))

# Baseline survival at 8 years
h0 <- 0.996233
# Design matrix of predictors
des_matr <- as.data.frame(model.matrix(~ I(CreatD > 190) + I(AgeR > 25) + I(Ngreffe > 2) 
                                       + Rejet + I(creat3/10) + I(sqrt(Creat12)) + SexeR + Prot12 
                                       + I(Prot12^2) + SexeR:Prot12 + SexeR:I(Prot12^2), 
                                       data = pool_EV))
des_matr$`(Intercept)` <- NULL

# initial coef
coef <- c(-0.76811, -0.99039,  1.07866, 0.25468, -0.03844, 0.44031,  -0.86668,  0.55057, -0.02110, 0.50685, -0.07623)

# Prognostic index (PI)
pool_EV$PI <- as.vector(as.matrix(des_matr) %*% cbind(coef))

# Estimated absolute risk at 5 years (1 - S(t), 1 - survival at time t)
pool_EV$pred_ini <-  as.vector(1 - h0**exp(pool_EV$PI))


# auc-optim coef
coef <- c(-0.75072, -1.02316, 1.17295, 0.22288, 0.01881, 0.41551, -0.88001, 0.61121, 0.04077, 0.48605, -0.06115)
 
# Prognostic index (PI)
pool_EV$PI <- as.vector(as.matrix(des_matr) %*% cbind(coef))

# Estimated absolute risk at 5 years (1 - S(t), 1 - survival at time t)
pool_EV$pred_optim <-  as.vector(1 - h0**exp(pool_EV$PI))


# prediction for competing risk (need the fit_csh object computed in the previous script)
pool_EV$pred_cc <- predictRisk(
  object = fit_csh, 
  cause = 1, 
  newdata = pool_EV, 
  times = 8
)




# Calibration curve
mycc_ini <- get_cc(pmcalibration(
  y = Surv(pool_EV$Tretour.annee, pool_EV$Eretour),
  p = pool_EV$pred_ini,
  smooth= "rcs",
  time=8,
  ci="pw",
  transf="cloglog"
))


mycc_optim <- get_cc(pmcalibration(
  y = Surv(pool_EV$Tretour.annee, pool_EV$Eretour),
  p = pool_EV$pred_optim,
  smooth= "rcs",
  time=8,
  ci="pw",
  transf="cloglog"
))

mycc_ini$Coefficients <- "Initial"
mycc_optim$Coefficients <- "Discrimination-optimized"




score_vdata <- Score(
  list("csh_validation" = fit_csh),
  formula = Hist(Tretour.annee, evt_num) ~ 1,
  cens.model = "km",
  data = pool_EV,
  conf.int = TRUE,
  times = 8,
  metrics = c("auc", "brier"),
  summary = c("ipa"),
  cause = 1,
  plots = "calibration"
)

pseudos <- data.frame(score_vdata$Calibration$plotframe)


# Use linear loess (weighted local regression with polynomial degree = 1) smoothing
smooth_pseudos <- predict(
  stats::loess(pseudovalue ~ risk, data = pseudos, degree = 1, span = 0.33), 
  se = TRUE
)

ymin <- smooth_pseudos$fit - qt(0.975, smooth_pseudos$df) * smooth_pseudos$se
ymax <- smooth_pseudos$fit + qt(0.975, smooth_pseudos$df) * smooth_pseudos$se


mycc_cc <- data.frame(p=pseudos$risk, p_c=smooth_pseudos$fit, lower=ymin, upper=ymax, Coefficients="Competing risk")

mycc <- rbind(mycc_ini[-1,], mycc_optim[-1,], mycc_cc[,])



ggplot(mycc, aes(x = p, y = p_c, ymin=lower, ymax=upper, col=Coefficients, fill=Coefficients)) +
  geom_abline(intercept = 0, slope = 1, lty=2) +
  geom_line() +
  geom_ribbon(alpha = 1/4, linetype=0) +
  ggokabeito::scale_okabe_ito(aes = c("colour", "fill"), order=c(3,1,9)) +
  xlab("Estimated graft-failure probability") +
  ylab("Actual graft-failure probability") +
  labs_pubr() +
  theme_pubr()
 



# decision curve analysis

pool_EV$`Treat the high risk individuals` <- 1*(pool_EV$PI>4.17)

#renaming for a beautiful plot
pool_EV$`Clinical prediction model` <- pool_EV$pred_ini

dca(Surv(Tretour.annee, Eretour) ~ `Clinical prediction model`  + `Treat the high risk individuals`, 
    data = pool_EV[pool_EV$country=="France",],
    time = 8,
    thresholds = 1:100 / 100) |>
  plot(smooth = TRUE)

pool_EV$`Clinical prediction model` <- pool_EV$pred_optim

dca(Surv(Tretour.annee, Eretour) ~ `Clinical prediction model`  + `Treat the high risk individuals`, 
    data = pool_EV[pool_EV$country=="France",],
    time = 8,
    thresholds = 1:100 / 100) |>
  plot(smooth = TRUE)


pool_EV$`Treat using the KTFS` <- ifelse(pool_EV$pred_cc>1, 1, pool_EV$pred_cc) #need to deal with values slightly greater than 1 due to comp risk
pool_EV$evt_name <- plyr::revalue(pool_EV$evt_name, c("cens" = "censor"))
pool_EV$evt_name <- forcats::fct_relevel(pool_EV$evt_name, c("return", "death"), after=1)

dca(Surv(Tretour.annee, evt_name) ~ `Treat using the KTFS` + `Treat the high risk individuals`, 
              data = pool_EV,
              time = 8,
              thresholds = 1:100 / 100) |>
  plot(smooth = TRUE)

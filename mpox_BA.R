library(readxl)
library(dplyr)
library(tidyr)
library(confintr) #cramersv
library(ggplot2)
library(plot.matrix)
library(viridis)
library(caret)
library(glmnet)
library(logisticRR)
setwd("/home/mo/IT/BA")


data <- read_excel("data/MPox_final_merged_for-modelling_08032024.xlsx") %>% as_tibble()
#sum(data$MPX_DIAGNOSE == 1)
# -> 47 Patienten (~6,5%) hatten bereits eine Mpox Diagnose
#restliche daten über R terminal erhoben 

newdata <- read_excel("newdata/newdataMay.xlsx") %>% as_tibbl
#new data contains some additional column for predicted cases

betterdata <- read_excel("/home/mo/IT/BA/betterdata/MPox_final_merged_for-modelling_27052024.xlsx")

#############cleaning up data#########

usable_data <- data %>%
  #select(! stat_reported:prediction_all) %>% 
  select(! ID:Result)%>%
  select(! c(MPX_NICHTBEREITSCHAFT_M5, MPX_NICHTBEREITSCHAFT_M6,SEX_LOCATIONS_M8,SEX_LOCATIONS_M9))%>%
  select(! c(HILF1,HILF2,HILF3)) %>%
  rowwise()%>%
  mutate(across(everything(), ~replace(.x, .x == 777, 9)))%>%
  mutate( VAC_HESITATION_NO_RISK = as.integer(any(c_across(18:21) == 1, na.rm = TRUE)),
          VAC_HESITATION_COMPLICATED = as.integer(any(c_across(18:21) == 2, na.rm = TRUE)),
          VAC_HESITATION_SIDE_EFFECTS = as.integer(any(c_across(18:21) == 3, na.rm = TRUE)),
          VAC_HESITATION_DIAGNOSED = as.integer(any(c_across(18:21) == 4, na.rm = TRUE)),
          VAC_HESITATION_ENDEMIC = as.integer(any(c_across(18:21) == 5, na.rm = TRUE)),
          VAC_HESITATION_IF_SPREADS = as.integer(any(c_across(18:21) == 6, na.rm = TRUE)),
          SEX_LOCATIONS_PRIVATE = as.integer(any(c_across(45:52) == 1, na.rm = TRUE)),
          SEX_LOCATIONS_CLUBS = as.integer(any(c_across(45:52) == 2, na.rm = TRUE)),
          SEX_LOCATIONS_SAUNAS = as.integer(any(c_across(45:52) == 3, na.rm = TRUE)),
          SEX_LOCATIONS_PRIV_PARTIES = as.integer(any(c_across(45:52) == 4, na.rm = TRUE)),
          SEX_LOCATIONS_PUBLIC_PARTIES = as.integer(any(c_across(45:52) == 5, na.rm = TRUE)),
          SEX_LOCATIONS_CRUISING = as.integer(any(c_across(45:52) == 7, na.rm = TRUE)),
          SEX_LOCATIONS_CINEMA = as.integer(any(c_across(45:52) == 9, na.rm = TRUE)),
          SEX_LOCATIONS_ELSWHERE = as.integer(any(c_across(45:52) == 8, na.rm = TRUE)),
          SEX_LOCATIONS__ = as.integer(any(c_across(45:52) == 6, na.rm = TRUE))) %>% 
  ungroup() %>% 
  mutate( SEROLOGICAL_PREDICTION = data$prediction_all,
          DIAGNOSED_CASE = data$MPX_DIAGNOSE,
          SUSPECTED_CASE = newdata$mpox_suspected,
          MPX_ALL_CASES = ((data$MPX_DIAGNOSE)==1 | replace_na(betterdata$mpox_suspected,0) == 1)*1) %>% 
  select(!c(MPX_NICHTBEREITSCHAFT, MPX_NICHTBEREITSCHAFT_M2, MPX_NICHTBEREITSCHAFT_M3, MPX_NICHTBEREITSCHAFT_M4, #was recoded
            RADIO_2,O_RADIO_2,         #free text -> useless
            SEX_LOCATIONS,SEX_LOCATIONS_M2, SEX_LOCATIONS_M3, SEX_LOCATIONS_M4, SEX_LOCATIONS_M5, SEX_LOCATIONS_M6, SEX_LOCATIONS_M7, #"was recoded
            stat_reported, stat_measured, prediction_thresh,serostatus.delta, predicted, prediction_all, #Not properly interpretable
            SEX_PARTNERS,SEX_PARTNERS_M2, SEX_PARTNERS_M3)) #implicitly in the amount of partners, would need recoding anyways

    
  
print_colnames <-  function(data){
  for(i in 1:ncol(data)){
    cat(i," ",colnames(data)[i], "\n")
  }
}
print_colnames(usable_data)
### some translating to english, some names just generally made more clear
colnames(usable_data)[c(4:29,33,34,35)] <- c("KNOWLEDGE_MOSTLY_MSM",
                                    "KNOWLEDGE_STI",
                                    "KNOWLEDGE_MPOX_VACC",
                                    "MPX_DIAGNOSIS",
                                    "MPX_DIAGOSIS_YEAR",
                                    "MPX_DIAGNOSIS_MONTH",#9
                                    "MPX_SYMPTOMS",
                                    "MPX_CONTACTS",#11
                                    "MPX_SORROW",
                                    "MPX_VACCINATION",
                                    "MPX_VACCINATION_AMOUNT",
                                    "MPX_VACCINATION_YEAR",
                                    "MPX_VACCINATION_MONTH",#16
                                    "MPX_VACCINATION_WILLING",
                                    "POX_VACCINATION",
                                    "POX_VACCINATION_AMOUNT",
                                    "OPINION_VACCINATION",#20
                                    "BEHAVIOUR_MODIFIED",
                                    "REDUCED_PARTNERS",
                                    "REDUCED_LOCATIONS",
                                    "REDUCED_ANAL",
                                    "REDUCED_ORAL",
                                    "REDUCED_VAGINAL",
                                    "MORE_CONDOMS",
                                    "MORE_CONVERSATION",
                                    "MORE_EXCHANGE_CONTACTS",
                                    "SEX_PARTNERS_FEMALE_ANAL",
                                    "SEX_PARTNERS_NB",
                                    "SEX_PARTNERS_NB_ANAL")

#infected_data <- usable_data[usable_data$MPX_ALL_CASES == 1,]

  

#####correlation array plot############
array_plot <- function(usable_data, title = "", mode = "cramers v", alpha = 0.05){
nparam <- ncol(usable_data)
usable_data_factor_df <- as.data.frame(usable_data) %>%
  lapply(function(x) ifelse(is.na(x),"NA", x)) %>%
  lapply(factor)%>%
  as.data.frame()
significant <- nparam
insignificant <- 0
cor_matrix <- matrix(rep(0, nparam**2), nrow = nparam)
if(mode == "cramers v"){
  for( i in 1:nparam){
    for(j in 1:nparam){
      if(i < j){
        cor_matrix[i,j] <- cramersv(usable_data_factor_df[,c(i,j)])
      }else if(i == j){
        cor_matrix[i,i] <- 1
      }else{
        cor_matrix[i,j] <- cor_matrix[j,i]
      }
    }
  }
}else if(mode == "chi squared test"){
  for( i in 1:nparam){
    for(j in 1:nparam){
      if(i < j){
        if(chisq.test(usable_data_factor_df[,i],usable_data_factor_df[,j])$p.value < 0.05) {
          cor_matrix[i,j] <- 1
          significant <- significant + 1
        } else{
          cor_matrix[i,j] <- 0
          insignificant <- insignificant + 1
          }
      }else if(i == j){
        cor_matrix[i,i] <- 1
      }else{
        cor_matrix[i,j] <- cor_matrix[j,i]
      }
    }
  }
} else{ cat("falsche eingabe")
  return
}
colnames(cor_matrix) <- colnames(usable_data)
rownames(cor_matrix) <- colnames(usable_data)
par(mar = c(8,8,4,6))
par(las=2)
par(cex.axis = 0.5)
plot(cor_matrix, breaks = 25, xlab = "", ylab = "", col = viridis(25), main = title)
if(mode == "chi squared test"){cat("significant:", significant, " insignificant:", insignificant)}
}


array_plot(usable_data)
array_plot(usable_data, mode = "chi squared test")


#######  prepearing data for regression  ###########

regression_data <- tibble(
  AGE_LOW = usable_data$AGE == 1,
  AGE_HIGH = usable_data$AGE >= 3,
  #gender and sex skipped, lets subset msm later :)
  KNOWLEDGE_MOSTLY_MSM_NEGATIVE = usable_data$KNOWLEDGE_MOSTLY_MSM != 1,
  KNOWLEDGE_STI_NEGATIVE = usable_data$KNOWLEDGE_STI != 1,
  KNOWLEDGE_MPOX_VACC_NEGATIVE = usable_data$KNOWLEDGE_MPOX_VACC != 1,
  #Wissen vereinfacht, Teilnehmer hat NICHT "I already knew that" angekreuzt.
  #mit drin gelassen auch wenn es vermutlich nix zu entdecken gibt
  MPX_CONTAKT_NO_SEX = data$MPX_KONTAKTE == 2,
  MPX_CONTAKT_SEX = data$MPX_KONTAKTE == 3,
  #MPX_SORGE_HIGH = data$MPX_SORGE == 3 | data$MPX_SORGE == 4,
  #MPX_SORGE_NONE = data$MPX_SORGE == 1,
  #sorge isn't a good predictor as its only for those without mpox diagnosis, maybe llater regression with diagnosed cases excluded?
  #VACCINATION STUFF SKIPPED,
  #inclusive NICHTBEREITSCHAFT
  #inclusive POX, wie hier imputieren (unreliable answers -> dont worry about it)?
  SEX_LOCATIONS_PRIVATE  = usable_data$SEX_LOCATIONS_PRIVATE,
  SEX_LOCATIONS_CRUISING = usable_data$SEX_LOCATIONS_CRUISING,
  SEX_LOCATIONS_PRIVATE_PARTIES = usable_data$SEX_LOCATIONS_PRIV_PARTIES,
  SEX_LOCATIONS_PUBLIC_PARTIES = usable_data$SEX_LOCATIONS_PUBLIC_PARTIES,
  SEX_LOCATIONS_CLUBS = usable_data$SEX_LOCATIONS_CLUBS,
  SEX_LOCATIONS_SAUNAS = usable_data$SEX_LOCATIONS_SAUNAS,
  SEX_LOCATIONS_ELSEWHERE = usable_data$SEX_LOCATIONS_ELSWHERE,
  SEX_LOCATIONS__ = usable_data$SEX_LOCATIONS__,
  #alternativ: < x bzw > x Orte besucht mit -1 / 1 codiert?
  SEX_PARTNERS_FEMALE = replace_na(usable_data$SEX_PARTNERS_FEMALE,0) > 1,
  SEX_PARTNERS_NB_TRANS = replace_na(usable_data$SEX_PARTNERS_NB,0) > 1,
  #bei weiblichen / nichtbinären / trans Partnern wird anzahl außen vor gelassen, und zunächst getestet, ob es überhaupt Zusammenhang gibt
  SEX_PARTNERS_MALE_HIGH = replace_na(usable_data$SEX_PARTNERS_MALE,0) > quantile(replace_na(data$SEX_PARTNERS_MALE,0),0.75),
  SEX_PARTNERS_MALE_LOW = replace_na(usable_data$SEX_PARTNERS_MALE,0) < quantile(replace_na(data$SEX_PARTNERS_MALE,0),0.25),
  SEX_PARTNERS_MALE_ANAL_HIGH = replace_na(usable_data$SEX_PARTNERS_MALE_ANAL,0) > quantile(replace_na(data$SEX_PARTNERS_MALE_ANAL,0),0.75),
  SEX_PARTNERS_MALE_ANAL_LOW = replace_na(usable_data$SEX_PARTNERS_MALE_ANAL,0) < quantile(replace_na(data$SEX_PARTNERS_MALE_ANAL,0),0.25),
  #Imputuation evtl nicht angabracht (na -> 0), codierung experimentell
  HIV_NO_TEST = replace_na(usable_data$HIV_TEST,0) == 2,
  HIV_DIAGNOSIS_POS = replace_na(usable_data$HIV_DIAGNOSIS,0) == 1,
  PREP_USAGE = replace_na(usable_data$PREP,0) > 1 &  replace_na(usable_data$PREP,0) < 5,
  #Verhaltensänderung no/notanymore vertauscht
  CHANGE_GENERAL_YES = usable_data$BEHAVIOUR_MODIFIED == 1,
  CHANGE_GENERAL_NOT_ANYMORE = usable_data$BEHAVIOUR_MODIFIED == 2,
  REDUCED_PARTNERS = replace_na(usable_data$REDUCED_PARTNERS,0) == 1 | replace_na(usable_data$REDUCED_PARTNERS,0) == 2,
  REDUCED_LOCATIONS = replace_na(usable_data$REDUCED_LOCATIONS,0) == 1 | replace_na(usable_data$REDUCED_LOCATIONS,0) == 2,
  REDUCED_ANAL = replace_na(usable_data$REDUCED_ANAL,0) == 1 | replace_na(usable_data$REDUCED_ANAL,0) == 2,
  REDUCED_ORAL = replace_na(usable_data$REDUCED_ORAL,0) == 1 | replace_na(usable_data$REDUCED_ORAL,0) == 2,
  REDUCED_VAGINAL = replace_na(usable_data$REDUCED_VAGINAL,0) == 1 | replace_na(usable_data$REDUCED_VAGINAL,0) == 2,
  MORE_CONDOMS = replace_na(usable_data$MORE_CONDOMS,0) == 1,
  MORE_CONVERSATION = replace_na(usable_data$MORE_CONVERSATION,0) == 1,
  MORE_CONTACTINFO = replace_na(usable_data$MORE_EXCHANGE_CONTACTS,0) == 1
  #,MANY_LOCATIONS = replace_na(data$SEX_LOCATIONS_M2,0) > 0
) %>% 
  filter(!is.na(data$SEX_LOCATIONS) & usable_data$ASSIGNED_SEX == 1 & usable_data$AGE < 4) %>% 
  mutate(across(.cols = everything(), .fns = as.integer)) %>%
  as.data.frame()

#MPOX_CASES <- (usable_data$MPX_ALL_CASES == 1 | replace_na(betterdata$mpox_suspected,0) == 1)*1
MPOX_CASES <- usable_data$MPX_ALL_CASES[!is.na(data$SEX_LOCATIONS) & usable_data$ASSIGNED_SEX == 1 & usable_data$AGE < 4]

array_plot(mutate(regression_data, MPOX_CASES = MPOX_CASES))
array_plot(mutate(regression_data, MPOX_CASES = MPOX_CASES), mode = "chi squared test")


high_sex_data <- regression_data[regression_data$SEX_PARTNERS_MALE_ANAL_HIGH == 1,] %>% 
  select(!c(SEX_PARTNERS_MALE_ANAL_HIGH, SEX_PARTNERS_MALE_ANAL_LOW, SEX_PARTNERS_MALE_HIGH,SEX_PARTNERS_MALE_LOW))

high_sex_cases <- MPOX_CASES[regression_data$SEX_PARTNERS_MALE_ANAL_HIGH == 1]

array_plot(high_sex_data, mode = "chi squared test")
array_plot(test, mode = "chi squared test")

test <- regression_data[sample(1:600, 100),]

#infected_data <- usable_data[usable_data$MPX_ALL_CASES == 1,]

riskRatioPlot <- function(data, response, nboot = 100, adjusted = FALSE){
    #rr_vector <- rep(0, ncol(data))
  rr_boot_df <-tibble(values = numeric(nboot * (ncol(data))), 
                      names = character(nboot * (ncol(data))))
  if(adjusted){
    for(i in 1:ncol(data)){
      formula <- as.formula(paste("y ~ ", colnames(data)[i])) 
      rr_object <- logisticRR(formula = formula, data = mutate(data, y = response), boot = TRUE, n.boot = nboot)
      #rr_vector[i] <- rr_object$RR
      rr_boot_df$values[(1:nboot)+ (i-1)*nboot] <- rr_object$boot.rr #saving bootstrapped riskratio in "flattened matrix" vector
      rr_boot_df$names[(1:nboot)+ (i-1)*nboot] <- rep(colnames(data)[i], nboot) #saving the name of the parameter tied to every bootstrapped value for the plotting function
      print(paste(i, colnames(data)[i],": ", rr_object$RR))
    }
  }else{ #cude risk ratios
    for (i in 1:ncol(data)) {
       data2col <- tibble(x = data[,i], y= response)
       for(j in 1:nboot){
         bootstrapped_data <- sample_n(data2col,nrow(data2col), replace = TRUE) #jth bootstrap
         exposed <- bootstrapped_data[bootstrapped_data$x == 1,]
         unexposed <- bootstrapped_data[bootstrapped_data$x == 0,]
         rr_boot_df$values[(i-1)*nboot + j] <- (sum(exposed$y)/nrow(exposed))/(sum(unexposed$y)/nrow(unexposed))  #rr formula
         rr_boot_df$names[(i-1)*nboot + j] <- colnames(data)[i]
       }
       print(paste(i, colnames(data)[i])) #supervision of process
    }
  }
  #ylabtitle = if(adjusted){"Ajusted Risk Ratio"} else {"Crude Risk Ratio"}
  ggplot(rr_boot_df, aes(x = names, y = values)) +
    geom_boxplot()+
    coord_flip()+
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    labs(x = "", y = "")+
    theme_minimal()+
    theme(axis.text.x = element_text(size = 12))
}
set.seed(666)
riskRatioPlot(regression_data, MPOX_CASES, adjusted = FALSE, nboot = 1000)

########MODEL CREATION############
trainControl <- trainControl(method = "repeatedcv",
                             number = 15,
                             repeats = 5,
                             summaryFunction = twoClassSummary,
                             classProbs = TRUE)  #mehr 


alphavec <- seq(0,1, by = 0.05)
lambdavec <- c(0.005, 0.0075, 0.01, 0.02, 0.03, 0.04,0.05,0.06, 0.07,0.08, 0.09, 0.1,0.11,0.12,0.13,0.14,0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5)
lambdavec2 <- seq(from = 0, to = 1.5, by = 0.005)
tuneGrid <- expand.grid(alpha = alphavec, lambda = lambdavec)

set.seed(1)

model <- train(regression_data , as.factor(c("gesund","krank")[MPOX_CASES +1]),
               method = "glmnet",
               metric = "ROC",
               tuneGrid = tuneGrid,
               trControl = trainControl
)
model$bestTune

rocMatrix <- matrix(model$results$ROC, byrow = TRUE, ncol = length(lambdavec), 
                    dimnames = list(paste(alphavec),
                                    paste(lambdavec)))

plot(rocMatrix, breaks = 50, xlab = "lambda", ylab = "alpha", col = viridis(50), main = "")




set.seed(1)
ridgemodel <- train(regression_data , as.factor(c("gesund","krank")[MPOX_CASES +1]),
                   method = "glmnet",
                   metric = "ROC",
                   tuneGrid = expand.grid(alpha = 0, lambda = lambdavec2),
                   trControl = trainControl
)
ridge_params <- ridgemodel$finalModel %>% coef(ridgemodel$bestTune$lambda)

relaxed_model <- train(regression_data , as.factor(c("gesund","krank")[MPOX_CASES +1]),
                       method = "glmnet",
                       metric = "ROC",
                       tuneGrid = expand.grid(alpha = 0.1, lambda = lambdavec2),
                       trControl = trainControl
)
relaxed_params <- relaxed_model$finalModel %>% coef(relaxed_model$bestTune$lambda)


halfmodel <- train(regression_data, as.factor(c("gesund","krank")[MPOX_CASES +1]),
                   method = "glmnet",
                   metric = "ROC",
                   tuneGrid = expand.grid(alpha = 0.5, lambda = lambdavec2),
                   trControl = trainControl
)
half_parms <- halfmodel$finalModel %>% coef(halfmodel$bestTune$lambda)

conservativemodel  <- train(regression_data, as.factor(c("gesund","krank")[MPOX_CASES +1]),
                            method = "glmnet",
                            metric = "ROC",
                            tuneGrid = expand.grid(alpha = 0.9, lambda = lambdavec2),
                            trControl = trainControl
)
conservative_params <- conservativemodel$finalModel %>% coef(conservativemodel$bestTune$lambda)



lassomodel  <- train(regression_data, as.factor(c("gesund","krank")[MPOX_CASES +1]),
                            method = "glmnet",
                            metric = "ROC",
                            tuneGrid = expand.grid(alpha = 1, lambda = lambdavec2),
                            trControl = trainControl
)
lassoparams <- lassomodel$finalModel %>% coef(lassomodel$bestTune$lambda)

####prints latex code 
ro <- 3
for(i in 1:(ncol(regression_data)+1)) {
  cat(gsub("_", "\\\\_",c("(intercept)",colnames(regression_data))[i]), "&$", round(ridge_params[i], ro), "$&$", round(relaxed_params[i],ro), "$&$", round(half_parms[i],ro), "$&$", round(conservative_params[i],ro), "$&$", round(lassoparams[i],ro), "$\\\\ \\hline\n")
}



relaxed_param_matrix <- relaxed_model$finalModel %>% coef(s = lambdavec2)
paramcounts <- rep(0, length(lambdavec2))
for(i in 1:length(lambdavec2)){
  paramcounts[i] = sum(relaxed_param_matrix[,i] != 0) -1
}

ggplot(tibble(lambda = lambdavec2, AUC = relaxed_model$results$ROC, nparams = paramcounts/100 +0.4))+
  geom_point( aes(lambda,AUC),color = "black") +
  geom_line(aes(lambda, nparams), color = "red")+
  scale_y_continuous(sec.axis = sec_axis(~ (. - 0.4) * 100, name = "Amount nonzero parameters"))+ 
  theme_minimal()+
  theme(axis.title = element_text(size = 25),   # Increase size of axis labels
        axis.text = element_text(size = 25)) 




ggplot(tibble(lambda = lambdavec2, AUC = relaxed_model$results$ROC, nparams = paramcounts), aes(lambda, AUC)) +
  geom_point(color = "black") +
  scale_y_continuous(sec.axis = sec_axis(~ . * 60, name = "Repurchase Premium", breaks = pretty(relaxed_model$results$ROC), labels = pretty(relaxed_model$results$ROC))) + 
  theme_minimal() +
  theme(axis.title = element_text(size = 25),   # Increase size of axis labels
        axis.text = element_text(size = 25)) 

palet <- rainbow(n = ncol(regression_data))
glmnet_lasso_model <-glmnet(regression_data, MPOX_CASES, alpha = 1, family = "binomial")
coef(glmnet_lasso_model,s = lassomodel$bestTune$lambda)
plot(glmnet_lasso_model, xvar = "lambda", label = TRUE, col = palet) 
legend("topleft", legend = colnames(regression_data), col = palet, lty = 1, cex = 0.4)




set.seed(1)
high_sex_model <- train(high_sex_data , as.factor(c("gesund","krank")[high_sex_cases +1]),
                     method = "glmnet",
                     metric = "ROC",
                     tuneGrid = tuneGrid,
                     trControl = trainControl
)

high_sex_rocMatrix <- matrix(high_sex_model$results$ROC, byrow = TRUE, ncol = length(lambdavec), 
                             dimnames = list(paste(alphavec),
                                             paste(lambdavec)))

plot(high_sex_rocMatrix, breaks = 50, xlab = "lambda", ylab = "alpha", col = viridis(50), main = "")




high_sex_model$finalModel %>% coef(high_sex_model$bestTune$lambda)


set.seed(1)
high_sex_ridge <- train(high_sex_data , as.factor(c("gesund","krank")[high_sex_cases +1]),
                        method = "glmnet",
                        metric = "ROC",
                        tuneGrid = expand.grid(alpha = 0.0, lambda= lambdavec2),
                        trControl = trainControl
)
high_sex_ridge_params <- high_sex_ridge$finalModel %>% coef(s = high_sex_ridge$bestTune$lambda)

high_sex_relaxed <-  train(high_sex_data , as.factor(c("gesund","krank")[high_sex_cases +1]),
                                          method = "glmnet",
                                          metric = "ROC",
                                          tuneGrid = expand.grid(alpha = 0.1, lambda= lambdavec2),
                                          trControl = trainControl
)
high_sex_relaxed_params <- high_sex_relaxed$finalModel %>% coef(s = high_sex_relaxed$bestTune$lambda)

high_sex_half <-  train(high_sex_data , as.factor(c("gesund","krank")[high_sex_cases +1]),
                           method = "glmnet",
                           metric = "ROC",
                           tuneGrid = expand.grid(alpha = 0.1, lambda= lambdavec2),
                           trControl = trainControl
)
high_sex_half_params <- high_sex_half$finalModel %>% coef(s = high_sex_half$bestTune$lambda)

high_sex_conservative <-  train(high_sex_data , as.factor(c("gesund","krank")[high_sex_cases +1]),
                        method = "glmnet",
                        metric = "ROC",
                        tuneGrid = expand.grid(alpha = 0.1, lambda= lambdavec2),
                        trControl = trainControl
)
high_sex_conservative_params <- high_sex_conservative$finalModel %>% coef(s = high_sex_conservative$bestTune$lambda)

high_sex_lasso <-  train(high_sex_data , as.factor(c("gesund","krank")[high_sex_cases +1]),
                        method = "glmnet",
                        metric = "ROC",
                        tuneGrid = expand.grid(alpha = 0.1, lambda= lambdavec2),
                        trControl = trainControl
)
high_sex_lasso_params <- high_sex_lasso$finalModel %>% coef(s = high_sex_lasso$bestTune$lambda)








high_sex_rocMatrix <- matrix(high_sex_model$results$ROC, byrow = TRUE, ncol = length(lambdavec), 
                          dimnames = list(paste(alphavec),
                                          paste(lambdavec)))

plot(high_sex_rocMatrix, breaks = 50, xlab = "lambda", ylab = "alpha", col = viridis(50), main = "")


ro <- 3
for(i in 1:(ncol(high_sex_data)+1)) {
  cat(gsub("_", "\\\\_",c("(intercept)",colnames(high_sex_data))[i]), "&$", round(high_sex_ridge_params[i], ro), "$&$", round(high_sex_relaxed_params[i],ro), "$&$", round(high_sex_half_params[i],ro), "$&$", round(high_sex_conservative_params[i],ro), "$&$", round(high_sex_lasso_params[i],ro), "$\\\\ \\hline\n")
}





#elbowplot
par(cex.lab = 1.5,  # Axis labels size
    cex.axis = 1.3) # Tick labels size
elbowmodel <- cv.glmnet(as.matrix(regression_data), MPOX_CASES, family="binomial", type.measure="auc", nfolds= 50, gamma = 1)
plot(elbowmodel)

set.seed(1)
elbowrelaxedmodel <- cv.glmnet(as.matrix(regression_data), MPOX_CASES, family="binomial", type.measure="auc", nfolds= 50, gamma = 0.1)


rocMatrix <- matrix(model$results$ROC, byrow = TRUE, ncol = length(lambdavec), 
                    dimnames = list(paste(alphavec),
                                    paste(lambdavec)))

plot(rocMatrix, breaks = 50, xlab = "lambda", ylab = "alpha", col = viridis(50), main = "")


################## MODEL PERFORMANCE ################
#2 violins visualisation
violin_data <- tibble(
  MPOX_CASES = as.factor(high_sex_cases),
  MPOX_PREDICTION = predict(high_sex_model, newdata = regression_data, type = "prob")$krank
)
ggplot(violin_data, aes(x = MPOX_CASES, y = MPOX_PREDICTION)) + 
  geom_violin() 

#contingency table + sensitivity + specificity 
contingency_table <- function(inputmodel, inputdata, truth, cutoff){
  prediction  <- predict(inputmodel, newdata = inputdata, type = "prob")$krank
  print(
  matrix(c(sum(truth == 1 & prediction > cutoff), sum(truth == 0 & prediction > cutoff),
           sum(truth == 1 & prediction<= cutoff), sum(truth == 0 & prediction <= cutoff)), 
         byrow = TRUE, nrow = 2,  dimnames= list(predicted = c("pred_case", "pred_nocase") , truth = c("case", "no_case"))))
  cat("Specififcity:", sum(truth == 0 & prediction <= cutoff)/sum(truth == 0), "\n",
      "Sensitivity:", sum(truth == 1 & prediction > cutoff)/ sum(truth == 1))
}
contingency_table(model, regression_data, MPOX_CASES, 0.2)

contingency_table(high_sex_model, regression_data %>%  select(!c(SEX_PARTNERS_MALE_ANAL_HIGH, SEX_PARTNERS_MALE_ANAL_LOW, SEX_PARTNERS_MALE_HIGH,SEX_PARTNERS_MALE_LOW)), MPOX_CASES, 0.38)


tibble(truth = MPOX_CASES,  prediction = qlogis(predict(model, newdata = regression_data, type = "prob")$krank)) %>%
  ggplot(aes(x= prediction, y = truth)) +
  geom_point(alpha = 0.1, color = "blue") +
  stat_function(fun = plogis, xlim = c(-4, 2))




library(proto)
library(dplyr)
library(splitstackshape)
library(grf)
library(ggplot2)

# source('marginal_learning.R')
# ################################## NOT NEEDED AMBLER ##########################
prepare_Ambler_data <- function() {
  setwd('/Users/dianamirzoyan/Desktop/Diana/DISSERTATION/data/Ambler/Remittances Dataverse Files/data')
  final_data <- read.csv(
    'Ambler_data.csv'
  )
  stratify_data <- data.frame(treatment=final_data$treatment, cid=final_data$dayXlocation) # $1:length(
  pretreatment_vars <- final_data[
    c('us_stu_insch', 'us_stu_edlvl', 'us_stu_grade', 'us_stu_edlvl_all',
      'us_stu_female', 'us_stu_mig_rel', 'us_stu_age', 'us_stu_grade_all',
      'us_stu_yearsofed', 'us_stu_mig_rel_child', 'us_stu_mig_rel_sib',
      'us_stu_mig_rel_niece', 'us_stu_mig_rel_cousin')
  ]

  list(
    pretreatment_vars=pretreatment_vars,
    Y2=final_data$costs,
    Y1=final_data$er_amount_total,
    cluster=final_data$dayXlocation,
    treatment=final_data$treatment,
    stratify_data=stratify_data
  )
}

prepare_Ambler_data_upd <- function() {
  setwd('/Users/dianamirzoyan/Desktop/Diana/DISSERTATION/data/Ambler/Remittances Dataverse Files/data')
  final_data <- read.csv(
    'Ambler_data_upd.csv'
  )
  final_data <- na.omit(final_data)
  stratify_data <- data.frame(treatment=final_data$treatment_tmp2, cid=final_data$dayXlocation) # $1:length(
  pretreatment_vars <- final_data[
    c('us_stu_insch', 'us_stu_edlvl', 'us_stu_grade', 'us_stu_edlvl_all',
      'us_stu_female', 'us_stu_mig_rel', 'us_stu_age', 'us_stu_grade_all',
      'us_stu_yearsofed', 'us_stu_mig_rel_child', 'us_stu_mig_rel_sib',
      'us_stu_mig_rel_niece', 'us_stu_mig_rel_cousin')
  ]

  list(
    pretreatment_vars=pretreatment_vars,
    Y2=final_data$costs,
    Y1=final_data$er_amount_total,
    cluster=final_data$dayXlocation,
    treatment=final_data$treatment_tmp2,
    stratify_data=stratify_data
  )
}
################################## NOT NEEDED AMBLER ##########################

cumulative_te <- function(data, order_var) { # weights! + naming + more vars. vs linear algebra?
  data %>% arrange(desc({{ order_var }})) %>%
    mutate(
      'cumulative_Y1_{{ order_var }}' := cumsum(Y1*W/sum(W) - Y1*(1-W)/sum(1-W)),
      'cumulative_Y2_{{ order_var }}' := cumsum(Y2*W/sum(W) - Y2*(1-W)/sum(1-W)))
}


############################# АЛЬТЕРНАТИВНЫЙ ВАРИАНТ CUM SUM ###############
prepare_for_cum_te <- function(data) {
  n_s <- nrow(data)
  print(n_s)
  n_out <- n_s #min(2000, n_s)
  beta <- dbeta( # MAKE IT LAZY with a cache + make the matrix sparse and compute lazilly
    # JUST COMPUTE EFFICIENTLY. no need to devide by BETA
    # however, probably need it to decide on a threshold for sparsity -- and this is a simple way!
    # integer beta is easy!
    # sparsity which works both ways!
    matrix(rep((1:n_s)/n_s, n_out), nrow=n_out, byrow=TRUE),
    matrix(rep((1:n_out), n_s), nrow=n_out),
    matrix(rep((n_out:1), n_s), nrow=n_out)
  )
  
  beta <- beta/rowSums(beta)
  
  data %>%
    mutate(
      transfromed_Y1=beta %*% (Y1*W/sum(W) - Y1*(1-W)/sum(1-W)),
      transfromed_Y2=beta %*% (Y2*W/sum(W) - Y2*(1-W)/sum(1-W))
    )
}

cumulative_te <- function(data, order_var) { # weights! + naming + more vars. vs linear algebra?
  # min rep?
  
  data %>% arrange(desc({{ order_var }})) %>% mutate(
    'cumulative_Y1_{{ order_var }}' := cumsum(transfromed_Y1), # 
    'cumulative_Y2_{{ order_var }}' := cumsum(transfromed_Y2) # 
  )
}
############################# АЛЬТЕРНАТИВНЫЙ ВАРИАНТ CUM SUM ###############

StatCumulative <- ggproto("StatCumulative", Stat,
                          required_aes=c('order_var'),
                          compute_group = function(data, scales) {
                            data <- cumulative_te(data, order_var)
                            data$x <- data$cumulative_Y1
                            data$y <- data$cumulative_Y2
                            data.frame(data)
                          }
)


stat_cumulative <- function(mapping = NULL, data = NULL, geom = "line",
                            position = "identity", na.rm = FALSE, show.legend = NA, 
                            inherit.aes = TRUE, ...) {
  layer(
    stat = StatCumulative, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


# MAKE A SEPARATE PREIDICT/FIT PIPELINE
build_models <- function(data, group=c('treatment'), size=0.5, best_factor=0.5, variable_importance_th = 0) {
  # ADD KALLUS, ADD DUFLO/CHERNOZUKOV
  attach(data)
  train_ind <- as.integer(stratified(distinct(stratify_data), group=group, size=size)$cid)
  # train_ind <- as.integer(stratified(stratify_data, group=group, size=size)$cid)
  forest <- instrumental_forest(
    pretreatment_vars[cluster %in% train_ind,], 
    Y1[cluster %in% train_ind], 
    Y2[cluster %in% train_ind], 
    treatment[cluster %in% train_ind], 
    cluster=cluster[cluster %in% train_ind]
  )
  if (variable_importance_th > 0) {
    pretreatment_vars <- pretreatment_vars[,variable_importance(forest) > variable_importance_th]
    forest <- instrumental_forest(
      pretreatment_vars[cluster %in% train_ind,],
      Y1[cluster %in% train_ind], 
      Y2[cluster %in% train_ind], 
      treatment[cluster %in% train_ind],
      cluster=cluster[cluster %in% train_ind]
    )
    print(variable_importance(forest))
  }
  
  cforest_1 <- causal_forest(
    pretreatment_vars[cluster %in% train_ind,], 
    Y1[cluster %in% train_ind], 
    treatment[cluster %in% train_ind], 
    cluster=cluster[cluster %in% train_ind])#, tune.parameters = 'all') #complete???
  
  cforest_2 <- causal_forest(
    pretreatment_vars[cluster %in% train_ind,], 
    Y2[cluster %in% train_ind], 
    treatment[cluster %in% train_ind], 
    cluster=cluster[cluster %in% train_ind])#, tune.parameters = 'all')
  
  cforest_regr <- causal_forest(
    pretreatment_vars[cluster %in% train_ind,], 
    Y1[cluster %in% train_ind] - best_factor*Y2[cluster %in% train_ind], 
    treatment[cluster %in% train_ind], 
    cluster=cluster[cluster %in% train_ind])#, tune.parameters = 'all')
  
  predictions <- predict(forest, pretreatment_vars[!cluster %in% train_ind,])$predictions
  cpredictions_1 <- predict(cforest_1, pretreatment_vars[!cluster %in% train_ind,])$predictions # pass data
  cpredictions_2 <- predict(cforest_2, pretreatment_vars[!cluster %in% train_ind,])$predictions
  cpredictions_regr <- predict(cforest_regr, pretreatment_vars[!cluster %in% train_ind,])$predictions
  
  cpredictions <- cpredictions_1/cpredictions_2
  retval <- data.frame(
    predictions=predictions, 
    cpredictions_2=cpredictions_2, 
    cpredictions_regr=cpredictions_regr, 
    cpredictions_1=cpredictions_1,
    cpredictions=cpredictions,
    Y1=Y1[!cluster %in% train_ind],
    Y2=Y2[!cluster %in% train_ind],
    W=treatment[!cluster %in% train_ind],
    random=runif(length(Y1[!cluster %in% train_ind]),0,1)
  )
  detach(data)
  retval
}

###############################   DUFFLO DATA #######################################
prepare_Duflo_data <- function() {
  setwd('/Users/dianamirzoyan/Desktop/Diana/DISSERTATION/data/Duflo_/ETP_JPubE_Dataverse_Files/heterogeneous_marginal_costs')
  final_data <- read.csv(
    'Dufflo_exp/Dufflo_data-2.csv'
  )
  print(nrow(final_data))
  # final_data$costs <- final_data$treatment * 21.68 ### кенийский шиллиг
  final_data$costs <- final_data$treatment * 0.2 ### доллар
  stratify_data <- data.frame(treatment=final_data$treatment, cid=final_data$schoolid) # $1:length(
  pretreatment_vars <- final_data[
    c('girl', 'std_mark',	'realpercentile',	'kcpe2001',	'kcpe2004',	
      'total_2004',	'rotation',	'total1_2005',	'streams1_2005',	
      'init_clsize',	'Nteachers',	'gradesize'	,'schoolsize')
  ]
  
  list(
    pretreatment_vars=pretreatment_vars, 
    Y2=final_data$costs, 
    Y1=final_data$Y1, 
    cluster=final_data$schoolid, 
    treatment=final_data$treatment,
    stratify_data=stratify_data
  )
}
prepare_Duflo_data_upd <- function() {
  setwd('/Users/dianamirzoyan/Desktop/Diana/DISSERTATION/data/Duflo_/ETP_JPubE_Dataverse_Files/heterogeneous_marginal_costs')
  final_data <- read.csv(
    'Dufflo_data_new.csv'
  )
  final_data['Y2'] <- final_data$costs * 0.013765 ## перевели в доллар по курсу 2006 года
  final_data$Y2 <- final_data$Y2 / 1.3235 # уровень инфляции доллара с янв 2006 по дек 2020 года
  stratify_data <- data.frame(treatment=final_data$treatment, cid=final_data$schoolid)
  pretreatment_vars <- final_data[
    c('girl', 'std_mark',	'realpercentile',	'kcpe2001',	'kcpe2004',	
      'total_2004',	'rotation',	'total1_2005',	'streams1_2005',	
      'init_clsize',	'Nteachers',	'gradesize'	,'schoolsize')
  ]
  
  list(
    pretreatment_vars=pretreatment_vars, 
    Y2=final_data$Y2, 
    Y1=final_data$Y1, 
    cluster=final_data$schoolid, 
    treatment=final_data$treatment,
    stratify_data=stratify_data
  )
}

# prepare_Duflo_data_upd_0 <- function() {
#   setwd('/Users/dianamirzoyan/Desktop/Diana/DISSERTATION/data/Duflo_/ETP_JPubE_Dataverse_Files/heterogeneous_marginal_costs')
#   final_data <- read.csv(
#     'Dufflo_data_new_0.csv'
#   )
#   final_data['Y2'] <- final_data$costs * 0.013765 ## перевели в доллар по курсу 2006 года
#   final_data$Y2 <- final_data$Y2 / 1.3235 # уровень инфляции доллара с янв 2006 по дек 2020 года
#   stratify_data <- data.frame(treatment=final_data$treatment, cid=final_data$schoolid)
#   pretreatment_vars <- final_data[
#     c('girl', 'std_mark',	'realpercentile',	'kcpe2001',	'kcpe2004',	
#       'total_2004',	'rotation',	'total1_2005',	'streams1_2005',	
#       'init_clsize',	'Nteachers',	'gradesize'	,'schoolsize')
#   ]
#   
#   list(
#     pretreatment_vars=pretreatment_vars, 
#     Y2=final_data$Y2, 
#     Y1=final_data$Y1, 
#     cluster=final_data$schoolid, 
#     treatment=final_data$treatment,
#     stratify_data=stratify_data
#   )
# }
##################################################################################################################
set.seed(0)
school_data_Dufflo <- prepare_Duflo_data_upd() #build_models(school_data_Dufflo,  size=0.7, variable_importance_th = 0.03)
df_Dufflo <- build_models(school_data_Dufflo)
df_Dufflo_transf <- prepare_for_cum_te(df_Dufflo)
df_Dufflo_transf %>%
  cumulative_te(cpredictions_1) %>%
  cumulative_te(random) %>%
  cumulative_te(predictions) %>%
  cumulative_te(cpredictions) %>%
  cumulative_te(cpredictions_regr) %>%
  ggplot() + 
  geom_line(aes(y=cumulative_Y1_random, x=cumulative_Y2_random, color='Random')) +
  geom_line(aes(y=cumulative_Y1_predictions, x=cumulative_Y2_predictions, color='RATIO')) +
  # geom_line(mapping=aes(y=cumulative_Y1_cpredictions_1, x=cumulative_Y2_cpredictions_1, color='model O only')) +
  geom_line(mapping=aes(y=cumulative_Y1_cpredictions, x=cumulative_Y2_cpredictions, color='SEP')) +
  geom_line(mapping=aes(y=cumulative_Y1_cpredictions_regr, x=cumulative_Y2_cpredictions_regr, color='LIN_COMB')) +
  labs(x='Потраченные доллары на 1 ученика ',
       y='Эффект на учеников (рост успеваемости в баллах)'
  )

RMSE(df_Dufflo$Y1, df_Dufflo$random) ## 'Random'
RMSE(df_Dufflo$Y1, df_Dufflo$predictions) # 'RATIO'
RMSE(df_Dufflo$Y1, df_Dufflo$cpredictions) ## 'SEP'
RMSE(df_Dufflo$Y1, df_Dufflo$cpredictions_regr) ## 'LIN_COMB'

MAE(df_Dufflo$Y1, df_Dufflo$random) ## 'Random'
MAE(df_Dufflo$Y1, df_Dufflo$predictions) # 'RATIO'
MAE(df_Dufflo$Y1, df_Dufflo$cpredictions) ## 'SEP'
MAE(df_Dufflo$Y1, df_Dufflo$cpredictions_regr) ## 'LIN_COMB'

df_Dufflo['exp_name'] <- 'Dufflo'

df_Dufflo_cumsum <- df_Dufflo_transf %>%
  cumulative_te(cpredictions_1) %>%
  cumulative_te(random) %>%
  cumulative_te(predictions) %>%
  cumulative_te(cpredictions) %>%
  cumulative_te(cpredictions_regr)
df_Dufflo_cumsum['exp_name'] <- 'Dufflo'
# ############## не гонять
# set.seed(0)
# rem_data <- prepare_Ambler_data()
# rem_data$cluster <- rem_data$cluster * 100000
# rem_data$stratify_data$cid <- rem_data$stratify_data$cid * 100000
# df <- build_models(rem_data)
# 
# df %>%
#   cumulative_te(cpredictions_1) %>%
#   cumulative_te(random) %>%
#   cumulative_te(predictions) %>%
#   cumulative_te(cpredictions) %>%
#   cumulative_te(cpredictions_regr) %>%
#   ggplot() + geom_line(aes(y=cumulative_Y1_random, x=cumulative_Y2_random, color='Random')) +
#   geom_line(aes(y=cumulative_Y1_predictions, x=cumulative_Y2_predictions, color='RATIO')) +
#   # geom_line(mapping=aes(y=cumulative_Y1_cpredictions_1, x=cumulative_Y2_cpredictions_1, color='model O only')) +
#   geom_line(mapping=aes(y=cumulative_Y1_cpredictions, x=cumulative_Y2_cpredictions, color='SEP')) +
#   geom_line(mapping=aes(y=cumulative_Y1_cpredictions_regr, x=cumulative_Y2_cpredictions_regr, color='LIN_COMB')) +
#   labs(x='Процент оказываемого воздействия (участие в ETP)', y='Эффект на учеников (рост успеваемости в баллах)') +
#   ggtitle("Графический анализ гетерогенности эффективности - сравнение моделей")
# 
# getwd()
# ###################
# # treatment_tmp1 # ГОНЯТЬ
# ###################
set.seed(0)
rem_data <- prepare_Ambler_data_upd()

rem_data$cluster <- rem_data$cluster * 100000
rem_data$stratify_data$cid <- rem_data$stratify_data$cid * 100000
df <- build_models(rem_data)
df_transf <- prepare_for_cum_te(df)
df_transf %>%
  cumulative_te(cpredictions_1) %>%
  cumulative_te(random) %>%
  cumulative_te(predictions) %>%
  cumulative_te(cpredictions) %>%
  cumulative_te(cpredictions_regr) %>%
  ggplot() + geom_line(aes(y=cumulative_Y1_random, x=cumulative_Y2_random, color='Random')) +
  geom_line(aes(y=cumulative_Y1_predictions, x=cumulative_Y2_predictions, color='RATIO')) +
  # geom_line(mapping=aes(y=cumulative_Y1_cpredictions_1, x=cumulative_Y2_cpredictions_1, color='model O only')) +
  geom_line(mapping=aes(y=cumulative_Y1_cpredictions, x=cumulative_Y2_cpredictions, color='SEP')) +
  geom_line(mapping=aes(y=cumulative_Y1_cpredictions_regr, x=cumulative_Y2_cpredictions_regr, color='LIN_COMB')) +
  labs(x='Потраченные доллары на 1 ученика - мигранта',
       y='Эффект на учеников (рост денежных трансферов на образование)')
# 
# ###################
# # treatment_tmp2 #
# ###################
# prepare_Ambler_data_upd <- function() {
#   setwd('/Users/dianamirzoyan/Desktop/Diana/DISSERTATION/data/Ambler/Remittances Dataverse Files/data')
#   final_data <- read.csv(
#     'Ambler_data_upd.csv'
#   )
#   print(length(final_data$id))
#   final_data <- na.omit(final_data)
#   print(length(final_data$id))
#   stratify_data <- data.frame(treatment=final_data$treatment_tmp2, cid=final_data$dayXlocation) # $1:length(
#   print(length(stratify_data$cid))
#   pretreatment_vars <- final_data[
#     c('us_stu_insch', 'us_stu_edlvl', 'us_stu_grade', 'us_stu_edlvl_all',
#       'us_stu_female', 'us_stu_mig_rel', 'us_stu_age', 'us_stu_grade_all',
#       'us_stu_yearsofed', 'us_stu_mig_rel_child', 'us_stu_mig_rel_sib',
#       'us_stu_mig_rel_niece', 'us_stu_mig_rel_cousin')
#   ]
# 
#   list(
#     pretreatment_vars=pretreatment_vars,
#     Y2=final_data$costs,
#     Y1=final_data$er_amount_total,
#     cluster=final_data$dayXlocation,
#     treatment=final_data$treatment_tmp2,
#     stratify_data=stratify_data
#   )
# }
# set.seed(0)
# rem_data <- prepare_Ambler_data_upd()
# # rem_data <- na.omit(rem_data)
# rem_data$cluster <- as.integer(rem_data$cluster * 100000)
# rem_data$stratify_data$cid <- as.integer(rem_data$stratify_data$cid * 100000)
# df <- build_models(rem_data)
# 
# df %>%
#   cumulative_te(cpredictions_1) %>%
#   cumulative_te(random) %>%
#   cumulative_te(predictions) %>%
#   cumulative_te(cpredictions) %>%
#   cumulative_te(cpredictions_regr) %>%
#   ggplot() + geom_line(aes(y=cumulative_Y1_random, x=cumulative_Y2_random, color='Random')) +
#   geom_line(aes(y=cumulative_Y1_predictions, x=cumulative_Y2_predictions, color='RATIO')) +
#   # geom_line(mapping=aes(y=cumulative_Y1_cpredictions_1, x=cumulative_Y2_cpredictions_1, color='model O only')) +
#   geom_line(mapping=aes(y=cumulative_Y1_cpredictions, x=cumulative_Y2_cpredictions, color='SEP')) +
#   geom_line(mapping=aes(y=cumulative_Y1_cpredictions_regr, x=cumulative_Y2_cpredictions_regr, color='LIN_COMB')) +
#   labs(x='Процент оказываемого воздействия (участие в ETP)', y='Эффект на учеников (рост успеваемости в баллах)') +
#   ggtitle("Графический анализ гетерогенности эффективности - сравнение моделей")
# 


########################### Giligan Data ##################
############## не гонять

prepare_Giligan_data_upd <- function() {
  setwd('/Users/dianamirzoyan/Desktop/Diana/DISSERTATION/data/Gilligan')
  pre_final_data <- read.csv(
    'Gilligan_data.csv'
  )
  final_data <- na.omit(pre_final_data)
  stratify_data <- data.frame(
    treatment=final_data$treatment, cid=final_data$a01
  ) 
  pretreatment_vars <- final_data[
    c('ple_tried', 'ple_pass','pre_theta_sd',
      'books_with', 
      'b_homework_pp', 'b_fem_pp', 'b_elec_pp', 'b_cement_pp', 'b_radio_pp',
      'b_edhigh_teacher', 'b_fem_teacher', 
      'decile', 'scores_top'
    )
  ]
  
  list(
    pretreatment_vars=pretreatment_vars, 
    Y2=final_data$decile_decile, 
    Y1=final_data$ple_pass, 
    cluster=final_data$a01, 
    treatment=final_data$treatment,
    stratify_data=stratify_data
  )
}

prepare_Giligan_data_upd <- function() {
  setwd('/Users/dianamirzoyan/Desktop/Diana/DISSERTATION/data/Gilligan')
  pre_final_data <- read.csv(
    'Gilligan_data.csv'
  )
  final_data <- na.omit(pre_final_data)
  # final_data$decile_decile <- final_data$decile_decile * 0.000276 ### переводим в доллары
  final_data$decile_decile <- final_data$decile_decile * 0.000291 ## перевели в доллар по курсу 2016 года
  final_data$decile_decile <- final_data$decile_decile / 1.0938 # уровень инфляции доллара с апр 2016 по дек 2020 года
  stratify_data <- data.frame(
    treatment=final_data$treatment, cid=final_data$a01
  ) 
  pretreatment_vars <- final_data[
    c(
      'books_with', 
      'decile',
      'b_homework_pp', 'b_fem_pp',
      'b_edhigh_teacher'#'b_fem_teacher'
    )
  ]
  
  list(
    pretreatment_vars=pretreatment_vars, 
    Y2=final_data$decile_decile, 
    Y1=final_data$ple_tried, #endline_score
    cluster=final_data$a01, 
    treatment=final_data$treatment,
    stratify_data=stratify_data
  )
}

set.seed(0)

school_data <- prepare_Giligan_data_upd()
df <- build_models(school_data)
# df_transf <- prepare_for_cum_te(df)
df %>%
  cumulative_te(cpredictions_1) %>%
  cumulative_te(random) %>%
  cumulative_te(predictions) %>%
  cumulative_te(cpredictions) %>%
  cumulative_te(cpredictions_regr) %>%
  ggplot() + geom_line(aes(y=cumulative_Y1_random, x=cumulative_Y2_random, color='Random')) +
  geom_line(aes(y=cumulative_Y1_predictions, x=cumulative_Y2_predictions, color='RATIO')) +
  # geom_line(mapping=aes(y=cumulative_Y1_cpredictions_1, x=cumulative_Y2_cpredictions_1, color='model O only')) +
  geom_line(mapping=aes(y=cumulative_Y1_cpredictions, x=cumulative_Y2_cpredictions, color='SEP')) +
  geom_line(mapping=aes(y=cumulative_Y1_cpredictions_regr, x=cumulative_Y2_cpredictions_regr, color='LIN_COMB')) +
  labs(x='Потраченные доллары на 1 ученика  ',
       y='Эффект на учеников (количество попыток сдать экзамен)')

random_RMSE <- c()
random_MAE <- c()
RATIO_RMSE <- c()
RATIO_MAE <- c()
SEP_RMSE <- c()
SEP_MAE <- c()
LIN_COMB_RMSE <- c()
LIN_COMB_MAE <- c()
for (i in 0:100){
  set.seed(i)
  print(i)
  school_data_Dufflo <- prepare_Giligan_data_upd()
  df <- build_models(school_data_Dufflo)
  
  random_RMSE[i+1] <- RMSE(df$Y1, df$random) ## 'Random'
  RATIO_RMSE[i+1] <- RMSE(df$Y1, df$predictions) # 'RATIO'
  SEP_RMSE[i+1] <- RMSE(df$Y1, df$cpredictions) ## 'SEP'
  LIN_COMB_RMSE[i+1] <- RMSE(df$Y1, df$cpredictions_regr) ## 'LIN_COMB'
  
  random_MAE[i+1] <- MAE(df$Y1, df$random) ## 'Random'
  RATIO_MAE[i+1] <- MAE(df$Y1, df$predictions) # 'RATIO'
  SEP_MAE[i+1] <- MAE(df$Y1, df$cpredictions) ## 'SEP'
  LIN_COMB_MAE[i+1] <- MAE(df$Y1, df$cpredictions_regr) ## 'LIN_COMB'
}

df_Giligan_cumsum <- df %>%
  cumulative_te(cpredictions_1) %>%
  cumulative_te(random) %>%
  cumulative_te(predictions) %>%
  cumulative_te(cpredictions) %>%
  cumulative_te(cpredictions_regr)

df_Giligan_cumsum['exp_name'] <- 'Giligan'

df_cumsum <- rbind(df_Dufflo_cumsum,df_Giligan_cumsum)

min(df_Dufflo_cumsum$cumulative_Y2_predictions)
max(df_Dufflo_cumsum$cumulative_Y2_predictions)

min(df_Giligan_cumsum$cumulative_Y2_predictions)
max(df_Giligan_cumsum$cumulative_Y2_predictions)

df_cumsum_sorted <- as.data.frame(
  sort(decreasing = TRUE, x = df_cumsum$cumulative_Y2_predictions)
)
colnames(df_cumsum_sorted) <- c('predicted_Y2')
df_cumsum_sorted['priority'] <- seq(from=1, to=nrow(df_cumsum_sorted), by=1)

df_tmp <- inner_join(df_cumsum, df_cumsum_sorted, by=c('cumulative_Y2_predictions' = 'predicted_Y2'))

RMSE(df$Y1, df$random)
MAE(df$Y1, df$cpredictions)

RMSE(df$Y1, df$random)/sqrt(var(df$Y1))
1 - sqrt(var(df$Y1))/RMSE(df$Y1, df$random)





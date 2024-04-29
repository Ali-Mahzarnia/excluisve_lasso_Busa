library(xlsx)
library(ggplot2)
library("readxl")
library("emmeans")
options(java.parameters = "-Xmx1000m")


path = "master_df_all_all_all_all_all_combined.xlsx"
df_original = read_excel(path)


df = df_original[as.numeric( df_original$risk_for_ad )<2 ,]

mean_names = colnames(df)[grepl("meanfa", colnames(df))]
sd_names = colnames(df)[grepl("sdfa", colnames(df))]
num_st_names = colnames(df)[  grepl("num_sl", colnames(df)) & ! grepl("num_sl_assym", colnames(df)) ]
vol_st_names = colnames(df)[  grepl("vol_sl", colnames(df)) & ! grepl("vol_sl_assym", colnames(df)) ]
len_st_names = colnames(df)[  grepl("len_sl", colnames(df)) & ! grepl("len_sl_assym", colnames(df)) ]
BUAN_st_names = colnames(df)[  grepl("BUAN", colnames(df))  ]
AI_st_names = colnames(df)[  grepl("assym", colnames(df)) & !grepl("sl", colnames(df))  ]


age = df["age"]


df_run_x =df_run_x = cbind(  df[sd_names], df[num_st_names], df[vol_st_names], df[len_st_names], df[BUAN_st_names], df[AI_st_names] ) #cbind( df[mean_names], df[sd_names], df[num_st_names], df[vol_st_names], df[len_st_names] )
all_names = c(sd_names, num_st_names, vol_st_names, len_st_names,BUAN_st_names , AI_st_names )
df_run_x = t(apply(df_run_x, 1, unlist))
df_run_x = t(apply(df_run_x, 1, as.numeric))
df_run_x = scale(df_run_x)
df_run_x[is.na(df_run_x)]  = 0 
df_run_y = t(age)
df_run_y = apply(df_run_y, 1, unlist)
df_run_y = apply(df_run_y, 1, as.numeric)
df_run = cbind(df_run_y,df_run_x)

library(ExclusiveLasso)
set.seed(123)
v.group = c( rep(1,length(sd_names)), rep(2,length(num_st_names)), rep(3,length(vol_st_names)), rep(4,length(len_st_names)), rep(5,length(BUAN_st_names)), rep(6,length(AI_st_names))  )

ex_cv <- cv.exclusive_lasso( df_run_x,df_run_y, groups = v.group,family="gaussian", 
                             intercept = T, nfolds=10)
plot(ex_cv)

ex_cv$lambda.min
ex_cv$lambda.1se


sqrt(ex_cv$cvm[ex_cv$lambda == ex_cv$lambda.min])
sd(df_run_y)


ex <- exclusive_lasso( df_run_x,df_run_y, groups = v.group,family="gaussian", 
                       intercept = T,  lambda = ex_cv$lambda.min )



coefs = ex$coef
index= which(coefs != 0)



results =cbind(all_names[index], coefs[index] )
write.xlsx2( results , paste0(length(mean_names),"_bundles_full_model_buan_ai.xlsx") )


# x=df_run_x[1,]
# ex$intercept + t(coefs)%*%x 


age_gap_whole = as.data.frame(predict(ex) - df_run_y   )
colnames(age_gap_whole) = "age_gap"
age_gap_whole$geno =  df$geno
age_gap_whole$age = df$age  

# mod = lm(age_gap ~ as.numeric(age)*geno, data = age_gap_whole)
# emm = emmeans(mod, "geno", at = list(df_run_y = 0))
# pval = pairs(emm)

ggplot(age_gap_whole, aes( x=as.numeric(age), y= age_gap, fill = geno, color = geno  ) )+
  geom_smooth( method='lm')#+geom_text(x=40, y=15, label=paste0("Intercept P Values" , ) )

ggsave(  paste0(length(mean_names),"_bundles_age_gap_vs_age_buan_ai.png"), plot= last_plot())


############## now bias correction

index3 = df$geno=="APOE33"

y_3_pred =  predict(ex, df_run_x[index3,])
y_3_true = df_run_y[index3]


df_sick = df_original[as.numeric( df_original$risk_for_ad )>1 ,]

df_sick_x = cbind(  df_sick[sd_names], df_sick[num_st_names], df_sick[vol_st_names], df_sick[len_st_names], df_sick[BUAN_st_names],  df_sick[AI_st_names]  )
df_sick_x = t(apply(df_sick_x, 1, unlist))
df_sick_x = t(apply(df_sick_x, 1, as.numeric))
df_sick_x = scale(df_sick_x)
df_sick_x[is.na(df_sick_x)]  = 0 
age_sick = df_sick["age"]
df_risk_y = t(age_sick)
df_risk_y = apply(df_risk_y, 1, unlist)
df_risk_y = apply(df_risk_y, 1, as.numeric)
df_risk_all = cbind(df_risk_y,df_sick_x)


y_scik_pred = predict(ex, df_sick_x)
y_scik_true = df_risk_y




###### alpha beta
lm = lm(y_3_pred~y_3_true)
alpha = lm$coefficients[1]
beta = lm$coefficients[2]

age_predicted_whole = as.data.frame( predict(ex))
age_predicted_whole = rbind(age_predicted_whole,y_scik_pred )
colnames(age_predicted_whole) = "age_pred"
age_predicted_whole$geno =  c( df$geno, df_sick$geno)
age_predicted_whole$age = c ( df$age , df_risk_y)

age_predicted_whole[,1] = (as.numeric(age_predicted_whole[,1])-alpha)/beta

age_predicted_whole$age_gap_corr= age_predicted_whole$age_pred - as.numeric(age_predicted_whole$age)
age_predicted_whole$risk = c(df$risk_for_ad, df_sick$risk_for_ad)

age_predicted_whole_plot = age_predicted_whole[age_predicted_whole$risk<2,]
age_gap_sick_df = age_predicted_whole[age_predicted_whole$risk>1, ]

ggplot(age_predicted_whole_plot, aes( x=as.numeric(age), y= age_gap_corr, fill = geno, color = geno  ) )+
  geom_point()+geom_smooth( method='lm') + #+geom_text(x=40, y=15, label=paste0("Intercept P Values" , ) )
  geom_point(data=age_gap_sick_df,  aes(x=as.numeric(age),y=age_gap_corr ),colour="green")+theme_bw()




ggsave(  paste0(length(mean_names),"_bundles_age_gap_vs_age_after_correction_buan_ai.png"), plot= last_plot())


ggplot(age_predicted_whole_plot, aes( x=as.numeric(age), y= age_pred, fill = geno, color = geno  ) )+
  geom_point()+geom_smooth( method='lm') + #+geom_text(x=40, y=15, label=paste0("Intercept P Values" , ) )
  geom_point(data=age_gap_sick_df,  aes(x=as.numeric(age),y=age_pred ),colour="green") + theme_bw()

ggsave(  paste0(length(mean_names),"_bundles_age_pred_vs_age_after_correction_buan_ai.png"), plot= last_plot())



# 
# ##############################
# ###
# ### only apoe3 
# ##############################
# 
# df3 = df_original[df_original$geno == "APOE33",] 
# df = df3
# 
# 
# mean_names = colnames(df)[grepl("meanfa", colnames(df))]
# sd_names = colnames(df)[grepl("sdfa", colnames(df))]
# num_st_names = colnames(df)[  grepl("num_sl", colnames(df)) & ! grepl("num_sl_assym", colnames(df)) ]
# vol_st_names = colnames(df)[  grepl("vol_sl", colnames(df)) & ! grepl("vol_sl_assym", colnames(df)) ]
# len_st_names = colnames(df)[  grepl("len_sl", colnames(df)) & ! grepl("len_sl_assym", colnames(df)) ]
# 
# age = df["age"]
# 
# 
# df_run_x = cbind( df[mean_names], df[sd_names], df[num_st_names], df[vol_st_names], df[len_st_names] )
# all_names = c(mean_names,sd_names, num_st_names, vol_st_names, len_st_names )
# df_run_x = t(apply(df_run_x, 1, unlist))
# df_run_x = t(apply(df_run_x, 1, as.numeric))
# df_run_x = scale(df_run_x)
# # df_run_x[is.na(df_run_x)]  = 0 
# df_run_y = t(age) 
# df_run_y = apply(df_run_y, 1, unlist)
# df_run_y = apply(df_run_y, 1, as.numeric)
# df_run = cbind(df_run_y,df_run_x)
# 
# library(ExclusiveLasso)
# v.group = c( rep(1,length(mean_names)) , rep(2,length(sd_names)), rep(3,length(num_st_names)), rep(4,length(vol_st_names)), rep(5,length(len_st_names))  )
# 
# ex_cv <- cv.exclusive_lasso( df_run_x,df_run_y, groups = v.group,family="gaussian", 
#                              intercept = T, nfolds=10)
# plot(ex_cv)
# 
# ex_cv$lambda.min
# ex_cv$lambda.1se
# 
# 
# sqrt(ex_cv$cvm[ex_cv$lambda == ex_cv$lambda.min])
# sd(df_run_y)
# 
# 
# ex <- exclusive_lasso( df_run_x,df_run_y, groups = v.group,family="gaussian", 
#                        intercept = T,  lambda = ex_cv$lambda.min )
# 
# 
# 
# coefs = ex$coef
# index= which(coefs != 0)
# 
# 
# 
# 
# results =cbind(all_names[index], coefs[index] )
# write.xlsx2( results , paste0(length(mean_names),"_bundles_apoe3_model.xlsx") )
# 
# 
# age_gap_3 = df_run_y - predict(ex)
# 
# 
# #########
# ########
# ## only apoe 4 
# #######
# 
# 
# 
# df4 = df_original[df_original$geno == "APOE44",] 
# df = df4
# 
# 
# mean_names = colnames(df)[grepl("meanfa", colnames(df))]
# sd_names = colnames(df)[grepl("sdfa", colnames(df))]
# num_st_names = colnames(df)[  grepl("num_sl", colnames(df)) & ! grepl("num_sl_assym", colnames(df)) ]
# vol_st_names = colnames(df)[  grepl("vol_sl", colnames(df)) & ! grepl("vol_sl_assym", colnames(df)) ]
# len_st_names = colnames(df)[  grepl("len_sl", colnames(df)) & ! grepl("len_sl_assym", colnames(df)) ]
# 
# age = df["age"]
# 
# 
# df_run_x = cbind( df[mean_names], df[sd_names], df[num_st_names], df[vol_st_names], df[len_st_names] )
# all_names = c(mean_names,sd_names, num_st_names, vol_st_names, len_st_names )
# df_run_x = t(apply(df_run_x, 1, unlist))
# df_run_x = t(apply(df_run_x, 1, as.numeric))
# df_run_x = scale(df_run_x)
# # df_run_x[is.na(df_run_x)]  = 0 
# df_run_y = t(age)
# df_run_y = apply(df_run_y, 1, unlist)
# df_run_y = apply(df_run_y, 1, as.numeric)
# df_run = cbind(df_run_y,df_run_x)
# 
# library(ExclusiveLasso)
# v.group = c( rep(1,length(mean_names)) , rep(2,length(sd_names)), rep(3,length(num_st_names)), rep(4,length(vol_st_names)), rep(5,length(len_st_names))  )
# 
# ex_cv <- cv.exclusive_lasso( df_run_x,df_run_y, groups = v.group,family="gaussian", 
#                              intercept = T, nfolds=10)
# plot(ex_cv)
# 
# ex_cv$lambda.min
# ex_cv$lambda.1se
# 
# 
# sqrt(ex_cv$cvm[ex_cv$lambda == ex_cv$lambda.min])
# sd(df_run_y)
# 
# 
# ex <- exclusive_lasso( df_run_x,df_run_y, groups = v.group,family="gaussian", 
#                        intercept = T,  lambda = ex_cv$lambda.min )
# 
# 
# 
# coefs = ex$coef
# index= which(coefs != 0)
# 
# 
# 
# results =cbind(all_names[index], coefs[index] )
# write.xlsx2( results , paste0(length(mean_names),"_bundles_apoe4_model.xlsx"))
# 
# 
# age_gap_4 = df_run_y - predict(ex)
# 
# age_gap = rbind(age_gap_3, age_gap_4)
# age_gap = as.numeric(age_gap)
# age_gap = cbind (age_gap, c( rep("APOE33", length(age_gap_3) ) , rep("APOE44", length(age_gap_4) )     ))
# colnames(age_gap) = c("Gap", "APOE" )
# 
# age_gap = as.data.frame(age_gap)
# 
# library(ggplot2)
# dodge <- position_dodge(width = 0.5 )
# mycolors <- c('blueviolet', 'chartreuse1', 'red', 'azure3')
# 
# ggplot(age_gap, aes(x =APOE, y=as.numeric(Gap) , fill = APOE )  ) + 
#   geom_violin(inherit.aes=TRUE,position=dodge) +
#   scale_color_manual(values=mycolors)+
#   scale_fill_manual(values=mycolors)+
#   geom_boxplot(color="black", outlier.color="black", width=0.4, alpha=.6, position=dodge) +
#   geom_dotplot(binaxis='y', stackdir='center', stackratio = 0.7, dotsize=1, alpha=0.6, position=dodge)+
#   labs(title = "Age Gap")+
#   theme_minimal()+
#   theme_bw()+
#   labs(x = "Genoptype", y = paste0( "Age Gap") ) +
#   theme_bw() 
# ggsave( paste0(length(mean_names),"_bundles_Age_gap.png" ), plot = last_plot(), dpi = 300)
# 
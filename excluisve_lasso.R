library(xlsx)
path = "master_df_all_all_all_all.xlsx"
df_original = read.xlsx2(path, sheetName = "Sheet1")


df = df_original

mean_names = colnames(df)[grepl("meanfa", colnames(df))]
sd_names = colnames(df)[grepl("sdfa", colnames(df))]
num_st_names = colnames(df)[  grepl("num_sl", colnames(df)) & ! grepl("num_sl_assym", colnames(df)) ]
vol_st_names = colnames(df)[  grepl("vol_sl", colnames(df)) & ! grepl("vol_sl_assym", colnames(df)) ]
len_st_names = colnames(df)[  grepl("len_sl", colnames(df)) & ! grepl("len_sl_assym", colnames(df)) ]

age = df["age"]


df_run_x = cbind( df[mean_names], df[sd_names], df[num_st_names], df[vol_st_names], df[len_st_names] )
all_names = c(mean_names,sd_names, num_st_names, vol_st_names, len_st_names )
df_run_x = t(apply(df_run_x, 1, unlist))
df_run_x = t(apply(df_run_x, 1, as.numeric))
df_run_x = scale(df_run_x)
# df_run_x[is.na(df_run_x)]  = 0 
df_run_y = t(age)
df_run_y = apply(df_run_y, 1, unlist)
df_run_y = apply(df_run_y, 1, as.numeric)
df_run = cbind(df_run_y,df_run_x)

library(ExclusiveLasso)
set.seed(123)
v.group = c( rep(1,length(mean_names)) , rep(2,length(sd_names)), rep(3,length(num_st_names)), rep(4,length(vol_st_names)), rep(5,length(len_st_names))  )

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
write.xlsx2( results , paste0(length(mean_names),"_bundles_full_model.xlsx") )
##############################
###
### only apoe3 
##############################

df3 = df_original[df_original$geno == "APOE33",] 
df = df3


mean_names = colnames(df)[grepl("meanfa", colnames(df))]
sd_names = colnames(df)[grepl("sdfa", colnames(df))]
num_st_names = colnames(df)[  grepl("num_sl", colnames(df)) & ! grepl("num_sl_assym", colnames(df)) ]
vol_st_names = colnames(df)[  grepl("vol_sl", colnames(df)) & ! grepl("vol_sl_assym", colnames(df)) ]
len_st_names = colnames(df)[  grepl("len_sl", colnames(df)) & ! grepl("len_sl_assym", colnames(df)) ]

age = df["age"]


df_run_x = cbind( df[mean_names], df[sd_names], df[num_st_names], df[vol_st_names], df[len_st_names] )
all_names = c(mean_names,sd_names, num_st_names, vol_st_names, len_st_names )
df_run_x = t(apply(df_run_x, 1, unlist))
df_run_x = t(apply(df_run_x, 1, as.numeric))
df_run_x = scale(df_run_x)
# df_run_x[is.na(df_run_x)]  = 0 
df_run_y = t(age) 
df_run_y = apply(df_run_y, 1, unlist)
df_run_y = apply(df_run_y, 1, as.numeric)
df_run = cbind(df_run_y,df_run_x)

library(ExclusiveLasso)
v.group = c( rep(1,length(mean_names)) , rep(2,length(sd_names)), rep(3,length(num_st_names)), rep(4,length(vol_st_names)), rep(5,length(len_st_names))  )

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
write.xlsx2( results , paste0(length(mean_names),"_bundles_apoe3_model.xlsx") )


age_gap_3 = df_run_y - predict(ex)


#########
########
## only apoe 4 
#######



df4 = df_original[df_original$geno == "APOE44",] 
df = df4


mean_names = colnames(df)[grepl("meanfa", colnames(df))]
sd_names = colnames(df)[grepl("sdfa", colnames(df))]
num_st_names = colnames(df)[  grepl("num_sl", colnames(df)) & ! grepl("num_sl_assym", colnames(df)) ]
vol_st_names = colnames(df)[  grepl("vol_sl", colnames(df)) & ! grepl("vol_sl_assym", colnames(df)) ]
len_st_names = colnames(df)[  grepl("len_sl", colnames(df)) & ! grepl("len_sl_assym", colnames(df)) ]

age = df["age"]


df_run_x = cbind( df[mean_names], df[sd_names], df[num_st_names], df[vol_st_names], df[len_st_names] )
all_names = c(mean_names,sd_names, num_st_names, vol_st_names, len_st_names )
df_run_x = t(apply(df_run_x, 1, unlist))
df_run_x = t(apply(df_run_x, 1, as.numeric))
df_run_x = scale(df_run_x)
# df_run_x[is.na(df_run_x)]  = 0 

df_run_y = t(age)
df_run_y = apply(df_run_y, 1, unlist)
df_run_y = apply(df_run_y, 1, as.numeric)
df_run = cbind(df_run_y,df_run_x)

library(ExclusiveLasso)
v.group = c( rep(1,length(mean_names)) , rep(2,length(sd_names)), rep(3,length(num_st_names)), rep(4,length(vol_st_names)), rep(5,length(len_st_names))  )

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
write.xlsx2( results , paste0(length(mean_names),"_bundles_apoe4_model.xlsx"))


age_gap_4 = df_run_y - predict(ex)

age_gap = rbind(age_gap_3, age_gap_4)
age_gap = as.numeric(age_gap)
age_gap = cbind (age_gap, c( rep("APOE33", length(age_gap_3) ) , rep("APOE44", length(age_gap_4) )     ))
colnames(age_gap) = c("Gap", "APOE" )

age_gap = as.data.frame(age_gap)

library(ggplot2)
dodge <- position_dodge(width = 0.5 )
mycolors <- c('blueviolet', 'chartreuse1', 'red', 'azure3')

ggplot(age_gap, aes(x =APOE, y=as.numeric(Gap) , fill = APOE )  ) + 
  geom_violin(inherit.aes=TRUE,position=dodge) +
  scale_color_manual(values=mycolors)+
  scale_fill_manual(values=mycolors)+
  geom_boxplot(color="black", outlier.color="black", width=0.4, alpha=.6, position=dodge) +
  geom_dotplot(binaxis='y', stackdir='center', stackratio = 0.7, dotsize=1, alpha=0.6, position=dodge)+
  labs(title = "Age Gap")+
  theme_minimal()+
  theme_bw()+
  labs(x = "Genoptype", y = paste0( "Age Gap") ) +
  theme_bw() 
ggsave( paste0(length(mean_names),"_bundles_Age_gap.png" ), plot = last_plot(), dpi = 300)


#setting "event" variables
df %>% 
  mutate(event = case_when(non_informative_censoring == 1 ~ 0,
                           stability == 1 ~ 1,
                           death == 1 ~ 2,
                           (stability == 0 & death == 0) ~ 3)) -> df
#Multiple imputation with 100 imputed datasets
m <- 100 
imp <- mice(df, m = m ,seed=12345 ,maxit=50 ,printFlag=FALSE)
#Calculating propensity score within each imputed dataset
imp_stack <- complete(imp,action="long") %>% 
  as_data_frame
imp_stack<-imp_stack %>%
  group_by(.imp)%>%
  nest() %>%
  mutate(ps=map(data,function(df){
    ps_model<-glm(steroid ~ age + gender + adl + hot + bun + rr + hospital, family=binomial,data=df)
    ps<-predict(ps_model,type="response")
    return(ps)
  }))%>%
  unnest()
#Propensity score matching
imp_stack<-imp_stack %>%
  group_by(.imp)%>%
  nest() %>%
  mutate(match.weight=map(data,function(df){
    m.out<-matchit(steroid ~ ps, data = df,replace=FALSE,caliper=0.2)
    match.weight <- m.out$weights
  }))%>%
  unnest()  
imp_matched<-imp_stack %>% filter(match.weight==1)
table(imp_matched$steroid)
#Histogram describing propensity score in each group
histogram_imp_matched <- 
  ggplot(data = imp_matched, aes(x = ps)) +
  geom_histogram() +
  facet_grid(. ~ steroid)
histogram_imp_matched
#Balance check 
tab_imp_matched <- CreateTableOne(vars=c("age","gender", "adl", "hot", "bun", "rr", "ams", "hr", "hospital"),strata="steroid", data=imp_matched, test=FALSE)
tab_imp_matched %>% 
  print(smd=TRUE) %>% 
  write.csv(file = "table.csv")
imp_matched %>% 
  group_by(.imp) %>% 
  nest() %>% 
  mutate(smd=map(data,function(df){
    smd <- ExtractSmd(CreateTableOne(vars=c("age","gender", "adl", "hot", "bun", "rr", "ams", "hr", "hospital"),strata="steroid", data=df, test=FALSE))
    df <- as_data_frame(smd) %>% 
      mutate(name=rownames(smd))
    names(df)<-c("smd","name")
    df
  })) %>% 
  select(.imp, smd) %>% 
  unnest() %>% 
  ggplot(mapping=aes(x=name,y=smd,group=.imp))+
  geom_line()+
  geom_hline(yintercept=0.1)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,vjust=0.5),
        legend.key=element_blank(),
        plot.title=element_text(hjust=0.5),
        strip.background=element_blank())
#Outcome analysis 
#Fine and Gray model within each dataset
imp_analysis <- imp_matched %>%
  group_by(.imp)%>%
  nest() 
model <- function(df){
  CI.multi <- model.matrix(object = ~ steroid + adl + hot + age + gender + bun + rr + ams + hr + hospital, data=df)
  CI.multi <- CI.multi[,-1]
  res_cimulti <- crr(ftime = df$hospitalization, fstatus = df$event, cov1 = CI.multi, failcode = 1, cencode = 0, variance=TRUE)
}
#Rubin's combining rule
list_of_models <- map(imp_analysis$data, model)
list_of_coefs <- list()
for(i in unique(imp_analysis$.imp)){
  list_of_coefs[[i]] <- list_of_models[[i]]$coef 
} 
list_of_vcovs <- list()
for(i in unique(imp_analysis$.imp)){
  list_of_vcovs[[i]] <- list_of_models[[i]]$var 
} 
summary(MIcombine(list_of_coefs, list_of_vcovs)) -> results
exp(results[1,1:4])

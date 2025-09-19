# Load packages 
  pacman::p_load(tidyverse,data.table,readxl,ggplot2,writexl,RColorBrewer,lubridate,tableone)
  # Load data
  dx <- readRDS("Diagnosis.RDS") # Diagnosis data 
  rx <- readRDS("Prescription_20002021_20092023.RDS") # Prescription data
  # Format prescription date and diagnosis date
  rx[,DispensingDate:=as.Date(`DispensingDate(yyyy-mm-dd)`)]
  rx[,Disenddate:=DispensingDate+as.numeric(DispensingDuration)]
  setorder(rx,DispensingDate)
  dx$date <- as_date(dx$date)
  setorder(dx,date)
  setnames(rx,"TherapeuticClassification(BNF,Principal)","BNFcode")  
  #Remove duplicate records
  dx <- distinct(dx)
  rx <- distinct(rx)
  
# 1 Generate drug user cohort --------
  # We have prescription records from 2009 to 2023, thus we use 2009 as reference year, any patients with 
  # prescription of biological and target synthetic disease modifying anti-rheumatic drugs (b/ts DMARDs)
  # in 2009 were removed to ascertain first time exposure to b/ts DMARDs
  
  druglist <- c("INFLIXIMAB","ADALIMUMAB","CERTOLIZUMAB","ETANERCEPT",
                "GOLIMUMAB","ABATACEPT","RITUXIMAB","SARILUMAB","TOCILIZUMAB",
                "BARICITINIB","TOFACITINIB","UPADACITINIB") # List of b/ts DMARDs
  alldrug<- rx[str_detect(DrugName,regex("Infliximab|Adalimumab|Certolizumab|Etanercept|Golimumab|Abatacept|Rituximab|Sarilumab|Tocilizumab|Baricitinib|Tofacitinib|Upadacitinib",ignore_case = T))]
  for (i in 1:length(druglist)) {
    alldrug[str_detect(DrugName,regex(druglist[i],ignore_case = T)),agent:=toupper(druglist[i])]
    print(i)} #Standardize drug name to upper case
  # alldrug[DispensingDate==Disenddate,Disenddate:=]
  preuser <- alldrug[DispensingDate<as.Date("2010-01-01")]#RA patients with prescription of b/ts DMARDs in 2009
  alldrug<- alldrug[!ReferenceKey%in%preuser$ReferenceKey]#Exclude these previous user's prescription information
  alldrug[,indexdate:=min(DispensingDate),ReferenceKey]#Identify the earliest date of using b/ts DMARDs
  initre <- alldrug[DispensingDate==indexdate]
  initre <-distinct(initre,indexdate,agent,ReferenceKey)#Identify the first b/ts DMARDs products patients used
  alldrug <- left_join(alldrug,initre[,.(indexdate,ReferenceKey,initre=agent)])
  alldrug <- distinct(alldrug,ReferenceKey,DispensingDate,Disenddate,agent,.keep_all=T)
  single <- alldrug[,.N,.(ReferenceKey,agent)][N==1]
  
  #Tidy up the prescription period of b/ts DMARDs, follow up until 1) Treatment discontinuation  
  #2)switch to other b/ts DMARDs 3)prescription gap greater than 180 days
  cc <- alldrug %>%
    arrange(ReferenceKey, DispensingDate) %>%
    mutate(DispensingDate = ymd(DispensingDate), 
           Disenddate = ymd(Disenddate)) %>%
    group_by(ReferenceKey) %>%
    mutate(indx = c(0, cumsum(as.numeric(lead(DispensingDate)) > 
                                cummax(as.numeric(Disenddate)))[-n()])) %>%
    group_by(ReferenceKey, indx, agent) %>%
    summarise(start = first(DispensingDate), end = last(Disenddate))
  
  calc_time <- function(refDate, endDate) {
    period <- as.period(interval(refDate, endDate), unit = "day")
    period$day
  }
  
  cc.start <- cc %>%
    select(ReferenceKey, start,agent)
  
  cc.end <- cc %>%
    select(ReferenceKey, start = end,agent)
  
  cc.startend <- rbind(cc.start, cc.end) %>%
    distinct()
  
  cc.periods <- cc.startend %>%
    group_by(ReferenceKey,agent) %>%
    mutate(end = lead(start, 1, order_by = start)) %>%mutate(end=case_when(is.na(end)~start,!is.na(end)~end))
  
  cc.overlapped <- cc.periods %>%#This process will automatcally delete these patients with only 1 prescrption records as it requires at least 2 records
    left_join(cc, by = "ReferenceKey", suffix = c("_1", "_2")) %>%
    filter(start_1 <= end_2,
           start_2 <= end_1) %>% mutate(agent=agent_1) %>% 
    select(ReferenceKey, agent, start = start_1, end = end_1) %>%
    arrange(ReferenceKey, start, agent) %>%
    group_by(ReferenceKey, start, end,agent) %>%
    group_by(ReferenceKey) %>%
    mutate(end = if_else(end == lead(start, default = as.Date("1970-01-01")), as.Date(end) - 1, end)) 
  
  cc.gap <- cc.overlapped %>%
    mutate(gap = as.numeric(start - lag(end, default = as.Date("1970-01-01"))),
           gap = if_else(gap > 180, 1, gap),
           change = if_else(gap > 180 | agent != lag(agent), 1, 0),#If gap greater than 180 days we consider it discontinued the therapy
           regimen = cumsum(change),
           duration = calc_time(start, end))
  setDT(cc.gap)
  # Keep prescription of first b/ts DMARDs for each patient  
  cc.new <- left_join(cc.gap,initre[,.(indexdate,ReferenceKey,initre=agent)])
  cc.new <- cc.new[agent==initre,.(start=min(start),end=max(end)),.(ReferenceKey,agent)]
  #Some prescription end date exceed 2023-12-31, standardize it to 2023-12-31 (end of follow up)
  cc.new <- cc.new[start<=as.Date("2022-12-31")]

  dur <- data.table(read_xlsx("BtsDMARDs_drug_effect.xlsx"))
  single_assign <- semi_join(cc.new,single,by=c("ReferenceKey","agent"))
  single_assign[,period:=as.numeric(end-start)]
  for (i in 1:length(dur$Agent)) {
    single_assign[str_detect(agent,regex(dur$Agent[i],ignore_case=T)),assign_period:=dur$Treatment_effect_duration[i]]
  }
  single_assign[,true_period:=pmax(period,assign_period,na.rm=T)]
  cc.new <- merge(cc.new,single_assign[,.(ReferenceKey,agent,true_period)],by=c("ReferenceKey","agent"),all.x=T)
  cc.new[!is.na(true_period),end:=start+true_period]
  cc.new[end>as.Date("2023-12-31"),end:=as.Date("2023-12-31")]
  
  ifcode <- setDT(readxl::read_excel("ICD9_opportunisitc_infection.xlsx",sheet = 1))
  infection <- dx[str_detect(code,"^01[0-8]|^031|^117.5|^115|^117.3|^114|^116|^112|^078.5|^070.2|^070.3|^070.41|^070.44|^070.51|^070.54|^070.71|^130|^136.3|^053|^054|^027.0|^48[0-6]|^513|^038|^790.7|^320|^323|^421|^422.92|^590|^711.0|^711.9|^730|^003|^996.66|^68[0-6]|^566|^567|^466|^572.0|^614.[035]|^604|^599.0|^465|^527.8|^380.1|^728.86|^041|^009|^079")]
  infection <-infection[!str_detect(type,regex("o",ignore_case = T))]#remove outpatient diagnosis
  infection <-distinct(infection,ReferenceKey,code,date)#remove duplicate records
  antiinfe <- rx[str_detect(BNFcode,"^5.[123]")]
  trueif <- merge(infection, antiinfe,all.x=T,allow.cartesian = T)
  trueif[,itv:=as.numeric(difftime(DispensingDate,date,units="days"))]
  trueif <- distinct(trueif[itv>=-14&itv<=14][,.(ReferenceKey,code,date)])
  infection <-semi_join(infection,trueif)

  cca <- left_join(cc.new,infection)
  #Case defined as infection occured within 30 days following end of prescription, 30 days is a grace period
  case <-cca[!is.na(code)&date>start&date<end+30&date-start<365,head(.SD,1),ReferenceKey]
  #The rest cohort were selected as control
  control <- cca[!ReferenceKey%in%case$ReferenceKey][,head(.SD,1),ReferenceKey]
  ccc <- rbind(case[,outcome:="Infected"],control[,outcome:="Not_infected"])
  

  
# 2 Generate the variables --------
  
  # Age and sex
  ide <- setDT(readRDS("Identification.RDS"))
  ccc <-merge(ccc,ide[,.(ReferenceKey,sex=Sex,dateofbirth=as.Date(`DateofBirth(yyyy-mm-dd)`))],all.x = T)
  ccc[,age:=as.integer(round((start-dateofbirth)/365,0))]

  # Baseline comorbidities 
  dxcode <- read_xlsx("Features.xlsx",sheet="Co-morbidities")
  dxcc <- merge(dx[ReferenceKey%in%ccc$ReferenceKey],ccc[,.(ReferenceKey,start)],by=c("ReferenceKey"),allow.cartesian = T)
  dxcc <- dxcc[date<start]
  for (i in 1:length(dxcode$Regex)) {
    l <- dxcc[str_detect(code,dxcode$Regex[i])]$ReferenceKey
    ccc[,dxcode$Name[i]:=as.factor(ifelse(ReferenceKey%in%l,1,0))]
    print(i)
  } 
  
  # Remove individuals with age below 18 
  ccc <- ccc[!age<18]
  
  # Lab results
  lab <- readRDS("lab.RDS")
  lab[,LISReferenceDatetime:=as_date(substr(LISReferenceDatetime,1,10))]
  lab[,LISResultNumericResult:=as.numeric(LISResultNumericResult)]
  labbox <- merge(ccc,lab,by = c("ReferenceKey"),all.x = T)
  # Lab test conducted within 90 days prior to b/ts DMARDs initiation
  labbox_90 <- labbox[start-90<=LISReferenceDatetime&LISReferenceDatetime<=start]
  # Standardized the name of lab test, as in CDARS a lab test may have different test names
  # These names were determined manually 
  labname <- c("Albumin","Alanine|ALT","APTT","Aspartate|AST","C3 Complement|Complement C3|C3|Complement 3","C4 Complement|Complement C4|C4|Complement 4","C-reactive protein|CRP","ESR|Erythrocyte sedimentation rate","Ferritin","Fibrinogen","Haematocrit","Hemoglobin|A1C","Indirect Bilirubin|Unconjugated Bilirubin|Bilirubin, Unconjugated|Bilirubin, Indirect","IgA","IgE","IgG","IgM","Lymphocytes|%Lymphocyte|Lymphocyte%|%Lymphocytes|MHA_LYMPHOCYTE, %_%","Platelet|PLT","^Prothrombin Time","Red blood cell|rbc","Creatinine","Total Bilirubin|Bilirubin, total|Bilirubin,Total|Bilirubin (Total)|MBI_Total Bilirubin_umol","^Thrombin Time","Urea","White blood cell|wbc","Anti-nuclear|ANA")
  labnames <- c("Albumin","Alanine transaminase","APTT","Aspartate","Complement 3","Complement 4","CRP","ESR","Ferritin","Fibrinogen","Haematocrit","Hemoglobin","Indirect Bilirubin","IgA","IgE","IgG","IgM","Lymphocytes","Platelet","Prothrombin Time","Red blood cell","Creatinine","Total Bilirubin","^Thrombin Time","Urea","White blood cell","Anti-nuclear")
  #This module reported the missing rate of different lab tests
  sublab <- data.table()
  for (i in 1:length(labname)) {
    z <- labbox_90[str_detect(LISTestDescription,regex(labname[i],ignore_case = T))]
    z$testname <- labnames[i]
    sublab <-rbind(sublab,z)
    print(i)
  } 
  #Remove lab test with missing rate greater than 30%
  completeness_all <- sublab[,round(uniqueN(ReferenceKey)/length(ccc$ReferenceKey),2),testname]
  completeness <- sublab[,round(uniqueN(ReferenceKey)/length(ccc$ReferenceKey),1),testname][V1>=0.7]
  sublab <- sublab[testname%in%completeness$testname]
  # The following code is unique to the Hong Kong data
  # We check if all the lab test have standard test unit
  sublab[,unique(LISTestUnit),testname]#CRP have two different test units due to the naming
  sublab <- sublab[!(str_detect(testname,"CRP")&LISTestUnit!="mg/dL")]#Remove lab test with unwantted unit
  #Now check again we see every test has one unique unit
  Unit <- sublab[,unique(LISTestUnit),testname]
  #Merge the data table
  labcombine <- dcast(sublab[,mean(LISResultNumericResult,na.rm=T),.(testname,ReferenceKey)], ReferenceKey ~ testname,value.var = "V1")
  ccc <- merge(ccc,labcombine,all.x = T)
  colnames(ccc) <- make.names(colnames(ccc))
  
  
  #Medication within 90 days before b/ts DMARDs initiation
  rxcc <- merge(ccc[,.(ReferenceKey,start,end)], rx[ReferenceKey%in%ccc$ReferenceKey],all.x=T,allow.cartesian = T)
  setorder(rxcc,DispensingDate,Disenddate)
  
  calc_time_2 <- function(A1,A2,B1,B2) {
  period <- max(0, min(A2, B2) - max(A1, B1) + 1) 
  return(period)
  }
  
  # #Create matching table
  # nsaid_table <- data.table(Frequency=unique(rx[!(is.na(Dosage)|is.na(DrugStrength)|is.na(DrugFrequency)|is.na(DispensingDate))][str_detect(Route,"ORAL")&!str_detect(DrugStrength,"ML")][str_detect(DrugName,"CELECOXIB|ETORICOXIB|DICLOFENAC|IBUPROFEN|MEFENAMIC|NAPROXEN|PIROXICAM|SULINDAC")]$DrugFrequency))
  # pred_table <- data.table(Frequency=unique(rx[!(is.na(Dosage)|is.na(DrugStrength)|is.na(DrugFrequency)|is.na(DispensingDate))][str_detect(Route,"ORAL")&!str_detect(DrugStrength,"ML")][str_detect(DrugName,"^PREDNISOLONE")]$DrugFrequency))
  # opioid_table <- data.table(Frequency=unique(rx[!(is.na(Dosage)|is.na(DrugStrength)|is.na(DrugFrequency)|is.na(DispensingDate))][str_detect(Route,"ORAL")&!str_detect(DrugStrength,"ML")][str_detect(DrugName,"CODEINE|DEXTROPROPOXYPHENE|DIHYDROCODEINE|METHADONE|MORPHINE|OXYCODONE|TRAMADOL")]$DrugFrequency))
  # frequency_table <- distinct(rbind(nsaid_table,pred_table,opioid_table))
  # write_xlsx(frequency_table,"Frequency_table_ML_Final.xlsx")
  
  #Remove missing value and remove unwanted indications
  rx_nsaid <-rxcc[str_detect(DrugName,"CELECOXIB|ETORICOXIB|DICLOFENAC|IBUPROFEN|MEFENAMIC|NAPROXEN|PIROXICAM|SULINDAC")]
  rx_nsaid <- rx_nsaid[!(is.na(Dosage)|is.na(DrugStrength)|is.na(DrugFrequency)|is.na(DispensingDate))][str_detect(Route,"ORAL")&!str_detect(DrugStrength,"ML")]    

  rx_pred <-rxcc[str_detect(DrugName,"^PREDNISOLONE")]
  rx_pred <-rx_pred[!(is.na(Dosage)|is.na(DrugStrength)|is.na(DrugFrequency)|is.na(DispensingDate))][Route=="ORAL"&!str_detect(DrugStrength,"ML")]

  rx_opioid <- rxcc[str_detect(DrugName,"CODEINE|DEXTROPROPOXYPHENE|DIHYDROCODEINE|METHADONE|MORPHINE|OXYCODONE|TRAMADOL")]
  rx_opioid <-rx_opioid[!(is.na(Dosage)|is.na(DrugStrength)|is.na(DrugFrequency)|is.na(DispensingDate))][Route=="ORAL"&!str_detect(DrugStrength,"ML")]
  
  

  frequency_table <- data.table(readxl::read_excel("frequency_table_ML_Final.xlsx"))
  for (i in 1:length(frequency_table$Frequency_text)) {
    rx_nsaid[DrugFrequency==frequency_table$Frequency_text[i],Frequency:=frequency_table$Frequency_numeric[i]]
    rx_opioid[DrugFrequency==frequency_table$Frequency_text[i],Frequency:=frequency_table$Frequency_numeric[i]]
    rx_pred[DrugFrequency==frequency_table$Frequency_text[i],Frequency:=frequency_table$Frequency_numeric[i]]
    print(i)
  }
  rx_nsaid <- rx_nsaid[!Frequency==0]
  rx_opioid <- rx_opioid[!Frequency==0]
  rx_pred <- rx_pred[!Frequency==0]

  ddd <- setDT(read_excel("Drug DDD.xlsx"))
    
  rx_pred[,T_enddate:=lead(DispensingDate),.(start,ReferenceKey)]
  rx_pred[Disenddate>=T_enddate,Disenddate:=T_enddate-1]
  rx_pred[,overlap:=mapply(calc_time_2,start-90,start,DispensingDate,Disenddate)]
  rx_pred[,strength:=parse_number(DrugStrength)]
  rx_pred[,dosage:=parse_number(Dosage)]
  rx_pred[,truedose:=strength*Frequency*dosage]
  rx_pred <- rx_pred[,.(prednisolone=sum(overlap*truedose,na.rm = T)),.(ReferenceKey,start,end)]
  rx_pred[,prednisolone:=round(prednisolone/90,2),.(ReferenceKey,start,end)]
  ddd_pred <- ddd[Type=="STERIOD"]$`WHO defined daily dose`
  rx_pred$prednisolone <- rx_pred$prednisolone/ddd_pred
  ccc <- merge(ccc,rx_pred,all.x=T,allow.cartesian = T,by=c("ReferenceKey","start","end"))
  ccc[is.na(prednisolone),prednisolone:=as.numeric(0)]
  
  nsaidname <- c("CELECOXIB","ETORICOXIB","DICLOFENAC","IBUPROFEN","MEFENAMIC","NAPROXEN","PIROXICAM","SULINDAC")
  for (i in 1:length(nsaidname)) {
    rx_nsaid[str_detect(DrugName,regex(nsaidname[i],ignore_case = T)),agent:=nsaidname[i]]
    print(i)    
  }    
 
  rx_nsaid[,T_enddate:=lead(DispensingDate),.(start,ReferenceKey,agent)]
  rx_nsaid[Disenddate>=T_enddate,Disenddate:=T_enddate-1]
  rx_nsaid[,overlap:=mapply(calc_time_2,start-90,start,DispensingDate,Disenddate)]
  rx_nsaid[,strength:=parse_number(DrugStrength)]
  rx_nsaid[,dosage:=parse_number(Dosage)]
  rx_nsaid[dosage>5,dosage:=1]#some dosage was written in combination of drug strength and dosage
  rx_nsaid[,truedose:=strength*Frequency*dosage]  
  for (i in 1:length(ddd[Type=="NSAID"]$`WHO defined daily dose`)) {
    rx_nsaid[str_detect(DrugName,regex(ddd[Type=="NSAID"]$Description[i],ignore_case = T)),ddd_nsaid:=ddd[Type=="NSAID"]$`WHO defined daily dose`[i]]
    print(i)
  }
  rx_nsaid <- rx_nsaid[,.(nsaid=sum(overlap*truedose/ddd_nsaid,na.rm = T)),.(ReferenceKey,start,end)]
  rx_nsaid[,nsaid:=round(nsaid/90,2),.(ReferenceKey,start,end)]  
  ccc <- merge(ccc,rx_nsaid,all.x=T,allow.cartesian = T,by=c("ReferenceKey","start","end"))
  ccc[is.na(nsaid),nsaid:=as.numeric(0)] 


  opioidname <- c("CODEINE","DEXTROPROPOXYPHENE","DIHYDROCODEINE","METHADONE","MORPHINE","OXYCODONE","TRAMADOL")
  for (i in 1:length(opioidname)) {
    rx_opioid[str_detect(DrugName,regex(opioidname[i],ignore_case = T)),agent:=opioidname[i]]
    print(i)    
  }  
  
  rx_opioid[,T_enddate:=lead(DispensingDate),.(start,ReferenceKey,agent)]
  rx_opioid[Disenddate>=T_enddate,Disenddate:=T_enddate-1]
  rx_opioid[,overlap:=mapply(calc_time_2,start-90,start,DispensingDate,Disenddate)]
  rx_opioid[,strength:=parse_number(DrugStrength)]
  rx_opioid[,dosage:=parse_number(Dosage)]  
  rx_opioid[,truedose:=strength*Frequency*dosage]  
  for (i in 1:length(ddd[Type=="OPIOID"]$`WHO defined daily dose`)) {
    rx_opioid[str_detect(DrugName,regex(ddd[Type=="OPIOID"]$Description[i],ignore_case = T)),ddd_opioid:=ddd[Type=="OPIOID"]$`WHO defined daily dose`[i]]
    print(i)
  }
  rx_opioid <- rx_opioid[,.(opioid=sum(overlap*truedose/ddd_opioid,na.rm = T)),.(ReferenceKey,start,end)]
  rx_opioid[,opioid:=round(opioid/90,2),.(ReferenceKey,start,end)]  
  ccc <- merge(ccc,rx_opioid,all.x=T,allow.cartesian = T,by=c("ReferenceKey","start","end"))
  ccc[is.na(opioid),opioid:=as.numeric(0)]
  
  #Binary variable, whether experienced infection within previous 3 years
  pif <- merge(ccc[,.(ReferenceKey,start)],infection,all.x=T,by=c("ReferenceKey"))
  ccc[,pastif:=as.factor(ifelse(ReferenceKey%in%pif[date<start&date>start-3*365]$ReferenceKey,1,0))]
  

  
  ifsum <- data.table()
  for (i in 1:length(ifcode$Description)) {
    ct <- copy(ccc)[outcome=="Infected",.N,code][str_detect(code,ifcode$ICD9[i]),.(Des=ifcode$Description[i],icd9=ifcode$ICD9[i],Count=sum(N))]
    ifsum <- rbind(ifsum,ct)
    print(i)
  }# Summarize the type of infection

  ife <- copy(ccc)[,c("start","end","code","date","dateofbirth","ReferenceKey","true_period"):=NULL]
  ife$outcome <- as.factor(ife$outcome)
  ife$outcome <- factor(ife$outcome,levels = c("Not_infected","Infected"))
  
  tableone <- copy(ife)
  s1 <- print(CreateTableOne(data=tableone[,.SD],strata = "outcome",addOverall=T),showAllLevels = F, quote = F, noSpaces = T,test = T,smd = T)  
  
  saveRDS(ife,"data_final.RDS")  
  write_xlsx(list("Infection summary"=ifsum,"Lab unit summary"=Unit,"Missing rate"=completeness_all, "Tableone"=data.table(s1,keep.rownames = T)),"Output_final.xlsx")
  

  
  # 3 Construct machine learning model----
  pacman::p_load(tidyverse,data.table,ggplot2,writexl,RColorBrewer,
                 lubridate,tictoc,doParallel,MLeval,recipes,
                 epitools,fmsb,caret)
  
  ife2 <- readRDS("data_final.RDS")
  
  #1 Remove near zero variance variables
  nzv <- nearZeroVar(ife2,98/2)
  ife3 <-ife2[,.SD,.SDcol=-nzv]
  summary(ife3)  

  # 2 K-means imputation, centering, scaling, and correlation
  imputation_k <-preProcess(ife3,method = c("center", "scale","corr","knnImpute"),K=5)
  pred_k <- predict(imputation_k, ife3)

  # 3 Split the data at a ratio 8:2
  set.seed(126)
  idx1 <- createDataPartition(
    y = pred_k$outcome,#Split according to study outcome
    p = .8,# Percentage 
    list = FALSE
  )
  
  train1 <- pred_k[idx1,]
  test1 <- pred_k[-idx1,]
  
  #Logistic regression
  set.seed(126)
  lg_model1 <- train(
    outcome ~ .,
    data = train1,
    method = "glm",
    family = "binomial",
    trControl = trainControl(
      method = "repeatedcv",
      number = 10,
      repeats = 10,
      sampling = "down",
      classProbs=T,
      savePredictions = T)
  )
  #Plot AUROC curve
  roc_lg <- evalm(lg_model1,plots='r',rlinethick=0.8,fsize=8,positive = "Infected",showplots = TRUE)
  test1$predict_lg1 <- factor(ifelse(predict(lg_model1,test1, type="prob",na.action = na.pass)[2]>0.5,"Infected","Not_infected"),levels = c("Not_infected","Infected"))
  confusionMatrix(test1$predict_lg1,test1$outcome,positive = "Infected" )

  # LASSO regression 
  set.seed(126)
  lasso_model1<-train(  
    outcome ~ ., 
    data = train1, 
    method = "glmnet",
    family="binomial",
    metric="ROC",
    tuneGrid = expand.grid(alpha = 1, lambda = seq(0.01,0.05,by=0.01)),
    trControl = trainControl(
      method = "repeatedcv",
      number = 10,
      repeats = 10,
      sampling = "down",
      classProbs=T,
      savePredictions = T,
      summaryFunction = twoClassSummary)
    ) 
  
  predictors(lasso_model1)
  roc_lasso <- evalm(lasso_model1,plots='r',rlinethick=0.8,fsize=8,positive = "Infected",showplots = TRUE)
  test1$predict_lasso1 <- factor(ifelse(predict(lasso_model1,test1, type="prob")[2]>0.5,"Infected","Not_infected"),levels = c("Not_infected","Infected"))
  confusionMatrix(test1$predict_lasso1,test1$outcome,positive = "Infected")
  
  #Ridge regression 
  set.seed(126)
  ridge_model1<-train(  
    outcome ~ ., 
    data = train1, 
    method = "glmnet",
    family="binomial",
    metric="ROC",
    tuneGrid = expand.grid(alpha = 0, lambda = seq(0.01,0.05,by=0.01)),
    trControl = trainControl(
      method = "repeatedcv",
      number = 10,
      repeats = 10,
      sampling = "down",
      classProbs=T,
      savePredictions = T,
      summaryFunction = twoClassSummary)
  ) 

  roc_ridge <- evalm(ridge_model1,plots='r',rlinethick=0.8,fsize=8,positive = "Infected",showplots = TRUE)
  test1$predict_ridge1 <- factor(ifelse(predict(ridge_model1,test1, type="prob")[2]>0.5,"Infected","Not_infected"),levels = c("Not_infected","Infected"))
  confusionMatrix(test1$predict_ridge1,test1$outcome,positive = "Infected" )

  # Random forest
  cl <- makePSOCKcluster(detectCores() - 4)#Set multiple-cores
  registerDoParallel(cl)
  getDoParWorkers()
  tic()
  set.seed(126)
  rf_model1 <- train(
    outcome ~ ., 
    data = train1, 
    method = "ranger",
    metric = "ROC",
    tuneLength=10,
    trControl = trainControl(
      method = "repeatedcv",
      number = 10,
      repeats = 10,
      sampling = "down",
      classProbs=T,
      savePredictions = T,
      summaryFunction = twoClassSummary)
  )
  toc()
  stopCluster(cl)
  registerDoSEQ()
  
  roc_rf <- evalm(rf_model1,plots='r',rlinethick=0.8,fsize=8,positive = "Infected")
  test1$predict_rf1 <- factor(ifelse(predict(rf_model1,test1, type="prob")[2]>0.5,"Infected","Not_infected"),levels = c("Not_infected","Infected"))
  confusionMatrix(test1$predict_rf1,test1$outcome,positive = "Infected")
  library(epiR)
  rval <- epi.tests(confusionMatrix(table(test1$predict_rf1,test1$outcome),positive = "Infected")$table, conf.level = 0.95,digits=3)
  print(rval)
  
  
  #XGBoost

  tune_grid_xgb <- expand.grid(nrounds = c(50,100,150,200,250),
                          max_depth = c(1,2,3),
                          eta = c(0.05,0.1,0.3,0.5),
                          gamma = 0,
                          min_child_weight=1,
                          colsample_bytree = c(0.6,0.8),
                          subsample = 1)
  
  cl <- makePSOCKcluster(detectCores() - 4)
  registerDoParallel(cl)
  getDoParWorkers()
  tic()
  set.seed(126)
  xgb_model1 <- train(outcome ~., 
      data = train1, 
      method = "xgbTree",
      metric="ROC",
      tuneGrid = tune_grid_xgb,
      trControl=trainControl(
      method = "repeatedcv",
      number = 10,
      repeats = 10,
      sampling = "down",
      classProbs=T,
      savePredictions = T,
      summaryFunction = twoClassSummary)
      )
  toc()
  stopCluster(cl)
  registerDoSEQ()
  
  roc_xgb <- evalm(xgb_model1,plots='r',rlinethick=0.8,fsize=8,positive = "Infected")
  test1$predict_xgb1 <- factor(ifelse(predict(xgb_model1,test1, type="prob")[2]>0.5,"Infected","Not_infected"),levels = c("Not_infected","Infected"))
  confusionMatrix(test1$predict_xgb1,test1$outcome,positive = "Infected" )

  #average Neural network

  tune_grid_nnn <- expand.grid(size=c(1,3,5,7,9),
                              decay=c(0,0.001,0.01,0.1,0.2),
                              bag=10)
  
   cl <- makePSOCKcluster(detectCores() - 4)
    registerDoParallel(cl)
    getDoParWorkers()
    tic()
    set.seed(126)
    nnn_model1 <- train(outcome ~., 
        data = train1, 
        method = "avNNet",
        metric="ROC",
        tuneGrid=tune_grid_nnn,
        trControl=trainControl(
        method = "repeatedcv",
        number = 10,
        repeats = 10,
        sampling = "down",
        classProbs=T,
        savePredictions = T,
      summaryFunction = twoClassSummary)
        )
    toc()
    stopCluster(cl)
    registerDoSEQ()

    roc_nnn <- evalm(nnn_model1,plots='r',rlinethick=0.8,fsize=8,positive = "Infected")
    test1$predict_nnn1 <- factor(ifelse(predict(nnn_model1,test1, type="prob")[2]>0.5,"Infected","Not_infected"),levels = c("Not_infected","Infected"))
    confusionMatrix(test1$predict_nnn1,test1$outcome,positive = "Infected" )

  #SVM
  
   cl <- makePSOCKcluster(detectCores() - 4)
    registerDoParallel(cl)
    getDoParWorkers()
    tic()
    set.seed(126)
    svm_model1 <- train(outcome ~., 
        data = train1, 
        method = "svmLinear",
        metric="ROC",
        trControl=trainControl(
        method = "repeatedcv",
        number = 10,
        repeats = 10,
        sampling = "down",
        classProbs=T,
        savePredictions = T,
      summaryFunction = twoClassSummary)
        )
    toc()
    stopCluster(cl)
    registerDoSEQ()

    roc_svm <- evalm(svm_model1,plots='r',rlinethick=0.8,fsize=8,positive = "Infected")
    test1$predict_svm1 <- factor(ifelse(predict(svm_model1,test1, type="prob")[2]>0.5,"Infected","Not_infected"),levels = c("Not_infected","Infected"))
    confusionMatrix(test1$predict_svm1,test1$outcome,positive = "Infected" )    
    
  #Save all the trained model for quick use in the next time
  save(lg_model1,rf_model1,xgb_model1,nnn_model1,ridge_model1,lasso_model1,svm_model1, file = "output_algorithm_final.RData")
    
 
  # Observe the ROC plots for each model
  library(pROC)
  dev_lg <- roc(test1$outcome,predict(lg_model1,test1, type="prob")$`Infected`,plot=T,legacy.axes=T,col="#EF476F",lwd=3,ci=T,print.auc=T)
  dev_lasso <- roc(test1$outcome,predict(lasso_model1,test1, type="prob")$`Infected`,plot=T,legacy.axes=T,col="#EF476F",lwd=3,ci=T,print.auc=T)
  dev_ridge <- roc(test1$outcome,predict(ridge_model1,test1, type="prob")$`Infected`,plot=T,legacy.axes=T,col="#EF476F",lwd=3,ci=T,print.auc=T)
  dev_rf <- roc(test1$outcome,predict(rf_model1,test1, type="prob")$`Infected`,plot=T,legacy.axes=T,col="#EF476F",lwd=3,ci=T,print.auc=T)
  dev_xgb <- roc(test1$outcome,predict(xgb_model1,test1, type="prob")$`Infected`,plot=T,legacy.axes=T,col="#EF476F",lwd=3,ci=T,print.auc=T)
  dev_nnn <- roc(test1$outcome,predict(nnn_model1,test1, type="prob")$`Infected`,plot=T,legacy.axes=T,col="#EF476F",lwd=3,ci=T,print.auc=T)
  dev_svm <- roc(test1$outcome,predict(svm_model1,test1, type="prob")$`Infected`,plot=T,legacy.axes=T,col="#EF476F",lwd=3,ci=T,print.auc=T)

  # #Plot ROC curve
  jpeg("ROC_final.jpg",width = 3600,height = 2400,res = 300)
  par(pty="s")
  roc(test1$outcome,predict(lg_model1,test1, type="prob")$`Infected`,plot=T,legacy.axes=T,col="#ef476f",lwd=3)
  plot.roc(test1$outcome,predict(lasso_model1,test1, type="prob")$`Infected`,col="#f78c6b",lwd=3,add=T)
  plot.roc(test1$outcome,predict(svm_model1,test1, type="prob")$`Infected`,col="#ffd166",lwd=3,add=T)
  plot.roc(test1$outcome,predict(xgb_model1,test1, type="prob")$`Infected`,col="#06d6a0",lwd=3,add=T)
  plot.roc(test1$outcome,predict(rf_model1,test1, type="prob")$`Infected`,col="#073b4c",lwd=3,add=T)
  plot.roc(test1$outcome,predict(nnn_model1,test1, type="prob")$`Infected`,col="#118ab2",lwd=3,add=T)
  legend("bottomright",legend=c("Logistic Regression: AUROC = 0.830 (0.778-0.882)","LASSO Regression: AUROC = 0.836 (0.786-0.886)",
                                "Support vector machines: AUROC = 0.833 (0.782-0.884)", "XGBoost: AUROC = 0.795 (0.741-0.849)",
                                "Random Forest: AUROC = 0.840 (0.793-0.888)","Neural Network: AUROC = 0.804 (0.747-0.862)"),
         col=c("#ef476f","#f78c6b","#ffd166","#06d6a0","#073b4c","#118ab2"),lwd=3)
  dev.off()
  # 
  #Generate Shap value----
  packageurl <- "https://cran.r-project.org/src/contrib/Archive/kernelshap/kernelshap_0.4.1.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
  library(kernelshap)
  library(shapviz)
  library(kernlab)  
  
  set.seed(420)
  explaindata <- as.data.frame(train1[,.SD,.SDcol=colnames(train1)[-2]])#Explain data
  bkdata <- as.data.frame(train1[sample(nrow(train1), 200)])#Background data
  #In previous step, random forest was selected as the optimal model
  #It was therefore used to generate the SHAP value here
  #It takes around 3 hours using the Hong Kong data
  cl <- makePSOCKcluster(detectCores() - 4)
  registerDoParallel(cl)
  getDoParWorkers()
  tic()  
  shap_v <- kernelshap(rf_model1, explaindata, bg_X = bkdata)
  stopCluster(cl)
  registerDoSEQ()  
  toc() 
  
  save(shap_v, file = "shapval_final.RData")
 
  shap_v2 <- shapviz(object = shap_v)
  colnames(shap_v2$Infected) <- c("B/ts DMARDs","Sex","Age","Coronary artery disease","Hyperlipidemia","Thyroid disorder","Osteoporosis","Hypertension",
                                  "Diabetes mellitus","Cerebrovascular disease","COPD",
                                  "Ulcers","Malignancy","Alanine transaminase","Albumin","CRP","Creatinine","ESR",
                                  "Lymphocytes","Platelet","RBC", "Bilirubin","Urea","WBC","Prednisolone","NSAID","Opioid","Previous infection")
  p <- sv_importance(shap_v2$Infected, kind = "beeswarm",show_numbers = F,max_display = 30)
  
  #Plot beeswarm plot for SHAP value, keep the top 15 important factors
  jpeg("SHAP_ridge_beeswarm_final.jpg",width = 3840,height = 2160,res = 300)
  p+theme_classic() + theme(
    axis.text = element_text(size = 14),       # Axis text font size
    axis.title = element_text(size =14),     # Axis title font size
    plot.title = element_text(size = 14),     # Plot title font size
    legend.text = element_text(size = 14),    # Legend text font size
    legend.title = element_text(size = 14)    # Legend title font size
  )
  dev.off()  
  
  #Shap value explore categorical variable - b/ts DMARDs
  btsshap <- data.table(agent=shap_v$X$agent,shap=shap_v$S$Infected[,1])
  btsshap[,median:=median(shap),agent]
  tiff("SHAP_agent_final.jpg",width = 3840,height = 2160,res = 300,compression = "jpeg")
  ggplot(aes(x = reorder(agent,median),y=shap,fill=agent),data = btsshap)+geom_violin()+stat_summary(fun=median, geom="point", size=2, color="gray7")+
    scale_fill_brewer(palette="Paired") + theme(legend.position="none")+xlab("B/ts DMARDs")+ylab("SHAP value")+ 
    theme(axis.text.x = element_text(angle=45,hjust = 1),panel.border = element_blank(),
          panel.grid.major = element_blank(),    axis.text = element_text(size = 12),       # Axis text font size
          axis.title = element_text(size =14),     # Axis title font size
          plot.title = element_text(size = 14),     # Plot title font size
          legend.text = element_text(size = 14),    # Legend text font size
          legend.title = element_text(size = 14), 
          panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white"))+geom_hline(yintercept = 0,linetype="dashed")
  dev.off()  
  
  #Draw spider plot for phenotype among infected and non-infected population ----
  #Get the top 15 features' names
  imp <- rownames(shapviz:::.get_imp(get_shap_values(shap_v2)))
  #Select the continuous variable manually
  imp_ct <- c("Albumin","Creatinine","CRP","ESR","Lymphocytes","Platelet","Red.blood.cell","Urea","White.blood.cell")
  imp_ct_rader <- ife3[,.SD,.SDcol=c("outcome",imp_ct)]
  
  datasp <- rbind(as.data.frame(imp_ct_rader[,-c("outcome")][,lapply(.SD,quantile,0.95,na.rm=T)][,lapply(.SD, round,0)]),
  rep(0,9),
  as.data.frame(imp_ct_rader[,lapply(.SD,mean,na.rm=T),outcome][,-c("outcome")]))
  rownames(datasp)[3:4] <- c("Not_infected","Infected")
  # Color vector
  colors_border=c( rgb(0.2,0.5,0.5,0.6), rgb(0.7,0.5,0.1,0.6) )
  colors_in=c( rgb(0.2,0.5,0.5,0.2) , rgb(0.7,0.5,0.1,0.2) )
  
  # plot with default options:
  tiff("Spider_plot_final.jpg",width = 3840,height = 2160,res = 300,compression = "jpeg")
  radarchart( datasp  , axistype=2 ,
      #custom polygon
      pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
      #custom the grid
      cglcol="black", cglty=1, axislabcol="black", cglwd=0.8,
      #custom labels
      vlcex=1.2,palcex = 1.2,vlabels=c("Albumin (g/L)","Creatinine	(umol/L)     ", "CRP (mg/dL)             ","ESR (mm/hr)        ", "Lymphocytes (%)","Platelet (10^9/L)","                       Red blood cell (10^12/L)","             Urea (mmol/L)","    White blood cell	(10^9/L)") 
      )
  legend(x=2, y=1, legend = rownames(datasp[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "black", cex=1.2, pt.cex=3)
  dev.off()

  # Find cutoff thresholds with maximum odds ratio----
  
  F1cutoff1 <- data.table()
  biomarkers1 <- c("CRP","Creatinine","ESR","Urea","White.blood.cell")

  for (i in 1:length(biomarkers1)) {
    data <- as.data.frame(ife3[,.(outcome,feature=get(biomarkers1[i]))][!is.na(feature)])
    thresholds <- seq(quantile(data$feature,0.05), quantile(data$feature,0.95),by=0.1)
    for (j in 1:length(thresholds)) {
    data$group <- ifelse(data$feature > thresholds[j],  "Infected", "Not_infected")#greater than threshold was considered as risk factors
    data$outcome <- factor(data$outcome, levels = c("Infected","Not_infected")) 
    data$group <- factor(data$group, levels = c("Infected","Not_infected")) 
    comat <- confusionMatrix(data=data$group,reference=data$outcome)
    F1 <- data.table(Description=biomarkers1[i],
                        Threshold=thresholds[j],
                        F1=round(comat$byClass[7],4))
    F1cutoff1 <- rbind(F1cutoff1,F1)
    } 
}
  
  F1cutoff2 <- data.table()
  biomarkers2 <- c("Albumin","Lymphocytes","Red.blood.cell","Platelet")

  for (i in 1:length(biomarkers2)) {
    data <- as.data.frame(ife3[,.(outcome,feature=get(biomarkers2[i]))][!is.na(feature)])
    thresholds <- seq(quantile(data$feature,0.05), quantile(data$feature,0.95),by=0.1)
    for (j in 1:length(thresholds)) {
    data$group <- ifelse(data$feature > thresholds[j],  "Infected", "Not_infected")#greater than threshold was considered as risk factors
    data$outcome <- factor(data$outcome, levels = c("Infected","Not_infected")) 
    data$group <- factor(data$group, levels = c("Infected","Not_infected")) 
    comat <- confusionMatrix(data=data$group,reference=data$outcome)
    F1 <- data.table(Description=biomarkers2[i],
                        Threshold=thresholds[j],
                        F1=round(comat$byClass[7],4))
    F1cutoff2 <- rbind(F1cutoff2,F1)
    } 
}
  
  F1sum1 <- F1cutoff1 %>%group_by(Description) %>%  filter(F1==max(F1)) %>% distinct(F1,.keep_all = T)
  F1sum2 <- F1cutoff2 %>%group_by(Description) %>%  filter(F1==max(F1)) %>% distinct(F1,.keep_all = T)
  write_xlsx(list("F1_1"=acsum1,"F1_2"=acsum2),"output_final_threshold.xlsx")
  
  
  
  
  #Sample size calculation
  library(pmsampsize)
  pmsampsize(type = "b", 
           cstatistic = 0.84, 
           parameters = 28, 
           prevalence = 277/3159)
  # library(pmvalsampsize)
  # pmvalsampsize(type = "b", 
  #          cstatistic = 0.73, 
  #          prevalence = 52/1845,lpnorm = c(-5,2.5), oeciwidth = 1, simobs = 1000)
  
  
#Minimum sample size required for new model development based on user inputs = 2372, 
#with 231 events (assuming an outcome prevalence = 0.097) and an EPP = 11.5 
  
  
  #Plot external validation----
  allofus <- read.csv("aou_rf.csv",header=T)  
  library(pROC)
  
  roc(allofus$outcome,allofus$`Infected`,plot=T,legacy.axes=T,col="#EF476F",lwd=3,ci=T,print.auc=T)
  
  allofus$predict <- factor(ifelse(allofus$Infected>0.5,"Infected","Not_infected"),levels = c("Not_infected","Infected"))
  allofus$outcome <- factor(allofus$outcome,levels = c("Not_infected","Infected"))
  confusionMatrix(allofus$predict,allofus$outcome,positive = "Infected" )
  
  test1$predict_rf1 <- factor(ifelse(predict(rf_model1,test1, type="prob")[2]>0.5,"Infected","Not_infected"),levels = c("Not_infected","Infected"))
  confusionMatrix(test1$predict_rf1,test1$outcome,positive = "Infected")
  
  library(epiR)
  rval <- epi.tests(confusionMatrix(table(allofus$predict,allofus$outcome),positive = "Infected")$table, conf.level = 0.95,digits=3)
  print(rval)
  
  jpeg("ROC_ex.jpg",width = 3600,height = 2400,res = 300)
  par(pty="s")
  roc(allofus$outcome,allofus$`Infected`,plot=T,legacy.axes=T,col="#1E1E2F",lwd=3)
  plot.roc(test1$outcome,predict(rf_model1,test1, type="prob")$`Infected`,col="#F9A620",lwd=3,add=T)
  legend("bottomright",legend=c("AllofUs: AUROC = 0.729 (0.665-0.793)","CDARS: AUROC = 0.840 (0.793-0.888)"),
         col=c("#1E1E2F","#F9A620"),lwd=3)
  dev.off()  
  
  
  jpeg("ROC_combine.jpg",width = 4800,height = 2400,res = 300)
  par(mfrow = c(1, 2))
  par(pty="s")
  f1 <- roc(test1$outcome,predict(lg_model1,test1, type="prob")$`Infected`,plot=T,legacy.axes=T,col="#ef476f",lwd=3,  cex.axis = 1.5,cex.lab = 1.5)
  plot.roc(test1$outcome,predict(lasso_model1,test1, type="prob")$`Infected`,col="#f78c6b",lwd=3,add=T)
  plot.roc(test1$outcome,predict(svm_model1,test1, type="prob")$`Infected`,col="#ffd166",lwd=3,add=T)
  plot.roc(test1$outcome,predict(xgb_model1,test1, type="prob")$`Infected`,col="#06d6a0",lwd=3,add=T)
  plot.roc(test1$outcome,predict(rf_model1,test1, type="prob")$`Infected`,col="#073b4c",lwd=3,add=T)
  plot.roc(test1$outcome,predict(nnn_model1,test1, type="prob")$`Infected`,col="#118ab2",lwd=3,add=T)
  legend("bottomright",legend=c("Logistic Regression: AUROC = 0.830 (0.778-0.882)","LASSO Regression: AUROC = 0.836 (0.786-0.886)",
                                "Support vector machines: AUROC = 0.833 (0.782-0.884)", "XGBoost: AUROC = 0.795 (0.741-0.849)",
                                "Random Forest: AUROC = 0.840 (0.793-0.888)","Neural Network: AUROC = 0.804 (0.747-0.862)"),
         col=c("#ef476f","#f78c6b","#ffd166","#06d6a0","#073b4c","#118ab2"),lwd=3,cex = 1)
  mtext("A", side = 3, line = 0.5,adj= 0.05, cex = 1.5, font = 2)  # Add Figure 1A to the top-left

  par(pty="s")
  f2 <- roc(allofus$outcome,allofus$`Infected`,plot=T,legacy.axes=T,col="#1E1E2F",lwd=3,  cex.axis = 1.5,cex.lab = 1.5)
  plot.roc(test1$outcome,predict(rf_model1,test1, type="prob")$`Infected`,col="#F9A620",lwd=3,add=T)
  legend("bottomright",legend=c("AllofUs: AUROC = 0.729 (0.665-0.793)","CDARS: AUROC = 0.840 (0.793-0.888)"),
         col=c("#1E1E2F","#F9A620"),lwd=3,cex = 1)
  mtext("B", side = 3, line = 0.5, adj = 0.05, cex = 1.5, font = 2)  # Add Figure 1B to the top-left  
  dev.off()  


#End---------------------  
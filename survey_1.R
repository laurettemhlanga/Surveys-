

Survey_1 <- function(n, 
                     #N_0,
                     prevalence,
                     standdev,
                     Enumaration_Areas
                     ,Prob_Recency_ptve
                     ){
  # takes a vector of enumertaion sizes, selects n using pps, assigns the respective prevalence 
  # produces a  data frame with sampling weights ea size and prevalence
  # browser()
Total_population <- sum(Enumaration_Areas)

Sampling_probabilities <- Enumaration_Areas / Total_population

#Enumaration_Areas_sample = (Enumaration_Areas/sum(Enumaration_Areas)) * N_0

Study_EAs <- sample(x = Enumaration_Areas, 
                    size = n, 
                    prob = Sampling_probabilities, 
                    replace = F)

Sampling_fraction = Study_EAs / Total_population

rsd_strata = standdev/prevalence

prevalences_strata = rnorm(length(Study_EAs), mean = prevalence, sd = standdev)
cof_var_strata= sd(prevalences_strata ) / mean(prevalences_strata)

prevalences_strata  = mean(prevalences_strata ) + (prevalences_strata  - mean(prevalences_strata ))*(rsd_strata/cof_var_strata)


factor_strata = sum(Study_EAs * prevalences_strata) / sum(Study_EAs)
prevalences_strata = (prevalences_strata / factor_strata) * prevalence

prevcheck = sum(Study_EAs * prevalences_strata) / sum(Study_EAs)

Survey_data = data.frame(EA_id = 1: n,
                         EA_size = Study_EAs,
                         Sampling_weights = 1 / Sampling_fraction,
                         Prevalence_EA = prevalences_strata
                         ,Prevalence_EA_Recency = prevalences_strata * Prob_Recency_ptve
                         )

return(list(check = prevcheck, Survey_data = Survey_data))
}




Survey_2 <- function(n, 
                     N_0,
                     prevalence,
                     standdev,
                     Enumaration_Areas_sizes
                     ,Prob_Recency_ptve
){
  # takes a vector of enumertaion sizes, selects n using pps, assigns the respective prevalence 
  # produces a  data frame with sampling weights ea size and prevalence
  # browser()
  Total_population <- sum(Enumaration_Areas_sizes)
  
  EA_Sampling_probabilities <- Enumaration_Areas_sizes / Total_population
  
  Study_EAs <- sample(x = Enumaration_Areas_sizes, 
                      size = n, 
                      prob = EA_Sampling_probabilities, 
                      replace = F)
  
  
  EA_Sampling_fraction = Study_EAs / Total_population
  #####here 
  Enumaration_Areas_sample = round((Study_EAs/sum(Study_EAs)) * N_0)
  
  
  
  Sample_Sampling_fraction = Enumaration_Areas_sample  / sum(Enumaration_Areas_sample)
  
  CoV_overall = standdev/prevalence
  
  prevalences_sampled_R = rnorm(length(Study_EAs), mean = prevalence, sd = standdev)
  CoV_Sampled_R = sd(prevalences_sampled_R ) / mean(prevalences_sampled_R)
  
  prevalences_sampled_R_Adj = mean( prevalences_sampled_R ) + (prevalences_sampled_R - mean( prevalences_sampled_R )) * (CoV_overall / CoV_Sampled_R)
  
  
  factor_sampled_R = sum(Enumaration_Areas_sample * prevalences_sampled_R_Adj) / sum(Enumaration_Areas_sample)
  prevalences_final = (prevalences_sampled_R_Adj / factor_sampled_R ) * prevalence
  
  prevcheck = sum(Enumaration_Areas_sample * prevalences_final) / sum(Enumaration_Areas_sample)
  
  EA_std_dev = sqrt(((Study_EAs - Enumaration_Areas_sample)/(Study_EAs*( Enumaration_Areas_sample - 1))) * (prevalences_final *(1-prevalences_final)))
    
  Survey_data = data.frame(EA_id = 1: n,
                           EA_size = Study_EAs,
                           EA_sample = Enumaration_Areas_sample,
                           EA_Sampling_weights = 1 / EA_Sampling_fraction,
                           Sample_Sampling_weights = 1 / Sample_Sampling_fraction,
                           Base_weight =  1/( EA_Sampling_fraction * Sample_Sampling_fraction), 
                           Prevalence_EA = prevalences_final,
                           EA_std_dev = EA_std_dev
                           ,Prevalence_EA_Recency = prevalences_final * Prob_Recency_ptve
  )
  
  return(list(check = prevcheck, Survey_data = Survey_data))
}



##########################Above embedded in one function###################################################
###########################################################################################################


Generate_Prevalence <- function(n, 
                                N_0,
                                prevalence,
                                standdev,
                                Enumaration_Areas
                                ,Prob_Recency_ptve
                                ){
  #genertaes prevalences for the respecive sampled areas 
  #has an option for when everyone in the given sampled region is interviewed vs 
  #when a subset of the sampled region is interviwed and the sample size should be equal to a designated samplesize
  #outputs dataframe with weights and  per EA and weights per indididual 
  Total_population <- sum(Enumaration_Areas)
  
  EA_Sampling_probabilities <- Enumaration_Areas / Total_population
  
  Study_EAs <- sample(x = Enumaration_Areas, 
                      size = n, 
                      prob = EA_Sampling_probabilities, 
                      replace = F)
  
  
  EA_Sampling_fraction = Study_EAs / Total_population
  if (is.null(N_0) == TRUE){
    
    CoV_overall = standdev/prevalence
    
    prevalences_sampled_R = rnorm(length(Study_EAs), mean = prevalence, sd = standdev)
    CoV_Sampled_R = sd(prevalences_sampled_R ) / mean(prevalences_sampled_R)
    
    prevalences_sampled_R_Adj  = mean(prevalences_sampled_R ) + (prevalences_sampled_R  - mean(prevalences_sampled_R  ))*(CoV_overall/CoV_Sampled_R)
    
    
    factor_sampled_R = sum(Study_EAs * prevalences_sampled_R_Adj) / sum(Study_EAs)
    prevalences_final = (prevalences_sampled_R_Adj/ factor_sampled_R) * prevalence
    
    prevcheck = sum(Study_EAs * prevalences_final) / sum(Study_EAs)
    
    Survey_data = data.frame(EA_id = 1: n,
                             EA_size = Study_EAs,
                             Sampling_weights = 1 / EA_Sampling_fraction,
                             Prevalence_EA = prevalences_final
                             ,Prevalence_EA_Recency = prevalences_final * Prob_Recency_ptve)
    
    survey = list(check = prevcheck, Survey_data = Survey_data)
  }else{
    
    
    Enumaration_Areas_sample = round((Study_EAs/sum(Study_EAs)) * N_0)
    
    
    
    Sample_Sampling_fraction = Enumaration_Areas_sample  / sum(Enumaration_Areas_sample)
    
    CoV_overall = standdev/prevalence
    
    prevalences_sampled_R = rnorm(length(Study_EAs), mean = prevalence, sd = standdev)
    CoV_Sampled_R = sd(prevalences_sampled_R ) / mean(prevalences_sampled_R)
    
    prevalences_sampled_R_Adj = mean( prevalences_sampled_R ) + (prevalences_sampled_R - mean( prevalences_sampled_R )) * (CoV_overall / CoV_Sampled_R)
    
    
    factor_sampled_R = sum(Enumaration_Areas_sample * prevalences_sampled_R_Adj) / sum(Enumaration_Areas_sample)
    prevalences_final = (prevalences_sampled_R_Adj / factor_sampled_R ) * prevalence
    
    prevcheck = sum(Enumaration_Areas_sample *  prevalences_final) / sum(Enumaration_Areas_sample)
    
    Survey_data = data.frame(EA_id = 1: n,
                             EA_sample = Enumaration_Areas_sample,
                             Primary_Sampling_weights = 1 / EA_Sampling_fraction,
                             Unit_Sampling_weights = 1 / Sample_Sampling_fraction,
                             Base_weight = (1 / EA_Sampling_fraction) * (1 / Sample_Sampling_fraction),
                             Prevalence_EA = prevalences_final
                             ,Prevalence_EA_Recency = prevalences_final * Prob_Recency_ptve
    )
    survey = list(check = prevcheck, Survey_data = Survey_data)
  }
  return(survey)
}


survey_data = Survey_2(n = 11, 
              N_0 = 20000,
              prevalence = 0.02,
              standdev = 0.01,
              Enumaration_Areas_sizes= sample(2000:5000, 20, T)
             ,Prob_Recency_ptve = 0.01
              )


survey_data = Generate_Prevalence(n = 11,
              N_0 = 20000,        
              prevalence = Prop_Prevalence(23),
              standdev =  0.04,
              Enumaration_Areas = sample(2000:5000, 20, T)
              ,Prob_Recency_ptve = Prop_Prev_Recency(23))





survey_data = Generate_Prevalence(n = 11,
                       N_0 = 20000,        
                       prevalence = Prop_Prevalence(23),
                       standdev =  0.04,
                       Enumaration_Areas = sample(2000:5000, 20, T)
                       ,Prob_Recency_ptve = Prop_Prev_Recency(23))






#population prevalence and standard error 

Prevalence_weighted <- function(dataset
                                ){
 
   weighted_prevalence = sum(dataset$Base_weight * dataset$Prevalence_EA) / sum(dataset$Base_weight)
   
   weighted_standard_error =  sum(((dataset$Sample_Sampling_weights )^ 2) * ((dataset$EA_size - dataset$EA_sample)/(dataset$EA_size*(dataset$EA_sample - 1))) * dataset$Prevalence_EA * (1 - dataset$Prevalence_EA))
   return(c(Pop_prev = weighted_prevalence,  Pop_prev_se = weighted_standard_error ))
}

dataset = survey_data$Survey_data

Prevalence_weighted(survey_data$Survey_data)


library(survey)


y = survey_data$Survey_data

data1 <- svydesign(id = ~EA_id, weights = ~Sample_Sampling_weights, data = y)

#individualing the data set 
svymean(~Prevalence_EA, data1)


Individual_level_data = function(Surveydata
                                 ){
  #takes in the survey data in 1 and disintergrates it into individual 
 data = NULL
  for (i in 1:nrow(survey_data$Survey_data)){
  Rawdata = data.frame(id = 1:survey_data$Survey_data$EA_size[i],
                     EA_id = rep(survey_data$Survey_data$EA_id[i],length(1:survey_data$Survey_data$EA_size[i])),
                     weight = rep(survey_data$Survey_data$Base_weight[1], length(1:survey_data$Survey_data$EA_size[i])),
                     HIV_status = sample(x=c(0,1), size = survey_data$Survey_data$EA_size[i], replace = T, prob = c((1-survey_data$Survey_data$Prevalence_EA[i]), survey_data$Survey_data$Prevalence_EA[i])))
  
  Rawdata$HIV_Recency_Ptve = ifelse(Rawdata$HIV_status == 1,
                                  sample(x = c(0,1), size = survey_data$Survey_data$EA_size[i], replace = T, prob = c((1-survey_data$Survey_data$Prevalence_EA_Recency[i]), survey_data$Survey_data$Prevalence_EA_Recency[i])),
                                  0)
  
  data = rbind(data, Rawdata)
  }
  data = cbind(Unique_id = 1:sum(survey_data$Survey_data$EA_size), data)
  return(data)
}


Individual_level_data(survey_data$Survey_data)


#############################################################################################################################################
#############################################################################################################################################
#using survey genertaed data in the "Midpoint_inc_from_prevalence" function 



survey_data1 = Survey_2(n = 11,
                        N_0 = 20000,        
                        prevalence = Prop_Prevalence(17),
                        standdev =  0.04,
                        Enumaration_Areas = sample(2000:5000, 20, T)
                        ,Prob_Recency_ptve = Prop_Prev_Recency(17))
  
survey_data2 = Survey_2(n = 11,
                        N_0 = 20000,        
                        prevalence = Prop_Prevalence(23),
                        standdev =  0.04,
                        Enumaration_Areas = sample(2000:5000, 20, T)
                        ,Prob_Recency_ptve = Prop_Prev_Recency(23))





data1 <- svydesign(id = ~EA_id,  weights = ~Base_weight,data = survey_data1$Survey_data)
data2 <- svydesign(id = ~EA_id,  weights = ~Base_weight,data = survey_data2$Survey_data)

#data2 <- svydesign(strata = ~strata_number, id = ~EA_Number,  weights = ~sampling_weights, data = z, nest=TRUE)

# survey_data$stratification <- sample(1:3, 20, replace = TRUE)
# survey_data$cluster <- sample(1:10, 20, replace = TRUE)



summary(data1)
summary(data2)


estimated_prvalence1 <- svymean(~Prevalence_EA, data1, deff = TRUE)
estimated_prvalence2 <- svymean(~Prevalence_EA, data2, deff = TRUE)



estimated_prvalence_Rec1 <- svymean(~Prevalence_EA_Recency, data2, deff = TRUE)
estimated_prvalence_Rec2 <- svymean(~Prevalence_EA_Recency, data2, deff = TRUE)





incidi = Midpoint_inc_from_prevalence(T1 = 17, T2 = 22, excess_mortality_estimate = 0.1, 
                                      pop_prevalence = c(0.044153, 0.194458),
                                      DE1 = 0.8636, DE2 = 1.418, sigma_prev_1 = 0.010438 , 
                                      sigma_prev_2 = 0.015587, sample_size_1 = 104063, 
                                      sample_size_2 = 123842)

incidi_2 = Midpoint_inc_from_prevalence(T1 = 17, T2 = 22, excess_mortality_estimate = 0.1, 
                                      pop_prevalence = Prop_Prevalence,
                                      DE1 = 1, DE2 = 1, sigma_prev_1 = 0 , 
                                      sigma_prev_2 = 0, sample_size_1 = 104063, 
                                      sample_size_2 = 123842)









 

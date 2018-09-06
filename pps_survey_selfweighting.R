

#written by laurette Mhlanga an adaptation on the sampling with pps 
#psu to achieve equal weights in a 
#for every individual 



partition_prevalence <- function(overall_prevalence1,
                                 overall_prevalence2,
                                 sigma_prevalence,
                                 cluster_number,
                                 cluster_size
 #function that divides/partitions  the overall prevalence and assigns to clusters
 #into given clusters                                                                  
                                 
){
  
  CoV_overall =  sigma_prevalence/overall_prevalence1
  
  cluster_prevalences = rnorm(length(cluster_number), 
                              mean = overall_prevalence1, 
                              sd = sigma_prevalence)
  # to avoid negative prevalences we check is there are any negatives in the prevalences vector and if they are we replace
  #the negative values by the mean of the simulated cluster_prevalences
  
  if (any(cluster_prevalences < 0)){
    
    x = replace(cluster_prevalences, which(cluster_prevalences<0), mean(cluster_prevalences))
    
  }else{
    
    cluster_prevalences = cluster_prevalences 
  }
  
  
  
  cluster_CoV = sd(cluster_prevalences) / mean(cluster_prevalences)
  
  cluster_prevalences_R_Adj = mean( cluster_prevalences ) + (cluster_prevalences - mean(cluster_prevalences )) * (CoV_overall / cluster_CoV)
  
  
  Cluster_factor = sum(cluster_size * cluster_prevalences_R_Adj) / sum(cluster_size)
  cluster_prevalences_final1 = (cluster_prevalences_R_Adj / Cluster_factor ) * overall_prevalence1
  
  cluster_prevalences_final2 = (cluster_prevalences_final1)*(overall_prevalence2 / overall_prevalence1)
  
  prevcheck = c(sum(cluster_size * cluster_prevalences_final1) / sum(cluster_size), sum(cluster_size * cluster_prevalences_final2) / sum(cluster_size))
   
  return(list(prevcheck, data.frame( cluster_number = cluster_number, 
                                     cluster_prevalence_t1 = cluster_prevalences_final1, 
                                     cluster_prevalence_t2 = cluster_prevalences_final2 )))
  
}







Survey_pps <- function(cluster_number,
                       cluster_size , 
                       num_cluster_sample,
                       ind_per_cluster, 
                       overall_prevalence1,
                       overall_prevalence2,
                       sigma_prevalence
){
  
  #Steps in appling Probability propotional to size(pps) and ways to ensure equal sampling weights 
  #individual (calculation of basic weights)
  
  survey_data = data.frame(cluster_id = cluster_number,
                           cluster_population = cluster_size)
  
  survey_data$cumulative_sum = cumsum(survey_data$cluster_population)
  
  
  survey_data = cbind(survey_data, partition_prevalence(overall_prevalence1 = overall_prevalence1,
                                                        overall_prevalence2 = overall_prevalence2,
                                                        sigma_prevalence = sigma_prevalence,
                                                        cluster_number = cluster_number,
                                                        cluster_size = cluster_size)[[2]][, -1])
  
  
  sampling_interval = survey_data$cumulative_sum[length(cluster_number)] / num_cluster_sample
  
  threshold = survey_data$cumulative_sum[max(which(survey_data$cumulative_sum < sampling_interval))]
  
  #threshold added to avoid getting NAs this ensures the first random number is not greater than
  # the sampling interval
  
  random_number = sample(1 : threshold, 1)
  
  cluster_series = cumsum(c(random_number, rep(sampling_interval, num_cluster_sample-1)))
  
  id_clusters_sampled = as.vector(rep(NA, num_cluster_sample))
  
  counter = 1
  for (tt in cluster_series){
    
    id_clusters_sampled[counter] = which(survey_data$cumulative_sum > tt)[1]
    
    counter = counter + 1
  }
  
  survey_data = survey_data[id_clusters_sampled, ] 
  
  survey_data$cluster_series = cluster_series
  
  survey_data$sampling_fraction_1 = (survey_data$cluster_population * num_cluster_sample) / sum(cluster_size)
  
  survey_data$ind_per_cluster = rep(ind_per_cluster, length(id_clusters_sampled))
  
  survey_data$sampling_fraction_2 = survey_data$ind_per_cluster / survey_data$cluster_population
  
  survey_data$overall_weight = 1 / (survey_data$sampling_fraction_1 * survey_data$sampling_fraction_2)
  
  return(survey_data)
}




Bootstrap_prevalence <- function(N_iterations,
                                 cluster_number,
                                 cluster_size,
                                 num_cluster_sample,
                                 ind_per_cluster,
                                 overall_prevalence1,
                                 overall_prevalence2,
                                 sigma_prevalence
                                 ){
#the function simulates the survey N_iterations times  to get a distribution of 
# Prevalence_t1 and Prevalence_t2 
Prevalence_t1 = as.vector(rep(NA, N_iterations))
Prevalence_t2 = as.vector(rep(NA, N_iterations))

for (ii in 1:N_iterations){
  
  survey_data <- Survey_pps(cluster_number = cluster_number,
                            cluster_size = cluster_size, 
                            num_cluster_sample = num_cluster_sample,
                            ind_per_cluster = ind_per_cluster, 
                            overall_prevalence1 = overall_prevalence1,
                            overall_prevalence2 = overall_prevalence2,
                            sigma_prevalence = sigma_prevalence)
  
  Prevalence_t1[ii]  = mean(survey_data$cluster_prevalence_t1)
  Prevalence_t2[ii]  = mean(survey_data$cluster_prevalence_t2)
}
return(data.frame( N_iterations = N_iterations,
                   Prevalence_t1 = Prevalence_t1,
                   Prevalence_t2 = Prevalence_t2
                  ))
}








library(survey)

suvery_object <- svydesign(id = ~cluster_id, 
                           data = survey_data, 
                           weights = ~overall_weight)

Prevalence_t1 <- svymean(~cluster_prevalence_t1, suvery_object, deff = T)
Prevalence_t2 <- svymean(~cluster_prevalence_t2, suvery_object, deff = T)


prevalence = c(Prevalence_t1 , Prevalence_t2)


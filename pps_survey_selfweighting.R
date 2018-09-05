#written by laurette Mhlanga an adaptation on the sampling with pps 
#psu to achieve equal weights in a 
#for every individual 



partition_prevalence <- function(overall_prevalence1,
                                 overall_prevalence2,
                                 Standard_error_prevalence,
                                 cluster_number,
                                 cluster_size
  #function that divides/partitions  the overall prevalence and assigns to clusters
  #into given clusters                                  
                                 
){
  
  CoV_overall =  Standard_error_prevalence/overall_prevalence1
  
  cluster_prevalences = rnorm(length(cluster_number), 
                              mean = overall_prevalence1, 
                              sd = Standard_error_prevalence)
  
  cluster_CoV = sd(cluster_prevalences) / mean(cluster_prevalences)
  
  cluster_prevalences_R_Adj = mean( cluster_prevalences ) + (cluster_prevalences - mean(cluster_prevalences )) * (CoV_overall / cluster_CoV)
  
  
  Cluster_factor = sum(cluster_size * cluster_prevalences_R_Adj) / sum(cluster_size)
  cluster_prevalences_final1 = (cluster_prevalences_R_Adj / Cluster_factor ) * overall_prevalence1
  
  cluster_prevalences_final2 = (cluster_prevalences_final1)*(overall_prevalence2 / overall_prevalence1)
  
  prevcheck =c(sum(cluster_size * cluster_prevalences_final1) / sum(cluster_size), sum(cluster_size * cluster_prevalences_final2) / sum(cluster_size))
  
  return(list(prevcheck, data.frame(cluster_number = cluster_number, 
                                    cluster_prevalence_t1 = cluster_prevalences_final1, 
                                    cluster_prevalence_t2 = cluster_prevalences_final2)))
  
}



prevalence = partition_prevalence(overall_prevalence1 = 0.2,
                                  overall_prevalence2 = 0.3 ,
                                  Standard_error_prevalence = 0.01,
                                  cluster_size = c(1028, 555, 390, 1309, 698, 907,
                                                   432, 897, 677, 501, 867, 867, 
                                                   1002, 1094, 668, 500, 835, 
                                                   396, 630, 483, 319, 569, 987, 598, 
                                                   375, 387, 465, 751, 365, 448),
                                  cluster_number = 1:30)[[2]]



Survey_pps <- function(cluster_number,
                       cluster_size , 
                       num_cluster_sample,
                       ind_per_cluster, 
                       overall_prevalence1,
                       overall_prevalence2,
                       Standard_error_prevalence
                      ){
  
  #Steps in appling Probability propotional to size(pps) and ways to ensure equal sampling weights 
  #individual (calculation of basic weights)
  
  survey_data = data.frame(cluster_id = cluster_number,
                           cluster_population = cluster_size)
  
  survey_data$cumulative_sum = cumsum(survey_data$cluster_population)
  
  
  survey_data = cbind(survey_data, partition_prevalence(overall_prevalence1 = overall_prevalence1,
                                          overall_prevalence2 = overall_prevalence2,
                                          Standard_error_prevalence = Standard_error_prevalence,
                                          cluster_number = cluster_number,
                                          cluster_size = cluster_size)[[2]][, -1])
  
  
  sampling_interval = survey_data$cumulative_sum[length(cluster_number)] / num_cluster_sample
  
  random_number = sample(1 : sampling_interval, 1)
  
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
  
survey_data <- Survey_pps(cluster_number = 1:30,
            cluster_size = c(1028, 555, 390, 1309, 698, 907,
                             432, 897, 677, 501, 867, 867, 
                             1002, 1094, 668, 500, 835, 
                             396, 630, 483, 319, 569, 987, 598, 
                             375, 387, 465, 751, 365, 448), 
            num_cluster_sample = 10,
            ind_per_cluster = 300, 
            overall_prevalence1 = 0.2,
            overall_prevalence2 = 0.3,
            Standard_error_prevalence = 0.01)




library(survey)

suvery_object <- svydesign(id = ~cluster_id, 
                           data = survey_data, 
                           weights = ~overall_weight)

Prevalence_t1 <- svymean(~cluster_prevalence_t1, suvery_object, deff = T)
Prevalence_t2 <- svymean(~cluster_prevalence_t2, suvery_object, deff = T)


prevalence = c(Prevalence_t1 , Prevalence_t2)


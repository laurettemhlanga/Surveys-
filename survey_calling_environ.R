


prevalence = partition_prevalence(overall_prevalence1 = 0.2,
                                  overall_prevalence2 = 0.3 ,
                                  sigma_prevalence = 0.01,
                                  cluster_size = c(1028, 555, 390, 1309, 698, 907,
                                                   432, 897, 677, 501, 867, 867, 
                                                   1002, 1094, 668, 500, 835, 
                                                   396, 630, 483, 319, 569, 987, 598, 
                                                   375, 387, 465, 751, 365, 448),
                                  cluster_number = 1:30)[[2]]





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
                          sigma_prevalence = 0.01)





Bootstrapp_Prevalence = Bootstrap_prevalence( N_iterations = 10000,
                                           cluster_number = 1:30,
                                           cluster_size = c(1028, 555, 390, 1309, 698, 907,
                                                            432, 897, 677, 501, 867, 867, 
                                                            1002, 1094, 668, 500, 835, 
                                                            396, 630, 483, 319, 569, 987, 598, 
                                                            375, 387, 465, 751, 365, 448), 
                                           num_cluster_sample = 10,
                                           ind_per_cluster = 300, 
                                           overall_prevalence1 = 0.2,
                                           overall_prevalence2 = 0.3,
                                           sigma_prevalence = 0.01)




library(survey)

suvery_object <- svydesign(id = ~cluster_id, 
                           data = survey_data, 
                           weights = ~overall_weight)

boot = as.svrepdesign(suvery_object, type = "JK1")

Prevalence <- svymean(~cluster_prevalence_t1+cluster_prevalence_t2, suvery_object, boot)





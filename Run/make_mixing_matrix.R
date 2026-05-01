
## script to generate the mixing matrix and aggregated Australian population

pop <- read.csv("Data/australian_population.csv")
pop_agg <- data.frame(age = seq(0, 80, by = 5), population = NA_integer_)
for(i in seq(0, 80, by = 5)) pop_agg[pop_agg$age == i,]$population <- sum(pop[pop$age >= i & pop$age < ifelse(i < 80, i + 5, 101),]$population)
mixing <- RSVModels::get_contact_matrix(pop = pop_agg)

write.csv(pop_agg, "Data/australian_population_aggregated.csv", row.names = FALSE)
write.csv(mixing, "Data/mixing.csv", row.names = FALSE)

rm(pop, pop_agg, i, mixing)

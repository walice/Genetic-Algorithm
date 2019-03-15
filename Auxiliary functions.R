plot.pop.evol <- function(result, filter = F, gen.sequence = NULL){
  pop.size <- result$pop.size
  pop.evolution <- as.data.frame(result$pop.evolution)
  pop.evolution <- gather(pop.evolution,
                          Gen, Individual)
  pop.evolution$Generation <- sapply(strsplit(pop.evolution$Gen, "_"), "[", 1)
  pop.evolution$Variable <- sapply(strsplit(pop.evolution$Gen, "_"), "[", 2)
  pop.evolution <- pop.evolution %>% 
    dplyr::select(Generation, Variable, Individual)
  
  if(filter == TRUE){
    plot.data <- pop.evolution %>%
      filter(Generation %in% gen.sequence)
  } else {
    plot.data <- pop.evolution
  }
  
  iter.num <- length(unique(plot.data$Generation))
  nvars <- length(unique(plot.data$Variable))
  
  if (nvars > 1){
    V1 <- plot.data %>%
      filter(Variable == "V1") %>%
      rename(Var1 = Individual) %>%
      dplyr::select(-Variable)
    V2 <- plot.data %>%
      filter(Variable == "V2") %>%
      rename(Var2 = Individual) %>%
      dplyr::select(-Variable)
    plot.data <- cbind(V1, V2 %>% dplyr::select(Var2))
  }
  
  if (nvars == 1){
    g <- ggplot(plot.data,
                aes(x = 1:(pop.size*iter.num),
                    y = Individual,
                    col = Generation))
  } else {
    g <- ggplot(plot.data,
                aes(x = Var1, y = Var2, col = Generation))
  }
  g <- g +
    geom_point() +
    labs(x = ifelse(nvars == 1, "", "Variable 1"),
         y = ifelse(nvars == 1, "Solution", "Variable 2"),
         title = "Evolution of solution in population") +
    scale_color_discrete(breaks = sort(as.numeric(plot.data$Generation))) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  if(nvars > 1){
    g <- g + geom_count(show.legend = F)
  }
  g  
}

plot.fitness.evol <- function(result, filter = F, gen.sequence = NULL){
  pop.size <- result$pop.size
  fitness.evolution <- as.data.frame(result$fitness.evolution)
  fitness.evolution <- gather(fitness.evolution,
                              Generation, Individual)
  
  if(filter == TRUE){
    plot.data <- fitness.evolution %>%
      filter(Generation %in% gen.sequence)
  } else {
    plot.data <- fitness.evolution
  }
  
  iter.num <- length(unique(plot.data$Generation))
  
  ggplot(plot.data, 
         aes(x = 1:(pop.size*iter.num), y = Individual, 
             col = Generation)) +
    geom_point() +
    labs(y = "Fitness",
         title = "Evolution of fitness in population") +
    scale_color_discrete(breaks = sort(as.numeric(plot.data$Generation))) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}

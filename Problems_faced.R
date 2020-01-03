## Problems we have faced up
library(OncoSimulR)
library(ggplot2)

## In this script we are showing one interesting problem we had 
## with our plot functions. This problem motivated us to test
## our functions (please go to "test_functions.R").

## Firstly we will show how low mutation rate (mu) affects the 
## results from the simulations, due to not all genotypes are 
## generated and "compositionPop2" can not work properly

################# "compositionPop2" can not work properly #################

## Plot box plot
simul_boxplot2 <- function(df, main = FALSE, xlab = "Genotype", ylab = "N") {
  ## Create box plot, title and axis parameters
  e <- ggplot(df, aes(x = Genotype, y = N)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11)
    )
  ## No title
  if (main ==FALSE) {
    e + geom_boxplot(aes(fill = Genotype)) +
      stat_summary(fun.y = mean, geom = "point",
                   shape = 18, size = 2.5, color = "#FC4E07") +
      xlab(xlab) + ylab(ylab)
  }
  ## Title
  else{
    e + geom_boxplot(aes(fill = Genotype)) +
      stat_summary(fun.y = mean, geom = "point",
                   shape = 18, size = 2.5, color = "#FC4E07") +
      labs(title = main) +
      xlab(xlab) + ylab(ylab)
  }
}

## Extract and create a data frame with results from several simulations
compositionPop2 <- function(objPop, ...) {
  clon_labels <- c("WT", objPop[[1]]$geneNames)
  listPop <- vapply(objPop, function(x) tail(x[[1]], 1)[1, -1], as.double(1:length(clon_labels)))
  dfPop <- data.frame("Genotype" = rep(clon_labels, length(listPop)/length(clon_labels)),
                      "N" = c(listPop))
  simul_boxplot2(dfPop, ...)
}

avc <- function (a, v, c) {
  data.frame(Genotype = c("WT", "G", "V", "A"),
             Fitness = c("1",
                         paste0("1 + ", a, " * (f_1 + 1)"),
                         paste0("1 + ", a, " * f_1 + ", v, " * (f_2 + 1) - ", c),
                         paste0("1 + ", a, " * f_1 + ", v, " * f_2")))
}

afavc <- allFitnessEffects(genotFitness = avc(37.5, 2, 1), # Parametros de la interaccion
                           frequencyDependentFitness = TRUE,
                           frequencyType = "rel")

simulation <- oncoSimulPop(5,
                           mc.cores = 6,
                           afavc,
                           model = "McFL",
                           onlyCancer = FALSE,
                           finalTime = 25,
                           mu = 1e-4,
                           initSize = 4000,
                           keepPhylog = TRUE,
                           seed = NULL,
                           errorHitMaxTries = FALSE,
                           errorHitWallTime = FALSE)

plot(simulation[[1]], show = "genotypes", type = "line")

compositionPop2(simulation)

## Error in vapply(objPop, function(x) tail(x[[1]], 1)[1, -1], as.double(1:length(clon_labels))) : 
## Los valores deben ser de longitud 4, 
## pero el resultado FUN(X [[3]]) es la longitud 2 

## Secondly we show an alternative fucntion to build up the box plot 
## with the results from different simulations to see agreement among them. 
## The alternative function is "compositionPop". This function does not break
## if listPop is a list instead of a matrix, and is able to pass a data frame
## to "simul_boxplot2" to create the box plot. 
## Nevertheless, as you will see, box plot make no sense. They show different
## results from the trajectories plot. 

################# "compositionPop" shows wrong results #################

compositionPop <- function(objPop) {
  condi <- c("WT", objPop[[1]]$geneNames)
  listPop <- lapply(objPop, function(x) tail(x[[1]], 1)[1, -1])
  dfPop <- data.frame(matrix(unlist(listPop), ncol = length(condi), byrow = TRUE))
  colnames(dfPop) <- condi
  simul_boxplot2(dfPop)
}


avc <- function (a, v, c) {
  data.frame(Genotype = c("WT", "G", "V", "A"),
             Fitness = c("1",
                         paste0("1 + ", a, " * (f_1 + 1)"),
                         paste0("1 + ", a, " * f_1 + ", v, " * (f_2 + 1) - ", c),
                         paste0("1 + ", a, " * f_1 + ", v, " * f_2")))
}

afavc <- allFitnessEffects(genotFitness = avc(37.5, 2, 1),
                           frequencyDependentFitness = TRUE,
                           frequencyType = "rel")

simulation <- oncoSimulPop(10,
                           mc.cores = 6,
                           afavc,
                           model = "McFL",
                           onlyCancer = FALSE,
                           finalTime = 25,
                           mu = 1e-4,
                           initSize = 4000,
                           keepPhylog = TRUE,
                           seed = NULL,
                           errorHitMaxTries = FALSE,
                           errorHitWallTime = FALSE)


## trajectorie plot
plot(simulation, show = "genotypes", type = "line")

## As you can chek from the plots above, "G" cell type
## achieves the maximum size. Moreover, sometimes is the 
## only mutant genotype. Now try box plot. 


## Boxplot
compositionPop(simulation)

## It makes no sense. First there is a warning message telling
## us that the data has different lengths. 
## Secondly boxplot gives completely different results, showing
## population sizes for mutant genotypes that was not present
## at the trajectories plots. 


## Why? Check the "test_functions.R". But also, let's check 
## the value of listPop and then the value of  dfPop
condi <- c("WT", simulation[[1]]$geneNames)
listPop <- lapply(simulation, function(x) tail(x[[1]], 1)[1, -1])
listPop
## As we check in "test_functions.R" not all genotypes have results in each 
## simulation (each vagon from listPop). listPop is not rectangular. 

## See dfPop
dfPop <- data.frame(matrix(unlist(listPop), ncol = length(condi), byrow = TRUE))
colnames(dfPop) <- condi
dfPop
## Awful! Values from different simulations are mixed.
## That is the reason why box plot are wrong. 

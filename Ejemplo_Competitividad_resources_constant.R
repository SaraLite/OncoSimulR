## Competitive situations considering resources as a constant
## r ->  intrinsic growth rates
## K -> carrying capacities (max capacity of the environment)
## a_12 -> competition coefficient. The effect of specie 2 on specie 1

## The signs of a_12 and a_21 therefore reflect the relationship
## between the two species:
##  a_12 a_21 Relationship
##    -   -   Mutualistic
##    -   0   Commensal
##    0   -   Commensal
##    +   -   Parasitic
##    -   +   Parasitic
##    +   +   Competitive
##    +   0   Amensalism
##    0   +   Amensalism

## We must introduce a WT specie in the model due to "to_genotFitness_std"
## which requires WT to appears once.

# https://books.google.es/books?id=HgAgMwsYOtkC&pg=PA167&lpg=PA167&dq=a21+a12+Lotka&source=bl&ots=SSTfj46eM6&sig=ACfU3U1aG5qdkQT56tmSrO0SC7SROkKVeQ&hl=es&sa=X&ved=2ahUKEwivovuP-9PmAhWRyoUKHRslAREQ6AEwA3oECAcQAQ#v=onepage&q=a21%20a12%20Lotka&f=false

## Hay que modificar el codigo para ponerle un contexto biologico
## pero lo que se va a ver aqui es la importancia de los parametros
## de competitividad en la dinamica de dos poblaciones

## Hay 4 pruebas:
## Prueba 1: No hay competitividad, los coeficientes de competitividad = 0
## Vemos que las poblaciones se acercan a su K (maxima poblacion soportada)
## sp2 tiene mas individuos que sp1

## Prueba 2: Ligera competitividad. Coef  a_12 = 0.25 a_21= 0.5
## Vemos que la sp 1 tiene mas efecto sobre sp2 y sp2 ahora tiene menos 
## individuos que sp1

## Prueba 3: Mayor competitividad, pero siguen pudiendo coexistir
## se comprueba viendo que el producto de los coef es menor que 1.
## a_12 * a_21 = 0.75 * 1 < 1. 
## La competencia ha aumentando y vemos que hay fluctuacion de laqs poblaciones
## Ninguna llega a la K

## Prueba 4: Competitividad muy alta, insostenible. a_12 = 0.9 a_21= 1.4
## La especie 1 afecta mucho a especie 2. La especie 2 se extingue y
## la especie 1, libre de competencia, llega a la K1. 

## Cosas que mejorar:
## 1- Graficas: sacar la leyenda y ponerla a mano
## 2- Intentar poner mas de 2 especies (una movida matematica)
## 3- Podemos cambiar los parametros para que sea mas exagerado el efecto


###############################################################################

## Lotka-Volterra expressions for competition models
## ( continuous-time differential equations)
G_fe_LV <- function(r1, r2, K1, K2, a_12, a_21, awt = 1,
                    gt = c("WT", "S1", "S2")) {
  data.frame(Genotype = gt,
             Fitness = c(
               paste0("max(0.1, 1 - ", awt, " * (n_2 + n_1))"),
               paste0("1 + ", r1,
                      " * ( 1 - (n_1 + ", a_12, " * n_2)/", K1,
                      ")"),
               paste0("1 + ", r2,
                      " * ( 1 - (n_2 + ", a_21, " * n_1)/", K2,
                      ")")
             ))
}

## Show expressions for birth rates
G_fe_LV("r1", "r2", "K1", "K2", "a_12", "a_21", "awt")


## Remember the above are the birth rates for a model with death rate = 1 (the "Exp"
## model). That is why we added a 1 to the birth rate. If you subtract 1 from the
## expressions for the birth rates of predators and prey above you get the standard
## expressions for the differential equations for Lotka-Volterra model of competition

## Test 1: No competition effects
## K2 > K1; r2 > r1 (specie 2 has greater capacity and growth rate)
## a21 > a12 (specie 1 has more effect on specie 2)
fe_competition <-
  allFitnessEffects(
    genotFitness =
      G_fe_LV(1.2, 1.7, 350, 400, 0, 0,
              gt = c("WT", "S1", "S2")),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs")

competition <- oncoSimulIndiv(fe_competition,
                              model = "Exp",
                              onlyCancer = FALSE,
                              finalTime = 500,
                              mu = 1e-1,
                              initSize = 4000,
                              keepPhylog = TRUE,
                              seed = NULL,
                              errorHitMaxTries = FALSE,
                              errorHitWallTime = FALSE)
competition

plot(competition, show = "genotypes", type = 'line',
     col = c("black", "green", "red", "blue"),
     xlim = c(20, 500))

## Test 2: Competition with Stable Coexistance
## K2 > K1; r2 > r1 (specie 2 has greater capacity and growth rate)
## a21 > a12 (specie 1 has more effect on specie 2)
## Vemos que ahora la competencia provoca que se lleguen a menor n? de indv
fe_competition <-
  allFitnessEffects(
    genotFitness =
      G_fe_LV(1.2, 1.7, 350, 400, 0.25, 0.5,
              gt = c("WT", "S1", "S2")),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs")

competition <- oncoSimulIndiv(fe_competition,
                              model = "Exp",
                              onlyCancer = FALSE,
                              finalTime = 500,
                              mu = 1e-1,
                              initSize = 4000,
                              keepPhylog = TRUE,
                              seed = NULL,
                              errorHitMaxTries = FALSE,
                              errorHitWallTime = FALSE)
competition

plot(competition, show = "genotypes", type = 'line',
     col = c("black", "green", "red", "blue"),
     xlim = c(20, 500),
     ylab = "Number of organisms")


## Test 3: Estable competition but more aggressive
# Vemos que las dos especies no son compatibles y el estado no es estable
## Ambas impactan bastante sobre la otra, pero pueden coexitir
## Se alejan del K original
fe_competition <-
  allFitnessEffects(
    genotFitness =
      G_fe_LV(1.2, 1.7, 350, 400, 0.75, 1,
              gt = c("WT", "S1", "S2")),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs")

competition <- oncoSimulIndiv(fe_competition,
                              model = "Exp",
                              onlyCancer = FALSE,
                              finalTime = 500,
                              mu = 1e-1,
                              initSize = 4000,
                              keepPhylog = TRUE,
                              seed = NULL,
                              errorHitMaxTries = FALSE,
                              errorHitWallTime = FALSE)
competition

plot(competition, show = "genotypes", type = 'line',
     col = c("black", "green", "red", "blue"),
     xlim = c(10, 500),
     ylab = "Number of organisms")


## Test 4: Unestable competition
# Vemos que las dos especies no son compatibles y el estado no es estable
## Una de ellas se extingue, y la otra llega a un valor cercano a la K.
## Sobrevive S1 que es la que m?s impacto tiene sobre la otra (a_21 > a_12)
fe_competition <-
  allFitnessEffects(
    genotFitness =
      G_fe_LV(1.2, 1.7, 350, 400, 0.9, 1.4,
              gt = c("WT", "S1", "S2")),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs")

competition <- oncoSimulIndiv(fe_competition,
                              model = "Exp",
                              onlyCancer = FALSE,
                              finalTime = 500,
                              mu = 1e-1,
                              initSize = 4000,
                              keepPhylog = TRUE,
                              seed = NULL,
                              errorHitMaxTries = FALSE,
                              errorHitWallTime = FALSE)
competition

plot(competition, show = "genotypes", type = 'line',
     col = c("black", "green", "red", "blue"),
     xlim = c(10, 500),
     ylab = "Number of organisms")

###############################################################################

## Now an example of "Amensalism" 

fe_ammensalism <-
  allFitnessEffects(
    genotFitness =
      G_fe_LV(1.2, 1.7, 350, 400, 0, 0,
              gt = c("WT", "S1", "S2")),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs")

ammensalism <- oncoSimulIndiv(fe_ammensalism,
                              model = "Exp",
                              onlyCancer = FALSE,
                              finalTime = 500,
                              mu = 1e-1,
                              initSize = 4000,
                              keepPhylog = TRUE,
                              seed = NULL,
                              errorHitMaxTries = FALSE,
                              errorHitWallTime = FALSE)
competition

plot(competition, show = "genotypes", type = 'line',
     col = c("black", "green", "red", "blue"),
     xlim = c(20, 500))
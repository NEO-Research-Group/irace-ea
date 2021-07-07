DifferentialEvolution <- R6::R6Class("DifferentialEvolution", cloneable=FALSE,
  public = list(
    type = "DE",
    num_parents = 4,
    CR = NULL,
    Fscale = NULL,

    initialize = function(CR, Fscale) {
      self$CR <- CR
      self$Fscale <- Fscale
    },
    next_children = function(...) invisible(),
      
    variation = function(parents, domain = NULL) {
      irace.assert(length(parents) == self$num_parents)
      X_base <- parents[1]
      X_r1 <- parents[2]
      X_r2 <- parents[3]
      X_r3 <- parents[4]
      value <- DE_best_1_bin(X_base, X_r1, X_r2, X_r3, self$CR, self$Fscale)
      return(value)
    }
))
# polynomial mutation proposed by Agrawal and Deb
# https://www.iitk.ac.in/kangal/papers/k2012016.pd
# distributionIndex is eta
mutation_polynomial <- function(y, yl, yu, distributionIndex)
{
  if (yl == yu) return(yl)

  irace.assert(yl < yu)
  distributionIndex <- distributionIndex + 1.

  delta1 <- (y - yl) / (yu - yl)
  delta2 <- (yu - y) / (yu - yl)
  rnd <- runif(1)
  mutPow <- 1. / distributionIndex
  if (rnd < 0.5) {
    xy <- 1. - delta1
    val <- 2. * rnd + (1. - 2. * rnd) * (xy^distributionIndex)
    deltaq <- (val^mutPow) - 1.
  } else {
    xy <- 1. - delta2
    val <- 2. * (1. - rnd) + 2. * (rnd - 0.5) * (xy^distributionIndex)
    deltaq <- 1. - (val^mutPow)
  }
  y <- y + deltaq * (yu - yl)
  y <- min(max(y, yl), yu)
  return(y)
}

int2binary <- function(x, n_bits = 128) {
  irace.assert(is.wholenumber(x))
  x <- head(intToBits(as.integer(x)), n_bits)
  as.integer(x)
}

# Precompute for speed.
.powers_of_two <- 2L^(0L:63L)
binary2int <- function(x) {
  irace.assert(length(x) <= length(.powers_of_two))
  sum(.powers_of_two[1:length(x)] * x)
}

binary_bitflip <- function(x, prob) ifelse(runif(length(x)) < prob, 1L - x, x)

mutation_intbitflip <- function(x, lower, upper, prob)
{
  x <- x - lower # to avoid negative numbers
  n_bits <- ceiling(log2(1 + upper - lower))
  bitval <- int2binary(x, n_bits)
  newval <- binary_bitflip(bitval, prob / n_bits)
  x <- binary2int(newval)
  x <- min(x + lower, upper) # It could happen that becomes larger than the maximum.
  x
}

Mutation <- R6::R6Class("Mutation", cloneable=FALSE,
  public = list(
    probMut = 1.,
    n_variables = 1.,
    initialize = function(prob) {
      irace.assert(prob >= 0 && prob <= 1)
      self$probMut <- prob
    },
    next_children = function(n_variables) {
      self$n_variables <- n_variables
    }
))

PolynomialMutation <- R6::R6Class("PolynomialMutation", inherit = Mutation, cloneable=FALSE,
  public = list(
    distributionIndex = 1., # MANUEL: What is the default?
    initialize = function(probMut, distributionIndex) {
      self$distributionIndex <- distributionIndex
      super$initialize(probMut)
    },
    
    do = function(x, lower, upper) {
      if (runif(1) >= self$probMut / self$n_variables) return(x)
      mutation_polynomial(x, lower, upper, self$distributionIndex)
    }
))

IntegerBitFlip <- R6::R6Class("IntegerBitFlip", inherit = Mutation, cloneable=FALSE,
  public = list(
    do = function(x, lower, upper) mutation_intbitflip(x, lower, upper, self$probMut / self$n_variables)
))
      
get_mutation <- function(x, probMut, mutParam)
  switch(x,
         polynomial = PolynomialMutation$new(probMut, distributionIndex = mutParam[["distributionIndex"]]),
         bitflip = IntegerBitFlip$new(probMut),
         irace.error("Unknown mutation"))


crossover_sbx <- function(p1, p2, distributionIndex, domain)
{
  if (is.null(distributionIndex))
    distributionIndex <- 20

  distributionIndex <- distributionIndex + 1
  EPS <- 1.0e-14
  if (domain[1] >= domain[2]) return(p1)
  if (abs(p1 - p2) < EPS) return(p1)
  if (runif(1) >= 0.5) return(p1)

  if (p1 < p2) {
    y1 <- p1
    y2 <- p2
  } else {
    y1 <- p2
    y2 <- p1
  }
      
  lowerBound <- domain[1]
  upperBound <- domain[2]
      
  rand <- runif(1)

  sbx_betaq <- function(beta, eta_c_p1, rand01) {
    alpha <- 2. - (beta ^ (-eta_c_p1))
    if (rand01 < (1. / alpha))
      return((rand01 * alpha)^(1. / eta_c_p1))
    return ((1. / (2. - rand01 * alpha))^(1. / eta_c_p1))
  }

  if (runif(1) < 0.5) {
    beta <- 1. + (2. * (y1 - lowerBound) / (y2 - y1))
    betaq <- sbx_betaq(beta, distributionIndex, rand)
    c1 <- 0.5 * (y1 + y2 - betaq * (y2 - y1))
    c1 <- min(max(c1, lowerBound), upperBound)
    return(c1)
  } else {
    beta <- 1. + (2. * (upperBound - y2) / (y2 - y1))
    betaq <- sbx_betaq(beta, distributionIndex, rand)
    c2 <- 0.5 * (y1 + y2 + betaq * (y2 - y1))
    c2 <- min(max(c2, lowerBound), upperBound)
    return (c2)
  }
}

# Process ONLY one parameter
crossover_uniform <- function(p1, p2)
{
  if (runif(1) <= 0.5) return(p1) else return(p2)
}



Crossover <- R6::R6Class("Crossover", cloneable=FALSE,
  public = list(
    probCross = 1,
    do_crossover = TRUE,
    initialize = function(probCross) {
      irace.assert(probCross >= 0 && probCross <= 1)
      self$probCross <- probCross
    },
    next_children = function(...) {
      self$do_crossover <- (runif(1) <= self$probCross)
    },
    do = function(p1, p2, domain) {
      if (self$do_crossover) return(self$apply(p1, p2, domain))
      return(p1)
    }    
))


SBX <- R6::R6Class("SBX", inherit = Crossover, cloneable=FALSE,
 public = list(
   distributionIndex = 20,
   initialize = function(probCross, distributionIndex) {
     self$distributionIndex <- distributionIndex
     super$initialize(probCross)
   },
   apply = function(p1, p2, domain) crossover_sbx(p1,p2, self$distributionIndex, domain)
))

UniformCX <- R6::R6Class("UniformCX", inherit = Crossover, cloneable=FALSE,
 public = list(
   apply = function(p1, p2, ...) crossover_uniform(p1,p2)
))

NoCX <- R6::R6Class("NoCX", inherit = Crossover, cloneable=FALSE,
 public = list(
   next_children = function(...) {
      self$do_crossover <- FALSE
    },
    apply = function(p1, p2, ...) stop("Should never be executed")
))

TwoPoint <- R6::R6Class("TwoPoint", inherit = Crossover, cloneable=FALSE, 
  public = list(
    cut_points = c(-1,-1),
    idx = 0,
    
    initialize = function(probCross) {
      super$initialize(probCross)
    },

    next_children = function(n_variables, ...) {
      self$cut_points <- sort(runif(2) * n_variables)
      self$idx <- 0
      super$next_children(...)
    },
        
    apply = function(p1, p2, ...) {
      if (self$idx <= self$cut_points[1] || self$idx >= self$cut_points[2])
        value <- p1
      else
        value <- p2
      self$idx <- self$idx + 1
      return(value)
    }
))

get_crossover <- function(x, probCross, crossParam)
  switch(x,
         uniform = UniformCX$new(probCross),
         sbx = SBX$new(probCross, distributionIndex = crossParam[["distributionIndex"]]),
         "2point" = TwoPoint$new(probCross),
         none = NoCX$new(probCross),
         irace.error("Unknown crossover"))


GeneticAlgorithm <- R6::R6Class("EvolutionaryOperator", cloneable=FALSE,
  public = list(
    type = "GA",
    num_parents = 4,
    crossover = NULL,
    mutation = NULL,

    initialize = function(probCross, crossParam, probMut, mutParam, crossover, mutation) {
      self$crossover <- get_crossover(crossover, probCross, crossParam)
      self$mutation <- get_mutation(mutation, probMut, mutParam)
    },
    next_children = function(n_variables) {
      self$crossover$next_children(n_variables = n_variables)
      self$mutation$next_children(n_variables = n_variables)
    },
  
    variation = function(parents, domain) {
      irace.assert(length(parents) == self$num_parents)
      p1 <- parents[1]
      p2 <- parents[2]
      # Crossover
      newVal <- self$crossover$do(p1, p2, domain)
      # Mutation
      newVal <- self$mutation$do(newVal, domain[1], domain[2])
      return(newVal)
    }
))


# Evolutionary Operators --------------------------------------------------

DE_best_1_bin <- function(X_base, X_r1, X_r2, X_r3, CR, Fscale)
{
  if (runif(1) <= CR)
    return(X_r1 + Fscale * (X_r2 - X_r3)) # base reproduction
  else
    return(X_base)
}

jade <- function(X_base, X_r1, X_r2, X_best, CR, Fscale)
{
  if (runif(1) <= CR)
    return (X_base + Fscale*(X_best - X_base) + Fscale*(X_r1 - X_r2))
  else
    return(X_base)
}

get_parents <- function(parameters, eliteConfigurations, numParents,
                        digits, forbidden, repair)
{
  # Random select from elites
  indexEliteParent <- sample.int (n = nrow(eliteConfigurations), 
                                  size = min(numParents, nrow(eliteConfigurations)),
                                  prob = eliteConfigurations[[".WEIGHT."]])
  #parents <- eliteConfigurations[indexEliteParent, !names(eliteConfigurations) %in% c(".ID.", ".ALIVE.", ".RANK.", ".WEIGHT."), drop = FALSE]
  parents <- eliteConfigurations[indexEliteParent, !names(eliteConfigurations) %in% c(".ALIVE.", ".RANK.", ".WEIGHT."), drop = FALSE]
  
  # Check if it is necessary generate new parent solutions
  extraSol <- max(numParents - nrow(eliteConfigurations), 0)
  if (extraSol > 0) {
    # MANUEL: If we don't have enough parents we sampleUniform, but why not
    # simply best configurations found previously?
    # Generate new parents configurations
    newConfigurationsAux <- sampleUniform(parameters, extraSol, digits, forbidden, repair)
    newConfigurationsAux[[".ID."]] <- 0
    print("names")
    print(names(eliteConfigurations))
    print(names(newConfigurationsAux))
    print(names(parents))
    parents <- rbind(parents, newConfigurationsAux)
  }
  return(parents)
}


# Sampling process with evolutionary operators
ea_generate <- function (parameters, eliteConfigurations,
                         nbNewConfigurations, digits, 
                         forbidden = NULL, repair = NULL)
{
  if (nbNewConfigurations <= 0) {
    irace.error ("The number of configurations to generate appears to be negative or zero.")
  }
  namesParameters <- names(parameters$conditions)
  newConfigurations  <-
    as.data.frame(matrix(nrow = nbNewConfigurations,
                         ncol = length(namesParameters) + 1,
                         dimnames = list(NULL, c(namesParameters, ".PARENT."))
                         ))
  empty_configuration <- new_empty_configuration(parameters)
  
  for (idxConfiguration in seq_len(nbNewConfigurations)) {
    forbidden.retries <- 0
    .irace$ea_variation$next_children(parameters$nbVariable)
    while (forbidden.retries < 100) {
      # Get the parents for the configuration
      numParents <- .irace$ea_variation$num_parents
      parents <- get_parents(parameters, eliteConfigurations, numParents, digits, forbidden, repair)
      configuration <- empty_configuration
      print(paste("configuration1: ", length(configuration), length(namesParameters) + 1, sep =' '))
      padre <- parents[1, ".ID."]
      if(is.null(padre) || is.na(padre)) padre <- 0
      configuration[[".PARENT."]] <- padre
      print(paste("configuration1: ", length(configuration), length(namesParameters) + 1, sep =' '))
      print(configuration)
      
      # Sample a value for every parameter of the new configuration.
      for (p in seq_along(namesParameters)) {
        # FIXME: We must be careful because parameters$types does not
        # have the same order as parameters$conditions. Ideally, we
        # should fix this or make it impossible to confuse them.
        currentParameter <- namesParameters[p]
        currentType <- parameters$types[[currentParameter]]
        if (!conditionsSatisfied(parameters, configuration, currentParameter)) {
          # Some conditions are unsatisfied.
          # Should be useless, NA is (always?) assigned when matrix created
          newVal <- NA
          
        } else if (isFixed(currentParameter, parameters)) {
          # We don't even need to sample, there is only one possible value !
          newVal <- get.fixed.value (currentParameter, parameters)
          # The parameter is not a fixed and should be sampled
        } else if (currentType %in% c("i", "r")) {
          domain <- getDependentBound(parameters, currentParameter, configuration)
          # Dependent domains could be not available because of inactivity of parameters
          # on which they are depedent. In this case, the dependent parameter becomes 
          # not active and we return NA.
          if (anyNA(domain)) {
            newVal <- NA
          } else {
            newVal <- .irace$ea_variation$variation(parents = parents[, currentParameter],
                                                    domain = domain)
            newVal <- numeric_value_round(currentType, newVal,
                                          lowerBound = domain[1], upperBound = domain[2], digits)
          }
        } else {
          irace.internal.error("Unexpected parameter type in ea_generate")
        }
        configuration[[p]] <- newVal
      }
      #print(paste("configuration1: ", length(configuration), sep =''))
      configuration <- as.data.frame(configuration, stringsAsFactors = FALSE)
      if (!is.null(repair)) {
        configuration <- repair(configuration, parameters, digits)
      }
      if (is.null(forbidden)
          || nrow(checkForbidden(configuration, forbidden)) == 1) {
        #print(paste("configuration2: ", length(configuration), sep =''))
        #print(configuration)
        #configuration[[".PARENT."]] <- parents[1, ".ID."]
        newConfigurations[idxConfiguration,] <- configuration
        break
      }
      forbidden.retries <- forbidden.retries + 1
    }
    if (forbidden.retries >= 100) {
      irace.error("irace tried 100 times to sample from the model a configuration not forbidden without success, perhaps your constraints are too strict?")
    }
  }
  return (newConfigurations)
}


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
  next_children = function() {  },
  
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

GeneticAlgorithm <- R6::R6Class("EvolutionaryOperator", cloneable=FALSE,
public = list(
  type = "GA",
  crossParam = NULL,
  probCross = NULL,
  probMut = NULL,
  do_crossover = TRUE,
  mutParam = NULL,
  num_parents = 4,

  initialize = function(probCross, crossParam = 0.5, probMut, mutParam) {
      self$probCross <- probCross
      self$crossParam <- crossParam
      self$probMut <- probMut
      self$mutParam <- mutParam
  },
  next_children = function() {
    self$do_crossover <- (runif(1) <= self$probCross)
    
  },
  variation = function(parents, domain) {
    irace.assert(length(parents) == self$num_parents)
    p1 <- parents[1]
    p2 <- parents[2]
    # Crossover
    newVal <- if (docrossover) crossover_uniform(p1, p2, self$crossParam) else p1
    # Mutation
    newVal <- if (runif(1) <= self$probMut) mutation_integerPolinomial(newVal, domain[1], domain[2], self$mutParam) else newVal
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

# Process ONLY one parameter
crossover_uniform <- function(p1, p2, crossParam)
{
  if (runif(1) <= crossParam) return(p1) else return(p2)
}

mutation_integerPolinomial <- function(x, yl, yu, distributionIndex)
{
  y <- x
  if (yl == yu) {
    y = yl
  } else {
    delta1 <- (y - yl) / (yu - yl)
    delta2 <- (yu - y) / (yu - yl)
    rnd <- runif(1)
    mutPow <- 1.0 / (distributionIndex + 1.0)
    if (rnd <= 0.5) {
      xy <- 1.0 - delta1
      val <- 2.0 * rnd + (1.0 - 2.0 * rnd) * (xy^(distributionIndex + 1.0))
      deltaq = (val^mutPow) - 1.0
    } else {
      xy <- 1.0 - delta2
      val <- 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (xy^(distributionIndex + 1.0))
      deltaq <- 1.0 - (val^mutPow)
    }
    y <- y + deltaq * (yu - yl)
  }
  return(y)
}

get_parents <- function(parameters, eliteConfigurations, numParents,
                        digits, forbidden, repair)
{
  # Random select from elites
  indexEliteParent <- sample.int (n = nrow(eliteConfigurations), 
                                  size = min(numParents, nrow(eliteConfigurations)),
                                  prob = eliteConfigurations[[".WEIGHT."]])
  parents <- eliteConfigurations[indexEliteParent, , drop = FALSE]
  
  # Check if it is necessary generate new parent solutions
  extraSol <- max(numParents - nrow(eliteConfigurations), 0)
  if (extraSol > 0) {
    # MANUEL: If we don't have enough parents we sampleUniform, but why not
    # simply best configurations found previously?
    # Generate new parents configurations
    newConfigurationsAux <- sampleUniform(parameters, extraSol, digits, forbidden, repair)
    parents <- rbind(parents, newConfigurationsAux)
  }
  return(parents)
}


# Sampling process with evolutionary operators
ea_generate <- function (parameters, eliteConfigurations,
                         nbConfigurations, digits, 
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
    .irace$ea_variation$next_children()
    while (forbidden.retries < 100) {
      # Get the parents for the configuration
      numParents <- .irace$ea_variation$num_parents
      parents <- get_parents(parameters, eliteConfigurations, numParents, digits, forbidden, repair)
      configuration <- empty_configuration
      configuration[[".PARENT."]] <- parents[1, ".ID."]
      
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
      configuration <- as.data.frame(configuration, stringsAsFactors = FALSE)
      if (!is.null(repair)) {
        configuration <- repair(configuration, parameters, digits)
      }
      if (is.null(forbidden)
          || nrow(checkForbidden(configuration, forbidden)) == 1) {
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

#######################################
## GENERATE CONFIGURATIONS
#######################################

## When called with an unconditional parameter, it
## must return TRUE
conditionsSatisfied <- function (parameters, partialConfiguration, paramName)
{
  condition <- parameters$conditions[[paramName]]
  # If there is no condition, do not waste time evaluating it.
  if (isTRUE(condition)) return(TRUE)

  v <- eval(condition, as.list(partialConfiguration))
  # Return TRUE if TRUE, FALSE if FALSE or NA
  ## FIXME: If we byte-compile the condition, then we should incorporate the
  ## following into the condition directly. See readForbiddenFile.
  v <- !is.na(v) && v
  return(v)
}

new.empty.configuration <- function(parameters)
{
  newConfigurationsColnames <- c(names(parameters$conditions), ".PARENT.")
  return(setNames(as.list(rep(NA, length(newConfigurationsColnames))),
                  newConfigurationsColnames))
}

get.fixed.value <- function(param, parameters)
{
  value <- parameters$domain[[param]][1]
  type <- parameters$types[[param]]
  if (type == "i") {
    return (as.integer(value))
  } else if (type == "c" || type == "o") {
    return (value)
  } else {
    irace.assert (type == "r")
    return (as.double(value))
  }
}

### Uniform sampling for the initial generation
sampleUniform <- function (parameters, nbConfigurations, digits,
                           forbidden = NULL, repair = NULL)
{
  namesParameters <- names(parameters$conditions)
  newConfigurations  <-
    as.data.frame(matrix(nrow = nbConfigurations,
                         ncol = length(namesParameters) + 1,
                         dimnames = list(NULL, c(namesParameters, ".PARENT."))
                         ))
  empty.configuration <- new.empty.configuration(parameters)

  for (idxConfiguration in seq_len(nbConfigurations)) {
    forbidden.retries <- 0
    while (forbidden.retries < 100) {
      configuration <- empty.configuration
      for (p in seq_along(namesParameters)) {
        currentParameter <- namesParameters[p]
        if (!conditionsSatisfied(parameters, configuration, currentParameter)) {
          configuration[[p]] <- NA
          next
        }
        # FIXME: We must be careful because parameters$types does not have the
        # same order as namesParameters, because we sample in the order of the
        # conditions.
        currentType <- parameters$types[[currentParameter]]
        if (isFixed(currentParameter, parameters)) {
          # We don't even need to sample, there is only one possible value !
          newVal <- get.fixed.value (currentParameter, parameters)
          # The parameter is not a fixed and should be sampled          
        } else if (currentType %in% c("i", "r")) {
          newVal <- sample.unif(currentParameter, parameters, currentType, digits)
        } else {
          irace.assert(currentType %in% c("c","o"))
          possibleValues <- parameters$domain[[currentParameter]]
          newVal <- sample(possibleValues, 1)
        }
        configuration[[p]] <- newVal
      }
      configuration <- as.data.frame(configuration, stringsAsFactors=FALSE)
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

# To be called the first time before the second race (with indexIter =
# 2) Nb configurations is the number of configurations at the end
# included the elite ones obtained from the previous iteration
sampleModel <- function (parameters, eliteConfigurations, model,
                         nbNewConfigurations, digits, forbidden = NULL,
                         repair = NULL)
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
  empty.configuration <- new.empty.configuration(parameters)
  
  for (idxConfiguration in seq_len(nbNewConfigurations)) {
    forbidden.retries <- 0
    while (forbidden.retries < 100) {
      # Choose the elite which will be the parent.
      indexEliteParent <- sample.int (n = nrow(eliteConfigurations), size = 1,
                                      prob = eliteConfigurations[[".WEIGHT."]])
      eliteParent <- eliteConfigurations[indexEliteParent, ]
      idEliteParent <- eliteParent[[".ID."]]
      configuration <- empty.configuration
      configuration[[".PARENT."]] <- idEliteParent
      
      # Sample a value for every parameter of the new configuration.
      for (p in seq_along(namesParameters)) {
        # FIXME: We must be careful because parameters$types does not
        # have the same order as parameters$conditions. Ideally, we
        # should fix this or make it impossible to confuse them.
        currentParameter <- namesParameters[p]
        currentType <- parameters$types[[currentParameter]]
        if (!conditionsSatisfied(parameters, configuration, currentParameter)) {
          # Some conditions are unsatisfied.
          # Should be useless, NA is ?always? assigned when matrix created
          newVal <- NA
          
        } else if (isFixed(currentParameter, parameters)) {
          # We don't even need to sample, there is only one possible value !
          newVal <- get.fixed.value (currentParameter, parameters)
          # The parameter is not a fixed and should be sampled
        } else if (currentType %in% c("i", "r")) {
          mean <- as.numeric(eliteParent[currentParameter])
          # If there is not value we obtain it from the model
          if (is.na(mean)) mean <- model[[currentParameter]][[as.character(idEliteParent)]][2]
          if (is.na(mean)) {
            # The elite parent does not have any value for this parameter,
            # let's sample uniformly.
            newVal <- sample.unif(currentParameter, parameters, currentType, digits)
                                                                                             
          } else {
            stdDev <- model[[currentParameter]][[as.character(idEliteParent)]][1]
            newVal <- sample.norm(mean, stdDev, currentParameter, parameters, currentType, digits)
          }
        } else if (currentType == "o") {
          possibleValues <- paramDomain(currentParameter, parameters)  
          value <- eliteParent[currentParameter]
          
          if (is.na(value)) {
            # The elite parent does not have any value for this
            # parameter, let's sample uniformly
            ## FIXME: We should save the last used parameter in the model and use it here.
            newVal <- sample(possibleValues, 1)
          } else {
            # Find the position within the vector of possible
            # values to determine the equivalent integer.
            mean <- match(value, possibleValues) # Return index of value in array
            stdDev <- model[[currentParameter]][[as.character(idEliteParent)]]

            # Sample with truncated normal distribution as an integer.
            # See sample.norm() for an explanation.
            newValAsInt <- floor(rtnorm(1, mean + 0.5, stdDev, lower = 1,
                                        upper = length(possibleValues) + 1L))

            # The probability of this happening is very small, but it can happen.
            if (newValAsInt == length(possibleValues) + 1L)
              newValAsInt <- length(possibleValues)
            
            irace.assert(newValAsInt >= 1L && newValAsInt <= length(possibleValues))
            # Get back to categorical values, find the one corresponding to the
            # newVal
            newVal <- possibleValues[newValAsInt]
          } 
        } else if (currentType == "c") {
          # FIXME: Why is idEliteParent character?
          # FIXME: Why the model is <parameter><Parent>? It makes more sense to be <Parent><parameter>.
          probVector <- model[[currentParameter]][[as.character(idEliteParent)]]
          possibleValues <- paramDomain(currentParameter, parameters)
          newVal <- sample(x = possibleValues, size = 1, prob = probVector)
        } else {
          irace.internal.error("Unexpected condition in sampleModel")
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

transform.from.log <- function(x, transf, lowerBound, upperBound)
{
  trLower <- attr(transf, "lower") 
  trUpper <- attr(transf, "upper")
  x <- exp(trLower + (trUpper - trLower) * x)
  return(x)
}

transform.to.log <- function(x, transf, lowerBound, upperBound)
{
  trLower <- attr(transf, "lower") 
  trUpper <- attr(transf, "upper")
  return((log(x) - trLower)/(trUpper - trLower))
}
## How to sample integer values?
#
# The problem: If we have an integer with domain [1,3] and we sample a real value
# and round, then there are more chances of getting 2 than 1 or 3:
# [1, 1,5) -> 1
# [1.5, 2,5) -> 2
# [2.5, 3) -> 3
#
# The solution: Sample in [lowerbound, upperbound + 1], that is, [1, 4], then floor():
# [1, 2) -> 1
# [2, 3) -> 2
# [3, 4) -> 3
#
# Why floor() and not trunc()?
# Because trunc(-1.5) -> -1, while floor(-1.5) -> -2, so for a domain [-3,-1]:
#
# [-3, -2) -> -3
# [-2, -1) -> -2
# [-1, 0)  -> -1
#
# Issue 1: We can sample 4 (upperbound + 1). In that case, we return 3.
#
# Issue 2: When sampling from a truncated normal distribution, the extremes are
# not symmetric.
#
# nsamples <- 100000
# table(floor(rtnorm(nsamples, mean=1, sd=1, lower=1,upper=4)))/nsamples
# table(floor(rtnorm(nsamples, mean=3, sd=1, lower=1,upper=4)))/nsamples
#
# To make them symmetric, we translate by 0.5, so that the mean is at the
# actual center of the interval that will produce the same value after
# truncation, e.g., given an integer value of 1, then mean=1.5, which is at the
# center of [1,2).
#
# nsamples <- 100000
# table(floor(rtnorm(nsamples, mean=1.5, sd=1, lower=1,upper=4)))/nsamples
# table(floor(rtnorm(nsamples, mean=3.5, sd=1, lower=1,upper=4)))/nsamples
#
# The above reasoning also works for log-transformed domains, because 
# floor() happens in the original domain, not in the log-transformed one,
# except for the case of log-transformed negative domains, where we have to
# translate by -0.5.
# 
numeric.value.round <- function(type, value, lowerBound, upperBound, digits)
{  
  irace.assert(is.finite(value))
  if (type == "i") {
    value <- floor(value)
    upperBound <- upperBound - 1L # undo the above for the assert
    # The probability of this happening is very small, but it could happen.
    if (value == upperBound + 1L)
      value <- upperBound
  } else
    value <- round(value, digits)

  irace.assert(value >= lowerBound && value <= upperBound)
  return (value)
}

# Sample value for a numerical parameter.
sample.unif <- function(param, parameters, type, digits = NULL)
{
  lowerBound <- paramLowerBound(param, parameters)
  upperBound <- paramUpperBound(param, parameters)
  transf <- parameters$transform[[param]]
  if (type == "i") {
    # +1 for correct rounding before floor()
    upperBound <- 1L + upperBound
  }
  if (transf == "log") {
    value <- runif(1, min = 0, max = 1)
    value <- transform.from.log(value, transf, lowerBound, upperBound)
  } else {
    value <- runif(1, min = lowerBound, max = upperBound)    
  }
  value <- numeric.value.round(type, value, lowerBound, upperBound, digits)
  return(value)
}

sample.norm <- function(mean, sd, param, parameters, type, digits = NULL)
{
  lowerBound <- paramLowerBound(param, parameters)
  upperBound <- paramUpperBound(param, parameters)
  transf <- parameters$transform[[param]]
  if (type == "i") {
    upperBound <- 1L + upperBound
    # Because negative domains are log-transformed to positive domains.
    mean <- mean + 0.5
  }
  
  if (transf == "log") {
    trMean <- transform.to.log(mean, transf, lowerBound, upperBound)
    value <- rtnorm(1, trMean, sd, lower = 0, upper = 1)
    value <- transform.from.log(value, transf, lowerBound, upperBound)
  } else {
    value <- rtnorm(1, mean, sd, lowerBound, upperBound)
  }

  value <- numeric.value.round(type, value, lowerBound, upperBound, digits)
  return(value)
}


# Evolutionary Operators --------------------------------------------------

DE_best_1_bin <- function(X.base, X.r1, X.r2, X.r3, CR, Fscale) {
  child <- if (runif(1) <= CR)
    X.r1 + Fscale*(X.r2 - X.r3) # base reproduction
  else X.base #+ 0.5*(Fscale + 1)*(X.r1 + X.r2 - 2*X.base)
  return(child)
}

jade <- function(X.base, X.r1, X.r2, X.best, CR, Fscale) {
  if (runif(1) <= CR)
    return (X.base + Fscale*(X.best - X.base) + Fscale*(X.r1 - X.r2))
  else return(X.base)
}

# Process ONLY one parameter
crossover_uniform <- function(p1, p2, crossParam) {
  if (runif(1) <= crossParam) return(p1) else return(p2)
}

mutation_integerPolinomial <- function(x, yl, yu, distributionIndex) {
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
      val <- 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (xy ^(distributionIndex + 1.0))
      deltaq <- 1.0 - (val^mutPow)
    }
    y <- y + deltaq * (yu - yl)
  }
  return(y)
}

get_parents <- function(parameters, eliteConfigurations, numParents, configurations.colnames, digits, forbidden, repair) {
  # MANUEL: Since these are not really parameters, they should be ".CR."
  # and ".Fscale.", but do we really need to store them per configuration?
  parents <- as.data.frame(matrix(ncol = length(configurations.colnames),
                                  nrow = 0,
                                  dimnames = list(NULL, configurations.colnames)))
  # Random select from elites
  indexEliteParent <- sample.int (n = nrow(eliteConfigurations), 
                                  size = min(numParents, nrow(eliteConfigurations)),
                                  prob = eliteConfigurations[[".WEIGHT."]])
  drops <- c(".ALIVE.", ".RANK.", ".WEIGHT.")
  eliteParentRemoved <- eliteConfigurations[indexEliteParent, !(names(eliteConfigurations) %in% drops)]
  parents <- rbind(parents, eliteParentRemoved)
  
  # Check if it is necessary generate new parent solutions
  extraSol <- if (numParents > nrow(eliteConfigurations)) (numParents - nrow(eliteConfigurations)) else 0
  if (extraSol > 0) {
    # MANUEL: If we don't have enough parents we sampleUniform, but why not
    # simply best configurations found previously?
    # Generate new parents configurations
    newConfigurationsAux <- sampleUniform(parameters, extraSol, digits, forbidden, repair)
    newConfigurationsAux <-
      cbind (.ID. = max(0, parents$.ID.) + 1:nrow(newConfigurationsAux),
             newConfigurationsAux)
    parents <- rbind(parents, newConfigurationsAux)
  }
  return(parents)
}

# Sampling process with evolutionary operators
sampleModelEA <- function (parameters, eliteConfigurations, model,
                           nbNewConfigurations, digits, scenario,
                           forbidden = NULL, repair = NULL, solSet = NULL)
{
  # message(paste('sampleModelEA', operation, toString(eaparameters), sep = " "))
  if (nbNewConfigurations <= 0) {
    irace.error ("The number of configurations to generate appears to be negative or zero.")
  }
  namesParameters <- names(parameters$conditions)
  newConfigurations  <-
    as.data.frame(matrix(nrow = nbNewConfigurations,
                         ncol = length(namesParameters) + 3,#length(namesParameters) + 1,
                         dimnames = list(NULL, c(namesParameters, ".PARENT.", "CR", "Fscale"))
    ))
  empty.configuration <- new.empty.configuration(parameters)
  
  # Column names to configuration sets (e.g. parents)
  configurations.colnames <- c(".ID.", namesParameters, ".PARENT.", "CR", "Fscale")
  
  for (idxConfiguration in seq_len(nbNewConfigurations)) {
    
    forbidden.retries <- 0
    while (forbidden.retries < 100) {
      configuration <- empty.configuration
      
      idEliteParent <- NA 
      # Choose the elite which will be the parent.zz
      numParents <- switch(scenario$ealgorithm, "DE" = 4, "GA" = 2, "JADE" = 4)
      docrossover <- 1
      if (scenario$ealgorithm == "GA") {
        probCross <- scenario$probCross
        docrossover <- (runif(1) <= probCross)
        probMut <- scenario$probMut
        mutParam <- scenario$mutParam # distributionIndex <- 20
      }
      if (scenario$ealgorithm == "JADE") {
        Fscale <- rcauchy(1, location = scenario$muF, scale = 0.1)
        Fscale <- clamp(Fscale, lower = 0, upper = 1)
        configuration[["Fscale"]] <- Fscale
        CR <- rnorm(1, mean = scenario$muCR, sd = 0.1)
        CR <- clamp(CR, 0, 1)
        configuration[["CR"]] <- CR
      }
      
      # Get the parents for the configuration
      parents <- get_parents(parameters, eliteConfigurations, numParents, configurations.colnames, digits, forbidden, repair)
      # FIXME: Always is the first parent that it is assigned as .PARENT.
      configuration[[".PARENT."]] <- parents[1,][[".ID."]]
      
      # MANUEL: Cual es la diferencia entre parents y parentsAll?
      parentsAll <-
        as.data.frame(matrix(ncol = length(configurations.colnames),
                             nrow = 0,
                             dimnames = list(NULL, configurations.colnames)))
      if (scenario$ealgorithm == "JADE") {
        drops <- c(".ALIVE.", ".RANK.", ".WEIGHT.")
        eliteParentRemoved <- eliteConfigurations[ , !(names(eliteConfigurations) %in% drops)]
        parentsAll <- rbind(parents, eliteParentRemoved)
      }
      
      config <- eliteConfigurations
      
      if (scenario$ealgorithm == "JADE") {
        # Calculation of a set of configurations joining the current parents, elites, and archive
        # PuA = parents U elites U archive
        drops <- c(".ALIVE.", ".RANK.", ".WEIGHT.")
        eliteParents <- eliteConfigurations[ , !(names(eliteConfigurations) %in% drops)]
        PuA <- rbind(parents, eliteParents, solSet[[1]])
        
        indexPuA <- sample.int (n = nrow(PuA), size = 1)
        parents[2,] <- PuA[indexPuA, ] # We modify this parent based to the original article
        
        eliteParents <- eliteParents[order(".WEIGHT."),]
        # Sample only from the 50% best ones
        irace.assert(nrow(eliteParents) > 0)
        indexBest <- sample.int(n = ceiling(0.5 * nrow(eliteParents)), size = 1)
        parents[3,] <- eliteParents[nrow(eliteParents) - indexBest + 1,] # change the parent to one in the 50% of the best elites
      }
      # Sample a value for every parameter of the new configuration.
      for (p in seq_along(namesParameters)) {
        # FIXME: We must be careful because parameters$types does not
        # have the same order as parameters$conditions. Ideally, we
        # should fix this or make it impossible to confuse them.
        currentParameter <- namesParameters[p]
        currentType <- parameters$types[[currentParameter]]
        if (!conditionsSatisfied(parameters, configuration, currentParameter)) {
          # Some conditions are unsatisfied.
          # Should be useless, NA is ?always? assigned when matrix created
          newVal <- NA
          
        } else if (isFixed(currentParameter, parameters)) {
          # We don't even need to sample, there is only one possible value !
          newVal <- get.fixed.value (currentParameter, parameters)
          # The parameter is not a fixed and should be sampled
        } else if (currentType %in% c("i", "r")) {
          if (scenario$ealgorithm == "DE") {
            irace.assert(nrow(parents) == 4)
            X.base <- as.numeric(parents[4,][currentParameter])
            X.r1 <- as.numeric(parents[1,][currentParameter])
            X.r2 <- as.numeric(parents[2,][currentParameter])
            X.r3 <- as.numeric(parents[3,][currentParameter])
            
            CR <- scenario$deCR
            Fscale <- scenario$deF # default 0.5
            
            newVal <- DE_best_1_bin(X.base, X.r1, X.r2, X.r3, CR, Fscale)
            lowerBound <- paramLowerBound(currentParameter, parameters)
            upperBound <- paramUpperBound(currentParameter, parameters)
            newVal <- clamp(newVal, lower = lowerBound, upper = upperBound)
          } 
          if (scenario$ealgorithm == "JADE") {
            irace.assert(nrow(parents) == 4)
            X.base <- as.numeric(parents[4,][currentParameter])
            X.r1 <- as.numeric(parents[1,][currentParameter])
            X.r2 <- as.numeric(parents[2,][currentParameter])
            X.r3 <- as.numeric(parents[3,][currentParameter])
            
            newVal <- jade(X.base, X.r1, X.r2, X.r3, CR, Fscale) # X.r3=Xbest
            lowerBound <- paramLowerBound(currentParameter, parameters)
            upperBound <- paramUpperBound(currentParameter, parameters)
            newVal <- clamp(newVal, lower = lowerBound, upper = upperBound)
          } else if (scenario$ealgorithm == "GA") {
            irace.assert(nrow(parents) == 2)
            p1 <- as.numeric(parents[1,][currentParameter])
            p2 <- as.numeric(parents[2,][currentParameter])
            # Crossover
            crossParam <- scenario$crossParam# default 0.5
            newVal <- if (docrossover) crossover_uniform(p1, p2, crossParam) else p1
            
            # Mutation
            domutation <- (runif(1) <= probMut)
            newVal <- if (domutation) mutation_integerPolinomial(newVal, lowerBound, upperBound, mutParam) else newVal
            
            lowerBound <- paramLowerBound(currentParameter, parameters)
            upperBound <- paramUpperBound(currentParameter, parameters)
            newVal <- clamp(newVal, lowerBound, upperBound)
          }
          newVal <- ifelse(currentType == "i", round(newVal),
                           round(newVal, digits))
          
        } else if (currentType == "o") {
          # MANUEL: Deberíamos reemplazar este código con un irace.assert() porque no deberíamos llegar nunca aquí.
          possibleValues <- paramDomain(currentParameter, parameters)  
          value <- eliteParent[currentParameter]
          
          if (is.na(value)) {
            # The elite parent does not have any value for this
            # parameter, let's sample uniformly
            ## FIXME: We should save the last used parameter in the model and use it here.
            newVal <- sample(possibleValues, 1)
          } else {
            # Find the position within the vector of possible
            # values to determine the equivalent integer.
            mean <- match(value, possibleValues) # Return index of value in array
            stdDev <- model[[currentParameter]][[as.character(idEliteParent)]]
            
            # Sample with truncated normal distribution as an integer.
            # See sample.norm() for an explanation.
            newValAsInt <- floor(rtnorm(1, mean + 0.5, stdDev, lower = 1,
                                        upper = length(possibleValues) + 1L))
            # The probability of this happening is very small, but it can happen.
            if (newValAsInt == length(possibleValues) + 1L)
              newValAsInt <- length(possibleValues)
            
            irace.assert(newValAsInt >= 1L && newValAsInt <= length(possibleValues))
            # Get back to categorical values, find the one corresponding to the
            # newVal
            newVal <- possibleValues[newValAsInt]
          } 
        } else if (currentType == "c") {
          # MANUEL: Deberíamos reemplazar este código con un irace.assert() porque no deberíamos llegar nunca aquí.
          # FIXME: Why is idEliteParent character?
          # FIXME: Why the model is <parameter><Parent>? It makes more sense to be <Parent><parameter>.
          probVector <- model[[currentParameter]][[as.character(idEliteParent)]]
          possibleValues <- paramDomain(currentParameter, parameters)
          newVal <- sample(x = possibleValues, size = 1, prob = probVector)
        } else {
          irace.internal.error("Unexpected condition in sampleModelEA")
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
  # drops <- c(".ID.")
  # newConfigurations <- newConfigurations[ , !(names(newConfigurations) %in% drops)]
  return (newConfigurations)
}
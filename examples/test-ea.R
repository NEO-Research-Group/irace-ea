library(irace)

instances <- 1:20

parameters <- readParameters(text='
r1 "" r (-5, 5)
r2 "" r (-5, 5)
r3 "" r (-5, 5)
r4 "" r (-5, 5)
i1 "" i (0, 10)
i2 "" i (0, 10)
i3 "" i (0, 10)
i4 "" i (0, 10)
')

target.runner <- function(experiment, scenario)
{
  debugLevel    <- scenario$debugLevel
  configuration.id  <- experiment$id.configuration
  instance.id   <- experiment$id.instance
  seed          <- experiment$seed
  configuration <- experiment$configuration
  instance      <- experiment$instance

  value <- 0
  for (p in names(configuration)) {
    if (startsWith(p, "r")) {
      value <- value - configuration[[p]]
    } else {
      value <- value + configuration[[p]]
    } 
  }
  result <- list(cost = value, call = toString(experiment))
  return(result)
}

run_irace <- function(...)
{
  args <- list(...)
  scenario <- list(targetRunner = target.runner, instances = instances,
                   maxExperiments = 1000, seed = 1234567, debugLevel = 1)
  scenario <- modifyList(scenario, args)
  printScenario(scenario)
  scenario <- checkScenario (scenario)
  confs <- irace(scenario = scenario, parameters = parameters)
  return(confs)
}

run_irace(ea_variation = "GA")


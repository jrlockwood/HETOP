condition_get <- function(xicc, xcv, xng, xcutscores, xreps) {
  #  get the reps of the condition specified by icc, cv, ng and cutscores
  #  from the file hopsim1-counts.dta
  
  d <- read.dta(file = 'hopsim1-counts.dta')
  d <- subset(d, icc == xicc & cv == xcv & ng == xng & (rep %in% xreps))
  
  if (xcutscores == 'mid') {
    d <- within(d, {
      level0 = nb_20
      level1 = nb_50 - nb_20
      level2 = nb_80 - nb_50
      level3 = ng - nb_80
    }
    )
  } else if (xcutscores == 'wide') {
    d <- within(d, {
      level0 = nb_05
      level1 = nb_50 - nb_05
      level2 = nb_95 - nb_50
      level3 = ng - nb_95
    }
    )
  } else if (xcutscores == 'skewed') {
    d <- within(d, {
      level0 = nb_05
      level1 = nb_30 - nb_05
      level2 = nb_55 - nb_30
      level3 = ng - nb_55
    }
    )
  } else if (xcutscores == 'many') {
    d <- within(d, {
      level0 = nb_5
      level1 = nb_25 - nb_05
      level2 = nb_50 - nb_25
      level3 = nb_75 - nb_50
      level4 = nb_95 - nb_75
      level5 = ng - nb_95
    }
    )
  }
  
  if (xcutscores != 'many') { 
    d <- d[ , c('id', 'rep', 'cv', 'icc', 'ng', 'mtrue', 'strue', 'samp_mean', 'samp_sd', 'level0', 'level1', 'level2', 'level3') ]
  } else {
    d <- d[ , c('id', 'rep', 'cv', 'icc', 'ng', 'mtrue', 'strue', 'samp_mean', 'samp_sd', 'level0', 'level1', 'level2', 'level3', 'level4', 'level5') ]
  }
  
  return(d)
  
}

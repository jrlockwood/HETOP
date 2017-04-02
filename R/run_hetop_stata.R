###########################################################
## FUNCTION TO CALL OUR HETOP.ADO STATA ROUTINE
###########################################################

# NOTE: this requires:
#     hetop.ado, oglm.ado, oglminit.ado, oglm_ll.ado, oglm_p.ado, oglm_svy_check.ado, vgap.ado
#     to be in the working directory
#     oglm could be install via ssc install, but vgap and hetop cannot...
# the function will call stata via system() and return estimates in a similar format
# to those returned by our R MLE functions
# will include both the prime estimates (in sums constraints)
# and star estimates which are "bias-corrected"

run_hetop_stata <- function(ngk, modtype = "hetop", maxit = 1000, noisy = FALSE, adodir=NULL, workdir=NULL) {
  ## ngk is a GxK matrix of counts
  ## modtype is passed to hetop in stata and can be hetop or homop
  ## maxit is maximum number of iterations to try
  ## noisy controls whether stata log file is printed in .Rout file
  ## adodir is the directory with the necessary ado files and will be the
  ##    working directory for stata commands
  ## workdir is the directory for storing .do and .dta files along the way (usually tempdir passed in)
  
  if (is.null(adodir)) {
    stop("ado dir must not be NULL")
  }
  if (is.null(workdir)) {
    stop("workdir dir must not be NULL")
  }
  
  .K <- ncol(ngk)
  .G <- nrow(ngk)
  
  statadir <- workdir

  statainfile  <- paste0(tempfile("statain", tmpdir = statadir), ".dta")
  statadofile  <- paste0(statadir, "/statascript.do")
  stataoutfile <- paste0(tempfile("stataout", tmpdir = statadir), ".dta")
  
  if(file.exists(statainfile)) file.remove(statainfile)
  if(file.exists(statadofile)) file.remove(statadofile)
  if(file.exists(stataoutfile)) file.remove(stataoutfile)
  
  ngk <- data.frame(id = (1:G), ngk)
  write.dta(ngk, file = statainfile)
  
  # some pieces of the do-file
  usefilecall <- paste0("use \"", statainfile, "\", clear")
  hetopcall <- paste0("noi cap hetop id X , numcats(", .K, ") modtype(", modtype, ") identify(sums) iterate(", maxit, ")")
  savefilecall <- paste0("saveold \"", stataoutfile, "\" , replace")
  
  # create do-file
  cat("
      clear all
      set more off
      ", file = statadofile)
  cat(paste0("cd \"", adodir, "\""), file = statadofile, append=TRUE, sep="\n")
  cat("
      pwd
      which hetop
      which oglminit
      ", file = statadofile, append = TRUE)
  cat(paste0("local K=", .K), file = statadofile, append=TRUE, sep = "\n")
  cat(paste0("local G=", .G), file = statadofile, append=TRUE, sep = "\n")
  cat(usefilecall, file = statadofile, append = TRUE, sep = "\n")
  cat(hetopcall, file = statadofile, append = TRUE)
  cat("
      if _rc == 0 & e(converged) == 1 {
      foreach m in mstar sstar mprime sprime {
      mat `m' = e(`m')
      svmat `m'
      rename `m'1 `m'
      }
      g icc_bc = e(icchat)
      mat cutstar = e(cutsstar)
      mat cutprime = e(cutsprime)
      forv i = 1/`=`K'-1' {
      g cut`i'_bc = cutstar[`i',1]
      g cut`i'_mle = cutprime[`i',1]
      }
      g loglik = e(ll)
      g iters = e(ic)
      g converged = e(converged)
      g rc_code = 0
      g hetop_prob = 0
      noi di \"### hetop converged normally in \" e(ic) \" iterations ###\"
      }
      else {
      g mstar = 0
      g mprime = 0
      g sstar = 1
      g sprime = 1
      g icc_bc = 0
      forv i = 1/`=`K'-1' {
      g cut`i'_bc = 0
      g cut`i'_mle = 0
      }
      g loglik = e(ll)
      g iters = e(ic)
      g converged = e(converged)
      g rc_code = _rc
      g hetop_prob = 1
      noi di \"### hetop failed to converge normally after \" e(ic) \" iterations with error code \" _rc \" ###\"
      }
      compress
      ", file = statadofile, append = TRUE)
  cat(savefilecall, file = statadofile, append = TRUE, sep = "\n")
  
  # call stata
  system(paste0("stata -b do \"", statadofile, "\""))
  
  # print output and clean up log file
  if (noisy) print(readLines("statascript.log"))
  if (file.exists("statascript.log")) file.remove("statascript.log")
  
  # get and return results; check for problems
  statares <- read.dta(file = stataoutfile)
  
  if (statares$hetop_prob[1] == 1) {
    cat(paste("\n############# problem with stata hetop function #############\n"))
  }
  
  if(file.exists(statadofile)) file.remove(statadofile)
  if(file.exists(stataoutfile)) file.remove(stataoutfile)
  if(file.exists(statainfile)) file.remove(statainfile)
  
  return(list(mug = statares$mprime,
              sigmag = statares$sprime,
              mug_star = statares$mstar,
              sigmag_star = statares$sstar,
              icc_bc = statares$icc_bc[1],
              cutpoints_star = statares[1,c('cut1_bc', 'cut2_bc', 'cut3_bc')],
              cutpoints = statares[1,c('cut1_mle', 'cut2_mle', 'cut3_mle')],
              hetopdetails = list(iter = statares$iters[1],
                                  conv = statares$converged[1],
                                  rc_code = statares$rc_code[1],
                                  loglik = statares$loglik[1]))
  )
  
}

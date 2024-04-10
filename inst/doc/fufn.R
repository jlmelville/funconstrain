fufn <- function(fnum=NULL){
  # return list with tfn=function, tgr=gradient given fn number and n
  if (is.null(fnum)) stop("ffn needs a function number fnum")
  if ((fnum < 1) || (fnum > 35)) stop("fnum must be in [1, 35]")
#  cat("entering ffn, fnum=",fnum,"\n")
  # select function
  funnam <- c("rosen", "freud_roth", "powell_bs", "brown_bs", "beale", 
              "jenn_samp", "helical", "bard", "gauss", "meyer", "gulf",
              "box_3d", "powell_s", "wood", "kow_osb", "brown_den", 
              "osborne_1", "biggs_exp6", "osborne_2", "watson", "ex_rosen", 
              "ex_powell", "penalty_1", "penalty_2", "var_dim", "trigon", 
              "brown_al", "disc_bv", "disc_ie", "broyden_tri", "broyden_band", 
              "linfun_fr", "linfun_r1", "linfun_r1z", "chebyquad")
#  print(str(funnam))
  fname <- funnam[as.integer(fnum)]
#  cat("fname:", fname,"\n")
  while (fnum %in% 1:35) {
    ameth <- optimx::ctrldefault(2)$bdmeth # Choose only bounded methods
    ameth <- ameth[ameth != "lbfgsb3c"] ## ?? Temporarily remove lbfgsb3c
    ameth <- c(ameth, "L-BFGS-B")
    # ?? may want to test allmeth to check that inappropriate methods are captured
    #    cat("in while, fnum=", fnum); tmp <- readline("cont.")
    mm <- 0 # in case m value needed
  if (fnum == 1) {
     n <- 2 # fixed
     mm <- 2
     tt <- rosen()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 2) {
     n <- 2 # fixed
     mm <- 2
     tt <- freud_roth()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 3) {
     n <- 2 # fixed
     mm <- 2
     tt <- powell_bs()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 4) {
     n <- 2 # fixed
     mm <- 3
     tt <- brown_bs()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
## BAD -- reset 20240323
#     lo <- rep((min(xx0)-0.1), n)
#     up <- rep((max(xx0)+0.1), n)
     lo <- -1e20
     up <- -lo
     break }

  if (fnum == 5) {
     n <- 2 # fixed
     mm <- 3
     tt <- beale()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 6) {
     n <- 2 # fixed
     mm <- 10
     tt <- jenn_samp()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 7) {
     n <- 3 # fixed
     tt <- helical()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 8) {
     n <- 3 # fixed
     mm <- 15
     tt <- bard()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 9) {
     n <- 3 # fixed
     mm <- 15
     tt <- gauss()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 10) {
     n <- 3 # fixed
     m <- 16 # ?? how to return
     tt <- meyer()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 11) {
     n <- 3
     mm <- 99
     tt <- gulf()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 12) {
     n <- 3
     mm <- 20
     tt <- box_3d()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 13) {
     n <- 4
     tt <- powell_s()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 14) {
     n <- 4
     tt <- wood()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 15) {
     mm <- 11
     n <- 4
     tt <- kow_osb()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 16) {
     mm <- 20
     n <- 4
     tt <- brown_den()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 17) {
     mm <- 33
     n <- 5
     tt <- osborne_1()
     ameth<-ameth[-which(ameth=="L-BFGS-B")] # remove L-BFGS-B from this case
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     lo[4] <- 0
     lo[5] <- 0
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 18) {
     mm <- 20
     n <- 6
     tt <- biggs_exp6()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 19) {
     mm <- 65
     n <- 11
     tt <- osborne_2()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 20) {
     n <-8
     mm <- 31
     tt <- watson()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 21) {
     n <- 10
     tt <- ex_rosen()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 22) {
     n <- 20
     tt <- ex_powell()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 23) {
     n <- 10
     mm <- n + 1
     tt <- penalty_1()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 24) {
     n <- 10
     mm <- n + 1
     tt <- penalty_2()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 25) {
     n <- 6
     mm <- n + 2
     tt <- var_dim()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 26) {
     n <- 8
     tt <- trigon()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 27) {
     n <- 8
     mm <- n
     tt <- brown_al()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 28) {
     n <- 6
     mm <- n
     tt <- disc_bv()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 29) {
     n <- 8
     mm <- n
     tt <- disc_ie()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 30) {
     n <- 8
     mm <- n
     tt <- broyden_tri()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 31) {
     n <- 8
     mm <- n
     tt <- broyden_band()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 32) {
     mm <- 10
     n <- 8
     tt <- linfun_fr()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }

  if (fnum == 33) {
    mm <- 10
    n <- 8
    tt <- linfun_r1()
    if (is.function(tt$x0)) {
      xx0<-tt$x0(n)
    }
    else xx0 <- tt$x0
    lo <- rep((min(xx0)-0.1), n)
    up <- rep((max(xx0)+0.1), n)
    break }

  if (fnum == 34) {
    mm <- 10
    n <- 8
    tt <- linfun_r1z()
    if (is.function(tt$x0)) {
      xx0<-tt$x0(n)
    }
    else xx0 <- tt$x0
    lo <- rep((min(xx0)-0.1), n)
    up <- rep((max(xx0)+0.1), n)
    break }

  if (fnum == 35) {
     n <- 8
     m <- n
     tt <- chebyquad()
     if (is.function(tt$x0)) {
       xx0<-tt$x0(n)
     }
     else xx0 <- tt$x0
     lo <- rep((min(xx0)-0.1), n)
     up <- rep((max(xx0)+0.1), n)
     break }
  }    
# NOTE: bounds are experimental only
  mask <- rep(1L, n) # masks set to "free" (not masked)
  val <- list(npar = n, fffn=tt$fn, ffgr=tt$gr, ffhe=tt$he, x0=xx0, lo=lo, up=up, 
              mask=mask, fname=fname, ameth=ameth)
#  cat("val:"); print(val); tmp<-readline('exit ffn')
  val
} # end fufn


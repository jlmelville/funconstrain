# fufnrun.R -- J C Nash 2024-4-8
## ?? fixing kkt
# RFO.txt is input file
source("./fufn.R") # ensure fufn() loaded
library(funconstrain)  # get the functions
library(optimx)
mycon<-file("RFO.txt", open="r", blocking = TRUE)
sfname<-readLines(mycon, n=1)
if (length(sfname) == 0) {
   cat("no sink file\n")
} else { 
  cat("opening sink file ",sfname,"\n")
  sink(sfname, split=TRUE) 
} # open sink file
cat("sink file name=",sfname,"\n")

lin2 <- readLines(mycon, n=1)
cat("probs =",lin2,"\n")
if (length(lin2) == 0) stop("Unexpected null probs")
txt<-paste("probc<-c(",lin2,")","")
tryparse<-eval(parse(text=txt))
# ?? should we check it worked?
cat("Problem numbers:\n"); print(probc)
# check loop
for (iprob in probc){ # loop over problems
  if ( (iprob < 1) || (iprob > 35) ) {
    stop('Problem number out of range. Stopping.')
  }
} # end check loop
meths <- readLines(mycon, n=1)
if (length(meths) == 0) stop("Unexpected null meths")
cat("Methods:\n")
cat(meths,"\n")
methvec<-paste("methc<-c(",meths,")","")
tryparse<-eval(parse(text=methvec))
# cat("methods in list form:"); print(methc)
tbounds<-readLines(mycon, n=1)
have.bounds<-FALSE
if (tbounds == "TRUE") have.bounds<-TRUE
cat("have.bounds:",have.bounds,"\n")
close(mycon)
for (iprob in probc){ # loop over problems
   tfun <- fufn(fnum=iprob)
   # print(tfun)
   cat("Problem:", tfun$fname,"\n")
   x0 <- tfun$x0
   if (have.bounds){
     lo <- tfun$lo
     up <- tfun$up
   } 
   else {
     lo <- -Inf
     up <- Inf
   }
   tfn <- tfun$fffn
   attr(tfn, "fname") <- tfun$fname
   tgr <- tfun$ffgr
   the <- tfun$ffhe
   nx0<-length(x0)
#   cat("about to call opm\n")
   if (have.bounds) {
     t21 <-opm(x0, tfn, tgr, hess=the, lower=lo, upper=up, method=methc, 
               contro=list(trace=0))
   } else {
     t21 <-opm(x0, tfn, tgr, hess=the, method=methc, contro=list(trace=0))
   }
   print(summary(t21, order=value, par.select=1:min(nx0,5)))
   cat("END :", tfun$fname,"\n\n")
}
sink()


# fufnrun.R -- J C Nash 2024-4-8
## ?? fixing kkt
# Assume fufn.R has been loaded
source("./fufn.R")
sfname <- readline("Sink name=")
# source("~/optimr.R")
sink(sfname, split=TRUE)
library(funconstrain)  # get the functions
library(optimx)
tmp <- readline("begin")
# iprob <- as.numeric(readline("Prob #"))
iprob <- 1
while (iprob %in% 1:35){
tfun <- ffn(fnum=iprob)
# print(tfun)
cat("Problem:", tfun$fname,"\n")
x0 <- tfun$x0
lo <- tfun$lo
up <- tfun$up
tfn <- tfun$fffn
attr(tfn, "fname") <- tfun$fname
tgr <- tfun$ffgr
the <- tfun$ffhe
ameth<-unlist(tfun$ameth)
# ameth<-"ALL"
# ?? masking?
cat("about to call opm\n")
nx0<-length(x0)
t21 <-opm(x0, tfn, tgr, hess=the, lower=lo, upper=up, method=ameth, contro=list(trace=0))
print(summary(t21, order=value, par.select=1:min(nx0,5)))
cat("END :", tfun$fname,"\n\n")

# iprob <- as.numeric(readline("Prob #")) # Must have as.numeric
 tmp <- readline("next")
iprob <- iprob + 1
}
sink()

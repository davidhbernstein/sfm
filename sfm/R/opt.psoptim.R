opt.psoptim   <- function(fn, start_v, lower.psoptim, upper.psoptim=NA, maxit.psoptim, psopt.TF, rand.order = TRUE){
start_feval   <-  fn(start_v)
if(isTRUE(psopt.TF ==TRUE)){  
set.seed(1234)
  
opt00 <- psoptim(par     = start_v, 
                 fn      = fn,
                 lower   = lower.psoptim , 
                 upper   = upper.psoptim,
                 control = list(trace          = 1,
                                REPORT         = maxit.psoptim/10,
                                trace.stats    = TRUE,
                                maxit          = maxit.psoptim,
                                rand.order     = rand.order))
  
if(isTRUE(start_feval > opt00$value )) {start_v <- opt00$par
start_feval   <-  fn(start_v)} else{print("no improvement from psoptim")}
}
assign("start_v",     start_v,     envir=parent.frame())
assign("start_feval", start_feval, envir=parent.frame())  
}
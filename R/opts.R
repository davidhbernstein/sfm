opt.bobyqa    <- function(fn, start_v, lower.bobyqa, maxit.bobyqa, bob.TF, rhobeg = NA, rhoend  =NA){
start_feval   <-  fn(start_v)
if(isTRUE(bob.TF==TRUE)){  

bob1   <- bobyqa(par = start_v, 
                 fn = fn,
                 lower   = lower.bobyqa, 
                 control = list(iprint = 2, 
                                maxfun = maxit.bobyqa,
                                rhobeg = rhobeg,
                                rhoend = rhoend))
  
if(isTRUE(start_feval > bob1$fval )) {start_v <- bob1$par
start_feval   <-  fn(start_v)}else{print("no improvement from bobyqa")}
}
assign("start_v",     start_v,     envir=parent.frame())
assign("start_feval", start_feval, envir=parent.frame())  
}

opt.optim     <- function(fn, start_v, lower.optim, upper.optim, maxit.optim, opt.TF, method, optHessian, trace=1){
  start_feval   <-  fn(start_v)
  if(isTRUE(opt.TF ==TRUE)){
    
    opt <- optim(par     = start_v, 
                 fn      = fn,
                 lower   = lower.optim, 
                 upper   = upper.optim,
                 hessian = optHessian, 
                 method  = method,
                 control = list(maxit   = maxit.optim, 
                                REPORT  = maxit.optim/10, 
                                trace   = trace))
    
    if(isTRUE(start_feval > opt$value )) {start_v <- opt$par
    start_feval   <-  fn(start_v)} else{print("no improvement from optim")}
  }
  assign("start_v",     start_v, envir=parent.frame())
  assign("opt",         opt,     envir=parent.frame())
  assign("start_feval", start_feval, envir=parent.frame())
}

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
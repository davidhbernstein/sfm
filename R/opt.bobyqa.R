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
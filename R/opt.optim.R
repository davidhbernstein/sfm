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
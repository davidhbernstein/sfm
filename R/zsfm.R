zsfm <- function(formula, 
                model_name    = c("ZISF","ZISF_Z"),
                data, 
                maxit.bobyqa  = 10000,
                maxit.psoptim = 1000,
                maxit.optim   = 1000, 
                REPORT        = 1, 
                trace         = 2,
                pgtol         = 0, 
                start_val     = FALSE,
                PSopt         = FALSE,
                bob           = TRUE,
                optHessian    = TRUE,
                inefdec       = TRUE,
                upper         = NA,
                Method        = "L-BFGS-B",
                logit         = TRUE){
  
data_proc(formula,   data, model_name, individual = NULL, inefdec)
start_cs( formula_x ,data_orig, x_vars_vec, intercept, model_name, n_x_vars, start_val) 
data_proc2(data, data_x, fancy_vars, fancy_vars_z, data_z, y_var, x_vars_vec, halton_num=NA, individual=NA, N, model_name)

if(model_name %in% c("ZISF","ZISF_Z") ){
like.fn = function(x){
      
      if(model_name %in% c("ZISF","ZISF_Z")){x_x_vec <- x[4:as.numeric(n_x_vars+3)]}
        
      eps     <- (inefdec_n*(Y  - as.matrix(data_i_vars)%*%x_x_vec))
      
      if(model_name == "ZISF"){
  
  	## Calculate logit/probit probabilities for p parameter 
  	## in ZISF
  	gamma <- x[1]
  	prob <- exp(-abs(gamma)) 
  	
  	## Now the variance parameters
  	sigvsq <- x[2]^2
  	sigusq <- x[3]^2

  	sigv <- sqrt(sigvsq)
  	sigu <- sqrt(sigusq)

  	## Reparametrize the log-likelihood function
  	lambda     <- sigu/sigv
  	sigsq      <- sigvsq+sigusq
  	sig        <- sqrt(sigsq)           

  	## Now the likelihood function
  	e <- eps

  	f1 <- -0.5*log(2*pi*sigvsq)-(0.5/sigvsq)*e^2
  	f2 <- log(2/sig)+log(dnorm(e/sig))+log(pnorm(e*lambda/sig))
  	f  <- prob*exp(f1)+(1-prob)*exp(f2)

  	like <- log(f+1e-10) }

      like[like==-Inf]         <-  -sqrt(.Machine$double.xmax/length(like))
      like[like== Inf]         <-  -sqrt(.Machine$double.xmax/length(like))
      like[is.nan(like)]       <-  -sqrt(.Machine$double.xmax/length(like))
      
      return(-sum(like[is.finite(like)]))}  
  
start.time()
   
opt.bobyqa(fn=like.fn, start_v=start_v, lower.bobyqa=lower_bob, maxit.bobyqa=maxit.bobyqa, bob.TF=bob) 

lower.start(start_v, model_name, differ=1)

opt.psoptim(fn=like.fn, start_v, lower.psoptim=lower1,
            upper.psoptim=lower1, maxit.psoptim, psopt.TF=PSopt, rand.order = FALSE)  

lower.start(start_v, model_name, differ=0.5)
opt.optim(fn = like.fn, start_v = start_v, lower.optim =lower1,
          upper.optim=upper1, maxit.optim=maxit.optim, opt.TF=TRUE, method=Method, optHessian= TRUE)

end.time(start_time)    
if(optHessian==FALSE){st_err  <- rep(NA,length(opt$par))}
if(optHessian==TRUE){ st_err  <- if (isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0 ) ){ rep(NA,length(opt$par)) }   else{sqrt(diag(solve(opt$hessian)))}}
t_val      <- opt$par/st_err
out[1,]    <- opt$par
out[2,]    <- st_err
out[3,]    <- t_val
print(t(out))
    
if(model_name %in% c("ZISF") ){

beta <- opt$par[-c(1:3)]
eps  <- (inefdec_n*(Y  - as.matrix(data_i_vars)%*%beta))  
  
z <- 1
  
  ## Calculate logit/probit probabilities for p parameter 
  ## in ZISF
  gamma <- opt$par[1]
  prob  <- exp(-gamma)
  
  ## Now the variance parameters
  sigvsq <- opt$par[2]^2
  sigusq <- opt$par[3]^2
  
  sigv <- sqrt(sigvsq)
  sigu <- sqrt(sigusq)
  
  ## Reparametrize the log-likelihood function
  lambda     <- sigu/sigv
  sigsq      <- sigvsq+sigusq
  sig        <- sqrt(sigsq)           
  
  ## Now the likelihood function
  e <- eps 
  
  f1 <- -0.5*log(2*pi*sigvsq)-(0.5/sigvsq)*e^2
  f2 <- log(2/sig)+log(dnorm(e/sig))+log(pnorm(e*lambda/sig))
  f   <- prob*exp(f1)+(1-prob)*exp(f2)
  
  # print.default(f)
  
  post.prob <- prob*exp(f1)/f
  
  mustar      <- e*sigusq/sigsq
  sigstarsq   <- sigusq*sigvsq/(sigusq+sigvsq)
  sigstar     <- sqrt(sigstarsq)
  
  zz     <- mustar/sigstar
  jlms   <- mustar+sigstar*dnorm(zz)/pnorm(zz) }

  if(model_name %in%  c("ZISF","ZISF_Z") ){
    ret_stuff <- list(t(out),c(opt),total_time,start_v,model_name,formula, jlms,post.prob)
    names(ret_stuff)  <- c("out","opt","total_time","start_v","model_name","formula","jlms","post.prob")}else{print("No model name")}
    
    return(ret_stuff)}

else {return(c("This is not a valid command"))}}




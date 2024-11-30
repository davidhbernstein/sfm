sfm <- function(formula, 
                model_name    = c("NHN","NHN-MDPD","NHN-PSI","NHN-MLQE","NHN_Z","NE","NE_Z","NR","THT","NTN"), 
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
                eta           = 0.01,   
                alpha         = 0.2){
  
data_proc(formula,   data, model_name, individual = NULL, inefdec)
start_cs( formula_x ,data_orig, x_vars_vec, intercept, model_name, n_x_vars, start_val,n_z_vars,z_vars) 
data_proc2(data, data_x, fancy_vars, fancy_vars_z, data_z, y_var, x_vars_vec, halton_num=NA, individual=NA, N, model_name)

if(model_name %in% c("NHN","NE","NR","NHN-MDPD","NHN-PSI","NHN-MLQE","THT","NTN","NHN_Z","NE_Z") ){
like.fn = function(x){
      
      if(model_name %in% c("NHN","NE","NR","NHN-MDPD","NHN-PSI","NHN-MLQE")){x_x_vec <- x[3:as.numeric(n_x_vars + 2)]}
      if(model_name ==     "THT"){                                      x_x_vec <- x[4:as.numeric(n_x_vars + 3)]}
      if(model_name ==     "NTN"){                                      x_x_vec <- x[4:as.numeric(n_x_vars + 3)]}

      if(model_name %in% c("NE_Z","NHN_Z")){data_z_vars <- as.matrix(data.frame(subset(data,select = z_vars)))
                                            x_x_vec     <- x[2:as.numeric(n_x_vars+1)]
                                            z_z_vec     <- x[as.numeric(n_x_vars+2):as.numeric(length(start_v))]}

      eps     <- (inefdec_n*(Y  - as.matrix(data_i_vars)%*%x_x_vec))
      
      if(model_name=="NHN_Z"){
        sigma_u_fun    <- exp(as.matrix(data_z_vars)%*%z_z_vec)
        sigma_v_fun    <- x[1]
        sigma_fun      <- sqrt(sigma_v_fun^2 + sigma_u_fun^2)
        lamb_fun       <- sigma_u_fun/sigma_v_fun
        like           <-       log(pmax(   (2/sigma_fun)  * 
                                              dnorm(eps/sigma_fun)*  
                                              pnorm(-eps*lamb_fun/sigma_fun)  , eps*0+.Machine$double.eps) )}
      
      if(model_name=="NE_Z"){
        sigma_u_fun    <- exp(data_z_vars%*%z_z_vec)
        sigv           <- x[1]
        l1             <- log(1/sigma_u_fun)
        l2             <- pnorm( -(eps/sigv) - (sigv    /sigma_u_fun), log.p = TRUE)
        l3             <- (eps/sigma_u_fun) + (sigv^2 /  (2*sigma_u_fun^2) )
        like           <-  l1+l2+l3}
      
      if(model_name == "NHN"){
      like  <-      as.numeric(log(     pmax(   (2/x[2])    * 
                                              dnorm(eps/x[2]) *
                                              pnorm(-eps*x[1]/x[2])  ,  eps*0+.Machine$double.eps )    ))}
      
      if(model_name == "NE"){
      l1    <- log(1/x[2])
      l2    <- pnorm( -(eps/x[1]) - (x[1]    /x[2]), log.p = TRUE)
      l3    <- (eps/x[2]) + (x[1]^2 /  (2*x[2]^2)  )
      
      like  <-  l1+l2+l3}

      if(model_name=="NR"){
        sigv           <- x[1]
        sigu          <- x[2]
        sigma          <- sqrt(2*sigv**2+sigu^2)
        z              <- (eps*sigu/sigv)/sigma
        like           <- (log(sigv)- 2*log(sigma) - 1/2*(eps/sigv)^2 + 1/2*z^2 + 
                          log(sqrt(2/pi)*exp(-1/2*z**2) - z*(1-erf(z/sqrt(2)))))}
      
      if(model_name == "NHN-MLQE"){
        NNN    <- length(data)
        QQ     <- 1 - (1/ (10*log(NNN+10)) )
        like0  <- (2/x[2])  *    dnorm(eps/x[2])*  pnorm(-eps*x[1]/x[2])    
        like   <- sum(  (like0^(1-QQ) - 1) /  (1-QQ)  )}
      
      if(model_name == "NHN-MDPD"){
        fn_inner  <- function(z){exp(-(1+alpha)*(z^2/2))*pnorm(-z*x[1])^(1+alpha)}
        like      <-  -as.numeric((1/ (x[2]^alpha))* 
                               (( (sqrt(2) / sqrt(pi))*integrate(f=fn_inner,lower = -Inf,upper = Inf)[[1]])-    
                                  ((1+(1/alpha))*(1/N)*sum(exp(-alpha*(eps^2/ (2*x[2]^2)))*pnorm(-(eps*x[1])/x[2])^alpha))))}
      
      if(model_name == "NHN-PSI"){
        fn_int  <- function(z){ 
          like0 <- (2 / x[2]) * dnorm(z/x[2]) * pnorm(-z*x[1]/x[2])
          like  <- like0^(eta+1) / (eta+1)
          return(like)}
        
        fn_sum  <- function(z){ 
          like0 <- (2 / x[2]) * dnorm(z/x[2]) * pnorm(-z*x[1]/x[2])
          like  <- like0^eta / eta
          return(like)}
        
        like  <-  mean(fn_sum(eps)) - integrate(f=fn_int, lower = -Inf, upper = Inf)[[1]] }
      
      if(model_name=="THT"){
        sig_u   <- x[1]
        sig_v   <- x[2]
        a       <- x[3]
        lamb    <- sig_u/sig_v
        sig     <- sqrt(sig_v^2 + sig_u^2)
        like    <- as.numeric(log(pmax(2*dt(eps, df=a)*
                   pt((-eps*lamb/sig)*sqrt((a+1)/(sig^{-2}*eps^2 + a))  ,df=a+1)  ,  eps*0+.Machine$double.eps))) }
      
      if(model_name=="NTN"){
        lam  <- x[1]
        sig  <- x[2]
        mu   <- x[3]
        
        l1   <- -log(sig^2)/2
        l2   <- -log(2*pi)/2
        l3   <- -(1/(2*sig^2))*(-eps-mu)^2  
        l4   <-  pnorm(((mu/lam)-eps*lam)/sig,   log.p=TRUE)  
        l5   <- -pnorm((mu/sig)*sqrt(1+lam^(-2)),log.p=TRUE)  
        like <- l1 + l2 + l3 + l4 + l5}
      
      like[like==-Inf]         <-  -sqrt(.Machine$double.xmax/length(like))
      like[like== Inf]         <-  -sqrt(.Machine$double.xmax/length(like))
      like[is.nan(like)]       <-  -sqrt(.Machine$double.xmax/length(like))
      
      return(-sum(like[is.finite(like)]) ) }  
  
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

## TE Measurements     
if(model_name %in% c("NHN","NHN-MDPD","NHN-PSI","NHN-MLQE") ){
beta  <- opt$par[-c(1:2)]
lamb  <- opt$par[1]
sig   <- opt$par[2]
sig_u <- (lamb*sig) / sqrt(1+lamb^2)
sig_v <- sig_u/lamb
eps_hat    <- inefdec_n*(Y - rowSums(t(t(data_i_vars)*beta))) 
sig_star   <- sig_u*sig_v/sig
inner      <- (lamb*eps_hat)/sig
exp_u_hat  <- ((1-pnorm((sig_u*sig_v/sig) + inner ))/  pmax( (1-pnorm(inner)), .Machine$double.eps)    )*exp( (sig_u^2/sig^2)*(eps_hat + 0.5*sig_v^2) )
exp_u_hat  <- pmax(exp_u_hat, 0)
exp_u_hat  <- pmin(exp_u_hat, 1)}

if(model_name=="NHN_Z"){
NX         <- n_x_vars + 1
NZ1        <- n_x_vars + 2
NZ2        <- n_x_vars + n_z_vars + 1 
beta       <- opt$par[c(2:NX)]
delta      <- opt$par[c(NZ1:NZ2)]
sig_v      <- opt$par[1]
sig_u      <- exp((as.matrix(as.matrix(data.frame(subset(data,select = z_vars)))))%*%delta)
lamb       <- sig_u/sig_v 
sig        <- sqrt(sig_u^2 + sig_v^2)
eps_hat    <- inefdec_n*(Y - rowSums(t(t(data_i_vars)*beta)))
sig_star   <- (sig_u*sig_v)  /sig
inner      <- (lamb*eps_hat) /sig
exp_u_hat  <- ((1-pnorm((sig_u*sig_v/sig) + inner ))/  pmax( (1-pnorm(inner)), .Machine$double.eps)  )*exp( (sig_u^2/sig^2)*(eps_hat + 0.5*sig_v^2) )
exp_u_hat  <- pmax(exp_u_hat, 0)
exp_u_hat  <- pmin(exp_u_hat, 1)}

    if(model_name %in%  c("NHN","NHN-MDPD","NHN-PSI","NHN-MLQE","NHN_Z") ){
    ret_stuff <- list(t(out),c(opt),total_time,start_v,model_name,formula, exp_u_hat)
    names(ret_stuff)  <- c("out","opt","total_time","start_v","model_name","formula","exp_u_hat")}else{
    ret_stuff <- list(t(out),c(opt),total_time,start_v,model_name,formula)
    names(ret_stuff)  <- c("out","opt","total_time","start_v","model_name","formula")}

    return(ret_stuff)}

else {return(c("This is not a valid command"))}}




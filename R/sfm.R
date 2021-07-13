sfm <- function(formula, 
                model_name = c("NHN","NHN_Z","NE","NE_Z","THT","NTN"), 
                data, 
                maxit=100000,
                maxit_2 =1000, 
                REPORT=1, 
                trace=2,
                pgtol=0, 
                start_val=FALSE,
                inefdec=TRUE, 
                bob=TRUE,
                PSopt=FALSE){
  
  data_orig  <- data  
  
  ## lm code for data   
  data_conform <- function (formula, data,na.action, method = "qr", 
                            model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
                            contrasts = NULL, offset, ...) 
  {
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
                 "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    if (method == "model.frame")return(mf) else if (method != "qr") 
      warning(gettextf("method = '%s' is not supported. Using 'qr'", method), domain = NA)
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf, contrasts)
    return(x)}
  
  form_parts <- base::strsplit(as.character(formula), "|", fixed = TRUE)
  if(isFALSE(unique(grepl( "|", deparse(formula), fixed = TRUE))) ){ formula_x  <- formula
  y_var     <- gsub(" ", "", noquote(as.character( unlist( strsplit( deparse(formula_x), "~", fixed = TRUE)[[1]][1]))))}else{
  formula_x <- paste(form_parts[[2]], "~",form_parts[[3]][1], sep = "")
  y_var     <- gsub(" ", "", noquote(as.character( unlist(strsplit( formula_x, "~", fixed = TRUE))[[1]])))}
  
  if(length(unlist(form_parts))>3){
    formula_z <- paste(form_parts[[2]], "~",form_parts[[3]][2], sep = "")
    z_vars    <- as.character(gsub(" ", "", noquote(as.character(unlist(form_parts)[[4]])))) 
    z_vars    <- noquote(gsub("+", " ", z_vars, fixed=TRUE))
    z_vars    <- unlist(strsplit(z_vars, " "))
    data_z    <- data_conform(formula = formula_z, data = data)
    if(model_name=="NHN"){model_name <- "NHN_Z"}
    if(model_name=="NE"){model_name <- "NE_Z"}
    if(model_name %in% c("THT","NTN","CHC","NU")){print("Currently building this functionality") 
    return(c("Currently building this functionality"))}}
  
  data_x    <- data_conform(formula = formula_x, data = data)
  
  method         <- "L-BFGS-B"
  
  intercept      <- if(isTRUE(grepl(-1, gsub("[[:space:]]", "",as.character(formula_x))))) {0} else{1}
  inefdec_n      <- if(isTRUE(inefdec) ) {1} else{-1}
  inefdec_TF     <- if(isTRUE(inefdec) ) {TRUE} else{FALSE}
  
  x_vars_vec     <- if(model_name %in% c("WMLE","FD") & intercept==1){colnames(data_x)[-c(1)]}else {colnames(data_x)}
  n_x_vars       <- length(x_vars_vec)
  x_vars         <- x_vars_vec
  x_x_vec        <- rep(0,length= n_x_vars)
  fancy_vars   <- setdiff(colnames(data_x),colnames(data))
  fancy_vars_z <- NULL 
  
  if(length(unlist(form_parts))>3){    
    intercept_z    <- if(isTRUE(grepl(-1, gsub("[[:space:]]", "",as.character(formula_z))))) {0} else{1}
    z_vars_vec     <- if(model_name %in% c("WMLE","FD") & intercept_z==1){colnames(data_z)[-c(1)]}else {colnames(data_z)}
    n_z_vars       <- length(z_vars_vec)
    z_vars         <- z_vars_vec
    z_z_vec        <- rep(0,length= n_z_vars)
    fancy_vars_z   <- setdiff(colnames(data_z),colnames(data))}
  
  ## Starting values  
  plm_lm         <- lm(formula_x ,data_orig)
  
  ## let lm code handle missing observations
  data           <- data[rownames(data_x),]
  N              <- nrow(unique(data))
  
  if (isTRUE(length(fancy_vars)>0)) { data_inter   <- cbind(  data, data_x[,fancy_vars]  )
  colnames(data_inter) <- c(colnames(data),fancy_vars)
  data  <- data_inter} else  {}
  
  if (isTRUE(length(fancy_vars_z)>0)) { data_inter   <- cbind(  data, data_z[,fancy_vars_z]  )
  colnames(data_inter) <- c(colnames(data),fancy_vars_z)
  data  <- data_inter} else  {}
  
  beta_hat       <- if(isTRUE(intercept==0)) {plm_lm$coefficients[x_vars_vec]} else{plm_lm$coefficients[x_vars_vec][-1]}
  epsilon_hat    <- plm_lm$residuals
  beta_0_st      <- if(isTRUE(intercept==0)) {NA} else{plm_lm$coefficients[c(1)]}
  # sfa_eps        <- sfa(epsilon_hat   ~1, ineffDecrease = inefdec_TF)
  # exp_u          <- sfa_eps$mleParam[c(1)]
  sigma_u        <- 0.1  #coef(sfa_eps,extraPar=TRUE)[c("sigmaU")]
  sigma_v        <- 0.1  #coef(sfa_eps,extraPar=TRUE)[c("sigmaV")]
  mu             <- 0.1 
  beta_0         <- beta_0_st  
  lambda         <- sigma_u/sigma_v
  sigma          <- sqrt(sigma_u^2 + sigma_v^2)
  
  start_v_ntn    <- if(is.na(beta_0_st)) {unname(c(lambda,sigma,mu,beta_hat))} else{unname(c(lambda,sigma,mu,beta_0,beta_hat)) }
  start_v_t      <- if(is.na(beta_0_st)) {unname(c(sigma_u,sigma_v,1,beta_hat))} else{unname(c(sigma_u,sigma_v,1,beta_0,beta_hat)) }
  start_v_ne     <- if(is.na(beta_0_st)) {unname(c(sigma_v,sigma_u,beta_hat))} else{unname(c(sigma_v,sigma_u,beta_0,beta_hat)) }
  start_v_nhn    <- if(is.na(beta_0_st)) {unname(c(lambda,sigma,beta_hat))} else{unname(c(lambda,sigma,beta_0,beta_hat)) }
  
  if(model_name=="NHN"){
    start_v <-  start_v_nhn
    out            <- matrix(0,nrow = 3,ncol = length(start_v))
    colnames(out)  <- c("lambda","sigma",c(names(plm_lm$coefficients)))
    lower_bob <- c(rep(.Machine$double.eps,2),rep(-Inf,n_x_vars))}
  if(model_name=="NE"){
    start_v <- start_v_ne
    out            <- matrix(0,nrow = 3,ncol = length(start_v))
    colnames(out)  <- c("sigv","sigu",c(names(plm_lm$coefficients)))
    lower_bob <- c(rep(.Machine$double.eps,2),rep(-Inf,n_x_vars))}
  if(model_name=="THT"){
    start_v <- start_v_t
    out            <- matrix(0,nrow = 3,ncol = length(start_v))
    colnames(out)  <- c("sigv","sigu","a",c(names(plm_lm$coefficients))) 
    lower_bob <- c(rep(.Machine$double.eps,3),rep(-Inf,n_x_vars))}
  if(model_name=="NTN"){
    start_v <-  start_v_ntn
    out            <- matrix(0,nrow = 3,ncol = length(start_v))
    colnames(out)  <- c("lambda","sigma","mu",c(names(plm_lm$coefficients))) 
    lower_bob <- c(rep(.Machine$double.eps,2),rep(-Inf,n_x_vars+1))}
  if(model_name %in% c("NHN_Z","NE_Z") ){
      start_v <- start_v_nhn
      out            <- matrix(0,nrow = 3,ncol = length(start_v))}
  
  if (isTRUE(is.numeric(start_val))) {start_v <- start_val} 
  
  rownames(out)  <- c("par","st_err","t-val") 
  
  Y             <-  as.matrix(subset(data,select = y_var))
  data_i_vars   <-  as.matrix(data.frame(subset(data,select = x_vars_vec)))
  
  upper <- NA
  
  ###################################################################

  if(model_name %in% c("NHN","NE","THT","NTN") ){

    print(start_time <- Sys.time())
    
    fn = function(x){
      
      if(model_name %in% c("NHN","NE")){x_x_vec <- x[3:as.numeric(n_x_vars + 2)]}
      if(model_name ==     "THT"){      x_x_vec <- x[4:as.numeric(n_x_vars + 3)]}
      if(model_name ==     "NTN"){      x_x_vec <- x[4:as.numeric(n_x_vars + 3)]}
      
      eps     <- (inefdec_n*(Y  - as.matrix(data_i_vars)%*%x_x_vec))
      
      if(model_name == "NHN"){
      like  <-      as.numeric(log(     pmax(   (2/x[2])    * 
                                              dnorm(eps/x[2]) *
                                              pnorm(-eps*x[1]/x[2])  ,  eps*0+.Machine$double.eps )    ))}
      if(model_name == "NE"){
      l1    <- log(1/x[2])
      l2    <- pnorm( -(eps/x[1]) - (x[1]    /x[2]), log.p = TRUE)
      l3    <- (eps/x[2]) + (x[1]^2 /  (2*x[2]^2)  )
      
      like  <-  l1+l2+l3}
      
      if(model_name=="THT"){
        sig_u   <- x[1]
        sig_v   <- x[2]
        a       <- x[3]
        lamb    <- sig_u/sig_v
        sig     <- sqrt(sig_v^2 + sig_u^2)
        
        like  <-  as.numeric(log(     pmax(2*dt(eps, df=a)*
                                             pt((-eps*lamb/sig)*sqrt((a+1)/(sig^{-2}*eps^2 + a))  ,df=a+1)  ,  eps*0+.Machine$double.eps)))  }
      
      if(model_name=="NTN"){
        lam   <- x[1]
        sig   <- x[2]
        mu    <- x[3]
        
        l1 <- -log(sig^2)/2
        l2 <- -log(2*pi)/2
        l3 <- -(1/(2*sig^2))*(-eps-mu)^2  
        l4 <-  pnorm(((mu/lam)-eps*lam)/sig,   log.p=TRUE)  
        l5 <- -pnorm((mu/sig)*sqrt(1+lam^(-2)),log.p=TRUE)  
        
        like <- l1 + l2 + l3 + l4 + l5}
      
      like[like==-Inf]         <-  -sqrt(.Machine$double.xmax/length(like))
      like[like== Inf]         <-  -sqrt(.Machine$double.xmax/length(like))
      like[is.nan(like)]       <-  -sqrt(.Machine$double.xmax/length(like))
      
      return(-sum(like[is.finite(like)]) ) }  
    
    ## Optimization time 
    print(start_fun   <- fn(start_v))
    
    if(isTRUE(bob == TRUE)){
      bob1   <- bobyqa(par = start_v, 
                       fn = fn,
                       lower = lower_bob ,
                       control = list(iprint = 2, maxfun= maxit)) }
    
    if (isTRUE(start_fun > bob1$fval) ) {start_v <- bob1$par}
    
    if (PSopt==TRUE){
    differ  <- 1
    
    if(model_name %in% c("NHN","NE","NTN") ){lower_pso <- c(rep(.00000000001,2), start_v[-c(1:2)] - differ )}
    if(model_name=="THT"){                   lower_pso <- c(rep(.00000000001,3), start_v[-c(1:3)] - differ )}

    opt00  <- psoptim(par = start_v,
                      fn=fn,
                      lower = lower_pso,
                      upper = c(start_v + differ ),
                      control = list(trace      = 1,
                                     REPORT     = maxit_2/10,
                                     maxit      = maxit_2,
                                     trace.stats= TRUE,
                                     rand.order = FALSE))
    
    if (bob1$fval > opt00$value) {start_v <-  opt00$par } }
    
    if(model_name %in% c("NHN","NE") ){lower_opt <- c(rep(.Machine$double.eps,2),rep(-Inf,n_x_vars))}
    if(model_name=="THT"){             lower_opt <- c(rep(.Machine$double.eps,3),rep(-Inf,n_x_vars))}
    if(model_name=="NTN"){             lower_opt <- c(rep(.Machine$double.eps,2),rep(-Inf,n_x_vars+1))}
    
    opt <- optim(par = start_v, 
                 fn = fn,
                 lower = lower_opt, upper=upper,
                 control = list(maxit = maxit, REPORT = REPORT, trace = trace, pgtol=pgtol),
                 hessian = TRUE, method = method)
    
    st_err     <- if (isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0 ) ){ rep(NA,length(opt$par)) }   else{sqrt(diag(solve(opt$hessian)))}
    t_val      <- opt$par/st_err
    out[1,]    <- opt$par
    out[2,]    <- st_err
    out[3,]    <- t_val
    print(t(out))
    
    if(model_name=="NHN"){
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
    
    print(total_time <- Sys.time() - start_time)
    
    if(model_name=="NHN"){
    ret_stuff <- list(t(out),c(opt),total_time,start_v,model_name,formula, exp_u_hat)
    names(ret_stuff)  <- c("out","opt","total_time","start_v","model_name","formula","exp_u_hat")}else{
    ret_stuff <- list(t(out),c(opt),total_time,start_v,model_name,formula)
    names(ret_stuff)  <- c("out","opt","total_time","start_v","model_name","formula")
    }
    
    return(ret_stuff)}
  
  ###################################################################
  
  if(model_name %in% c("NE_Z","NHN_Z")){
    print(start_time <- Sys.time())
    start_time <- Sys.time()
    
    n_z_vars       <- length(z_vars)
    z_z_vec        <- rep(0,length= n_z_vars)
    
    ## Starting values  
    plm_pcs        <- lm(formula_x ,data)
    beta_hat       <- if(isTRUE(intercept==0)){plm(formula_x ,data,effect = "individual")$coefficients[c(x_vars_vec)]} else{plm_pcs$coefficients[-c(1)]}
    epsilon_hat    <- plm_pcs$residuals
    beta_0_st      <- if(isTRUE(intercept==0)) {NA} else{plm_pcs$coefficients[c(1)]}
    beta_se        <- summary(plm_pcs)[[4]][c("(Intercept)"),2]
    sigma_v        <- 0.1 
    beta_0         <- beta_0_st  
    delta          <- rep(0.1,length(z_vars))
    
    start_v        <- if(is.na(beta_0_st)) {unname(c(sigma_v,beta_hat,delta))} else{unname(c(sigma_v,beta_0,beta_hat,delta))}
    
    if (isTRUE(is.numeric(start_val))) {start_v <- start_val} 
    
    out            <- matrix(0,nrow = 3,ncol = length(start_v))
    rownames(out)  <- c("par","st_err","t-val") 
    colnames(out)  <- c("sigma_v",c(names(plm_pcs$coefficients)),z_vars) 
    
    Y             <-  as.matrix(subset(data,select = y_var))
    data_i_vars   <-  as.matrix(data.frame(subset(data,select = x_vars_vec)))
    data_z_vars   <-  as.matrix(data.frame(subset(data,select = z_vars)))
    
    
    fn = function(x){
      for (q in 1:n_x_vars){
        v             <- q + 1
        x_x_vec[q]    <- x[v]}
      
      for (qq in 1:n_z_vars){
        v              <- qq  + 1 + n_x_vars
        z_z_vec[qq]    <- x[v]}
      
      eps     <- inefdec_n*(Y  - as.matrix(data_i_vars)%*%x_x_vec)
      
      if(model_name=="NHN_Z"){
      sigma_u_fun    <- exp(as.matrix(data_z_vars)%*%z_z_vec)
      sigma_v_fun    <- x[1]
      sigma_fun      <- sqrt(sigma_v_fun^2 + sigma_u_fun^2)
      lamb_fun       <- sigma_u_fun/sigma_v_fun
      
      like  <-       log(pmax(   (2/sigma_fun)  * 
                                   dnorm(eps/sigma_fun)*  
                                   pnorm(-eps*lamb_fun/sigma_fun)  , eps*0+.Machine$double.eps) )}
      
      if(model_name=="NE_Z"){
        sigma_u_fun    <- exp(data_z_vars%*%z_z_vec)
        sigv           <- x[1]
        
        l1    <- log(1/sigma_u_fun)
        l2    <- pnorm( -(eps/sigv) - (sigv    /sigma_u_fun), log.p = TRUE)
        l3    <- (eps/sigma_u_fun) + (sigv^2 /  (2*sigma_u_fun^2) )
        
        like  <-  l1+l2+l3
      }
      
      
      like[like==-Inf]         <-  -sqrt(.Machine$double.xmax/length(like))
      like[like== Inf]         <-  -sqrt(.Machine$double.xmax/length(like))
      like[is.nan(like)]       <-  -sqrt(.Machine$double.xmax/length(like))
      
      return(-sum(like[is.finite(like)]) ) }  
    
    print("function:  ")    
    start_feval   <-  fn(start_v)
    print(start_feval   <-  fn(start_v))
    
    if(isTRUE(bob==TRUE)){
      bob1   <- bobyqa(par = start_v, fn = fn,
                       lower = c(.Machine$double.eps , rep(-.Machine$double.xmax^.1,length(start_v[-c(1)]))   ),
                       control = list(iprint = 2, 
                                      maxfun = maxit,
                                      rhobeg = 0.01,
                                      rhoend  =1e-12))
      
      if(isTRUE(start_feval > bob1$fval )) {start_v <- bob1$par} else{print("no improvement from bobyqa")}
    }
    
    
    start_feval   <-  fn(start_v)
    differ  <- 1
    
    if(isTRUE(PSopt==TRUE)){
      
      differ  <- 1
      
      opt00 <- psoptim(par     = start_v, 
                       fn      = fn,
                       lower = c(0.00001 , rep(start_v[-c(1)]) -differ  ),
                       upper   = c(start_v+differ),
                       control = list(trace          = 1,
                                      REPORT         = maxit_2/10,
                                      trace.stats    = TRUE,
                                      maxit          = maxit_2))
      
      if(isTRUE(start_feval > opt00$value )) {start_v <- opt00$par} else{print("no improvement from psoptim")}
      
    }
    
    differ  <- 1
    
    opt <- optim(par = start_v, fn = fn,
                 lower = c(0.00001 , start_v[-c(1)] -differ  ),
                 control = list(maxit = maxit_2, 
                                REPORT = REPORT, 
                                trace = trace, 
                                pgtol=pgtol,
                                factr  = 1e9),   
                 hessian = TRUE,
                 method = method)
    
    start_v    <- opt$par
    
    st_err     <- if (isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0 ) ){ rep(NA,length(opt$par)) }   else{sqrt(diag(solve(opt$hessian)))}
    t_val      <- opt$par/st_err
    out[1,]    <- opt$par
    out[2,]    <- st_err
    out[3,]    <- t_val
    print(t(out))
    
    if(model_name=="NHN_Z"){
    ## TE Measurements 
    NX    <- n_x_vars + 1
    NZ1   <- n_x_vars + 2
    NZ2   <- n_x_vars + n_z_vars + 1 
    
    beta   <- opt$par[c(2:NX)]
    delta  <- opt$par[c(NZ1:NZ2)]
    
    sig_v  <- opt$par[1]
    sig_u  <- exp((as.matrix(data_z_vars))%*%delta)
    lamb   <- sig_u/sig_v 
    sig    <- sqrt(sig_u^2 + sig_v^2)
    
    eps_hat    <- inefdec_n*(Y - rowSums(t(t(data_i_vars)*beta)))
    
    sig_star   <- (sig_u*sig_v)  /sig
    inner      <- (lamb*eps_hat) /sig
    exp_u_hat  <- ((1-pnorm((sig_u*sig_v/sig) + inner ))/  pmax( (1-pnorm(inner)), .Machine$double.eps)  )*exp( (sig_u^2/sig^2)*(eps_hat + 0.5*sig_v^2) )
    
    exp_u_hat  <- pmax(exp_u_hat, 0)
    exp_u_hat  <- pmin(exp_u_hat, 1)}
    
    print(total_time <- Sys.time() - start_time)

    if(model_name=="NHN_Z"){
      ret_stuff <- list(t(out),c(opt),total_time,start_v,model_name,formula, exp_u_hat)
      names(ret_stuff)  <- c("out","opt","total_time","start_v","model_name","formula","exp_u_hat")}else{
        ret_stuff <- list(t(out),c(opt),total_time,start_v,model_name,formula)
        names(ret_stuff)  <- c("out","opt","total_time","start_v","model_name","formula")
      }
    
    return(ret_stuff)} 
  
  ###################################################################

  else {return(c("This is not a valid command"))}}




psfm <- function(formula, 
                 model_name   = c("TRE_Z","GTRE_Z","TRE","GTRE","WMLE","FD","GTRE_SEQ1","GTRE_SEQ2"), 
                 data, 
                 maxit.bobyqa = 100,
                 maxit.psoptim= 10,
                 maxit.optim  = 10,
                 REPORT       = 1,
                 trace        = 3, 
                 pgtol        = 0, 
                 individual,
                 halton_num   = NULL,
                 start_val    = FALSE,
                 gamma        = FALSE,
                 PSopt        = FALSE,
                 bob          = TRUE,
                 optHessian   = TRUE,
                 inefdec      = TRUE, 
                 Method       = "L-BFGS-B"){
    
data_proc(  formula,   data, model_name, individual, inefdec)
start_panel(formula_x, data, model_name, start_val, intercept, x_vars_vec)
data_proc2(data, data_x, fancy_vars, fancy_vars_z, data_z, y_var, x_vars_vec, halton_num, individual, N, model_name)

if(model_name %in% c("GTRE","TRE") ){
  
fn_1 = function(x){
          
   if(model_name == "GTRE"){x_x_vec       <- x[5:as.numeric(n_x_vars + 4)]}
   if(model_name == "TRE"){ x_x_vec       <- x[4:as.numeric(n_x_vars + 3)]}
          
 fn1 = function(ii){
            
   if(model_name == "GTRE"){eps[[ii]]     <- Y[[ii]]  - x[3]*R_h1[[ii]] + x[4]*R_h2[[ii]] * inefdec_n}
   if(model_name == "TRE"){ eps[[ii]]     <- Y[[ii]]  - x[3]*R_h1[[ii]]}
    
            for (qq in 1:n_x_vars) {
              eps[[ii]] <- eps[[ii]] - x_x_vec[qq]*matrix(rep(data_i_vars[[ii]][,qq],R),t[[ii]],R)  
            }
   
            eps[[ii]] <- inefdec_n*eps[[ii]]
            z0          <-  2/x[2]
            z1          <-  eps[[ii]]/x[2]
            z2          <- -eps[[ii]]*x[1]/x[2]
            z1[z1>  37] <-  37 
            z1[z1< -37] <- -37 
            z2[z2>   8] <-   8 
            z2[z2< -37] <- -37 
            
            prod_vec_n  <- log(max(   mean(   colProds(    
                         z0* dnorm(z1)* pmax(pnorm(z2), eps[[ii]]*0+.Machine$double.xmin)    
                       )), .Machine$double.xmin))
            return(-prod_vec_n)}
          
fn1_apply                         <- unlist(lapply(1:N, fn1))
fn1_apply[is.nan(fn1_apply)]      <- sqrt(.Machine$double.xmax/length(x))
fn1_apply[is.infinite(fn1_apply)] <- sqrt(.Machine$double.xmax/length(x))
          
return( sum( fn1_apply[is.finite(fn1_apply)] ) ) }

start.time()      
lower.start(start_v, model_name, differ=3)
opt.bobyqa(fn=fn_1, start_v=start_v, lower.bobyqa=lower1, maxit.bobyqa=maxit.bobyqa, bob.TF=bob)
      
lower.start(start_v, model_name, differ=2)
opt.psoptim(fn=fn_1, start_v=start_v, lower.psoptim=lower1,upper.psoptim=upper1, maxit.psoptim, psopt.TF=PSopt)

lower.start(start_v, model_name, differ=0.5)
opt.optim(fn = fn_1, start_v = start_v, lower.optim =lower1 ,upper.optim=upper1, 
          maxit.optim=maxit.optim, opt.TF=TRUE, method=Method, optHessian= optHessian)

end.time(start_time)      
st_err     <- if (isTRUE(any(opt$hessian==0))  | optHessian ==FALSE ){ rep(NA,length(opt$par)) }   else{sqrt(pmax(diag(solve(opt$hessian)),0))}
t_val      <- opt$par/st_err
out[1,]    <- opt$par
out[2,]    <- st_err
out[3,]    <- t_val
print(t(out))
      
      ## TE Measurements : GTRE       
      if(model_name == "GTRE"){ 
        
        beta  <- opt$par[-c(1:4)]
        lamb  <- opt$par[1]
        sig   <- opt$par[2]
        sig_u <- (lamb*sig) / sqrt(1+lamb^2)
        sig_v <- sig_u/lamb
        sig_r <- opt$par[3]
        sig_h <- opt$par[4]
        
        e_i  <-   as.list(rep(0,N))
        A    <-   as.list(rep(0,N))
        SIG  <-   as.list(rep(0,N))
        VEE  <-   as.list(rep(0,N))
        LAM  <-   as.list(rep(0,N))
        ARR  <-   as.list(rep(0,N))
        n    <-   sum(t)
        U    <-   rep(0,n)
        
        for(i in 1:N){
          e_i[[i]]     <- pmin(inefdec_n*(Y[[i]][,1] - rowSums(t(t(data_i_vars[[i]]) * beta))),Y[[i]][,1]*0)
          A[[i]]       <- -cbind(rep(1,t[i]),diag(t[i]))
          SIG[[i]]     <- sig_v^2*diag(t[i]) + sig_r^2*rep(1,t[i])%*%t(rep(1,t[i]))
          VEE[[i]]     <- rbind( c(sig_h^2,rep(0,t[i]))   ,  cbind( rep(0,t[i]) , sig_u^2*diag(t[i])) )
          LAM[[i]]     <- solve( solve(VEE[[i]])  +  t(A[[i]]) %*% solve(SIG[[i]]) %*% A[[i]]    )
          ARR[[i]]     <- LAM[[i]] %*% t(A[[i]]) %*% solve(SIG[[i]])
        }
        
        res_d_fn <- function(i){ptmvnorm(lowerx = rep(0, t[i]+1), upperx = rep(Inf, t[i]+1), 
                                         mean = as.numeric(ARR[[i]] %*% e_i[[i]] ),
                                         sigma=LAM[[i]])[1]}
        
        res_n_fn <- function(i){ptmvnorm(lowerx = rep(0, t[i]+1), upperx=rep(Inf, t[i]+1), 
                                         mean =  as.numeric(ARR[[i]] %*% e_i[[i]] + LAM[[i]]%*% c(-1,rep(0,t[i])) ), 
                                         sigma=LAM[[i]])[1]}
        
        res_d <- lapply(seq(1,N,1), res_d_fn)
        res_n <- lapply(seq(1,N,1), res_n_fn)
        
        H_fn    <- function(i){(max(res_n[[i]], .Machine$double.eps)/max(res_d[[i]],.Machine$double.eps) )*
            exp(t(c(-1,rep(0,t[i])))%*%ARR[[i]]%*%e_i[[i]] + 0.5*t(c(-1,rep(0,t[i])))%*%LAM[[i]]%*%c(-1,rep(0,t[i])) )}
        
        H    <- unlist(lapply(seq(1,N,1), H_fn))
        H    <- pmin(H,rep(1,length(H)))
        print(paste("Mean TE of H: ",round(mean(H),3), sep = ""))
        
        
        new_t_exp <- as.list(rep(0,n))
        
        for(i in 1:N){
          for (j in 1:t[i]) {
            h <- cumsum(t)[i]-t[i] + j
            new_t      <- rep(0,t[i] +1)
            new_t[j+1] <- -1 
            new_t_exp[[h]] <- new_t
            
          }}
        
        
t_cum      <-   c(cumsum(t))
t_exp      <-   e_i_exp    <-   A_exp  <- SIG_exp <- VEE_exp <- LAM_exp <- ARR_exp <- res_d_exp  <-   as.list(rep(0,n))
        
for(m in 1:N){
B  <- t_cum[m]
A  <- B +1 - t[m] 
t_exp[A:B]     <- rep(t[m],  t[m])
e_i_exp[A:B]   <- rep(e_i[m],t[m])
A_exp[A:B]     <- rep(A[m],  t[m])
SIG_exp[A:B]   <- rep(SIG[m],t[m])
VEE_exp[A:B]   <- rep(VEE[m],t[m])
LAM_exp[A:B]   <- rep(LAM[m],t[m])
ARR_exp[A:B]   <- rep(ARR[m],t[m])
res_d_exp[A:B] <- rep(res_d[m],t[m])
}
        
        res_n_t_fn <- function(i){ptmvnorm(  lowerx=rep(0, t_exp[[i]]+1), upperx=rep(Inf, t_exp[[i]]+1),
                                             mean= as.numeric(ARR_exp[[i]] %*% e_i_exp[[i]] + LAM_exp[[i]]%*% new_t_exp[[i]] ), 
                                             sigma=LAM_exp[[i]])[1]}
        
        res_n_t <- lapply(seq(1,n,1), res_n_t_fn)
        
        
        U_fn    <- function(i){(max(res_n_t[[i]], .Machine$double.eps) /  max(res_d_exp[[i]], .Machine$double.eps))*
            exp(t(new_t_exp[[i]])%*%ARR_exp[[i]]%*%e_i_exp[[i]] + 0.5*t(new_t_exp[[i]])%*%LAM_exp[[i]]%*%new_t_exp[[i]] )}
        
        U     <- unlist(lapply(seq(1,n,1), U_fn))
        U     <- pmin(U,rep(1,length(U)))
        print(paste("Mean TE of U: ",round(mean(U),3), sep = ""))  }
      
      ## TE Measurements : TRE       
      if(model_name == "TRE"){
        beta  <- opt$par[-c(1:3)]
        lamb  <- opt$par[1]
        sig   <- opt$par[2]
        sig_u <- (lamb*sig) / sqrt(1+lamb^2)
        sig_v <- sig_u/lamb
        
        Y_mean           <- rep(0,N)
        X_mean           <- matrix(0,N,ncol = length(x_vars_vec))
        colnames(X_mean) <- x_vars_vec
        
        for (ii in 1:N) {
          data_i[[ii]]  <- data[which(data[,c(individual)]==indiv[ii]),]
          Y_mean[ii]    <- mean(as.numeric(data_i[[ii]][,y_var]))
          X_mean[ii,]   <- colMeans(data.frame(data_i[[ii]][,c(x_vars_vec)] ))
        }
        
        r_hat_m      <- Y_mean - rowSums(t(beta*t(X_mean))) + inefdec_n*sqrt(2/pi)*sig_u
        r_hat_m_exp  <- rep(0,sum(t))
        t_cum        <- c(cumsum(t))
        
        for(m in 1: length(t)){
          B  <- t_cum[m]
          A  <- B +1 - t[m] 
          r_hat_m_exp[A:B] <- rep(r_hat_m[m],t[m])}
        
eps_hat    <- pmin(inefdec_n*(data[,y_var] - rowSums(t(t(data[,c(x_vars_vec)])*beta)) -r_hat_m_exp ),data[,y_var]*0)
sig_star   <- sig_u*sig_v/sig
inner      <- (lamb*eps_hat)/sig
exp_u_hat  <- ((1-pnorm((sig_u*sig_v/sig) + inner ))/(1-pnorm(inner)))*exp( (sig_u^2/sig^2)*(eps_hat + 0.5*sig_v^2) )
U          <- exp_u_hat}
      
#cor(U,exp(-p_data_trial$u))
#cor(H,exp(-unique(p_data_trial$h) ))

if(model_name == "GTRE"){      
ret_stuff <- list(t(out),c(opt),total_time,start_v,model_name,formula, U, H)
names(ret_stuff)  <- c("out","opt","total_time","start_v","model_name","formula","U","H")}
if(model_name == "TRE"){      
ret_stuff <- list(t(out),c(opt),total_time,start_v,model_name,formula, U)
names(ret_stuff)  <- c("out","opt","total_time","start_v","model_name","formula","U")}
      
return(ret_stuff)}   
  
if(model_name %in% c("GTRE_Z","TRE_Z") ){

delta          <- rep(0.1,length(z_vars))
      
if(isTRUE(is.numeric(start_val))){start_v <- start_val}
if(isTRUE(model_name == "GTRE_Z") & isFALSE(is.numeric(start_val))){start_v  <- if(is.na(beta_0_st)) {unname(c(sigma_v,sigma_r,sigma_h,beta_hat,delta))} else{unname(c(sigma_v,sigma_r,sigma_h,beta_0,beta_hat,delta))} }
if(isTRUE(model_name == "TRE_Z") & isFALSE(is.numeric(start_val))){start_v   <- if(is.na(beta_0_st)) {unname(c(sigma_v,sigma_r,beta_hat,delta))} else{unname(c(sigma_v,sigma_r,beta_0,beta_hat,delta))}}

out            <- matrix(0,nrow = 3,ncol = length(start_v))
rownames(out)  <- c("par","st_err","t-val") 
colnames(out)  <- if(model_name=="GTRE_Z"){c("sigv","sigr","sigh",colnames(data_x),z_vars)} else{c("sigv","sigr",colnames(data_x),z_vars)} 
      
      if (isTRUE(is.numeric(start_val))) {start_v <- start_val}
      
      if (isTRUE(is.numeric(halton_num))) {R <- halton_num}else{
      R              <- ceiling(sqrt(nrow(data)))+100 }          ## Integral reps   
      R_H            <- randtoolbox::halton(R,8,normal = TRUE)   ## halton seqs
      
      indiv           <- noquote(as.vector(unique(data[,c(individual)])))
      t               <- rep(0, N)
      data_i <- Y <- eps <- data_i_vars <- data_z_vars <- R_h1 <- R_h2 <- as.list(rep(0,N))
      
      for (ii in 1:N) {
        data_i[[ii]]        <-  data[which(data[,c(individual)]==indiv[ii]),]
        t[ii]               <-  nrow(data_i[[ii]])
        R_h1[[ii]]          <-  t(matrix(rep(R_H[,7],t[[ii]]),R,t[[ii]]))
        R_h2[[ii]]          <-  abs(t(matrix(rep(R_H[,8],t[[ii]]),R,t[[ii]])))
        Y[[ii]]             <-  matrix(rep(data_i[[ii]][,y_var],R),t[[ii]],R)
        data_i_vars[[ii]]   <- data.frame(data_i[[ii]][,c(x_vars_vec)] )
        data_z_vars[[ii]]   <- data.frame(data_i[[ii]][,c(z_vars)]     )} 
      
fn <-  function(x){
if(model_name == "GTRE_Z"){      
x_x_vec    <- x[4:as.numeric(n_x_vars + 3)]

for (qq in 1:n_z_vars){
v              <- qq  + 3 + n_x_vars
z_z_vec[qq]    <- x[v]}
}
      
if(model_name == "TRE_Z"){  
x_x_vec    <- x[3:as.numeric(n_x_vars + 2)]

for (qq in 1:n_z_vars){
v              <- qq  + 2 + n_x_vars
z_z_vec[qq]    <- x[v]} }
          
fn1 = function(ii){ 
  if(model_name == "GTRE_Z"){
   eps[[ii]]     <- Y[[ii]]  - x[2]*R_h1[[ii]]  + x[3]*R_h2[[ii]] * inefdec_n}
  
  if(model_name == "TRE_Z"){   
   eps[[ii]]     <- Y[[ii]]  - x[2]*R_h1[[ii]]   }
  
  for (qq in 1:n_x_vars) {
    eps[[ii]] <- eps[[ii]] - x_x_vec[qq]*matrix(rep(data_i_vars[[ii]][,qq],R),t[[ii]],R)  
  }
  eps[[ii]] <- inefdec_n*eps[[ii]]
  
  sigma_u_fun    <- sqrt(exp(as.matrix(data_z_vars[[ii]])%*%z_z_vec))
  sigma_v_fun    <- x[1]
  sigma_fun      <- sqrt(sigma_v_fun^2 + sigma_u_fun^2)
  lamb_fun       <- sigma_u_fun/sigma_v_fun     
  
  prod_vec_n  <- log(mean(colProds((2/matrix(rep(sigma_fun,R),t[[ii]],R))* 
                 dnorm(eps[[ii]]/matrix(rep(sigma_fun,R),t[[ii]],R))*   
       pmax(pnorm(-eps[[ii]]*matrix(rep(lamb_fun,R),t[[ii]],R)/matrix(rep(sigma_fun,R),t[[ii]],R)),eps[[ii]]*0+.Machine$double.eps))))
  
  return(-prod_vec_n)
}

fn1_apply  <- unlist(lapply(1:N, fn1))

fn1_apply[which(fn1_apply==Inf)]   <-  (.Machine$double.xmax)^.1
fn1_apply[which(fn1_apply==-Inf)]  <- -(.Machine$double.xmax)^.1

return( sum( fn1_apply[is.finite(fn1_apply)] ) )} 

start.time()
differ<- 1.5
if(model_name == "TRE_Z"){lower1   <- c(rep(.Machine$double.eps,2) , start_v[-c(1:2)] -differ)}
if(model_name == "GTRE_Z"){lower1   <- c(rep(.Machine$double.eps,3) , start_v[-c(1:3)] -differ)}    
        
opt.bobyqa(fn=fn, start_v=start_v, lower.bobyqa=lower1, maxit.bobyqa=maxit.bobyqa, bob.TF=bob)      

differ<- 10
if(model_name == "TRE_Z"){lower1   <- c(rep(.Machine$double.eps,2) , start_v[-c(1:2)] -differ)}
if(model_name == "GTRE_Z"){lower1   <- c(rep(.Machine$double.eps,3) , start_v[-c(1:3)] -differ)}  

opt.psoptim(fn=fn, start_v, lower.psoptim=lower1, upper.psoptim=c(start_v+differ), 
            maxit.psoptim=maxit.psoptim, psopt.TF=PSopt)
      
differ  <- 1
if(model_name == "TRE_Z"){lower1   <- c(rep(.Machine$double.eps,2) , start_v[-c(1:2)] -differ )}
if(model_name == "GTRE_Z"){lower1  <- c(rep(.Machine$double.eps,3) , start_v[-c(1:3)] -differ)}

opt.optim(fn = fn, start_v = start_v, lower.optim =lower1 ,upper.optim=c(start_v+differ), 
          maxit.optim=maxit.optim, opt.TF=TRUE, method=Method, optHessian= optHessian)

end.time(start_time)      
st_err     <- if (isTRUE(any(opt$hessian==0) | optHessian ==FALSE ) ){ rep(NA,length(opt$par)) }   else{sqrt(pmax(diag(solve(opt$hessian)),0))}
t_val      <- opt$par/st_err
out[1,]    <- opt$par
out[2,]    <- st_err
out[3,]    <- t_val
print(t(out))
      
      ## TE Measurements       
      if(model_name == "GTRE_Z"){         
        
        NX    <- n_x_vars + 3
        NZ1   <- n_x_vars + 4
        NZ2   <- n_x_vars + n_z_vars + 3 
        
        beta   <- opt$par[c(4:NX)]
        delta  <- opt$par[c(NZ1:NZ2)]
        
        sig_v  <- max(opt$par[1],0.00001)
        sig_r  <- max(opt$par[2],0.00001)
        sig_h  <- max(opt$par[3],0.00001)
        
        sig_u  <- sqrt(exp(as.matrix(data.frame(subset(data,select = z_vars)))%*%delta))
        # lamb   <- sig_u/sig_v 
        # sig    <- sqrt(sig_u^2 + sig_v^2)
        
        e_i  <-   as.list(rep(0,N))
        A    <-   as.list(rep(0,N))
        SIG  <-   as.list(rep(0,N))
        VEE  <-   as.list(rep(0,N))
        LAM  <-   as.list(rep(0,N))
        ARR  <-   as.list(rep(0,N))
        n    <-   sum(t)
        U    <-   rep(0,n)
        
        for(i in 1:N){
          e_i[[i]]     <- pmin(inefdec_n*(Y[[i]][,1] - rowSums(t(t(data_i_vars[[i]]) * beta))),Y[[i]][,1]*0)
          A[[i]]       <- -cbind(rep(1,t[i]),diag(t[i]))
          SIG[[i]]     <- sig_v^2*diag(t[i]) + sig_r^2*rep(1,t[i])%*%t(rep(1,t[i]))
          VEE[[i]]     <- rbind( c(sig_h^2,rep(0,t[i]))   ,  cbind( rep(0,t[i]) , sig_u[i]^2*diag(t[i])) )
          LAM[[i]]     <- solve( solve(VEE[[i]])  +  t(A[[i]]) %*% solve(SIG[[i]]) %*% A[[i]]    )
          ARR[[i]]     <- LAM[[i]] %*% t(A[[i]]) %*% solve(SIG[[i]])
        }
        
        res_d_fn <- function(i){ptmvnorm(lowerx = rep(0, t[i]+1), upperx = rep(Inf, t[i]+1), 
                                         mean = as.numeric(ARR[[i]] %*% e_i[[i]] ),
                                         sigma=LAM[[i]])[1]}
        
        res_n_fn <- function(i){ptmvnorm(lowerx = rep(0, t[i]+1), upperx=rep(Inf, t[i]+1), 
                                         mean =  as.numeric(ARR[[i]] %*% e_i[[i]] + LAM[[i]]%*% c(-1,rep(0,t[i])) ), 
                                         sigma=LAM[[i]])[1]}
        
        res_d <- lapply(seq(1,N,1), res_d_fn)
        res_n <- lapply(seq(1,N,1), res_n_fn)
        
        H_fn    <- function(i){(max(res_n[[i]], .Machine$double.eps)/max(res_d[[i]],.Machine$double.eps) )*
            exp(t(c(-1,rep(0,t[i])))%*%ARR[[i]]%*%e_i[[i]] + 0.5*t(c(-1,rep(0,t[i])))%*%LAM[[i]]%*%c(-1,rep(0,t[i])) )}
        
        H    <- unlist(lapply(seq(1,N,1), H_fn))
        H    <- pmin(H,rep(1,length(H)))
        # print(paste("Mean TE of H: ",round(mean(H),3), sep = ""))
        
        
        new_t_exp <- as.list(rep(0,n))
        
        for(i in 1:N){
          for (j in 1:t[i]) {
            h <- cumsum(t)[i]-t[i] + j
            new_t      <- rep(0,t[i] +1)
            new_t[j+1] <- -1 
            new_t_exp[[h]] <- new_t
            
          }}
        
        t_cum      <-   c(cumsum(t))
        t_exp      <-   as.list(rep(0,n))
        e_i_exp    <-   as.list(rep(0,n))
        A_exp      <-   as.list(rep(0,n))
        SIG_exp    <-   as.list(rep(0,n))
        VEE_exp    <-   as.list(rep(0,n))
        LAM_exp    <-   as.list(rep(0,n))
        ARR_exp    <-   as.list(rep(0,n))
        res_d_exp  <-   as.list(rep(0,n))
        
        for(m in 1:N){
          B  <- t_cum[m]
          A  <- B +1 - t[m] 
          t_exp[A:B]     <- rep(t[m],  t[m])
          e_i_exp[A:B]   <- rep(e_i[m],t[m])
          A_exp[A:B]     <- rep(A[m],  t[m])
          SIG_exp[A:B]   <- rep(SIG[m],t[m])
          VEE_exp[A:B]   <- rep(VEE[m],t[m])
          LAM_exp[A:B]   <- rep(LAM[m],t[m])
          ARR_exp[A:B]   <- rep(ARR[m],t[m])
          res_d_exp[A:B] <- rep(res_d[m],t[m])
        }
        
        res_n_t_fn <- function(i){ptmvnorm(  lowerx=rep(0, t_exp[[i]]+1), upperx=rep(Inf, t_exp[[i]]+1),
                                             mean= as.numeric(ARR_exp[[i]] %*% e_i_exp[[i]] + LAM_exp[[i]]%*% new_t_exp[[i]] ), 
                                             sigma=LAM_exp[[i]])[1]}
        
        res_n_t <- lapply(seq(1,n,1), res_n_t_fn)
        
        
        U_fn    <- function(i){(max(res_n_t[[i]], .Machine$double.eps) /  max(res_d_exp[[i]], .Machine$double.eps))*
            exp(t(new_t_exp[[i]])%*%ARR_exp[[i]]%*%e_i_exp[[i]] + 0.5*t(new_t_exp[[i]])%*%LAM_exp[[i]]%*%new_t_exp[[i]] )}
        
        U <- unlist(lapply(seq(1,n,1), U_fn))
        
        U <- pmin(U,rep(1,length(U)))   
      }
      
      
      ## TE Measurements            
      if(model_name == "TRE_Z"){
        NX    <- n_x_vars + 2
        NZ1   <- n_x_vars + 3
        NZ2   <- n_x_vars + n_z_vars + 2
        
        beta   <- opt$par[c(3:NX)]
        delta  <- opt$par[c(NZ1:NZ2)]
        
        sig_v  <- opt$par[1]
        sig_u  <- sqrt(exp((as.matrix(data.frame(data[,c(z_vars)] )))%*%delta))
        lamb   <- sig_u/sig_v
        sig    <- sqrt(sig_u^2 + sig_v^2)
        
        eps_hat    <- pmin( inefdec_n*(data[,y_var] - rowSums(t(t(data.frame(data[,c(x_vars_vec)] ))*beta))), 5)
        
        sig_star   <- (sig_u*sig_v)  /sig
        inner      <- (lamb*eps_hat) /sig
        U          <- ((1-pnorm((sig_u*sig_v/sig) + inner ))/ pmax(1-pnorm(inner),.Machine$double.eps) )*exp( (sig_u^2/sig^2)*(eps_hat + 0.5*sig_v^2) )
        U          <- pmax(U, 0)
        U          <- pmin(U, 1)
        
      }      

if(model_name == "GTRE_Z"){      
ret_stuff <- list(t(out),c(opt),total_time,start_v,model_name,formula, U, H)
names(ret_stuff)  <- c("out","opt","total_time","start_v","model_name","formula","U","H")}
if(model_name == "TRE_Z"){      
ret_stuff <- list(t(out),c(opt),total_time,start_v,model_name,formula, U)
names(ret_stuff)  <- c("out","opt","total_time","start_v","model_name","formula","U")}
return(ret_stuff)}    

if(model_name == "WMLE"){
start.wmle(formula_x, data, model_name, start_val, intercept, x_vars_vec, gamma, individual, N, y_var, plm_wmle)

like.wmle  <-  function(x){ 
x_x_vec    <- x[3:as.numeric(n_x_vars + 2)]  
    
fn1 = function(i){ 
      
eps_t    <- Y[[i]]  
      
for (qq in 1:n_x_vars) {
eps_t  <- eps_t - x_x_vec[qq]*demean(as.numeric(data_i_vars[[i]][,qq])) }
eps_t  <- eps_t*inefdec_n
eps_t1 <- eps_t[1:t[[i]]-1]
      
if(gamma==FALSE)   {E_t      <- ((x[1]*x[1])/t[[i]])*one_t[[i]]%*%t(one_t[[i]])}
else               {E_t      <- ( (x[1]/(1-x[1])) /t[[i]] )*one_t[[i]]%*%t(one_t[[i]])}
E_t1     <- (1/t[[i]])*one_t1[[i]]%*%t(one_t1[[i]])
l1       <- dmnorm(x= eps_t1,                  mean = rep(0,t[[i]]-1), varcov = x[2]*x[2]*(I_t1[[i]] - E_t1 ))
if(gamma==FALSE)   {l2       <- pmnorm(x= -(x[1]/x[2])*eps_t,  mean = rep(0,t[[i]]),       varcov =  I_t[[i]] + E_t )}
else               {l2       <- pmnorm(x= -(sqrt(x[1]/(1-x[1]))/x[2])*eps_t,  mean = rep(0,t[[i]]),  varcov =  I_t[[i]] + E_t )}
 
prod_vec_n  <- log(l1) + log(l2[1])  ## Log likelihood   
return(-prod_vec_n)
}
fn1_apply  <- unlist(lapply(1:N, fn1))
return(sum( fn1_apply[is.finite(fn1_apply)] ))}  

start.time()

differ <- 2
opt.bobyqa(fn=like.wmle, start_v=start_v, lower.bobyqa=c(rep(.Machine$double.eps,2),rep(-Inf,n_x_vars)), 
           maxit.bobyqa=maxit.bobyqa, bob.TF=bob)

opt.psoptim(fn=like.wmle, start_v, lower.psoptim=c(rep(.Machine$double.eps,2), start_v[-c(1:2)]-differ),
            upper.psoptim=c(start_v+differ), maxit.psoptim=maxit.psoptim, psopt.TF=PSopt)

opt.optim(fn = like.wmle, start_v = start_v, lower.optim =c(rep(.Machine$double.eps,2),rep(-Inf,n_x_vars)),
          upper.optim=c(start_v+differ), maxit.optim=maxit.optim, opt.TF=TRUE, method=Method, optHessian= optHessian)

end.time(start_time)      
beta             <- opt$par[-c(1:2)]  ## estimate the r_i's and then e(exp(u)|eps)'s 
lamb             <- opt$par[1]
sig              <- opt$par[2]
sig_u            <- (lamb*sig) / sqrt(1+lamb^2)
sig_v            <- sig_u/lamb
Y_mean           <- rep(0,N)
X_mean           <- matrix(0,N,ncol = length(x_vars_vec))
colnames(X_mean) <- x_vars_vec
      
for (ii in 1:N) {
        data_i[[ii]]  <- data[which(data[,c(individual)]==indiv[ii]),]
        Y_mean[ii]    <- mean(as.numeric(data_i[[ii]][,y_var]))
        X_mean[ii,]   <- colMeans(data.frame(data_i[[ii]][,c(x_vars_vec)] ))
      }
      
      r_hat_m      <- Y_mean - rowSums(t(beta*t(X_mean))) + inefdec_n*sqrt(2/pi)*sig_u
      r_hat_m_exp  <- rep(0,sum(t))
      t_cum        <- c(cumsum(t))
      
      for(m in 1: length(t)){
        B  <- t_cum[m]
        A  <- B +1 - t[m] 
        r_hat_m_exp[A:B] <- rep(r_hat_m[m],t[m])}
      
eps_hat    <- pmin(inefdec_n*(data[,y_var] - rowSums(t(t(data[,c(x_vars_vec)])*beta))- r_hat_m_exp),data[,y_var]*0)
sig_star   <- sig_u*sig_v/sig
inner      <- (lamb*eps_hat)/sig
exp_u_hat  <- ((1-pnorm((sig_u*sig_v/sig) + inner ))/(1-pnorm(inner)))*exp( (sig_u^2/sig^2)*(eps_hat + 0.5*sig_v^2) )
      
#cor(exp_u_hat,exp(-data$u))   
      
st_err     <- if (isTRUE(any(opt$hessian==0))  | optHessian ==FALSE ){ rep(NA,length(opt$par)) }   else{sqrt(pmax(diag(solve(opt$hessian)),0))}
t_val      <- opt$par/st_err
out[1,]    <- opt$par
out[2,]    <- st_err
out[3,]    <- t_val
print(t(out))
ret_stuff <- list(t(out),c(opt),total_time,start_v,r_hat_m,exp_u_hat ,model_name,formula, data)
names(ret_stuff) <- c("out","opt","total_time","start_v","r_hat_m","exp_u_hat","model_name","formula","data")
return(ret_stuff)
}
    
if(model_name == "FD"){
start.time()
      
if (isTRUE(is.numeric(start_val))) {start_v <- start_val} else{
        
        beta_hat         <- plm_fd$coefficients[c(x_vars_vec)]
        epsilon_hat      <- plm_fd$residuals
        beta_se          <- as.data.frame(summary(plm_fd)[1])$coefficients.Std..Error
        sfa_eps        <- pcs_c(Y  = as.numeric(epsilon_hat))[[1]]$par
        
        exp_u          <- sfa_eps[3]
        sigma_v        <- sqrt(sfa_eps[2]^2/ (1 + sfa_eps[1]^2))    #coef(sfa_eps,extraPar=TRUE)[c("sigmaV")]
        sigma_u        <- sigma_v*sfa_eps[1]                        #coef(sfa_eps,extraPar=TRUE)[c("sigmaU")]
        
        sigmaSq_u        <- sigma_u^2
        sigmaSq_v        <- sigma_v^2
        mu               <- 0.1
        delta            <- rep(0.1,length(z_vars))
        start_v          <- unname(c(sigmaSq_u,sigmaSq_v,mu,beta_hat,delta))}
      
      
      out              <- matrix(0,nrow = 3,ncol = length(start_v))
      rownames(out)    <- c("par","st_err","t-val") 
      colnames(out)    <- c("sig_u2","sig_v2","mu",x_vars_vec,z_vars) 
      
      
      indiv           <- noquote(as.vector(unique(data[,c(individual)])))
      t               <- rep(0, N)
      data_i          <- as.list(rep(0,N))
      Y               <- as.list(rep(0,N))
      eps             <- as.list(rep(0,N))
      data_i_vars     <- as.list(rep(0,N))
      SIGMA           <- as.list(rep(0,N))
      data_z_vars     <- as.list(rep(0,N))
      
      for (ii in 1:N) {
        data_i[[ii]]              <- data[which(data[,c(individual)]==indiv[ii]),]
        t[ii]                     <- nrow(data_i[[ii]])
        Y[[ii]]                   <- as.numeric(data_i[[ii]][,y_var])
        data_i_vars[[ii]]         <- data.frame(data_i[[ii]][,c(x_vars_vec)] )
        data_z_vars[[ii]]         <- data.frame(data_i[[ii]][,c(z_vars)] )
        
        if(t[ii] == 2) {SIGMA[[ii]] = 2}
        if(t[ii] == 3) {SIGMA[[ii]] = matrix(c(2,-1,-1,2),2,2)} 
        if(t[ii] > 3) {
          SIGMA[[ii]]               <- matrix(0, nrow = t[ii]-1, ncol = t[ii]-1)
          diag(SIGMA[[ii]])         <-  2
          diag(SIGMA[[ii]][-1,  ])  <- -1
          diag(SIGMA[[ii]][  ,-1])  <- -1}}
      
like.fd <- function(x){ 
        
        for (q in 1:n_x_vars){
          v             <- q + 3
          x_x_vec[q]    <- x[v]}
        
        for (qq in 1:n_z_vars){
          m <- v + qq
          z_z_vec[qq]    <- x[m] }
        
        fn1 = function(i){
          
          eps_t    <- diff(Y[[i]], lag = 1) 
          eps_h    <- rep(0, t[i])
          
          for (qq in 1:n_x_vars) {eps_t <- eps_t - x_x_vec[qq] * diff(as.numeric(data_i_vars[[i]][,qq]), lag = 1)}
          for (qq in 1:n_z_vars) {eps_h <- eps_h + z_z_vec[qq] * as.numeric(data_z_vars[[i]][,qq])  }
          
          eps_h     <- diff(exp(eps_h), lag = 1)    
          eps_t     <- eps_t*inefdec_n            ## not exactly sure if this is sufficient for cost
          SIG       <- x[2]*SIGMA[[i]]
          sig_star2 <- (t(eps_h) %*% qr.solve(SIG) %*% eps_h  + (1/x[1]))^-1
          mu_star   <- (  (x[3]/x[1]) - t(eps_t) %*% qr.solve(SIG) %*% eps_h  )*sig_star2
          
          l1     <-  - 0.5 * (t[i]-1) * log(2*pi) 
          l2     <-  - 0.5 * log(t[i])
          l3     <-  - 0.5 * (t[i]-1) * log(x[2])
          l4     <-  - 0.5 * t(eps_t) %*% qr.solve(SIG) %*% eps_t
          l5     <-    0.5 * ( (mu_star^2/sig_star2) -  ((x[3]^2) / x[1]) )
          l6     <-    log( sqrt(sig_star2)* max(pnorm(mu_star/ sqrt(sig_star2) ),.Machine$double.eps) )
          l7     <-  - log( sqrt(x[1]) *     max(pnorm( x[3] /  sqrt(x[1])      ),.Machine$double.eps) )
          
          
          prod_vec_n  <- sum(l1, l2, l3, l4, l5, l6, l7)      
          return(-prod_vec_n)}
        
        fn1_apply  <- unlist(lapply(1:N, fn1))
        
        return( sum( fn1_apply[is.finite(fn1_apply)] ) )}

differ  <- 2            

opt.bobyqa(fn=like.fd,  start_v=start_v, lower.bobyqa=c(rep(0.000001,2) , rep(-Inf,n_x_vars+1), rep(-Inf,length(z_vars))), 
           maxit.bobyqa=maxit.bobyqa, bob.TF=bob)      

opt.psoptim(fn=like.fd, start_v=start_v, lower.psoptim=c(rep(.Machine$double.eps,2), start_v[-c(1:2)]- differ),
            upper.psoptim=c(start_v+differ), maxit.psoptim, psopt.TF=PSopt)

opt.optim(fn = like.fd, start_v = start_v, lower.optim =c(rep(0.000001,2) , rep(-Inf,n_x_vars+1), rep(-Inf,length(z_vars))),
          upper.optim=c(start_v+differ),maxit.optim=maxit.optim, opt.TF=TRUE, method=Method, optHessian= optHessian)

end.time(start_time)
## TE Measurements: ui's
NX    <- n_x_vars + 3
NZ1   <- n_x_vars + 4
NZ2   <- n_x_vars + n_z_vars + 3
beta  <- opt$par[c(4:NX)]
delta <- opt$par[c(NZ1:NZ2)]
sigu2 <- opt$par[1]
sigv2 <- opt$par[2]
mu    <- opt$par[3]
u_hat <- rep(0, sum(t))
h_hat <- exp(as.matrix(data.frame(subset(data,select = z_vars)))%*%delta)

for (i in 1:N) {
for (tt in 1:t[i]) {
eps_t     <- diff(Y[[i]], lag = 1) 
eps_h     <- rep(0, t[i]-1)

for (qq in 1:n_x_vars) {
m         <- 3 + qq
eps_t     <- eps_t - opt$par[m] * diff(as.numeric(data_i_vars[[i]][,qq]), lag = 1)}

for (qq in 1:n_z_vars) {
m         <- 3 + length(n_x_vars) + qq 
eps_h     <- eps_h + opt$par[m] * diff(as.numeric(data_z_vars[[i]][,qq]), lag = 1)}

eps_t     <- eps_t*inefdec_n
SIG       <- sigv2*SIGMA[[i]]
sig_star2 <- (t(eps_h) %*% qr.solve(SIG) %*% eps_h  + (1/sigu2))^-1
mu_star   <- (  (mu/sigu2) - t(eps_t) %*% qr.solve(SIG) %*% eps_h  )*sig_star2 

num        <- if(i>1) { cumsum(t)[i-1]  + tt } else {tt}
u_hat[num]   <- h_hat[num] * (mu_star + (sqrt(sig_star2)*dnorm(mu_star/sqrt(sig_star2)) / max(pnorm(mu_star/sqrt(sig_star2)),.Machine$double.eps)   ) ) 
}}
      
exp_u_hat  <- exp(-u_hat)
exp_u_hat  <- pmax(exp_u_hat, 0)
exp_u_hat  <- pmin(exp_u_hat, 1)
      
st_err     <- if (isTRUE(any(opt$hessian==0))  | optHessian ==FALSE){ rep(NA,length(opt$par)) }   else{sqrt(diag(qr.solve(opt$hessian)))} 
t_val      <- opt$par/st_err 
out[1,]        <- opt$par
out[2,]        <- st_err
out[3,]        <- t_val
print(t(out))
ret_stuff <- list(t(out),c(opt),total_time,start_v,model_name,formula, u_hat, h_hat, exp_u_hat, data)
names(ret_stuff) <- c("out","opt","total_time","start_v","model_name","formula", "u_hat", "h_hat", "exp_u_hat", "data")
return(ret_stuff) 
}

if(model_name == "GTRE_SEQ1") {
start.time()
      
## Sequential Method  
      sfa_eps       <- sfa(epsilon_hat   ~1, ineffDecrease = inefdec_TF)
      sfa_alp       <- sfa(alpha_hat     ~1, ineffDecrease = inefdec_TF)
      exp_eta       <- sfa_alp$mleParam[c(1)]
      exp_u         <- sfa_eps$mleParam[c(1)]
      sigma_u       <- coef(sfa_eps,extraPar=TRUE)[c("sigmaU")]
      sigma_v       <- coef(sfa_eps,extraPar=TRUE)[c("sigmaV")]
      sigma_r       <- coef(sfa_alp,extraPar=TRUE)[c("sigmaU")]
      sigma_h       <- coef(sfa_alp,extraPar=TRUE)[c("sigmaV")]
      beta_0        <- beta_0_st + exp_u +exp_eta
      Lambda        <- sigma_u/sigma_v
      Sigma         <- sqrt(sigma_u^2 + sigma_v^2)
      sigma_se_r    <- NA
      sigma_r_se    <- NA
      sigma_h_se    <- NA
      lambda_se_r   <- NA
      beta_0_se     <- NA 
      gamma_uv      <- coef(sfa_eps,extraPar=TRUE)[c("gamma")]
      sigmaSq_uv    <- coef(sfa_eps,extraPar=TRUE)[c("sigmaSq")]
      gamma_hr      <- coef(sfa_alp,extraPar=TRUE)[c("gamma")]
      sigmaSq_hr    <- coef(sfa_alp,extraPar=TRUE)[c("sigmaSq")]
      gamma_uv_se   <- summary(sfa_eps)[[25]][c("gamma"),2]
      sigmaSq_uv_se <- summary(sfa_eps)[[25]][c("sigmaSq"),2]
      gamma_hr_se   <- summary(sfa_alp)[[25]][c("gamma"),2]
      sigmaSq_hr_se <- summary(sfa_alp)[[25]][c("sigmaSq"),2]
      
      other_parms            <- as.matrix(c(Lambda,Sigma,beta_0))
      rownames(other_parms)  <- c("lambda","sigma","beta_0")
      start_v                <- unname(c(gamma_uv,sigmaSq_uv,gamma_hr,sigmaSq_hr,beta_hat))
      out                    <- matrix(0,nrow = 3,ncol = length(start_v))
      rownames(out)          <- c("par","st_err","t-val") 
      colnames(out) <- if(isTRUE(intercept==0)) {c("gamma_uv","sigmaSq_uv","gamma_hr","sigmaSq_hr",colnames(data_x) )} else{c("gamma_uv","sigmaSq_uv","gamma_hr","sigmaSq_hr",colnames(data_x)[-c(1)] )}
      
      end.time(start_time)
      st_err     <- rep(NA,ncol(out)) 
      st_err     <- if(isTRUE(intercept==0)) {c(gamma_uv_se,sigmaSq_uv_se,gamma_hr_se,sigmaSq_hr_se,summary(plm_gtre)$coefficients[,2])} else{c(gamma_uv_se,sigmaSq_uv_se,gamma_hr_se,sigmaSq_hr_se,summary(plm_gtre)$coefficients[-c(1),2])}
      t_val      <- start_v/st_err
      out[1,]    <- start_v
      out[2,]    <- st_err
      out[3,]    <- t_val
print(t(out))
return(list(t(out),total_time,other_parms,model_name,formula, data))}

if(model_name == "GTRE_SEQ2") {
start.time()
## Sequential Method following 1995 paper  
## take second and third moments of alpha_hat and epsilon_hat
      alp_2m <- mean(alpha_hat^2)
      alp_3m <- min(0,mean(alpha_hat^3))
      eps_2m <- mean(epsilon_hat^2)
      eps_3m <- min(0,mean(epsilon_hat^3))
      
      gamma_uv   <- min(1, 1/(eps_2m*(sqrt(pi/2)*(pi/(pi-4))*eps_3m)^(-2/3) + (2/pi)  ))
      gamma_hr   <- min(1, 1/(alp_2m*(sqrt(pi/2)*(pi/(pi-4))*alp_3m)^(-2/3) + (2/pi)  ))
      sigmaSq_uv <- eps_2m + (2/pi)*(sqrt(pi/2)*(pi/(pi-4))*eps_3m)^(2/3)
      sigmaSq_hr <- alp_2m + (2/pi)*(sqrt(pi/2)*(pi/(pi-4))*alp_3m)^(2/3)
      beta_0     <- beta_0_st + sqrt(2/pi)*(sqrt(pi/2)*(pi/(pi-4))*eps_3m)^(1/3) + sqrt(2/pi)*(sqrt(pi/2)*(pi/(pi-4))*alp_3m)^(1/3)
      
      ## calculate the ten needed central moments 
      mu_2_eps   <- sigmaSq_uv*((1-gamma_uv) + gamma_uv*((pi -2)/pi) )
      mu_3_eps   <- sigmaSq_uv^(3/2)*(sqrt(2/pi)*(1-(4/pi))*gamma_uv^(3/2))
      mu_4_eps   <- sigmaSq_uv^2*(3*(1-gamma_uv)^2  + ((6*(pi-2)*gamma_uv*(1-gamma_uv))/pi) + gamma_uv^2*(3- (4/pi) - (12/pi^2)) )
      mu_5_eps   <- sigmaSq_uv^(5/2)*gamma_uv^(3/2)*sqrt(2/pi)*(10*(1-(4/pi))*(1-gamma_uv) + (7-(20/pi)-(16/pi^2))*gamma_uv )
      mu_6_eps   <- sigmaSq_uv^3*( 15*(1-gamma_uv)^3   +  
                                     (45*(pi-2)*(1-gamma_uv)^2*gamma_uv/pi)  +
                                     15*(3 - (4/pi) - (12/pi^2))*(1-gamma_uv)*gamma_uv^2  +
                                     (15  -  (6/pi) - (100/pi^2)  - (40/pi^3)  )*gamma_uv^3         )
      
      mu_2_alp   <- sigmaSq_hr*((1-gamma_hr) + gamma_hr*((pi -2)/pi) )
      mu_3_alp   <- sigmaSq_hr^(3/2)*(sqrt(2/pi)*(1-(4/pi))*gamma_hr^(3/2))
      mu_4_alp   <- sigmaSq_hr^2*(3*(1-gamma_hr)^2  + ((6*(pi-2)*gamma_hr*(1-gamma_hr))/pi) + gamma_hr^2*(3- (4/pi) - (12/pi^2)) )
      mu_5_alp   <- sigmaSq_hr^(5/2)*gamma_hr^(3/2)*sqrt(2/pi)*(10*(1-(4/pi))*(1-gamma_hr) + (7-(20/pi)-(16/pi^2))*gamma_hr )
      mu_6_alp   <- sigmaSq_hr^3*( 15*(1-gamma_hr)^3   +  
                                     (45*(pi-2)*(1-gamma_hr)^2*gamma_hr/pi)  +
                                     15*(3 - (4/pi) - (12/pi^2))*(1-gamma_hr)*gamma_hr^2  +
                                     (15  -  (6/pi) - (100/pi^2)  - (40/pi^3)  )*gamma_hr^3         )
      
      var_2m_eps  <-  (1/nrow(data))*(mu_4_eps - mu_2_eps^2)
      var_3m_eps  <-  (1/nrow(data))*(mu_6_eps - mu_3_eps^3 -6*mu_2_eps*mu_4_eps + 9*mu_2_eps) 
      cov_23m_eps <-  (1/nrow(data))*(mu_5_eps - 4*mu_2_eps*mu_3_eps)
      
      var_2m_alp  <-  (1/N)*(mu_4_alp - mu_2_alp^2)
      var_3m_alp  <-  (1/N)*(mu_6_alp - mu_3_alp^3 -6*mu_2_alp*mu_4_alp + 9*mu_2_alp) 
      cov_23m_alp <-  (1/N)*(mu_5_alp - 4*mu_2_alp*mu_3_alp)
      
      ## define 8 needed derivatives 
      d_beta_0_d_m3_eps  <- (pi/(pi-4))*(1/3)    *       (sqrt(pi/2)*(pi/(pi-4))*eps_3m)^(-2/3)
      d_sigma2_d_m3_eps  <- sqrt(2/pi)*(pi/(pi-4))*(2/3)*(sqrt(pi/2)*(pi/(pi-4))*eps_3m)^(-1/3)
      d_gamma_d_m2_eps   <- -(sqrt(pi/2)*(pi/(pi-4))*eps_3m)^(-2/3)*(eps_3m*(  sqrt(pi/2)*(pi/(pi-4))*eps_3m    )  + (2/pi))
      d_gamma_d_m3_eps   <- eps_2m*sqrt(pi/2)*(pi/(pi-4))*(2/3)*(sqrt(pi/2)*(pi/(pi-4))*eps_3m)^(-5/3)*(eps_3m*(sqrt(pi/2)*(pi/(pi-4))*eps_3m )^(-2/3) + (2/pi) )^(-2)
      
      d_beta_0_d_m3_alp  <- (pi/(pi-4))*(1/3)    *       (sqrt(pi/2)*(pi/(pi-4))*alp_3m)^(-2/3)
      d_sigma2_d_m3_alp  <- sqrt(2/pi)*(pi/(pi-4))*(2/3)*(sqrt(pi/2)*(pi/(pi-4))*alp_3m)^(-1/3)
      d_gamma_d_m2_alp   <- -(sqrt(pi/2)*(pi/(pi-4))*alp_3m)^(-2/3)*(alp_3m*(  sqrt(pi/2)*(pi/(pi-4))*alp_3m    )  + (2/pi))
      d_gamma_d_m3_alp   <- alp_2m*sqrt(pi/2)*(pi/(pi-4))*(2/3)*(sqrt(pi/2)*(pi/(pi-4))*alp_3m)^(-5/3)*(alp_3m*(sqrt(pi/2)*(pi/(pi-4))*alp_3m )^(-2/3) + (2/pi) )^(-2)
      
      beta_0_se     <- sqrt(beta_se[1]  + d_beta_0_d_m3_eps^2*var_3m_eps + d_beta_0_d_m3_alp^2*var_3m_alp )
      
      sigmaSq_uv_se <- sqrt(var_2m_eps + d_sigma2_d_m3_eps^2*var_3m_eps  + d_sigma2_d_m3_eps * cov_23m_eps)
      gamma_uv_se   <- sqrt(d_gamma_d_m2_eps^2*var_2m_eps  + d_gamma_d_m3_eps^2* var_3m_eps + d_gamma_d_m2_eps*d_gamma_d_m3_eps*cov_23m_eps     )
      
      sigmaSq_hr_se <- sqrt(var_2m_alp + d_sigma2_d_m3_alp^2*var_3m_alp  + d_sigma2_d_m3_alp * cov_23m_alp)
      gamma_hr_se   <- sqrt(d_gamma_d_m2_alp^2*var_2m_alp  + d_gamma_d_m3_alp^2* var_3m_alp + d_gamma_d_m2_alp*d_gamma_d_m3_alp*cov_23m_alp     )
      
      start_v        <- if(is.na(beta_0_st)) {unname(c(gamma_uv,sigmaSq_uv,gamma_hr,sigmaSq_hr,beta_hat))} else {unname(c(gamma_uv,sigmaSq_uv,gamma_hr,sigmaSq_hr,beta_0,beta_hat))}
      end.time(start_time)
      out            <- matrix(0,nrow = 3,ncol = length(start_v))
      rownames(out)  <- c("par","st_err","t-val") 
      colnames(out)  <- if(isTRUE(intercept==0)) {c("gamma_uv","sigmaSq_uv","gamma_hr","sigmaSq_hr",colnames(data_x) )} else{c("gamma_uv","sigmaSq_uv","gamma_hr","sigmaSq_hr",colnames(data_x) )}
      
      st_err     <- rep(NA,ncol(out)) 
      st_err     <- if(is.na(beta_0_st)){c(gamma_uv_se,sigmaSq_uv_se,gamma_hr_se,sigmaSq_hr_se,summary(plm_gtre)$coefficients[c(x_vars_vec),2])} else {c(gamma_uv_se,sigmaSq_uv_se,gamma_hr_se,sigmaSq_hr_se,beta_0_se,summary(plm_gtre)$coefficients[-c(1),2])}
      t_val      <- start_v/st_err
      out[1,]    <- start_v
      out[2,]    <- st_err
      out[3,]    <- t_val

      ## look at individual sigmas 
      Sigma_u <- sqrt(gamma_uv*sigmaSq_uv)    
      Sigma_v <- sqrt((1-gamma_uv)*sigmaSq_uv) 
      Sigma_h <- sqrt(gamma_hr*sigmaSq_hr)   
      Sigma_r <- sqrt((1-gamma_hr)*sigmaSq_hr) 
      Lambda  <- sigma_u/sigma_v
      Sigma   <- sqrt(sigmaSq_uv)
      
      other_parms    <- as.matrix(c(Sigma_u,Sigma_v,Sigma_h, Sigma_r, Lambda, Sigma))
      rownames(other_parms) <- c("sigma_u","sigma_v","sigma_h", "sigma_r","lambda","sigma")
      print(t(out))
      return(list(t(out),total_time,other_parms,model_name,formula, data))
    }
    
else {return(c("This is not a valid command"))}}


library(Rcpp)

# general sample estimates
require("compiler")

#sourceCpp( paste( getwd(),"/srcresidualcokurtosisSF.cpp",sep=""))
#sourceCpp( paste( getwd(),"/srcresidualcokurtosisMF.cpp",sep=""))

`%k%` <- function (X, Y) {
  # kronecker product when X has 2 rows and Y is a vector
  # out= [ Xy[1] | Xy[2] | ... | Xy[n] ]
  t(as.vector(X) %*% Y)
}

coskewnes = function (U){
  cAssets = ncol(U)
  T = nrow(U)
  M3 = matrix(rep(0, cAssets^3), nrow = cAssets, ncol = cAssets^2)
  for (t in 1:T) {
    centret = as.numeric(U[t, ])
    M3 = M3 + tcrossprod(centret) %k% t(centret)
  }
  return(1/T * M3)
}

coskewness <- cmpfun(coskewnes)

cokurtosi = function (U)
{
  cAssets = ncol(U)
  T = nrow(U)
  M4 = matrix(rep(0, cAssets^4), nrow = cAssets, ncol = cAssets^3)
  for (t in 1:T) {
    centret = as.numeric(U[t, ] )
    #M4 = M4 + (centret %*% t(centret)) %x% t(centret) %x% t(centret)
    M4 = M4 + tcrossprod(centret) %k% t(centret) %k% t(centret)
  }
  return(1/T * M4)
}

cokurtosis <- cmpfun(cokurtosi)


# single factor model


covarianceS <- function( beta , stockM2 , factorM2 ){
  beta = as.numeric(beta); # beta of the stock with the factor index
  stockM2 = as.numeric(stockM2); # idiosyncratic variance of the stock
  factorM2 = as.numeric(factorM2); #variance of the factor
  beta = matrix( beta , ncol = 1)
  S = (beta%*%t(beta))*factorM2
  D = diag( stockM2 )
  return( S + D )
}

covarianceSF = cmpfun(covarianceS)

coskewnessS <- function( beta , stockM3 , factorM3 ){
  N = length(beta);
  beta = as.numeric(beta); # beta of the stock with the factor index
  stockM3 = as.numeric(stockM3); # idiosyncratic third moment of the stock
  factorM3 = as.numeric(factorM3); #third moment of the factor
  beta = matrix( beta , ncol = 1)
  S = ((beta%*%t(beta))%x%t(beta))*factorM3
  D = matrix( rep(0,N^3) , ncol=N^2)
  for( i in 1:N ){
    col = (i-1)*N+i
    D[i,col] = stockM3[i]
  }
  return( S + D )
}

coskewnessSF = cmpfun(coskewnessS)

cokurtosisS <- function( beta , stockM2 , stockM4 , factorM2 , factorM4 ){
  N = length(beta);
  beta = as.numeric(beta); # beta of the stock with the factor index
  stockM2 = as.numeric(stockM2); # idiosyncratic second moment of the stock
  stockM4 = as.numeric(stockM4); # idiosyncratic fourth moment of the stock  
  factorM2 = as.numeric(factorM2); #second moment of the factor
  factorM4 = as.numeric(factorM4); #fourth moment of the factor  
  beta = matrix( beta , ncol = 1)
  S = ((beta%*%t(beta))%x%t(beta)%x%t(beta))*factorM4
  D = residualcokurtosisSF(sstockM2=stockM2, sstockM4=stockM4, mfactorM2 = factorM2, bbeta = beta )
  return( S + D )
}

cokurtosisSF = cmpfun(cokurtosisS)


# multi factor model

covarianceM <- function( beta , stockM2 , factorM2 ){
  beta = as.matrix(beta); # beta of the stock with the factor index
  stockM2 = as.numeric(stockM2); # idiosyncratic variance of the stock
  factorM2 = as.matrix(factorM2); #variance of the factor
  S = beta%*%factorM2%*%t(beta)
  D = diag( stockM2 )
  return( S + D )
}

covarianceMF = cmpfun(covarianceM)

coskewnessM <- function( beta , stockM3 , factorM3 ){
  beta = as.matrix(beta); # beta of the stock with the factor index
  N = nrow(beta) 
  stockM3 = as.numeric(stockM3); # idiosyncratic third moment of the stock
  factorM3 = as.matrix(factorM3); #third moment of the factor
  S = (beta%*%factorM3)%*%(t(beta)%x%t(beta))
  D = matrix( rep(0,N^3) , ncol=N^2)
  for( i in 1:N ){
    col = (i-1)*N+i
    D[i,col] = stockM3[i]
  }
  return( S + D )
}

coskewnessMF = cmpfun(coskewnessM)

cokurtosisM <- function( beta , stockM2 , stockM4 , factorM2 , factorM4 ){
  beta = as.matrix(beta); # element i, j is the beta of the stock i with the factor j 
  N = nrow(beta) # number of stocks
  stockM2 = as.numeric(stockM2); # idiosyncratic second moment of the stock
  stockM4 = as.numeric(stockM4); # idiosyncratic fourth moment of the stock  
  factorM2 = as.matrix(factorM2); #second moment of the factor
  factorM4 = as.matrix(factorM4); #fourth moment of the factor  
  S = (beta%*%factorM4)%*%(t(beta)%x%t(beta)%x%t(beta))
  # betacov
  betacov = as.numeric(beta%*%factorM2%*%t(beta))
  D = residualcokurtosisMF(NN=N, sstockM2=stockM2, sstockM4=stockM4, bbetacov=betacov)
  return( S + D )
}

cokurtosisMF = cmpfun(cokurtosisM)


## portfolio moments and the derivatives

m2 = function(w,M2)
{
  return(t(w)%*%M2%*%w)
}
derm2 = function(w,M2)
{
  return(2*M2%*%w)
}

m3 = function(w,M3)
{
  return(t(w)%*%M3%*%(w%x%w))
}
derm3 = function(w,M3)
{
  return(3*M3%*%(w%x%w))
}
m4 = function(w,M4)
{
  return(t(w)%*%M4%*%(w%x%w%x%w))
}
derm4 = function(w,M4)
{
  return(4*M4%*%(w%x%w%x%w))
}

StdDevfun = function(w,sigma){ return(  sqrt( t(w)%*%sigma%*%w  )) }

GVaRfun = function(w,alpha,mu,sigma){ return (- (t(w)%*%mu) - qnorm(alpha)*sqrt( t(w)%*%sigma%*%w  ) ) }

mVaRfun = function(w,alpha,mu,sigma,M3,M4){
  pm4 = t(w)%*%M4%*%(w%x%w%x%w) ; pm3 = t(w)%*%M3%*%(w%x%w) ; pm2 =  t(w)%*%sigma%*%w ; 
  skew = pm3 / pm2^(3/2);
  exkurt = pm4 / pm2^(2) - 3; z = qnorm(alpha);
  h = z + (1/6)*(z^2 -1)*skew
  h = h + (1/24)*(z^3 - 3*z)*exkurt - (1/36)*(2*z^3 - 5*z)*skew^2;
  return (- (t(w)%*%mu) - h*sqrt( pm2  ) ) }

resmVaRfun = function(w,alpha,mu,sigma,ressigma,M3,M4){
  pm4 = t(w)%*%M4%*%(w%x%w%x%w) ; pm3 = t(w)%*%M3%*%(w%x%w) ; pm2 =  t(w)%*%sigma%*%w ; respm2 =  t(w)%*%resSigma%*%w ;
  skew = pm3 / respm2^(3/2);
  exkurt = pm4 / respm2^(2) - 3; z = qnorm(alpha);
  h = z + (1/6)*(z^2 -1)*skew
  h = h + (1/24)*(z^3 - 3*z)*exkurt - (1/36)*(2*z^3 - 5*z)*skew^2;
  return (- (t(w)%*%mu) - h*sqrt( pm2  ) ) }

GESfun = function(w,alpha,mu,sigma,M3,M4){
  return (- (t(w)%*%mu) + dnorm(qnorm(alpha))*sqrt(t(w)%*%sigma%*%w)/alpha ) }

operMESfun = function(w,alpha,mu,sigma,M3,M4){
  pm4 = t(w)%*%M4%*%(w%x%w%x%w) ; pm3 = t(w)%*%M3%*%(w%x%w) ; pm2 =  t(w)%*%sigma%*%w ; 
  skew = pm3 / pm2^(3/2);
  exkurt = pm4 / pm2^(2) - 3; z = qnorm(alpha);
  h = z + (1/6)*(z^2 -1)*skew
  h = h + (1/24)*(z^3 - 3*z)*exkurt - (1/36)*(2*z^3 - 5*z)*skew^2;
  E = dnorm(h)
  E = E + (1/24)*(   Ipower(4,h) - 6*Ipower(2,h) + 3*dnorm(h)   )*exkurt
  E = E +  (1/6)*(   Ipower(3,h) - 3*Ipower(1,h)   )*skew;
  E = E + (1/72)*(  Ipower(6,h) -15*Ipower(4,h)+ 45*Ipower(2,h) - 15*dnorm(h) )*(skew^2)
  E = E/alpha
  return (- (t(w)%*%mu) - sqrt(pm2)*min(-E,h) )
}


Portmean = function(w,mu,precision=4)
{
  return( list(  round( t(w)%*%mu , precision) , round ( as.vector(w)*as.vector(mu) , precision ) , round( as.vector(w)*as.vector(mu)/t(w)%*%mu) , precision) )
}

Portsd =  function(w,sigma,precision=4)
{
  pm2 = m2(w,sigma)
  dpm2 = derm2(w,sigma)
  dersd = (0.5*as.vector(dpm2))/sqrt(pm2);
  contrib = dersd*as.vector(w)
  return(list(  round( sqrt(pm2) , precision ) , round( contrib , precision ) , round ( contrib/sqrt(pm2) , precision) )) 
}


PortgausVaR =  function(alpha,w,mu,sigma,precision=4){
  location = t(w)%*%mu
  pm2 = m2(w,sigma)
  dpm2 = derm2(w,sigma)
  VaR = - location - qnorm(alpha)*sqrt(pm2)
  derVaR = - as.vector(mu)- qnorm(alpha)*(0.5*as.vector(dpm2))/sqrt(pm2);
  contrib = derVaR*as.vector(w) 
  return(list( round( VaR , precision ) , round ( contrib , precision ) , round( contrib/VaR , precision) )) 
}

PortgausES =  function(alpha,w,mu,sigma,precision=4){
  location = t(w)%*%mu
  pm2 = m2(w,sigma)
  dpm2 = derm2(w,sigma)
  ES = - location + dnorm(qnorm(alpha))*sqrt(pm2)/alpha
  derES = - as.vector(mu) + (1/alpha)*dnorm(qnorm(alpha))*(0.5*as.vector(dpm2))/sqrt(pm2);
  contrib = as.vector(w)*derES;
  return(list( round( ES , precision ) , round( contrib , precision) , round( contrib/ES , precision) )) 
}

PortSkew =  function(w,sigma,M3)
{
  pm2 = m2(w,sigma)
  pm3 = m3(w,M3)
  skew = pm3 / pm2^(3/2);
  return( skew )
}

PortKurt =  function(w,sigma,M4)
{
  pm2 = m2(w,sigma)
  pm4 = m4(w,M4)
  kurt = pm4 / pm2^(2) ;
  return( kurt )
}

PortMVaR =  function(alpha,w,mu,sigma,M3,M4,precision=4)
{
  z = qnorm(alpha)
  location = t(w)%*%mu
  pm2 = m2(w,sigma)
  dpm2 = as.vector( derm2(w,sigma) )
  pm3 = m3(w,M3)
  dpm3 = as.vector( derm3(w,M3) )
  pm4 = m4(w,M4)
  dpm4 = as.vector( derm4(w,M4) )
  
  skew = pm3 / pm2^(3/2);
  exkurt = pm4 / pm2^(2) - 3;
  
  derskew = ( 2*(pm2^(3/2))*dpm3 - 3*pm3*sqrt(pm2)*dpm2 )/(2*pm2^3)
  derexkurt = ( (pm2)*dpm4 - 2*pm4*dpm2    )/(pm2^3)
  
  h = z + (1/6)*(z^2 -1)*skew 
  h = h + (1/24)*(z^3 - 3*z)*exkurt - (1/36)*(2*z^3 - 5*z)*skew^2; 
  
  MVaR = - location - h*sqrt(pm2)
  
  derGausVaR = - as.vector(mu)- qnorm(alpha)*(0.5*as.vector(dpm2))/sqrt(pm2);
  derMVaR = derGausVaR + (0.5*dpm2/sqrt(pm2))*( -(1/6)*(z^2 -1)*skew  - (1/24)*(z^3 - 3*z)*exkurt + (1/36)*(2*z^3 - 5*z)*skew^2 )
  derMVaR = derMVaR + sqrt(pm2)*( -(1/6)*(z^2 -1)*derskew  - (1/24)*(z^3 - 3*z)*derexkurt + (1/36)*(2*z^3 - 5*z)*2*skew*derskew  )
  contrib = as.vector(w)*as.vector(derMVaR)
  return(list(  round( MVaR , precision) , round( contrib , precision ), round (contrib/MVaR , precision ) ) ) 
}

derIpower = function(power,h)
{
  
  fullprod = 1;
  
  if( (power%%2)==0 ) #even number: number mod is zero
  {
    pstar = power/2;
    for(j in c(1:pstar))
    {
      fullprod = fullprod*(2*j)
    }
    I = -fullprod*h*dnorm(h);
    
    for(i in c(1:pstar) )
    { 
      prod = 1;
      for(j in c(1:i) )
      {
        prod = prod*(2*j)
      }
      I = I + (fullprod/prod)*(h^(2*i-1))*(2*i-h^2)*dnorm(h)
    }
  }else{
    pstar = (power-1)/2
    for(j in c(0:pstar) )
    {
      fullprod = fullprod*( (2*j)+1 )
    }
    I = -fullprod*dnorm(h);
    
    for(i in c(0:pstar) )
    { 
      prod = 1;
      for(j in c(0:i) )
      {
        prod = prod*( (2*j) + 1 )
      }
      I = I + (fullprod/prod)*(h^(2*i)*(2*i+1-h^2) )*dnorm(h)
    }
  }
  return(I)
}


PortMES =  function(alpha,w,mu,sigma,M3,M4,precision=4)
{
  z = qnorm(alpha)
  location = t(w)%*%mu
  pm2 = m2(w,sigma)
  dpm2 = as.vector( derm2(w,sigma) )
  pm3 = m3(w,M3)
  dpm3 = as.vector( derm3(w,M3) )
  pm4 = m4(w,M4)
  dpm4 = as.vector( derm4(w,M4) )
  
  skew = pm3 / pm2^(3/2);
  exkurt = pm4 / pm2^(2) - 3;
  
  derskew = ( 2*(pm2^(3/2))*dpm3 - 3*pm3*sqrt(pm2)*dpm2 )/(2*pm2^3)
  derexkurt = ( (pm2)*dpm4 - 2*pm4*dpm2    )/(pm2^3)
  
  h = z + (1/6)*(z^2 -1)*skew 
  h = h + (1/24)*(z^3 - 3*z)*exkurt - (1/36)*(2*z^3 - 5*z)*skew^2; 
  
  derh = (1/6)*(z^2 -1)*derskew + (1/24)*(z^3 - 3*z)*derexkurt - (1/18)*(2*z^3 - 5*z)*skew*derskew
  
  E = dnorm(h)
  E = E + (1/24)*(   Ipower(4,h) - 6*Ipower(2,h) + 3*dnorm(h)   )*exkurt
  E = E +  (1/6)*(   Ipower(3,h) - 3*Ipower(1,h)   )*skew;
  E = E + (1/72)*(  Ipower(6,h) -15*Ipower(4,h)+ 45*Ipower(2,h) - 15*dnorm(h) )*(skew^2)
  E = E/alpha
  MES = - location + sqrt(pm2)*E
  
  derMES = -mu + 0.5*(dpm2/sqrt(pm2))*E
  derE = (1/24)*(   Ipower(4,h) - 6*Ipower(2,h) + 3*dnorm(h)   )*derexkurt 
  derE = derE +  (1/6)*(   Ipower(3,h) - 3*Ipower(1,h)   )*derskew
  derE = derE + (1/36)*(  Ipower(6,h) -15*Ipower(4,h)+ 45*Ipower(2,h) - 15*dnorm(h) )*skew*derskew
  X = -h*dnorm(h) + (1/24)*(  derIpower(4,h) - 6*derIpower(2,h) -3*h*dnorm(h)  )*exkurt 
  X = X + (1/6)*( derIpower(3,h) - 3*derIpower(1,h) )*skew 
  X = X + (1/72)*( derIpower(6,h) - 15*derIpower(4,h) + 45*derIpower(2,h) + 15*h*dnorm(h)  )*skew^2
  derE = derE+derh*X  # X is a scalar
  derE = derE/alpha
  derMES = derMES + sqrt(pm2)*derE
  contrib = as.vector(w)*as.vector(derMES)
  return(list(  round( MES , precision ) , round( contrib , precision ), round( contrib/MES, precision )) ) 
}

derMES =  function(alpha,w,mu,sigma,M3,M4,precision=4)
{
  z = qnorm(alpha)
  location = t(w)%*%mu
  pm2 = m2(w,sigma)
  dpm2 = as.vector( derm2(w,sigma) )
  pm3 = m3(w,M3)
  dpm3 = as.vector( derm3(w,M3) )
  pm4 = m4(w,M4)
  dpm4 = as.vector( derm4(w,M4) )
  
  skew = pm3 / pm2^(3/2);
  exkurt = pm4 / pm2^(2) - 3;
  
  derskew = ( 2*(pm2^(3/2))*dpm3 - 3*pm3*sqrt(pm2)*dpm2 )/(2*pm2^3)
  derexkurt = ( (pm2)*dpm4 - 2*pm4*dpm2    )/(pm2^3)
  
  h = z + (1/6)*(z^2 -1)*skew 
  h = h + (1/24)*(z^3 - 3*z)*exkurt - (1/36)*(2*z^3 - 5*z)*skew^2; 
  
  derh = (1/6)*(z^2 -1)*derskew + (1/24)*(z^3 - 3*z)*derexkurt - (1/18)*(2*z^3 - 5*z)*skew*derskew
  
  E = dnorm(h)
  E = E + (1/24)*(   Ipower(4,h) - 6*Ipower(2,h) + 3*dnorm(h)   )*exkurt
  E = E +  (1/6)*(   Ipower(3,h) - 3*Ipower(1,h)   )*skew;
  E = E + (1/72)*(  Ipower(6,h) -15*Ipower(4,h)+ 45*Ipower(2,h) - 15*dnorm(h) )*(skew^2)
  E = E/alpha
  MES = - location + sqrt(pm2)*E
  
  derMES = -mu + 0.5*(dpm2/sqrt(pm2))*E
  derE = (1/24)*(   Ipower(4,h) - 6*Ipower(2,h) + 3*dnorm(h)   )*derexkurt 
  derE = derE +  (1/6)*(   Ipower(3,h) - 3*Ipower(1,h)   )*derskew
  derE = derE + (1/36)*(  Ipower(6,h) -15*Ipower(4,h)+ 45*Ipower(2,h) - 15*dnorm(h) )*skew*derskew
  X = -h*dnorm(h) + (1/24)*(  derIpower(4,h) - 6*derIpower(2,h) -3*h*dnorm(h)  )*exkurt 
  X = X + (1/6)*( derIpower(3,h) - 3*derIpower(1,h) )*skew 
  X = X + (1/72)*( derIpower(6,h) - 15*derIpower(4,h) + 45*derIpower(2,h) + 15*h*dnorm(h)  )*skew^2
  derE = derE+derh*X  # X is a scalar
  derE = derE/alpha
  derMES = derMES + sqrt(pm2)*derE
  return(derMES)
}


operPortMES =  function(alpha,w,mu,sigma,M3,M4)
{
  z = qnorm(alpha)
  location = t(w)%*%mu
  pm2 = m2(w,sigma)
  dpm2 = as.vector( derm2(w,sigma) )
  pm3 = m3(w,M3)
  dpm3 = as.vector( derm3(w,M3) )
  pm4 = m4(w,M4)
  dpm4 = as.vector( derm4(w,M4) )
  
  skew = pm3 / pm2^(3/2);
  exkurt = pm4 / pm2^(2) - 3;
  
  derskew = ( 2*(pm2^(3/2))*dpm3 - 3*pm3*sqrt(pm2)*dpm2 )/(2*pm2^3)
  derexkurt = ( (pm2)*dpm4 - 2*pm4*dpm2    )/(pm2^3)
  
  h = z + (1/6)*(z^2 -1)*skew 
  h = h + (1/24)*(z^3 - 3*z)*exkurt - (1/36)*(2*z^3 - 5*z)*skew^2; 
  I1 = Ipower(1,h); I2 = Ipower(2,h); I3 = Ipower(3,h); I4 = Ipower(4,h);  I6 = Ipower(6,h); 
  
  derh = (1/6)*(z^2 -1)*derskew + (1/24)*(z^3 - 3*z)*derexkurt - (1/18)*(2*z^3 - 5*z)*skew*derskew
  
  E = dnorm(h)
  E = E + (1/24)*(   I4 - 6*I2 + 3*dnorm(h)   )*exkurt
  E = E +  (1/6)*(   I3 - 3*I1   )*skew;
  E = E + (1/72)*( I6 -15*I4+ 45*I2 - 15*dnorm(h) )*(skew^2)
  E = E/alpha
  
  MES = - location - sqrt(pm2)*min(-E,h)
  
  if(-E<=h){
    derMES = -mu + 0.5*(dpm2/sqrt(pm2))*E
    derE = (1/24)*(   I4 - 6*I2 + 3*dnorm(h)   )*derexkurt 
    derE = derE +  (1/6)*(   I3 - 3*I1   )*derskew
    derE = derE + (1/36)*(  I6 -15*I4 + 45*I2 - 15*dnorm(h) )*skew*derskew
    X = -h*dnorm(h) + (1/24)*(  derIpower(4,h) - 6*derIpower(2,h) -3*h*dnorm(h)  )*exkurt 
    X = X + (1/6)*( derIpower(3,h) - 3*derIpower(1,h) )*skew 
    X = X + (1/72)*( derIpower(6,h) - 15*derIpower(4,h) + 45*derIpower(2,h) + 15*h*dnorm(h)  )*skew^2
    derE = derE+derh*X  # X is a scalar
    derE = derE/alpha
    derMES = derMES + sqrt(pm2)*derE
  }else{
    derMES = -mu - 0.5*(dpm2/sqrt(pm2))*h - sqrt(pm2)*derh ; 
  }
  contrib = as.vector(w)*as.vector(derMES)
  return(list(  MES ,  contrib , contrib/MES ) ) 
}

deroperPortMES =  function(alpha,w,mu,sigma,M3,M4)
{
  z = qnorm(alpha)
  location = t(w)%*%mu
  pm2 = m2(w,sigma)
  dpm2 = as.vector( derm2(w,sigma) )
  pm3 = m3(w,M3)
  dpm3 = as.vector( derm3(w,M3) )
  pm4 = m4(w,M4)
  dpm4 = as.vector( derm4(w,M4) )
  
  skew = pm3 / pm2^(3/2);
  exkurt = pm4 / pm2^(2) - 3;
  
  derskew = ( 2*(pm2^(3/2))*dpm3 - 3*pm3*sqrt(pm2)*dpm2 )/(2*pm2^3)
  derexkurt = ( (pm2)*dpm4 - 2*pm4*dpm2    )/(pm2^3)
  
  h = z + (1/6)*(z^2 -1)*skew 
  h = h + (1/24)*(z^3 - 3*z)*exkurt - (1/36)*(2*z^3 - 5*z)*skew^2; 
  I1 = Ipower(1,h); I2 = Ipower(2,h); I3 = Ipower(3,h); I4 = Ipower(4,h);  I6 = Ipower(6,h); 
  
  derh = (1/6)*(z^2 -1)*derskew + (1/24)*(z^3 - 3*z)*derexkurt - (1/18)*(2*z^3 - 5*z)*skew*derskew
  
  E = dnorm(h)
  E = E + (1/24)*(   I4 - 6*I2 + 3*dnorm(h)   )*exkurt
  E = E +  (1/6)*(   I3 - 3*I1   )*skew;
  E = E + (1/72)*( I6 -15*I4+ 45*I2 - 15*dnorm(h) )*(skew^2)
  E = E/alpha
  
  MES = - location - sqrt(pm2)*min(-E,h)
  
  if(-E<=h){
    derMES = -mu + 0.5*(dpm2/sqrt(pm2))*E
    derE = (1/24)*(   I4 - 6*I2 + 3*dnorm(h)   )*derexkurt 
    derE = derE +  (1/6)*(   I3 - 3*I1   )*derskew
    derE = derE + (1/36)*(  I6 -15*I4 + 45*I2 - 15*dnorm(h) )*skew*derskew
    X = -h*dnorm(h) + (1/24)*(  derIpower(4,h) - 6*derIpower(2,h) -3*h*dnorm(h)  )*exkurt 
    X = X + (1/6)*( derIpower(3,h) - 3*derIpower(1,h) )*skew 
    X = X + (1/72)*( derIpower(6,h) - 15*derIpower(4,h) + 45*derIpower(2,h) + 15*h*dnorm(h)  )*skew^2
    derE = derE+derh*X  # X is a scalar
    derE = derE/alpha
    derMES = derMES + sqrt(pm2)*derE
  }else{
    derMES = -mu - 0.5*(dpm2/sqrt(pm2))*h - sqrt(pm2)*derh ; 
  }
  return( derMES )
}

centeredmoment = function(series,power)
{
  location = mean(series);
  out = sum( (series-location)^power  )/length(series);
  return(out);
}

operMES =  function(series,alpha,r)
{
  z = qnorm(alpha)
  location = mean(series);
  m2 = centeredmoment(series,2)
  m3 = centeredmoment(series,3)
  m4 = centeredmoment(series,4)
  skew = m3 / m2^(3/2);
  exkurt = m4 / m2^(2) - 3;
  
  h = z + (1/6)*(z^2 -1)*skew 
  if(r==2){ h = h + (1/24)*(z^3 - 3*z)*exkurt - (1/36)*(2*z^3 - 5*z)*skew^2};
  
  MES = dnorm(h)
  MES = MES + (1/24)*(   Ipower(4,h) - 6*Ipower(2,h) + 3*dnorm(h)   )*exkurt
  MES = MES +  (1/6)*(   Ipower(3,h) - 3*Ipower(1,h)   )*skew;
  MES = MES + (1/72)*(  Ipower(6,h) -15*Ipower(4,h)+ 45*Ipower(2,h) - 15*dnorm(h) )*(skew^2)
  MES = - location - (sqrt(m2))*min( -MES/alpha , h )
  return(MES)
}



Ipower = function(power,h)
{
  
  fullprod = 1;
  
  if( (power%%2)==0 ) #even number: number mod is zero
  {
    pstar = power/2;
    for(j in c(1:pstar))
    {
      fullprod = fullprod*(2*j)
    }
    I = fullprod*dnorm(h);
    
    for(i in c(1:pstar) )
    { 
      prod = 1;
      for(j in c(1:i) )
      {
        prod = prod*(2*j)
      }
      I = I + (fullprod/prod)*(h^(2*i))*dnorm(h)
    }
  }
  else{
    pstar = (power-1)/2
    for(j in c(0:pstar) )
    {
      fullprod = fullprod*( (2*j)+1 )
    }
    I = -fullprod*pnorm(h);
    
    for(i in c(0:pstar) )
    { 
      prod = 1;
      for(j in c(0:i) )
      {
        prod = prod*( (2*j) + 1 )
      }
      I = I + (fullprod/prod)*(h^(  (2*i) + 1))*dnorm(h)
    }
  }
  return(I)
}

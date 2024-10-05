# Sebastian Fossati
# 1/2013


ur.test <- function(yt,trend,method,penalty="BIC",kmax=NULL,kmin=NULL){

#	Implements unit root tests discussed in Ng, S. 
#	and P. Perron (2001): "Lag Length Selection and 
#	the Construction of Unit Root Tests with Good 
#	Size and Power", Econometrica, 69, 6, 1519-1554.
#	
#	Inputs: 
#			yt -> univariate time series
#           trend -> the order of the polynomial trend
#               a) trend=c for constant
#               b) trend=ct for constant and linear trend
#			method -> unit root test to be constructed
#				a) method=adf.ols for augmented Dickey-Fuller test
#				b) method=adf.gls for Dickey-Fuller GLS test
#				c) method=pt for feasible point optimal test
#           penalty -> penalty for lag length selection
#				a) penalty -> MAIC for MAIC 
#				b) penalty -> AIC for AIC 
#				c) penalty -> BIC for BIC 
#           kmax -> maximum number of lags
#           kmin -> minimum number of lags
#	
#	Output: 
#			stat: 	the test statistic
#			cv: 	the corresponding .01, .05, and .10 critical values
#   		krule:	the order of the autoregression selected by 
#			the data dependent rule (used in constructing 
#			the spectral density estimator at frequency zero)
#			phi:	estimate of phi
	
	# set data as matrix!
	yt <- as.matrix(yt)
	nt <- nrow(yt)
	if (is.null(kmax)){kmax <- floor(12*(nt/100)^(.25))} 
	if (is.null(kmin)){kmin <- 0} 
	# get deterministic components
	if(trend=="c"){z <- matrix(1,nr=nt); cbar <- -7}
	if(trend=="ct"){
		z <- cbind(matrix(1,nr=nt),matrix(seq(1:nt)))
		cbar <- -13.5
		}
	# transform the data
	tmp <- glsd(yt,z,cbar)
	yd <- tmp$yd
	phi <- solve(t(yd[1:(nt-1)])%*%yd[1:(nt-1)])%*%
					t(yd[1:(nt-1)])%*%yd[2:nt]
	ssra <- tmp$ssr
	# get optimal number of lags
	if(kmax==kmin){krule <- kmax}
	else{krule <- s2ar(yd,penalty,kmax,kmin)}
	
	# adf-ols
	if(method=="adf.ols"){
		# construct adf-ols
		adf.ols <- adf.ols(yt,z,krule)
		# 5% critical values
		if(trend=="c"){cv <- c(-3.43,-2.86,-2.57)}
		if(trend=="ct"){cv <- c(-3.96,-3.41,-3.13)}
		listss <- list(method=method,stat=adf.ols$tstat,cv=cv,
								penalty=penalty,k=krule,phi=phi)
		}
	# adf-gls of ERS (1996)
	if(method=="adf.gls"){
		# construct adf-ols
		adf.gls <- adfk(yd,krule)
		# 5% critical values
		if(trend=="c"){cv <- c(-2.58,-1.98,-1.62)}
		if(trend=="ct"){cv <- c(-3.42,-2.91,-2.62)}
		listss <- list(method=method,stat=adf.gls$tstat,cv=cv,
								penalty=penalty,k=krule,phi=phi)
		}
	# point optimal test of ERS (1996)
	if(method=="pt"){
		# construct adf-gls for sar
		adf.gls <- adfk(yd,krule)
		sar <- adf.gls$s2vec
		# construct pt
		tmp <- glsd(yt,z,0)
		ssr1 <- tmp$ssr
		pt <- (ssra - (1+cbar/nt)*ssr1)/sar
		# 5% critical values
		if(trend=="c"){cv <- c(1.78,3.17,4.45)}
		if(trend=="ct"){cv <- c(4.03,5.48,6.67)}
		listss <- list(method=method,stat=pt,cv=cv,
							penalty=penalty,k=krule,phi=phi)
		}		

	return( listss )
	
	}

print.ur.test <- function(test){
	method <- test$method
	# print results
	cat(sprintf("\n Test for Unit Root: "))
	if(method=="adf.ols"){cat(sprintf("ADF Test"))}
	if(method=="adf.gls"){cat(sprintf("ADF Test with GLS Detrending"))}
	if(method=="pt"){cat(sprintf("ERS Test"))}
	cat(sprintf("\n\n Null Hypothesis: there is a unit root"))
	cat(sprintf("\n Test Statistic: %.3f",test$stat))
	cat(sprintf("\n Critical Values (.01,.05,.10): %.2f %.2f %.2f",
									test$cv[1],test$cv[2],test$cv[3]))		
	cat(sprintf("\n\n Lag Order Selection Rule: %s",test$penalty))
	cat(sprintf("\n Selected Lag Order: %1.0f",test$k))
	cat(sprintf("\n Estimated Coefficient: %.4f \n\n",test$phi))
}

adfk <- function(yt,kstar){
	reg <- lagn(yt,1)
	dyt <- diffk(yt,1)
	if(kstar>0){for(i in 1:kstar){reg <- cbind(reg,lagn(dyt,i))}}
	dyt <- trimr(dyt,kstar+1,0)
	reg <- trimr(reg,kstar+1,0)
	b <- solve(t(reg)%*%reg)%*%t(reg)%*%dyt
	ee <- dyt - reg%*%b
	nef <- nrow(dyt)-kstar-1
	s2e <- t(ee)%*%ee/nef
	xx <- solve(t(reg)%*%reg)
	sre <- xx[1,1]*s2e
	sumb <- 0
	if(kstar>0){sumb <- sum(b[2:(kstar+1)])}
	return(list(tstat=b[1]/sqrt(sre),rho=b[1]+1,s2vec=s2e/(1-sumb)^2))
	}

adf.ols <- function(yt,z,kstar){
	reg <- cbind(lagn(yt,1),z)
	dyt <- diffk(yt,1)
	if(kstar>0){for(i in 1:kstar){reg <- cbind(reg,lagn(dyt,i))}}
	dyt <- trimr(dyt,kstar+1,0)
	reg <- trimr(reg,kstar+1,0)
	b <- solve(t(reg)%*%reg)%*%t(reg)%*%dyt
	ee <- dyt - reg%*%b
	nef <- nrow(dyt)-kstar-1-ncol(z)
	s2e <- t(ee)%*%ee/nef
	xx <- solve(t(reg)%*%reg)
	sre <- xx[1,1]*s2e
	return(list(tstat=b[1]/sqrt(sre)))
	}

s2ar <- function(yt,penalty,kmax,kmin){
	tau <- rep(0,kmax+1)
	s2e <- 999999*rep(1,kmax+1)
	dyt <- diffk(yt,1)
	reg <- lagn(yt,1)
	if(kmax>0){for(i in 1:kmax){reg <- cbind(reg,lagn(dyt,i))}}
	dyt <- trimr(dyt,kmax+1,0)
	reg <- trimr(reg,kmax+1,0)
	sumy <- t(reg[,1])%*%reg[,1]
	nef <- nrow(dyt)
	for(k in kmin:kmax){
		b <- solve(t(reg[,1:(k+1)])%*%reg[,1:(k+1)])%*%
									t(reg[,1:(k+1)])%*%dyt
		e <- dyt - reg[,1:(k+1)]%*%b
		s2e[k+1] <- t(e)%*%e/nef
		tau[k+1] <- b[1]*b[1]*sumy/s2e[k+1]
		}
	kk <- seq(0,kmax)
	if(penalty=="AIC"){ic <- log(s2e)+2*kk/nef}
	else if(penalty=="BIC"){ic <- log(s2e)+log(nef)*kk/nef}
	else {ic <- log(s2e)+2*(kk+tau)/nef}
	tmp <- which.min(ic)
	return(tmp-1)
	}

glsd <- function(yt,z,cbar){
	nt <- nrow(yt)
	abar <- 1 + cbar/nt
	ya <- matrix(0,nt,1)
	za <- matrix(0,nt,ncol(z))
	ya[1,1] <- yt[1,1]
	za[1,] <- z[1,]
	ya[2:nt,1] <- yt[2:nt,1]-abar*yt[1:(nt-1),1]
	za[2:nt,] <- z[2:nt,]-abar*z[1:(nt-1),]
	# construct gls detrended series
	bhat <- solve(t(za)%*%za)%*%t(za)%*%ya
	yd <- yt - z%*%bhat
	ssr <- t(ya-za%*%bhat)%*%(ya-za%*%bhat)
	return(list(yd=yd,ssr=ssr))
	}
	
olsd <- function(yt,z){
	# construct ols detrended series
	bhat <- solve(t(z)%*%z)%*%t(z)%*%yt
	yd <- yt - z%*%bhat
	return(yd=yd)
	}

diffk <- function(x,k){
	if(k==0){return(x)}
	else{return(rbind(matrix(0,k,ncol(x)),diff(x,lag=k)))}
	}

lagn <- function(x,n){
	if(n<0){return(rbind(trimr(x,abs(n),0),
						matrix(0,abs(n),ncol(x))))}
	else{return(rbind(matrix(0,n,ncol(x)),trimr(x,0,n)))}
	}

trimr <- function(x,n1,n2){
	return(as.matrix(x[(n1+1):(nrow(x)-n2),],nc=ncol(x)))
	}


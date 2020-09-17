
##############  FUNCTIONS FOR COMPOSITIONAL DATA ANALYSIS ###########
### most function names are those of "compositions" package with an "m"
### as initial letter, thus avoiding overlapping.
### 
### this minipackage was initially written for substituting "compositions"
### in a period in which "compositions" was not available in CRAN.
### Now it can be used to avoid the " burden of the object directed
### programming" and intrincacies of "compositions" in simple
### scripts
###    This version was completed in May 2017 by JJ Egozcue
###    Revised in August 2017
#####################################################################
# This mini-library contains the following functions:
#   mbuildcontrast(W=c(1,-1))  build contrasts coded as +1,-1,0 code
#                    its main use is constructig basis contrast matrices
#                    from coded sequential binary partitions
#   mcls(x)     closes to 1 the rows of x
#   milr(x, V)  assign ilr coordinates to the compositions in the
#               rows of x associated with V, the contrast matrix basis
#   milrInv(z, V) computes the ilr inverse following the contrast matrix
#               basis V.
#   mclr(x)  computes the clr of x containing compositions by rows
#   mvariation(x,optnor="none")  computes the variation matrix of x 
#               containingcompositions by rows. Normalization allowed.
#   mcenmat(X,optcen="row") arithmetically centers the matrix X by rows
#                 optcen="row" or columns optcen="col"
#   mBPPOP(x,xpop=1,choice=c(1,2),biscale=1,punch=1,colist=1:10,
#             optpdf=0,
#             colray="red",cextext=1,lwdray=1,pchpoint=1)
#                 computes the centered  clr of x and its svd
#                 plots a biplot in which possible populations
#                 are differently colored
#   mcountsprop(X) estimate multinomial probabilities or proportions 
#                from observed counts, possibly containing zeros.
#   mvariatzeros(X,toleig=0.000001) computes variation matrix using all
#                non zero data
#   mPBclustvar(x) compute the contrast matrix of Principal Balances
#                approximated by hierachical clustering of variables 
#                (Wards method)
#   mmerge2sign(Merge) from a $merge object of a cluster analysis
#                builds asign code for a SBP
#   mTablediscrimErrors(ogroup,predgroup) computes the table of errors
#               in a discriminant analysis. Center probabilities against
#               observed groups.
#
# function mcheckfunctions() carries out check of the previous functions
#
#################################################################
# If desired these functions can be inserted in a workspace
# inserting the next source statement:
#source("~/JJE/projects/R-lib-JJE/proper-libs/mini_compositions_lib.R")
####################################################################

####################################
##### function mbuildcontrast ###
##### Builds the contrast matrix ###
##### from a SBP code (signs) ######
##### for an ilr-transform    ######
# contrasts are defined by columns
# intput W is a matrix with D rows
#     and a undefined number of columns
#     containing 1,-1,0
# W also can be a vector e.g. c(1,1,-1,0)
# The matrix W does not need to correspond
# to a sequential binary partition,
# or to have strictly D-1 columns. Each
# row define a balance and it is computed
# independently of other columns.
####################################
## J. J. Egozcue, May 2014, May 2017
####################################
mbuildcontrast <-function(W=c(1,-1)){
  if(length(W)<2){
  	print("improper dimension") 
    return(W)
  }
	
   W = as.matrix(W)
   D = nrow(W)
   nc = ncol(W)
   isPos = (W>0)
   isNeg = (W<0)
   nPos=colSums(isPos)
   nNeg=colSums(isNeg)
   valPos=sqrt(nNeg/(nPos*(nPos+nNeg)))
   valNeg=-sqrt(nPos/(nNeg*(nPos+nNeg)))
   WW = isPos * outer(rep(1,D),valPos) + 
   	    isNeg * outer(rep(1,D),valNeg)
   return(WW)
}


###### function mcls #######
###### closes rows of positive
# components to sum 1 #####
###########################
# JJ Egozcue, May 2017
###########################

mcls <- function(x){
	if(is.null(dim(x)) ){x=t(as.matrix(x))}
  csx = rowSums(x)
  xx = sweep(x, 1, csx, "/")
  return(xx)
}
####################################



###### function milr ##############
## computes the ilr transformation
## of x  (n,D)
## using the contrast matrix V
###########################
# JJ Egozcue, May 2016
###########################
milr <- function(x, V){
	if(is.null(dim(x)) ){x=t(as.matrix(x))}
  xilr = log(as.matrix(x)) %*% V
  return(xilr)
}
####################################

##### function milrInv ############
## computes ilr inverse of rows of z
## using the contrast matrix V
###########################
# JJ Egozcue, May 2017
###########################
milrInv <- function(z, V){
  if(is.null(dim(z)) ){z=t(as.matrix(z))}
  zcls = mcls(exp(z %*% t(V)))
  return(zcls)
}

##### function mclr ##############
## computes the clr of x by rows
###########################
# JJ Egozcue, May 2017
###########################
mclr <- function(x){
	if(is.null(dim(x)) ){x=t(as.matrix(x))}
  logx = log(x)
  rs = outer(rowSums(logx),rep((1/ncol(x)),length=(ncol(x))))
  xclr=logx-rs
  return(xclr)
}
###############################

###### function mvariation #######
## computes the variation matrix
## of a compositional data matrix x
## There are three options of normalization
## optnor="none"  variation matrix as is
## optnor="minassoc" variation over minimum
##        associated matrix
## optnor="linear01", rational 0,1 transform
## These normalizations are
## "minassoc" : t_ij = (d-1) v_ij / 2 totvar
## "linear01" : t_ij = totvar / [totvar + ((d-1)/2) v_ij]
## The value of v_ij which distributes uniformly total
## variance is v_u = 2 totvar/(d-1)
## characteristic values
##  "no normalization"   "minassoc"  "linear01"
##     0                    0             1
##     v_u                  1             0.5
##     Inf                  Inf           0
#########################################################
# JJ Egozcue, May 2017
###########################
mvariation <- function(x,optnor="none"){
	if(is.null(dim(x)) ){x=t(as.matrix(x))}
  d = ncol(x)
  xclr = mclr(x)
  co = var(xclr)
  va = diag(co)
  co1 = matrix(rep(va,each=d),ncol = d)
  co2 = matrix(rep(va,times=d),ncol = d)
  varia = -2 * co + co1 + co2
  if(optnor=="none"){  return(varia) }
  totvar = sum(varia)/(2*d)
  varian= varia*((d-1)/(2*totvar))
  if(optnor=="minassoc"){
  	return(varian)
  }
  if(optnor=="linear01"){
  	varian1 = totvar/(varia*((d-1)/2)+totvar)
  	return(varian1)
  }
}
##################################

#### mcenmat function ############
# arithmetic centering of a matrix 
# by rows or columns
# options "row" default, "col" 
##################################
# JJ Egozcue May 2017
##################################
mcenmat <- function(X,optcen="row"){
	if(optcen=="row"){
	 Xcen = X - outer( rowMeans(X),rep(1,length=ncol(X)) )  
	}
  if(optcen=="col"){
   Xcen = X - outer( rep(1,length=nrow(X)), colMeans(X))
  }
	return(Xcen)
}
##################################

######## FUNCTION mBPPOP ####################
# BiPlot, colored by POPulation 
# (not in "compositions" package)
# Programmed by J.J. Egozcue (2014)
##### draws a CoDa-biplot with data coming from 
# coded populations
# carries out clr of data set
# centres the clr
# carries out svd of centred clr
# plots biplot with data colored by pop (population)
##### input:
# x compositional data by rows (matrix)
# xpop factor indicating the population of each row in x
# biscale    = 1 covariance biplot
#            = 0 form biplot
# punch      = 0 do not plot data points
#            = 1 plot symbols for data points
#            = 2 plot numbers for data points
# choice[1:2] = PC's to be plotted, eg c(1,2), c(1,3)
# colist a sequence for colors of populations
# optpdf = 1  prints on pdf
#### output: a list containing 
# the svd matrices: U, V, and singular values in D
# explined variance in the biplot explvar
# additionally
# a pdf file contining the biplot is printed in wd
###################################################

mBPPOP <- function(x,xpop=1,choice=c(1,2),
                  biscale=1,punch=1,colist=1:10,
                  optpdf=0,
                  colray="red",cextext=1,lwdray=1,pchpoint=1){
  numxpop = as.numeric(xpop)
  colpoint = colist[numxpop]
  # clr of x
  logx=log(x)
  rs = outer(rowSums(logx)/ncol(logx),rep(1,length=ncol(logx)))
  xclr= logx-rs
  # centring xclr
  cxclr = xclr-
    outer(rep(1,length=nrow(xclr)),(apply(xclr,2,"sum")))/nrow(xclr)  
  
  # svd (cxclr)
  SVDxclr = svd(cxclr)
  U = SVDxclr$u
  V = SVDxclr$v
  D = SVDxclr$d
  # scores and loadings
  ## covariance biplot
  if(biscale==1){
    ld=t(diag(D)%*%t(V))/sqrt(nrow(cxclr))
    mainT="covariance biplot"
    fileT="covBPPOP"
    mld=max(abs(ld[,choice]))
    msc=max(abs(U[,choice]))
    sc=U *(mld/msc)  # scaling scores
    xylimit=c(-mld,mld)
  }
  
  ## form biplot
  ## scaling: unit norm of V-vectors
  if(biscale==0){
    sc=U%*%diag(D) 
    ld=V
    mainT="form biplot"
    fileT="formBPPOP"
    mld = max(abs(ld[,choice]))     # scaling basis vectors
    msc = max(abs(sc[,choice]))
    sc = sc*(mld/msc)
    xylimit = c(-mld,mld)
  }
  
  # numeric output
  variances = D^2
  totvar=sum(variances)/(nrow(x)-1)
  explvar = (variances[choice[1]]+variances[choice[2]])/sum(variances)
  # names
  #  clrnames=paste("clr.",colnames(x),sep="")
  clrnames=colnames(x)
  if(choice[1] == 1){
    expl=100*variances[choice[1]]/sum(variances)
    xlabel=paste("first axis,", " var% ",format(expl,nsmall=1,digits=3),sep="")}
  if(choice[1] == 2){
    expl=100*variances[choice[1]]/sum(variances)
    xlabel=paste("second axis,", " var% ",format(expl,nsmall=1,digits=3),sep="")}
  if(choice[1] == 3){
    expl=100*variances[choice[1]]/sum(variances)
    xlabel=paste("third axis,", " var% ",format(expl,nsmall=1,digits=3),sep="")}
  if(choice[1] >= 4){
    expl=100*variances[choice[1]]/sum(variances)
    xlabel=paste(paste(choice[1],"th axis",sep=""), " var% ",format(expl,nsmall=1,digits=3),sep="")}
  if(choice[2] == 1){
    expl=100*variances[choice[2]]/sum(variances)
    ylabel=paste("first axis,"," var% ",format(expl,nsmall=1,digits=3),sep="")}
  if(choice[2] == 2){
    expl=100*variances[choice[2]]/sum(variances)
    ylabel=paste("second axis,"," var% ",format(expl,nsmall=1,digits=3),sep="")}
  if(choice[2] == 3){
    expl=100*variances[choice[2]]/sum(variances)
    ylabel=paste("third axis,"," var% ",format(expl,nsmall=1,digits=3),sep="")}
  if(choice[2] >= 4){
    expl=100*variances[choice[2]]/sum(variances)
    ylabel=paste(paste(choice[2],"th axis",sep=""), " var% ",format(expl,nsmall=1,digits=3),sep="")}
  
  if(punch==0){pun="n"}
  if(punch==1){pun="p"}
  if(punch==2){pun="n"}
  
  # pdf output
  filenam=paste(fileT,".pdf",sep="")
  if(optpdf==1){
    pdf(filenam, width=5, height=5, fam="Times")}
  
  plot(sc[,choice],col=colpoint,type=pun,cex=0.8,asp=1,
       xlim=xylimit,ylim=xylimit,main=mainT,
       xlab=xlabel,ylab=ylabel,pch=pchpoint)
  #  this is for changing punch in place of color
  #  plot(sc[,choice],col="black",type=pun,cex=0.8,asp=1,
  #       pch=colpoint,
  #       xlim=xylimit,ylim=xylimit,main=mainT,
  #       xlab=xlabel,ylab=ylabel)
  
  if(punch==2){
#    text(sc[,choice],labels=(1:nrow(sc)),col=colpoint,cex=0.8)	    
  	text(sc[,choice],labels=pchpoint,col=colpoint,cex=0.8)
  }
  for(i in 1:ncol(x)){
    xx=rbind(c(0,0),ld[i,choice])
    lines(xx, col=colray,type="l",lwd=lwdray)
    xtext = ld[i,choice[1]]
    ytext = ld[i,choice[2]]
    text(xtext,ytext,labels=clrnames[i],col=colray,
         pos=2,offset=0.3,cex=cextext)
  }
  if(optpdf==1){
    dev.off() }
  
  lout = list("U",U,"V",V,"D",D,"explvar",explvar,"totvar",totvar)
  return(lout)
} 
#############################################


###### mcountsprop function ##################
# estimates proportions of multinomial samples
# including zeros or low counts.
# The number of multinomial trials is different 
# from row to row of the samples in X: D columns
# for categories; nsam rows number of multinomial
# samples.
# The procedure is based on 
# GMB estimation of proportions in
# GMB zero count replacement by Martin-Fernandez et al.2015 
# Statistical Modelling
# It consists of a point Bayesian estimation of 
# multinomial proportions (posterior mean value). The prior
# proportions are assess on using leave-one-out. The strength
# is computed for each row as the geometric mean of the
# observed counts.
# On the output, the estimated matrix of proportions. 
##############################################
#  JJ Egozcue, May 2017, amended May 2018
##############################################
mcountsprop<-function(X){
	D = ncol(X)
	n = nrow(X)
	# number of multinomial trials per row
	# and number of counts per category
	numtrials = rowSums(X)
	numcat = colSums(X)
	# compute m : leave-one-out proportions along sample
  alpha = outer(rep(1,n),numcat) - X
  salpha = rowSums(alpha)
  # prior proportions
  m = alpha / outer(salpha,rep(1,D))
  # strengths
  ss = exp(-rowMeans(log(m)))
  # counts+prior proportions
  Xpr = X + ss * m
  # estimated proportions 
  Xpr1 = Xpr / outer(rowSums(Xpr),rep(1,D) )
  return(Xpr1)
}
################################################

####### mvariatzeros ###########################
# estimates the variation matrix of a compositional
# data set (by rows) possibly containing zeros.
# Zero values are ignored when computing variances 
# of simple log-ratios.
# At the end of estimation the negative definition 
# of the matrix is checked and corrected if necessary
################################################
# JJ Egozcue, May 2017
################################################
mvariatzeros<-function(X,toleig=0.000001){
	# put negative numbers in place of 0's
	D = ncol(X)
	mvariatzeros=matrix(0,ncol=D,nrow=D)
  Xna = X -(X<=0)	
  # take log's producing NA's
  Xlog = log(Xna)
  #indices combining two columns of Xlog
  indcomb1=rep(1:D,D)
  indcomb2=as.vector(t(matrix(indcomb1,nrow=D,ncol=D)))
  # log differences = log ratios
  for(i in 1:D){
  	for(j in 1:i){
      xx = Xlog[,i]-Xlog[,j]
      mvariatzeros[i,j]=var(xx,na.rm=TRUE)
  	}
  }
  mvariatzeros = mvariatzeros + t(mvariatzeros)
  # check positive definiteness and correct to toleig
  Xsvd=svd(mvariatzeros,nu=0,nv=0)
  if(any(Xsvd[[1]]<=0)){
  	 print("variation is not positive definite")
  	 print("correction applied, non zero diagonal elements")
  	 print("give an idea of the size of error")
  	 Xsvd=svd(mvariatzeros)
  	 Lamb=((Xsvd$d>0)*Xsvd$d) + (Xsvd$d<=0)*toleig
  	 mvariatzeros = as.matrix(Xsvd$u) %*% diag(Lamb) %*% t(as.matrix(Xsvd$v))
  }
  return(mvariatzeros)
}
###############################################

#### function checks of functions #####
mcheckfunctions <- function(){
	# some compositional data in X
print("Data")
	X = rbind(  c(1,2,3,4,5),
							c(1,1,2,3,4),
							c(7,4,3,1,1),
							c(1,1,1,1,1) )
	X1=X[1,]
# check mcls
print("check mcls")
  print(X1)
	cX = mcls(X)
	print(cX)
	print(rowSums(cX))
	cX1=mcls(X1)
	print(cX1)

# check  mbuildcontrast
print("mbuildcontrast")
	W = cbind(  c(1,1,-1,-1,-1),
							c(1,-1,0,0,0),
							c(0,0,1,-1,-1),
							c(0,0,0, 1,-1)  )
	W1 =  c(1,1,-1,-1,-1)
	
	V=mbuildcontrast(W)
	print(W)
  print(V)
  V1=mbuildcontrast(W1)
  print(W1)
  print(V1)
  
#  check milr(x, V) and milrInv(z, V)
print("milr and milrInv")
  balX=milr(X,V)
  print(balX)
  XX = milrInv(balX,V)
  print((X/XX))
  print((mcls(X)/XX))

# check mclr
print("check mclr")
    Xclr= mclr(X)
  print(Xclr)
  XX = mcls(exp(Xclr))
  print(mcls(X)/XX)
  
# check mvariation
print("check mvariation")
  Xvar = mvariation(X)
  print(Xvar)  
  covbalX= cov(balX)
  XXvar=  (t(V) %*% Xvar %*% V)*(-0.5) 
  print(covbalX/XXvar)

# check mBPPOP 
print("check mBPPOP")
  colnames(X)=c("one","two","three","four","five")
  Xpop=c(1,1,2,2)
  
  mBPPOP(x=X,xpop=Xpop,choice=c(1,2),biscale=1,punch=1,
  			     colist=c("blue","springgreen4"),
             optpdf=0,
             colray="red",cextext=1,lwdray=1,pchpoint=1)
  # position of the center
  # near to (0,0,0,0,0) but a little shift
  print(colMeans(Xclr))
  
  # check mvariatzeros
  Xmvar=mvariatzeros(X)
  print("comparison variation matrices")
  print(Xmvar)
  print(Xvar)
  print("compare ilr covariance")
  ilrcov = ( t(V) %*% Xmvar %*% V)/(-2)
  print(ilrcov)
  print(var(balX))
}  
######################################################

###################################################
# function mmerge2sign transform the tree code $merge
# obtained in a hyerachical cluster into a sign
# code of a SBP. The result can be used to construct
# a contrast matrix V using mbuildConstrast
# It is based on the function gsi.merge2signary in
# package compositions.
# Merge a (D,2) matrix of type $merge
# V on the output a SBP sign code (D,D-1) matrix
# containing the SBP code.
##############################################
# JJ Egozcue, August 2017
##############################################
mmerge2sign<-function(Merge){
    V = matrix(0, ncol = nrow(Merge) + 1, nrow = nrow(Merge))
    for (i in 1:nrow(Merge)) {
        for (j in 1:2) {
            weight = (-1)^j
            k = Merge[i, j]
            if (k < 0) {
                V[i, abs(k)] = weight
            }
            if (k > 0) {
                take = as.logical(V[k, ])
                V[i, take] = rep(weight, sum(take))
            }
        }
    }
    revV = V[nrow(V):1,]
    return(t(revV))
}
################################################

###############################################
# function mPBclustvar(x)
# given a CoDa data set x[,1:D] a cluster analysis
# of the variables (columns) is carried out using
# Ward's method in stats{hclust}
# no zeros or negatives allowed
# The contrast matrix is computed and returned
#################################################
# JJ Egozcue, August 2017
#################################################
mPBclustvar<-function(x){
	# distances from variation
  hdist = sqrt(mvariation(x))
  hdist=as.dist(hdist)
  # cluster
  hmerge=hclust(d=hdist,method="ward.D2")
  # signs from $merge; contrast matrix
  Vsigns=mmerge2sign(hmerge$merge)
  V = mbuildcontrast(Vsigns)
  rownames(V)=colnames(x)
  return(V)
}
#################################################

##################################################
# function mTablediscrimErrors(ogroup,predgroup)
# Builds a table of centers of probabilities 
# in the observed groups after a discriminant
# analysis (probabilities as compositional data)
# it substitutes the traditional table of counts
# of predicted against observed.
################################################
# JJ Egozcue, August 2017
##################################################
mTablediscrimErrors<-function(ogroup,predgroup,namegroups=NULL){
	ngroup = ncol(predgroup)
	if(is.null(namegroups)){namegroups=1:ngroup}
	Cenp=matrix(0,ncol=ngroup,nrow=ngroup)
	colnames(Cenp)=paste("p-",namegroups,sep="")
	rownames(Cenp)=paste("o-",namegroups,sep="")
	for(i in 1:ngroup){
	 Cenp[i,]=mcls(  exp(colMeans(log(predgroup[ogroup==i,]))) )
	}
  return(Cenp)
}
#######################################################





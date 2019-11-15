################################################################################
#
#                                                                                     
#   Filename:	  multiset_data_generator.R    												  
#                                                                                     
#   Project :   BiomJ article "Multiset sparse redundancy analysis for high 
#               dimensional omics data"
#   Date    :   24-01-2018
#
#
################################################################################


multiset_data_generator <- function(N = 500, k = 4, m = 2,
                                    p0=c(1000,500,200,10),
                                    p1=c(10,8,8,2)){
  #GENERATE DATA####
  # input
  ## size of the data
  # N = 500    # number of individuals
  # k = 4      # number of datasets
  # m = 2      # number of latent variables (LV's) per dataset
  # p0         # number of noise variables per dataset
  # p1         # number of highly correlated variables per dataset
  
  # association parameters between the LV's of the k datasets
  # e.g. ksi(dataset4) = b1*ksi(dataset3) + b2*ksi(dataset2) + ...
  b=array(0,dim=c(m,k,k))
  b[1,1,2]=0.8
  b[1,2,3]=0.7
  b[1,3,4]=0.4
  b[2,1,2]=0.6
  b[2,2,3]=0.6
  b[2,3,4]=0.2

  
  # specify the regression coefficients of the relevant variables per dataset on the associated LVs
  # X = a1*ksi1 + a2*ksi2 + ....
  a=array(0,dim=c(max(p1),m,k))
  a[1:p1[1],1:2,1]=
    matrix(
      c(.5,0,
        .5,0,
        .5,0,
        .3,0,
        .3,0,
        .2,0.2,
        0.2,0.3,
        0.2,0.7,
        0.2,0.6,
        0.2,0.6),nrow=p1[1],ncol=m,byrow=TRUE)
  a[1:p1[2],1:2,2]=
    matrix(
      c(.5,0,
        .5,0,
        .5,0,
        .4,0,
        0.4,0,
        0.3,0.7,
        0.2,0.6,
        0.2,0.6),nrow=p1[2],ncol=m,byrow=TRUE)
  a[1:p1[3],1:2,3]=
    matrix(
      c(.5,0,
        .5,0,
        .5,0,
        .3,0,
        .3,0,
        .2,0.2,
        0.2,0.3,
        0.2,0.6),nrow=p1[3],ncol=m,byrow=TRUE)
  a[1:p1[4],1:2,4]=
    matrix(
      c(.5,.1,
        .1,0.6),nrow=p1[4],ncol=m,byrow=TRUE)
  
  
  # generate ksi's
  ksi=array(NA,dim=c(N,m,k))
  for (j in 1:m) {                                     # loop over the number of LV's
    ksi[,j,1] = rnorm(N,0,1)                          # generate values for the LV of the first set of variables
    for (el in 2:k) {                                 # loop over the other sets of variables
      meanx=0
      sumb2=0
      for (ell in (el-1):1) {                        # calculate per person the mean and sd of the LVs
        meanx=meanx+b[j,ell,el]*ksi[,j,ell]
        sumb2=sumb2+b[j,ell,el]
      }
      sdx=max(0.0001,sqrt(1-sumb2))
      ksi[,j,el] = rnorm(N,meanx,sdx)                # sample LV values from the normal distribution
    }
    ksi[,j,]=scale(ksi[,j,])
  }
  
  # generate manifest data
  X=array(NA,dim=c(N,max(p0+p1),k))
  for (el in 1:k) {
    for (j in 1:p0[el]) {                            # sample values for the irrelevant manifest variables
      X[,j,el]=rnorm(N,0,1)
    }
    meanx=ksi[,1:m,el]%*%t(a[1:p1[el],1:m,el])       # calculate per person means and sds of the relevant manifest variables
    suma2=apply(a[,,el],1,sum)
    for (j in (p0[el]+1):(p0[el]+p1[el])) {
      sdx=min(0.0001,sqrt(1-suma2[(j-p0[el])]))
      X[,j,el]=rnorm(N,meanx[,(j-p0[el])],sdx)      # sample values for the manifest variables
    }
  }
  
  X1=X[,,1]
  X2=X[,,2]
  X3=X[,,3]
  X4=X[,,4]
  if (length(which(is.na(apply(X1,2,sd,na.rm=TRUE))))>0) {X1=X1[,-which(is.na(apply(X1,2,sd,na.rm=TRUE)))]}
  if (length(which(is.na(apply(X2,2,sd,na.rm=TRUE))))>0) {X2=X2[,-which(is.na(apply(X2,2,sd,na.rm=TRUE)))]}
  if (length(which(is.na(apply(X3,2,sd,na.rm=TRUE))))>0) {X3=X3[,-which(is.na(apply(X3,2,sd,na.rm=TRUE)))]}
  if (length(which(is.na(apply(X4,2,sd,na.rm=TRUE))))>0) {X4=X4[,-which(is.na(apply(X4,2,sd,na.rm=TRUE)))]}
  
  c(dim(X1),dim(X2),dim(X3),dim(X4))
  
  X1=scale(X1)
  X2=scale(X2)
  X3=scale(X3)
  X4=scale(X4)
  
  result <- list(X1,
                 X2,
                 X3,
                 X4
  )
  
  names(result) <- c("X1",
                     "X2",
                     "X3",
                     "X4"
  )
  
  result
}

#end of function **********
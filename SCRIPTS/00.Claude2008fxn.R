# GEOMETRIC MORPHOMETRICS
# from Morphometrics with R (Claude 2008)

# Function for centering a configuration on the centroid (version 2)
trans1<-function(M){
  scale(M,scale=F)
}


# Function to calculate centroid coordinates of a landmark configuration
centcoord<-function(M){
  apply(M,2,mean)
}


# Function to calculate centroid size and scale a configuration to unit centroid size
centsiz<-function(M){
  p<-dim(M)[1]
  size<-sqrt(sum(apply(M,2,var))*(p-1))
  list("centroid_size"=size, "scaled"=M/size)
}

# Function for calculating the p interlandmark distances between two configurations
ild2<-function(M1,M2){
  sqrt(apply((M1-M2)^2,1,sum))
}

# Function for calculating average shape from an array
# This is the same function as used in the geomorph package
mshape<-function(A){
  apply(A,c(1,2),mean)
}


# Function for placing two configurations in partial Procrustes superimposition
# Requires centsiz, trans1, and ild2 functions (above)
pPsup<-function(M1,M2){
  k<-ncol(M1)
  Z1<-trans1(centsiz(M1)[[2]])
  Z2<-trans1(centsiz(M2)[[2]])
  sv<-svd(t(Z2)%*%Z1)
  U<-sv$v; V<-sv$u; Delt<-sv$d
  sig<-sign(det(t(Z2)%*%Z1))
  Delt[k]<-sig*abs(Delt[k]); V[,k]<-sig*V[,k] # to remove any reflection
  Gam<-U%*%t(V)
  beta<-sum(Delt)
  list(Mp1=Z1%*%Gam,Mp2=Z2,rotation=Gam,DP=sqrt(sum(ild2(Z1%*%Gam,Z2)^2)),rho=acos(beta))
}

# Function for placing two configurations in full Procrustes superimposition
# Requires centsiz and trans1 functions (above)
fPsup<-function(M1,M2){
  k<-ncol(M1)
  Z1<-trans1(centsiz(M1)[[2]])
  Z2<-trans1(centsiz(M2)[[2]])
  sv<-svd(t(Z2)%*%Z1)
  U<-sv$v; V<-sv$u; Delt<-sv$d
  sig<-sign(det(t(Z2)%*%Z1))
  Delt[k]<-sig*abs(Delt[k]); V[,k]<-sig*V[,k] # to remove any reflection
  Gam<-U%*%t(V)
  beta<-sum(Delt)
  list(Mp1=beta*Z1%*%Gam,Mp2=Z2,rotation=Gam,scale=beta,DF=sqrt(1-beta^2))
}


# Partial General Procrustes Analysis
# Requires trans1, centsiz, ild2, mshape, and pPsup functions (above)
pgpa<-function(A){
  p<-dim(A)[1]; k<-dim(A)[2]; n<-dim(A)[3]
  temp2<-temp1<-array(NA,dim=c(p,k,n)) # empty array for storing configurations
  Siz<-numeric(n) # zero vector to contain centroid size
  for (i in 1:n){ # translate and scale configurations to unit centroid size
    Acs<-centsiz(A[,,i])
    Siz[i]<-Acs[[1]]
    temp1[,,i]<-trans1(Acs[[2]])
  }
  Qm1<-dist(t(matrix(temp1,k*p,n))) # quantity that must be minimized during rotation
  Q<-sum(Qm1); iter<-0
  while (abs(Q)>0.00001){ # loop until quantity is minimized
    for (i in 1:n){
      M<-mshape(temp1[,,-i]) # define mean shape ignoring configuration that is going to be rotated
      temp2[,,i]<-pPsup(temp1[,,i],M)[[1]]
    }
    Qm2<-dist(t(matrix(temp2,k*p,n)))
    Q<-sum(Qm1)-sum(Qm2)
    Qm1<-Qm2
    iter=iter+1
    temp1<-temp2
  }
  list("rotated"=temp2,"it.number"=iter,"Q"=Q,"intereuclid.dist"=Qm2,"mshape"=centsiz(mshape(temp2))[[2]],"cent.size"=Siz)
}


# Function for rotating configurations along their major axes, checking for eventual reflection
aligne<-function(A){
  B<-A
  n<-dim(A)[3]; k<-dim(A)[2]
  for (i in 1:n){
    Ms<-scale(A[,,i],scale=F)
    sv<-eigen(var(Ms))
    M<-Ms%*%sv$vectors
    B[,,i]<-M
  }
  B
}


# Function for projecting from shape space into tangent space (sterographic projection)
# Configurations must have already been put in Procrustes superimposition
# Requires mshape function (above)
stp<-function(A){
  p<-dim(A)[1]; k<-dim(A)[2]; n<-dim(A)[3]
  Yn<-mshape(A)
  B<-array(NA,dim=c(p,k,n))
  for (i in 1:n){
    rho<-2*asin((sqrt(sum((A[,,i]-Yn)^2)))/2)
    B[,,i]<-A[,,i]/(cos(rho))
  }
  return(B)
}


# Function for projecting from shape space into tangent space (orthogonal projection)
# Configurations must have already been put in Procrustes superimposition
# Requires mshape and centsiz functions (above)
orp<-function(A){
  p<-dim(A)[1]; k<-dim(A)[2]; n<-dim(A)[3]
  Y1<-as.vector(centsiz(mshape(A))[[2]])
  oo<-as.matrix(rep(1,n))%*%Y1
  I<-diag(1,k*p)
  mat<-matrix(NA,n,k*p)
  for (i in 1:n){
    mat[i,]<-as.vector(A[,,i])
  }
  Xp<-mat%*%(I-(Y1%*%t(Y1)))
  Xp1<-Xp+oo
  array(t(Xp1),dim=c(p,k,n))
}


# Function for partial Procrustes superimposition with reference form aligned along principal axis
# Data are also projected into tangent space with in orthogonal projection
# Requires pgpa, trans1, centsiz, mshape, and pPsup functions
procalign<-function(A){
  pA<-pgpa(A); n<-dim(A)[3]; k<-dim(A)[2]
  A<-pA$rotated; msh<-pA$mshape
  A1<-A
  sv<-eigen(var(msh))
  V<-sv$vectors;
  rotmsh<-msh%*%V
  for (i in 1:n){
    A1[,,i]<-pPsup(A[,,i],rotmsh)$Mp1
  }
  list("rotated"=orp(A1), "meansh"=rotmsh)
}



# Function for calculating uniform terms following Bookstein 1996
# This function yields uniform scores that are essentially identical to those produced by PGAgen (IMP)
# Requires procalign, pgpa, trans1, centsiz, mshape, and pPsup functions
uniform2D<-function(A){
  n<-dim(A)[3]
  kp<-dim(A)[1]*dim(A)[2]
  temp<-procalign(A)
  msh<-temp$meansh
  proc<-temp$rotated
  X<-t(matrix(proc,kp,n))
  V<-X-rep(1,n)%*%t(as.vector(msh))
  alph<-sum(msh[,1]^2)
  gam<-sum(msh[,2]^2)
  U1<-c(sqrt(alph/gam)*msh[,2],sqrt(gam/alph)*msh[,1])
  U2<-c(-sqrt(gam/alph)*msh[,1],sqrt(alph/gam)*msh[,2])
  score<-V%*%cbind(U1,U2)
  list("scores"=score,"uniform"=cbind(U1,U2),"meanshape"=msh,"rotated"=proc)
}




# THIN-PLATE SPLINES
# Function for calculating the position of interpolated coordinates
tps2d<-function(M,matr,matt){
  p<-dim(matr)[1]; q<-dim(M)[1]; n1<-p+3
  P<-matrix(NA,p,p) # compute matrices P, Q, and L
  for (i in 1:p){
    for (j in 1:p){
      r2<-sum((matr[i,]-matr[j,])^2)
      P[i,j]<-r2*log(r2) # calculate U-functions
    }
  }
  P[which(is.na(P))]<-0
  Q<-cbind(1,matr)
  L<-rbind(cbind(P,Q),cbind(t(Q),matrix(0,3,3)))
  m2<-rbind(matt,matrix(0,3,2))
  coefx<-solve(L)%*%m2[,1]
  coefy<-solve(L)%*%m2[,2]
  fx<-function(matr,M,coef){
    Xn<-numeric(q)
    for (i in 1:q){
      Z<-apply((matr-matrix(M[i,],p,2,byrow=T))^2,1,sum)
      Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*(Z*log(Z)))
    }
    Xn
  }
  matg<-matrix(NA,q,2)
  matg[,1]<-fx(matr,M,coefx)
  matg[,2]<-fx(matr,M,coefy)
  matg
}


# Function for drawing a TPS
# Requires tps2d function (above)
tps<-function(matr,matt,n){ # n is the desired number of grid columns in the TPS
  xm<-min(matt[,1])
  ym<-min(matt[,2])
  xM<-max(matt[,1])
  yM<-max(matt[,2])
  rX<-xM-xm; rY<-yM-ym # define the range of the graph to estimate grid size
  a<-seq(xm-1/5*rX,xM+1/5*rX,length=n)
  b<-seq(ym-1/5*rX,yM+1/5*rX,by=(xM-xm)*7/(5*(n-1)))
  m<-round(0.5+(n-1)*(2/5*rX+yM-ym)/(2/5*rX+xM-xm))
  M<-as.matrix(expand.grid(a,b))
  ngrid<-tps2d(M,matr,matt) # calculate coordinates at gridline intersections
  plot(ngrid,cex=0.2,asp=1,axes=F,xlab="",ylab="")
  for (i in 1:m){
    lines(ngrid[(1:n)+(i-1)*n,])
  }
  for (i in 1:n){
    lines(ngrid[(1:m)*n-i+1,])
  }
}
rm(list=ls())
load("C:/Users/zhu2/Documents/getpathway/model20170215/disease.rda")
load("C:/Users/zhu2/Documents/getpathway/model20170215/phenet.rda")
phe <- phenet[[1]]
disease[,5] <- ifelse(disease[,5]%in%c(4,5),1,0); disease <- disease[,-1]
colnames(disease) <- strsplit('hypertension,diabetes,lewy_body,AD,obesity',',')[[1]]

qpca <- function(A,rank=0,ifscale=TRUE){
  if(ifscale){A <- scale(as.matrix(A))[,]}
  A.svd <- svd(A)
  if(rank==0){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
  }
  d <- d[d > 1e-8]
  r <- length(d)
  prop <- d^2; info <- sum(prop)/sum(A.svd$d^2);prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop,info=info)
  return(rlt)
}

qq_plot <- function(p_value){
  n = length(p_value);
  exp = -log10((c(1:n)-0.5)/n);
  rgen = -log10(sort(p_value));
  plot(exp,rgen,xlab="-log10(Expect)",ylab="-log10(Real)");
  abline(0,1,col="red")
}

p_ginv_sq <- function(X,p){
  X.eigen = eigen(X);
  X.rank = sum(X.eigen$values>1e-8);
  X.value = X.eigen$values[1:X.rank]^(-1*p);
  if (length(X.value)==1){
    D = as.matrix(X.value);
  }else{
    D = diag(X.value);
  }
  rlt = X.eigen$vectors[,1:X.rank] %*% D %*% t(X.eigen$vectors[,1:X.rank]);
  return(rlt);
}
mrank <- function(X){
  X.svd = svd(X);
  X.rank = sum(X.svd$d>1e-6);
  return(X.rank);
}
mrank_sq <- function(X){
  X.eigen = eigen(X);
  X.rank = sum(Re(X.eigen$values)>1e-6);
  return(X.rank);
}
CCA_chisq_test <- function(rho,n,p,q){
  tstat = -1*n*sum(log(1-rho^2));
  p_value = pchisq(tstat,(p*q),lower.tail=FALSE);
  return(p_value);          
}
cca <- function(A,B){
  A <- scale(A); B <- scale(B)
  n = nrow(A);
  p = mrank(A);
  q = mrank(B);
  if (p <= q){
    X = A;
    Y = B;
  }else{
    X = B;
    Y = A;
  }
  R = p_ginv_sq(cov(X),0.5) %*% cov(X,Y) %*% p_ginv_sq(cov(Y),1) %*% cov(Y,X) %*% p_ginv_sq(cov(X),0.5);
  k = mrank_sq(R);
  d = Re(eigen(R)$values);
  rho = d[1:k]^(0.5);
  rho[rho >= 0.9999]=0.9;
  chisq_p = CCA_chisq_test(rho,n,p,q);
  return(c("chisq_p"=chisq_p,"df"=p*q));
}

################################

d2d <- sapply(1:ncol(disease),function(i){
  sapply(1:ncol(disease),function(j){
    cca(disease[,i,drop=F],disease[,j,drop=F])[1]
  })
})
dimnames(d2d) <- list(colnames(disease),colnames(disease))

#################################
# Conditional CCA
#################################

ccca <- function(x,y,con=rep(1,length(x))){
  if(all(x==con)){
    con <- rep(1,length(x))
  }
  conp <- cca(cbind(x[con==1]),cbind(y[con==1]))[1]
  tp <- t.test(y,y[con==1])$p.value
  return(c(n=sum(con==1),ccap=conp,tp=tp))
}

rlt.ccca <- do.call(rbind,
  lapply(1:ncol(disease),function(i){
    d <- disease[,i]
    rlt <- lapply(1:ncol(phe),function(j){
      p <- phe[,j]
      t(sapply(1:ncol(disease),function(k){
        d2 <- disease[,k]
        c(d=i,p=j,d2=k,ccca(d,p,d2))
      }))
    })
    do.call(rbind,rlt)
  })
)

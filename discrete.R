rm(list=ls())
library(Matrix)
library(MASS)
library(survival)
library(glmnet)
library(Ckmeans.1d.dp)
library(igraph)
library(ncvreg)
setwd("e:/heye/subgroup analysis in cox model/revision_code/discrete/")

######data generating process######
genedata<-function(sam,n1,n2,n3,p1,p2,p3,c)
{
  set.seed(sam)
  label<-c(rep("1",n1),rep("2",n2),rep("3",n3))####three prognostic subgroup with the subgroup size are n1,n2,n3 respectively. 
  resam<-c(sample(rep(5:30),n,replace=TRUE)) ####cluster size m_i randomly sampled from 5 to 30
  obs<-rep(1:n,resam)
  labels<-rep(label,resam)
  df<-data.frame(obs,labels)
  N<-sum(resam)
  x<-matrix(rbinom(p*N,1,0.5),N,p)
  w<-rbinom(N,1,0.5)####the treatment
  u<-rexp(N,1)
  z<-matrix(0,N,q)####the confounder 
  z[,1]<-rbinom(N,1,0.5)
  z[,2]<-rnorm(N,0,1)
  index<-aggregate(resam,by=list(type=label),sum)[,2]
  treat_v<-rep(alpha3,times=index)####the treatment effect 
  eta_v<-rep(alpha2,times=index)######the frailty effect 
  inter_v<-rep(alpha4,times=c(p1,p2,p3))####three predictive subgroups with the subgroup size are p1,p2 and p3. 
  bet<-numeric() 
  bet[1:index[1]]<-x[1:index[1],]%*%alpha11
  bet[(index[1]+1):sum(index[1:2])]<-x[(index[1]+1):sum(index[1:2]),]%*%alpha12
  bet[(sum(index[1:2])+1):N]<-x[(sum(index[1:2])+1):N,]%*%alpha13
  
  t<-u/exp(eta_v+treat_v*w+w*(x%*%inter_v)+bet+z%*%gamma)
  censor<-runif(N,0,c)  #####the censoring time 
  status<-(t<=censor); tau<-c(pmin(t,censor))
  time<-t[status==1];n0<-length(time)
  risk<-(outer(tau,rep(1,n0))>=outer(rep(1,N),time))
  return(list(time=tau,status=status,Risk=risk,X=x,W=w,Z=z,Label=resam,N=N,df=df))
}

######Update the individual prognostic effect######

fun_heter<-function(step,grad,eta0,v0,g0,gamma0,alpha10,alpha20,alpha30,beta0,cind0)
{
  betanew<-matrix(0,n,p)
  betanew[,ind0]<-0
  etanew<-eta0-step*grad[,1]
  vnew<-v0-step-step*grad[,2]
  betanew[,cind0]<-beta0[,cind0]-step*grad[,-c(1,2)]
  ex<-newx%*%as.vector(t(betanew))+ind%*%etanew+ind_treat%*%vnew+wx%*%g0+z%*%gamma0
  temp<-t(exp(ex))%*%risk####1*n0
  dis<-matrix(0,n,K0)
  for (k in 1:K0)
  {
    dis[,k]<-apply((betanew-kronecker(t(alpha10[k,]),rep(1,n)))^2,1,sum)+(etanew-alpha20[k])^2+(vnew-alpha30[k])^2
  }
  avedis<-(apply(dis^s,1,mean))^(1/s)
  fun<--sum(ex[status==1,])+sum(log(temp))+lam1*sum(avedis)
  return(fun)
}

######update the center of the prognostic subgroup ########

fun_homo<-function(thetaest,alphaest,lam3)
{ 
  theta0<-thetaest;alpha0<-alphaest
  wei<-matrix(0,n,K0)
  for(j in 1:n)
  {
    dis<-apply((t(theta0[j,]%o%rep(1,K0))-alpha0)^2,1,sum) ###################may be a mistake 
    W<-(mean(dis^s))^(1/s-1)
    wei[j,]<-W*dis^(s-1)
  }
  
  ###update alpha20 and alpha30###
  
  alp2<-apply(theta0[,1]%o%rep(1,K0)*wei,2,sum)/apply(wei,2,sum)
  alp3<-apply(theta0[,2]%o%rep(1,K0)*wei,2,sum)/apply(wei,2,sum)
  
  ###updat alpah10 with adaptive lasso####
  p0<-dim(alpha0)[2]-2
  alp1<-matrix(0,K0,p0)
  for (k in 1:K0)
  {
    newy<-as.vector(t(theta0[,-c(1,2)]*(sqrt(wei[,k])%o%rep(1,p0))))
    newx<-kronecker(sqrt(wei[,k]),diag(p0))
    fit0<-lm(newy~newx-1)
    fit<-ncvfit(newx,newy,alpha=1,lambda=lam3,penalty="lasso",penalty.factor=1/abs(fit0$coef))
    alp1[k,]<-fit$beta
  }
  return(list(alp1=alp1,alp2=alp2,alp3=alp3))
}

######update the predictive effect#####
fun1_heter<-function(step2,grad2,g0,betaest,etaest,vest,gamma0,alpha40)
{
  gnew<-g0-step2*grad2
  ex<-newx%*%as.vector(t(betaest))+ind%*%etaest+ind_treat%*%vest+wx%*%gnew+z%*%gamma0
  temp<-t(exp(ex))%*%risk 
  dis<-(gnew%o%rep(1,K10)-rep(1,p)%o%alpha40)^2
  avedis<-(apply(dis^s,1,mean))^(1/s)
  fun<--sum(ex[status==1,])+sum(log(temp))+lam2*sum(avedis)
  return(fun)
}

########update the center of the predictive effcet######

fun1_homo<-function(gest,alp4est)
{
  g0<-gest
  alpha40<-alp4est
  wei<-matrix(0,p,K10)
  dis<-(g0%o%rep(1,K10)-rep(1,p)%o%alpha40)^2
  W<-(apply(dis^s,1,mean))^(1/s-1)
  wei<-(W%o%rep(1,K10))*dis^(s-1)
  gchat<-apply(g0%o%rep(1,K10)*wei,2,sum)/apply(wei,2,sum)
  return(gchat=gchat)
}



######the main function #######

n=30;s=-1;n00<-200
p=50;q=2
K0=3;K10=3;n1=15;n2=7;n3=8
p1=p3=15;p2=20
lam1=3###tuning parameter for prognostic subgroup 
lam2=1####tuning parameter for predictive subgroup 
lam3=0.001###tuning parameter for variable selection 
realp1<-realp2<-list()
#####the prognostic effect of the biomarker######
alpha11<-matrix(c(-2,-2.5,0,rep(0,p-3)),p,1)
alpha12<-matrix(c(0,-1,-1.5,rep(0,p-3)),p,1)
alpha13<-matrix(c(2,1,1.5,rep(0,p-3)),p,1)
####the fraility #######
alpha2<-c(-1,0,1)
####the treatment effect####
alpha3<-c(-0.5,0,0.5)
####the interaction effect between treatment and biomarker######
alpha4<-c(-0.75,0,0.75)
####the confounder effect#####
gamma<-c(0,-0.5)
c=exp(4)

RI<-matrix(0,n00,5)
for (sam in 1:n00)
{
  data<-genedata(sam,n1,n2,n3,p1,p2,p3,c)
  head(data)
  time<-data$time;status<-data$status;x<-data$X;w<-data$W;z<-data$Z;label<-data$Label;N<-data$N;df<-data$df
  y<-Surv(time,status);risk<-data$Risk
  ind<-ind_treat<-matrix(0,N,n)
  newx<-matrix(0,N,n*p)
  for(j in 1:n)
  {
    m0<-label[j]
    if(j==1)
    {
      ind[1:m0,1]<-1
      ind_treat[1:m0,1]<-w[1:m0]
      newx[1:m0,1:p]<-x[1:m0,1:p]
    }
    if(j>1)
    {
      m1<-sum(label[1:(j-1)])
      ind[(m1+1):(m1+m0),j]<-1
      ind_treat[(m1+1):(m1+m0),j]<-w[(m1+1):(m1+m0)]
      newx[(m1+1):(m1+m0),((j-1)*p+1):(j*p)]<-x[(m1+1):(m1+m0),1:p]
    }
    
  }
  
  wx<-(w%o%rep(1,p))*x
  
  newxzw<-cbind(ind,ind_treat,newx,wx,z)

  #####initial value#######
  fit=glmnet(newxzw,y,family="cox",alpha=0,lambda=1,penalty.factor=c(rep(1,(2+p)*n),rep(0,p),rep(0,q)))
  
  eta0<-fit$beta[c(1:n)]
  v0<-fit$beta[(n+1):(2*n)]
  beta0<-matrix(fit$beta[(2*n+1):(n*(p+2))],n,p,byrow=TRUE)
  g0<-fit$beta[(n*(p+2)+1):(n*(p+2)+p)]
  clu0<-Ckmeans.1d.dp(g0,k=K10)
  index0<-clu0$cluster
  alpha40<-clu0$center
  gamma0<-fit$beta[-c(1:(n*(p+2)+p))]
  clu<-kmeans(cbind(eta0,v0,beta0),K0,nstart=25)
  cen<-clu$centers[order(clu$centers[,1]),]
  clu1<-kmeans(cbind(eta0,v0,beta0),K0,nstart=25,centers=cen)
  alpha20<-clu1$center[,1]
  alpha30<-clu1$center[,2]
  alpha10<-clu1$center[,-c(1,2)]
  index1<-clu1$cluster
  RI[sam,1]<-compare(index1,c(rep(1,n1),rep(0,n2),rep(-1,n3)),method="rand")####the RI for initial prognostic subgroup####
  RI[sam,2]<-compare(index0,c(rep(1,p1),rep(2,p2),rep(3,p3)),method="rand")####the RI for initial predictive subgroup#####
  RI[sam,3]<-sum(status)/N

#####loop######
  step_eps<-0.01
  e1=1;q1=1;eps=10^(-6);emax=500
  while(e1>eps&q1<=emax)
  {
    ind0<-which(apply(alpha10,2,function(x)(all(x==0)))==TRUE)
    cind0<-setdiff(1:p,ind0);p0<-length(ind0);alpha10[,ind0]<-0
    if(length(cind0)==1)
    {
      q1=emax+1
      break
    }
    beta0[,ind0]<-0
    ex<-newx%*%as.vector(t(beta0))+ind%*%eta0+ind_treat%*%v0+wx%*%g0+z%*%gamma0
    temp<-t(exp(ex))%*%risk#####1*n0
    ln1<-matrix(0,n,p-p0+2)
    theta0<-cbind(eta0,v0,beta0[,cind0])
    alpha0<-cbind(alpha20,alpha30,alpha10[,cind0])
    ln2<-matrix(0,n,dim(alpha0)[2])
    for (j in 1:n)
    {
      indj<-which(df$obs==j)
      xi<-cbind(w[indj],x[indj,cind0])
      deltaj<-status[indj]
      riskj<-risk[indj,]
      exj<-exp(eta0[j]+w[indj]*v0[j]+x[indj,cind0]%*%beta0[j,cind0]+wx[indj,]%*%g0+z[indj,]%*%gamma0)
      temp20<-t(riskj)%*%((exj[,1]%o%rep(1,p-p0+1))*xi)/(temp[1,]%o%rep(1,p-p0+1))
      temp30<-t(riskj)%*%(exj[,1])/temp[1,]
      ln1[j,1]<-sum(deltaj)-sum(temp30)
      if (sum(deltaj)==1)
      {
        ln1[j,-c(1)]<-xi[deltaj==1]-apply(temp20,2,sum)
      }
      if (sum(deltaj)>1)
      {
        ln1[j,-c(1)]<-apply(xi[deltaj==1,],2,sum)-apply(temp20,2,sum)
      }
      dis<-apply((t(theta0[j,]%o%rep(1,K0))-alpha0)^2,1,sum)
      dis[which(dis==0)]<-eps
      W<-(mean(dis^s))^(1/s-1)
      ln2[j,]<-apply(W*2*diag(dis^(s-1))%*%(t(theta0[j,]%o%rep(1,K0))-alpha0)/K0,2,sum)
    }
    
    grad<--ln1+lam1*ln2
  step_q1<-tryCatch({opt<-optimize(f=fun_heter,lower=0, 
  upper =1,grad=grad,eta0=eta0,v0=v0,g0=g0,gamma0=gamma0,alpha10=alpha10,alpha20=alpha20,alpha30=alpha30,beta0=beta0,cind0=cind0)
  opt$minimum},warning=function(w){
  print("warning")
  step_eps
  },finally={opt$minimum})

    betaest<-matrix(0,n,p)
    betaest[,ind0]<-0
    etaest<-eta0-step_q1*grad[,c(1)]
    vest<-v0-step_q1*grad[,c(2)]
    betaest[,cind0]<-beta0[,cind0]-step_q1*grad[,-c(1,2)]
    thetaest<-cbind(etaest,vest,betaest[,cind0])
    
    ####update the center of the prognostic effect#####

    alp<-fun_homo(thetaest,alpha0,lam3)####lam3
    alp1est<-matrix(0,K0,p)
    alp1est[,cind0]<-alp[[1]]
    alp2est<-alp[[2]]
    alp2est[order(alp2est)[2]]<-0
    alp3est<-alp[[3]]
    
    ####update the predictive effect####
    ex<-newx%*%as.vector(t(betaest))+ind%*%etaest+ind_treat%*%vest+wx%*%g0+z%*%gamma0
    temp<-t(exp(ex))%*%risk####1*n0
    temp20<-t(risk)%*%((exp(ex[,1])%o%rep(1,p))*wx)/(temp[1,]%o%rep(1,p))
    ln12<-apply(wx[status==1,],2,sum)-apply(temp20,2,sum)### the first order of the g
    
    dis<-(g0%o%rep(1,K10)-rep(1,p)%o%alpha40)
    dis[which(dis==0)]<-eps
    W<-(apply(dis^(2*s),1,mean))^(1/s-1)
    
    ln22<- W*apply(dis^(2*(s-1))*2*dis/K10,1,sum)
    grad2<--ln12+lam2*ln22

     step_q2<-tryCatch({opt<-optimize(f=fun1_heter,lower=0,upper =1,grad2=grad2,g0=g0,betaest=betaest,
     etaest=etaest,vest=vest,gamma0=gamma0,alpha40=alpha40)
     opt$minimum},warning=function(w){
     print("warning")
      step_eps},finally={opt$minimum})
    #print(step_q2)
    gest<-g0-step_q2*grad2
######update the center of the predictive effect#####
    alp4est<-fun1_homo(gest,alpha40)
    exw<-(newx%*%as.vector(t(betaest))+ind%*%etaest+ind_treat%*%vest+wx%*%gest)[,1]
    ####update the confounder effect####
    gamest<-coxph(y~z+offset(exw))$coef
    e1<-max((alp1est-alpha10)^2,(alp2est-alpha20)^2,(alp3est-alpha30)^2,(alp4est-alpha40)^2,(gest-g0)^2,(gamest-gamma0)^2)
    alpha10<-alp1est ####the center estimators biomarker 
    alpha20<-alp2est####fraility estimators 
    alpha30<-alp3est####treatment 
    alpha40<-alp4est ####predictive estimators 
    beta0<-betaest
    eta0<-etaest
    v0<-vest
    g0<-gest
    gamma0<-gamest####confounder 
    q1<-q1+1
  }

  km<-kmeans(cbind(eta0,v0,beta0),K0,nstart=25,centers=cbind(alpha20,alpha30,alpha10))
  index0<-km$cluster
  RI[sam,4]<-compare(index0,c(rep(1,n1),rep(0,n2),rep(-1,n3)),method="rand")
  print(alpha10)
  clu0<-Ckmeans.1d.dp(g0,k=K10)$cluster
  RI[sam,5]<-compare(clu0,c(rep(1,p1),rep(2,p2),rep(3,p3)),method="rand")
  ####refine estimators#######
  if((all((alpha10[1,]==alpha10[2,])==TRUE)|all((alpha10[1,]==alpha10[3,])==TRUE)|all((alpha10[3,]==alpha10[2,])==TRUE))==TRUE)
  {
    realp1[[sam]]<-alp1est[order(alp1est[,1]),]
    realp2[[sam]]<-c(sort(alp2est),sort(alp3est),sort(alp4est),gamest)
  }
  if((all((alpha10[1,]==alpha10[2,])==TRUE)|all((alpha10[1,]==alpha10[3,])==TRUE)|all((alpha10[3,]==alpha10[2,])==TRUE))==FALSE)
  {
    cind<-list()
    indexnew<-rep(0,n)
    indexnew[index0==order(alpha10[,1])[1]]=1
    cind[[1]]<-which(alpha10[order(alpha10[,1])[1],]!=0)
    indexnew[index0==order(alpha10[,1])[2]]=2
    cind[[2]]<-which(alpha10[order(alpha10[,1])[2],]!=0)
    indexnew[index0==order(alpha10[,1])[3]]=3
    cind[[3]]<-which(alpha10[order(alpha10[,1])[3],]!=0)
    
    p10<-length(cind[[1]])
    p20<-length(cind[[2]])
    p30<-length(cind[[3]])
    ppp<-rep(1:K0,c(p10,p20,p30))
    if (p10==0|p20==0|p30==0)
    {
      realp1[[sam]]<-alp1est[order(alp1est[,1]),]
      realp2[[sam]]<-c(sort(alp2est),sort(alp3est),sort(alp4est),gamest)
    }
    
    if(p10!=0&p20!=0&p30!=0)
    {
      
      p00<-p10+p20+p30
      cnew<-matrix(0,N,K0)
      wnew<-matrix(0,N,K0)
      xnew<-matrix(0,N,p00)
      for (k in 1:K0)
      {
        ind1<-which(indexnew==k)
        cnew[df$obs%in%ind1,k]=1
        wnew[which(df$obs%in%ind1==TRUE),k]=w[which(df$obs%in%ind1==TRUE)]
        xnew[which(df$obs%in%ind1==TRUE),which(ppp==k)]=x[which(df$obs%in%ind1==TRUE),cind[[k]]]
        
      }
      
      wx1<-apply(wx[,clu0==1],1,sum)
      wx2<-apply(wx[,clu0==2],1,sum)
      wx3<-apply(wx[,clu0==3],1,sum)
      
      fit<-coxph(y~xnew+cnew[,c(1,3)]+wnew+wx1+wx2+wx3+z)
      
      for (k in 1:K0)
      {alp1est[k,c(cind[[k]])]<-fit$coef[which(ppp==k)]}
      
      alp2est[c(1,3)]<-fit$coef[(p00+1):(p00+2)]
      alp2est[c(2)]<-0
      alp3est<-cbind(fit$coef[(p00+K0):(p00+2*K0-1)])
      alp4est<-cbind(fit$coef[c((p00+2*K0):(p00+2*K0+K10-1))])
      alp5est<-cbind(fit$coef[-c(1:(p00+2*K0+K10-1))])
      realp1[[sam]]<-alp1est
      realp2[[sam]]<-c(alp2est,alp3est,alp4est,alp5est)
    }
  }
  print(sam)
}


save(realp1,file="(5,30)robustlam1=3,lam2=1,lam3=0.001alp1.RData")
save(realp2,file="(5,30)robustlam1=3,lam2=1,lam3=0.001alp2.RData")
write.table(RI,"(5,30)robustlam1=3,lam2=1,lam3=0.001RI.r")


re1<-re2<-re3<-matrix(0,n00,p)
re4<-matrix(0,n00,K0)
re5<-matrix(0,n00,K0)
re6<-matrix(0,n00,K10)
re7<-matrix(0,n00,q)
for (l in 1:n00)
{
  re1[l,]<-realp1[[l]][1,]
  re2[l,]<-realp1[[l]][2,]
  re3[l,]<-realp1[[l]][3,]
  re4[l,]<-realp2[[l]][1:K0]
  re5[l,]<-realp2[[l]][(K0+1):(2*K0)]
  re6[l,]<-realp2[[l]][(2*K0+1):(2*K0+K10)]
  re7[l,]<-realp2[[l]][(2*K0+K10+1):(2*K0+K10+q)]
  print(l)
}

re<-matrix(0,21,5)
re[20,1]<-1-mean(apply(re1[,c(1,2)],1,function(x)sum(x==0)))/2
re[20,2]<-1-mean(apply(re2[,c(2,3)],1,function(x)sum(x==0)))/2
re[20,3]<-1-mean(apply(re3[,c(1,2,3)],1,function(x)sum(x==0)))/3

delte<-which(apply(cbind(re1[,c(1,2)]==0,re2[,c(2,3)]==0,re3[,c(1,2,3)]==0),1,function(x)(any(x)))==TRUE)
retain<-setdiff(c(1:n00),delte)
bias1<-apply(re1[retain,c(1,2)],2,mean)-c(-2,-2.5)
sd1<-apply(re1[retain,c(1,2)],2,sd)
bias1/sd1

bias2<-apply(re2[retain,c(2,3)],2,mean)-c(-1,-1.5)
sd2<-apply(re2[retain,c(2,3)],2,sd)

bias2/sd2

bias3<-apply(re3[retain,c(1,2,3)],2,mean)-c(2,1,1.5)
sd3<-apply(re3[retain,c(1,2,3)],2,sd)

bias3/sd3

bias4<-apply(re4[retain,],2,mean)-c(-1,0,1)
sd4<-apply(re4[retain,],2,sd)
bias4/sd4


bias5<-apply(re5[retain,],2,mean)-c(-0.5,0,0.5)
sd5<-apply(re5[retain,],2,sd)
bias5/sd5

bias6<-apply(re6[retain,],2,mean)-c(-0.75,0,0.75)
sd6<-apply(re6[retain,],2,sd)
bias6/sd6

bias7<-apply(re7[retain,],2,mean)-c(0,-0.5)
sd7<-apply(re7[retain,],2,sd)
bias7/sd7

res<-matrix(0,5,3)
res[1,1]<-sqrt((sum(bias1^2)+sum(bias2^2)+sum(bias3^2)+sum(bias4^2)+sum(bias5^2))/12)
res[1,2]<-sqrt((sum(sd1^2)+sum(sd2^2)+sum(sd3^2)+sum(sd4^2)+sum(sd5^2))/12)
res[2,1]<-sqrt(mean(bias6^2))
res[2,2]<-sqrt(mean(sd6^2))
res[3,1]<-sqrt(mean(bias7^2))
res[3,2]<-sqrt(mean(sd7^2))
res[1:3,3]<-sqrt(res[1:3,1]^2+res[1:3,2]^2)

res[4,1:2]<-apply(RI[retain,-c(1:3)],2,mean)
res[5,1]<-1-(mean(apply(re1[,c(1,2)],1,function(x)sum(x==0))+apply(re2[,c(2,3)],1,function(x)sum(x==0))+apply(re3[,c(1,2,3)],1,function(x)sum(x==0))))/7
res[5,2]<-1-(mean(apply(re1[retain,-c(1,2)],1,function(x)sum(x!=0))+apply(re2[retain,-c(2,3)],1,function(x)sum(x!=0))+apply(re3[retain,-c(1,2,3)],1,function(x)sum(x!=0))))/143



round(res,digit=3)




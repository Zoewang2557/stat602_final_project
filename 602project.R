#########################cooking time

y1<-c(10.5,13,8.5,10,17.75,15.75,10.75,15.25,16.75,11.25,11,10,16.5,12.5,15,15.5,
      11.75,7.25,11,7.75,10,10.25,9.5,7.75,11.25,7,9.5,7.25,16.75,8.25,9.75,8)

A = rep(c(-1,1),16)
B = rep(c(-1,-1,1,1), 8)
C = rep(c(rep(-1,4),rep(1,4)), 4)
D = rep(c(rep(-1,8),rep(1,8)), 2) 
E = c(rep(-1,16), rep(1,16))
Block<-c('1','3','3','1','4','2','2','4','2','4','4','2','3','1','1','3',
         '2','4','4','2','3','1','1','3','1','3','3','1','4','2','2','4')
time<-data.frame(y1,A,B,C,D,E,Block)
fit_time<-lm(y1~A*B*C*D*E+Block)
summary(fit_time)
anova(fit_time)

est.effects<-2*fit_time$coefficients[-c(1,7,8,9,20,29,32)]
a.est<-abs(est.effects)
g<-length(est.effects)

########Lenth
s0<-1.5*median(a.est)
a1<-a.est[a.est<2.5*s0]
pse<-1.5*median(a1)
v<-g/3
alpha<-0.1
ga<-0.5*(1-(1-alpha)^(1/g))
a.est[a.est>qt(ga,v,lower.tail = F)*pse]

#########Dong
m1<-length(a1)
s12<-sum(a1^2)/m1
s1<-sqrt(s12)
a2<-a.est[a.est<2.5*s1]
m2<-length(a2)
s22<-sum(a2^2)/m2
s2<-sqrt(s22)
a.est>s2*qt(ga,m2,lower.tail = F)
a.est[a.est>s2*qt(ga,m2,lower.tail = F)]

par(mfrow=c(1,3),pty="s",cex=1.2)
rd.1<-lm(y1~E+Block)#residual plots of dong&lenth method
summary(rd.1)
res1<-rd.1$residuals
fit1<-rd.1$fitted.values
plot(res1 ~ fit1,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y1 ~ fit1,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(res1)
qqline(res1)

f1<-lm(y1~A+B+C+D+E+Block)
summary(f1)
resi <- f1$residuals#residual plots of transformation
fitt<- f1$fitted.values
plot(resi ~ fitt,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y1~ fitt,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resi)
qqline(resi)
########transformation

#notes method
gm <- exp(mean(log(y1))) ## geometric mean
ssr <- NULL
lambda.seq <- seq(from=-2,to=3,length.out=20) 
for(lambda in lambda.seq){
  if(lambda == 0){
    y <- gm*log(y1)
  } else {
    y <- (y1^lambda-1)/(lambda * gm^{lambda-1})
  }
  fit <- lm(y ~ A+B+C+D+E+Block)
  ssr <- c(ssr,sum(fit$resid^2))
} 
plot(lambda.seq,ssr,type="b",xlab=expression(lambda),
       ylab=expression(S[lambda]),main = "Select lambda")
y<-y1^(-0.5)
y<-y1^ll
model<-lm(y~A+B+C+D+E+Block)
summary(model)
anova(model)

resid <- model$residuals#residual plots of transformation
fitted <- model$fitted.values
plot(resid ~ fitted,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y ~ fitted,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resid)
qqline(resid)
#seems A,B,C,E


####selection
library(MASS)
time.t1<-lm(y1~A*B*C*D*E+Block-A:B:C:D:E)
time.bc1<-boxcox(time.t1,lambda = seq(-5,5-0.05))
time.bc1$x[which.max(time.bc1$y)]

time.t2<-lm(y1~(A+B+C+D+E)^3+Block)
time.bc2<-boxcox(time.t2)
time.bc2$x[which.max(time.bc2$y)]

time.t3<-lm(y1~(A+B+C+D+E)^2+Block)
resi<-time.t3$residuals
fitt<-time.t3$fitted.values
plot(resi ~ fitt,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y1 ~ fitt,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resi)
qqline(resi)

time.bc3<-boxcox(time.t3)
l<-time.bc3$x[which.max(time.bc3$y)]
y.t<-y1^l
model.t<-lm(y.t~(A+B+C+D+E)^2+Block)
summary(model.t)
anova(model.t)
resid.t<- model.t$residuals#residual plots of transformation
fitted.t <- model.t$fitted.values
fit.t<-fitted.t^(1/l)
plot(resid.t ~ fitted.t,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y1 ~ fit.t,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resid.t)
qqline(resid.t)



#tranformation studentized mm
coef.t<-model.t$coefficients[-c(1,7,8,9)]
N<-32
g.t<-15
v.t<-13
M<-3.1
se<-0.0014770
coef.t[abs(coef.t)>M*se]
#A,C,E,AE



######time dependence?
res.b1<-resid[c(1, 4, 14, 15, 22, 23, 25, 28)]
res.b2<-resid[c(6, 7,  9, 12, 17, 20, 30, 31)]
res.b3<-resid[c(2, 3, 13, 16, 21, 24, 26, 27)]
res.b4<-resid[c(5, 8, 10, 11, 18, 19, 29, 32)]
run1<-c(2,7,5,1,4,6,3,8)
run2<-c(2,6,3,4,1,8,7,5)
run3<-c(6,8,4,3,1,2,5,7)
run4<-c(8,4,1,3,2,7,5,6)
run<-c(run1,run2+8,run3+16,run4+24)
plot(run,resid.t,xlab = "run orders",ylab = "residuals",main = "run orders vs. residuals")
abline(v=8.5,col="darkorange",lty=3)
abline(v=16.5,col="darkorange",lty=3)
abline(v=24.5,col="darkorange",lty=3)
plot(run1,res.b1)
plot(run2,res.b2)
plot(run3,res.b3)
plot(run4,res.b4)










###########################silky
y2x<-c(6,5,7,7,6,7,8,6,4,6,6,6,6,7,8,7,6,7,5,6,8,6,7,7,7,7,7,7,5,5,7,8)
y2x<-c(6,5,7,5,6,7,8,6,4,6,6,6,8,7,8,7,7,6,6,6,8,6,7,7,7,7,7,8,6,5,7,8)
y2w<-c(7,5,8,6,6,6,8,7,6,5,7,6,8,7,9,7,8,5,7,5,9,6,9,7,7,7,7,7,6,6,8,8)
y2<-c(y2x,y2w)
y3x<-c(4,5,5,7,3,7,4,6,6,7,5,8,4,7,8,8,3,6,6,7,3,6,4,7,4,6,6,7,4,4,5,7)
y3w<-c(4,5,5,6,4,7,4,7,6,8,5,7,5,7,6,8,3,6,4,7,3,5,4,6,3,6,5,6,3,7,4,7)
y3<-c(y3x,y3w)

kg<-data.frame(y2,a,b,c,d,e,block)
eg<-data.frame(y3,a,b,c,d,e,block)
t<-data.frame(y4,a,b,c,d,e,block)

##taster1
fit_kg1<-lm(y2x~A*B*C*D*E+Block)
summary(fit_kg1)
anova(fit_kg1)
est.t1<-2*fit_kg1$coefficients[-c(1,7,8,9,20,29,32)]
a.t1<-abs(est.t1)
g.t1<-length(est.t1)
#Lenth
s0.t1<-1.5*median(a.t1)
a1.t1<-a.t1[a.t1<2.5*s0.t1]
pse.t1<-1.5*median(a1.t1)
v.t1<-g.t1/3
alpha<-0.1
ga.t1<-0.5*(1-(1-alpha)^(1/g.t1))
a.t1[a.t1>qt(ga.t1,v.t1,lower.tail = F)*pse.t1]
#Dong
m1.t1<-length(a1.t1)
s12.t1<-sum(a1.t1^2)/m1.t1
s1.t1<-sqrt(s12.t1)
a2.t1<-a.t1[a.t1<2.5*s1.t1]
m2.t1<-length(a2.t1)
s22.t1<-sum(a2.t1^2)/m2.t1
s2.t1<-sqrt(s22.t1)
a.t1>s2.t1*qt(ga.t1,m2.t1,lower.tail = F)
a.t1[a.t1>s2.t1*qt(ga.t1,m2.t1,lower.tail = F)]
#trans
gm.t1 <- exp(mean(log(y2x))) ## geometric mean
ssr <- NULL
lambda.seq <- seq(from=-2,to=6,length.out=30) 
for(lambda in lambda.seq){
  if(lambda == 0){
    y <- gm.t1*log(y2x)
  } else {
    y <- (y2x^lambda-1)/(lambda * gm.t1^{lambda-1})
  }
  fit <- lm(y ~ A+B+C+D+E+Block)
  ssr <- c(ssr,sum(fit$resid^2))
} 
plot(lambda.seq,ssr,type="b",xlab=expression(lambda),
     ylab=expression(S[lambda]))
y.t1<-y2x^2
model.t1<-lm(y.t1~A+B+C+D+E+Block)
summary(model.t1)
anova(model.t1)
coef.t1<-model.t1$coefficients[-c(1,7,8,9)]
g.t1<-5
v.t1<-23
M.t1<-2.46
se.t1<-2.05577
coef.t1[abs(coef.t1)>M.t1*se.t1]


##taster2
fit_kg2<-lm(y2w~A*B*C*D*E+Block)
summary(fit_kg2)
anova(fit_kg2)
est.t2<-2*fit_kg2$coefficients[-c(1,7,8,9,20,29,32)]
a.t2<-abs(est.t2)
g.t2<-length(est.t2)
#Lenth
s0.t2<-1.5*median(a.t2)
a1.t2<-a.t2[a.t2<2.5*s0.t2]
pse.t2<-1.5*median(a1.t2)
v.t2<-g.t2/3
alpha<-0.1
ga.t2<-0.5*(1-(1-alpha)^(1/g.t2))
a.t2[a.t2>qt(ga.t2,v.t2,lower.tail = F)*pse.t2]

rd.t2<-lm(y2w~A+B+C+Block)#residual plots of dong&lenth method
summary(rd.t2)
res.t2<-rd.t2$residuals
fit.t2<-rd.t2$fitted.values
plot(res.t2 ~ fit.t2,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y2w ~ fit.t2,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(res.t2)
qqline(res.t2)


#Dong
m1.t2<-length(a1.t2)
s12.t2<-sum(a1.t2^2)/m1.t2
s1.t2<-sqrt(s12.t2)
a2.t2<-a.t2[a.t2<2.5*s1.t2]
m2.t2<-length(a2.t2)
s22.t2<-sum(a2.t2^2)/m2.t2
s2.t2<-sqrt(s22.t2)
a.t2>s2.t2*qt(ga.t2,m2.t2,lower.tail = F)
a.t2[a.t2>s2.t2*qt(ga.t2,m2.t2,lower.tail = F)]

rd.t2<-lm(y2w~A+C+Block)#residual plots of dong&lenth method
summary(rd.t2)
res.t2<-rd.t2$residuals
fit.t2<-rd.t2$fitted.values
plot(res.t2 ~ fit.t2,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y2w ~ fit.t2,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(res.t2)
qqline(res.t2)


#trans
gm.t2 <- exp(mean(log(y2w))) ## geometric mean
ssr <- NULL
lambda.seq <- seq(from=-2,to=3,length.out=20) 
for(lambda in lambda.seq){
  if(lambda == 0){
    y <- gm.t2*log(y2w)
  } else {
    y <- (y2w^lambda-1)/(lambda * gm.t2^{lambda-1})
  }
  fit <- lm(y ~ A+B+C+D+E+Block)
  ssr <- c(ssr,sum(fit$resid^2))
} 
plot(lambda.seq,ssr,type="b",xlab=expression(lambda),
     ylab=expression(S[lambda]))
y.t2<-y2w^(1/3)
model.t2<-lm(y.t2~A+B+C+D+E+Block)
summary(model.t2)
anova(model.t2)

resid.t2 <- model.t2$residuals#residual plots of transformation
fitted.t2 <- model.t2$fitted.values
plot(resid.t2 ~ fitted.t2,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y2w ~ fitted.t2,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resid.t2)
qqline(resid.t2)
coef.t2<-model.t2$coefficients[-c(1,7,8,9)]
g.t2<-5
v.t2<-23
M.t2<-2.46
se.t2<-0.0122954
coef.t2[abs(coef.t2)>M.t2*se.t2]




###as replication
a<-rep(A,2)
b<-rep(B,2)
c<-rep(C,2)
d<-rep(D,2)
e<-rep(E,2)
block<-rep(Block,2)

dd<-data.frame(y2,a,b,c,d,e,block,taster)
fit_kg<-lm(y2~a*b*c*d*e+block,data = dd)
summary(fit_kg)
anova(fit_kg)
est.effects2<-2*fit_kg$coefficients[-c(1,7,8,9,20,29,32)]
res2 <- fit_kg$residuals#residual plots of transformation
fit2 <- fit_kg$fitted.values
plot(res2~fit2,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y2 ~ fit2,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(res2)
qqline(res2)
#smm
v2<-32
g2<-length(est.effects2)
se2.est<-0.09504*2
M30<-3.1
M35<-3.07
est.effects2[abs(est.effects2)>M30*se2.est]
#no need to do transformation



###################################egg

##taster1 A,B
fit_eg1<-lm(y3x~A*B*C*D*E+Block)
summary(fit_eg1)
anova(fit_eg1)
est.e1<-2*fit_eg1$coefficients[-c(1,7,8,9,20,29,32)]
a.e1<-abs(est.e1)
g.e1<-length(est.e1)
#Lenth
s0.e1<-1.5*median(a.e1)
a1.e1<-a.t1[a.e1<2.5*s0.e1]
pse.e1<-1.5*median(a1.e1)
v.e1<-g.e1/3
alpha<-0.1
ga.e1<-0.5*(1-(1-alpha)^(1/g.e1))
a.e1[a.e1>qt(ga.e1,v.e1,lower.tail = F)*pse.e1]

rd.e1<-lm(y3x~A+Block)#residual plots of dong&lenth method
summary(rd.e1)
res.e1<-rd.e1$residuals
fit.e1<-rd.e1$fitted.values
plot(res.e1 ~ fit.e1,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3x ~ fit.e1,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(res.e1)
qqline(res.e1)


#Dong
m1.e1<-length(a1.e1)
s12.e1<-sum(a1.e1^2)/m1.e1
s1.e1<-sqrt(s12.e1)
a2.e1<-a.t1[a.e1<2.5*s1.e1]
m2.e1<-length(a2.e1)
s22.e1<-sum(a2.e1^2)/m2.e1
s2.e1<-sqrt(s22.e1)
a.e1>s2.e1*qt(ga.e1,m2.e1,lower.tail = F)
a.e1[a.e1>s2.e1*qt(ga.e1,m2.e1,lower.tail = F)]

rd.e1<-lm(y3x~A+B+Block)#residual plots of dong&lenth method
summary(rd.e1)
res.e1<-rd.e1$residuals
fit.e1<-rd.e1$fitted.values
plot(res.e1 ~ fit.e1,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3x ~ fit.e1,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(res.e1)
qqline(res.e1)


#trans
gm.e1 <- exp(mean(log(y3x))) ## geometric mean
ssr <- NULL
lambda.seq <- seq(from=-2,to=6,length.out=30) 
for(lambda in lambda.seq){
  if(lambda == 0){
    y <- gm.e1*log(y3x)
  } else {
    y <- (y3x^lambda-1)/(lambda * gm.e1^{lambda-1})
  }
  fit <- lm(y ~ A+B+C+D+E+Block)
  ssr <- c(ssr,sum(fit$resid^2))
} 
plot(lambda.seq,ssr,type="b",xlab=expression(lambda),
     ylab=expression(S[lambda]))
y.e1<-y3x^(0.8)
model.e1<-lm(y.e1~A+B+C+D+E+Block)
summary(model.e1)
anova(model.e1)

resid.e1 <- model.e1$residuals#residual plots of transformation
fitted.e1 <- model.e1$fitted.values
fit.e1<-fitted.e1^(1/0.8)
plot(resid.e1 ~ fitted.e1,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3x ~ fit.e1,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resid.e1)
qqline(resid.e1)



ta1.t1<-lm(y3x~A*B*C*D*E+Block-A:B:C:D:E)
ta1.bc1<-boxcox(ta1.t1)
ta1.bc1$x[which.max(ta1.bc1$y)]#no way

ta1.t2<-lm(y3x~(A+B+C+D+E)^3+Block)
summary(ta1.t2)
resi<-ta1.t2$residuals
fitt<-ta1.t2$fitted.values
plot(resi ~ fitt,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3x ~ fitt,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resi)
qqline(resi)
ta1.bc2<-boxcox(ta1.t2)
(l2<-ta1.bc2$x[which.max(ta1.bc2$y)])
y.ta1<-y3x^l2
model.ta1<-lm(y.ta1~(A+B+C+D+E)^3+Block)
summary(model.ta1)
anova(model.ta1)
resid.ta1 <- model.ta1$residuals#residual plots of transformation
fitted.ta1 <- model.ta1$fitted.values
fit.ta1<-fitted.ta1^(1/l2)
plot(resid.ta1 ~ fitted.ta1,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3x ~ fit.ta1,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resid.ta1)
qqline(resid.ta1)

coef.e1<-ta1.t2$coefficients[-c(1,7,8,9,20,29)]
g.e1<-23
v.e1<-5
M.e1<-4.1
se.e1<-0.20010
coef.e1[abs(coef.e1)>M.e1*se.e1]

ta1.t3<-lm(y3x~(A+B+C+D+E)^2+Block)
resi<-ta1.t3$residuals
fitt<-ta1.t3$fitted.values
plot(resi ~ fitt,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3x ~ fitt,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resi)
qqline(resi)
ta1.bc3<-boxcox(ta1.t3)
(l2<-ta1.bc3$x[which.max(ta1.bc3$y)])
y.ta1<-y3x^l2
model.ta1<-lm(y.ta1~(A+B+C+D+E)^2+Block)
summary(model.ta1)
anova(model.ta1)
resid.ta1 <- model.ta1$residuals#residual plots of transformation
fitted.ta1 <- model.ta1$fitted.values
fit.ta1<-fitted.ta1^(1/l2)
plot(resid.ta1 ~ fitted.ta1,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3x ~ fit.ta1,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resid.ta1)
qqline(resid.ta1)
coef.e2<-ta1.t3$coefficients[-c(1,7,8,9)]
g.e2<-15
v.e2<-13
M.e2<-3.1
se.e2<-0.17442
coef.e2[abs(coef.e2)>M.e2*se.e2]

ta1.t4<-lm(y3x~A+B+C+D+E+Block)
ta1.bc4<-boxcox(ta1.t4)
(l2<-ta1.bc4$x[which.max(ta1.bc4$y)])
y.ta1<-y3x^l2
model.ta1<-lm(y.ta1~A+B+C+D+E+Block)
summary(model.ta1)
anova(model.ta1)
resid.ta1 <- model.ta1$residuals#residual plots of transformation
fitted.ta1 <- model.ta1$fitted.values
fit.ta1<-fitted.ta1^(1/l2)
plot(resid.ta1 ~ fitted.ta1,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3x ~ fit.ta1,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resid.ta1)
qqline(resid.ta1)






##taster2 A,D,E
fit_eg2<-lm(y3w~A*B*C*D*E+Block)
summary(fit_eg2)
anova(fit_eg2)
est.e2<-2*fit_eg2$coefficients[-c(1,7,8,9,20,29,32)]
a.e2<-abs(est.e2)
g.e2<-length(est.e2)
#Lenth
s0.e2<-1.5*median(a.e2)
a1.e2<-a.e2[a.e2<2.5*s0.e2]
pse.e2<-1.5*median(a1.e2)
v.e2<-g.e2/3
alpha<-0.1
ga.e2<-0.5*(1-(1-alpha)^(1/g.e2))
a.e2[a.e2>qt(ga.e2,v.e2,lower.tail = F)*pse.e2]

rd.e2<-lm(y3w~A+Block)#residual plots of dong&lenth method
summary(rd.e2)
res.e2<-rd.e2$residuals
fit.e2<-rd.e2$fitted.values
plot(res.e2 ~ fit.e2,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3w ~ fit.e2,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(res.e2)
qqline(res.e2)


#Dong
m1.e2<-length(a1.e2)
s12.e2<-sum(a1.e2^2)/m1.e2
s1.e2<-sqrt(s12.e2)
a2.e2<-a.e2[a.e2<2.5*s1.e2]
m2.e2<-length(a2.e2)
s22.e2<-sum(a2.e2^2)/m2.e2
s2.e2<-sqrt(s22.e2)
a.e2>s2.e2*qt(ga.e2,m2.e2,lower.tail = F)
a.e2[a.e2>s2.e2*qt(ga.e2,m2.e2,lower.tail = F)]

rd.e2<-lm(y3w~A+E+Block)#residual plots of dong&lenth method
summary(rd.e2)
res.e2<-rd.e2$residuals
fit.e2<-rd.e2$fitted.values
plot(res.e2 ~ fit.e2,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3w ~ fit.e2,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(res.e2)
qqline(res.e2)


#trans
gm.e2 <- exp(mean(log(y3w))) ## geometric mean
ssr <- NULL
lambda.seq <- seq(from=-2,to=3,length.out=20) 
for(lambda in lambda.seq){
  if(lambda == 0){
    y <- gm.e2*log(y3w)
  } else {
    y <- (y3w^lambda-1)/(lambda * gm.e2^{lambda-1})
  }
  fit <- lm(y ~ A+B+C+D+E+Block)
  ssr <- c(ssr,sum(fit$resid^2))
} 
plot(lambda.seq,ssr,type="b",xlab=expression(lambda),
     ylab=expression(S[lambda]))
(lambda<-1)
y.e2<-y3w
model.e2<-lm(y.e2~A+B+C+D+E+Block)
summary(model.e2)
anova(model.e2)

resid.e2 <- model.e2$residuals#residual plots of transformation
fitted.e2 <- model.e2$fitted.values
plot(resid.e2 ~ fitted.e2,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3w ~ fitted.e2,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resid.e2)
qqline(resid.e2)



ta2.t1<-lm(y3w~A*B*C*D*E+Block-A:B:C:D:E)
ta2.bc1<-boxcox(ta2.t1)
ta2.bc1$x[which.max(ta2.bc1$y)]#no way

ta2.t2<-lm(y3w~(A+B+C+D+E)^3+Block)
resi<-ta2.t2$residuals
fitt<-ta2.t2$fitted.values
plot(resi ~ fitt,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3x ~ fitt,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resi)
qqline(resi)
ta2.bc2<-boxcox(ta2.t2)
(l2<-ta2.bc2$x[which.max(ta2.bc2$y)])
y.ta2<-y3w^l2
model.ta2<-lm(y.ta2~(A+B+C+D+E)^3+Block)
summary(model.ta2)
anova(model.ta2)
resid.ta2 <- model.ta2$residuals#residual plots of transformation
fitted.ta2 <- model.ta2$fitted.values
fit.ta2<-fitted.ta2^(1/l2)
plot(resid.ta2 ~ fitted.ta2,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3w ~ fit.ta2,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resid.ta2)
qqline(resid.ta2)



ta2.t3<-lm(y3w~(A+B+C+D+E)^2+Block)
resi<-ta2.t3$residuals
fitt<-ta2.t3$fitted.values
plot(resi ~ fitt,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3x ~ fitt,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resi)
qqline(resi)
ta2.bc3<-boxcox(ta2.t3)
(l2<-ta2.bc3$x[which.max(ta2.bc3$y)])
y.ta2<-y3w^l2
model.ta2<-lm(y.ta2~(A+B+C+D+E)^2+Block)
summary(model.ta2)
anova(model.ta2)
resid.ta2 <- model.ta2$residuals#residual plots of transformation
fitted.ta2 <- model.ta2$fitted.values
fit.ta2<-fitted.ta2^(1/l2)
plot(resid.ta2 ~ fitted.ta2,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3w ~ fit.ta2,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resid.ta2)
qqline(resid.ta2)#this

coef.e2<-model.ta2$coefficients[-c(1,7,8,9)]
g.e2<-15
v.e2<-13
M.e2<-3.1
se.e2<-0.0092236
coef.e2[abs(coef.e2)>M.e2*se.e2]


ta2.t4<-lm(y3w~A+B+C+D+E+Block)
resi<-ta2.t4$residuals
fitt<-ta2.t4$fitted.values
plot(resi ~ fitt,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3x ~ fitt,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resi)
qqline(resi)
ta2.bc4<-boxcox(ta2.t4)
(l2<-ta2.bc4$x[which.max(ta2.bc4$y)])
y.ta2<-y3w^l2
model.ta2<-lm(y.ta2~A+B+C+D+E+Block)
summary(model.ta2)
anova(model.ta2)
resid.ta2 <- model.ta2$residuals#residual plots of transformation
fitted.ta2 <- model.ta2$fitted.values
fit.ta2<-fitted.ta2^(1/l2)
plot(resid.ta2 ~ fitted.ta2,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3w ~ fit.ta2,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resid.ta2)
qqline(resid.ta2)
coef.t<-ta2.t4$coefficients[-c(1,7,8,9)]
N<-32
g.t<-5
v.t<-23
M<-2.46
se<-0.12208
coef.t[abs(coef.t)>M*se]


###as replication
a<-rep(A,2)
b<-rep(B,2)
c<-rep(C,2)
d<-rep(D,2)
e<-rep(E,2)
block<-rep(Block,2)
taster<-c(rep('1',32),rep('2',32))
d<-data.frame(y3,a,b,c,d,e,block,taster)

fit_eg<-lm(y3~a*b*c*d*e+block+taster,data=d)
summary(fit_eg)
anova(fit_eg)
est.effects3<-2*fit_eg$coefficients[-c(1,7,8,9,10,21,30,33)]
res3 <- fit_eg$residuals#residual plots of transformation
fit3 <- fit_eg$fitted.values
plot(res3~fit3,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3 ~ fit3,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(res3)
qqline(res3)
#smm
v3<-31
g3<-length(est.effects3)
se3.est<-0.08531*2
M30<-3.1
est.effects3[abs(est.effects3)>M30*se3.est]


qqnorm(y3)
qqnorm(est.effects3)
qqline(est.effects3)
#transformation



rep.bc1<-boxcox(fit_eg,data=d)
(l3<-rep.bc1$x[which.max(rep.bc1$y)])
y.e<-y3^l3
da<-data.frame(y.e,a,b,c,d,e,block,taster)
model.e<-lm(y.e~a*b*c*d*e+block+taster,data = da)
summary(model.e)
anova(model.e)
resid.e <- model.e$residuals#residual plots of transformation
fitted.e <- model.e$fitted.values
fit.e<-fitted.e^(1/l3)
plot(resid.e ~ fitted.e,xlab="Fitted values",ylab="Residuals",main="Residuals vs fitted values")
plot(y3 ~ fit.e,xlab="Fitted values",ylab="Observed values",main="Observed values vs fitted values")
qqnorm(resid.e)
qqline(resid.e)


run<-c(run1,run2+8,run3+16,run4+24)
plot(run,res3[1:32],xlab = "run orders",ylab = "residuals",main = "taster1")
plot(run,res3[33:64],xlab = "run orders",ylab = "residuals",main = "taster2")
abline(v=8.5,col="darkorange",lty=3)
abline(v=16.5,col="darkorange",lty=3)
abline(v=24.5,col="darkorange",lty=3)

#Romain HATON et Leonardo MARTINS BIANCO
set.seed(999)
options(digits = 22)
library(bspec)
library(aod)
library(car)
#Partie 2
#Question 1)
n = 50
a = 5
b = 1
sigma2 = 2
XX = rnorm(n)
epsi = rnorm(50,mean = 0, sd = sigma2)
Z = a + b*XX + epsi 
Y = ((Z * 0.3) + 1)^(1/0.3)
modelZ=lm(Z~XX)
summary(modelZ)
modelZ.stdres = rstudent(modelZ) 
plot(modelZ$fitted.values, modelZ.stdres, ylab="Studentised Residuals",  xlab="Fitted") 
abline(0, 0) 
qqnorm(modelZ.stdres, main = "Normal Q-Q Plot", col = "darkgrey")
qqline(modelZ.stdres, col = "dodgerblue", lwd = 2)
modelY = lm(Y~XX)
summary(modelY)
modelY.stdres = rstudent(modelY) 
plot(modelY$fitted.values, modelY.stdres, ylab="Studentised Residuals",  xlab="Fitted") 
abline(0, 0) 
qqnorm(modelY.stdres, main = "Normal Q-Q Plot", col = "darkgrey")
qqline(modelY.stdres, col = "dodgerblue", lwd = 2)
#Question 2)
X = as.matrix(cbind(XX,a))

Q = diag(1,n) - X%*%solve(t(X)%*%X)%*%t(X) 
Lmle = function(Z){ 
  n = length(Z)
  sig2 = ( t(Z)%*%Q%*%Z )/n 
  return(-(n/2)*log(sig2))
}


lmin <- function(lambda, Y){
  n = length(Y)
  term1 <-Lmle(((Y^lambda) -1)/lambda)
  print(term1)
  term2 <- (lambda - 1) * sum(log(abs(Y)))
  print(term2)
  term3 <- -(n/2) *log(2*pi) - (n/2)
  print(term3)
  return (-(term1 + term2 + term3))
}
lambda = seq(0,2,0.001)
Vlmin = Vectorize(lmin,"lambda")
plot(lambda, Vlmin(lambda,Y), xlab= "lambda" , ylab = "-Lmax",type = "l")
#Question 3)
resopt = nlm(lmin,Y=Y,p=c(0.4),hessian=TRUE)
resopt
resopt$estimate
variance_lambda = solve(resopt$hessian)

#Question 4
borne_inf = resopt$estimate - sqrt(solve(resopt$hessian)*solve(t(XX)%*%XX))*qt(0.975, df = 49)
borne_sup = resopt$estimate + sqrt(solve(resopt$hessian)*solve(t(XX)%*%XX))*qt(0.975, df = 49)
# test 1
W_test1 = t(coef(modelY)-c(5,1))%*%solve(vcov(modelY))%*%(coef(modelY)-c(5,1))
pval1 = 1-pchisq(W_test1,2)
#test 2
Y_square = sign(Y)*sqrt(abs(Y))
modelY_square = lm(Y_square~XX)
W_test2 = t(coef(modelY_square)-c(5,1))%*%solve(vcov(modelY_square))%*%(coef(modelY_square)-c(5,1))
pval2 = 1-pchisq(W_test2,2)
# test3
Y_0.3 = ((Y^0.3) -1)/0.3
modelY_0.3 = lm(Y_0.3~XX)
W_test3 = t(coef(modelY_0.3)-c(5,1))%*%solve(vcov(modelY_0.3))%*%(coef(modelY_0.3)-c(5,1))
pval3 = 1-pchisq(W_test3,2)
#test 4
Y_log = log(Y)
modelY_log = lm(Y_log~XX)
W_test4 = t(coef(modelY_log)-c(5,1))%*%solve(vcov(modelY_log))%*%(coef(modelY_log)-c(5,1))
pval4 = 1-pchisq(W_test4,2)



#Question 5
# test 1
M = Y - a -b*XX
L_test1 = 2*(logLik(lm(M~XX),REML =TRUE)-logLik(lm(M ~ 1),REML = FALSE))
pval_lr1 = 1-pchisq(L_test1,2)
#test 2
M = Y_square - a -b*XX
L_test2 = 2*(logLik(lm(M~XX),REML =TRUE)-logLik(lm(M ~ 1),REML = FALSE))
pval_lr2 = 1-pchisq(L_test2,2)
# test3
M = Y_0.3 - a -b*XX
L_test3 = 2*(logLik(lm(M~XX),REML =TRUE)-logLik(lm(M ~ 1),REML = FALSE))
pval_lr3 = 1-pchisq(L_test3,2)
#test 4
M = Y_log - a -b*XX
L_test4 = 2*(logLik(lm(M~XX),REML =TRUE)-logLik(lm(M ~ 1),REML = FALSE))
pval_lr4 = 1-pchisq(L_test4,1)
#Question 6
valopt = powerTransform(modelY,family = "bcPower")
valopt



#Partie 3
NbCycleRupture <- read.csv("C:/Users/romain/NbCycleRupture.csv", sep=";")
dim(NbCycleRupture)

# Pour simplifier notation:
x1 <- NbCycleRupture$x1
x2 <- NbCycleRupture$x2
x3 <- NbCycleRupture$x3
y  <- NbCycleRupture$y

# Faisons le fit
model_1 <- lm (y ~ x1+x2+x3)
summary(model_1)
plot(model_1)



# On applique Box-Cox au modele:
transf_boxcox = boxcox(model_1, lambda = seq(-3,3))

# On trouve la valeur exacte de lambda donnée
var_lambda <- 0
model_2 <- lm(log(y) ~ x1 + x2 + x3)
summary(model_2)
plot(model_2)


# Maintenant analysons le modele avec des variables d'ordre plus haute
model_3 <- lm(y ~ x1 + x2 + x3 + I(x1^2) + I(x2^2) + I(x3^2) + I(x1*x2) + I(x1*x3) + I(x2*x3))
summary(model_3)
plot(model_3)

# COMPARAISON DES MODELES M1 ET M2
table_anova <- anova(model_1, model_3)


# On applique la meme Box-Cox que pour le modele 1 dans le modele 2:
model_4 <- lm(log(y) ~ x1 + x2 + x3 + I(x1^2) + I(x2^2) + I(x3^2) + I(x1*x2) + I(x1*x3) + I(x2*x3))
summary(model_4)
plot(model_4)

anova(model_2,model_4)

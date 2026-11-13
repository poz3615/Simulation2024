################################################
## R code to produce results in the Supplemental
## materials, Section S2 in the paper 
## "Best Practices for Modeling Temperature-Dependent 
## Arthropod Lifetimes in a Bayesian Framework: 
## A simulation study" by Zimmerman and Johnson
##
## code by Leah R Johnson, Virginia Tech, 2026
##
################################################



## setting up, numbers of samples
N<-100

## First simulations
## median around 25
set.seed(5432) # for reproducing the figures in the manuscript
my.shp<-1.75
my.med<-25

LTs25exp<-rexp(N, rate=log(2)/my.med)
LTs25w<-rweibull(N, 
                 scale= my.med/(log(2))^(1/my.shp), 
                 shape=my.shp)

## reducing the median 
my.med2<-5
LTs5exp<-rexp(N, rate=log(2)/my.med2)
LTs5w<-rweibull(N, 
                 scale= my.med2/(log(2))^(1/my.shp), 
                 shape=my.shp)


png(filename = "supmat-sim-example1.png")
par(mfrow=c(2,2), mar=c(2.5, 3, 3, 0.5) + 0.1)
hist(LTs25exp, main="Exponential: M=25",
     ylab="")
#abline(v=mean(LTs25exp), col=2)
hist(LTs25w, main="Webull: M=25",
     ylab="")
#abline(v=mean(LTs25w), col=2)

hist(LTs5exp, main="M=5",
     ylab="")
#abline(v=mean(LTs5exp), col=2)
hist(LTs5w, main="M=5",
     ylab="")
#abline(v=mean(LTs5w), col=2)
dev.off()

library(kableExtra)
kable(
  cbind(c("Exp: M=25", 
      "Weibull: M=25",
      "Exp: M=5", 
      "Weibull: M=5"),
      rbind(
        round(quantile(LTs25exp, probs=c(0.025, 0.5, 0.9725)),2),
        round(quantile(LTs25w, probs=c(0.025, 0.5, 0.9725)),2),
        round(quantile(LTs5exp, probs=c(0.025, 0.5, 0.9725)),2),
        round(quantile(LTs5w, probs=c(0.025, 0.5, 0.9725)),2)
        )
      )
  , format = "latex")

## simulation about transformation
## of data

n<-1000
my.med<-5
imLT.e<-miLT.e<-imLT.w<-miLT.w<-NULL
for(i in 1:n){
  #temp<-rnorm(N, mean=10, sd=2)#rate=1/10)
  temp.e<-rexp(N, rate=log(2)/my.med)
  temp.w<- rweibull(N, 
                  scale= my.med/(log(2))^(1/my.shp), 
                  shape=my.shp)
  if(i==1){
    imLT.e<-1/mean(temp.e)
    miLT.e<-mean(1/temp.e)
    imLT.w<-1/mean(temp.w)
    miLT.w<-mean(1/temp.w)
  }else{
    imLT.e <- c(imLT.e, 1/mean(temp.e))
    miLT.e <- c(miLT.e, mean(1/temp.e))
    imLT.w <- c(imLT.w, 1/mean(temp.w))
    miLT.w <- c(miLT.w, mean(1/temp.w))
  }
}

## estimates of the rate param
## across simulations
c(mean(imLT.e), mean(miLT.e),
  mean(imLT.w), mean(miLT.w))
c(sd(imLT.e), sd(miLT.e),
  sd(imLT.w), sd(miLT.w))

## estimated lifetimes derived across 
## many rate estimates
c(mean(1/imLT.e), mean(1/miLT.e),
  mean(1/imLT.w), mean(1/miLT.w))
c(sd(1/imLT.e), sd(1/miLT.e),
  sd(1/imLT.w), sd(1/miLT.w))

### comparison between mortality rate

png(filename = "supmat-sim-sim1.png")
par(mfrow=c(2,2), mar=c(2.5, 3, 3, 0.5) + 0.1)
hist(imLT.e, main="IoM - Exp",
     ylab="")
hist(miLT.e, main="MoI - Exp",
     ylab="")
hist(imLT.w, main="IoM - W",
     ylab="")
hist(miLT.w, main="MoI - W",
     ylab="")
dev.off()


png(filename = "supmat-sim-sim3.png")
par(mfrow=c(2,2), mar=c(2.5, 3, 3, 0.5) + 0.1)
hist(imLT.e, main="IoM - Exp",
     ylab="")
hist(miLT.e, main="MoI - Exp",
     ylab="")
hist(1/imLT.e, main="1/IoM - Exp",
     ylab="")
hist(1/miLT.e, main="1/MoI - Exp",
     ylab="")
dev.off()

par(mfrow=c(2,2))
hist(imLT.w, main="imLT - W",
     ylab="")
hist(miLT.w, main="miLT - W",
     ylab="")
hist(1/imLT.w, main="1/imLT - W",
     ylab="")
hist(1/miLT.w, main="1/miLT - W",
     ylab="")

## simulation about transformation
## of data

my.shp2<-5
imLT.w2<-miLT.w2<-NULL
for(i in 1:n){
  temp.w<- rweibull(N, 
                    scale= my.med/(log(2))^(1/my.shp), 
                    shape=my.shp2)
  if(i==1){

    imLT.w2<-1/mean(temp.w)
    miLT.w2<-mean(1/temp.w)
  }else{
    imLT.w2 <- c(imLT.w2, 1/mean(temp.w))
    miLT.w2 <- c(miLT.w2, mean(1/temp.w))
  }
}

## estimates of the rate param
## across simulations
c(mean(imLT.e), mean(miLT.e),
  mean(imLT.w), mean(miLT.w),
  mean(imLT.w2), mean(miLT.w2))
c(sd(imLT.e), sd(miLT.e),
  sd(imLT.w), sd(miLT.w),
  sd(imLT.w2), sd(miLT.w2))

## estimated lifetimes derived across 
## many rate estimates
c(mean(1/imLT.e), mean(1/miLT.e),
  mean(1/imLT.w), mean(1/miLT.w),
  mean(1/imLT.w2), mean(1/miLT.w2))
c(sd(1/imLT.e), sd(1/miLT.e),
  sd(1/imLT.w), sd(1/miLT.w),
  sd(1/imLT.w2), sd(1/miLT.w2))


png(filename = "supmat-sim-sim2.png")
par(mfrow=c(2,2), mar=c(3, 3, 4, 0.5) + 0.1)
hist(imLT.w, main="IoM - k=1.75",
     ylab="")
hist(miLT.w, main="MoI - k=1.75",
     ylab="")
hist(imLT.w2, main="IoM - k=3",
     ylab="", xlab="")
hist(miLT.w2, main="MoI - k=3",
     ylab="", xlab="")
dev.off()



hist(1/imLT.w2, main="1/imLT - W, sym",
     ylab="", xlab="")
hist(1/miLT.w2, main="1/miLT - W, sym",
     ylab="", xlab="")


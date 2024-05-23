require(deSolve)
require(ggplot2)
require(patchwork)
cols <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#66CC99")

### The aim of this code is to make a figure that illustrates
### the isolation-only strategy in terms of Hansen & Day
### in terms of the public health terminology that I defined
#### elimination, suppression, and circuit breaker.


##### Parameters
# max daily rate of community member isolation
u1max <- 0.7
# max daily rate of traveler isolation
u2max <- 0
# constraint on community member isolation
C1max <- 500
# constraint on traveller isolation
C2max <- 100
# transmission rate
beta<-.0002
# importation rate
theta<-2
theta<-0
mu<-0.334
# relative transmissibility of travellers
c<-1
T<-70
# daily vaccination rate
rho<-0
# start time for community isolation
tstart1<-0
# in case we wanted to do two switches
toff1 <- 1000 # currently set to not turn-off
# start time for restrictions on travellers
tstart2<-0


# community isolation function
u1fun<-function(t){
  u1 = u1max
  if(t<tstart1){
  u1 = 0}
  if(t<(toff1+10) & t>toff1){
    u1 = 0
  }
  return(u1)
}

# travel restrictions function
u2fun<-function(t){
  if(t>tstart2){
    u2 = u2max}
  else{
    u2 = 0
  }
  return(u2)
}

# This is designed to trigger termination of the simulation when I1>0.5
rootfun <- function(t, y, parameters){
  I1 = y[2]
  # However, with importations, it doesn't make sense to have this as
  # the endpoint of the simulation, so instead I just set Imin = 0 and have
  # the simulation end at T=100.
  y1 = 1 - I1 
  return(y1)
}


SI<-function(t, y, parameters){
  S = y[1]
  I1 = y[2]
  I2 = y[3]
  C1 = y[4]
  C2 = y[5]
  iso = y[6]
  u1 = u1fun(t)
  u2 = u2fun(t)
  
  if(C1 > C1max){
    u1 = 0
  }
  if(C2 > C2max){
    u2 = 0
  }
  
  dS = - beta*S*(I1 + c*I2) - rho*S
  dI1 = beta*S*(I1 + c*I2) - mu*I1 - u1*I1
  dI2 = theta - 2*mu*I2 - u2*I2
  dC1 = u1*I1
  dC2 = u2*I2
  diso = u1*I1 - mu*iso
  return(list(c(dS,dI1,dI2,dC1,dC2,diso)))
}

### ELIMINATION
out <- ode(y = c(S=5000,I1=10,I2=0,C1=0,C2=0,iso=0), parms = NULL, times = seq(0,T,.05), func = SI, events = list(func = NULL, root = TRUE, terminalroot = 1), rootfun=rootfun)
out<- data.frame(out)
J = cumsum(beta*out$S*(out$I1+c*out$I2)*diff(c(0,out$time)))
J.elim <- J
u = rep(0,length(J))
u[which(diff(c(out$C1[1],out$C1))*diff(c(0,out$time))>0.001)] = max(c(out$I1,out$iso))
out.elim <- data.frame(out,u=u)

g.elim = ggplot(out.elim, aes(x = time, y = I1)) +
  geom_ribbon(aes(ymin = 0, ymax = u), fill = cols[7], alpha = 0.2)+
  geom_line(color = cols[7], size = 2)+
  geom_line(aes(y=iso))+
  ggtitle("Elimination")+ylab("Community prevalence") + xlim(0,T)+theme_classic()

##### SUPPRESSION (same parameters except a lower u1max)
u1max <- 0.6
out <- ode(y = c(S=5000,I1=10,I2=0,C1=0,C2=0,iso=0), parms = NULL, times = seq(0,T,.05), func = SI, events = list(func = NULL, root = TRUE, terminalroot = 1), rootfun=rootfun)
out <- data.frame(out)
out.supp<-data.frame(out)
J = cumsum(beta*out$S*(out$I1+c*out$I2)*diff(c(0,out$time)))
J.supp <- J
u = rep(0,length(J))
u[which(diff(c(out$C1[1],out$C1))*diff(c(0,out$time))>0.001)] = max(c(out$I1,out$iso))
out.supp <- data.frame(out,u=u)

g.supp = ggplot(out.supp, aes(x = time, y = I1)) +
  geom_ribbon(aes(ymin = 0, ymax = u), fill = cols[2], alpha = 0.2)+
  geom_line(color = cols[2], size = 2)+
  geom_line(aes(y=iso))+ggtitle("Suppression")+ylab("") + xlim(0,T)+theme_classic()

###### CIRCUIT BREAKER (a lower value of the constraint)
### (also times of on and off to illustrate that timing doesn't matter)
C1max <- 400
toff1<-10
out <- ode(y = c(S=5000,I1=10,I2=0,C1=0,C2=0,iso=0), parms = NULL, times = seq(0,T,.05), func = SI, events = list(func = NULL, root = TRUE, terminalroot = 1), rootfun=rootfun)
out <- data.frame(out)
out.circ<-out
J = cumsum(beta*out$S*(out$I1+c*out$I2)*diff(c(0,out$time)))
u = rep(0,length(J))
u[which(diff(c(out$C1[1],out$C1))*diff(c(0,out$time))>0.001)] = max(c(out$I1,out$iso))
out.circ <- data.frame(out,u=u,J=J)
print(tail(out$C1,1))

g.circ = ggplot(out.circ, aes(x = time, y = I1)) +
  geom_ribbon(aes(ymin = 0, ymax = u), fill = cols[4], alpha = 0.2)+
  geom_line(color = cols[4], size = 2)+geom_line(aes(y=iso))+ggtitle("Circuit breaker 1")+ylab("Community prevalence") + xlim(0,T) + theme_classic()


#### A different circuit-breaker with the same parameters as above
tstart1<-5
toff1<-1000
out <- ode(y = c(S=5000,I1=10,I2=0,C1=0,C2=0,iso=0), parms = NULL, times = seq(0,T,.05), func = SI, events = list(func = NULL, root = TRUE, terminalroot = 1), rootfun=rootfun)
out <- data.frame(out)
out.circ.2<-out
print(tail(out$C1,1))
J = cumsum(beta*out$S*(out$I1+c*out$I2)*diff(c(0,out$time)))
u = rep(0,length(J))
u[which(diff(c(out$C1[1],out$C1))*diff(c(0,out$time))>0.001)] = max(c(out$I1,out$iso))
out.circ2 <- data.frame(out,u=u,J=J)
print(tail(out$C1,1))

g.circ.2 = ggplot(out.circ.2, aes(x = time, y = I1)) +
  geom_ribbon(aes(ymin = 0, ymax = u), fill = cols[6], alpha = 0.2)+
  geom_line(color = cols[6], size = 2)+
  geom_line(aes(y=iso))+
  ggtitle("Circuit breaker 2")+ylab("") + xlim(0,T) + theme_classic()



J.circ = ggplot()+
  geom_ribbon(data=out.circ,aes(x=time,ymin = 0, ymax = u*(max(J)/max(I1))), fill = cols[4], alpha = 0.2)+
  geom_ribbon(data=out.circ.2,aes(x=time,ymin = 0, ymax = u*(max(J)/max(I1))), fill = cols[6], alpha = 0.2)+
  geom_line(data=out.circ, aes(x = time, y = J),color = cols[4], size = 2)+
  geom_line(data=out.circ.2, aes(x = time, y = J),color = cols[6], size = 2)+ylab("Cumulative cases")+ggtitle("Circuit breaker 1 & 2 equivlance")+theme_classic()

g1=g.elim+g.supp
g2=g.circ+g.circ.2
fig = g1/g2/J.circ+
  plot_annotation(tag_levels = 'A')
ggsave("~/Desktop/Work/Students/MSc/Adu-Boahen/Thesis_2024-Updates/Amys codes/figures/Isolation_only_no_importations.png", width = 20, height = 15, units = "cm")

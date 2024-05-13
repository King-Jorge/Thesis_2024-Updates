require(deSolve)
# max daily rate of community member isolation
u1max <- .3
# max daily rate of traveller isolation
u2max <- 0
# constraint on community member isolation (this is a large constraint for an
# essentially unconstrained problem)
C1max <- 4000
# constraint on traveller isolation
C2max <- 500
# transmission rate
beta<-.0003
# importation rate (high for illustrative purposes)
theta<-3
# While the mu's aren't evenly spaced, these mu's give good
# spacing on the x-axis variables in the plot.
muvec <- 1/seq(5,0.1,-0.1)
# relative transmissibility of travellers
c<-1

T<-100

# daily vaccination rate
rho<-0 # no vaccination

# start time for community isolation
tstart1<-0

# community isolation function - assumes that u1 = umax until the constraint is
# used up. We know that this is not a unique optimal control, but it is one
# optimal control and therefore sufficient to demonstrate the relationships.
u1fun<-function(t){
  if(t>tstart1){
  u1 = u1max}
  else{
    u1 = 0
  }
  return(u1)
}

# travel restrictions function
u2fun<-function(t){
    u2 = 0
  return(u2)
}


# This is designed to trigger termination of the simulation when I1>0.5
rootfun <- function(t, y, parameters){
  I1 = y[2]
  C1 = y[4]
  y1 = 0.5 - I1 
  return(y1)
}

# An event function that does nothing because termination is the event
eventfun <- function(t, y, pars) {
  return(y)
}

SI<-function(t, y, parameters){
  S = y[1]
  I1 = y[2]
  I2 = y[3]
  C1 = y[4]
  C2 = y[5]
  u1 = u1fun(t)
  u2 = u2fun(t)
  
  if(C1 > C1max){
    u1 = 0
  }
  
dS = - beta*S*(I1 + c*I2) - rho*S
dI1 = beta*S*(I1 + c*I2) - mu*I1 - u1*I1
dI2 = theta - 2*mu*I2 - u2*I2
dC1 = u1*I1
dC2 = u2*I2
  return(list(c(dS,dI1,dI2,dC1,dC2)))
}

output = NULL
for(i in seq(1,length(muvec))){
mu<-muvec[i]
out <- ode(y = c(S=5000,I1=10,I2=0,C1=0,C2=0), parms = NULL, times = seq(0,T,.05), func = SI, events = list(func = NULL, root = TRUE, terminalroot = 1), rootfun = rootfun)
out <- data.frame(out)
# This is the correct objective function with importations
J = sum(beta*out$S*(out$I1+c*out$I2)*diff(c(0,out$time)))
output = data.frame(rbind(output,c(mu=mu,J=J,w_T = tail(out$C1,1), S_T = tail(out$S,1), I_T = tail(out$I1,1), T = tail(out$time,1))))
}

plot(output$w_T,output$S, typ="l", ylim = c(0,5000), ylab="S(T)", xlab="number isolated, w_max")

# repeat with  constraint
C1max = 1450

output2 = NULL
for(i in seq(1,length(muvec))){
  mu<-muvec[i]
  out <- ode(y = c(S=5000,I1=10,I2=0,C1=0,C2=0), parms = NULL, times = seq(0,T,.05), func = SI, events = list(func = NULL, root = TRUE, terminalroot = 1), rootfun = rootfun)
  out <- data.frame(out)
  # This is the correct objective function with importations
  J = sum(beta*out$S*(out$I1+c*out$I2)*diff(c(0,out$time)))
  output2 = data.frame(rbind(output2,c(mu=mu,J=J,w_T = tail(out$C1,1), S_T = tail(out$S,1), I_T = tail(out$I1,1), T = tail(out$time,1))))
}

points(output2$w_T,output2$S, pch=19, ylim = c(0,5000))

# Parameters from the lower panel of Hansen and Day, Figure 1.
u1max <-1
C1max <- 4000 # unconstrained.
output = NULL
for(i in seq(1,length(muvec))){
  mu<-muvec[i]
  out <- ode(y = c(S=5000,I1=10,I2=0,C1=0,C2=0), parms = NULL, times = seq(0,T,.05), func = SI, events = list(func = NULL, root = TRUE, terminalroot = 1), rootfun = rootfun)
  out <- data.frame(out)
  # This is the correct objective function with importations
  J = sum(beta*out$S*(out$I1+c*out$I2)*diff(c(0,out$time)))
  output = data.frame(rbind(output,c(mu=mu,J=J,w_T = tail(out$C1,1), S_T = tail(out$S,1), I_T = tail(out$I1,1), T = tail(out$time,1))))
}

plot(output$w_T,output$S_T, typ="l", ylim = c(0,5000), ylab="S(T)", xlab="number isolated, w_max")

# repeat with  constraint
C1max <- 1450

output2 = NULL
for(i in seq(1,length(muvec))){
  mu<-muvec[i]
  out <- ode(y = c(S=5000,I1=10,I2=0,C1=0,C2=0), parms = NULL, times = seq(0,T,.05), func = SI, events = list(func = NULL, root = TRUE, terminalroot = 1), rootfun = rootfun)
  out <- data.frame(out)
  # This is the correct objective function with importations
  J = sum(beta*out$S*(out$I1+c*out$I2)*diff(c(0,out$time)))
  output2 = data.frame(rbind(output2,c(mu=mu,J=J,w_T = tail(out$C1,1), S_T = tail(out$S,1), I_T = tail(out$I1,1), T = tail(out$time,1))))
  if(i==1){
    out2 = out
  }
}

points(output2$w_T,output2$S_T, ylim = c(0,5000),pch=19)

plot(out$time, out$I1, typ="l")
lines(out$time, 10*out$I2, typ="l", col="red")
plot(out2$time, out2$I1, typ="l")
lines(out2$time, 10*out2$I2, typ="l", col = "red")
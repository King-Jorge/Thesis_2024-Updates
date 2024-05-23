require(deSolve)
# max daily rate of community member isolation
u1max <-0.3
u1max <- 0.9

# max daily rate of traveler isolation
u2max <- 10

# constraint on community member isolation
C1max <- 600

# constraint on traveller isolation
C2max <- 100

# transmission rate
beta<-.0002
# importation rate
theta<-2
mu<-0.334
# relative transmissibility of travellers
c<-1
T<-100

# daily vaccination rate
rho<-0

# start time for community isolation
tstart1<-0
# in case we wanted to do two switches
toff1 <- 1000 # currently set to not turn-off
# an example of another 
ton <- toff1+5

# start time for restrictions on travellers
tstart2<-0


# community isolation function
u1fun<-function(t){
  if(toff1>t & t>tstart1){
  u1 = u1max}
  else if (t>ton){
    u1 = u1max
  }
  else{
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
  return(list(c(dS,dI1,dI2,dC1,dC2)))
}

out <- ode(y = c(S=5000,I1=10,I2=3,C1=0,C2=0), parms = NULL, times = seq(0,T,.05), func = SI, events = list(func = NULL, root = TRUE, terminalroot = 1), rootfun=rootfun)
out <- data.frame(out)

plot(out$time, out$I1, typ = "l")
# Travellers are x50 to allow for visualization
lines(out$time, 50*out$I2, typ = "l", col="red")
# visualization to see that all resources are used up
plot(out$time, out$C1, typ = "l")
plot(out$time, out$C2, typ = "l")

# This is the correct objective function with importations
J = sum(beta*out$S*(out$I1+c*out$I2)*diff(c(0,out$time)))
# Print the objective function. Try different tstart1 and tstart2.
# When there is no vaccintion, the objective function value is the
# same, provided all resources are used.
print(c(tstart1,tstart2,J,tail(out$C1,1),tail(out$C2,1)))


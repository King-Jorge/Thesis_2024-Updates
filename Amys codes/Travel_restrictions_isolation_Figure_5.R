require(deSolve)
require(ggplot2)
require(patchwork)
require(ggpubr)
cols <- c("#D55E00","#E69F00","#F0E442","#66CC99","#009E73","white", "#CC79A7","#999999","#56B4E9")

##### Parameters
# constraint on community member isolation
C1max <- 1000
# constraint on traveller isolation
C2max <- 30
# transmission rate
beta<-.0002
# importation rate
theta<-1
mu<-0.334
# relative transmissibility of travellers
c<-1
T<-1000

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

  if(C1 > C1max){
    u1 = 0
  }
  if(C2 > C2max){
    u2 = 0
  }

  dS = - beta*S*(I1 + c*I2)
  dI1 = beta*S*(I1 + c*I2) - mu*I1 - u1*I1
  dI2 = theta - 2*mu*I2 - u2*I2
  dC1 = u1*I1
  dC2 = u2*I2
  return(list(c(dS,dI1,dI2,dC1,dC2)))
}
stratMx = NULL
u1vec <- seq(2,0,-.025)
u2vec <- seq(3,0,-.025)

for(i in seq(1,length(u1vec))){
  u1 = u1vec[i]
  for(j in seq(1,length(u2vec))){
    u2 = u2vec[j]
  out <- ode(y = c(S=5000,I1=10,I2=1,C1=0,C2=0), parms = NULL, times = seq(0,T,1), func = SI, events = list(func = NULL, root = TRUE, terminalroot = 1), rootfun=rootfun, method = "radau")
  out<- data.frame(out)
if(max(out$I1)<=10 & max(out$time)<T & max(out$C1)<C1max & max(out$C2)<C2max){
  outcome = "elimination"
} else if(max(out$I1)>10 & max(out$C1)<C1max & max(out$C2)>=C2max & max(out$time)<T){
  outcome = "supp/circ"
} else if(max(out$I1)>10 & max(out$C1)>=C1max & max(out$C2)<C2max& max(out$time)<T){
  outcome = "circ/supp"
} else if(max(out$I1)>10 & max(out$C1)<C1max & max(out$C2)<C2max& max(out$time)<T){
  outcome = "suppression"
} else if(max(out$time)==T){
  outcome = "no end"
} else if(max(out$I1)>10 & max(out$C1)>=C1max & max(out$C2)>=C2max& max(out$time)<T){
  outcome = "circuit breaker"
  }

J = cumsum(beta*out$S*(out$I1+c*out$I2)*diff(c(0,out$time)))
stratMx = rbind(stratMx, data.frame(u1 = u1, u2 = u2, outcome = outcome, J=J, T = max(out$time)))
}}

g1=ggplot(stratMx, aes(u1, u2, fill= outcome)) +
  geom_tile()+xlab("community daily max, u1max")+ylab("traveler daily max, u2max")+
  ggtitle("Low importations")+
  scale_fill_manual(values=cols, breaks = c("elimination", "suppression", "supp/circ", "circ/supp", "circuit breaker"))+theme_classic()+theme(legend.title=element_blank())

theta <- 2
stratMx2 = NULL
for(i in seq(1,length(u1vec))){
  u1 = u1vec[i]
  for(j in seq(1,length(u2vec))){
    u2 = u2vec[j]
    out <- ode(y = c(S=5000,I1=10,I2=1,C1=0,C2=0), parms = NULL, times = seq(0,T,1), func = SI, events = list(func = NULL, root = TRUE, terminalroot = 1), rootfun=rootfun, method = "radau")
    out<- data.frame(out)
    if(max(out$I1)<=10 & max(out$time)<T & max(out$C1)<C1max & max(out$C2)<C2max){
      outcome = "elimination"
    } else if(max(out$I1)>10 & max(out$C1)<C1max & max(out$C2)>=C2max & max(out$time)<T){
      outcome = "supp/circ"
    } else if(max(out$I1)>10 & max(out$C1)>=C1max & max(out$C2)<C2max& max(out$time)<T){
      outcome = "circ/supp"
    } else if(max(out$I1)>10 & max(out$C1)<C1max & max(out$C2)<C2max& max(out$time)<T){
      outcome = "suppression"
    } else if(max(out$time)==T){
      outcome = "no end"
    } else if(max(out$I1)>10 & max(out$C1)>=C1max & max(out$C2)>=C2max& max(out$time)<T){
      outcome = "circuit breaker"
    }

    stratMx2 = rbind(stratMx2, data.frame(u1 = u1, u2 = u2, outcome = outcome))
  }}

g2=ggplot(stratMx2, aes(u1, u2, fill= outcome)) +
  geom_tile()+xlab("community daily max, u1max")+ylab("traveler daily max, u2max")+
  ggtitle("High importations")+
  scale_fill_manual(values=cols, breaks = c("elimination", "suppression", "supp/circ", "circ/supp", "circuit breaker"))+theme_classic()+theme(legend.title=element_blank())

g3=ggarrange(g1,g2,common.legend = TRUE, legend="right")

g4=ggplot(stratMx, aes(u1, u2, fill= T)) +
  geom_tile()+xlab("community daily max, u1max")+ylab("traveler daily max, u2max")+
  ggtitle("Duration of outbreak (low importations)")+
  scale_fill_gradient(low = cols[7], high = "black")+theme_classic()+theme(legend.title=element_blank())

g5=ggplot(stratMx, aes(u1, u2, fill= J)) +
  geom_tile()+xlab("community daily max, u1max")+ylab("traveler daily max, u2max")+
  ggtitle("New cases (low importations)")+
  scale_fill_gradient(low = cols[9], high = "black")+theme_classic()+theme(legend.title=element_blank())

g7=g3/(g4+g5)+plot_annotation(tag_levels = 'A')+plot_layout(heights = c(1.5, 1))

ggsave("~/Desktop/Work/Students/MSc/Adu-Boahen/Thesis_2024-Updates/Amys codes/figures/Heat_map.png", width = 30, height = 20, units = "cm")

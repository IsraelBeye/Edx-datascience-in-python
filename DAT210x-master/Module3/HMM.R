
# HMM implementation for the data found in
#https://www.12000.org/my_notes/HMM_program/report/report.htm#x1-60003

len <- 5; 
states<- c("S", "T")
observations<-c("A","C","T")
#-----------------------------------------------------------
#function to calculate transition matrix
#-----------------------------------------------------------
tm <- matrix(data=0, ncol=2, nrow=2);
rownames(tm)<-states
colnames(tm)<-states

#from S
tm[1,1] <- 0.7;tm[1,2] <- 0.3;
#from T
tm[2,1] <- 0.4;  tm[2,2] <- 0.6;
#tm[2,1] <- 0.3;  tm[2,2] <- 0.7
tm

#-----------------------------------------------------------
#function to calculate emission probabilities
#-----------------------------------------------------------

em <- matrix(0, nrow=2, ncol=3);
rownames(em)<-states
colnames(em)<-observations

#em[1,1] <- 0.5; em[1,2] <- 0.4;  em[1,3] <- 0.1;
#em[2,1] <- 0.1; em[2,2] <- 0.3; em[2,3] <- 0.6;
em[1,1] <- 0.4; em[1,2] <- 0.4;  em[1,3] <- 0.2;
em[2,1] <- 0.25; em[2,2] <- 0.55; em[2,3] <- 0.2;
em


#-----------------------------------------------------------
#prepare storage
sta <- numeric(len);
obs <- numeric(len);


#choose random starting position
#prior <- c(0.6, 0.4);
prior <- c(0.4, 0.6);
sta[1] <- sample(1:2, 1, prob=prior);
obs[1] <- sample(1:3, 1, prob=em[sta[1],]);


#run MC
for(i in 2:len){
  sta[i] <- sample(1:2, 1, prob=tm[sta[i-1],]);
  obs[i] <- sample(1:3, 1, prob=em[sta[i],]);
}
#return
data.frame(sta=sta, obs=obs);

#-----------------------------------------------------------
#Function to calculate LL
#-----------------------------------------------------------
forward <- matrix(data=0, nrow=2, ncol=len);

forward[,1] <- em[,obs[1]] * prior;
s <- sum(forward[,1]);
forward[,1] <- forward[,1] / s;
scale <- log10(s);

#now make recursion    
for(i in 2:len){
  #swap entries
  #calc forward
  forward[,i] <- (forward[,i-1] %*% tm) * em[,obs[i]];
  
  s <- sum(forward[,i]);
  forward[,i] <- forward[,i] / s; #calculates the probability of being in state si at time t given observing observation o1to ot inclusive
  #p(St=si,O1:t(forward[,i] or alpha)/s(sum over all the forward prob at time i)
  scale <- scale + log10(s);    
}
list("forward"= forward, "scale"=scale, "em"=em, "tm"=tm)

backward <- matrix(data=0, nrow=2, ncol=len);
gamma<-matrix(data=0, nrow=2, ncol=len)

#Start with initial  

backward[,len] <- 1;  
#now make recursion    
for(i in (len-1):1){
  
  #calc backward
  backward[,i] <- (backward[,i+1]*em[,obs[i+1]])%*% tm ;
  s=sum(backward[,i])
  backward[,i] <- backward[,i] / s;
  gamma[,i]<-forward[,i]*backward[,i]
}

list("forward"= forward, "backward"=backward, "gamma"=gamma, "scale"=scale, "em"=em, "tm"=tm)



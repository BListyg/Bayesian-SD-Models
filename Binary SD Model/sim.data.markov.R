sim.data.markov <- function(trans.mat,chain.length){
#Generates a markov chain with transition matrix, trans.mat, and length, chain.length
#See sim.data.markov.example for an example of usage
#Author:Jeffrey Annis
#Date:3.29.14  

  num.states = length(trans.mat[,1])
  init.probs = rep(1,num.states)/num.states #all states have eq. prob at t=1
  state = NULL
  state <- sample(1:num.states,1,prob=init.probs) #get initial state
  
  for(i in 1:(chain.length-1)){
    trans.probs <- trans.mat[state[i],] #get the transition probabilities
    state[i+1] <- sample(1:num.states,1,prob=trans.probs) #sample next state based on transition probs
  }

  return(state)
  
}


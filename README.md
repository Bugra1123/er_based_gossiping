# Randomized Gossiping with Effective Resistance Weights: Asynchronous Gossiping Experiments

This repository contains the codes to obtain the comparison results on the asynchronous gossiping algorithms presented in the Numerical Experiments section of the paper [[Can et al, 2019]](https://arxiv.org/abs/1907.13110) which is to appear in IEEE TCNS at the moment this document was prepared.  The experiments compare the performance of effective resistance based asynchronous gossiping algorithm with the Fastest Mixing Markov Chain probabilities (FMMC), Fastest Quantum Gossiping (FQG), asynhronous gossiping with Metropolis based probabilities, and asynchronous randomized uniform gossiping on barbell graphs, small-world graphs, and stochastic block models  (SBM) (see the paper for further details on the graphs and the gossiping algorithms). 

# ER based gossiping vs FMMC on barbell graph and small-world graph

We first compute the effective resistance based probabilities using decentralized normalized Kaczmarz algorithm proposed in [[Aybat and Gurbuzbalaban, 2017]](https://ieeexplore.ieee.org/document/8308701) and collect wake-up and communication  probabilities together with computation times and number of iterations to achieve pre-determined accuracy (_TimeNormKacmarz_Barbell.m_,_TimeNormKacmarz_SmallWorld.m_). Similarly, we find the same quantities for decentralized subgradient method proposed in [[Boyd et al., 2006]](Ju1qrYW8_LYAAAAA:zotP08JQaiRltEt3Hj1R6ySmsxbxBdVwiGsHSCWtinY5m7_4t7iENAz5DXCriYMQ8Ds20bn1lGM) (_TimeFMC_v3_Barbell.m_ and _TimeFMC_v3_smallWorld.m_). 
Upon computing these probabilities, we did asynchronous gossiping on barbell graphs and small-word graphs then store average wake-up times required for algorithm to achieve a relative error less than 0.01 (_ComparisonConsensus_v3_Barbell.m_ and _ComparisonConsensus_v3_SmallWorld.m_).

# Wake-up time comparisons between ER based gossiping with randomized uniform gossiping, and Metropolis weights based gossiping on barbell graphs, small-world graphs, and SBM.

The code _KaczmarzERProb.m_ computes the effective resistance based wake-up and communication probabilities using normalized Kaczmarz method on given graph and  _ERProb.m_ computes the exact ER-based probabilities using pseudo inverse of Laplacian. The wake-up and communication probabilities of Metropolis based gossiping are computed using _METProb.m_ and for the probabilites of uniform gossiping we directly used the graph adjecancy and diagonal degree matrix. The codes _GossipingComparisonSmallWorld.m_, _GossipingComparisonSBMk.m_, and _GossipingComparisonBarbell.m_ calculate the mean and the standard deviation of number of wake-ups required to achieve the 10^{-6} accuracy over 250 scenarios on small-world, SBM, and barbell graphs. 

# Comparison of spectral gaps of expected iteration matrices of ER based gossiping, Metropolis based gossiping, and FQG on barbell graphs and SBM
In this part, we compute the spectral gaps of expected communication matrices implied by ER based probabilities, Metropolis based probabilities, and the probabilities of FQG on barbell graphs and SBM. The algorithm FQGopt.m finds the communication probabilities of FQG by solving SDP problem given in [[Jafarizadeh, 2019]](Ju1qrYW8_LYAAAAA:zotP08JQaiRltEt3Hj1R6ySmsxbxBdVwiGsHSCWtinY5m7_4t7iENAz5DXCriYMQ8Ds20bn1lGM), we used ERProb.m and METProb.m to compute ER based probabilities and Metropolis based probabilities. _FQGEigComparisonSBM.m_ and _FQGEigComparisonBarbell.m_ compare the spectral gap of expected iteration matrices of these randomized gossiping algorithms on SBM and barbell graphs, respectively. 

# Miscellaneous
We give more information on the functions included in the repository below. 
## SBMk
Generates a stochastic block model. 
### _Syntax_
[L,D,A,Lnorm,Connected]=SBMk(n,k,p,q)
### _Description_
[L,D,A,Lnorm,Connected]=SBMk(n,k,p,q) generates a stochastic block model consists of n nodes and k clusters where probability that two nodes belonging to samere cluster is connected is p, and the probability that two nodes belonging to seperate clusters are connected is q. Returns the Laplacian (L), diagonal degree matrix (D), adjacency matrix (A), normalized Laplacian (Lnorm), and boolean variable (Connected) that checks whether the network is connnected or not. 

## small_world
Generates a small-world graph. This code was provided by Dr. Necdet Serhat Aybat.
### _Syntax_
[L,A,d]=smallworld_graph(N,M)
### _Description_
[L,A,d]=smallworld_graph(N,M) generates a small-world graph consists of N nodes and m edges. Returns the Laplacian (L), diagonal degree matrix (d), and the adjacency matrix (A).

## barbell_graph
Generates a barbell graph. This code was provided by Dr. Necdet Serhat Aybat.
### Syntax
[L,A,d]=graph_barbell(N)
### Description
[L,A,d]=graph_barbell(N) generates a barbell graph consists of 2N nodes. Returns the Laplacian (L), diagonal degree matrix (d), and the adjacency matrix (A).

## Gossip
Asynchronous randomized gossiping algorithm. 
### Syntax
[WakeUps,Std]=Gossip(Sim,P_WakeUp,P_Comm,Dimension,Acc)
### Description
[WakeUps,Std]=Gossip(Sim,P_WakeUp,P_Comm,Dimension,Acc) runs an asynchronous randomized gossiping where each node i wakes up with probability, P_WakeUp(i), and communicates with its neighbor j with probability, P_Comm(i,j), Initially a random vector x(i,:) of the size (1,Dimension) is assigned to node i, then at each iteration node i wakes up and communicates with its neighbor j to receive variable x(j,:) and both agents make the updates 

x(i,:) = 1/2*(x(i,:)+x(j,:)), & x(j,:)= 1/2*(x(i,:)+x(j,:))

The gossiping continues until the relative error between the matrix x and the matrix whose rows are the average of the initialized values \{x(i,:)\}_{i=1,..,N} at the beginning (with respect to Frobenius norm) reaches the accuracy level Acc. This process is repeated for (Sim)-number of iterations and the mean (WakeUps) and the standard deviation (Std) of the number of wake-ups observed over these simulations are returned by the algorithm. 
   

%% Comparison of ER, Metropolis, and classic gossiping on stochastic block matrices
% Simulates the randomized gossiping algorithms on SBM with provided communication and
% wake up probabilities for given number of iterations. 

rng(15)           % Fix random seed generator for simulations.
%% Parameters
n=120;            % The number of agents in the network

ncluster=3;       % The number of clusters                       

innerprob=0.9;    % The probability two nodes belonging to same cluster are connected.

outerprob=0.05;   % The probability that two nodes not belonging to same cluster are connected.

Dimension=5;      % The dimension of the model, i.e. number of features stored at agents.

Accuracy=1e-6;    % Accuracy on the Frobenius distance between the matrix X 
                  % where x[i,:] corresponds to local variable, and the average vector x_mean=1/N sum_{i=1,..,N}x[i,:]

Sim=250;          % The total number of simulations that are going to be generated.

% Define the network
[L,D,A,Lnorm,Connected]=SBMk(n,ncluster,innerprob,outerprob); % Calculate the Laplacian, Adjecency, and degree matrices 

G=graph(A);
close all
plot(G)      % Display the network

% Define the set of nodes.
Nodes=1:n;


%% Kacmarz ER gossiping
% Calculate the relevant probabilities
[PRes_Kac_Com,PRes_Kac_WakeUp,IterKac]=KaczmarzERProb(L,A,0.01);
% Do the gossiping
[WakeUp_KacER,Std_KacER]=Gossip(Sim,PRes_Kac_WakeUp,PRes_Kac_Com,Dimension,Accuracy);

%% ER gossiping (with exact ER probabilities)
% Calculate the relevant probabilities
[PRes_Com,PRes_WakeUp]=ERProb(L,A);
% Do the gossiping
[WakeUp_ER,Std_ER]=Gossip(Sim,PRes_WakeUp,PRes_Com,Dimension,Accuracy);

%% Uniform gossiping
% Calculate the relevant probabilities
PUni_Com= inv(diag(D))*A;
PUni_WakeUp=1/n*(1:n);
% Do the gossiping
[WakeUp_Uni,Std_Uni]=Gossip(Sim,PUni_WakeUp,PUni_Com,Dimension,Accuracy);

%% Metropolis gossiping
% Calculate the relevant probabilities
[PMet_Com,PMet_WakeUp]= METProb(A);
% Do the gossiping
[WakeUp_Met,Std_Met]=Gossip(Sim,PMet_WakeUp,PMet_Com,Dimension,Accuracy);

%% Output 
Data={'Method'    , 'WakeUps'    ,'Std','Num of Nodes','Inner Prob', 'Outer Prob','IterKac';
      'KacER'     , WakeUp_KacER , Std_KacER, n            , innerprob  ,  outerprob  , IterKac;
      'ER'        , WakeUp_ER    , Std_ER   , n            , innerprob  ,  outerprob  , 'NA';
      'Uniform'   , WakeUp_Uni   , Std_Uni  , n            , innerprob  ,  outerprob  , 'NA';
      'Metropolis', WakeUp_Met   , Std_Met  , n            , innerprob  ,  outerprob  , 'NA'};
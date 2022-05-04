%% Comparison of ER, Metropolis, and classic gossiping on small world graph
% Simulates the randomized gossiping algorithms on small world graphs with provided communication and
% wake up probabilities for given number of iterations. 


rng(15)             % Fix random seed generator for simulations.
%% Parameters 
n=100;                  % The number of agents in the network

m=floor(0.05*(n^2-n));  % The average degree of the network.    
% m=floor(0.4*(n^2-n));

Dimension=5;      % The dimension of the model, i.e. number of features stored at agents.

Accuracy=1e-6;    % Accuracy on the Frobenius distance between the matrix X 
                  % where x[i,:] corresponds to local variable, and the average vector x_mean=1/N sum_{i=1,..,N}x[i,:]

Sim=250;          % The total number of simulations that are going to be generated.


%% Define the network
[L,A,D]=smallworld_graph(n,m);

G=graph(A);
close all
plot(G,'NodeLabel',{})  % Plot the graph

% Define the set of nodes
Nodes=1:n;

%% Kacmarz ER gossiping
% Calculate the relevant probabilities
[PRes_Kac_Com,PRes_Kac_WakeUp,IterKac]=KaczmarzERProb(L,A,0.01);
% Do the gossiping
[WakeUp_KacER, Std_KacER]=Gossip(Sim,PRes_Kac_WakeUp,PRes_Kac_Com,Dimension,Accuracy);

%% ER gossiping (with exact ER probabilities)
% Calculate the relevant probabilities
[PRes_Com,PRes_WakeUp]=ERProb(L,A);
[WakeUp_ER, Std_ER]=Gossip(Sim,PRes_WakeUp,PRes_Com,Dimension,Accuracy);

%% Uniform gossiping
% Calculate the relevant probabilities
PUni_Com= inv(diag(D))*A;
PUni_WakeUp=1/n*(1:n);
% Do the gossiping
[WakeUp_Uni, Std_Uni]=Gossip(Sim,PUni_WakeUp,PUni_Com,Dimension,Accuracy);

%% Metropolis gossiping
% Calculate the relevant probabilities
[PMet_Com,PMet_WakeUp]= METProb(A);
% Do the gossiping
[WakeUp_Met, Std_Met]=Gossip(Sim,PMet_WakeUp,PMet_Com,Dimension,Accuracy);

%% Output 
Data={'Method'    , 'WakeUps'    , 'Std' ,'Num of Nodes','m', "IterKac";
      'KacER'     , WakeUp_KacER , Std_KacER  ,n            , m , IterKac;
      'ER'        , WakeUp_ER    , Std_ER     ,n            , m , 'NA';
      'Uniform'   , WakeUp_Uni   , Std_Uni    ,n            , m , 'NA';
      'Metropolis', WakeUp_Met   , Std_Met    ,n            , m , 'NA'};
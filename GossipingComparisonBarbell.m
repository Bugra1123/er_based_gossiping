%% Comparison of ER, Metropolis, and classic gossiping on barbell graph
% Simulates the randomized gossiping algorithms on barbell graph with provided communication and
% wake up probabilities for given number of iterations. 

rng(15)             % Fix random seed generator for simulations.
%% Parameters 
n=20;               % The number of agents at the complete graph part,i.e. K_n, of the barbell graph
                    % so that the total number of agents is 2n

Dimension=5;        % The dimension of the model, i.e. number of features stored at agents.

Accuracy=1e-6;      % Accuracy on the Frobenius distance between the matrix X 
                    % where x[i,:] corresponds to local variable, and the average vector x_mean=1/N sum_{i=1,..,N}x[i,:]

Sim=250;            % The total number of simulations that are going to be generated.

% Define the network
[L,A,D]=graph_barbell(n);   % Calculate the Laplacian, Adjecency, and degree matrices 
G=graph(A);
close all
plot(G)                     % Display the network

% Define the set of nodes.
Nodes=1:(2*n);      

%% Kacmarz ER gossiping
% Calculate the relevant probabilities
[PRes_Kac_Com,PRes_Kac_WakeUp,IterKac]=KaczmarzERProb(L,A,0.01);
% Do the gossiping
[WakeUp_KacER, Std_KacER]=Gossip(Sim,PRes_Kac_WakeUp,PRes_Kac_Com,Dimension,Accuracy);

%% ER gossiping (with exact ER probabilities)
% Calculate the relevant probabilities
[PRes_Com,PRes_WakeUp]=ERProb(L,A);
% Do the gossiping
[WakeUp_ER,Std_ER]=Gossip(Sim,PRes_WakeUp,PRes_Com,Dimension,Accuracy);

%% Uniform gossiping
% Calculate the relevant probabilities
PUni_Com= inv(diag(D))*A;
PUni_WakeUp=1/(2*n)*Nodes;
% Do the gossiping
[WakeUp_Uni,Std_Uni]=Gossip(Sim,PUni_WakeUp,PUni_Com,Dimension,Accuracy);

%% Metropolis gossiping
% Calculate the relevant probabilities
[PMet_Com,PMet_WakeUp]= METProb(A);
% Do the gossiping
[WakeUp_Met,Std_Met]=Gossip(Sim,PMet_WakeUp,PMet_Com,Dimension,Accuracy);

%% Output 
Data={'Method'    , 'WakeUps'    , "Std",'Num of Nodes','IterKac';
      'KacER'     , WakeUp_KacER , Std_KacER ,2*n          , IterKac;
      'ER'        , WakeUp_ER    , Std_ER,2*n          , 'NA';
      'Uniform'   , WakeUp_Uni   , Std_Uni,2*n          , 'NA';
      'Metropolis', WakeUp_Met   , Std_Met,2*n          , 'NA'};
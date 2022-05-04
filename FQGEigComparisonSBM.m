%% Comparison of ERGossip with FQG on SBM
FG_Data_SBM={"n", "Lambda ER","Eigen_MET","Lambda FG","CPUTime ER","CPUTime_FG","Status"};
% Parameters 
k=3;      % Number of clusters
p=0.9;    % Inner cluster probability
q=0.01;   % Outer cluster probability
for  n=8:2:56
    % Compute the garph Laplacian and adjacency matrices
    [L,D,Adj,Lnorm,Connected]=SBMk(n,k,p,q);
    [n,m]=size(Adj);
    sprintf("Computing eigenvalues of gossiping algorithms on SBM(%g,%g,%g,%g)",n,k,p,q)
    G=graph(Adj);
    % Plot the graph for reference. 
    if n==30
        plot(G,"Layout","auto")
    end
    %% Retrieve the probabilities and the eigenvalues of FQG
    sprintf("Fastest quantum gossiping...")
    time_start=cputime;
    [Eigen_FG,PComm, PWakeUp,WBar, Pall, Status]=FQGopt(Adj);
    CPUtime_FG=cputime-time_start;
    %% Retrieve the probabilities and the eigenvalues of ER based gossiping
    sprintf("ER gossiping...")
    time_start=cputime;
    [PRes_Comm,PRes_WakeUp]=ERProb(L, Adj);
    CPUtime_ER=cputime-time_start;
    Ptemp=diag(PRes_WakeUp)*PRes_Comm;
    Ptemp2=Ptemp+Ptemp';
    DiagER=sum(Ptemp2);
    W=eye(n)-1/2*diag(DiagER)+1/2*Ptemp2;
    EigER=eig(W);
    Eigen_ER=EigER(n-1);
    
    %% Retrieve the probabilities and the eigenvalues of Metropolis based gossiping
    sprintf("Metropolis gossiping...")
    [PMet_Comm,PMet_WakeUp]=METProb(Adj);
    PMET=diag(PMet_WakeUp)*PMet_Comm;
    Ps=PMET+PMET';
    DiagMET=diag(sum(Ps));
    WMET=eye(n)-1/2*DiagMET+1/2*Ps;
    eigens=eig(WMET);
    Eigen_MET=eigens(n-1);
    FG_Data_SBM=[FG_Data_SBM;{n,-log(1-Eigen_ER),-log(1-Eigen_MET),-log(1-Eigen_FG),CPUtime_ER,CPUtime_FG},Status];
    
end

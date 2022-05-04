% Comparison of ERGossip with FQG on barbell graphs
FG_Data_Barbell={"n", "Lambda ER","lambda Met","Lambda FG","CPUTime ER","CPUTime_FG","Status"};
for n=10:2:50
    %% Compute the graph Laplacian and adjacency matrices
    [L,Adj,degree]=graph_barbell(n);
    [n,m]=size(Adj);
    sprintf("Computing eigenvalues for n=%g on barbell graph",n)
    %% Retrieve the probabilities and the eigenvalues of FQG
    sprintf("Fastestquantum gossiping...")
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
    FG_Data_Barbell=[FG_Data_Barbell;{n,-log(1-Eigen_ER), -log(1-Eigen_MET),-log(1-Eigen_FG),CPUtime_ER,CPUtime_FG}, Status];
    
end
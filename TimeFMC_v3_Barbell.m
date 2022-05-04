%% Subgradient method to find the communication weights of FMMC. 
% Computes the average computation to find the communication
% probabilities of FMMC. Saves them to file named ComparisonResults_FMC.
AverageTime=[];
AccuracyData={};
IterationData={};
global vect
global eps
index=0;
for n=vect
    [EDGE,L,degree]=graph_barbell(n);
    Adj=diag(diag(L))-L;
    m=0.5*sum(sum(Adj));
    edges=[];
    for i=1:2*n
        for j=i:2*n
            if Adj(i,j)~=0
                edges=[edges;i,j];
            end
        end
    end
    al=[];
    for i=1:length(edges)
        al(edges(i,1),i)=1;
        al(edges(i,2),i)=-1;
    end
    % Comparison of Computation Times for FMMC
    cvx_begin sdp
    variable s
    variable P(2*n,2*n)
    minimize(s)
    subject to
    Ps=(P+P');
    Ds=diag(sum(P+P'));
    W = eye(2*n)-1/(4*n)*Ds+1/(4*n)*Ps;
    for i=1:2*n
        for j=1:2*n
            if EDGE(i,j)==1
                P(i,j) >= 0;
            else
                P(i,j)==0;
            end
        end
    end
    W-1/(2*n)*ones(2*n,2*n) <= s*eye(2*n,2*n)
    P*ones(2*n,1)==ones(2*n,1)
    cvx_end
    P_fmc=full(P)
    W_fmc=full(W)
    %% Decentralized Subgradient Method
    communication=0;
    edges=[];
    for i=1:2*n
        for j=i:2*n
            if Adj(i,j)~=0
                edges=[edges;i,j];
            end
        end
    end
    for i=(m+1):2*m
        edges(i,1)=edges(i-m,2);
        edges(i,2)=edges(i-m,1);
    end
    node_probabilities=ones(1,2*n)/(2*n);
    p=[];
    % Start with uniform distribution probabilities Pij=1/di;
    for i=1:2*n
        Mi=edges(:,1)==i; % This collects the locations of Ni in edges
    end
    Psub=inv((diag(diag(L))))*Adj;
    Wsub=eye(2*n)-1/(4*n)*(diag((Psub+Psub')*ones(2*n,1)))+1/(4*n)*(Psub+Psub');
    error= norm(P_fmc-Psub,"fro")/norm(P_fmc,'fro');
    AccuracySub=[error];
    CalctimeSub=zeros(1,2*n);
    convergence=0;
    calc_time=0;
    iter=1;
    Iteration_gen=0;
    % Error on n=20 starts at 0.0584 but so I changed the error threshold
    % to 0.05 only at that because for the rest the algorithm does not
    % converge for that threshold. 
    all_count=0;
    gk_count=0;
    
    while error>eps && iter<10^4    % Iteration of Subgradient step k
        vk=rand(2*n,1);
        vk=vk/norm(vk);
        [Ueig,eigs]=eig(Wsub);
        eigs=diag(eigs);
        evalue=sort(eigs,'descend');
        findvec=eigs==evalue(2);
        evector=Ueig(:,findvec==1);
        error=norm(evector-vk)/(2*n);
        while error>eps  % Second eigenvector    
            vk=Wsub*vk;
            Iteration_gen=Iteration_gen+sum(degree);
            % Find sumvk for each node via consensus
            consns=vk;
            % Orthogonalization
            err_cons1=norm(consns-mean(vk)*ones(2*n,1))/sqrt(2*n);
            while err_cons1> eps
                Iteration_gen=Iteration_gen+1;
                node_cons=randsample(1:2*n,1,true,1/(2*n)*ones(2*n,1));
                Mi=edges(:,1)==node_cons;
                pi=Psub(node_cons,edges(Mi,2));
                edge_cons=randsample(edges(Mi,2),1,true,pi);
                consns(node_cons)=0.5*(consns(node_cons)+consns(edge_cons));
                consns(edge_cons)=consns(node_cons);
                err_cons1=norm(consns-mean(vk)*ones(2*n,1))/sqrt(2*n);
            end
            sumvk=mean(consns);
            vk=vk-sumvk*ones(2*n,1);
            consns2=vk.^2;
            err_cons2=norm(consns2-mean(vk.^2)*ones(2*n,1))/sqrt(2*n);
            while err_cons2>eps
                Iteration_gen=Iteration_gen+1;
                node_cons2=randsample(1:2*n,1,true,1/(2*n)*ones(1,2*n));
                Mi=edges(:,1)==node_cons2;
                pi=Psub(node_cons2,edges(Mi,2));
                edge_cons2=randsample(edges(Mi,2),1,true,pi);
                communication=communication+1;
                consns2(node_cons2)=0.5*(consns2(node_cons2)+consns2(edge_cons2));
                consns2(edge_cons2)=consns2(node_cons2);
                err_cons2=norm(consns2-mean(vk.^2)*ones(2*n,1))/sqrt(2*n);
            end
            sumsqr=sqrt(2*n*mean(consns2));
            vk=vk/sumsqr;
            error=min(norm(evector-vk),norm(evector+vk));
            error=error/sqrt(2*n);
        end
        u=vk;
% % Real Eigenvector
          alpha=24*n^(2)/(iter); % alpha at barbell 
%           [Usub,Vsub]=eig(Wsub);
%           u=Usub(:,2*n-1);
         for i = 1:2*n
            Mi=edges(:,1)==i; % This collects the locations of Ni in edges
            Ni=edges(Mi,2);
            ni=length(Ni);
            for j=1:ni
                Iteration_gen=Iteration_gen+1;
                gk=-1/(4*n)*(u(i)-u(Ni(j)))^2;
                Psub(i,Ni(j))=Psub(i,Ni(j))-alpha*gk;
                if abs(alpha*gk)<0.00001
                   gk_count=gk_count+1; 
                end
                all_count=all_count+1;
            end  %Creating gradient
            % Orthagonal Projection
            pi=Psub(i,Ni)';
            sumpi=sum(max(0,pi));
            if sumpi<=1
                pi=max(0,pi);
            else
                sumbisect= @(x) (ones(1,ni)*max(0,(pi-x)))-1;
                xroot=bisection(sumbisect,0,max(pi),10^-3);
                pi=max(0,pi-xroot);
            end
            for j=1:ni
                Psub(i,Ni(j))=pi(j,1);
            end
        end
        Dsub=diag((Psub+Psub')*ones(2*n,1));
        Wsub=eye(2*n)-1/(4*n)*Dsub+1/(4*n)*(Psub+Psub');
        eigen_sub=sort(eig(Wsub),'descend');
        iter=iter+1;
        error= norm(P_fmc-Psub,"fro")/(norm(P_fmc,'fro'));
        AccuracySub=[AccuracySub,error];
    end
    IterationData=[IterationData;Iteration_gen];
    AccuracyData=[AccuracyData;AccuracySub];
    index=index+1;
    if index==1
        Psub1=Psub;
    elseif index==2
        Psub2=Psub;
    elseif index==3 
        Psub3=Psub;
    elseif index==4
        Psub4=Psub;
    elseif index==5
        Psub5=Psub;
    end
end



%AverageTime_FMC=AverageTime/60;
AccuracyData_FMC=AccuracyData;
IterationData_FMC=IterationData;
% save('ComparisonResults_FMC','AverageTime_FMC')
save('ComparisonResults_FMC','AccuracyData_FMC');
save('ComparisonResults_FMC','IterationData_FMC','-append');
save('ComparisonResults_FMC','Psub1','-append');
save('ComparisonResults_FMC','Psub2','-append');
save('ComparisonResults_FMC','Psub3','-append');
save('ComparisonResults_FMC','Psub4','-append');
save('ComparisonResults_FMC','Psub5','-append');




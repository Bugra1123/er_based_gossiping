%% Decentralized normalized Kaczmarz to compute ER weights
% Computes the average computation to find the ER based communication
% probabilities together with the effective resistances. Saves them to file
% named ComparisonResults_Kac.
AverageTime_Kac=[];
IterationData_Kac={};
AccuracyData_Kac={};
Data_Distance_Kac=zeros(4,500);
global vect
% vect=[5,10,15,20,25];
global eps
% eps=0.01;
index=0;
for n=vect
    %% Graph Create
    %M=floor(0.80*2*n);
    %[L,degree]=smallworld_graph(2*n,M);
    [EDGE,L,degree]=graph_barbell(n);
    Adj=diag(diag(L))-L;
    m=0.5*sum(sum(Adj));
    Accrcy=0.01;
    iter_max=10^4;
    % Collect edges
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
    %%Real Values
    Linv=pinv(L);
    R=ones(2*n,1)*diag(Linv)'+diag(Linv)*ones(1,2*n)-2*Linv;
    R=R.*EDGE;
    Dr=R*ones(2*n,1);
    P_opt=inv(diag(Dr))*R;
    %% Kacmarz Method
    iter_sim=0;
    % We obtain the effective resistance rows in distributed matter
    I=eye(2*n);
    I=I-1/(2*n);
    w=degree.^2+degree;
    D=diag(1./sqrt(w));
    Lb=D*L;
    N=2*n;
    pb=ones(N,1)/N;
    cdf_pb = cumsum(pb);
    Ib=D*I;
    X=zeros(2*n);
    iter_kac=0;
    error=1;
    Accuracy=[];
    %
    while error>eps && iter_kac<6*10^8 
        for i=1:N
            iter_kac=iter_kac+degree(i);
            for j=1:N
                X(:,j)=X(:,j)+(Ib(i,j)-Lb(i,:)*X(:,j))*Lb(i,:)';
            end
        end
        R=ones(2*n,1)*diag(X)'+diag(X)*ones(1,2*n)-2*X;
        R=R.*EDGE;
        Dr=R*ones(2*n,1);
        P=pinv(diag(Dr))*R;
        error=norm(P-P_opt,'fro')/(norm(P_opt,'fro'));
        Accuracy=[Accuracy, error];
    end  
    
    
    AccuracyData_Kac=[AccuracyData_Kac;Accuracy];
    %     AverageTime_Kac=[AverageTime_Kac,calc_time_kac/(2*n)];
    IterationData_Kac=[IterationData_Kac;iter_kac];
    index=index+1;
    if index==1
        Pi_res1=P;
        R1=Dr;
    elseif index==2
        Pi_res2=P;
        R2=Dr;
    elseif index==3 
        Pi_res3=P;
        R3=Dr;
    end
end
AverageTime_Kac=AverageTime_Kac/60;
save('ComparisonResults_Kac','AccuracyData_Kac');
save('ComparisonResults_Kac','AverageTime_Kac','-append');
save('ComparisonResults_Kac','IterationData_Kac','-append');
save('ComparisonResults_Kac','Pi_res1','-append')
save('ComparisonResults_Kac','R1','-append')
% % save('ComparisonResults_Kac','Pi_res2','-append')
% % save('ComparisonResults_Kac','R2','-append')

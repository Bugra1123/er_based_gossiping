Data=["n","c","LowerBound","SpecGap","UpperBound","LowerRate","UpperRate"];
for n=5:5:5
    for c=5:1:10
        C=ones(n)-eye(n);
        A=C;
        for i=1:(c-1)
            A=[A,zeros(i*n,n);
                zeros(n,i*n),C];
            A(i*n,i*n+1)=1;
            A(i*n+1,i*n)=1;
        end
        Deg=sum(A');
        P=1/(n*c)*inv(diag(Deg))*A;
        D=sum(P'+P);
        WUni=eye(length(D))-0.5*diag(D)+0.5*(P+P');
        eigen=sort(eig(WUni),'descend');
        eigen2=eigen(2);
        % Conductance
        c_s=1/floor(c/2);
        CUni=c_s/(c*n^3);
        % Comparison
        eigen2;
        Lowerbound=1-2*CUni;
        Upperbound=1-CUni^2;
        LowerRate=Lowerbound/eigen2;
        UpperRate=Upperbound/eigen2;
        Data=[Data; 
            c*n, c, -log(Lowerbound), -log(eigen2), -log(Upperbound),LowerRate,UpperRate];
    end
end
%% PLot the graph
G=graph(A)
plot(G)
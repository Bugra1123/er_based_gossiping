function [Eigen,PComm,PWakeUp,WBar, Pall, Status]=FQGopt(Edges)
% Edges: The edges of the required network structure.
[n,m]=size(Edges);
cvx_begin sdp
variables s P(n,n) 
expression W0(n,n)
minimize(s)
subject to
W0=zeros(n,n);
for i=1:n
    for j=1:n
        ei=zeros(n,1);
        ej=zeros(n,1);
        ei(i,1)=1;
        ej(j,1)=1;
        W0=W0+P(i,j)*(eye(n,n)-0.5*(ei-ej)*(ei-ej)');
    end
end
W=W0;
for i=1:n
    for j=1:n
        if Edges(i,j)==1
            P(i,j) >=0;
        else
            P(i,j) ==0;
        end
    end
end
W-1/(n)*ones(n,n) <= s*eye(n,n);
sum(sum(P))==1;
cvx_end
Status=cvx_status
Eigen=s; 
PWakeUp=sum(P');
PWakeUp=full(PWakeUp);
PComm=inv(diag(PWakeUp))*P;
PComm=full(PComm);
Pall=P;
WBar=W;
end
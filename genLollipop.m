function [L,A,d]=genLollipop(n,k)
% n the complete graph with n+1 vertices
% k the number of vertices in the gpath

% First create the complete graph
Acomp=ones(n+1)-eye(n+1);
Apath=zeros(k);

for i=1:(k-1)
    Apath(i,i+1)=1;
    Apath(i+1,i)=1;
end
A=[Acomp, zeros(n+1,k);zeros(k,n+1),Apath];
A(n+1,n+2)=1;
A(n+2,n+1)=1;
d=sum(A');
L=diag(d)-A;
end

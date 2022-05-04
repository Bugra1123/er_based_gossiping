function [L,A,Deg]=c_barbell(n,c)

C=ones(n)-eye(n);
A=C;
for i=1:(c-1)
    A=[A,zeros(i*n,n);
        zeros(n,i*n),C];
    A(i*n,i*n+1)=1;
    A(i*n+1,i*n)=1;
end
Deg=sum(A');
L=diag(Deg)-A;
end
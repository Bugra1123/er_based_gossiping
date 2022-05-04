function [L,A,d]=graph_barbell(N)
% The function is provided by Necdet Serhat Aybat
A = ones(N);
A = A-eye(N);
A = [A, zeros(N); zeros(N), A];
A(N, N+1) = 1;
A(N+1, N) = 1;

d=sum(A,2);
L=diag(d)-A;
end

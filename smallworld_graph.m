function [L,A,d]=smallworld_graph(N,M)
% The small world function is provided by Necdet Serhat Aybat.
E = zeros(N);
for i=1:N-2
    for j=i+2:N
        E(i,j)=1;
    end
end
ind = find(E==1);
E = zeros(N);

for i=1:N-1
    E(i,i+1)=1;
end
E(1,N)=1;

rng(1);
rand_ind = ind(randperm(length(ind),M));
E(rand_ind)=1;

E=E+E';
d=sum(E,2);
A=E;
L=diag(d)-E;
end


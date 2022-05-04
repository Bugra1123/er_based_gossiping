function [L,D,A,Lnorm,Connected]=SBMk(n,k,p,q)
%% Parameters
% n: the number of the edges
% k: The number of clusters in the group (each of them has equal population)
% p: The inner cluster probability
% q: The outer cluster probability
%% Create the adjecency matrix
A=zeros(n);
for i=1:(k-1)
    clustsize= floor(n/k);
    beg_ind=clustsize*(i-1)+1;
    end_ind=beg_ind+clustsize-1;
    A(beg_ind:end_ind,beg_ind:end_ind)=RandCompGraph(clustsize,p);
    for i=beg_ind:end_ind
        for j= end_ind+1:n
            if rand()<=q
                A(i,j)=1;
                A(j,i)=1;
            end
        end
    end
end
clustsize= n-floor(n/k)*(k-1);
beg_ind=n-clustsize+1;
end_ind=n;
A(beg_ind:end_ind,beg_ind:end_ind)=RandCompGraph(clustsize,p);
% for i=beg_ind:end_ind
%     for j= i+1:n
%         if rand()<=q
%             A(i,j)=1;
%             A(j,i)=1;
%         end
%     end
% end

% Compute the degree matrix
D=sum(A);
D=D';
Dmat=diag(D);
% Compute the Laplacian
L=Dmat-A;
G=graph(A);
% plot(G);
% Check if the network is connnected or not
if sum(D==0)>0
    sprintf('Not feasible')
    Connected=0;
    Lnorm=zeros(n);
    return
    % Compute the normalized matrix
else
    Lnorm=inv(sqrt(Dmat))*L*inv(sqrt(Dmat));
    Eigens=sort(eig(Lnorm));
    if Eigens(2)<=0
        sprintf('Not feasible')
        Connected=0;
    else
        Connected=1;
    end
    
end


end

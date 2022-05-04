function [PMet_Com,PMet_WakeUp]=METProb(A)
[n,m]=size(A);
D=sum(A');
for i=1:n
    for j=1:n
        if A(i,j)==1
            PM(i,j)=1/max(D(i),D(j));
        end
    end
end
PM_Diag=sum(PM');
PMet_Com=(eye(n)-diag(PM_Diag))+PM;
% PMet_Com= 0.5*eye(n)+0.5*PMet_Com;
PMet_Com=PMet_Com;
PMet_WakeUp=1/n*ones(1,n);
end
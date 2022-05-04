function [lambda_ER,lambda_met,lambda_uni,lambda_opt]=EigCompLol(n)
k=floor(sqrt(sqrt(2)*n-1))-2;
[L,A,d]=genLollipop(n,k);
%% Effective Resistance weights
[PRes_Comm,PRes_WakeUp]=ERProb(L, A);
PER=diag(PRes_WakeUp)*PRes_Comm;
Ps=PER+PER';
DiagER=diag(sum(Ps));
WER=eye(k+n+1)-1/2*DiagER+1/2*Ps;
eigens=eig(WER);
lambda_ER=eigens(n+k);
%% Metropolis weights
[PMet_Comm,PMet_WakeUp]=METProb(A);
PMET=diag(PMet_WakeUp)*PMet_Comm;
Ps=PMET+PMET';
DiagMET=diag(sum(Ps));
WMET=eye(k+n+1)-1/2*DiagMET+1/2*Ps;
eigens=eig(WMET);
lambda_met=eigens(n+k);
%% Uniform Probs
Puni=1/(n+k+1)*inv(diag(d))*A;
G=graph(0.5*(Puni+Puni'));
plot(G)
Ps=Puni+Puni';
DiagUni=diag(sum(Ps));
WUni=eye(n+k+1)-1/2*DiagUni+1/2*Ps;
eigens=eig(WUni);
lambda_uni=eigens(n+k);

%% Optimal graph
A=6*(n-1)*(n+k+1)+(k+1)*(6*k*sqrt(2*n*(n+1))+(n+1)*(6+k*(k+2))+k^2*(3*n+k+2));
lambda_opt=1-6*(n+k+1)/A;
end


% Compute the upper bound on the diameter of SBM.
n=100; 
innerprob=0.9;
outerprob=0.01;
ncluster=2;
ExpL=zeros(n);
ExpA=zeros(n);
for i=1:1e3
    [L,D,A,Lnorm,Connected]=SBMk(n,ncluster,innerprob,outerprob);
    ExpA=ExpA+A;
    ExpL=ExpL+L;
end
ExpL=1e-3*ExpL;
ExpA=1e-3*ExpA;
eig(ExpA)


    
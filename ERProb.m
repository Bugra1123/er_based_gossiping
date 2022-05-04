function [PRes_Comm,PRes_WakeUp]=ERProb(L, Adj)
[n,m]=size(Adj);
Linv=pinv(L);
R=ones(n,1)*diag(Linv)'+diag(Linv)*ones(1,n)-2*Linv;
R=R.*Adj;
Ri=sum(R');
PRes_Comm=inv(diag(Ri))*R;
PRes_WakeUp=1/(sum(Ri))*Ri;
end
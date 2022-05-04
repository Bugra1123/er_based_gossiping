function [PRes_Com,PRes_WakeUp,KacIter, err_diff]=KaczmarzERProb(L, Adj, eps)
[n,m]=size(Adj);
Linv=pinv(L);
R=ones(n,1)*diag(Linv)'+diag(Linv)*ones(1,n)-2*Linv;
R=R.*Adj;
Ri=sum(R');
P_opt=inv(diag(Ri))*R;
iter_sim=0;
% We obtain the effective resistance using Kaczmarz
I=eye(n);
I=I-1/n;
degree=sum(Adj');
w=degree.^2+degree;
D=diag(1./sqrt(w));
Lb=D*L;
Ib=D*I;
X=zeros(n);
iter_kac=0;
error=1;
err_diff=1;
P=zeros(n);
while error>eps && err_diff>1e-06
% while error>eps
    Pold=P;
    for i=1:n
        err_diff=error;
        iter_kac=iter_kac+degree(i);
        for j=1:n
            X(:,j)=X(:,j)+(Ib(i,j)-Lb(i,:)*X(:,j))*Lb(i,:)';
        end
        Rtemp=ones(n,1)*diag(X)'+diag(X)*ones(1,n)-2*X;
        Rtemp=Rtemp.*Adj;
        Dr=Rtemp*ones(n,1);
        P=pinv(diag(Dr))*Rtemp;
        error=norm(P-P_opt,'fro')/(norm(P_opt,'fro'));
    end
    err_diff=norm(P-Pold,'fro')/(norm(Pold,'fro'));
end
Ri=sum(Rtemp');
PRes_Com=P;
PRes_WakeUp=1/(sum(Ri))*Ri;
KacIter=iter_kac;
end
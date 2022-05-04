%% Comparison of Consensus
% Code first runs the Time FMC_v2 and TimeKaczmarz in order to obtain convergence data anda iterations 
% tha accuracy threshold at eps=0.05;Time_FMC is 0.075 for n=5,10 however due to
% computational issues and scaling we chose 0.05 for n=20;

global vect 
vect=[10,15,20,25];
global eps 
eps=0.05;
%%Setup Data
 run TimeKaczmarz_Smallworld
 run TimeFMC_v3_smallWorld
% 
load('ComparisonResults_Kac')
load('ComparisonResults_FMC')
%  
 
% Dimension of the data
p=1;
index=0;
iter_max=10^4;
% Data initialization
BigDataSub=[];
DataSubDiff=[];
Data_Distance_Sub={};
Data_Distance_Kac={};
Communication_Kac={};
Communication_Sub={};
for n=vect
    index=index+1;
    Communication_Kac{index,1}=[];
    Communication_Sub{index,1}=[];
    %% Graph Create
    n=n/2;
    L=Lap_Kac{index};
    Adj=diag(diag(L))-L;
    m=0.5*sum(sum(Adj));
    Accrcy=0.01;
    eps=0.01;
    edges=[];
    for i=1:2*n
        for j=i:2*n
            if Adj(i,j)~=0
                edges=[edges;i,j];
            end
        end
    end
    for i=(m+1):2*m
        edges(i,1)=edges(i-m,2);
        edges(i,2)=edges(i-m,1);
    end
    Al=[];
    for l=1:m
        al=zeros(2*n,1);
        al(edges(l,1))=1;
        al(edges(l,2))=-1;
        Al(:,l)=al;
    end
    %% Kacmarz Method
    iter_kac=1;
    if index==1
        Pi_res=Pi_res1;
        Ri=R1;
    elseif index==2
        Pi_res=Pi_res2;
        Ri=R2;
    elseif index==3
        Pi_res=Pi_res3;
        Ri=R3;
            elseif index==4
        Pi_res=Pi_res4;
        Ri=R4;
    end
    %%Asynchronous Consensus Kacmarz
    wakeprob=1/(sum(Ri))*Ri; %2n-1
    totalwakesres=0;
    sim=10^2;
    for l=1:sim
        xinitial=100*randn(2*n,p);
        xaverage=mean(xinitial);
        xdiffuse=xinitial;
        difference=xdiffuse-xaverage;
        distance=norm(difference)/(sqrt(p)*norm(xaverage));
        wakesres=0;
        Away=[];
        comm_times=0;
        while distance>0.01
            comm_times=comm_times+1;
            nodewake=randsample(1:2*n,1,true,wakeprob);
            wakesres=wakesres+1;
             Mi=Pi_res(nodewake,:)>0;
            pi=Pi_res(nodewake,:);
            communicate=randsample(1:2*n,1,true,pi);
            xdiffuse(nodewake)=0.5*(xdiffuse(nodewake)+xdiffuse(communicate));
            xdiffuse(communicate)=xdiffuse(nodewake);
            difference=xdiffuse-xaverage;
            distance=norm(difference)/(sqrt(p)*norm(xaverage));
            Away=[Away,distance];
        end
        
        if n==vect(1)/2
            Communication_Kac{1,1}=[Communication_Kac{1,1},comm_times];
        elseif n==vect(2)/2
            Communication_Kac{2,1}=[Communication_Kac{2,1},comm_times];
        elseif n==vect(3)/2
            Communication_Kac{3,1}=[Communication_Kac{3,1},comm_times];
        elseif  n==vect(4)/2
            Communication_Kac{4,1}=[Communication_Kac{4,1},comm_times];
        end
        totalwakesres=wakesres+totalwakesres;
    end
    
    %% Asynchronous Gossip for Subdiff
    totalwakes=0;
    if index==1
        P=Psub1;
    elseif index==2
        P=Psub2;
    elseif index==3
        P=Psub3;
    elseif index==4
        P=Psub4;
    end
    for l=1:sim
        comm_times=0;
        xinitial=100*randn(2*n,p);
        xaverage=mean(xinitial);
        xdiffuse=xinitial;
        difference=xdiffuse-xaverage;
        distance=norm(difference)/(sqrt(p)*norm(xaverage));
        wakes=0;
        uniformdist=1/(2*n)*ones(2*n,1);
        away=[];
        while distance>0.01
            comm_times=comm_times+1;
            nodewake=randsample(1:2*n,1,true,uniformdist);
            wakes=wakes+1;
            Mi=P(nodewake,:)>0;
            pi=P(nodewake,Mi);
            communicate=randsample(edges(Mi,2),1,true,pi);
            xdiffuse(nodewake)=0.5*(xdiffuse(nodewake)+xdiffuse(communicate));
            xdiffuse(communicate)=xdiffuse(nodewake);
            difference=xdiffuse-xaverage;
            distance=norm(difference)/(sqrt(p)*norm(xaverage));
            away=[away,distance];          
        end
         if n==vect(1)/2
            Communication_Sub{1,1}=[Communication_Sub{1,1},comm_times];
        elseif n==vect(2)/2
            Communication_Sub{2,1}=[Communication_Sub{2,1},comm_times];
        elseif n==vect(3)/2
            Communication_Sub{3,1}=[Communication_Sub{3,1},comm_times];
        elseif n==vect(4)/2
            Communication_Sub{4,1}=[Communication_Sub{4,1},comm_times];
         end
        totalwakes=totalwakes+wakes;
    end
    DataSubDiff=[DataSubDiff;2*n,m,totalwakes/sim];
    BigDataSub=[BigDataSub;DataSubDiff];
end
%% Table 
% A=zeros(3,2); 
% for i=1:3 
%      A(i,1)=mean(Communication_Kac{i,1});
%      A(i,2)=mean(Communication_Sub{i,1});
% end
% Table=mat2dataset(A);
% Table.Properties.VarNames={'ERA','FMCA'}

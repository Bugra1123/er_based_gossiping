%% Comparison of Consensus
% Code first runs the Time FMC_v3 and TimeNormKaczmarz in order to obtain computation time
% and number of iterations to obtain the communication and wake-up probabilities. 
% The accuracy threshold at eps=0.05;Time_FMC is 0.075 for n=5,10 however due to
% computational issues and scaling we chose 0.05 for n=20.

global vect 
vect=[50];
global eps 
eps=0.01;
%% Calculate the weight computation times.
run TimeNormKaczmarz_Barbell
run TimeFMC_v3_Barbell
% Load the results  
load('ComparisonResults_Kac')
load('ComparisonResults_FMC')


% Problem setup
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
    left=1:n;
    [EDGE,L,degree]=graph_barbell(n);
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
    %% Asynchronous Consensus ER weights computed with Kaczmarz
    I=eye(2*n);
    I=I-1/(2*n);
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
    elseif index==5
        Pi_res=Pi_res5;
        Ri=R5;
    end
    wakeprob=1/(sum(Ri))*Ri;    % Set wake-up probabilities
    totalwakesres=0;
    sim=10^2;                   % Number of simulations.
    for l=1:sim
        % Initialize the local data x[i,:] at each node i
        xinitial=10*randn(2*n,p)+500*ones(2*n,p);
        xinitial((n+1):2*n,:)=10*randn(n,p)-500*ones(n,p);
        % The average vector x_mean=1/N \sum_{i=1,..,N}x[i,:]
        xaverage=mean(xinitial);
        xdiffuse=xinitial;
        % Compute the relative error
        difference=xdiffuse-xaverage;
        distance=norm(difference)/(sqrt(p)*norm(xaverage));
        wakesres=0;
        % Collect the relative error in a vector named Away
        Away=[];
        comm_times=0;
        while distance>0.01
            comm_times=comm_times+1;
            % Each node wakes up
            nodewake=randsample(1:2*n,1,true,wakeprob);
            wakesres=wakesres+1;
            % Communicates with one of its neighbors with probability pi
            pi=Pi_res(nodewake,:);
            communicate=randsample(1:2*n,1,true,pi);
            xdiffuse(nodewake)=0.5*(xdiffuse(nodewake)+xdiffuse(communicate));
            xdiffuse(communicate)=xdiffuse(nodewake);
            % Compute the relative error
            difference=xdiffuse-xaverage;
            distance=norm(difference)/(sqrt(p)*norm(xaverage));
            % Collect the relative error at Away
            Away=[Away,distance];
        end
        % Save data
        if n==vect(1)
            Communication_Kac{1,1}=[Communication_Kac{1,1},comm_times];
        elseif n==vect(2)
            Communication_Kac{2,1}=[Communication_Kac{2,1},comm_times];
        elseif n==vect(3)
            Communication_Kac{3,1}=[Communication_Kac{3,1},comm_times];
        end
        totalwakesres=wakesres+totalwakesres;
    end
    
    %% Asynchronous Gossip with FMMC weights computed with subgradient method.
    totalwakes=0;
    if index==1
        P=Psub1;
    elseif index==2
        P=Psub2;
    elseif index==3
        P=Psub3;
    end
    for l=1:sim
        comm_times=0;
        xinitial=10*randn(2*n,p)+500*ones(2*n,p);
        xinitial((n+1):2*n,:)=10*randn(n,p)-500*ones(n,p);
        xaverage=mean(xinitial);
        xdiffuse=xinitial;
        difference=xdiffuse-xaverage;
        distance=norm(difference)/(sqrt(p)*norm(xaverage));
        wakes=0;
        uniformdist=1/(2*n)*ones(2*n,1);
        away=[];
        while distance>0.01
            comm_times=comm_times+1;
            % Each node wakes up
            nodewake=randsample(1:2*n,1,true,uniformdist);
            wakes=wakes+1;
            Mi=Adj(nodewake,:)>0;
            pi=P(nodewake,Mi);
            % Communicates with one of its neighbors with probability pi
            communicate=randsample(edges(Mi,2),1,true,pi);
            xdiffuse(nodewake)=0.5*(xdiffuse(nodewake)+xdiffuse(communicate));
            xdiffuse(communicate)=xdiffuse(nodewake);
            % Compute the relative error
            difference=xdiffuse-xaverage;
            distance=norm(difference)/(sqrt(p)*norm(xaverage));
            % Collect the relative error at away
            away=[away,distance];          
        end
        % Store data 
        if n==vect(1)
            Communication_Sub{1,1}=[Communication_Sub{1,1},comm_times];
        elseif n==vect(2)
            Communication_Sub{2,1}=[Communication_Sub{2,1},comm_times];
        elseif n==vect(3)
            Communication_Sub{3,1}=[Communication_Sub{3,1},comm_times];
        end
        totalwakes=totalwakes+wakes;
    end
    DataSubDiff=[DataSubDiff;2*n,m,totalwakes/sim];
    BigDataSub=[BigDataSub;DataSubDiff];
end
%% Table 
A=zeros(3,2); 
for i=1:3 
     A(i,1)=mean(Communication_Kac{i,1});
     A(i,2)=mean(Communication_Sub{i,1});
end
Table=mat2dataset(A);
Table.Properties.VarNames={'ERA','FMCA'}

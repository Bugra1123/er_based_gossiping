function [WakeUps,Std]=Gossip(Sim,P_WakeUp,P_Comm,Dimension,Acc)
    % Asynchronous randomized gossiping algorithm
    % Parameters: 
    %   Sim: The number of simulations.
    %   P_Wakeup:   Wake-up probabilities for the agents.
    %   P_Comm:     Communication probabilities for each agent.
    %   Dimension:  The dimnesion of the problem, i.e. the dimension of the
    %               local iterates x[i,:]
    %   Acc:        The accuracy for the relative Frobenius norm of
    %               distance between X=x[i,:] and average vector x_mean=1/N
    %               \sum_{i=1,..N}x[i.:]
[n,m]=size(P_Comm);
Nodes=1:n;
WakeUps=0;
WakeUp_vec=[];
for sample=1:Sim
    xinitial=1000*rand(n,Dimension); 
    xmean=mean(xinitial);
    xdiffuse=xinitial;
    diff=xdiffuse-xmean;
%     error=norm(diff,'fro')/(sqrt(Dimension)*norm(xmean));
    error=norm(diff,'fro')^2/(Dimension*n*norm(xmean)^2);
    WakeUps_Sample=0;
    while error>Acc
        % Node wakes up
        WakeUps_Sample=WakeUps_Sample+1;
        NodeWake=randsample(1:n,1,true,P_WakeUp);
        % Chooses a neighbor to communicate
        Neighbor=Nodes(P_Comm(NodeWake,:)>0);
        Pi=P_Comm(NodeWake,P_Comm(NodeWake,:)>0);
        Node_Comm=randsample(Neighbor,1,true,Pi);
        % NodeWake gossips with neighbor. 
        xdiffuse(NodeWake,:)=0.5*(xdiffuse(NodeWake,:)+xdiffuse(Node_Comm,:));
        xdiffuse(Node_Comm,:)=xdiffuse(NodeWake,:);
        diff=xdiffuse-xmean;
        error=(norm(diff,'fro'))^2/(Dimension*n*norm(xmean))^2;
    end
%     WakeUps=WakeUps+WakeUps_Sample;
    WakeUp_vec=[WakeUp_vec,WakeUps_Sample];
end
% WakeUps=WakeUps/Sim;
WakeUps=mean(WakeUp_vec);
Std=std(WakeUp_vec);
end

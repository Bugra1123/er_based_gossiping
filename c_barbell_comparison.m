close all
Data=["n","c","SpecGapUni","SpecGapER", "Ratio"];
% DataUni=["n", "c", "SpecGapUni","LowerUni", "UpperUni"];
DataUni=["n", "c", "AvTimeUni","LowerUni", "UpperUni"];
% DataUni=[];
% DataER=["n", "c", "SpecGapER","LowerER", "UpperER"];
DataER=["n", "c", "AvTimeER","LowerER", "UpperER"];
% DataER=[];
for c=10
    for n=[50, 100,500,1000]
        % Adjecency matrix of c_Barbell graph
        [L,A,Degree]=c_barbell(n,c);
        % PLot the graph
        if n==10 && c==5
            figure(1)
            G=graph(A)
            plot(G)
        end
        %% W_uniform
        P=1/(n*c)*inv(diag(Degree))*A;
        D=sum(P'+P);
        WUni=eye(length(D))-0.5*diag(D)+0.5*(P+P');
%         eigen=sort(eig(WUni),'descend');
        opts.tol=1e-50;
        eigen=sort(eigs(WUni,2,'LM',opts),'descend');
        eigen2=eigen(2);
        %%% Conductance for uniform
        c_s=1/floor(c/2);
        CUni=c_s/(c*n^3);
        %%% Comparison
        Lowerbound=1-2*CUni;
        Upperbound=1-CUni^2;
%         DataUni=[DataUni; 
%             n, c, -1/log(eigen2), -1/log(Lowerbound), -1/log(Upperbound)];
        DataUni=[DataUni; 
            n, c, -log(1-eigen2), log(2*CUni)/log(1-eigen2), 2*log(CUni)/log(1-eigen2)];
        
%         DataUni=[DataUni; 
%             n, c,eigen2, Lowerbound, Upperbound];
        %% W_er        
        [PRes_Comm,PRes_WakeUp]=ERProb(L, A);
        PER=diag(PRes_WakeUp)*PRes_Comm;
        Ps=PER+PER';
        DiagER=diag(sum(Ps));
        Wer=eye(length(D))-1/2*DiagER+1/2*Ps;
%         eigen=sort(eig(Wer),'descend');
        eigen=sort(eigs(Wer,2,'LM',opts),'descend');
        eigen2er=eigen(2);
       
        %%% Conductance for ER
        c_s=1/floor(c/2);
        CER=c_s/(2*n*(c*n-1));
        %%% Comparison
        Lowerbounder=1-2*CER;
        Upperbounder=1-CER^2;
%         DataER=[DataER; 
%             n, c, -1/log(eigen2er), -1/log(Lowerbounder), -1/log(Upperbounder)];
        DataER=[DataER; 
             n, c, -log(1-eigen2er), log(2*CER)/log(1-eigen2er), 2*log(CER)/log(1-eigen2er)];
%         DataER=[DataER; 
%              n, c, -1/log(eigen2er), log(eigen2er)/log(Lowerbounder), log(eigen2er)/log(Upperbounder)];
        % Comparison Data 
        Data=[Data; 
            n, c, -log(1-eigen2),-log(1-eigen2er), log(eigen2er)/log(eigen2)];
    end
end
% %% Plots 
% close all
% figure(2)
% i=1;
% l=length(DataUni);
% DataUniNum=str2double(DataUni(2:l,:))
% DataERNum=str2double(DataER(2:l,:))
% for c=10
%     idx=DataUniNum(:,2)==10;
%     x=DataUniNum(idx,1)-0.1;
%     y=DataUniNum(idx,3);
%     yneg=DataUniNum(idx,4);
%     ypos=DataUniNum(idx,5);
% %     errorbar(log(x),zeros(length(x),1),log(yneg),log(ypos),'s','markerfacecolor','blue','markersize',0.1,'linewidth',1)
%     plot(x,ypos,'b--','LineWidth',2.5)
%     hold on
%     plot(x,yneg,'b-','LineWidth',2.5)
%     hold on
%     x=DataERNum(idx,1)+0.1;
%     y=DataERNum(idx,3);
%     yneg=DataERNum(idx,4);
%     ypos=DataERNum(idx,5);
% %     errorbar(log(x),zeros(length(x),1),log(yneg),log(ypos), 'o','markerfacecolor','red','markersize',0.1, 'linewidth',1)
%     plot(x,ypos,'r--','LineWidth',2.5)
%     hold on 
%     plot(x,yneg,'r-','LineWidth',2.5)
%     hold on 
%     i=i+1;
% end
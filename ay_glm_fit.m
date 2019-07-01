% load data
clear
close all
clc
% Data = generateCovariates('jn029-11');
datmat = 'pr002-4';
[Data, D, W] = generateCovariates(datmat);
% important indexces
len = 350;
ind= [1 46  84  124 164 189 228 267 333 len]; % Ali original samples
% ind= [1 43  83  129 164 189 214 254 300 len]; % Corrected samples
% Xt
Xt = Data(1).timevec;
% No of trials and cells
nt = length(Data);
nc = size(Data(1).Ca1_spikes,1);
ns = size(Data(1).Ca1_spikes,2);
% I assume the result are correct
Bs = [];
Is = [];
X = [];
Y = [];
i = 1;    
TrialInfo= [];
for j=1:nt
        % Observed spikees
        Y = [Y;Data(j).Ca1_spikes(i,:)'];
        % Desing matrix
        Temp = zeros(ns,22);
        % First element
        Temp(ind(1):ind(2)-1,1)=1;
        Temp(ind(2):ind(3)-1,2)=1;
        % CW and CCW
        dCW  = find(diff(Data(j).CW)==1);
        dCCW = find(diff(Data(j).CCW)==1);
        % CW-CW
        if  length(dCW)==2 
            % The first half is the same for CW
            Temp(ind(3):ind(4)-1,3)= 1;
            Temp(ind(4):ind(5)-1,4)= 1;
            Temp(ind(5):ind(6)-1,5)= 1;
            Temp(ind(6):ind(7)-1,6)= 1;
            % The second half will be unique for CW-CW
            Temp(ind(7):ind(8)-1,11) = 1;
            Temp(ind(8):ind(9)-1,12) = 1;
            Temp(ind(9):ind(10),13)= 1;
            % Type of trial
            TrialInfo= [TrialInfo;1];
        end
        % CW-CCW
        if length(dCW)==1 &&  dCW(1) < dCCW(1) 
            % The first half is the same for CW
            Temp(ind(3):ind(4)-1,3)= 1;
            Temp(ind(4):ind(5)-1,4)= 1;
            Temp(ind(5):ind(6)-1,5)= 1;
            Temp(ind(6):ind(7)-1,6)= 1;
            % The second half will be different for CW-CCW
            Temp(ind(7):ind(8)-1,14)= 1;
            Temp(ind(8):ind(9)-1,15)= 1;
            Temp(ind(9):ind(10),16)= 1;
            % Type of trial
            TrialInfo= [TrialInfo;2];
        end
        % CCW-CW
        if length(dCCW)==1 && dCCW(1) < dCW(1)
            % The first half is the same for CCW
            Temp(ind(3):ind(4)-1,7)= 1;
            Temp(ind(4):ind(5)-1,8)= 1;
            Temp(ind(5):ind(6)-1,9)= 1;
            Temp(ind(6):ind(7)-1,10)= 1;
            % The second half is unique for CCW-CW
            Temp(ind(7):ind(8)-1,17)= 1;
            Temp(ind(8):ind(9)-1,18)= 1;
            Temp(ind(9):ind(10),19)= 1;
            % Type of trial 
            TrialInfo= [TrialInfo;3];
        end
        % CCW-CCW
        if length(dCCW)==2
            % The first half is the same for CCW
            Temp(ind(3):ind(4)-1,7)= 1;
            Temp(ind(4):ind(5)-1,8)= 1;
            Temp(ind(5):ind(6)-1,9)= 1;
            Temp(ind(6):ind(7)-1,10)= 1;
            % The second half is CCW and CCW
            Temp(ind(7):ind(8)-1,20)= 1;
            Temp(ind(8):ind(9)-1,21)= 1;
            Temp(ind(9):ind(10),22)= 1;
            % Trial type
            TrialInfo= [TrialInfo;4];
        end
        X = [X; Temp];
end
% By this point, we have built X and Y

Xz = [];
for s=1:4
        Temp = zeros(ns,22);
        Temp(ind(1):ind(2)-1,1) = 1;
        Temp(ind(2):ind(3)-1,2) = 1;
        if  s==1 
            % the first half is the same
            Temp(ind(3):ind(4)-1,3)  = 1;
            Temp(ind(4):ind(5)-1,4)  = 1;
            Temp(ind(5):ind(6)-1,5)  = 1;
            Temp(ind(6):ind(7)-1,6)  = 1;
            
            Temp(ind(7):ind(8)-1,11) = 1;
            Temp(ind(8):ind(9)-1,12) = 1;
            Temp(ind(9):ind(10),13)= 1;
        end
        if s==2 
            % the first half is the same
            Temp(ind(3):ind(4)-1,3)  = 1;
            Temp(ind(4):ind(5)-1,4)  = 1;
            Temp(ind(5):ind(6)-1,5)  = 1;
            Temp(ind(6):ind(7)-1,6)  = 1;
            
            Temp(ind(7):ind(8)-1,14) = 1;
            Temp(ind(8):ind(9)-1,15) = 1;
            Temp(ind(9):ind(10),16)= 1;
            
        end
        if s==3
            % the first half is the same
            Temp(ind(3):ind(4)-1,7)  = 1;
            Temp(ind(4):ind(5)-1,8)  = 1;
            Temp(ind(5):ind(6)-1,9)  = 1;
            Temp(ind(6):ind(7)-1,10) = 1;
            
            
            Temp(ind(7):ind(8)-1,17) = 1;
            Temp(ind(8):ind(9)-1,18) = 1;
            Temp(ind(9):ind(10),19)= 1;
        end
        if s==4
            % the first half is the same
            Temp(ind(3):ind(4)-1,7)  = 1;
            Temp(ind(4):ind(5)-1,8)  = 1;
            Temp(ind(5):ind(6)-1,9)  = 1;
            Temp(ind(6):ind(7)-1,10) = 1;
            
            Temp(ind(7):ind(8)-1,20) = 1;
            Temp(ind(8):ind(9)-1,21) = 1;
            Temp(ind(9):ind(10),22)= 1;
        end
        Xz = [Xz; Temp];
end
%%
YS = zeros(nc,1600);
NT = [length(find(TrialInfo==1)),length(find(TrialInfo==2)),length(find(TrialInfo==3)),length(find(TrialInfo==4))];
iTable = [];
for i=1:nc
    Y = [];
    disp(i)
    for j=1:nt
        Y = [Y;Data(j).Ca1_spikes(i,:)'];
        if TrialInfo(j)==1
            YS(i,1:len)=YS(i,1:len)+Data(j).Ca1_spikes(i,:);
        end
        if TrialInfo(j)==2
            YS(i,len+(1:len))=YS(i,len+(1:len))+Data(j).Ca1_spikes(i,:);
        end
        if TrialInfo(j)==3
             YS(i,2*len+(1:len))=YS(i,2*len+(1:len))+Data(j).Ca1_spikes(i,:);
        end
        if TrialInfo(j)==4
             YS(i,3*len+(1:len))=YS(i,3*len+(1:len))+Data(j).Ca1_spikes(i,:);
        end
    end
    
    iTable = [iTable; i sum(Y)];
    if sum(Y)> nt*1.25
        [B,~,stat] = glmfit(X,Y,'poisson','link','log','constant','off');
        Bs = [Bs; B'];
        Is = [Is; i];
        [Yt,Yt_low,Yt_high] = glmval(B,Xz,'log',stat,'constant','off');
        figure('units','normalized','outerposition',[0 0 1 1])
        subplot(2,2,1)
        title(num2str(i))
%         e = errorbar(1:length(B),B,2*stat.se);
%         hold on
%         e.Marker     = '*';
%         e.MarkerSize = 10;
%         e.Color   = 'red';
%         e.CapSize = 15;
%         e.LineWidth = 0.2;
        
        e = errorbar(3:6,B(3:6),2*stat.se(3:6));
        hold on
        e.Marker     = '*';
        e.MarkerSize = 10;
        e.Color   = 'cyan';
        e.CapSize = 15;
        e.LineWidth = 2.5;
        
        e = errorbar(3:6,B(7:10),2*stat.se(7:10));
        e.Marker     = 'o';
        e.MarkerSize = 10;
        e.Color      = 'magenta';
        e.CapSize    = 15;
        e.LineWidth  = 2.5;
        
        e = errorbar(8:2:12,B(11:13),2*stat.se(11:13));
        hold on
        e.Marker     = '*';
        e.MarkerSize = 10;
        e.Color      = 'red';
        e.CapSize    = 15;
        e.LineWidth  = 2.5;
        
        e = errorbar(0.1+(8:2:12),B(14:16),2*stat.se(14:16));
        hold on
        e.Marker     = '^';
        e.MarkerSize = 10;
        e.Color      = 'green';
        e.CapSize    = 15;
        e.LineWidth  = 2.5;
        
        e = errorbar(0.2+(8:2:12),B(17:19),2*stat.se(17:19));
        hold on
        e.Marker     = 'o';
        e.MarkerSize = 10;
        e.Color      = 'blue';
        e.CapSize    = 15;
        e.LineWidth  = 2.5;
        
        e = errorbar(0.3+(8:2:12),B(20:22),2*stat.se(20:22));
        hold on
        e.Marker     = '+';
        e.MarkerSize = 10;
        e.Color      = 'black';
        e.CapSize    = 15;
        e.LineWidth  = 2.5;
        
        axis tight
        title(['Cell ' num2str(i)])
        xlabel('Index')
        ylabel('Coef.')
        legend('cw','ccw','cw-cw','cw-ccw','ccw-cw','ccw-ccw','location','best')
        hold off
        xticks([3 4 5 6 8 10 12])
        xticklabels({'2.5-3.7','3.7-5','5-5.7','5.7-6.9','6.9-8','8-10.1','10.9-12.1'})
        y_min = max(min( B(3:end)-2*stat.se(3:end)),-11);
        y_max = min(max( B(3:end) +2*stat.se(3:end)),-2);
        ylim([y_min  y_max])
    
        subplot(4,2,2)
        plot(((1:len)/33),YS(i,1:len)/NT(1),'r','LineWidth',2);
        xlim([0 12])
        ylim([0 0.1])
        ylabel('norm-count')
        legend('cw-cw')
        
        
        subplot(4,2,4)
        plot(((1:len)/33),YS(i,len+(1:len))/NT(2),'g','LineWidth',2);
        xlim([0 12])
        ylim([0 0.1])
        ylabel('norm- ount')
        legend('cw-ccw')

        subplot(4,2,6)
        plot(((1:len)/33),YS(i,2*len+(1:len))/NT(3),'b','LineWidth',2);
        xlim([0 12])
        ylim([0 0.1])
        ylabel('norm-count')
        legend('ccw-cw')

        subplot(4,2,8)
        plot(((1:len)/33),YS(i,3*len+(1:len))/NT(4),'k','LineWidth',2);
        xlim([0 12])
        ylim([0 0.1])
        ylabel('norm-count')
        legend('ccw-ccw')
        xlabel('Time (sec)')
        
        subplot(2,2,3)
%         patch([1:len len:-1:1]',[Yt_low(1:len); Yt_high(len:-1:1)],'r')
%         hold on
%         alpha 0.5
%         patch([1:len len:-1:1]',[Yt_low(len+(1:len));  Yt_high(len+(len:-1:1))],'g')
%         patch([1:len len:-1:1]',[Yt_low(2*len+(1:len));  Yt_high(2*len+(len:-1:1))],'b')
%         patch([1:len len:-1:1]',[Yt_low(3*len+(1:len)); Yt_high(3*len+(len:-1:1))],'k')
        
         plot(((1:len))/33,Yt(1:len),'r','LineWidth',2);
         hold on
         plot(((1:len))/33,Yt(len+(1:len)),'g','LineWidth',2);
         plot(((1:len))/33,Yt(2*len+(1:len)),'b','LineWidth',2);
         plot(((1:len))/33,Yt(3*len+(1:len)),'k','LineWidth',2);
        
        
        hold off
        xlabel('Time (sec)')
        ylabel('Rate')
        legend('cw-cw','cw-ccw','ccw-cw','ccw-ccw','location','best')
        axis tight
        y_min = max(min(Yt_low),1e-6);
        y_max = min(max(Yt_high),1e2);
        ylim([y_min  Inf])
    
        saveas(gcf,['./pr002 Images/' datmat '/Cell ' num2str(i)],'tiff')
        close all
        clc
%         pause()
    end
    
end


subplot(2,1,1)
imagesc(Bs)
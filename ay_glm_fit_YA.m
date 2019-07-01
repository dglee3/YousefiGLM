clear
close all
clc

[Data, D, W, samps] = generateCovariates('jn029-11');
% No of trials and cells
nt = length(Data);
nc = size(Data(1).Ca1_spikes,1);
ns = size(Data(1).Ca1_spikes,2);
% I assume the result are correct
Bs = [];
Is = [];
cell_name = [];
X = [];
Y = [];
i = 1;
TrialInfo= [];

for j=1:nt
    % Observed spikees
    Y = [Y;Data(j).Ca1_spikes(i,:)'];
    [~,nz_cw,~] = find(Data(j).CW);
    nz_cw_mi = min(nz_cw);
    nz_cw_mx = max(nz_cw);
    [~,nz_ccw,~] = find(Data(j).CCW);
    nz_ccw_mi = min(nz_ccw);
    nz_ccw_mx = max(nz_ccw);
    
    if (isempty(nz_cw))
        num1 = find(diff(nz_ccw)>1);
        mat = [nz_ccw_mi nz_ccw_mx nz_ccw(num1) nz_ccw(num1+1)];
    elseif ( isempty(nz_ccw))
        num1 = find(diff(nz_cw)>1);
        mat = [nz_cw_mi nz_cw_mx nz_cw(num1) nz_cw(num1+1)];
    else
        mat = [nz_cw_mi nz_cw_mx nz_ccw_mi nz_ccw_mx];
    end
    ind = [1 sort(mat) samps];
    % Desing matrix
    Temp = zeros(ns,13);
    % First element
    Temp(ind(1):ind(2)-1,1)=1;
    % CW and CCW
    dCW  = find(diff(Data(j).CW)==1);
    dCCW = find(diff(Data(j).CCW)==1);
    % CW-CW
    if  length(dCW)==2
        % The first half is the same for CW
        Temp(ind(2):ind(3)-1,2)= 1;
        Temp(ind(3):ind(4)-1,3)= 1;
        % The second half will be unique for CW-CW
        Temp(ind(4):ind(5)-1,6) = 1;
        Temp(ind(5):ind(6)-1,7) = 1;
        % Type of trial
        TrialInfo= [TrialInfo;1];
    end
    % CW-CCW
    if length(dCW)==1 &&  dCW(1) < dCCW(1)
        % The first half is the same for CW
        Temp(ind(2):ind(3)-1,2)= 1;
        Temp(ind(3):ind(4)-1,3)= 1;
        % The second half will be different for CW-CCW
        Temp(ind(4):ind(5)-1,8)= 1;
        Temp(ind(5):ind(6)-1,9)= 1;
        % Type of trial
        TrialInfo= [TrialInfo;2];
    end
    % CCW-CW
    if length(dCCW)==1 && dCCW(1) < dCW(1)
        % The first half is the same for CCW
        Temp(ind(2):ind(3)-1,4)= 1;
        Temp(ind(3):ind(4)-1,5)= 1;
        % The second half is unique for CCW-CW
        Temp(ind(4):ind(5)-1,10)= 1;
        Temp(ind(5):ind(6)-1,11)= 1;
        % Type of trial
        TrialInfo= [TrialInfo;3];
    end
    % CCW-CCW
    if length(dCCW)==2
        % The first half is the same for CCW
        Temp(ind(2):ind(3)-1,4)= 1;
        Temp(ind(3):ind(4)-1,5)= 1;
        % The second half is CCW and CCW
        Temp(ind(4):ind(5)-1,12)= 1;
        Temp(ind(5):ind(6)-1,13)= 1;
        % Trial type
        TrialInfo= [TrialInfo;4];
    end
    X = [X; Temp];
end
% By this point, we have built X and Y

YS = zeros(nc,4*samps);
NT = [length(find(TrialInfo==1)),length(find(TrialInfo==2)),length(find(TrialInfo==3)),length(find(TrialInfo==4))];
for i=1:nc
    Y = [];
    for j=1:nt
        Y = [Y;Data(j).Ca1_spikes(i,:)'];
        if TrialInfo(j)==1
            YS(i,1:samps)=YS(i,1:samps)+Data(j).Ca1_spikes(i,:);
        end
        if TrialInfo(j)==2
            YS(i,samps+(1:samps))=YS(i,samps+(1:samps))+Data(j).Ca1_spikes(i,:);
        end
        if TrialInfo(j)==3
            YS(i,2*samps+(1:samps))=YS(i,2*samps+(1:samps))+Data(j).Ca1_spikes(i,:);
        end
        if TrialInfo(j)==4
            YS(i,3*samps+(1:samps))=YS(i,3*samps+(1:samps))+Data(j).Ca1_spikes(i,:);
        end
    end
    
    if i == 86%sum(Y)> 373
        disp(['Cell ' num2str(i)])
        cell_name = [cell_name;i];
        [B,dev,stat] = glmfit(X,Y,'poisson','link','log','constant','off');
        %% Time rescaling theorem:
        Lambda_ = exp(X*B);
        spike_v = zeros(length(Y),1);
        [ind_s,~] = find(Y);
        spike_v(ind_s) = 1;
        figure(1)
        subplot(2,4,5)
        [acf,lags,bounds] = autocorr(spike_v-Lambda_,10);%% Yalda:ACF
        plot(lags,acf,'Linewidth',2)
        hold on
        plot(lags,zeros(11,1)+bounds(1),'r--')
        plot(lags,zeros(11,1)+bounds(2),'r--')
        axis tight
        xlabel('Lag');ylabel('ACF')
        title(['ACF Plot: Cell ' num2str(i)])
        hold off
        figure(1)
        subplot(2,4,6)
        g = zeros(1,length(ind_s));
        h = zeros(1,length(ind_s));
        for k=1:length(ind_s)
            if k>1
                Tao = sum(Lambda_(ind_s(k-1):ind_s(k)));
            else
                Tao = sum(Lambda_(1:ind_s(k)));
            end
            g(k)=1-exp(-Tao);
            h(k)=(k-0.5)/length(ind_s);
            
        end
        plot(h,sort(g),'LineWidth',2)
        xlim([0,1]);
        ylim([0,1]);hold on
        plot([0,1],[0,1],'-.k')
        cb = 1.36/sqrt(length(find(Y)));
        hold on;plot([0 1],[0 1],'k')
        rand_ = linspace(0,1,length(find(Y)));
        plot(rand_,rand_-cb,'--k')
        plot(rand_,rand_+cb,'--k')
        xlim([0,1])
        ylim([0,1])
        ylabel('Model CDF');xlabel('Empiricial CDF')
        title(['KS Plot: Cell ' num2str(i)])
        hold off
        %%
        Bs = [Bs; B'];
        Is = [Is; i];
        figure(1)
        subplot(2,2,1)
        title(num2str(i))
        
        e = errorbar(2:3,B(2:3),2*stat.se(2:3));
        hold on
        e.Marker     = '*';
        e.MarkerSize = 10;
        e.Color   = 'cyan';
        e.CapSize = 15;
        e.LineWidth = 2.5;
        
        e = errorbar(2:3,B(4:5),2*stat.se(4:5));
        e.Marker     = 'o';
        e.MarkerSize = 10;
        e.Color      = 'magenta';
        e.CapSize    = 15;
        e.LineWidth  = 2.5;
        
        e = errorbar(4:2:6,B(6:7),2*stat.se(6:7));
        hold on
        e.Marker     = '*';
        e.MarkerSize = 10;
        e.Color      = 'red';
        e.CapSize    = 15;
        e.LineWidth  = 2.5;
        
        e = errorbar(0.1+(4:2:6),B(8:9),2*stat.se(8:9));
        hold on
        e.Marker     = '^';
        e.MarkerSize = 10;
        e.Color      = 'green';
        e.CapSize    = 15;
        e.LineWidth  = 2.5;
        
        e = errorbar(0.2+(4:2:6),B(10:11),2*stat.se(10:11));
        hold on
        e.Marker     = 'o';
        e.MarkerSize = 10;
        e.Color      = 'blue';
        e.CapSize    = 15;
        e.LineWidth  = 2.5;
        
        e = errorbar(0.3+(4:2:6),B(12:13),2*stat.se(12:13));
        hold on
        e.Marker     = '+';
        e.MarkerSize = 10;
        e.Color      = 'black';
        e.CapSize    = 15;
        e.LineWidth  = 2.5;
        
        axis tight
        title(['Cell ' num2str(i)])
        ylabel('Coef.')
        legend('cw','ccw','cw-cw','cw-ccw','ccw-cw','ccw-ccw','location','best')
        hold off
        set(gca,'xTick',[2 4],'xTickLabel',{'Direction 1 Start','Direction 2 Start'});
        
        y_min = max(min( B(3:end)-2*stat.se(3:end)),-11);
        y_max = min(max( B(3:end) +2*stat.se(3:end)),-2);
        ylim([y_min  y_max])
        
        subplot(4,2,2)
        plot(((1:samps)/33),YS(i,1:samps)/NT(1),'r','LineWidth',2);
        xlim([0 12])
        ylabel('norm-count')
        legend('cw-cw')
        
        
        subplot(4,2,4)
        plot(((1:samps)/33),YS(i,samps+(1:samps))/NT(2),'g','LineWidth',2);
        xlim([0 12])
        ylabel('norm-count')
        legend('cw-ccw')
        
        subplot(4,2,6)
        plot(((1:samps)/33),YS(i,2*samps+(1:samps))/NT(3),'b','LineWidth',2);
        xlim([0 12])
        ylabel('norm-count')
        legend('ccw-cw')
        
        subplot(4,2,8)
        plot(((1:samps)/33),YS(i,3*samps+(1:samps))/NT(4),'k','LineWidth',2);
        xlim([0 12])
        ylabel('norm-count')
        legend('ccw-ccw')
        xlabel('Time (sec)')
        
        
        pause()
    end
    
end


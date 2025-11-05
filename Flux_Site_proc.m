
%% ========================================================================
%  Comparison of DLEM-SIF and DLEM-Ag simulations with flux observations
%  Sites: US-Ne2 (Corn–Soybean Rotation) and US-Ne3 (Continuous Corn)
%  Variables: GPP, ET, and WUE
%  ------------------------------------------------------------------------
%  This script:
%   1. Loads model outputs and flux tower observations
%   2. Calculates ET and WUE
%   3. Filters unrealistic values
%   4. Compares DLEM-SIF vs DLEM-Ag with observations (2018–2023)
%   5. Plots GPP, ET, and WUE seasonal dynamics with evaluation metrics
%  ------------------------------------------------------------------------
%% ========================================================================
close all
clear
clc
%% ----------------------------- Load Data --------------------------------
load .\Data\Neflux.mat

% Read DLEM-SIF simulation outputs
simtable_Ne2 = DLEM_AG_read('.\Data\Ne2_DLEM_SIF.out');
simtable_Ne3 = DLEM_AG_read('.\Data\Ne3_DLEM_SIF.out');

% Read DLEM-Ag (original) simulation outputs
simtable_Ne2_ori = DLEM_AG_read('.\Data\Ne2_DLEM_Ag.out');
simtable_Ne3_ori = DLEM_AG_read('.\Data\Ne3_DLEM_Ag.out');

%% ----------------------------- Align Dates ------------------------------
% Convert simulation year/day to MATLAB datetime and align with observations

% --- Site US-Ne2 ---
simtable_Ne2.date = datetime(simtable_Ne2.year,1,1) + simtable_Ne2.day;
[~, loc] = ismember(simtable_Ne2.date, Neselct_s2.date);
simtable_Ne2.GPPobs = Neselct_s2.GPP(loc);
simtable_Ne2.NEPobs = -Neselct_s2.NEE(loc);
simtable_Ne2.Reobs  = Neselct_s2.Re(loc);
simtable_Ne2.ETobs  = Neselct_s2.fun1_ET(loc);

simtable_Ne2_ori.date = datetime(simtable_Ne2_ori.year,1,1) + simtable_Ne2_ori.day;
[~, loc] = ismember(simtable_Ne2_ori.date, Neselct_s2.date);
simtable_Ne2_ori.GPPobs = Neselct_s2.GPP(loc);
simtable_Ne2_ori.NEPobs = -Neselct_s2.NEE(loc);
simtable_Ne2_ori.Reobs  = Neselct_s2.Re(loc);
simtable_Ne2_ori.ETobs  = Neselct_s2.fun1_ET(loc);

% --- Site US-Ne3 ---
simtable_Ne3.date = datetime(simtable_Ne3.year,1,1) + simtable_Ne3.day;
[~, loc] = ismember(simtable_Ne3.date, Neselct_s2_NE3.date);
simtable_Ne3.GPPobs = Neselct_s2_NE3.GPP(loc);
simtable_Ne3.NEPobs = -Neselct_s2_NE3.NEE(loc);
simtable_Ne3.Reobs  = Neselct_s2_NE3.Re(loc);
simtable_Ne3.ETobs  = Neselct_s2_NE3.fun1_ET(loc);

simtable_Ne3_ori.date = datetime(simtable_Ne3_ori.year,1,1) + simtable_Ne3_ori.day;
[~, loc] = ismember(simtable_Ne3_ori.date, Neselct_s2_NE3.date);
simtable_Ne3_ori.GPPobs = Neselct_s2_NE3.GPP(loc);
simtable_Ne3_ori.NEPobs = -Neselct_s2_NE3.NEE(loc);
simtable_Ne3_ori.Reobs  = Neselct_s2_NE3.Re(loc);
simtable_Ne3_ori.ETobs  = Neselct_s2_NE3.fun1_ET(loc);

%% ----------------------------- Process Data -----------------------------
% Compute ET and WUE, and filter outliers

% Observed WUE
simtable_Ne2.WUEobs      = simtable_Ne2.GPPobs ./ simtable_Ne2.ETobs;
simtable_Ne3.WUEobs      = simtable_Ne3.GPPobs ./ simtable_Ne3.ETobs;
simtable_Ne2_ori.WUEobs  = simtable_Ne2_ori.GPPobs ./ simtable_Ne2_ori.ETobs;
simtable_Ne3_ori.WUEobs  = simtable_Ne3_ori.GPPobs ./ simtable_Ne3_ori.ETobs;

% Simulated ET (transpiration + evaporation)
simtable_Ne2.ET      = simtable_Ne2.g_trans + simtable_Ne2.g_evap;
simtable_Ne3.ET      = simtable_Ne3.g_trans + simtable_Ne3.g_evap;
simtable_Ne2_ori.ET  = simtable_Ne2_ori.g_trans + simtable_Ne2_ori.g_evap;
simtable_Ne3_ori.ET  = simtable_Ne3_ori.g_trans + simtable_Ne3_ori.g_evap;

% Simulated WUE (GPP/ET)
simtable_Ne2.WUE      = simtable_Ne2.g_gpp ./ simtable_Ne2.ET;
simtable_Ne3.WUE      = simtable_Ne3.g_gpp ./ simtable_Ne3.ET;
simtable_Ne2_ori.WUE  = simtable_Ne2_ori.g_gpp ./ simtable_Ne2_ori.ET;
simtable_Ne3_ori.WUE  = simtable_Ne3_ori.g_gpp ./ simtable_Ne3_ori.ET;

% Remove unrealistic WUE values (WUE < 0 or > 9)
datasets = {'simtable_Ne2', 'simtable_Ne3', 'simtable_Ne2_ori', 'simtable_Ne3_ori'};
for i = 1:numel(datasets)
    eval([datasets{i} '.WUE(' datasets{i} '.WUE>9 | ' datasets{i} '.WUE<0) = NaN;']);
    eval([datasets{i} '.WUEobs(' datasets{i} '.WUEobs>9 | ' datasets{i} '.WUEobs<0) = NaN;']);
end

%% ========================================================================
%  Plot GPP (Gross Primary Production)
%% ========================================================================
figure
ha = tight_subplot(2,6,[0.1,0.01],[0.08,0.05],[0.04,0.01]);

% ---------------- US-Ne2 ----------------
for k = 1:6
    axes(ha(k))
    starttime = datetime(2017+k,5,1);
    endtime   = datetime(2017+k,11,1);
    idx = simtable_Ne2.date >= starttime & simtable_Ne2.date <= endtime;

    plot(simtable_Ne2.day(idx)+1, simtable_Ne2.g_gpp(idx),'r-','LineWidth',1.5); hold on
    plot(simtable_Ne2_ori.day(idx)+1, simtable_Ne2_ori.g_gpp(idx),'b-','LineWidth',1.5);
    plot(simtable_Ne2.day(idx)+1, simtable_Ne2.GPPobs(idx),'ko','MarkerFaceColor','k','MarkerSize',3);
    xlim([120,310]); ylim([0,50]);
    xlabel('DOY','FontSize',12)

    if k==1
        ylabel('GPP (g C m^{-2} d^{-1})','FontSize',12)
        le=legend('DLEM-SIF','DLEM-Ag'); le.Position=[0.044 0.7943 0.0677 0.0521]; le.Box='off';
    else
        set(gca,'YTickLabel',[])
    end

    title(sprintf('US-Ne2 (%d)',k+2017),'FontSize',16)
    res1 = evaluate_regression_metrics(simtable_Ne2.g_gpp(idx),simtable_Ne2.GPPobs(idx));
    res2 = evaluate_regression_metrics(simtable_Ne2_ori.g_gpp(idx),simtable_Ne2_ori.GPPobs(idx));
    text(121,43,sprintf('{\\itR}^2 = %2.2f (%2.2f)\nRMSE = %2.2f (%2.2f)\nNRMSE = %2.2f%% (%2.2f%%)',...
        res1.R2,res2.R2,res1.RMSE,res2.RMSE,res1.NRMSE*100,res2.NRMSE*100));

    % Alternate background: corn vs soybean
    if mod(k,2)==0
        set(gca,'color',[0.9020 0.7608 0.4314]) % corn
    else
        set(gca,'color',[0.8235 0.9294 0.6863]) % soybean
    end
end

% ---------------- US-Ne3 ----------------
for k = 1:6
    axes(ha(k+6))
    starttime = datetime(2017+k,5,1);
    endtime   = datetime(2017+k,11,1);
    idx = simtable_Ne3.date >= starttime & simtable_Ne3.date <= endtime;

    plot(simtable_Ne3.day(idx)+1, simtable_Ne3.g_gpp(idx),'r-','LineWidth',1.5); hold on
    plot(simtable_Ne3_ori.day(idx)+1, simtable_Ne3_ori.g_gpp(idx),'b-','LineWidth',1.5);
    plot(simtable_Ne3.day(idx)+1, simtable_Ne3.GPPobs(idx),'ko','MarkerFaceColor','k','MarkerSize',3);
    xlim([120,310]); ylim([0,50]);
    xlabel('DOY','FontSize',12)

    if k==1
        ylabel('GPP (g C m^{-2} d^{-1})','FontSize',12)
    else
        set(gca,'YTickLabel',[])
    end

    title(sprintf('US-Ne3 (%d)',k+2017),'FontSize',16)
    res1 = evaluate_regression_metrics(simtable_Ne3.g_gpp(idx),simtable_Ne3.GPPobs(idx));
    res2 = evaluate_regression_metrics(simtable_Ne3_ori.g_gpp(idx),simtable_Ne3_ori.GPPobs(idx));
    text(121,43,sprintf('{\\itR}^2 = %2.2f (%2.2f)\nRMSE = %2.2f (%2.2f)\nNRMSE = %2.2f%% (%2.2f%%)',...
        res1.R2,res2.R2,res1.RMSE,res2.RMSE,res1.NRMSE*100,res2.NRMSE*100));

    if mod(k,2)==0
        set(gca,'color',[0.9020 0.7608 0.4314])
    else
        set(gca,'color',[0.8235 0.9294 0.6863])
    end
end

% Panel labels (a–l)
for i = 1:12
    axes(ha(i))
    text(284,44,sprintf('(%s)','a'+i-1),'FontSize',16)
end
set(gcf,'Position',[86 227 1551 662])
exportgraphics(gcf,'.\Fig\Fig_GPP.jpg','Resolution',300)

%% ========================================================================
%  Plot ET (Evapotranspiration)
%% ========================================================================
figure
ha = tight_subplot(2,6,[0.1,0.01],[0.08,0.05],[0.04,0.01]);

% ---------------- US-Ne2 ----------------
for k = 1:6
    axes(ha(k))
    starttime = datetime(2017+k,5,1);
    endtime   = datetime(2017+k,11,1);
    idx = simtable_Ne2.date >= starttime & simtable_Ne2.date <= endtime;

    plot(simtable_Ne2.day(idx)+1, simtable_Ne2.ET(idx),'r-','LineWidth',1.5); hold on
    plot(simtable_Ne2_ori.day(idx)+1, simtable_Ne2_ori.ET(idx),'b-','LineWidth',1.5);
    plot(simtable_Ne2.day(idx)+1, simtable_Ne2.ETobs(idx),'ko','MarkerFaceColor','k','MarkerSize',3);
    xlim([120,310]); ylim([0,12]); xlabel('DOY','FontSize',12)

    if k==1
        ylabel('ET (mm d^{-1})','FontSize',12)
        le=legend('DLEM-SIF','DLEM-Ag'); le.Position=[0.042 0.7913 0.0677 0.0521]; le.Box='off';
    else
        set(gca,'YTickLabel',[])
    end

    title(sprintf('US-Ne2 (%d)',k+2017),'FontSize',16)
    res1 = evaluate_regression_metrics(simtable_Ne2.ET(idx),simtable_Ne2.ETobs(idx));
    res2 = evaluate_regression_metrics(simtable_Ne2_ori.ET(idx),simtable_Ne2_ori.ETobs(idx));
    text(121,10,sprintf('{\\itR}^2 = %2.2f (%2.2f)\nRMSE = %2.2f (%2.2f)\nNRMSE = %2.2f%% (%2.2f%%)',...
        res1.R2,res2.R2,res1.RMSE,res2.RMSE,res1.NRMSE*100,res2.NRMSE*100));

    if mod(k,2)==0
        set(gca,'color',[0.9020 0.7608 0.4314])
    else
        set(gca,'color',[0.8235 0.9294 0.6863])
    end
end

% ---------------- US-Ne3 ----------------
for k = 1:6
    axes(ha(k+6))
    starttime = datetime(2017+k,5,1);
    endtime   = datetime(2017+k,11,1);
    idx = simtable_Ne3.date >= starttime & simtable_Ne3.date <= endtime;

    plot(simtable_Ne3.day(idx)+1, simtable_Ne3.ET(idx),'r-','LineWidth',1.5); hold on
    plot(simtable_Ne3_ori.day(idx)+1, simtable_Ne3_ori.ET(idx),'b-','LineWidth',1.5);
    plot(simtable_Ne3.day(idx)+1, simtable_Ne3.ETobs(idx),'ko','MarkerFaceColor','k','MarkerSize',3);
    xlim([120,310]); ylim([0,12]); xlabel('DOY','FontSize',12)

    if k==1
        ylabel('ET (mm d^{-1})','FontSize',12)
    else
        set(gca,'YTickLabel',[])
    end

    title(sprintf('US-Ne3 (%d)',k+2017),'FontSize',16)
    res1 = evaluate_regression_metrics(simtable_Ne3.ET(idx),simtable_Ne3.ETobs(idx));
    res2 = evaluate_regression_metrics(simtable_Ne3_ori.ET(idx),simtable_Ne3_ori.ETobs(idx));
    text(121,10,sprintf('{\\itR}^2 = %2.2f (%2.2f)\nRMSE = %2.2f (%2.2f)\nNRMSE = %2.2f%% (%2.2f%%)',...
        res1.R2,res2.R2,res1.RMSE,res2.RMSE,res1.NRMSE*100,res2.NRMSE*100));

    if mod(k,2)==0
        set(gca,'color',[0.9020 0.7608 0.4314])
    else
        set(gca,'color',[0.8235 0.9294 0.6863])
    end
end

% Panel labels
for i = 1:12
    axes(ha(i))
    text(284,11,sprintf('(%s)','a'+i-1),'FontSize',16)
end
set(gcf,'Position',[86 227 1551 662])
exportgraphics(gcf,'.\Fig\Fig_ET.jpg','Resolution',300)

%% ========================================================================
%  Plot WUE (Water Use Efficiency)
%% ========================================================================
figure
ha = tight_subplot(2,6,[0.1,0.01],[0.08,0.05],[0.04,0.01]);

% ---------------- US-Ne2 ----------------
for k = 1:6
    axes(ha(k))
    starttime = datetime(2017+k,5,1);
    endtime   = datetime(2017+k,11,1);
    idx = simtable_Ne2.date >= starttime & simtable_Ne2.date <= endtime;

    plot(simtable_Ne2.day(idx)+1, simtable_Ne2.WUE(idx),'r-','LineWidth',1.5); hold on
    plot(simtable_Ne2_ori.day(idx)+1, simtable_Ne2_ori.WUE(idx),'b-','LineWidth',1.5);
    plot(simtable_Ne2.day(idx)+1, simtable_Ne2.WUEobs(idx),'ko','MarkerFaceColor','k','MarkerSize',3);
    xlim([120,310]); ylim([0,12]); xlabel('DOY','FontSize',12)

    if k==1
        ylabel('WUE (g C m^{-2} mm^{-1})','FontSize',12)
        le=legend('DLEM-SIF','DLEM-Ag'); le.Position=[0.042 0.7913 0.0677 0.0521]; le.Box='off';
    else
        set(gca,'YTickLabel',[])
    end

    title(sprintf('US-Ne2 (%d)',k+2017),'FontSize',16)
    res1 = evaluate_regression_metrics(simtable_Ne2.WUE(idx),simtable_Ne2.WUEobs(idx));
    res2 = evaluate_regression_metrics(simtable_Ne2_ori.WUE(idx),simtable_Ne2_ori.WUEobs(idx));
       text(121,10.5,sprintf('{\\itR}^2 = %2.2f (%2.2f)\nRMSE = %2.2f (%2.2f)\nNRMSE = %2.2f%% (%2.2f%%)',...
        res1.R2,res2.R2,res1.RMSE,res2.RMSE,res1.NRMSE*100,res2.NRMSE*100));

    % Alternate background color: corn vs soybean
    if mod(k,2)==0
        set(gca,'color',[0.9020 0.7608 0.4314]) % corn
    else
        set(gca,'color',[0.8235 0.9294 0.6863]) % soybean
    end
end

% ---------------- US-Ne3 ----------------
for k = 1:6
    axes(ha(k+6))
    starttime = datetime(2017+k,5,1);
    endtime   = datetime(2017+k,11,1);
    idx = simtable_Ne3.date >= starttime & simtable_Ne3.date <= endtime;

    plot(simtable_Ne3.day(idx)+1, simtable_Ne3.WUE(idx),'r-','LineWidth',1.5); hold on
    plot(simtable_Ne3_ori.day(idx)+1, simtable_Ne3_ori.WUE(idx),'b-','LineWidth',1.5);
    plot(simtable_Ne3.day(idx)+1, simtable_Ne3.WUEobs(idx),'ko','MarkerFaceColor','k','MarkerSize',3);
    xlim([120,310]); ylim([0,12]); xlabel('DOY','FontSize',12)

    if k==1
        ylabel('WUE (g C m^{-2} mm^{-1})','FontSize',12)
    else
        set(gca,'YTickLabel',[])
    end

    title(sprintf('US-Ne3 (%d)',k+2017),'FontSize',16)
    res1 = evaluate_regression_metrics(simtable_Ne3.WUE(idx),simtable_Ne3.WUEobs(idx));
    res2 = evaluate_regression_metrics(simtable_Ne3_ori.WUE(idx),simtable_Ne3_ori.WUEobs(idx));
    text(121,10.5,sprintf('{\\itR}^2 = %2.2f (%2.2f)\nRMSE = %2.2f (%2.2f)\nNRMSE = %2.2f%% (%2.2f%%)',...
        res1.R2,res2.R2,res1.RMSE,res2.RMSE,res1.NRMSE*100,res2.NRMSE*100));

    if mod(k,2)==0
        set(gca,'color',[0.9020 0.7608 0.4314])
    else
        set(gca,'color',[0.8235 0.9294 0.6863])
    end
end

% Panel labels (a–l)
for i = 1:12
    axes(ha(i))
    text(284,11,sprintf('(%s)','a'+i-1),'FontSize',16)
end

% Set figure position and export
set(gcf,'Position',[86 227 1551 662])
exportgraphics(gcf,'.\Fig\Fig_WUE.jpg','Resolution',300)

%% ========================================================================
%  End of Script
%% ========================================================================



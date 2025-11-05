
%% ========================================================================
%  Comparison of DLEM-SIF and DLEM-Ag Simulations in the U.S. Midwest (2018–2023)
%  ------------------------------------------------------------------------
%  This script:
%   1. Loads DLEM-SIF and DLEM-Ag simulation outputs
%   2. Computes regional GPP, ET, and WUE for each year (2018–2023)
%   3. Produces spatial maps for each variable and model
%   4. Calculates statistics (Mean, SD, CV)
%   5. Analyzes correlations between SPEI and model outputs
%   6. Generates multi-year mean and statistical comparison plots
%
%  ------------------------------------------------------------------------
%% ========================================================================
close all
clear
clc

%% -------------------- Load Regional Data ---------------------------------
load .\Data\regionaldata.mat   % Contains names, sitesall_t, states, spei_m, etc.

%% -------------------- Initialize Variables -------------------------------
GPP_sif = nan(12069,6);
ET_sif  = nan(12069,6);
WUE_sif = nan(12069,6);

GPP_ag  = nan(12069,6);
ET_ag   = nan(12069,6);
WUE_ag  = nan(12069,6);

%% ========================================================================
%  1. Read DLEM-SIF Output Files (2018–2023)
%% ========================================================================
for yeard = 1:6
    filename = sprintf('.\\Data\\Midwest_DLEM_SIF_y%d.dlem', 2017 + yeard);

    % Read integer block
    fileid = fopen(filename, 'r');
    ini_data = fread(fileid, inf, "int");
    fclose(fileid);
    iniall_data = reshape(ini_data, [], 12069);

    % Read float block
    fileid = fopen(filename, 'r');
    ini_data = fread(fileid, inf, "float32");
    fclose(fileid);
    float_data = reshape(ini_data, [], 12069);

    % Replace first 5 integer rows
    float_data(1:5,:) = iniall_data(1:5,:);

    % Convert to table and extract variables
    restable = array2table(float_data', 'VariableNames', names.name);
    GPP_sif(:, yeard) = restable.g_gpp;
    ET_sif(:, yeard)  = restable.g_evap + restable.g_trans;
    WUE_sif(:, yeard) = restable.g_gpp ./ (restable.g_evap + restable.g_trans);
end

%% ========================================================================
%  2. Read DLEM-Ag Output Files (2018–2023)
%% ========================================================================
for yeard = 1:6
    filename = sprintf('.\\Data\\Midwest_DLEM_Ag_y%d.dlem', 2017 + yeard);

    % Read integer block
    fileid = fopen(filename, 'r');
    ini_data = fread(fileid, inf, "int");
    fclose(fileid);
    iniall_data = reshape(ini_data, [], 12069);

    % Read float block
    fileid = fopen(filename, 'r');
    ini_data = fread(fileid, inf, "float32");
    fclose(fileid);
    float_data = reshape(ini_data, [], 12069);

    % Replace first 5 rows
    float_data(1:5,:) = iniall_data(1:5,:);

    % Convert to table and extract variables
    restable = array2table(float_data', 'VariableNames', names.name);
    GPP_ag(:, yeard) = restable.g_gpp;
    ET_ag(:, yeard)  = restable.g_evap + restable.g_trans;
    WUE_ag(:, yeard) = restable.g_gpp ./ (restable.g_evap + restable.g_trans);
end

%% ========================================================================
%  3. Define Grid and Regional Extent
%% ========================================================================
SR = 0.1;                         % Spatial resolution (degrees)
xlon = 360 / SR + 1;              % Global grid columns
xlat = 180 / SR + 1;              % Global grid rows
lonnet = floor((sitesall_t.lon + 180) / SR) + 1;
latnet = floor((sitesall_t.lat + 90) / SR) + 1;
latlim = [35.5, 50];
lonlim = [-105, -80];

%% ========================================================================
%  4. Annual Maps of WUE (2018–2023)
%% ========================================================================
lim = [0.5, 4.5];
mycolor = brewermap(8, "Spectral");
figure
ha = tight_subplot(2,6,[0.05,0.0],[0.12,0.08],[0.0,0.0]);
set(gcf,'Position',[204 313 1550 512])

for yeard = 1:6
    % ------------------- DLEM-SIF -------------------
    gridded = rot90(accumarray([lonnet,latnet],WUE_sif(:,yeard),[xlon,xlat],@mean));
    rasterSize = size(gridded);
    R = georefcells([-90 90],[-180 180],rasterSize,"ColumnsStartFrom","north");
    [xmin,xmax,ymin,ymax] = deal(-104,-80.51,35.99,49.38);
    [B,RB] = geocrop(gridded,R,[ymin ymax],[xmin xmax]);
    B(B==0)=nan;
    [lat,lon] = geographicGrid(RB);
    mean_value = nanmean(WUE_sif(:,yeard));
    std_value  = nanstd(WUE_sif(:,yeard));
    CV_value   = std_value/mean_value*100;

    axes(ha(yeard))
    worldmap(latlim,lonlim); hold on
    pcolorm(lat,lon,B)
    colormap(mycolor)
    framem off; gridm off; mlabel off; plabel off
    clim(lim)
    textm(49,-87,sprintf('Mean=%.2f\nSD=%.2f\nCV=%.2f%%',mean_value,std_value,CV_value),...
        'FontSize',10,'Fontname','Times New Roman','FontWeight','bold')
    for jj=1:length(states)
        plotm(states(jj).Y,states(jj).X,'k')
    end
    title(sprintf('%d (DLEM-SIF)',yeard+2017),'FontSize',14)
end

for yeard = 1:6
    % ------------------- DLEM-Ag -------------------
    gridded = rot90(accumarray([lonnet,latnet],WUE_ag(:,yeard),[xlon,xlat],@mean));
    rasterSize = size(gridded);
    R = georefcells([-90 90],[-180 180],rasterSize,"ColumnsStartFrom","north");
    [xmin,xmax,ymin,ymax] = deal(-104,-80.51,35.99,49.38);
    [B,RB] = geocrop(gridded,R,[ymin ymax],[xmin xmax]);
    B(B==0)=nan;
    [lat,lon] = geographicGrid(RB);
    mean_value = nanmean(WUE_ag(:,yeard));
    std_value  = nanstd(WUE_ag(:,yeard));
    CV_value   = std_value/mean_value*100;

    axes(ha(yeard+6))
    worldmap(latlim,lonlim); hold on
    pcolorm(lat,lon,B)
    colormap(mycolor)
    framem off; gridm off; mlabel off; plabel off
    clim(lim)
    textm(49,-87,sprintf('Mean=%.2f\nSD=%.2f\nCV=%.2f%%',mean_value,std_value,CV_value),...
        'FontSize',10,'Fontname','Times New Roman','FontWeight','bold')
    for jj=1:length(states)
        plotm(states(jj).Y,states(jj).X,'k')
    end
    title(sprintf('%d (DLEM-Ag)',yeard+2017),'FontSize',14)
end

% Add colorbar
h = colorbar;
h.Location = 'south';
h.Position = [0.15 0.05 0.7 0.0374];
h.FontSize = 12; h.FontWeight = 'bold';
h.Label.String = 'WUE (g C m^{-2} mm^{-1})';
h.Label.Position = [2.3633 2.3 0];
h.Label.FontWeight = 'bold';
exportgraphics(gcf,'.\Fig\Fig_Midwest_WUE.jpg','Resolution',300)

%% ========================================================================
%  5. Annual Maps of ET (2018–2023)
%% ========================================================================
lim = [200,800];
mycolor = brewermap(6,"RdYlBu");
figure
ha = tight_subplot(2,6,[0.05,0.0],[0.12,0.08],[0.0,0.0]);
set(gcf,'Position',[204 313 1550 512])

for yeard = 1:6
    gridded = rot90(accumarray([lonnet,latnet],ET_sif(:,yeard),[xlon,xlat],@mean));
    rasterSize = size(gridded);
    R = georefcells([-90 90],[-180 180],rasterSize,"ColumnsStartFrom","north");
    [xmin,xmax,ymin,ymax] = deal(-104,-80.51,35.99,49.38);
    [B,RB] = geocrop(gridded,R,[ymin ymax],[xmin xmax]);
    B(B==0)=nan;
    [lat,lon] = geographicGrid(RB);
    mean_value = nanmean(ET_sif(:,yeard));
    std_value  = nanstd(ET_sif(:,yeard));
    CV_value   = std_value/mean_value*100;

    axes(ha(yeard))
    worldmap(latlim,lonlim)
    hold on
    pcolorm(lat,lon,B)
    colormap(mycolor)
    framem off; gridm off; mlabel off; plabel off
    clim(lim)
    textm(49,-87,sprintf('Mean=%.2f\nSD=%.2f\nCV=%.2f%%',mean_value,std_value,CV_value),...
        'FontSize',10,'Fontname','Times New Roman','FontWeight','bold')
    for jj=1:length(states)
        plotm(states(jj).Y,states(jj).X,'k')
    end
    title(sprintf('%d (DLEM-SIF)',yeard+2017),'FontSize',14)
end

for yeard = 1:6
    gridded = rot90(accumarray([lonnet,latnet],ET_ag(:,yeard),[xlon,xlat],@mean));
    rasterSize = size(gridded);
    R = georefcells([-90 90],[-180 180],rasterSize,"ColumnsStartFrom","north");
    [xmin,xmax,ymin,ymax] = deal(-104,-80.51,35.99,49.38);
    [B,RB] = geocrop(gridded,R,[ymin ymax],[xmin xmax]);
    B(B==0)=nan;
    [lat,lon] = geographicGrid(RB);
    mean_value = nanmean(ET_ag(:,yeard));
    std_value  = nanstd(ET_ag(:,yeard));
    CV_value   = std_value/mean_value*100;

    axes(ha(yeard+6))
    worldmap(latlim,lonlim)
    hold on
    pcolorm(lat,lon,B)
    colormap(mycolor)
    framem off; gridm off; mlabel off; plabel off
    clim(lim)
    textm(49,-87,sprintf('Mean=%.2f\nSD=%.2f\nCV=%.2f%%',mean_value,std_value,CV_value),...
        'FontSize',10,'Fontname','Times New Roman','FontWeight','bold')
    for jj=1:length(states)
        plotm(states(jj).Y,states(jj).X,'k')
    end
    title(sprintf('%d (DLEM-Ag)',yeard+2017),'FontSize',14)
end

h=colorbar;
h.Location='south';
h.Position=[0.15 0.05 0.7 0.0374];
h.FontSize=12; h.FontWeight='bold';
h.Label.String='ET (mm yr^{-1})';
h.Label.Position=[500.0003   2.5         0];
h.Label.FontWeight='bold';
exportgraphics(gcf,'.\Fig\Fig_Midwest_ET.jpg','Resolution',300)

%% ========================================================================
%  6. Annual Maps of GPP (2018–2023)
%% ========================================================================
lim=[0,2500];
mycolor=brewermap(5,"RdYlGn");
figure
ha=tight_subplot(2,6,[0.05,0.0],[0.12,0.08],[0.0,0.001]);
set(gcf,'Position',[204 313 1550 512])

for yeard=1:6
    gridded=rot90(accumarray([lonnet,latnet],GPP_sif(:,yeard),[xlon,xlat],@mean));
    rasterSize=size(gridded);
    R=georefcells([-90 90],[-180 180],rasterSize,"ColumnsStartFrom","north");
    [xmin,xmax,ymin,ymax]=deal(-104,-80.51,35.99,49.38);
    [B,RB]=geocrop(gridded,R,[ymin ymax],[xmin xmax]);
    B(B==0)=nan;
    [lat,lon]=geographicGrid(RB);
    mean_value=nanmean(GPP_sif(:,yeard));
    std_value=nanstd(GPP_sif(:,yeard));
    CV_value=std_value/mean_value*100;

    axes(ha(yeard))
    worldmap(latlim,lonlim)
    hold on
    pcolorm(lat,lon,B)
    colormap(mycolor)
    framem off; gridm off; mlabel off; plabel off
    clim(lim)
    textm(49,-87.5,sprintf('Mean=%.2f\nSD=%.2f\nCV=%.2f%%',mean_value,std_value,CV_value),...
        'FontSize',10,'Fontname','Times New Roman','FontWeight','bold')
    for jj=1:length(states)
        plotm(states(jj).Y,states(jj).X,'k')
    end
    title(sprintf('%d (DLEM-SIF)',yeard+2017),'FontSize',14)
end

for yeard=1:6
    gridded=rot90(accumarray([lonnet,latnet],GPP_ag(:,yeard),[xlon,xlat],@mean));
    rasterSize=size(gridded);
    R=georefcells([-90 90],[-180 180],rasterSize,"ColumnsStartFrom","north");
    [xmin,xmax,ymin,ymax]=deal(-104,-80.51,35.99,49.38);
    [B,RB]=geocrop(gridded,R,[ymin ymax],[xmin xmax]);
    B(B==0)=nan;
    [lat,lon]=geographicGrid(RB);
    mean_value=nanmean(GPP_ag(:,yeard));
    std_value=nanstd(GPP_ag(:,yeard));
    CV_value=std_value/mean_value*100;

    axes(ha(yeard+6))
    worldmap(latlim,lonlim)
    hold on
    pcolorm(lat,lon,B)
    colormap(mycolor)
    framem off; gridm off; mlabel off; plabel off
    clim(lim)
    textm(49,-87.5,sprintf('Mean=%.2f\nSD=%.2f\nCV=%.2f%%',mean_value,std_value,CV_value),...
        'FontSize',10,'Fontname','Times New Roman','FontWeight','bold')
    for jj=1:length(states)
        plotm(states(jj).Y,states(jj).X,'k')
    end
    title(sprintf('%d (DLEM-Ag)',yeard+2017),'FontSize',14)
end

h=colorbar;
h.Location='south';
h.Position=[0.15 0.05 0.7 0.0374];
h.FontSize=12; h.FontWeight='bold';
h.Label.Position=[1250   2.5         0];
h.Label.String='GPP (g C m^{-2} yr^{-1})';
h.Label.FontWeight='bold';
exportgraphics(gcf,'.\Fig\Fig_Midwest_GPP.jpg','Resolution',300)


%% ========================================================================
%  7. Multi-Year Mean Maps for GPP, ET, and WUE
%% ========================================================================

%% ---------------- Multi-Year Mean GPP ------------------------------------
lim = [0, 2500];
mycolor = brewermap(5, "RdYlGn");
figure
ha = tight_subplot(1, 2, [0.05, 0.0], [0.12, 0.08], [0.0, 0.0]);
set(gcf, 'Position', [655 313 1099 512])

% ----- DLEM-SIF -----
gridded = rot90(accumarray([lonnet, latnet], nanmean(GPP_sif, 2), [xlon, xlat], @mean));
rasterSize = size(gridded);
R = georefcells([-90 90], [-180 180], rasterSize, "ColumnsStartFrom", "north");
[xmin, xmax, ymin, ymax] = deal(-104, -80.51, 35.99, 49.38); % Midwest extent
[B, RB] = geocrop(gridded, R, [ymin ymax], [xmin xmax]);
B(B == 0) = nan;
[lat, lon] = geographicGrid(RB);

axes(ha(1))
worldmap(latlim, lonlim)
hold on
pcolorm(lat, lon, B)
colormap(mycolor)
framem off; gridm off; mlabel off; plabel off
clim(lim)
for jj = 1:length(states)
    plotm(states(jj).Y, states(jj).X, 'k') % State boundaries
end
title('DLEM-SIF', 'FontSize', 14)

% ----- DLEM-Ag -----
gridded = rot90(accumarray([lonnet, latnet], nanmean(GPP_ag, 2), [xlon, xlat], @mean));
[B, RB] = geocrop(gridded, R, [ymin ymax], [xmin xmax]);
B(B == 0) = nan;
[lat, lon] = geographicGrid(RB);

axes(ha(2))
worldmap(latlim, lonlim)
hold on
pcolorm(lat, lon, B)
colormap(mycolor)
framem off; gridm off; mlabel off; plabel off
clim(lim)
for jj = 1:length(states)
    plotm(states(jj).Y, states(jj).X, 'k')
end
title('DLEM-Ag', 'FontSize', 14)

% Colorbar
h = colorbar;
h.Location = 'south';
h.Position = [0.15 0.05 0.7 0.0374];
h.FontSize = 12;
h.FontWeight = 'bold';
h.Label.Position=[1250   2.5         0];
h.Label.String = 'GPP (g C m^{-2} yr^{-1})';
h.Label.FontWeight = 'bold';
exportgraphics(gcf, '.\Fig\Fig_Midwest_GPP_all.jpg', 'Resolution', 300)


%% ---------------- Multi-Year Mean ET -------------------------------------
lim = [200, 800];
mycolor = brewermap(6, "RdYlBu");
figure
ha = tight_subplot(1, 2, [0.05, 0.0], [0.12, 0.08], [0.0, 0.0]);
set(gcf, 'Position', [655 313 1099 512])

% ----- DLEM-SIF -----
gridded = rot90(accumarray([lonnet, latnet], nanmean(ET_sif, 2), [xlon, xlat], @mean));
[B, RB] = geocrop(gridded, R, [ymin ymax], [xmin xmax]);
B(B == 0) = nan;
[lat, lon] = geographicGrid(RB);

axes(ha(1))
worldmap(latlim, lonlim)
hold on
pcolorm(lat, lon, B)
colormap(mycolor)
framem off; gridm off; mlabel off; plabel off
clim(lim)
for jj = 1:length(states)
    plotm(states(jj).Y, states(jj).X, 'k')
end
title('DLEM-SIF', 'FontSize', 14)

% ----- DLEM-Ag -----
gridded = rot90(accumarray([lonnet, latnet], nanmean(ET_ag, 2), [xlon, xlat], @mean));
[B, RB] = geocrop(gridded, R, [ymin ymax], [xmin xmax]);
B(B == 0) = nan;
[lat, lon] = geographicGrid(RB);

axes(ha(2))
worldmap(latlim, lonlim)
hold on
pcolorm(lat, lon, B)
colormap(mycolor)
framem off; gridm off; mlabel off; plabel off
clim(lim)
for jj = 1:length(states)
    plotm(states(jj).Y, states(jj).X, 'k')
end
title('DLEM-Ag', 'FontSize', 14)

% Colorbar
h = colorbar;
h.Location = 'south';
h.Position = [0.15 0.05 0.7 0.0374];
h.FontSize = 12;
h.FontWeight = 'bold';
h.Label.Position=[500   2.5         0];
h.Label.String = 'ET (mm yr^{-1})';
h.Label.FontWeight = 'bold';
exportgraphics(gcf, '.\Fig\Fig_Midwest_ET_all.jpg', 'Resolution', 300)


%% ---------------- Multi-Year Mean WUE ------------------------------------
lim = [0.5, 4.5];
mycolor = brewermap(8, "Spectral");
figure
ha = tight_subplot(1, 2, [0.05, 0.0], [0.12, 0.08], [0.0, 0.0]);
set(gcf, 'Position', [655 313 1099 512])

% ----- DLEM-SIF -----
gridded = rot90(accumarray([lonnet, latnet], nanmean(WUE_sif, 2), [xlon, xlat], @mean));
[B, RB] = geocrop(gridded, R, [ymin ymax], [xmin xmax]);
B(B == 0) = nan;
[lat, lon] = geographicGrid(RB);

axes(ha(1))
worldmap(latlim, lonlim)
hold on
pcolorm(lat, lon, B)
colormap(mycolor)
framem off; gridm off; mlabel off; plabel off
clim(lim)
for jj = 1:length(states)
    plotm(states(jj).Y, states(jj).X, 'k')
end
title('DLEM-SIF', 'FontSize', 14)

% ----- DLEM-Ag -----
gridded = rot90(accumarray([lonnet, latnet], nanmean(WUE_ag, 2), [xlon, xlat], @mean));
[B, RB] = geocrop(gridded, R, [ymin ymax], [xmin xmax]);
B(B == 0) = nan;
[lat, lon] = geographicGrid(RB);

axes(ha(2))
worldmap(latlim, lonlim)
hold on
pcolorm(lat, lon, B)
colormap(mycolor)
framem off; gridm off; mlabel off; plabel off
clim(lim)
for jj = 1:length(states)
    plotm(states(jj).Y, states(jj).X, 'k')
end
title('DLEM-Ag', 'FontSize', 14)

% Colorbar
h = colorbar;
h.Location = 'south';
h.Position = [0.15 0.05 0.7 0.0374];
h.FontSize = 12;
h.FontWeight = 'bold';
h.Label.Position=[2.5   2.5         0];
h.Label.String = 'WUE (g C m^{-2} mm^{-1})';
h.Label.FontWeight = 'bold';
exportgraphics(gcf, '.\Fig\Fig_Midwest_WUE_all.jpg', 'Resolution', 300)



%% ========================================================================
%  8. Annual SPEI (2018–2023)
%% ========================================================================
close all
lim = [-1.5, 1.5];
mycolor = flipud(jet);
figure
ha = tight_subplot(2,3,[0.05,0.0],[0.12,0.08],[0.0,0.0]);
set(gcf,'Position',[204 279 817 546])

for yeard = 1:6
    gridded = rot90(accumarray([lonnet,latnet],spei_m(:,yeard),[xlon,xlat],@mean));
    rasterSize = size(gridded);
    R = georefcells([-90 90],[-180 180],rasterSize,"ColumnsStartFrom","north");
    [xmin,xmax,ymin,ymax] = deal(-104,-80.51,35.99,49.38);
    [B,RB] = geocrop(gridded,R,[ymin ymax],[xmin xmax]);
    B(B==0)=nan;
    [lat,lon] = geographicGrid(RB);
    mean_value = nanmean(spei_m(:,yeard));
    std_value  = nanstd(spei_m(:,yeard));

    axes(ha(yeard))
    worldmap(latlim,lonlim)
    hold on
    pcolorm(lat,lon,B)
    colormap(mycolor)
    framem off; gridm off; mlabel off; plabel off
    clim(lim)
    textm(49,-87.5,sprintf('Mean=%.2f\nSD=%.2f',mean_value,std_value),...
        'FontSize',10,'Fontname','Times New Roman','FontWeight','bold')
    for jj=1:length(states)
        plotm(states(jj).Y,states(jj).X,'k')
    end
    title(sprintf('%d',yeard+2017),'FontSize',14)
end

h=colorbar;
h.Location='south';
h.Position=[0.15 0.05 0.7 0.0374];
h.FontSize=12; h.FontWeight='bold';
h.Label.Position=[0,2,0];
h.Label.String='SPEI';
h.Label.FontWeight='bold';
exportgraphics(gcf,'.\Fig\Fig_Midwest_SPEI.jpg','Resolution',300)

%% ========================================================================
%  9. Correlation Between SPEI and WUE/ET/GPP
%% ========================================================================
for i = 1:12069
    disp(i)
    [rho_sif(i,1),pval_sif(i,1)] = corr(spei_m(i,:)',WUE_sif(i,:)');
    [rho_sif(i,2),pval_sif(i,2)] = corr(spei_m(i,:)',ET_sif(i,:)');
    [rho_sif(i,3),pval_sif(i,3)] = corr(spei_m(i,:)',GPP_sif(i,:)');
    [rho_ag(i,1),pval_ag(i,1)]   = corr(spei_m(i,:)',WUE_ag(i,:)');
    [rho_ag(i,2),pval_ag(i,2)]   = corr(spei_m(i,:)',ET_ag(i,:)');
    [rho_ag(i,3),pval_ag(i,3)]   = corr(spei_m(i,:)',GPP_ag(i,:)');
end

labels={'WUE','ET','GPP'};
figure
ha=tight_subplot(2,3,[0.05,0.0],[0.12,0.08],[0.0,0.0]);
set(gcf,'Position',[396 133 1160 740])

for k=1:3
    % ---------------- DLEM-SIF ----------------
    axes(ha(k))
    lim=[-1,1];
    mycolor=brewermap(10,"RdYlGn");
    gridded=rot90(accumarray([lonnet,latnet],rho_sif(:,k),[xlon,xlat],@mean));
    R=georefcells([-90 90],[-180 180],size(gridded),"ColumnsStartFrom","north");
    [B,RB]=geocrop(gridded,R,[35.99 49.38],[-104 -80.51]);
    B(B==0)=nan;
    [lat,lon]=geographicGrid(RB);
    worldmap(latlim,lonlim)
    hold on
    pcolorm(lat,lon,B)
    colormap(mycolor)
    framem off; gridm off; mlabel off; plabel off
    clim(lim)
    for jj=1:length(states)
        plotm(states(jj).Y,states(jj).X,'k')
    end
    % Mark p<0.05
    gridded_p = rot90(accumarray([lonnet,latnet],pval_sif(:,k),[xlon,xlat],@mean));
    [B_p,RB_p]=geocrop(gridded_p,R,[35.99 49.38],[-104 -80.51]);
    B_p(B_p==0)=nan;
    pindex = B_p<0.05;
    plotm(lat(pindex)-0.05,lon(pindex)+0.05,"k*","MarkerSize",0.5)
    title(strcat(labels{k},' (DLEM-SIF)'),'FontSize',14)
    textm(49,-83,sprintf('(%s)','a'+k-1),'FontSize',16,'FontWeight','Bold')
end

for k=1:3
    % ---------------- DLEM-Ag ----------------
    axes(ha(k+3))
    lim=[-1,1];
    mycolor=brewermap(10,"RdYlGn");
    gridded=rot90(accumarray([lonnet,latnet],rho_ag(:,k),[xlon,xlat],@mean));
    R=georefcells([-90 90],[-180 180],size(gridded),"ColumnsStartFrom","north");
    [B,RB]=geocrop(gridded,R,[35.99 49.38],[-104 -80.51]);
    B(B==0)=nan;
    [lat,lon]=geographicGrid(RB);
    worldmap(latlim,lonlim)
    hold on
    pcolorm(lat,lon,B)
    colormap(mycolor)
    framem off; gridm off; mlabel off; plabel off
    clim(lim)
    for jj=1:length(states)
        plotm(states(jj).Y,states(jj).X,'k')
    end
    % Mark p<0.05
    gridded_p = rot90(accumarray([lonnet,latnet],pval_ag(:,k),[xlon,xlat],@mean));
    [B_p,RB_p]=geocrop(gridded_p,R,[35.99 49.38],[-104 -80.51]);
    B_p(B_p==0)=nan;
    pindex = B_p<0.05;
    plotm(lat(pindex)-0.05,lon(pindex)+0.05,"k*","MarkerSize",0.5)
    title(strcat(labels{k},' (DLEM-Ag)'),'FontSize',14)
    textm(49,-83,sprintf('(%s)','a'+k-1+3),'FontSize',16,'FontWeight','Bold')
end

h=colorbar;
h.Location='south';
h.Position=[0.1 0.05 0.8 0.0175];
h.FontSize=16; h.FontWeight='bold';
h.Label.String='\itr';
h.Label.FontWeight='bold';
h.Label.Position=[0,3,0];
exportgraphics(gcf,'.\Fig\Fig_Midwest_SPEI_relationship.jpg','Resolution',300)

%% ========================================================================
%  10. WUE Distribution Histogram
%% ========================================================================
figure
hold on
h1 = histogram(WUE_sif(:),'Normalization','probability','BinWidth',0.2,...
    'FaceColor',[0.2 0.6 0.8],'EdgeColor','k','FaceAlpha',0.6);
h2 = histogram(WUE_ag(:),'Normalization','probability','BinWidth',0.2,...
    'FaceColor',[1.0 0.5 0.1],'EdgeColor','k','FaceAlpha',0.6);
xlim([0,5])
ylim([0,max([h1.Values,h2.Values])*1.1])
xlabel('WUE (g C m^{-2} mm^{-1})')
ylabel('Relative Frequency')
legend({'DLEM-SIF','DLEM-Ag'},'Location','northeast')
set(gca,'FontSize',12,'LineWidth',1)
grid on; box on; hold off
exportgraphics(gcf,'.\Fig\Fig_frequency.jpg','Resolution',300)

%% ========================================================================
%  11. Binned SPEI–WUE Relationships (2018–2023)
%% ========================================================================
spei_wue = [spei_m(:),WUE_sif(:),WUE_ag(:)];
spei_wue_t = array2table(spei_wue,"VariableNames",{'SPEI','WUE_SIF','WUE_AG'});
spei_wue_g = groupsummary(spei_wue_t,"SPEI",[-1.5:0.25:1.5],["mean","std"]);

for i=1:height(spei_wue_g)
    str1=strsplit(erase(string(spei_wue_g.disc_SPEI(i)),{'[','(',')',']'}),',');
    spei_wue_g.spei(i)=str2double(str1(1));
end
spei_wue_g(13,:) = [];

p_values = nan(height(spei_wue_g),1);
for i=1:height(spei_wue_g)
    idx = spei_wue_t.SPEI >= spei_wue_g.spei(i) & spei_wue_t.SPEI < spei_wue_g.spei(i)+0.25;
    data_SIF = spei_wue_t.WUE_SIF(idx);
    data_AG  = spei_wue_t.WUE_AG(idx);
    if numel(data_SIF)>=3 && numel(data_AG)>=3
        [~,p]=ttest2(data_SIF,data_AG);
        p_values(i)=p;
    end
end

x = spei_wue_g.spei;
mean_SIF = spei_wue_g.mean_WUE_SIF;
mean_AG  = spei_wue_g.mean_WUE_AG;
std_SIF  = spei_wue_g.std_WUE_SIF;
std_AG   = spei_wue_g.std_WUE_AG;
x_labels = string(spei_wue_g.disc_SPEI);

figure
hb = bar(x,[mean_SIF,mean_AG],'grouped','BarWidth',0.8);
hold on
hb(1).FaceColor=[0,0,0]; hb(1).FaceAlpha=0.7;
hb(2).FaceColor=[0.8,0.8,0.8]; hb(2).FaceAlpha=0.7;
xOffset = [-0.033,0.033];

mean_mat = [mean_SIF,mean_AG];
std_mat  = [std_SIF,std_AG];

for i=1:2
    errorbar(x+xOffset(i),mean_mat(:,i),std_mat(:,i),...
        'k','LineStyle','none','LineWidth',1.2,'CapSize',4);
end

for i=1:length(p_values)
    if ~isnan(p_values(i)) && p_values(i)<0.05
        x1 = x(i)+xOffset(1);
        x2 = x(i)+xOffset(2);
        y  = max(mean_mat(i,:)+std_mat(i,:))+0.1;
        plot([x1,x1,x2,x2],[y,y+0.05,y+0.05,y],'k','LineWidth',1.2)
        if p_values(i)<0.001
            stars='***';
        elseif p_values(i)<0.01
            stars='**';
        else
            stars='*';
        end
        text(mean([x1,x2]),y+0.07,stars,'HorizontalAlignment','center','FontSize',12)
    end
end

xlabel('SPEI Range')
ylabel('WUE (g C m^{-2} mm^{-1})')
legend({'DLEM-SIF','DLEM-Ag'},'Location','northeast')
set(gca,'FontSize',12,'LineWidth',1)
grid on; box on
ylim([0,5])
set(gcf,'Position',[680 375 790 503])
exportgraphics(gcf,'.\Fig\Fig_SPEI_WUE_allyear_p3.jpg','Resolution',300)

%% ========================================================================
%  End of Script


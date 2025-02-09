clc;clear;
%% Bootstrapping 
nrow=1000;
low_randlist=zeros(nrow,10);
high_randlist=zeros(nrow,10);
for i=1:nrow
    low_randlist(i,:) = randi([1 39],1,10) ;
    high_randlist(i,:) = randi([1 39],1,10) ;
end
%%
load("Adult_vital_1980_2018.mat"); %load adult time series 1965-2020
 Ms_mean = invlogit(logit(M_s_mean) + M_s_epsilon(16:54));
%%
uwind = ncread('ERA5_1959_2022_wind_t.nc','u',[1 1 1 1],[Inf Inf 1 Inf]);
 vwind = ncread('ERA5_1959_2022_wind_t.nc','v',[1 1 1 1],[Inf Inf 1 Inf]); 
 lon = ncread('ERA5_1959_2022_wind_t.nc','longitude');
 lat = ncread('ERA5_1959_2022_wind_t.nc','latitude'); 
 lat = lat(360:end);
 uwind = squeeze(uwind);  uwind = uwind(:,360:end,:);
 vwind = squeeze(vwind);  vwind = vwind(:,360:end,:);
%%
uwind_survival = zeros(1440,362,63);
for i = 1:63
     uwind_survival (:,:,i) = mean(uwind(:,:,((i-1)*12+1):((i-1)*12+12)),3);
end
uwind_survival = uwind_survival(:,:,22:60);
clear uwind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vwind_survival = zeros(1440,362,63);
for i = 1:63
     vwind_survival (:,:,i) = mean(vwind(:,:,((i-1)*12+1):((i-1)*12+12)),3);
end
vwind_survival = vwind_survival(:,:,22:60);
%%%%%%%%%%%%%%%%%%
windspeed_survival = sqrt(uwind_survival.^2 + vwind_survival.^2);
[LN,LT]=meshgrid(lon,lat);
clear  vwind
%% load sea level pressure data 
 %ncdisp('D:\Ruijiao\WA_personality_project\ENV analysis\spatial_correlation\climate_data\ERA5_1959_2022_SLP.nc');
 slp = ncread('climate_data/ERA5_1959_2022_SLP.nc','msl',[1 1 1 1],[Inf Inf 1 Inf]);
      slp = squeeze(slp);slp=slp(:,360:end,:);
%Jan 1 in year t-1 - December 1 in year t-1 for survival in year 1959-2021
slp_survival = zeros(1440,362,63);
for i = 1:63
    slp_survival (:,:,i) = mean(slp(:,:,((i-1)*12+1):((i-1)*12+12)),3);
end
slp_survival = slp_survival(:,:,22:60)/100; %take years from 1980-2020
clear slp
 %%
mycolormap = customcolormap([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1], {'#19196b','#3d549a','#638fca','#87caf7','#acdbf7','#fffceb','#feecce','#fbdab0','#f6a690','#eb626a','#de123d'},15);
mycolormap2 = customcolormap([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1], {'#19196b','#3d549a','#638fca','#87caf7','#acdbf7','#ffffff','#feecce','#fbdab0','#f6a690','#eb626a','#de123d'},15);
mean_colormap = customcolormap(linspace(0,1,8), {'#55499A','#5F7DAF','#6BC5B5','#B4E387','#EFED71','#F6AF5E','#E8615D','#A1353C'},15);
%%
fig=figure
fig.Position=[560,240,320,220]
m_proj('miller','lat',[-75 -5],'lon',[0 180]);
m_pcolor(lon,lat,diff_slp');shading interp
m_coast('color',[0.1 0.1 0.1]);hold on

m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),diff_uw_positive(2:30:end, 2:25:end)',diff_vw_positive(2:30:end, 2:25:end)',2,'color','k',"linewidth",1);
m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),diff_uw_negative(2:30:end, 2:25:end)',diff_vw_negative(2:30:end, 2:25:end)',2,'color',[.7 .7 .7],"linewidth",1);
%m_text(38,-53,"Crozet",'fontsize',10);
m_scatter(51.59,-46.25,'LineWidth',1.5,'MarkerEdgeColor',[0 0 0],"MarkerFaceColor",[1 0.8 0])
%mycolormap = customcolormap([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1], {'#19196b','#3d549a','#638fca','#87caf7','#acdbf7','#ffffff','#feecce','#fbdab0','#f6a690','#eb626a','#de123d'},15);
colormap(flipud(mycolormap2))
cb=colorbar('southoutside');
m_grid('linewi',1,'tickdir','out','ticklen',.005,'linestyl','none','color','k');
set(gca, 'Layer', 'top',"FontSize",14,'FontName',"Helvetica"); 
ax = gca;
axpos = ax.Position;
cb.Label.String = "SLP anomaly (mb)";
clim([-1 1]);
cb.Position(4) = 0.5*cb.Position(4); cb.Position(2) = 1.2*cb.Position(2);
ax.Position = axpos;
exportgraphics(fig,'anomaly_composite\wind_Fbs.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'anomaly_composite\wind_Fbs.png','BackgroundColor','none','Resolution',600)

%%
[sd,r]=sort(Ms_mean,'descend');
%% plot mean condition
fig=figure
fig.Position=[560,240,320,220]
uwind_survival_mean=mean(uwind_survival,3);
vwind_survival_mean=mean(vwind_survival,3);
m_proj('miller','lat',[-75 -5],'lon',[0 180]);
m_pcolor(lon,lat,mean(slp_survival,3)');shading interp
m_coast('color',[0.1 0.1 0.1]);hold on
v = [1,0.1];
m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),uwind_survival_mean(2:30:end, 2:25:end)',vwind_survival_mean(2:30:end, 2:25:end)',1.7,'color','k',"linewidth",1.2);
%m_text(38,-53,"Crozet",'fontsize',10);
m_scatter(51.59,-46.25,'LineWidth',1.5,'MarkerEdgeColor',[0 0 0],"MarkerFaceColor",[1 0.8 0])
%mycolormap = customcolormap([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1], {'#19196b','#3d549a','#638fca','#87caf7','#acdbf7','#ffffff','#feecce','#fbdab0','#f6a690','#eb626a','#de123d'},15);
colormap(flipud(mean_colormap))
cb=colorbar('southoutside');
m_grid('linewi',1,'tickdir','out','ticklen',.005,'linestyl','none','color','k');
set(gca, 'Layer', 'top',"FontSize",14,'FontName',"Helvetica"); 
ax = gca;
axpos = ax.Position;
cb.Label.String = "SLP (mb)";
caxis([980 1025]);
cb.Position(4) = 0.5*cb.Position(4); cb.Position(2) = 1.2*cb.Position(2);
ax.Position = axpos;
exportgraphics(fig,'anomaly_composite\wind_Ms_mean.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'anomaly_composite\wind_Ms_mean.png','BackgroundColor','none','Resolution',600)
%% 
uwind_survival_lowmean=mean(uwind_survival(:,:,r(30:39)),3);
vwind_survival_lowmean=mean(vwind_survival(:,:,r(30:39)),3);
uwind_survival_highmean=mean(uwind_survival(:,:,r(1:10)),3);
vwind_survival_highmean=mean(vwind_survival(:,:,r(1:10)),3);
slp_survival_mean=mean(slp_survival,3);

low_anomaly_uw = uwind_survival_lowmean - uwind_survival_mean;
low_anomaly_vw = vwind_survival_lowmean - vwind_survival_mean;
low_anomaly_slp = mean(slp_survival(:,:,r(30:39)),3) - slp_survival_mean;
low_anomaly_windspeed = mean(windspeed_survival(:,:,r(30:39)),3) - mean(windspeed_survival,3);
low_anomaly_windspeed_positive = real(low_anomaly_windspeed>0);
low_anomaly_uw_positive = low_anomaly_uw.*low_anomaly_windspeed_positive; 
low_anomaly_vw_positive = low_anomaly_vw.*low_anomaly_windspeed_positive; 
low_anomaly_windspeed_negative = real(low_anomaly_windspeed<=0);
low_anomaly_uw_negative = low_anomaly_uw.*low_anomaly_windspeed_negative; 
low_anomaly_vw_negative = low_anomaly_vw.*low_anomaly_windspeed_negative;

high_anomaly_uw = uwind_survival_highmean - uwind_survival_mean;
high_anomaly_vw = vwind_survival_highmean - vwind_survival_mean;
high_anomaly_slp = mean(slp_survival(:,:,r(1:10)),3) - slp_survival_mean;
high_anomaly_windspeed = mean(windspeed_survival(:,:,r(1:10)),3) - mean(windspeed_survival,3);
high_anomaly_windspeed_positive = real(high_anomaly_windspeed>0);
high_anomaly_uw_positive = high_anomaly_uw.*high_anomaly_windspeed_positive; 
high_anomaly_vw_positive = high_anomaly_vw.*high_anomaly_windspeed_positive; 
high_anomaly_windspeed_negative = real(high_anomaly_windspeed<=0);
high_anomaly_uw_negative = high_anomaly_uw.*high_anomaly_windspeed_negative; 
high_anomaly_vw_negative = high_anomaly_vw.*high_anomaly_windspeed_negative;

diff_uw = uwind_survival_highmean - uwind_survival_lowmean;
diff_vw = vwind_survival_highmean - vwind_survival_lowmean;
diff_slp = mean(slp_survival(:,:,r(1:10)),3) - mean(slp_survival(:,:,r(30:39)),3);
diff_windspeed = mean(windspeed_survival(:,:,r(1:10)),3) - mean(windspeed_survival(:,:,r(30:39)),3);
diff_windspeed_positive = real(diff_windspeed>=0);
diff_uw_positive = diff_uw.*diff_windspeed_positive; 
diff_vw_positive = diff_vw.*diff_windspeed_positive; 
diff_windspeed_negative = real(diff_windspeed<=0);
diff_uw_negative = diff_uw.*diff_windspeed_negative; 
diff_vw_negative = diff_vw.*diff_windspeed_negative;

    for j = 1:362
        low_rand_uw = zeros(1000,1440);
        low_rand_vw = zeros(1000,1440);
        low_rand_slp = zeros(1000,1440);

        high_rand_uw = zeros(1000,1440);
        high_rand_vw = zeros(1000,1440);
        high_rand_slp = zeros(1000,1440);

        diff_rand_uw = zeros(1000,1440);
        diff_rand_vw = zeros(1000,1440);
        diff_rand_slp = zeros(1000,1440);

        for t = 1:1000
            low_rand_uw(t,:)=mean(uwind_survival(:,j,low_randlist(t,:)),3) - uwind_survival_mean(:,j);
            low_rand_vw(t,:)=mean(vwind_survival(:,j,low_randlist(t,:)),3) - vwind_survival_mean(:,j);
            low_rand_slp(t,:)=mean(slp_survival(:,j,low_randlist(t,:)),3) - slp_survival_mean(:,j);

            high_rand_uw(t,:)=mean(uwind_survival(:,j,high_randlist(t,:)),3) - uwind_survival_mean(:,j);
            high_rand_vw(t,:)=mean(vwind_survival(:,j,high_randlist(t,:)),3) - vwind_survival_mean(:,j);
            high_rand_slp(t,:)=mean(slp_survival(:,j,high_randlist(t,:)),3) - slp_survival_mean(:,j);

            diff_rand_uw(t,:) = mean(uwind_survival(:,j,high_randlist(t,:)),3) - mean(uwind_survival(:,j,low_randlist(t,:)),3);
            diff_rand_vw(t,:) = mean(vwind_survival(:,j,high_randlist(t,:)),3) - mean(vwind_survival(:,j,low_randlist(t,:)),3);
            diff_rand_slp(t,:)= mean(slp_survival(:,j,high_randlist(t,:)),3) - mean(slp_survival(:,j,low_randlist(t,:)),3);
        end
        
        sig_low_anomaly_uw = (low_anomaly_uw(:,j) <= quantile(low_rand_uw, 0.05)'|low_anomaly_uw(:,j) >= quantile(low_rand_uw, 0.95)');
        sig_low_anomaly_vw = (low_anomaly_vw(:,j) <= quantile(low_rand_vw, 0.05)'|low_anomaly_vw(:,j) >= quantile(low_rand_vw, 0.95)');
        sig_low_anomaly_wind(:,j) =double( sig_low_anomaly_uw | sig_low_anomaly_vw);
        sig_low_anomaly_slp(:,j) = double((low_anomaly_slp(:,j) <= quantile(low_rand_slp, 0.05)')|(low_anomaly_slp(:,j) >= quantile(low_rand_slp, 0.95)'));

        sig_high_anomaly_uw = (high_anomaly_uw(:,j) >= quantile(high_rand_uw, 0.95)')|(high_anomaly_uw(:,j) <= quantile(high_rand_uw, 0.05)');
        sig_high_anomaly_vw = (high_anomaly_vw(:,j) >= quantile(high_rand_vw, 0.95)')|(high_anomaly_vw(:,j) <= quantile(high_rand_vw, 0.05)');
        sig_high_anomaly_wind(:,j) =double( sig_high_anomaly_uw | sig_high_anomaly_vw);
        sig_high_anomaly_slp(:,j) = double((high_anomaly_slp(:,j) >= quantile(high_rand_slp, 0.95)')|(high_anomaly_slp(:,j) <= quantile(high_rand_slp, 0.05)'));

        sig_diff_uw = (diff_uw(:,j) <= quantile(diff_rand_uw, 0.05)'|diff_uw(:,j) >= quantile(diff_rand_uw, 0.95)');
        sig_diff_vw = (diff_vw(:,j) <= quantile(diff_rand_vw, 0.05)'|diff_vw(:,j) >= quantile(diff_rand_vw, 0.95)');
        sig_diff_wind(:,j) =double( sig_diff_uw | sig_diff_vw);
        sig_diff_slp(:,j) = double(diff_slp(:,j) <= quantile(diff_rand_slp, 0.05)'|diff_slp(:,j) >= quantile(diff_rand_slp, 0.95)');
    end
clear low_rand_uw low_rand_vw low_rand_slp high_rand_uw high_rand_vw high_rand_slp diff_rand_uw diff_rand_vw diff_rand_slp
clear sig_diff_uw  sig_diff_vw sig_high_anomaly_uw sig_high_anomaly_vw sig_low_anomaly_uw sig_low_anomaly_vw
%%
scale_factor=2/max([max(diff_windspeed,[],'all'),max(high_anomaly_windspeed,[],'all'),max(low_anomaly_windspeed,[],'all')]);
%%
fig=figure
fig.Position=[560,240,320,220]
m_proj('miller','lat',[-75 -5],'lon',[0 180]);
%m_pcolor(lon,lat,(low_anomaly_slp.*sig_low_anomaly_slp)');shading interp
m_pcolor(lon,lat,low_anomaly_slp');shading interp;hold on
sig_lon=sig_low_anomaly_slp(2:10:end, 2:10:end).*lon(2:10:end); sig_lon(sig_lon==0)=NaN;
sig_lat=sig_low_anomaly_slp(2:10:end, 2:10:end).*lat(2:10:end)';sig_lat(sig_lat==0)=NaN;
m_scatter(sig_lon,sig_lat,1.5,'LineWidth',.001,'MarkerEdgeColor','none','MarkerFaceColor',[.2 .2 .2]);
m_coast('color',[0.1 0.1 0.1]);
v = [1,0.1];
%m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),(low_anomaly_uw_positive(2:30:end, 2:25:end).*sig_low_anomaly_wind(2:30:end, 2:25:end))',(low_anomaly_vw_positive(2:30:end, 2:25:end).*sig_low_anomaly_wind(2:30:end, 2:25:end))',1.7,'color',[6/256, 184/256, 129/256],"linewidth",1.5);
%m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),(low_anomaly_uw_negative(2:30:end, 2:25:end).*sig_low_anomaly_wind(2:30:end, 2:25:end))',(low_anomaly_vw_negative(2:30:end, 2:25:end).*sig_low_anomaly_wind(2:30:end, 2:25:end))',1.7,'color',[.4 .4 .4],"linewidth",1);
%m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),low_anomaly_uw_positive(2:30:end, 2:25:end)',low_anomaly_vw_positive(2:30:end, 2:25:end)',1.7,'color',[6/256, 184/256, 129/256],"linewidth",1.5);
%m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),low_anomaly_uw_negative(2:30:end, 2:25:end)',low_anomaly_vw_negative(2:30:end, 2:25:end)',1.7,'color',[.4 .4 .4],"linewidth",1);
m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),(low_anomaly_uw(2:30:end, 2:25:end).*sig_low_anomaly_wind(2:30:end, 2:25:end))',(low_anomaly_vw(2:30:end, 2:25:end).*sig_low_anomaly_wind(2:30:end, 2:25:end))',scale_factor,'color',[.0 .0 .0],"linewidth",1,'AutoScale','off')
m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),(low_anomaly_uw(2:30:end, 2:25:end).*(1-sig_low_anomaly_wind(2:30:end, 2:25:end)))',(low_anomaly_vw(2:30:end, 2:25:end).*(1-sig_low_anomaly_wind(2:30:end, 2:25:end)))',scale_factor,'color',[.7 .7 .7],"linewidth",1,'AutoScale','off')
m_scatter(51.59,-46.25,'LineWidth',1.5,'MarkerEdgeColor',[0 0 0],"MarkerFaceColor",[1 0.8 0])
colormap(flipud(mycolormap2))
cb=colorbar('southoutside');
m_grid('linewi',1,'tickdir','out','ticklen',.005,'linestyl','none','color','k');
set(gca, 'Layer', 'top',"FontSize",14,'FontName',"Helvetica"); 
ax = gca;
axpos = ax.Position;
cb.Label.String = "SLP anomaly (mb)";
caxis([-2 2]);
cb.Position(4) = 0.5*cb.Position(4); cb.Position(2) = 1.2*cb.Position(2);
ax.Position = axpos;
exportgraphics(fig,'anomaly_composite\wind_Ms_lowanomaly.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'anomaly_composite\wind_Ms_lowanomaly.png','BackgroundColor','none','Resolution',600)

%% plot high anomaly
fig=figure
fig.Position=[560,240,320,220]
m_proj('miller','lat',[-75 -5],'lon',[0 180]);
m_pcolor(lon,lat,high_anomaly_slp');shading interp;hold on
sig_lon=sig_high_anomaly_slp(2:10:end, 2:10:end).*lon(2:10:end); sig_lon(sig_lon==0)=NaN;
sig_lat=sig_high_anomaly_slp(2:10:end, 2:10:end).*lat(2:10:end)';sig_lat(sig_lat==0)=NaN;
m_scatter(sig_lon,sig_lat,1.5,'LineWidth',.001,'MarkerEdgeColor','none','MarkerFaceColor',[.2 .2 .2]);
m_coast('color',[0.1 0.1 0.1]);hold on
v = [1,0.1];
%m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),(high_anomaly_uw_positive(2:30:end, 2:25:end).*sig_high_anomaly_wind(2:30:end, 2:25:end))',(high_anomaly_vw_positive(2:30:end, 2:25:end).*sig_high_anomaly_wind(2:30:end, 2:25:end))',1.7,'color',[6/256, 184/256, 129/256],"linewidth",1.5);
%m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),(high_anomaly_uw_negative(2:30:end, 2:25:end).*sig_high_anomaly_wind(2:30:end, 2:25:end))',(high_anomaly_vw_negative(2:30:end, 2:25:end).*sig_high_anomaly_wind(2:30:end, 2:25:end))',1.7,'color',[.4 .4 .4],"linewidth",1);
m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),(high_anomaly_uw(2:30:end, 2:25:end).*sig_high_anomaly_wind(2:30:end, 2:25:end))',(high_anomaly_vw(2:30:end, 2:25:end).*sig_high_anomaly_wind(2:30:end, 2:25:end))',scale_factor,'color',[.0 .0 .0],"linewidth",1,'AutoScale','off')
m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),(high_anomaly_uw(2:30:end, 2:25:end).*(1-sig_high_anomaly_wind(2:30:end, 2:25:end)))',(high_anomaly_vw(2:30:end, 2:25:end).*(1-sig_high_anomaly_wind(2:30:end, 2:25:end)))',scale_factor,'color',[.7 .7 .7],"linewidth",1,'AutoScale','off')
m_scatter(51.59,-46.25,'LineWidth',1.5,'MarkerEdgeColor',[0 0 0],"MarkerFaceColor",[1 0.8 0])
colormap(flipud(mycolormap2))
cb=colorbar('southoutside');
m_grid('linewi',1,'tickdir','out','ticklen',.005,'linestyl','none','color','k');
set(gca, 'Layer', 'top',"FontSize",14,'FontName',"Helvetica"); 
ax = gca;
axpos = ax.Position;
cb.Label.String = "SLP anomaly (mb)";
clim([-2 2]);
cb.Position(4) = 0.5*cb.Position(4); cb.Position(2) = 1.2*cb.Position(2);
ax.Position = axpos;
exportgraphics(fig,'anomaly_composite\wind_Ms_highanomaly.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'anomaly_composite\wind_Ms_highanomaly.png','BackgroundColor','none','Resolution',600)
%%
fig=figure
fig.Position=[560,240,320,220]
m_proj('miller','lat',[-75 -5],'lon',[0 180]);
m_coast('color',[0.1 0.1 0.1]);hold on
%m_pcolor(lon,lat,(diff_slp.*sig_diff_slp)');shading interp;hold on
m_pcolor(lon,lat,diff_slp');shading interp;hold on
sig_lon=sig_diff_slp(2:10:end, 2:10:end).*lon(2:10:end); sig_lon(sig_lon==0)=NaN;
sig_lat=sig_diff_slp(2:10:end, 2:10:end).*lat(2:10:end)';sig_lat(sig_lat==0)=NaN;
m_scatter(sig_lon,sig_lat,1.5,'LineWidth',.001,'MarkerEdgeColor','none','MarkerFaceColor',[.2 .2 .2]);
m_coast('color',[0.1 0.1 0.1]);hold on
m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),(diff_uw_positive(2:30:end, 2:25:end).*sig_diff_wind(2:30:end, 2:25:end))',(diff_vw_positive(2:30:end, 2:25:end).*sig_diff_wind(2:30:end, 2:25:end))',scale_factor,'color',[40/256, 74/256, 247/256],"linewidth",1,'AutoScale','off');
m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),(diff_uw_negative(2:30:end, 2:25:end).*sig_diff_wind(2:30:end, 2:25:end))',(diff_vw_negative(2:30:end, 2:25:end).*sig_diff_wind(2:30:end, 2:25:end))',scale_factor,'color',[.0 .0 .0],"linewidth",1,'AutoScale','off');
m_quiver(LN(2:25:end, 2:30:end),LT(2:25:end, 2:30:end),(diff_uw(2:30:end, 2:25:end).*(1-sig_diff_wind(2:30:end, 2:25:end)))',(diff_vw(2:30:end, 2:25:end).*(1-sig_diff_wind(2:30:end, 2:25:end)))',scale_factor,'color',[.7 .7 .7],"linewidth",1,'AutoScale','off');
%m_text(38,-53,"Crozet",'fontsize',10);
m_scatter(51.59,-46.25,'LineWidth',1.5,'MarkerEdgeColor',[0 0 0],"MarkerFaceColor",[1 0.8 0])
%mycolormap = customcolormap([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1], {'#19196b','#3d549a','#638fca','#87caf7','#acdbf7','#ffffff','#feecce','#fbdab0','#f6a690','#eb626a','#de123d'},15);
colormap(flipud(mycolormap2))
cb=colorbar('southoutside');
m_grid('linewi',1,'tickdir','out','ticklen',.005,'linestyl','none','color','k');
ax = gca;
axpos = ax.Position;
cb.Label.String = "SLP anomaly (mb)";
clim([-2 2]);
cb.Position(4) = 0.5*cb.Position(4); cb.Position(2) = 1.2*cb.Position(2);
ax.Position = axpos;
set(gca, 'Layer', 'top',"FontSize",14,'FontName',"Helvetica"); 
exportgraphics(fig,'anomaly_composite\wind_Ms.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'anomaly_composite\wind_Ms.png','BackgroundColor','none','Resolution',600) 
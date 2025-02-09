clc;clear;
%generate covariate for mascarene heights
ncdisp('ERA5_1959_2022_SLP.nc');
%%
%time period cover 1959-2022
slp = ncread('ERA5_1959_2022_SLP.nc',"msl",[1 1 1 1],[Inf Inf 1 Inf]);
lon = ncread('ERA5_1959_2022_SLP.nc','longitude');
lat = ncread('ERA5_1959_2022_SLP.nc','latitude');
slp = squeeze(slp);
% slp during months for survival probabiity
slp_survival = zeros(1440,721,63);
for i = 1:63
    slp_survival (:,:,i) = mean( slp(:,:,((i-1)*12+1):((i-1)*12+12)),3);
end
% slp during months for breeding probability
slp_breed = zeros(1440,721,63);
for i = 1:63
    slp_breed (:,:,i) = mean( slp(:,:,((i-1)*12+9):((i-1)*12+12)),3);
end
% slp during months for breeding success probability
slp_success = zeros(1440,721,63);
for i = 1:63
    slp_success (:,:,i) = mean( slp(:,:,((i-1)*12+1):((i-1)*12+10)),3);
end
%%
% set range for mascarene high contour
mask_lat = find((lat>=-50).*(lat<=-20)==1);
mask_lon = find( (lon>=30).*(lon<=120)==1);
slp_survival=slp_survival(mask_lon,mask_lat,:);
slp_breed=slp_breed(mask_lon,mask_lat,:);
slp_success=slp_success(mask_lon,mask_lat,:);

%% Mh indices for survival probability 
s_mh_strength = [];
s_mh_lat = [];
s_mh_lon = [];
s_mh_max = [];
for i = 7:62
    temp = slp_survival(:,:,i)>=102000; % MH are area with SLP above 1020.0 millibars contour
    slp_contoured = slp_survival(:,:,i).*temp;
    slp_contoured(slp_contoured==0)=NaN;
    value = nanmean(slp_contoured,"all"); 
    s_mh_strength = [s_mh_strength;value];
    peak = max(slp_contoured,[],"all");
    s_mh_max = [s_mh_max;peak];
    [mh_lon,mh_lat] = find(slp_contoured==peak);
    mh_lat = lat(mask_lat(mh_lat(1)));
    mh_lon = lon(mask_lon(mh_lon(1)));
    s_mh_lat = [s_mh_lat;mh_lat];
    s_mh_lon = [s_mh_lon;mh_lon];
end
%% Mh indices for breeding probability
b_mh_strength = [];
b_mh_lat = [];
b_mh_lon = [];
b_mh_max = [];
for i = 7:62
    temp = slp_breed(:,:,i)>=102000;
    slp_contoured = slp_breed(:,:,i).*temp; % MH are area with SLP above 1020.0 millibars contour
    slp_contoured(slp_contoured==0)=NaN;
    value = nanmean(slp_contoured,"all"); 
    b_mh_strength = [b_mh_strength;value];
    peak = max(slp_contoured,[],"all");
    b_mh_max = [b_mh_max;peak];
    [mh_lon,mh_lat] = find(slp_contoured==peak);
    mh_lat = lat(mask_lat(mh_lat(1)));
    mh_lon = lon(mask_lon(mh_lon(1)));
    b_mh_lat = [b_mh_lat;mh_lat];
    b_mh_lon = [b_mh_lon;mh_lon];
end
%% Mh indices for breeding success
bs_mh_strength = []; 
bs_mh_lat = [];
bs_mh_lon = [];
bs_mh_max = [];
for i = 8:63
    temp = slp_success(:,:,i)>=102000;
    slp_contoured = slp_success(:,:,i).*temp;% MH are area with SLP above 1020.0 millibars contour
    slp_contoured(slp_contoured==0)=NaN;
    value = nanmean(slp_contoured,"all"); 
    bs_mh_strength = [bs_mh_strength;value];
    peak = max(slp_contoured,[],"all");
    bs_mh_max = [bs_mh_max;peak];
    [mh_lon,mh_lat] = find(slp_contoured==peak);
    mh_lat = lat(mask_lat(mh_lat(1)));
    mh_lon = lon(mask_lon(mh_lon(1)));
    bs_mh_lat = [bs_mh_lat;mh_lat];
    bs_mh_lon = [bs_mh_lon;mh_lon];
end
%%
clearvars -except b_mh_lon b_mh_lat b_mh_strength s_mh_lon s_mh_lat s_mh_strength bs_mh_lat bs_mh_lon bs_mh_strength
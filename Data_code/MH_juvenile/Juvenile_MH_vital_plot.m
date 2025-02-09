%% 
clc;clear;
%%
nrow=200;
load posterior_juvenile_ENV.mat
post_dist=zeros(nrow,7500);
for i=1:nrow
    post_dist(i,:) = randi([1 7500],1,7500) ;
end
cov = linspace(-3,3,nrow)';
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fs_ENSO_slope,"all");hold on
mu=Fs_ENSO;
points = invlogit(mu(post_dist)+Fs_ENSO_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[245 132 66]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[245 132 66]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[245 132 66]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[250 211 147]/256)
box off;
ylim([0 1]);
xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("ENSO")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/FJs_ENSO.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/FJs_ENSO.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fs_IOD_slope,"all");hold on
mu=Fs_IOD;
points = invlogit(mu(post_dist)+Fs_IOD_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[245 132 66]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[245 132 66]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[245 132 66]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[250 211 147]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("IOD")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/FJs_IOD.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/FJs_IOD.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fs_mhm_slope,"all");hold on
mu=Fs_mhm;
points = invlogit(mu(post_dist)+Fs_mhm_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[245 132 66]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[245 132 66]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[245 132 66]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[250 211 147]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("MH strength")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/FJs_mhm.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/FJs_mhm.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fs_mhlon_slope,"all");hold on
mu=Fs_mhlon;
points = invlogit(mu(post_dist)+Fs_mhlon_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[245 132 66]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[245 132 66]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[245 132 66]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[250 211 147]/256)
box off;ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("MH longitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/FJs_mhlon.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/FJs_mhlon.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fs_mhlat_slope,"all");hold on
mu=Fs_mhlat;
points = invlogit(mu(post_dist)+Fs_mhlat_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[245 132 66]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[245 132 66]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[245 132 66]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[250 211 147]/256)
box off;ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("MH latitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/FJs_mhlat.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/FJs_mhlat.png','BackgroundColor','none','Resolution',600)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fb_ENSO_slope,"all");hold on
mu=Fb_ENSO;
points = invlogit(mu(post_dist)+Fb_ENSO_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[245 132 66]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[245 132 66]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[245 132 66]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[250 211 147]/256)
box off;
ylim([0 1]);
xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("ENSO")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/FJb_ENSO.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/FJb_ENSO.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fb_IOD_slope,"all");hold on
mu=Fb_IOD;
points = invlogit(mu(post_dist)+Fb_IOD_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[245 132 66]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[245 132 66]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[245 132 66]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[250 211 147]/256)
box off;
ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("IOD")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/FJb_IOD.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/FJb_IOD.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fb_mhm_slope,"all");hold on
mu=Fb_mhm;
points = invlogit(mu(post_dist)+Fb_mhm_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[245 132 66]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[245 132 66]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[245 132 66]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[250 211 147]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("MH strength")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/FJb_mhm.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/FJb_mhm.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fb_mhlon_slope,"all");hold on
mu=Fb_mhlon;
points = invlogit(mu(post_dist)+Fb_mhlon_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[245 132 66]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[245 132 66]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[245 132 66]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[250 211 147]/256)
box off;ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("MH longitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/FJb_mhlon.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/FJb_mhlon.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fb_mhlat_slope,"all");hold on
mu=Fb_mhlat;
points = invlogit(mu(post_dist)+Fb_mhlat_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[245 132 66]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[245 132 66]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[245 132 66]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[250 211 147]/256)
box off; ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("MH latitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/FJb_mhlat.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/FJb_mhlat.png','BackgroundColor','none','Resolution',600)
%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fbs_ENSO_slope,"all");hold on
mu=Fbs_ENSO;
points = invlogit(mu(post_dist)+Fbs_ENSO_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[245 132 66]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[245 132 66]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[245 132 66]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[250 211 147]/256)
box off; ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("ENSO")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/FJbs_ENSO.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/FJbs_ENSO.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fbs_IOD_slope,"all");hold on
mu=Fbs_IOD;
points = invlogit(mu(post_dist)+Fbs_IOD_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[245 132 66]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[245 132 66]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[245 132 66]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[250 211 147]/256)
box off; ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("IOD")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/FJbs_IOD.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/FJbs_IOD.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fbs_mhm_slope,"all");hold on
mu=Fbs_mhm;
points = invlogit(mu(post_dist)+Fbs_mhm_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[245 132 66]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[245 132 66]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[245 132 66]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[250 211 147]/256)
box off; ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("MH strength")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/FJbs_mhm.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/FJbs_mhm.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fbs_mhlon_slope,"all");hold on
mu=Fbs_mhlon;
points = invlogit(mu(post_dist)+Fbs_mhlon_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[245 132 66]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[245 132 66]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[245 132 66]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[250 211 147]/256)
box off; ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("MH longitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/FJbs_mhlon.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/FJbs_mhlon.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fbs_mhlat_slope,"all");hold on
mu=Fbs_mhlat;
points = invlogit(mu(post_dist)+Fbs_mhlat_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[245 132 66]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[245 132 66]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[245 132 66]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[250 211 147]/256)
box off; ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("MH latitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/FJbs_mhlat.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/FJbs_mhlat.png','BackgroundColor','none','Resolution',600)
%%%%%%
%%%%%%
%%
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Ms_ENSO_slope,"all");hold on
mu=Ms_ENSO;
points = invlogit(mu(post_dist)+Ms_ENSO_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[2 156 107]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[2 156 107]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[2 156 107]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[141 252 195]/256)
box off; ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("ENSO")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/MJs_ENSO.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/MJs_ENSO.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Ms_IOD_slope,"all");hold on
mu=Ms_IOD;
points = invlogit(mu(post_dist)+Ms_IOD_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[2 156 107]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[2 156 107]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[2 156 107]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[141 252 195]/256)
box off; ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("IOD")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/MJs_IOD.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/MJs_IOD.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Ms_mhm_slope,"all");hold on
mu=Ms_mhm;
points = invlogit(mu(post_dist)+Ms_mhm_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[2 156 107]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[2 156 107]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[2 156 107]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[141 252 195]/256)
box off;ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("MH strength")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/MJs_mhm.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/MJs_mhm.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Ms_mhlon_slope,"all");hold on
mu=Ms_mhlon;
points = invlogit(mu(post_dist)+Ms_mhlon_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[2 156 107]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[2 156 107]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[2 156 107]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[141 252 195]/256)
box off; ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("MH longitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/MJs_mhlon.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/MJs_mhlon.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Ms_mhlat_slope,"all");hold on
mu=Ms_mhlat;
points = invlogit(mu(post_dist)+Ms_mhlat_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[2 156 107]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[2 156 107]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[2 156 107]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[141 252 195]/256)
box off; ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("MH latitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/MJs_mhlat.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/MJs_mhlat.png','BackgroundColor','none','Resolution',600)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mb_ENSO_slope,"all");hold on
mu=Mb_ENSO;
points = invlogit(mu(post_dist)+Mb_ENSO_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[2 156 107]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[2 156 107]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[2 156 107]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[141 252 195]/256)
box off; ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("ENSO")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/MJb_ENSO.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/MJb_ENSO.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mb_IOD_slope,"all");hold on
mu=Mb_IOD;
points = invlogit(mu(post_dist)+Mb_IOD_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[2 156 107]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[2 156 107]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[2 156 107]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[141 252 195]/256)
box off;ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("IOD")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/MJb_IOD.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/MJb_IOD.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mb_mhm_slope,"all");hold on
mu=Mb_mhm;
points = invlogit(mu(post_dist)+Mb_mhm_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[2 156 107]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[2 156 107]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[2 156 107]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[141 252 195]/256)
box off;ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("MH strength")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/MJb_mhm.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/MJb_mhm.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mb_mhlon_slope,"all");hold on
mu=Mb_mhlon;
points = invlogit(mu(post_dist)+Mb_mhlon_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[2 156 107]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[2 156 107]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[2 156 107]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[141 252 195]/256)
box off;ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("MH longitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/MJb_mhlon.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/MJb_mhlon.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mb_mhlat_slope,"all");hold on
mu=Mb_mhlat;
points = invlogit(mu(post_dist)+Mb_mhlat_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[2 156 107]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[2 156 107]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[2 156 107]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[141 252 195]/256)
box off;ylim([0 1]); xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("MH latitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/MJb_mhlat.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/MJb_mhlat.png','BackgroundColor','none','Resolution',600)
%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mbs_ENSO_slope,"all");hold on
mu=Mbs_ENSO;
points = invlogit(mu(post_dist)+Mbs_ENSO_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[2 156 107]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[2 156 107]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[2 156 107]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[141 252 195]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("ENSO")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/MJbs_ENSO.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/MJbs_ENSO.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mbs_IOD_slope,"all");hold on
mu=Mbs_IOD;
points = invlogit(mu(post_dist)+Mbs_IOD_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[2 156 107]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[2 156 107]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[2 156 107]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[141 252 195]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("IOD")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/MJbs_IOD.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/MJbs_IOD.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mbs_mhm_slope,"all");hold on
mu=Mbs_mhm;
points = invlogit(mu(post_dist)+Mbs_mhm_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[2 156 107]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[2 156 107]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[2 156 107]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[141 252 195]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("MH strength")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/MJbs_mhm.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/MJbs_mhm.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mbs_mhlon_slope,"all");hold on
mu=Mbs_mhlon;
points = invlogit(mu(post_dist)+Mbs_mhlon_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[2 156 107]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[2 156 107]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[2 156 107]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[141 252 195]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("MH longitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/MJbs_mhlon.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/MJbs_mhlon.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mbs_mhlat_slope,"all");hold on
mu=Mbs_mhlat;
points = invlogit(mu(post_dist)+Mbs_mhlat_slope(post_dist).*cov);
mean_line = invlogit(mean(mu)+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[2 156 107]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[2 156 107]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[2 156 107]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[141 252 195]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("MH latitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/MJbs_mhlat.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/MJbs_mhlat.png','BackgroundColor','none','Resolution',600)


%%
function [ x ] = logit(x)

x = log(x./(1-x));

end
function [ x ] = invlogit(x)

x = 1./(1+exp(-x));

end



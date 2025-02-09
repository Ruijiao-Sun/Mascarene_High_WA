%% 
clc;clear;
%%
nrow=200;
load posterior_adult_ENV.mat
post_dist=zeros(nrow,10000);
for i=1:nrow
    post_dist(i,:) = randi([1 18000],1,10000) ;
end
cov = linspace(-3,3,nrow)';
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fs_ENSO_slope,"all");hold on
mu=mean(Fs_ENSO,2);
points = invlogit(logit(mu(post_dist))+Fs_ENSO_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[153 51 51]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[153 51 51]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[153 51 51]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[225 175 175]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("ENSO")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Fs_ENSO.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Fs_ENSO.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fs_IOD_slope,"all");hold on
mu=mean(Fs_IOD,2);
points = invlogit(logit(mu(post_dist))+Fs_IOD_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[153 51 51]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[153 51 51]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[153 51 51]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[225 175 175]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("IOD")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Fs_IOD.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Fs_IOD.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fs_mhm_slope,"all");hold on
mu=mean(Fs_mhm,2);
points = invlogit(logit(mu(post_dist))+Fs_mhm_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[153 51 51]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[153 51 51]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[153 51 51]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[225 175 175]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("MH strength")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Fs_mhm.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Fs_mhm.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fs_mhlon_slope,"all");hold on
mu=mean(Fs_mhlon,2);
points = invlogit(logit(mu(post_dist))+Fs_mhlon_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[153 51 51]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[153 51 51]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[153 51 51]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[225 175 175]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("MH longitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Fs_mhlon.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Fs_mhlon.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fs_mhlat_slope,"all");hold on
mu=mean(Fs_mhlat,2);
points = invlogit(logit(mu(post_dist))+Fs_mhlat_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[153 51 51]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[153 51 51]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[153 51 51]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[225 175 175]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("MH latitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Fs_mhlat.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Fs_mhlat.png','BackgroundColor','none','Resolution',600)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fb_ENSO_slope,"all");hold on
mu=mean(Fb_ENSO,2);
points = invlogit(logit(mu(post_dist))+Fb_ENSO_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[153 51 51]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[153 51 51]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[153 51 51]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[225 175 175]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("ENSO")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Fb_ENSO.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Fb_ENSO.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fb_IOD_slope,"all");hold on
mu=mean(Fb_IOD,2);
points = invlogit(logit(mu(post_dist))+Fb_IOD_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[153 51 51]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[153 51 51]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[153 51 51]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[225 175 175]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("IOD")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Fb_IOD.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Fb_IOD.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fb_mhm_slope,"all");hold on
mu=mean(Fb_mhm,2);
points = invlogit(logit(mu(post_dist))+Fb_mhm_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[153 51 51]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[153 51 51]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[153 51 51]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[225 175 175]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("MH strength")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Fb_mhm.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Fb_mhm.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fb_mhlon_slope,"all");hold on
mu=mean(Fb_mhlon,2);
points = invlogit(logit(mu(post_dist))+Fb_mhlon_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[153 51 51]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[153 51 51]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[153 51 51]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[225 175 175]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("MH longitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Fb_mhlon.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Fb_mhlon.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fb_mhlat_slope,"all");hold on
mu=mean(Fb_mhlat,2);
points = invlogit(logit(mu(post_dist))+Fb_mhlat_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[153 51 51]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[153 51 51]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[153 51 51]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[225 175 175]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("MH latitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Fb_mhlat.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Fb_mhlat.png','BackgroundColor','none','Resolution',600)
%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fbs_ENSO_slope,"all");hold on
mu=mean(Fbs_ENSO,2);
points = invlogit(logit(mu(post_dist))+Fbs_ENSO_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[153 51 51]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[153 51 51]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[153 51 51]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[225 175 175]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("ENSO")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Fbs_ENSO.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Fbs_ENSO.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fbs_IOD_slope,"all");hold on
mu=mean(Fbs_IOD,2);
points = invlogit(logit(mu(post_dist))+Fbs_IOD_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[153 51 51]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[153 51 51]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[153 51 51]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[225 175 175]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("IOD")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Fbs_IOD.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Fbs_IOD.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fbs_mhm_slope,"all");hold on
mu=mean(Fbs_mhm,2);
points = invlogit(logit(mu(post_dist))+Fbs_mhm_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[153 51 51]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[153 51 51]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[153 51 51]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[225 175 175]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("MH strength")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Fbs_mhm.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Fbs_mhm.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fbs_mhlon_slope,"all");hold on
mu=mean(Fbs_mhlon,2);
points = invlogit(logit(mu(post_dist))+Fbs_mhlon_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[153 51 51]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[153 51 51]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[153 51 51]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[225 175 175]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("MH longitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Fbs_mhlon.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Fbs_mhlon.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Fbs_mhlat_slope,"all");hold on
mu=mean(Fbs_mhlat,2);
points = invlogit(logit(mu(post_dist))+Fbs_mhlat_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[153 51 51]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[153 51 51]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[153 51 51]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[225 175 175]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("MH latitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Fbs_mhlat.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Fbs_mhlat.png','BackgroundColor','none','Resolution',600)
%%%%%%
%%%%%%
%%
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Ms_ENSO_slope,"all");hold on
mu=mean(Ms_ENSO,2);
points = invlogit(logit(mu(post_dist))+Ms_ENSO_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[0 102 204]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[0 102 204]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[0 102 204]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[153, 204, 255]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("ENSO")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Ms_ENSO.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Ms_ENSO.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Ms_IOD_slope,"all");hold on
mu=mean(Ms_IOD,2);
points = invlogit(logit(mu(post_dist))+Ms_IOD_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[0 102 204]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[0 102 204]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[0 102 204]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[153, 204, 255]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("IOD")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Ms_IOD.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Ms_IOD.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Ms_mhm_slope,"all");hold on
mu=mean(Ms_mhm,2);
points = invlogit(logit(mu(post_dist))+Ms_mhm_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[0 102 204]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[0 102 204]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[0 102 204]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[153, 204, 255]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("MH strength")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Ms_mhm.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Ms_mhm.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Ms_mhlon_slope,"all");hold on
mu=mean(Ms_mhlon,2);
points = invlogit(logit(mu(post_dist))+Ms_mhlon_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[0 102 204]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[0 102 204]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[0 102 204]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[153, 204, 255]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("MH longitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Ms_mhlon.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Ms_mhlon.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Ms_mhlat_slope,"all");hold on
mu=mean(Ms_mhlat,2);
points = invlogit(logit(mu(post_dist))+Ms_mhlat_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[0 102 204]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[0 102 204]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[0 102 204]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[153, 204, 255]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Survival probability");xlabel("MH latitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Ms_mhlat.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Ms_mhlat.png','BackgroundColor','none','Resolution',600)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mb_ENSO_slope,"all");hold on
mu=mean(Mb_ENSO,2);
points = invlogit(logit(mu(post_dist))+Mb_ENSO_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[0 102 204]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[0 102 204]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[0 102 204]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[153, 204, 255]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("ENSO")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Mb_ENSO.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Mb_ENSO.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mb_IOD_slope,"all");hold on
mu=mean(Mb_IOD,2);
points = invlogit(logit(mu(post_dist))+Mb_IOD_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[0 102 204]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[0 102 204]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[0 102 204]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[153, 204, 255]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("IOD")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Mb_IOD.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Mb_IOD.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mb_mhm_slope,"all");hold on
mu=mean(Mb_mhm,2);
points = invlogit(logit(mu(post_dist))+Mb_mhm_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[0 102 204]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[0 102 204]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[0 102 204]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[153, 204, 255]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("MH strength")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Mb_mhm.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Mb_mhm.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mb_mhlon_slope,"all");hold on
mu=mean(Mb_mhlon,2);
points = invlogit(logit(mu(post_dist))+Mb_mhlon_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[0 102 204]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[0 102 204]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[0 102 204]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[153, 204, 255]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("MH longitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Mb_mhlon.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Mb_mhlon.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mb_mhlat_slope,"all");hold on
mu=mean(Mb_mhlat,2);
points = invlogit(logit(mu(post_dist))+Mb_mhlat_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[0 102 204]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[0 102 204]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[0 102 204]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[153, 204, 255]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding probability");xlabel("MH latitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Mb_mhlat.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Mb_mhlat.png','BackgroundColor','none','Resolution',600)
%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mbs_ENSO_slope,"all");hold on
mu=mean(Mbs_ENSO,2);
points = invlogit(logit(mu(post_dist))+Mbs_ENSO_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[0 102 204]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[0 102 204]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[0 102 204]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[153, 204, 255]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("ENSO")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Mbs_ENSO.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Mbs_ENSO.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mbs_IOD_slope,"all");hold on
mu=mean(Mbs_IOD,2);
points = invlogit(logit(mu(post_dist))+Mbs_IOD_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[0 102 204]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[0 102 204]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[0 102 204]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[153, 204, 255]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("IOD")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Mbs_IOD.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Mbs_IOD.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mbs_mhm_slope,"all");hold on
mu=mean(Mbs_mhm,2);
points = invlogit(logit(mu(post_dist))+Mbs_mhm_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[0 102 204]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[0 102 204]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[0 102 204]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[153, 204, 255]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("MH strength")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Mbs_mhm.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Mbs_mhm.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mbs_mhlon_slope,"all");hold on
mu=mean(Mbs_mhlon,2);
points = invlogit(logit(mu(post_dist))+Mbs_mhlon_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[0 102 204]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[0 102 204]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[0 102 204]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[153, 204, 255]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("MH longitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Mbs_mhlon.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Mbs_mhlon.png','BackgroundColor','none','Resolution',600)
%%
fig=figure("Units","inches","Position",[0 0 3 3])
hold on;
mean_slope=mean(Mbs_mhlat_slope,"all");hold on
mu=mean(Mbs_mhlat,2);
points = invlogit(logit(mu(post_dist))+Mbs_mhlat_slope(post_dist).*cov);
mean_line = invlogit(logit(mean(mu))+mean_slope*cov);
SEM = std(points,0,2);               % Standard Error
ts = tinv([0.025  0.975],nrow-1);      % T-Score
CI = mean_line + ts.*SEM;                      % Confidence Intervals
plot (cov, mean_line,"LineWidth",2,"Color",[0 102 204]/256);
plot(cov(1:3:end),CI((1:3:end),1),"--","Color",[0 102 204]/256,"LineWidth",0.3)
plot(cov(1:3:end),CI((1:3:end),2),"--","Color",[0 102 204]/256,"LineWidth",0.3)
ciplot(CI(:,1), CI(:,2), cov,[153, 204, 255]/256)
box off;ylim([0 1]);xlim([-3 3]);xticks([-3 0 3])
pbaspect([1 1 1])
ylabel("Breeding success");xlabel("MH latitude")
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica');
set(gca,'TickDir','out'); 
exportgraphics(fig,'figure/Mbs_mhlat.pdf','BackgroundColor','none','Resolution',600)
exportgraphics(fig,'figure/Mbs_mhlat.png','BackgroundColor','none','Resolution',600)
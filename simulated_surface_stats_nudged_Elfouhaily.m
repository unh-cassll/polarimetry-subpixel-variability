% Artificially "nudges" the gravity-capillary peak in Elfouhaily et al.'s
% wavenumber spectrum to be higher or lower, then uses those modified
% spectra to create simulated sea surfaces from which we model light
% reflection and compute our usual error metric of RMSE angle
%
% N. Laxague 2024
%
function simulated_surface_stats_nudged_Elfouhaily(fignum)

addpath _codes/

load('simulated_surface_comparison_stats_wavenumber_mult_U10_10_scale2048.mat')

% Set up limits and ticks/ticklabels
RMSE_lims = 2*[1e-1 1e0];
blocksize_lims = [1e-1 1e2]*2;
block_ticks = [1 10 100];
block_ticklabels = {'1','10','100'};

RMSE_ticks = [(1:9)*0.1 1 2];
RMSE_ticklabels = {'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1','2'};

% Prepare wavenumber spectra
Bk = k.^3.*Fk;
Bk_copy = Bk;
Bk_copy(k<50,:) = NaN;

% Preallocate variables of interest
num_mult = size(Bk,2);
kp = NaN*ones(num_mult,1);
lp_rmse = [kp kp];
p_rmse = lp_rmse;

% Some plotting options: label location and curve transparency value
text_x = 0.05;
text_y = 0.96;
alpha_val = 0.7;
load('coolwarm.mat')
cmap = interp1(1:size(coolwarm,1), coolwarm, linspace(1,size(coolwarm,1),num_mult), 'linear');

% Loop over the number of "nudges" to grab the RMSE and peak scale
for m = 1:num_mult

    ind = find(Bk(:,m)==max(Bk_copy(:,m),[],'omitnan'),1,'first');

    inds = ind-5:ind+5;
    p = polyfit(log10(k(inds)),log10(Bk(inds,m)),2);
    kp(m) = 10^(-p(2)/(2*p(1)));

    rmse_vec = squeeze(RMSE_mat(5,m,:));

    ind = find(rmse_vec==max(rmse_vec),1,'first');
    inds = ind-5:ind+5;
    p = polyfit(log10(full_resolution_m_px*block_vec(inds)),log10(rmse_vec(inds)),2);
    lp_rmse_bit = -p(2)/(2*p(1));
    p_rmse_bit = lp_rmse_bit^2*p(1)+lp_rmse_bit*p(2)+p(3);
    lp_rmse(m,1) = 10^(lp_rmse_bit);
    p_rmse(m,1) = 10^p_rmse_bit;

    rmse_vec = squeeze(RMSE_mat(6,m,:));

    ind = find(rmse_vec==max(rmse_vec),1,'first');
    inds = ind-5:ind+5;
    p = polyfit(log10(full_resolution_m_px*block_vec(inds)),log10(rmse_vec(inds)),2);
    lp_rmse_bit = -p(2)/(2*p(1));
    p_rmse_bit = lp_rmse_bit^2*p(1)+lp_rmse_bit*p(2)+p(3);
    lp_rmse(m,2) = 10^(lp_rmse_bit);
    p_rmse(m,2) = 10^p_rmse_bit;

end

% Construct cell array containing the lengthscale of the short wave peak
lp = 2*pi./kp;
lp_strings = {};
for m = 1:length(lp)
    lp_strings{m} = ['\lambda_p = ' sprintf('%0.1f',1000*lp(m)) ' mm'];
end

% First plot: wavenumber spectra
figure(fignum);clf
hold on
loglog(k,k.^3.*Fk,'k-','linewidth',4)
H_color = loglog(k,k.^3.*fliplr(Fk),'linewidth',3);colororder(flipud(cmap));
hold off
legend(H_color,flip(lp_strings),'Location','southwest')
box on
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
xlim([1e0 1e4])
ylim([5e-5 5e-2])
pbaspect([1 1 1])
xlabel('k [rad m^{-1}]')
ylabel('B(k) [rad]')
text(text_x,text_y,'(a)','FontSize',20,'HorizontalAlignment','center','Units','normalized')

% Second plot: RMSE as a function of block size
figure(fignum*10);clf
hold on
for i = 1:length(kp)
    plot(2*pi/kp(i)*1000*[1 1],[RMSE_lims(2) p_rmse(i,2)],'-','Color',[0 0 0 alpha_val*0.5],'linewidth',3)
    plot(2*pi/kp(i)*1000*[1 1],[RMSE_lims(2) p_rmse(i,2)],'-','Color',[cmap(i,:) alpha_val],'linewidth',2)
    plot(lp_rmse(i,2)*1000*[1 1],[p_rmse(i,2) RMSE_lims(1)],'-','Color',[0 0 0 alpha_val*0.5],'linewidth',3)
    plot(lp_rmse(i,2)*1000*[1 1],[p_rmse(i,2) RMSE_lims(1)],'-','Color',[cmap(i,:) alpha_val],'linewidth',2)
end
plot(1000*full_resolution_m_px*block_vec,squeeze(RMSE_mat(6,:,:)),'k.','markersize',30)
plot(1000*full_resolution_m_px*block_vec,squeeze(RMSE_mat(6,:,:)),'.','markersize',25)
hold off
box on
text(26,1.8,'peak B(k)','FontSize',20,'HorizontalAlignment','left')
text(22,0.55,'peak RMSE','FontSize',20,'HorizontalAlignment','left')
colororder(cmap)
ax=gca;
ax.XTick = block_ticks;
ax.XTickLabels = block_ticklabels;
ax.YTick = RMSE_ticks;
ax.YTickLabels = RMSE_ticklabels;
ax.XScale='log';
ax.YScale='log';
ylim(RMSE_lims)
xlim(blocksize_lims)
pbaspect([1 1 1])
xlabel('block size [mm]')
ylabel('RMSE, \theta [\circ]')
text(text_x,text_y,'(b)','FontSize',20,'HorizontalAlignment','center','Units','normalized')

% Third plot: scatter comparing peak length scales of spectra and RMSE
figure(fignum*100);clf
hold on
for n = 1:length(lp)
    plot(1000*lp_rmse(n,1),1000*lp(n),'o','markerfacecolor',cmap(n,:),'markeredgecolor','k','linewidth',1.5)
    plot(1000*lp_rmse(n,2),1000*lp(n),'s','markerfacecolor',cmap(n,:),'markeredgecolor','k','linewidth',1.5)
end
plot([0 30],[0 30],'k--','linewidth',2)
box on
hold off
xlim([0 30])
ylim([0 30])
pbaspect([1 1 1])
ylabel('\lambda_p, saturation spectrum [mm]')
xlabel('\lambda_p, RMSE [mm]')
legend(char(966),'\theta','Location','southeast')
ax=gca;
ax.XTick = 0:5:50;
ax.YTick = 0:5:50;
text(text_x,text_y,'(c)','FontSize',20,'HorizontalAlignment','center','Units','normalized')
% Plot observed spectra computed from block-averaged intensity fields
%
% Also plot percent difference in spectral energy density as a function of
% normalized wavenumber over a range of block sizes and wind forcing levels
%
% L. Hogan and N. Laxague 2024
%
function plot_block_averaged_spectra_comparison(fignum,labelsize)

% Pluck out our spectra
load("ASIT_Subpixel_Spect.mat")
S_all = S_all(:,:,2:end);
S_all = S_all./(S_all(:,:,1)./S_all(:,1,1));
k_all = k_all(:,:,2:end);
B_all = S_all.*k_all;
n_levels = size(S_all,2);

% Compute k*, the normalized reference wavenumber
k_all(k_all==0) = NaN;
kmin = min(k_all,[],3);
kmax = maxk(k_all,2,3);kmax = kmax(:,:,2);
kstar = (k_all-kmin)./(kmax-kmin);
kvec = squeeze(k_all(1,1,:));

% Pluck out run to show as example
ind = 36;
B_particular = squeeze(B_all(ind,:,:));

% Prepare colormap, figure label positions, and colorbar ticks/ticklabels
cmap = flipud(spectral(n_levels));
text_x = 0.95;
text_y = 0.95;
kstar_lims = [0 1];%+[-0.05 0.05];
cbar_ticklabels = {};
for i = 1:n_levels
    cbar_ticklabels{i} = sprintf('%0.1f',1000*pi./kmax(1,i));
end

% First plot: wavenumber saturation spectra
figure(fignum);clf
hold on
plot(kvec,B_particular,'k-','linewidth',3)
plot(kvec,B_particular,'linewidth',2)
hold off
box on
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
colororder(cmap)
xlim([4e0 4e3])
ylim([2e-5 2e-2])
pbaspect([1 1 1])
xlabel('k [rad m^{-1}]')
ylabel('B(k) [rad]')
colormap(cmap)
clim([0 1]+[-1 1]/(n_levels-1)/2)
cbar = colorbar;
cbar.Location = 'northoutside';
cbar.Ticks = (0:n_levels-1)/(n_levels-1);
cbar.TickLabels = cbar_ticklabels;
set(get(cbar,'label'),'string','Pixel Scale [mm]');

% Second plot
figure(fignum+1)
tlayout = tiledlayout(2,1);

original = repmat(S_all(:,1,:),[1 5 1]);
percentdiff = 100*(S_all(:,3:end,:)-original)./original;
percentdiff(percentdiff<-99.9) = NaN;
percentdiff_mean = squeeze(mean(percentdiff,1));
kstar_mean = squeeze(mean(kstar(:,3:end,:)));

% Second plot, first panel: percent difference with respect to pixel size
nexttile(1)
hold on
plot([0 1],[0 0],'k--','linewidth',1.5)
for i = 1:size(percentdiff_mean,1)
plot(kstar_mean(i,:),percentdiff_mean(i,:),'k-','linewidth',3)
plot(kstar_mean(i,:),percentdiff_mean(i,:),'color',cmap(i+2,:),'linewidth',2)
end
hold off
box on
ax_struc(1).ax = gca;
xlim(kstar_lims)
ylim([-100 20])
xlabel('k*')
ylabel('percent difference [%]')
colormap(cmap)
clim([0 1]+[-1 1]/(n_levels-1)/2)
cbar = colorbar;
cbar.Location = 'eastoutside';
cbar.Ticks = (0:n_levels-1)/(n_levels-1);
cbar.TickLabels = cbar_ticklabels;
set(get(cbar,'label'),'string','Pixel Scale [mm]');
text(text_x,text_y,'(a)','FontSize',labelsize,'HorizontalAlignment','center','Units','normalized')

% Second plot, second panel: percent difference with respect to wind forcing
ustar_lims = 0:0.1:0.5;
nbins = length(ustar_lims) - 1;

bin_cmap = viridis(nbins);

ind = 2;
k_mean = kstar_mean(ind,:);

nexttile(2)
hold on
plot([0 1],[0 0],'k--','linewidth',1.5)
for i = 1:nbins
    ustar_inds = find(ustar>ustar_lims(i) & ustar<=ustar_lims(i+1));
    percentdiff_mean_slice = squeeze(mean(percentdiff(ustar_inds,ind,:),1,'omitnan'))';
    plot(k_mean,percentdiff_mean_slice,'k-','LineWidth',3);
    plot(k_mean,percentdiff_mean_slice,'Color',bin_cmap(i,:),'LineWidth',2);
end
hold off
box on
xlim(kstar_lims)
ylim([-100 20])
xlabel('k*')
ylabel('percent difference [%]')
ax_struc(2).ax = gca;
colormap(ax_struc(2).ax,bin_cmap)

cbar = colorbar;
cbar.Location="eastoutside";
cbar.Ticks = 0:0.1:0.5;
clim([0 0.5])
set(get(cbar,'label'),'string','u_* [m s^{-1}]');
text(text_x,text_y,'(b)','FontSize',labelsize,'HorizontalAlignment','center','Units','normalized')

tile_cleaner(ax_struc,tlayout)

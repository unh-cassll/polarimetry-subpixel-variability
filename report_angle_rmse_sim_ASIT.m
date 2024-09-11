%
function report_angle_rmse_sim_ASIT(fignum)

cmap = flipud(spectral(256));

blocksize_lims = [1e-1 1e2]*2;
block_ticks = [1 10 100];
block_ticklabels = {'1','10','100'};
RMSE_lims = 10.^[-3 1];
RMSE_ticks = 10.^(log10(RMSE_lims(1)):1:log10(RMSE_lims(end)));
RMSE_ticklabels = {};
for n = 1:length(RMSE_ticks)
    num_decimal_points = -1*(floor(log10(RMSE_ticks(n))));
    if num_decimal_points < 0
        num_decimal_points = 0;
    end
    RMSE_ticklabels{n} = sprintf(['%0.' num2str(num_decimal_points) 'f'],RMSE_ticks(n));
end

mult_val = 2;

slope_lims = log10(mult_val*[1e-3 1e-1]);
cbar_ticks = log10(mult_val*[1e-3 3.16e-3 1e-2 3.16e-2 1e-1]);
slope_lims = log10([0.005 0.2]);
cbar_ticks = log10(logspace(slope_lims(1),slope_lims(end),5));
cbar_ticklabels = {};
for n = 1:length(cbar_ticks)
    cbar_ticklabels{n} = sprintf(['%0.' num2str(abs(floor(cbar_ticks(n)))+1) 'f'],10^cbar_ticks(n));
end

load('simulated_surface_comparison_stats_full_resolution.mat')
load('ASIT_Pyxis_stats.mat')

dx = 1.3e-3;
blocksize_vec = ASIT_Stats.blockscale*dx;
rms_slope = ASIT_Stats.rms_slope;
RMSE_ori = ASIT_Stats.RMSE_ori;
RMSE_aoi = ASIT_Stats.RMSE_aoi;

% block_vec(block_vec*full_resolution_m_px<2.6e-3) = NaN;

figure(fignum);clf
tlayout = tiledlayout(2,2);

ax_struc = struct();

h(1) = nexttile(1);
hold on
plot(1000*block_vec*full_resolution_m_px,squeeze(RMSE_mat(6,:,:)),'o','markerfacecolor','k','markeredgecolor','k','markersize',9)
for n = 1:12
    scatter(1000*block_vec*full_resolution_m_px,reshape(RMSE_mat(6,n,:),[],1),60,log10(reshape(SUB_BLOCK_MSS_mat(n,:).^0.5,[],1)),'filled')
end
f = fill([0.22 0.9 0.9 0.22],[5 5 9 9],'w');
f.FaceAlpha = 0.5;
hold off
box on
colormap(cmap)
clim(slope_lims)
ax_struc(1).ax=gca;
ax_struc(1).ax.XTick = block_ticks;
ax_struc(1).ax.XTickLabels = block_ticklabels;
ax_struc(1).ax.YTick = RMSE_ticks;
ax_struc(1).ax.YTickLabels = RMSE_ticklabels;
ax_struc(1).ax.XScale='log';
ax_struc(1).ax.YScale='log';
ylim(RMSE_lims)
xlim(blocksize_lims)
title('\theta')
xlabel('block size [mm]')
ylabel('rmse [\circ]')
text(0.44,7,'modeled','FontSize',16,'HorizontalAlignment','center')

h(2) = nexttile(2);
hold on
plot(1000*block_vec*full_resolution_m_px,squeeze(RMSE_mat(5,:,:)),'o','markerfacecolor','k','markeredgecolor','k','markersize',9)
for n = 1:12
    scatter(1000*block_vec*full_resolution_m_px,reshape(RMSE_mat(5,n,:),[],1),60,log10(reshape(SUB_BLOCK_MSS_mat(n,:).^0.5,[],1)),'filled')
end
ah = annotation('arrow','Color','k','HeadStyle','vback2');
ah.Position = [5.52e-1 7.6e-1 0 6e-2];
ah.LineWidth = 2.5;
text(4.05e-1,8.5e-2,'U_{10}','FontSize',20,'HorizontalAlignment','center')
hold off
box on
colormap(cmap)
clim(slope_lims)
ax_struc(2).ax=gca;
ax_struc(2).ax.XTick = block_ticks;
ax_struc(2).ax.XTickLabels = block_ticklabels;
ax_struc(2).ax.XScale='log';
ax_struc(2).ax.YScale='log';
ylim(RMSE_lims)
xlim(blocksize_lims)
title(char(966))
xlabel('block size [mm]')
ylabel('RMSE')

numbins = 1;
quantiles = [10 25 50 75 90];

rms_slope_binned = NaN*ones(length(quantiles),size(RMSE_aoi,2));
RMSE_aoi_binned = rms_slope_binned;
RMSE_ori_binned = rms_slope_binned;

for n = 1:size(blocksize_vec,2)

    slice_rms_slope = rms_slope(:,n);
    slice_aoi = RMSE_aoi(:,n);
    slice_ori = RMSE_ori(:,n);

    [rms_slope_quantiles,aoi_quantiles,~] = compute_quantiles_fixed_binsize(slice_rms_slope,slice_aoi,numbins,quantiles);
    [~,ori_quantiles,~] = compute_quantiles_fixed_binsize(slice_rms_slope,slice_ori,numbins,quantiles);

    rms_slope_binned(:,n) = rms_slope_quantiles;
    RMSE_aoi_binned(:,n) = aoi_quantiles;
    RMSE_ori_binned(:,n) = ori_quantiles;

end

msize = 13;

h(3) = nexttile(3);
hold on
plot(1000*blocksize_vec(1,:),RMSE_aoi_binned','s','markerfacecolor','k','markeredgecolor','k','markersize',msize)
for n = 1:length(quantiles)
    s = scatter(1000*blocksize_vec(1,:),RMSE_aoi_binned(n,:),0.9*msize^2,log10(rms_slope_binned(n,:)),'filled');
    s.Marker = 's';
end
f = fill([0.22 0.9 0.9 0.22],[5 5 9 9],'w');
f.FaceAlpha = 0.5;
hold off
box on
colormap(cmap)
clim(slope_lims)
ax_struc(3).ax=gca;
ax_struc(3).ax.XTick = block_ticks;
ax_struc(3).ax.XTickLabels = block_ticklabels;
ax_struc(3).ax.YTick = RMSE_ticks;
ax_struc(3).ax.YTickLabels = RMSE_ticklabels;
ax_struc(3).ax.XScale='log';
ax_struc(3).ax.YScale='log';
ylim(RMSE_lims)
xlim(blocksize_lims)
xlabel('block size [mm]')
ylabel('rmse [\circ]')
text(0.44,7,'observed','FontSize',16,'HorizontalAlignment','center')
for m = 1:length(quantiles)
    text(1.1,RMSE_aoi_binned(m,1),[sprintf('%0.0f',quantiles(m)) '^{th}'],'FontSize',16)
end

h(4) = nexttile(4);
hold on
plot(1000*blocksize_vec(1,:),RMSE_ori_binned','s','markerfacecolor','k','markeredgecolor','k','markersize',msize)
for n = 1:length(quantiles)
    s = scatter(1000*blocksize_vec(1,:),RMSE_ori_binned(n,:),0.9*msize^2,log10(rms_slope_binned(n,:)),'filled');
    s.Marker = 's';
end
hold off
box on
colormap(cmap)
clim(slope_lims)
ax_struc(4).ax=gca;
ax_struc(4).ax.XTick = block_ticks;
ax_struc(4).ax.XTickLabels = block_ticklabels;
ax_struc(4).ax.XScale='log';
ax_struc(4).ax.YScale='log';
ylim(RMSE_lims)
xlim(blocksize_lims)
xlabel('block size [mm]')
ylabel('RMSE')

cbar = colorbar(h(end));
clim(slope_lims)
cbar.Ticks = cbar_ticks;
cbar.TickLabels = cbar_ticklabels;
set(get(cbar,'Label'),'String','sub-block rms slope')
cbar.Layout.Tile = 'east';

tile_cleaner(ax_struc,tlayout)

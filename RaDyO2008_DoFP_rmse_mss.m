%
function RaDyO2008_DoFP_rmse_mss(fignum)

load('RaDyO2008_reprocess_stats.mat')
load('DoFP_intensity_avg_RMSE_struc.mat')

% Define colormap and panel label cell array
cmap = flipud(spectral(256));

text_x = 0.05;
text_y = 0.94;

% Define variable ticks/ticklabels
blocksize_ticks = [0.1 1 10 100];
blocksize_ticklabels = {'0.1','1','10','100'};

rmse_ticks = [0.1 1 10 100];
rmse_ticklabels = {'0.1','1','10','100'};

mss_percent_error_ticks = -40:10:40;

rms_mult = 3.6;
rms_slope_ticks = log10(logspace(log10(rms_mult*1e-2),log10(rms_mult/2*1e-1),5));
slope_lims = log10([0.02 0.2]);
rms_slope_ticks = log10(logspace(slope_lims(1),slope_lims(end),5));
rms_slope_ticklabels = {};
for n = 1:length(rms_slope_ticks)
    num_decimal_points = -1*(floor(rms_slope_ticks(n)));
    if num_decimal_points < 0
        num_decimal_points = 0;
    end
    rms_slope_ticklabels{n} = sprintf(['%0.' num2str(num_decimal_points+1) 'f'],10^rms_slope_ticks(n));
end

blocksize_vec = [1:1:10 12:2:20 24:4:40];

num_blocks = find(blocksize_vec==20,1,'first');

% Preallocate arrays of interest
ang_rmse = NaN*ones(length(running_struc),num_blocks);
mss_avg = ang_rmse;
mss_dofp = ang_rmse;
mss_tail = ang_rmse;

% Grab mean square slope (total and tail) and RMSE angle
for n = 1:length(running_struc)

    if n ~= 10 && n ~= 13

        for m = 1:num_blocks

            block_dofp = ['block' leading_zeroes(blocksize_vec(m),2)];
            dofp_ind = find(blocksize_vec==blocksize_vec(m)*2);
            block_avg = ['block' leading_zeroes(blocksize_vec(dofp_ind),2)];

            mss_dofp(n,m) = trapz(running_struc(n).(block_dofp).dofp_2x2.k,running_struc(n).(block_dofp).dofp_2x2.Sk);
            mss_avg(n,m) = trapz(running_struc(n).(block_avg).dofp_none.k,running_struc(n).(block_avg).dofp_none.Sk);

            ind_cut = length(running_struc(n).(block_dofp).dofp_2x2.k)+1;
            k_full = running_struc(n).block01.dofp_none.block_averaged_slopefield_spectra(1).k;
            Sk_full = running_struc(n).block01.dofp_none.block_averaged_slopefield_spectra(1).Sk;
            mss_tail(n,m) = trapz(k_full(ind_cut:end),Sk_full(ind_cut:end));

        end

            ang_rmse(n,:) = sum(RMSE_struc(n).ang,2);

    end

end

% Prepare arrays for scatterplots
ang_rmse_plot = reshape(ang_rmse,[],1);
mss_tail_plot = reshape(mss_tail,[],1);
mss_dofp_plot = reshape(mss_dofp,[],1);
mss_avg_plot = reshape(mss_avg,[],1);
blocksize_plot = reshape(repmat(blocksize_vec(1:num_blocks),length(running_struc),1),[],1);

msize = 7;
figure(fignum);clf
tlayout = tiledlayout(2,1);

ax_struc = struct();

% First panel: RMSE angle between DoFP imitation and intensity
% block-average
nexttile(1)
hold on
plot([0 60],[0 0],'--','Color',[0.6 0 0],'linewidth',2)
plot(blocksize_plot*2*1.3,ang_rmse_plot,'o','markerfacecolor','k','markeredgecolor','k','markersize',msize)
scatter(blocksize_plot*2*1.3,ang_rmse_plot,0.8*msize^2,log10(mss_tail_plot.^0.5),'filled')
hold off
colormap(cmap)
cbar=colorbar;
xlim([2e-1 2e2])
ylim(3*[1e-1 1e1])
clim([rms_slope_ticks(1) rms_slope_ticks(end)])
cbar.Ticks = rms_slope_ticks;
cbar.TickLabels = rms_slope_ticklabels;
box on
ax_struc(1).ax = gca;
ax_struc(1).ax.YScale = 'log';
ax_struc(1).ax.XScale = 'log';
ax_struc(1).ax.XTick = blocksize_ticks;
ax_struc(1).ax.XTickLabel = blocksize_ticklabels;
ax_struc(1).ax.YTick = rmse_ticks;
ax_struc(1).ax.YTickLabel = rmse_ticklabels;
xlabel('blocksize [mm]')
ylabel('rmse: DoFP vs. mean intensity [\circ]')
set(get(cbar,'Label'),'String','sub-block rms slope')
cbar.Location = 'northoutside';
text(text_x,text_y,'(a)','FontSize',20,'Units','normalized','HorizontalAlignment','center')

% Second panel: percent difference in mss between DofP imitation and
% intensity block-average
nexttile(2)
hold on
plot([0.2 200],[0 0],'k--','linewidth',1.5)
plot(blocksize_plot*2*1.3,100*(mss_dofp_plot./mss_avg_plot-1),'o','markerfacecolor','k','markeredgecolor','k','markersize',msize)
scatter(blocksize_plot*2*1.3,100*(mss_dofp_plot./mss_avg_plot-1),0.8*msize^2,log10(mss_tail_plot.^0.5),'filled')
hold off
colormap(cmap)
xlim([0.2 200])
ylim([-1 1]*40)
clim([rms_slope_ticks(1) rms_slope_ticks(end)])
box on
xlabel('blocksize [mm]')
ylabel('mss % difference')
ax_struc(2).ax = gca;
ax_struc(2).ax.XScale = 'log';
ax_struc(2).ax.XTick = blocksize_ticks;
ax_struc(2).ax.XTickLabel = blocksize_ticklabels;
ax_struc(2).ax.YTick = mss_percent_error_ticks;
text(text_x,text_y,'(b)','FontSize',20,'Units','normalized','HorizontalAlignment','center')

tile_cleaner(ax_struc,tlayout)
% Compares spectra computed from DoAm polarimeter data collected during the 
% 2008 RaDyO campaign with spectra computed from a version of those data
% which have been degraded to imitate a DoFP polarimeter
%
% N. Laxague 2024
% 
function RaDyO2008_DoFP_sim(fignum,labelsize)

% Load in field data
s = load('RaDyO_2008_omnidirectional_slope_spectra.mat');
RaDyO_2008_omnidirectional_slope_spectra = s.RaDyO_2008_omnidirectional_slope_spectra;
k = RaDyO_2008_omnidirectional_slope_spectra(1).k_avg;
s = load('RaDyO_wind_speed_stress.mat');
ustar = s.ustar;
ustar = ustar(1:length(RaDyO_2008_omnidirectional_slope_spectra));
[~,order] = sort(ustar);

% Plot limits and colorbars
cmap = viridis(length(ustar));
ustar_lims = [5 40]/100;
wind_ind = floor(255*((ustar-min(ustar)))/(ustar_lims(2)-min(ustar)))+1;
cmap_full = viridis(256);
klims = [3e0 3e3];
plims = [-1 1]*60;
lw_thick = 3;
lw_thin = 2;

% Construct panel label cell array and plot positions
labels = {'(a)','(b)','(c)'};
text_x = 0.08;
text_y = 0.93;

% Generate figure
figure(fignum);clf
tlayout = tiledlayout(1,3);

% Lay out panel labels and horizontal reference line
for m = 1:3
    nexttile(m)
    text(text_x,text_y,labels{m},'FontSize',labelsize,'HorizontalAlignment','center','Units','normalized')
    hold on
    if ~mod(m,3)
        plot(klims,0*klims,'k--','linewidth',2)
    end
end

% Plot black under-line
for ind = 1:length(RaDyO_2008_omnidirectional_slope_spectra)

    n = order(ind);

    if n ~= 10 && n ~=13

        k_avg = RaDyO_2008_omnidirectional_slope_spectra(n).k_avg;
        Sk_avg = RaDyO_2008_omnidirectional_slope_spectra(n).Sk_avg;

        nexttile(1)
        plot(k_avg,k_avg.*Sk_avg,'k-','linewidth',lw_thick)
        ax_struc(1).ax = gca;

        k_dofp = RaDyO_2008_omnidirectional_slope_spectra(n).k_dofp;
        Sk_dofp = RaDyO_2008_omnidirectional_slope_spectra(n).Sk_dofp;

        nexttile(2)
        plot(k_dofp,k_dofp.*Sk_dofp,'k-','linewidth',lw_thick)
        ax_struc(2).ax = gca;

        nexttile(3)
        plot(k_avg,100*(Sk_dofp-Sk_avg)./Sk_avg,'k-','linewidth',lw_thick)
        ax_struc(3).ax = gca;

    end

end

% Plot curves colored by wind forcing
for ind = 1:length(RaDyO_2008_omnidirectional_slope_spectra)

    n = order(ind);

    if n~=10 && n~=13

        k_avg = RaDyO_2008_omnidirectional_slope_spectra(n).k_avg;
        Sk_avg = RaDyO_2008_omnidirectional_slope_spectra(n).Sk_avg;

        nexttile(1)
        plot(k_avg,k_avg.*Sk_avg,'linewidth',lw_thin,'Color',cmap_full(wind_ind(n),:))
        ax_struc(1).ax = gca;

        k_dofp = RaDyO_2008_omnidirectional_slope_spectra(n).k_dofp;
        Sk_dofp = RaDyO_2008_omnidirectional_slope_spectra(n).Sk_dofp;

        nexttile(2)
        plot(k_dofp,k_dofp.*Sk_dofp,'linewidth',lw_thin,'Color',cmap_full(wind_ind(n),:))
        ax_struc(2).ax = gca;

        nexttile(3)
        plot(k_avg,100*(Sk_dofp-Sk_avg)./Sk_avg,'linewidth',lw_thin,'Color',cmap_full(wind_ind(n),:))
        ax_struc(3).ax = gca;

    end

end

hold off

% Enforce axis scaling, limits, ticks, etc. for spectra panels
B_inds = [1 2];
for ind = 1:length(B_inds)

    n = B_inds(ind);

    nexttile(n)
    box on
    ax_struc(n).ax.XScale = 'log';
    ax_struc(n).ax.YScale = 'log';
    ax_struc(n).ax.XTick = 10.^(0:3);
    ax_struc(n).ax.YTick = 10.^(-5:-1);
    xlim(klims)
    ylim([1e-4 1e-1])

    if n == 1 || n == 4
        ylabel('B(k) [rad]')
        cbar = colorbar;
        clim(ustar_lims*100)
        cbar.Ticks = 5:5:40;
        cbar.Location = 'south';
        set(get(cbar,'Label'),'String','u_* [cm s^{-1}]')

    else
        ax_struc(n).ax.YTickLabel = '';
    end
    xlabel('k [rad m^{-1}]')

end

% Enforce axis scaling, limits, ticks, etc. for ratio panel
ratio_inds = 3;
for ind = 1:length(ratio_inds)

    n = ratio_inds(ind);

    nexttile(n)
    box on
    colororder(cmap)

    ax_struc(n).ax.XScale = 'log';
    ax_struc(n).ax.XTick = 10.^(0:3);
    ax_struc(n).ax.YAxisLocation = 'right';

    xlim(klims)
    ylim(plims)
    xlabel('k [rad m^{-1}]')

    if mod(n,3)
        ax_struc(n).ax.YScale = 'log';
    end
    ylabel('% difference')
    hold on
    plot(max(k)*[1 1]*1,plims,'r-','linewidth',2)
    plot(max(k)*[1 1]*0.5,plims,'r-.','linewidth',2)
    plot(max(k)*[1 1]*0.25,plims,'r:','linewidth',2)
    hold off
    if n == 3
        text(max(k),plims(2)*1.12,'k_{Nyq}','FontSize',18,'Color','r','HorizontalAlignment','center')
        text(max(k)/2,plims(2)*1.12,'$\frac{1}{2}$','Interpreter','LaTeX','FontSize',20,'Color','r','HorizontalAlignment','center')
        text(max(k)/4,plims(2)*1.12,'$\frac{1}{4}$','Interpreter','LaTeX','FontSize',20,'Color','r','HorizontalAlignment','center')
    end

end

tlayout.TileSpacing = 'compact';

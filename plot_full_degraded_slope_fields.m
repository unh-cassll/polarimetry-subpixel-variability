% Loads in raw DoAm polarization data and degrades it through
% * block averaging of intensity or slope
% * mimicry of DoFP array
%
% N. Laxague 2024
%
function plot_full_degraded_slope_fields(fignum,labelsize)

% Load in raw data and set up full-size array
load('RaDyO_2008_example_polarized_intensities.mat')

dx_cm = 100*m_per_px;
N_full = size(I0,1);
x_full = 0:dx_cm:dx_cm*(N_full-1);

dofp_sim_type = 'none';
blocksize = 8;

x_avg = 0:dx_cm*blocksize:dx_cm*N_full;

% Compute along-look slope component, full resolution
[~,Sy] = compute_slope_field_from_polarization_intensities(I0,I45,I90,I135,dofp_sim_type,1);
Sy = Sy - mean(Sy,'all','omitnan');

slopefield_holder_struc(1).Sy = Sy;

% Compute along-look slope component, averaged at level of slope
[Sy_avg_slopefield] = block_average_and_subsample(Sy,blocksize*[1 1]);

slopefield_holder_struc(2).Sy = Sy_avg_slopefield;

% Compute along-look slope component, averaged at level of intensity
[~,Sy_avg_intensities] = compute_slope_field_from_polarization_intensities(I0,I45,I90,I135,dofp_sim_type,blocksize);
Sy_avg_intensities = Sy_avg_intensities - mean(Sy_avg_intensities,'all','omitnan');

slopefield_holder_struc(3).Sy = Sy_avg_intensities;

% Compute along-look slope component, imitation of DoFP array & bilinear
% interpolation
[~,Sy_DoFP_2x2] = compute_slope_field_from_polarization_intensities(I0,I45,I90,I135,'2x2',blocksize/4);
Sy_DoFP_2x2 = Sy_DoFP_2x2 - mean(Sy_DoFP_2x2,'all','omitnan');

Sy_DoFP_2x2 = Sy_DoFP_2x2(1:2:end,1:2:end);

slopefield_holder_struc(4).Sy = Sy_DoFP_2x2;

% Compute along-look slope component, imitation of DoFP array & 12-pixel
% kernel approach of Ratliff et al. (2009)
[~,Sy_DoFP_4x4] = compute_slope_field_from_polarization_intensities(I0,I45,I90,I135,'4x4',blocksize/4);
Sy_DoFP_4x4 = Sy_DoFP_4x4 - mean(Sy_DoFP_4x4,'all','omitnan');

slopefield_holder_struc(5).Sy = Sy_DoFP_4x4;


X_center = dx_cm*N_full*0.9;
Y_center = X_center*1.02;
DX = dx_cm*N_full*0.08;
DY = dx_cm*N_full*0.06;

xticks = 0:10:60;
cticks = -0.5:0.1:0.5;

label_cell = {'(a)','(b)','(c)','(d)','(e)'};

position_vec = [1 3 4 7 8];
span_vec = [2 1 1 1 1];

% Loop to create multi-panel figure
figure(fignum);clf
tlayout = tiledlayout(2,4);

ax_struc = struct();

for n = 1:5

    if n == 1
        x_plot = x_full;
        X = [-1 1 1 -1]*DX*0.5 + X_center*1.05;
        Y = [-1 -1 1 1]*DY*0.5 + Y_center*1.04;
    else
        x_plot = x_avg;
        X = [-1 1 1 -1]*DX + X_center;
        Y = [-1 -1 1 1]*DY + Y_center;
    end

    nexttile(position_vec(n),[1 1]*span_vec(n))
    hold on
    imagesc(x_plot,x_plot,slopefield_holder_struc(n).Sy);colormap('gray');
    f = fill(X,Y,'w');
    hold off
    box on
    xlim([0 dx_cm*N_full])
    ylim([0 dx_cm*N_full])
    clim([-1 1]*0.5)
    pbaspect([1 1 1])
    ax_struc(n).ax = gca;
    ax_struc(n).ax.XTick = xticks;
    ax_struc(n).ax.YTick = xticks;
    grid off
    f.FaceAlpha = 0.7;
    if n == 1
        text(X_center*1.05,Y_center*1.04,label_cell{n},'FontSize',labelsize,'HorizontalAlignment','center')
    else
        text(X_center,Y_center,label_cell{n},'FontSize',labelsize,'HorizontalAlignment','center')
    end
    

    if n == 1
        cbar = colorbar;
        set(get(cbar,'Label'),'String','S_y')
        cbar.Location = 'northoutside';
        cbar.Ticks = cticks;
        xlabel('x [cm]')
        ylabel('y [cm]')
    else
        ax_struc(n).ax.XTickLabel = '';
        ax_struc(n).ax.YTickLabel = '';
    end

end


tlayout.TileSpacing = 'compact';

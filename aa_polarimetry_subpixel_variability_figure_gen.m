% Figure generation script for
% "The Effects of Subpixel Variability on Polarimetric Sensing of Ocean Waves"
% 
% N. J. M. Laxague and co-authors, 2024
%

addpath _codes
addpath _codes/ElfouhailyEtAl1997
addpath _data
addpath _outputs

figure_style

close all;clear;clc

corner_x = 50;
corner_y = 50;

full_width = 1500;
full_height = 600;

fsize = 22;

n = 1;
print_options = {'none','svg','png'};
print_option = print_options{n};

dpi_val = 100;
dpi_string = ['-r' num2str(dpi_val)];

figpos = [corner_x corner_y full_width full_height];

out_path = '_outputs/manuscript_figures/';


%% Example light reflection from gravity-capillary wave

fignum = 1;
synthetic_slice_DoLP_example(fignum)
set(fignum,'Position',figpos.*[1 1 0.5 1])

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'MelvilleFedorov2015_surface_DoLP_example.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'MelvilleFedorov2015_surface_DoLP_example.png'],'-dpng',dpi_string);
end

%% Draw microgrid polarizers

fignum = 2;

draw_superpixel_elements(fignum)
set(fignum,'Position',figpos.*[1 1 1 1*0.8])

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'DoFP_superpixel_and_interpolation.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'DoFP_superpixel_and_interpolation.png'],'-dpng',dpi_string);
end

%% Example slope fields: full resolution and degraded

fignum = 3;

plot_full_degraded_slope_fields(fignum)
set(fignum,'Position',figpos.*[1 1 1 1.4])

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'slopefield_full_degraded_examples.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'slopefield_full_degraded_examples.png'],'-dpng',dpi_string);
end

%% Synthetic surface simulation

fignum = 4;

synthetic_surface_DoLP_ORI(fignum)
set(fignum,'Position',figpos)

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'synthetic_surface_DoLP_ORI.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'synthetic_surface_DoLP_ORI.png'],'-dpng',dpi_string);
end

%% Light reflection calculation

fignum = 5;

polarized_ray_tracing(fignum,fsize)
set(fignum,'Position',figpos.*[1 1 0.5 1])

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'light_reflection_demonstration.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'light_reflection_demonstration.png'],'-dpng',dpi_string);
end

%% Flowchart

% fignum = 6;

%% Surface stats: synthetic and ASIT

fignum = 7;

report_angle_rmse_sim_ASIT(fignum)
set(fignum,'Position',figpos.*[1 1 1 2])

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'simulated_and_ASIT_surface_angle_RMSE.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'simulated_and_ASIT_surface_angle_RMSE.png'],'-dpng',dpi_string);
end

%% Spectral nudge

fignum = 8;

simulated_surface_stats_nudged_Elfouhaily(fignum)

set(fignum,'Position',figpos.*[1 1 0.5 1.2])
set(fignum*10,'Position',figpos.*[1 1 0.5 1.2])
set(fignum*100,'Position',figpos.*[1 1 0.5 1.2])

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);pause(0.5);print('_outputs/manuscript_figures/spectral_nudge_saturation_spect.svg','-dsvg')
        figure(fignum*10);pause(0.5);print('_outputs/manuscript_figures/spectral_nudge_RMSE.svg','-dsvg')
        figure(fignum*100);pause(0.5);print('_outputs/manuscript_figures/spectral_nudge_lambda_p.svg','-dsvg')
    case 'png'
        figure(fignum);pause(0.5);print('_outputs/manuscript_figures/spectral_nudge_saturation_spect.png','-dpng',dpi_string)
        figure(fignum*10);pause(0.5);print('_outputs/manuscript_figures/spectral_nudge_RMSE.png','-dpng',dpi_string)
        figure(fignum*100);pause(0.5);print('_outputs/manuscript_figures/spectral_nudge_lambda_p.png','-dpng',dpi_string)
end

%% Effect of block-averaging on spectra, ASIT data

fignum = 9;

plot_block_averaged_spectra_comparison(fignum)
set(fignum,'Position',figpos.*[1 1 0.5 1.3])
set(fignum+1,'Position',figpos.*[1 1 0.5 2])

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'block_averaged_ASIT_spectra.svg'],'-dsvg')
        figure(fignum+1);print([out_path 'block_averaged_ASIT_spectra_normalized.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'block_averaged_ASIT_spectra.png'],'-dpng',dpi_string);
        figure(fignum+1);print([out_path 'block_averaged_ASIT_spectra_normalized.png'],'-dpng',dpi_string);
end


%% Degradation of RaDyO data: spectra

fignum = 11;

RaDyO2008_DoFP_sim(fignum)

set(fignum,'Position',figpos.*[1 1 1 1])

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'RaDyO2008_DoFP_sim_spectra_ratio.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'RaDyO2008_DoFP_sim_spectra_ratio.png'],'-dpng',dpi_string);
end

%% Degradation of RaDyO data: rmse angle

fignum = 12;

RaDyO2008_DoFP_rmse_mss(fignum)

set(fignum,'Position',figpos.*[1 0.1 0.5 2.2])

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'RaDyO2008_DoFP_rmse_mss.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'RaDyO2008_DoFP_rmse_mss.png'],'-dpng',dpi_string);
end

%% Subpixel RMS slope

fignum = 13;

subpixel_RMS_slope(fignum)

set(fignum,'Position',figpos.*[1 1 0.5 1.4])

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'subpixel_RMS_slope.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'subpixel_RMS_slope.png'],'-dpng',dpi_string);
end
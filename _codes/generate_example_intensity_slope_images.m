file_ind = 9;
t_on = tic;
in_dir = '../_data/RaDyO_2008_snippets/';
dir_info = dir([in_dir '*.nc']);
dofp_sim_strings = {'none','2x2'};

nc_name = [in_dir dir_info(file_ind).name];
info_var = ncinfo(nc_name);
lengths = [info_var.Dimensions.Length];
I0_stack = ncread(nc_name,'I0');
I45_stack = ncread(nc_name,'I45');
I90_stack = ncread(nc_name,'I90');
I135_stack = ncread(nc_name,'I135');


%%

% Load in data reduction matrix and DoLP:theta relationship
load('DRMSB.mat')
load('dolp_theta_vecs.mat')

% Only consider 0-90 degree portion of DoLP:theta relationship
ind = find(DOLP_full==max(DOLP_full),1,'first');
DOLP_unique = (0.0001:0.0001:1);
theta_unique = interp1(DOLP_full(1:ind),theta_full(1:ind),DOLP_unique);

[s1,~,~] = size(I0_stack);
if mod(s1,2)
    I0_stack = I0_stack(1:s1-1,1:s1-1,:);
    I45_stack = I45_stack(1:s1-1,1:s1-1,:);
    I90_stack = I90_stack(1:s1-1,1:s1-1,:);
    I135_stack = I135_stack(1:s1-1,1:s1-1,:);
end

% Compute average values of data reduction matrix
mults = NaN*ones(4,4);
for i = 1:4
    for j = 1:4
        mults(i,j) = mean(DRM(i,j,:,:),'all','omitnan');
    end
end

% Compute Stokes vector using data reduction matrix
S0 = mults(1,1)*double(I135_stack) + mults(1,2)*double(I45_stack) + mults(1,3)*double(I90_stack) + mults(1,4)*double(I0_stack);
S1 = mults(2,1)*double(I135_stack) + mults(2,2)*double(I45_stack) + mults(2,3)*double(I90_stack) + mults(2,4)*double(I0_stack);
S2 = mults(3,1)*double(I135_stack) + mults(3,2)*double(I45_stack) + mults(3,3)*double(I90_stack) + mults(3,4)*double(I0_stack);


% Compute DoLP and orientation angle
DOLP = sqrt(S1.^2+S2.^2)./S0;
ORI = 0.5*atan2(S2,S1)*180/pi;

clear S1 S2

% Estimate angle of incidence from DoLP
DOLP_int = floor(DOLP*10000);

clear DOLP

DOLP_int(DOLP_int<1) = 1;
DOLP_int(DOLP_int>10000) = 10000;
AOI = theta_unique(DOLP_int);


%%

i = 1000;

blocksize = 8;

in_array = I0_stack(:,:,i);

[out_array] = block_average_and_subsample(squeeze(in_array),[1 1]*blocksize);

figure(1);clf
set(gcf,'Position',[100 100 800 800])
imagesc(log10(double(in_array)));colormap('gray');grid off;box on;pbaspect([1 1 1])
ax = gca;
ax.XTick = [];
ax.YTick = [];

figure(2);clf
set(gcf,'Position',[100 100 800 800])
imagesc(log10(out_array));colormap('gray');grid off;box on;pbaspect([1 1 1])
ax = gca;
ax.XTick = [];
ax.YTick = [];

figure(3);clf
set(gcf,'Position',[100 100 800 800])
imagesc(AOI(:,:,i));colormap('gray');grid off;box on;pbaspect([1 1 1])
ax = gca;
ax.XTick = [];
ax.YTick = [];

cmap = gray(256);

obs_intensity_full = generate_im(log10(double(in_array)),cmap,min(log10(double(in_array)),[],'all'),max(log10(double(in_array)),[],'all'));
imwrite(obs_intensity_full,'../_outputs/manuscript_figures/obs_intensity_full.tiff')

obs_intensity_reduced = generate_im(log10(double(out_array)),cmap,min(log10(double(out_array)),[],'all'),max(log10(double(out_array)),[],'all'));
imwrite(obs_intensity_reduced,'../_outputs/manuscript_figures/obs_intensity_reduced.tiff')

aoi_slice = squeeze(AOI(:,:,i));
aoi_slice = inpaint_nans(aoi_slice);

obs_theta_full = generate_im(aoi_slice,cmap,min(AOI(:,:,i),[],'all'),max(AOI(:,:,i),[],'all'));
imwrite(obs_theta_full,'../_outputs/manuscript_figures/obs_theta_full.tiff')

%%

close all;clear;clc

load('example_wave_elevation_and_angles.mat')

inds = 1:4:2000;

wse_m = wse_m(inds,inds);
AOI = AOI(inds,inds);

[gx,gy] = gradient(wse_m);

blocksize = 8;

[out_array] = block_average_and_subsample(gy,[1 1]*blocksize);

figure(1);clf
set(gcf,'Position',[100 100 800 800])
imagesc(gy);colormap('gray');grid off;box on;pbaspect([1 1 1])
ax = gca;
ax.XTick = [];
ax.YTick = [];

figure(2);clf
set(gcf,'Position',[100 100 800 800])
imagesc(out_array);colormap('gray');grid off;box on;pbaspect([1 1 1])
ax = gca;
ax.XTick = [];
ax.YTick = [];

figure(3);clf
set(gcf,'Position',[100 100 800 800])
imagesc(AOI);colormap('gray');grid off;box on;pbaspect([1 1 1])
ax = gca;
ax.XTick = [];
ax.YTick = [];

cmap = gray(256);

sim_slope_full = generate_im(gy,cmap,min(gy,[],'all'),max(gy,[],'all'));
imwrite(sim_slope_full,'../_outputs/manuscript_figures/sim_slope_full.tiff')

sim_slope_reduced = generate_im(out_array,cmap,min(out_array,[],'all'),max(out_array,[],'all'));
imwrite(sim_slope_reduced,'../_outputs/manuscript_figures/sim_slope_reduced.tiff')

sim_theta_full = generate_im(AOI,cmap,min(AOI,[],'all'),max(AOI,[],'all'));
imwrite(sim_theta_full,'../_outputs/manuscript_figures/sim_theta_full.tiff')
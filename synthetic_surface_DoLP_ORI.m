%
function synthetic_surface_DoLP_ORI(fignum)

% grid_size = 2^13;
% 
% % Set up dimensions to loop over: wind speed, frames, blocks
% U10 = 5;
% nframes = 5;
% fps = 0.1;
% 
% % Set option that ignores sub-pixel mss variability
% subpixel_slope_var = 0;
% 
% dx = 5e-4;
% full_size_m = grid_size*dx;

% Comparison spectrum: 2 km to 0.1 mm wavelengths
% min_k = pi*1e-3;
% max_k = pi/dx;
% num_k = 1000;
% fetch_m = 1e6;

% % Compute total mean square slope and Hs for given wind speed
% [k,Fk] = Elfouhaily_omni(U10,min_k,max_k,num_k,fetch_m);
% mss = trapz(k,k.^2.*Fk);
% Hs = 4*sqrt(trapz(k,Fk));
% 
% % Simulate array of water surface elevation fields
% [x,y,eta] = ocean_simulator_compact_elfouhaily(U10,grid_size,dx,nframes,fps);
% 
% frame_ind = 5;
% 
% % Pluck out slice of simulated sea surface, model reflected Stokes
% wse_m = squeeze(eta(:,:,frame_ind));
% [S0,S1,S2] = model_reflected_Stokes_vector_compact(wse_m,dx,subpixel_slope_var);

% Compute slope field

% [gx,gy] = gradient(wse_m);
% 
% sx = gx/dx;
% sy = gy/dx;

% load('example_wave_elevation_and_angles.mat')
load('example_wave_angles.mat')
load('example_wave_elevation.mat')
load('example_wave_scale_info.mat')

wse_m = wse_m - mean(wse_m,'all');

% Make plots

grid_size = 2000;
num_ticks = 5;
tick_increment = floor(grid_size*dx/num_ticks*10)/10;

xlims = [0 max(x)];
ylims = [0 max(y)];
% eta_lims = [-1 1]*(ceil(100*2*max(abs(wse_m),[],'all'))-1);
eta_lims = [-1 1]*4;
eta_diff = abs(diff(eta_lims));
if eta_diff == 0
    eta_lims = [-1 1];
    eta_diff = 2;
end
eta_increment = eta_diff/4;
eta_ticks = eta_lims(1):eta_increment:eta_lims(2);

S_lims = [-1 1]*0.5;
DoLP_lims = [0 1];
AOI_lims = [-1 1]*30;
ORI_lims = [-1 1]*90;

AOI_ticks = -30:15:30;

figticks = 0:tick_increment:num_ticks*tick_increment;

max_array_size = 256;
subnum = floor(length(x)/max_array_size);

figure(fignum);clf
tlayout = tiledlayout(1,3);

nexttile()
% surf(x(1:subnum:end),y(1:subnum:end),100*wse_m(1:subnum:end,1:subnum:end),log10(S0(1:subnum:end,1:subnum:end)));colormap('gray');shading('flat')
% pbaspect([max(x) max(y) 2*max(abs(wse_m),[],'all')].*[1 1 10])
% xlim(xlims)
% ylim(ylims)
% zlim([-1 1]*100*2*max(abs(wse_m),[],'all'))
% % title('S_0')
% % xlabel('x [m]')
% % ylabel('y [m]')
% % zlabel('z [cm]')
% set(gcf,'Position',figpos)%.*[1 1 1.25 1])
% view([45 30])
% ax = gca;
% ax.XTick = figticks;
% ax.YTick = figticks;
% ax.XTickLabel = '';
% ax.YTickLabel = '';
% ax.ZTickLabel = '';
imagesc(x,y,100*wse_m);colormap('gray');shading('flat');view([0 -90])
pbaspect([1 1 1])
xlim(xlims)
ylim(ylims)
clim(eta_lims)
xlabel('x [m]')
ylabel('y [m]')
ax_struc(1).ax = gca;
ax_struc(1).ax.XTick = figticks;
ax_struc(1).ax.YTick = figticks;
grid('off')
cbar = colorbar;
cbar.Ticks = eta_ticks;
cbar.Location = 'northoutside';
set(get(cbar,'Title'),'String','\eta [cm]')

% figure(1);clf
% imagesc(x,y,sx);colormap('gray');colorbar;view([0 -90])
% pbaspect([1 1 1])
% xlim(xlims)
% ylim(ylims)
% clim(S_lims)
% title('S_x')
% xlabel('x [m]')
% ylabel('y [m]')
% set(gcf,'Position',figpos)
% ax = gca;
% ax.XTick = figticks;
% ax.YTick = figticks;
% grid('off')
% 
% figure(2);clf
% imagesc(x,y,sy);colormap('gray');colorbar;view([0 -90])
% pbaspect([1 1 1])
% xlim(xlims)
% ylim(ylims)
% clim(S_lims)
% title('S_y')
% xlabel('x [m]')
% ylabel('y [m]')
% set(gcf,'Position',figpos)
% ax = gca;
% ax.XTick = figticks;
% ax.YTick = figticks;
% grid('off')

% DoLP = sqrt(S1.^2+S2.^2)./S0;
% ORI = 0.5*atan2(S2,S1)*180/pi;
% ORI(ORI<0) = ORI(ORI<0) + 180;
% ORI = ORI - 90;
% inds = 4801:6800;
% [phase_unwrap,~]=Unwrap_TIE_DCT_Iter(ORI(inds,inds)*pi/180*4);
% ORI = phase_unwrap*180/pi/4;
% x = x(inds) - 2.4;
% y = y(inds) - 2.4;
% DoLP = DoLP(inds,inds);
% wse_m = wse_m(inds,inds);
% 
% s = load('dolp_theta_vecs_low.mat');
% DOLP_vec = s.DOLP_vec;
% theta_vec = s.theta_vec;
% 
% % Estimate angle of incidence from DoLP
% DOLP_int = floor(DoLP*10000);
% DOLP_int(DOLP_int<1) = 1;
% DOLP_int(DOLP_int>10000) = 10000;
% AOI = theta_vec(DOLP_int);

nexttile()
imagesc(x,y,AOI-35);colormap('gray');shading('flat');view([0 -90])
pbaspect([1 1 1])
xlim(xlims)
ylim(ylims)
clim(AOI_lims)
% title('DoLP')
xlabel('x [m]')
ylabel('y [m]')
ax_struc(2).ax = gca;
ax_struc(2).ax.XTick = figticks;
ax_struc(2).ax.YTick = figticks;
grid('off')
cbar = colorbar;
cbar.Ticks = AOI_ticks;
cbar.Location = 'northoutside';
set(get(cbar,'Title'),'String','\theta-\theta_i [\circ]')

nexttile()
imagesc(x,y,ORI);colormap('gray');colorbar;shading('flat');view([0 -90])
pbaspect([1 1 1])
xlim(xlims)
ylim(ylims)
clim(ORI_lims)
% title('orientation angle [\circ]')
xlabel('x [m]')
ylabel('y [m]')
ax_struc(3).ax = gca;
ax_struc(3).ax.XTick = figticks;
ax_struc(3).ax.YTick = figticks;
% ax.YAxisLocation = 'right';
grid('off')
cbar = colorbar;
cbar.Ticks = -90:45:90;
cbar.Location = 'northoutside';
% set(get(cbar,'Title'),'String','$\varphi [\circ]$','Interpreter','LaTeX')
set(get(cbar,'Title'),'String',[char(966) ' [\circ]'])

tile_cleaner(ax_struc,tlayout)

% %%
% 
% % Given inputs of
% % * incident Stokes vector
% % * water surface elevation field
% % * spatial resolution in meters
% % * azimuth angle of view in degrees
% % * incidence angle of view in degrees
% % * observation height in meters
% % * subpixel slope variance
% %
% % ... simulates the effect of the sea surface to modify the polarization state of reflected light
% 
% % Process following Mobley (2015)
% %
% % Code by N. Laxague 2024
% %
% function [S0,S1,S2] = model_reflected_Stokes_vector_compact(wse_m,dx,subpixel_slope_var)
% 
% % Virtual camera parameters
% incidence_view_deg = 35;
% azimuth_view_deg = 0;
% height_view_m = 100;
% S_inc = [1 0 0 0];
% 
% % Incident Stokes vector
% S_inc1 = S_inc(1);
% S_inc2 = S_inc(2);
% S_inc3 = S_inc(3);
% S_inc4 = S_inc(4);
% % camera_location = [sind(azimuth_view_deg+180)*tand(incidence_view_deg) cosd(azimuth_view_deg+180)*tand(incidence_view_deg) 1]'*height_view_m;
% camera_location = [0 -tand(incidence_view_deg)*height_view_m height_view_m]';
% 
% % Solar angles
% sun_polar_ang = 0;
% sun_az_ang = 0;
% 
% % Real index of refraction
% n = 1.34;
% 
% % Compute the surface normal field
% [gx,gy] = gradient(wse_m);
% 
% N_mat = NaN*repmat(wse_m,[1 1 3]);
% N_mat(:,:,1) = gx/dx;
% N_mat(:,:,2) = gy/dx;
% N_mat(:,:,3) = 0*gy + 1;
% 
% N_mat = N_mat./sqrt(sum(N_mat.^2,3));
% nx = N_mat(:,:,1);
% ny = N_mat(:,:,2);
% nz = N_mat(:,:,3);
% 
% % Run simulation
% [total_rows,total_cols] = size(wse_m);
% x = dx:dx:total_rows*dx;
% y = dx:dx:total_cols*dx;
% 
% % Produce subpixel slope distributions
% if subpixel_slope_var > 0
%     slope_var_x = subpixel_slope_var/2;
%     slope_var_y = subpixel_slope_var/2;
%     num_standard_dev = 3;
%     [Pxy,sx,sy] = produce_slope_distribution(slope_var_x,slope_var_y,num_standard_dev);
%     Pxy_copy = floor(Pxy);
%     Pxy_copy = Pxy_copy/min(Pxy_copy,[],'all');
%     slope_PDF_struc = struct();
%     increment = 0;
%     for i = 1:size(sx,1)
%         for j = 1:size(sx,2)
%             Pxy_count = Pxy_copy(i,j);
%             for counter = 1:Pxy_count
%                 increment = increment + 1;
%                 slope_PDF_struc(increment).sx = sx(i,j);
%                 slope_PDF_struc(increment).sy = sy(i,j);
%             end
%         end
%     end
%     sx_vec = [slope_PDF_struc.sx];
%     sy_vec = [slope_PDF_struc.sy];
%     num_rays_per_ensemble = 100;
%     rand_ind = rand(num_rays_per_ensemble,1);
%     rand_ind = floor((rand_ind-min(rand_ind))/(max(rand_ind)-min(rand_ind))*(length(sx_vec)-1))+1;
% else
%     sx_vec = 0;
%     sy_vec = 0;
%     num_rays_per_ensemble = 1;
%     rand_ind = 1;
% end
% 
% % Set up and run loop over all surface facets
% nx_vec = reshape(nx,[],1);
% ny_vec = reshape(ny,[],1);
% nz_vec = reshape(nz,[],1);
% X_vec = reshape(repmat(x,length(y),1),[],1);
% Y_vec = reshape(repmat(y',1,length(x)),[],1);
% Z_vec = reshape(wse_m,[],1);
% 
% total_inds = total_rows*total_cols;
% 
% 
% parfor ind = 1:total_inds
% 
%     % Grab surface normal vector and facet position
%     N0 = [nx_vec(ind) ny_vec(ind) nz_vec(ind)]';
%     X = X_vec(ind);
%     Y = Y_vec(ind);
%     Z = Z_vec(ind);
% 
%     % Compute reflected ray
%     ray_refl = camera_location - [X Y Z]';
%     ray_refl = ray_refl/norm(ray_refl);
% 
%     wiggle_x = sx_vec(rand_ind);
%     wiggle_y = sy_vec(rand_ind); %#ok<PFBNS>
%     wiggle_z = 0*wiggle_x;
% 
%     N = N0 + [wiggle_x;wiggle_y;wiggle_z];
%     N = N./sqrt(sum(N.^2));
% 
%     % Compute incident ray
%     ray_inc = ray_refl - (2*ray_refl'*N).*N;
% 
%     % Compute incidence angle with respect to surface facet
%     inc_dot_N = sum(ray_inc.*N,1);
%     theta_i = abs(acosd(inc_dot_N*-1));
%     theta_t = asind(sind(theta_i)./n);
% 
%     % Enforce sky-leaving radiance conditions
%     polar_ang = acosd([0 0 1]*N);
%     az_ang = atan2(N(1,:),N(2,:))*180/pi;
%     weight = (cosd(polar_ang-sun_polar_ang).*cosd(az_ang-sun_az_ang)).^2;
%     % weight = 0*N(1,:) + 1;
% 
%     % Ensure that ray is coming from above horizon
%     % NEED TO REPLACE THIS WITH MULTIPLE REFLECTION CALCULATION
%     pass_cond = ray_inc(3,:) < 0;
%     N = N(:,pass_cond);
%     ray_inc = ray_inc(:,pass_cond);
%     theta_i = theta_i(pass_cond);
%     theta_t = theta_t(pass_cond);
%     weight = weight(pass_cond);
% 
%     if ~isempty(N)
% 
%         % Compute h and s vectors (Mobley 2015)
%         % and from them, Stokes vector rotation angle alpha
%         h_inc = cross(repmat([0;0;1],1,length(theta_i)),ray_inc);
%         h_inc = h_inc./sqrt(sum(h_inc.^2));
%         s = cross(ray_inc,N);
%         s = s./sqrt(sum(s.^2));
%         h_cross_s = cross(h_inc,s);
% 
%         h_inc_dot_s = sum(h_inc.*s,1);
%         alpha_ang = abs(acosd(h_inc_dot_s));
%         antiparallel_inds = h_cross_s(3,:).*ray_inc(3,:) < 0;
%         alpha_ang(antiparallel_inds) = 360 - alpha_ang(antiparallel_inds);
%         c2a = cosd(2*alpha_ang);
%         s2a = sind(2*alpha_ang);
% 
%         % Rotation matrix, Stokes vector
%         % R = [1,0,            0,             0;
%         %      0,cosd(2*alpha),-sind(2*alpha),0;
%         %      0,sind(2*alpha),cosd(2*alpha), 0;
%         %      0,0,            0,             1];
% 
%         % Mueller matrix for reflection off surface, KA89
%         alpha = 1/2*(tand(theta_i-theta_t)./tand(theta_i+theta_t)).^2;
%         eta = 1/2*(sind(theta_i-theta_t)./sind(theta_i+theta_t)).^2;
%         gamma_Re = (tand(theta_i-theta_t).*sind(theta_i-theta_t))./(tand(theta_i+theta_t).*sind(theta_i+theta_t));
% 
%         % R_AM = [alpha+eta,alpha-eta,0,       0;
%         %         alpha-eta,alpha+eta,0,       0;
%         %         0,        0,        gamma_Re,0;
%         %         0,        0,        0,       gamma_Re];
% 
%         % Compute reflected Stokes vector
%         %S_refl = R*R_AM*S_inc;
%         S_refl_0 = S_inc1*(alpha+eta)+S_inc2*(alpha-eta);
%         S_refl_1 = S_inc1*c2a.*(alpha-eta)+S_inc2*c2a.*(alpha+eta)-S_inc3*s2a.*gamma_Re;
%         S_refl_2 = S_inc1*s2a.*(alpha-eta)+S_inc2*s2a.*(alpha+eta)+S_inc3*c2a.*gamma_Re;
%         S_refl_3 = S_inc4*gamma_Re;
% 
%         S_refl = median([S_refl_0; S_refl_1; S_refl_2; S_refl_3].*weight,2,'omitnan')./median(weight,'omitnan');
% 
% 
%         % Save Stokes vector in running array
%         image_holder_struc(ind).S_refl = S_refl;
%         image_holder_struc(ind).theta_i = median(theta_i.*weight)./median(weight);
% 
%     else
% 
%         % Save Stokes vector in running array
%         image_holder_struc(ind).S_refl = NaN*ones(4,1);
%         image_holder_struc(ind).theta_i = NaN;
% 
%     end
% 
% end
% 
% Stokes_array = permute(reshape([image_holder_struc.S_refl],4,total_rows,total_cols),[2 3 1]);
% 
% % Extract computed Stokes vector field
% S0 = squeeze(Stokes_array(:,:,1));
% S1 = squeeze(Stokes_array(:,:,2));
% S2 = real(squeeze(Stokes_array(:,:,3)));
% 
% end

end

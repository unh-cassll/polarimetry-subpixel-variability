% Simulates light reflection from a 1D gravity-capillary wave
% as modeled by Melville & Fedorov (2015)
%
% N. Laxague 2024
%
function synthetic_slice_DoLP_example(fignum)

text_x = 0.04;
text_y = 0.91;

in_data = load('MelvilleFedorov2015_length_10cm_steepness_0-375.mat');
x_cm = in_data.x_cm;
eta_cm = in_data.eta_cm;

% in_cmap = load('coolwarm.mat');
% coolwarm = in_cmap.coolwarm;

cmap = cividis(256);
cmap = pink(256);

x = (1e-4:1e-4:20e-2)' - 0.1;
x_cm = x_cm - 0.5;
wse_m = eta_cm/100;
wse_m_filt = smooth(wse_m,101);

[S,I,R] = quickly_get_Stokes_inc_refl(wse_m,x);
[S_filt,I_filt,R_filt] = quickly_get_Stokes_inc_refl(wse_m_filt,x);

% Extract computed Stokes vector field
S0 = squeeze(S(:,1));
S1 = squeeze(S(:,2));
S2 = real(squeeze(S(:,3)));

S0 = interp1(x,S0,x,'pchip');
S1 = interp1(x,S1,x,'pchip');
S2 = interp1(x,S2,x,'pchip');

DoLP = sqrt(S1.^2+S2.^2)./S0;

% Extract computed Stokes vector field
S0_filt = squeeze(S_filt(:,1));
S1_filt = squeeze(S_filt(:,2));
S2_filt = real(squeeze(S_filt(:,3)));

S0_filt = interp1(x,S0_filt,x,'pchip');
S1_filt = interp1(x,S1_filt,x,'pchip');
S2_filt = interp1(x,S2_filt,x,'pchip');

DoLP_filt = sqrt(S1_filt.^2+S2_filt.^2)./S0_filt;

figure(fignum);clf
tlayout = tiledlayout(2,1);
% set(gcf,'Position',figpos)

ax_struc = struct();

nexttile(1)
hold on
% REMOVED: traced rays
% sub_num = 30;
% refl_color = [0.7 0 0];
% inc_color = 0.25*[1 1 1];
% for i = 1:sub_num:length(x)
%     ray_inc_x = [x_cm(i) x_cm(i)+1000*I(i,1)];
%     ray_inc_z = [eta_cm(i) eta_cm(i)+-1000*I(i,3)];
%     plot(ray_inc_x,ray_inc_z,'-','Color',inc_color)
%     ray_refl_x = [x_cm(i) x_cm(i)+1000*R(i,1)];
%     ray_refl_z = [eta_cm(i) eta_cm(i)+1000*R(i,3)];
%     plot(ray_refl_x,ray_refl_z,'-','Color',refl_color)
% end
plot(x_cm,eta_cm,'k-','linewidth',6)
scatter(x_cm,eta_cm,15,DoLP,'filled')
hold off
xlim([0 10])
ylim([-1 2])
clim([0 1])
pbaspect([10 3 1])
colormap(cmap)
xlabel('x [cm]')
ylabel('\eta [cm]')
text(text_x,text_y,'(a)','FontSize',18,'Units','normalized','HorizontalAlignment','center')
grid off
box on
ax_struc(1).ax = gca;
cbar = colorbar;
set(get(cbar,'Title'),'String','DoLP')
cbar.Location = 'northoutside';

nexttile(2)
hold on
% for i = 1:sub_num:length(x)
%     ray_refl_x = [x_cm(i) x_cm(i)+1000*R_filt(i,1)];
%     ray_refl_z = [100*wse_m_filt(i) 100*wse_m_filt(i)+1000*R_filt(i,3)];
%     plot(ray_refl_x,ray_refl_z,'-','Color',refl_color)
%     ray_inc_x = [x_cm(i) x_cm(i)+1000*I_filt(i,1)];
%     ray_inc_z = [100*wse_m_filt(i) 100*wse_m_filt(i)+-1000*I_filt(i,3)];
%     plot(ray_inc_x,ray_inc_z,'-','Color',inc_color)
% end
plot(x_cm,wse_m_filt*100,'k-','linewidth',6)
scatter(x_cm,wse_m_filt*100,15,DoLP_filt,'filled')
hold off
xlim([0 10])
ylim([-1 2])
clim([0 1])
pbaspect([10 3 1])
colormap(cmap)
xlabel('x [cm]')
ylabel('\eta [cm]')
text(text_x,text_y,'(b)','FontSize',18,'Units','normalized','HorizontalAlignment','center')
grid off
box on
ax_struc(2).ax = gca;

tile_cleaner(ax_struc,tlayout)

tlayout.TileSpacing = 'none';

    function [S,I,R] = quickly_get_Stokes_inc_refl(wse_m,x)

        y = 0*x;

        % Virtual camera parameters
        incidence_view_deg = 20;
        azimuth_view_deg = 0;
        height_view_m = 100;
        S_inc = [1 0 0 0];

        % Incident Stokes vector
        S_inc1 = S_inc(1);
        S_inc2 = S_inc(2);
        S_inc3 = S_inc(3);
        S_inc4 = S_inc(4);
        camera_location = [tand(incidence_view_deg)*height_view_m 0 height_view_m]';

        % Solar angles
        sun_polar_ang = 0;
        sun_az_ang = 0;

        % Real index of refraction
        n = 1.34;

        % Compute the surface normal field
        dx = median(diff(x));
        [gx] = gradient(smooth(wse_m,11));
        nx = gx/dx;

        N_mat = NaN*ones(length(nx),3);
        N_mat(:,1) = gx/dx;
        N_mat(:,2) = 0*gx;
        N_mat(:,3) = 0*gx+1;

        N_mat = N_mat./sqrt(sum(N_mat.^2,2));

        nx = N_mat(:,1);
        ny = N_mat(:,2);
        nz = N_mat(:,3);

        S = NaN*ones(length(x),4);
        I = NaN*ones(length(x),3);
        R = I;

        for ind = 1:length(x)

            % Grab surface normal vector and facet position
            N = [nx(ind) ny(ind) nz(ind)]';
            X = x(ind);
            Y = y(ind);
            Z = wse_m(ind);

            % Compute reflected ray
            ray_refl = camera_location - [X Y Z]';
            ray_refl = ray_refl/norm(ray_refl);

            R(ind,:) = ray_refl;

            % Compute incident ray
            ray_inc = ray_refl - (2*ray_refl'*N).*N;

            I(ind,:) = ray_inc;

            % Compute incidence angle with respect to surface facet
            inc_dot_N = sum(ray_inc.*N,1);
            theta_i = abs(acosd(inc_dot_N*-1));
            theta_t = asind(sind(theta_i)./n);

            % Enforce sky-leaving radiance conditions
            polar_ang = acosd([0 0 1]*N);
            az_ang = atan2(N(1,:),N(2,:))*180/pi;
            weight = (cosd(polar_ang-sun_polar_ang).*cosd(az_ang-sun_az_ang)).^2;
            % weight = 0*N(1,:) + 1;

            % Ensure that ray is coming from above horizon
            % TODO: REPLACE THIS WITH MULTIPLE REFLECTION CALCULATION
            pass_cond = ray_inc(3,:) < 0;
            N = N(:,pass_cond);
            ray_inc = ray_inc(:,pass_cond);
            theta_i = theta_i(pass_cond);
            theta_t = theta_t(pass_cond);
            weight = weight(pass_cond);

            if ~isempty(N)

                % Compute h and s vectors (Mobley 2015)
                % and from them, Stokes vector rotation angle alpha
                h_inc = cross(repmat([0;0;1],1,length(theta_i)),ray_inc);
                h_inc = h_inc./sqrt(sum(h_inc.^2));
                s = cross(ray_inc,N);
                s = s./sqrt(sum(s.^2));
                h_cross_s = cross(h_inc,s);

                h_inc_dot_s = sum(h_inc.*s,1);
                alpha_ang = abs(acosd(h_inc_dot_s));
                antiparallel_inds = h_cross_s(3,:).*ray_inc(3,:) < 0;
                alpha_ang(antiparallel_inds) = 360 - alpha_ang(antiparallel_inds);
                c2a = cosd(2*alpha_ang);
                s2a = sind(2*alpha_ang);

                % Rotation matrix, Stokes vector
                % R = [1,0,            0,             0;
                %      0,cosd(2*alpha),-sind(2*alpha),0;
                %      0,sind(2*alpha),cosd(2*alpha), 0;
                %      0,0,            0,             1];

                % Mueller matrix for reflection off surface, KA89
                alpha = 1/2*(tand(theta_i-theta_t)./tand(theta_i+theta_t)).^2;
                eta = 1/2*(sind(theta_i-theta_t)./sind(theta_i+theta_t)).^2;
                gamma_Re = (tand(theta_i-theta_t).*sind(theta_i-theta_t))./(tand(theta_i+theta_t).*sind(theta_i+theta_t));

                % R_AM = [alpha+eta,alpha-eta,0,       0;
                %         alpha-eta,alpha+eta,0,       0;
                %         0,        0,        gamma_Re,0;
                %         0,        0,        0,       gamma_Re];

                % Compute reflected Stokes vector
                %S_refl = R*R_AM*S_inc;
                S_refl_0 = S_inc1*(alpha+eta)+S_inc2*(alpha-eta);
                S_refl_1 = S_inc1*c2a.*(alpha-eta)+S_inc2*c2a.*(alpha+eta)-S_inc3*s2a.*gamma_Re;
                S_refl_2 = S_inc1*s2a.*(alpha-eta)+S_inc2*s2a.*(alpha+eta)+S_inc3*c2a.*gamma_Re;
                S_refl_3 = S_inc4*gamma_Re;

                S(ind,:) = [S_refl_0 S_refl_1 S_refl_2 S_refl_3];

            else

                % Save Stokes vector in running array
                S(ind,:) = NaN*ones(1,4);

            end

        end

    end

end

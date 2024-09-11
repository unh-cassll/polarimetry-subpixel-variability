
% Constants
g = 9.81;       % acceleration due to gravity in m/s^2
sigma = 0.072;  % surface tension in N/m
rho_w = 1020;   % water density in kg/m^3

% Wind forcing
U10_vec = 1:0.5:12;
[~,u_star_a] = logistic_fit_drag(U10_vec,'U10');
rho_a = 1.226;  % density of air, kg/m^3
tau_total = rho_a*u_star_a.^2;

% Prepare wave scale and direction arrays
wave_dir_vec = linspace(-pi/2,pi/2,180)*180/pi;
min_k = 1e-3;
max_k = 1e3;
num_k = 500;
num_winds = length(U10_vec);
num_wave_dir = length(wave_dir_vec);

fetch_m = 100000;

[k,Fk] = Elfouhaily_omni(U10_vec,min_k,max_k,num_k,fetch_m);

B = repmat(k,1,num_winds).^3.*Fk;

B(B<0) = 0;

% Repeat blocks that vary in wavenumber, wind forcing, and wave direction
k_block = permute(repmat(k,[1 num_winds num_wave_dir]),[3 1 2]);
B_block = permute(repmat(B,[1 1 num_wave_dir]),[3 1 2]);
ustar_block = permute(repmat(u_star_a,[num_k 1 num_wave_dir]),[3 1 2]);
U10_block = permute(repmat(U10_vec,[num_k 1 num_wave_dir]),[3 1 2]);
wave_dir_block = permute(repmat(wave_dir_vec,[num_k 1 num_winds]),[2 1 3]);

% Compute wave frequency, celerity, and wave age
k_p = NaN*U10_vec;
for n = 1:num_winds
    ind = find(Fk==max(Fk),1,'first');
    k_p(n) = k(ind);
end
omega_p = sqrt(g*k_p+sigma/rho_w*k_p.^3);
cp = omega_p./k_p;
waveage_block = permute(repmat(cp,[num_k 1 num_wave_dir]),[3 1 2])./U10_block;

% Compute directional spreading function
[deltak] = Elfouhaily_spread(k_block,U10_block,waveage_block,ustar_block);
s_block = atanh(deltak)*2/log(2);
spreading_block = cosd(wave_dir_block).^(2*s_block)/(2);

% Produce 'F_block', the directional wavenumber elevation variance spectrum
F_block = k_block.^-4.*B_block.*spreading_block;

% Compute the upwind mean square slope, mss_u
mss_u_Elf = squeeze(trapz(k,trapz(wave_dir_vec*pi/180,cosd(wave_dir_block).^2.*k_block.^3.*F_block)));

% Compute the wave form stress
tau_w = 0.04*rho_w*u_star_a(:).^2.*mss_u_Elf(:);

% Cox & Munk [1954] mss
mss_u_CM_clean = 0.000+3.16e-3*U10_vec;
mss_c_CM_clean = 0.003+1.92e-3*U10_vec;
mss_u_CM_slick = 0.005+0.78e-3*U10_vec;
mss_c_CM_slick = 0.003+1.56e-3*U10_vec;

figure(1);clf
semilogy(U10_vec,tau_w,'.-')
xlabel('U_{10N} [m s^{-1}]')
ylabel('\tau_w [N m^{-2}]')

figure(2);clf
plot(U10N_m_s,mss_u,'.',U10_vec,mss_u_Elf,'.',U10_vec,trapz(k,k.^-1.*B),'.',U10_vec,mss_u_CM_clean,'--',U10_vec,mss_u_CM_slick,'--','linewidth',3)
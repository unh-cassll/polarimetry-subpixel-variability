% Constants
g = 9.81;       % acceleration due to gravity in m/s^2
sigma = 0.072;  % surface tension in N/m
rho_w = 1020;   % water density in kg/m^3
rho_a = 1.200;  % air density in kg/m^3

fetch_m = 1000; % fetch in meters

U10_vec = 2.5; % wind speed, meters per second

[~,u_star_a] = logistic_fit_drag(U10_vec,'U10');
tau_total = rho_a*u_star_a.^2;

% Prepare wave scale and direction arrays
wave_dir_vec = linspace(-pi/2,pi/2,180)*180/pi;
min_k = 1e-3;
max_k = 2e3;
num_k = 500;
num_winds = length(U10_vec);
num_wave_dir = length(wave_dir_vec);

% Compute omnidirectional wavenumber spectra
[k,Fk] = Elfouhaily_omni(U10_vec,min_k,max_k,num_k,fetch_m);
B = repmat(k,1,num_winds).^3.*Fk;
B(B<0) = 0;

% Repeat blocks that vary in wavenumber, wind forcing, and wave direction
k_block = permute(repmat(k,[1 num_winds num_wave_dir]),[3 1 2]);
B_block = permute(repmat(B,[1 1 num_wave_dir]),[3 1 2]);
ustar_block = permute(repmat(u_star_a,[num_k 1 num_wave_dir]),[3 1 2]);
U10_block = permute(repmat(U10_vec,[num_k 1 num_wave_dir]),[3 1 2]);
wave_dir_block = permute(repmat(wave_dir_vec,[num_k 1 num_winds]),[2 1 3]);

% Compute wave frequency and celerity given wavenumber and current profile
omega_block = sqrt(g*k_block+sigma/rho_w*k_block.^3);%
cp_block = omega_block./k_block;
waveage_block = cp_block./U10_block;

% Compute directional spreading function
[deltak] = Elfouhaily_spread(k_block,U10_block,waveage_block,ustar_block);
s_block = atanh(deltak)*2/log(2);
spreading_block = cosd(wave_dir_block).^(2*s_block)/(2*pi);

% Produce 'F_block', the directional wavenumber elevation variance spectrum
F_block = k_block.^-4.*B_block.*spreading_block;
% Computes slope field components from Ursa Division of Amplitude camera
%
% Inputs:
%   * I0, I45, I90, I135: raw polarized light intensity fields
%   * DoFP simulation flag (true/false): do we simulate division of focal
%   plane array?
%   * blocksize: do we perform discrete spatial averaging and downsampling?
%
% N. Laxague 2024
%
function [Sx,Sy] = compute_slope_field_from_polarization_intensities(I0_stack,I45_stack,I90_stack,I135_stack,dofp_sim_type,blocksize)

% Load in data reduction matrix and DoLP:theta relationship
load('DRMSB.mat')
load('dolp_theta_vecs.mat')

% Only consider 0-90 degree portion of DoLP:theta relationship
ind = find(DOLP_full==max(DOLP_full),1,'first');
DOLP_unique = (0.0001:0.0001:1);
theta_unique = interp1(DOLP_full(1:ind),theta_full(1:ind),DOLP_unique);

% Block spatial averaging
if blocksize(1) > 1
    blockvec = [blocksize blocksize];
    if length(size(I0_stack)) == 3
        blockvec = [blockvec 1];
    end
    [I0_stack] = block_average_and_subsample(I0_stack,blockvec);
    [I45_stack] = block_average_and_subsample(I45_stack,blockvec);
    [I90_stack] = block_average_and_subsample(I90_stack,blockvec);
    [I135_stack] = block_average_and_subsample(I135_stack,blockvec);
end

[s1,~,~] = size(I0_stack);
if mod(s1,2)
    I0_stack = I0_stack(1:s1-1,1:s1-1,:);
    I45_stack = I45_stack(1:s1-1,1:s1-1,:);
    I90_stack = I90_stack(1:s1-1,1:s1-1,:);
    I135_stack = I135_stack(1:s1-1,1:s1-1,:);
end

% DoFP simulation
if ~strcmp(dofp_sim_type,'none')
    [I0_stack,I45_stack,I90_stack,I135_stack] = compute_sparse_intensities(I0_stack,I45_stack,I90_stack,I135_stack,dofp_sim_type);
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

clear I0_stack I45_stack I90_stack I135_stack

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

clear DOLP_int

% Find extraordinarily bright or dark spots
S0_vec = sort(reshape(S0,[],1));
L = length(S0_vec);
if floor(0.00001*L) >= 1
    ind_low = find(S0_vec>S0_vec(floor(0.00001*L)),1,'first');
    S0_low = S0_vec(ind_low);
else
    S0_low = 0;
end
ind_high = find(S0_vec>S0_vec(floor(0.99999*L)),1,'first');
S0_high = S0_vec(ind_high);

% Remove them
inds_remove = S0 > S0_high | S0 < S0_low;
mask = 0*S0+1;
mask(inds_remove) = NaN;

clear S0

% Compute slope field components
Sx = -sind(ORI).*tand(AOI).*mask;
Sy = -cosd(ORI).*tand(AOI).*mask;

Sx = Sx - mean(Sx,'all','omitnan');
Sy = Sy - mean(Sy,'all','omitnan');

Sx(isnan(Sx)) = 0;
Sy(isnan(Sy)) = 0;

addpath ..

fetch_m = 1e6;
U10 = 10;
min_k = 1e-2;
max_k = 1e5;
num_k = 1e4;

% Define constants
g = 9.81; %m/s/s
sigma = 0.072;
rho = 1030;

% Compute drag/friction velocity
%CD = 1.2*10^-3;
%CD = 1e-3*(0.8+0.065*U10);
%u_star = sqrt(CD).*U10;
[~,u_star] = logistic_fit_drag(U10,'U10');

% Compute wavenumber, frequency, and phase speed
k = logspace(log10(min_k),log10(max_k),num_k)';     % rad/m
ang_f = sqrt(g*k+sigma/rho*k.^3);                   % rad/s
c = ang_f./k;                                       % phase speed

s = struct();

km = 370;
cm = sqrt(g/km+0.072/1020*km);

% Capillarity effects
multvec = [0.612 0.75 1 1.38 1.82];
for m = 1:length(multvec)

    % Retrieve wavenumber peak
    k0 = g./(U10.*U10);                                  % rad/m
    X0 = 2.2e4;                                         % "empirical constant"
    omega_c = 0.84*tanh((fetch_m*k0/X0).^0.4).^-0.75;     % critical inverse wave age
    kp = k0.*omega_c.^2;                                  % peak wavenumber
    ang_fp = sqrt(g*kp);                                % peak angular frequency
    cp = ang_fp./kp;                                     % phase speed at peak wavenumber
    omega = U10./cp;                                     % inverse wave age
    gamma = 1.7;                                        % JONSWAP parameter
    alpha_p = 0.006*omega_c.^0.55;
    sigma = 0.08*(1+4*omega_c.^-3);

    % Long wave portion
    Lpm = exp(-5/4*(kp./k).^2);                         % Pierson-Moskowitz shape spectrum
    big_gamma = exp(-(sqrt(k./kp)-1).^2./(2*sigma.^2));
    Jp = gamma.^big_gamma;                              % JONSWAP peak enhancement
    Fp = Lpm.*Jp.*exp(-omega/sqrt(10).*(sqrt(k./kp)-1));  % Long wave side effect function
    Bl = 0.5*alpha_p.*cp./c.*Fp;

    % Short wave portion
    if u_star < cm
        alpha_m = 1 + log(u_star/cm);
    else
        alpha_m = 1 + 3*log(u_star/cm);
    end
    alpha_m = alpha_m*10^-2;
    Fm = exp(-1/4*(k/(km*multvec(m))-1).^2);
    Bh = 0.5*alpha_m.*cm./c.*Fm;

    % Combine portions of spectra
    S = k.^-3.*(Bl+Bh);

    out_K = k;
    out_S = S;

    for n = 1:length(U10)

        ind = find(Bl(:,n) == max(Bl(:,n)),1,'first');

        out_S(:,n) = [k(1:ind).^-3.*Bl(1:ind,n)+S(ind,n)-k(ind).^-3.*Bl(ind,n); S(ind+1:end,n)];

    end

    s(m).F = out_S;
end

Fk = [s.F];

Bk = k.^3.*Fk;
Bk_copy = Bk;
Bk_copy(k<50,:) = NaN;

kp = NaN*ones(length(s),1);

for m = 1:length(s)


    ind = find(Bk(:,m)==max(Bk_copy(:,m),[],'omitnan'),1,'first');
    kp(m) = k(ind);

end

lp = 2*pi./kp;

lp_strings = {};
for m = 1:length(lp)
    lp_strings{m} = ['\lambda_p = ' sprintf('%0.1f',1000*lp(m)) ' mm'];
end

load('coolwarm.mat')

cmap = interp1(1:size(coolwarm,1), coolwarm, linspace(1,size(coolwarm,1),length(s)), 'linear');

figure(111);clf
hold on
loglog(k,k.^3.*Fk,'k-','linewidth',4)

H_color = loglog(k,k.^3.*fliplr(Fk),'linewidth',3);colororder(flipud(cmap));

hold off
legend(H_color,flip(lp_strings),'Location','southwest')
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
xlim([1e0 1e4])
ylim([1e-5 1e-1])
xlabel('k [rad m^{-1}]')
ylabel('B(k) [rad]')

mss = trapz(k,k.^2.*Fk);

figure(112);plot(kp,100*(mss-mss(3))/mss(3),'.')

kp/kp(3)



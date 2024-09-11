%
function subpixel_RMS_slope(fignum)

% Comparison spectrum: 2 km to 0.2 mm wavelengths
min_k = pi*1e-3;
max_k = pi*1e3;
num_k = 1000;
fetch_m = 1e6;

klims = [3e0 3e3];
rmss_lims = [3e-4 3e-1];

dU = 1;
U10 = 3:dU:14;

[k,Fk] = Elfouhaily_omni(U10,min_k,max_k,num_k,fetch_m);

Fk(Fk<0) = 0;

cumu_mss = cumtrapz(k,k.^2.*Fk);
total_mss = trapz(k,k.^2.*Fk);
sub_mss = total_mss - cumu_mss;
sub_mss(sub_mss<0) = 0;
sub_rms_slope = sqrt(sub_mss);

figure(fignum);clf

% axis for wavelength
axis_wavelength=axes('Position',[.1 .1 .8 1e-12]);
set(axis_wavelength,'Units','normalized');
set(axis_wavelength,'Color','none');

% axis for wavenumber with plot
axis_wavenumber=axes('Position',[.1 .24 .8 .7]);
set(axis_wavenumber,'Units','normalized');

hold on
plot(k,sub_rms_slope,'k-','linewidth',3)
plot(k,sub_rms_slope,'linewidth',2);colororder(viridis(length(U10)))
hold off
ylim(rmss_lims)
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
ax.YTick = 10.^(-3:1:0);
ax.YTickLabel = {'0.001','0.01','0.1','1'};
box on

colormap(viridis(length(U10)))
clim([U10(1) U10(end)]+[-1 1]*dU/2)
cbar = colorbar;
cbar.Ticks = U10;
cbar.Location = 'northoutside';
set(get(cbar,'Label'),'String','U_{10} [m s^{-1}]')

% set limits and labels
set(axis_wavenumber,'xlim',klims);
set(axis_wavelength,'xlim',2*pi./flip(klims));
axis_wavelength.XTick = 2*pi./[1e3 1e2 1e1];
axis_wavelength.XScale = 'log';
axis_wavelength.XTickLabel = {'6.3','62.8','628.3'};
axis_wavelength.XDir = 'reverse';
axis_wavelength.FontSize = axis_wavenumber.FontSize;
axis_wavelength.XMinorTick = "off";
xlabel(axis_wavenumber,'k_{cut} [rad m^{-1}]')
xlabel(axis_wavelength,'\lambda_{cut} [mm]')
ylabel(axis_wavenumber,'subpixel RMS slope')


%
function RaDyO2008_DoFP_sim(fignum)

% DoFP synthesis

s = load('RaDyO2008_reprocess.mat');
running_struc = s.running_struc;

s = load('RaDyO_wind_speed_stress.mat');
ustar = s.ustar;
ustar = ustar(1:length(running_struc));
[~,order] = sort(ustar);

cmap = viridis(length(ustar));

ustar_lims = [5 40]/100;

wind_ind = floor(255*((ustar-min(ustar)))/(ustar_lims(2)-min(ustar)))+1;

cmap_full = viridis(256);

dofp_sim_num = 2;
dofp_sim_string = ['dofp_' num2str(dofp_sim_num) 'x' num2str(dofp_sim_num)];

blocknums = [1 2 4 8];
block_ind = 1;

klims = [3e0 3e3];
plims = [-1 1]*100;

labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};
text_x = 0.08;
text_y = 0.93;

figure(fignum);clf
tlayout = tiledlayout(1,3);

for m = 1:3
    nexttile(m)
    text(text_x,text_y,labels{m},'FontSize',22,'HorizontalAlignment','center','Units','normalized')
    hold on
    if ~mod(m,3)
        plot(klims,0*klims,'k--','linewidth',2)
    end
end



k_ind_start = 1;
k_ind_end = 256/blocknums(block_ind)/dofp_sim_num;
k_ind_end = k_ind_end - 2;

lw_thick = 3;
lw_thin = 2;

for ind = 1:length(running_struc)

    n = order(ind);

    if n ~= 10 && n ~=13

        k=running_struc(n).(['block' num2str(blocknums(block_ind+1))]).dofp_none.k;
        Sk=running_struc(n).(['block' num2str(blocknums(block_ind+1))]).dofp_none.Sk;
        mss_off = trapz(k,Sk);      

        k=running_struc(n).(['block' num2str(blocknums(block_ind))]).(dofp_sim_string).k;
        Sk=running_struc(n).(['block' num2str(blocknums(block_ind))]).(dofp_sim_string).Sk;
        mss_on = trapz(k,Sk);     

        k=running_struc(n).(['block' num2str(blocknums(block_ind+1))]).dofp_none.k;
        Sk=running_struc(n).(['block' num2str(blocknums(block_ind+1))]).dofp_none.omniSk;
        Sk = Sk*mss_off/trapz(k,Sk);
        k_tail = k(k_ind_end:end);
        Sk_tail = Sk(k_ind_end:end);
        k = k(k_ind_start:k_ind_end);
        Sk = Sk(k_ind_start:k_ind_end);        

        nexttile(1)
        plot(k_tail,k_tail.*Sk_tail,'-','Color',[0 0 0 0.1],'linewidth',lw_thick)
        plot(k,k.*Sk,'k-','linewidth',lw_thick)
        ax_struc(1).ax = gca;

        k=running_struc(n).(['block' num2str(blocknums(block_ind))]).(dofp_sim_string).k;
        Sk=running_struc(n).(['block' num2str(blocknums(block_ind))]).(dofp_sim_string).omniSk;
        k = k(k_ind_start:k_ind_end);
        Sk = Sk(k_ind_start:k_ind_end);
        Sk = Sk*mss_on/trapz(k,Sk);

        nexttile(2)
        plot(k,k.*Sk,'k-','linewidth',lw_thick)
        ax_struc(2).ax = gca;

        k=running_struc(n).(['block' num2str(blocknums(block_ind))]).(dofp_sim_string).k;
        Sk=running_struc(n).(['block' num2str(blocknums(block_ind))]).(dofp_sim_string).omniSk;
        k = k(k_ind_start:k_ind_end);
        Sk = Sk(k_ind_start:k_ind_end);

        Sk_none = running_struc(n).(['block' num2str(blocknums(block_ind+1))]).dofp_none.omniSk(k_ind_start:k_ind_end);

        nexttile(3)
        plot(k,100*(Sk-Sk_none)./Sk_none,'k-','linewidth',lw_thick)
        ax_struc(3).ax = gca;

    end

end

for ind = 1:length(running_struc)

    n = order(ind);

    if n~=10 && n~=13

        k=running_struc(n).(['block' num2str(blocknums(block_ind+1))]).dofp_none.k;
        Sk=running_struc(n).(['block' num2str(blocknums(block_ind+1))]).dofp_none.Sk;
        mss_off = trapz(k,Sk);      

        k=running_struc(n).(['block' num2str(blocknums(block_ind))]).(dofp_sim_string).k;
        Sk=running_struc(n).(['block' num2str(blocknums(block_ind))]).(dofp_sim_string).Sk;
        mss_on = trapz(k,Sk);     

        k=running_struc(n).(['block' num2str(blocknums(block_ind+1))]).dofp_none.k;
        Sk=running_struc(n).(['block' num2str(blocknums(block_ind+1))]).dofp_none.omniSk;
        Sk = Sk*mss_off/trapz(k,Sk);
        k_tail = k(k_ind_end:end);
        Sk_tail = Sk(k_ind_end:end);
        k = k(k_ind_start:k_ind_end);
        Sk = Sk(k_ind_start:k_ind_end);    

        nexttile(1)
        plot(k_tail,k_tail.*Sk_tail,'linewidth',lw_thin,'Color',[cmap_full(wind_ind(n),:) 0.5])
        plot(k,k.*Sk,'linewidth',lw_thin,'Color',cmap_full(wind_ind(n),:))
        ax_struc(1).ax = gca;


        k=running_struc(n).(['block' num2str(blocknums(block_ind))]).(dofp_sim_string).k;
        Sk=running_struc(n).(['block' num2str(blocknums(block_ind))]).(dofp_sim_string).omniSk;
        k = k(k_ind_start:k_ind_end);
        Sk = Sk(k_ind_start:k_ind_end);
        Sk = Sk*mss_on/trapz(k,Sk);

        nexttile(2)
        plot(k,k.*Sk,'linewidth',lw_thin,'Color',cmap_full(wind_ind(n),:))
        ax_struc(2).ax = gca;

        k=running_struc(n).(['block' num2str(blocknums(block_ind))]).(dofp_sim_string).k;
        Sk=running_struc(n).(['block' num2str(blocknums(block_ind))]).(dofp_sim_string).omniSk;
        k = k(k_ind_start:k_ind_end);
        Sk = Sk(k_ind_start:k_ind_end);

        Sk_none = running_struc(n).(['block' num2str(blocknums(block_ind+1))]).dofp_none.omniSk(k_ind_start:k_ind_end);

        nexttile(3)
        plot(k,100*(Sk-Sk_none)./Sk_none,'linewidth',lw_thin,'Color',cmap_full(wind_ind(n),:))
        ax_struc(3).ax = gca;

    end

end

hold off

B_inds = [1 2];

for ind = 1:length(B_inds)

    n = B_inds(ind);

    nexttile(n)
    box on
    % colororder(cmap)
    ax_struc(n).ax.XScale = 'log';
    ax_struc(n).ax.YScale = 'log';
    ax_struc(n).ax.XTick = 10.^(0:3);
    ax_struc(n).ax.YTick = 10.^(-5:-1);
    xlim(klims)
    ylim([1e-4 1e-1])

    if n == 1 || n == 4
        ylabel('B(k) [rad]')
        cbar = colorbar;
        clim(ustar_lims*100)
        cbar.Ticks = 5:5:40;
        cbar.Location = 'south';
        set(get(cbar,'Label'),'String','u_* [cm s^{-1}]')
        
    else
        ax_struc(n).ax.YTickLabel = '';
    end
        xlabel('k [rad m^{-1}]')

end

ratio_inds = 3;

for ind = 1:length(ratio_inds)

    n = ratio_inds(ind);

    nexttile(n)
    box on
    colororder(cmap)

    ax_struc(n).ax.XScale = 'log';
    ax_struc(n).ax.XTick = 10.^(0:3);
    ax_struc(n).ax.YAxisLocation = 'right';

    xlim(klims)
    ylim(plims)
        xlabel('k [rad m^{-1}]')

    if mod(n,3)
        ax_struc(n).ax.YScale = 'log';
    end
    ylabel('% difference')
    hold on
    plot(max(k)*dofp_sim_num/2*[1 1]*1,plims,'r-','linewidth',2)
    plot(max(k)*dofp_sim_num/2*[1 1]*0.5,plims,'r-.','linewidth',2)
    plot(max(k)*dofp_sim_num/2*[1 1]*0.25,plims,'r:','linewidth',2)
    hold off
    if n == 3
        text(max(k)*dofp_sim_num/2,plims(2)*1.12,'k_{Nyq}','FontSize',18,'Color','r','HorizontalAlignment','center')
        text(max(k)*dofp_sim_num/4,plims(2)*1.12,'$\frac{1}{2}$','Interpreter','LaTeX','FontSize',20,'Color','r','HorizontalAlignment','center')
        text(max(k)*dofp_sim_num/8,plims(2)*1.12,'$\frac{1}{4}$','Interpreter','LaTeX','FontSize',20,'Color','r','HorizontalAlignment','center')
    end

end

tlayout.TileSpacing = 'compact';

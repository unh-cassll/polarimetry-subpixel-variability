
% Draws 'superpixel' microgrid polarizers
%
% N. Laxague 2024
%
function draw_superpixel_elements(fignum)

% Create first component of plot: full focal plane array with disparate
% microgrid polarizers
tilenum = 1;

titlesize = 20;
labelsize = 10;
I_size = 16;

figure(fignum);clf
tlayout = tiledlayout(2,6);

nexttile(tilenum,[2 2])
plot(NaN,NaN);box on;grid off;ax=gca;ax.XTick=[];ax.YTick=[];

X = 1;
Y = 1;

num_pixel = 8;

pos_starts = 0:1:num_pixel-1;

line_spacing = 0.25;

linewidth = 1;

for n = 1:length(pos_starts)

    x0 = pos_starts(n);

    for m = 1:length(pos_starts)

        y0 = pos_starts(m);

        if mod(n,2) && ~mod(m,2) % upper left of superpixel

            plot_vert_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing)

        end

        if ~mod(n,2) && mod(m,2) % lower right of superpixel

            plot_horiz_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing)

        end

        if mod(n,2) && mod(m,2)

            plot_diag_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing) % lower left of superpixel

        end

        if ~mod(n,2) && ~mod(m,2)

            plot_antidiag_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing) % upper right of superpixel

        end

    end

end

figure(fignum)
nexttile(tilenum)
xlim([0 num_pixel]);ylim([0 num_pixel])
pbaspect([1 1 1])
text(num_pixel*X/2,num_pixel*Y*1.05,'DoFP Sensing Array','FontSize',titlesize,'HorizontalAlignment','center')

% Plot pixel-by-pixel blocks
plot_horiz_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)
plot_vert_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)

% Plot superpixel blocks
plot_horiz_lines(fignum,tilenum,linewidth*2,0,0,num_pixel,num_pixel,2*X)
plot_vert_lines(fignum,tilenum,linewidth*2,0,0,num_pixel,num_pixel,2*X)

% Add bounding box to lower-right superpixel
x0 = num_pixel - 2*X;
y0 = 0;

figure(fignum)
hold on
f = fill(x0*[1 1 1 1]+2*X*[0 0 1 1],y0*[1 1 1 1]+2*Y*[0 1 1 0],'w');
f.FaceAlpha = 0.8;
plot([x0 x0+2*X],[y0 y0],'r','linewidth',2*linewidth)
plot([x0 x0+2*X],[y0 y0]+2*Y,'r','linewidth',2*linewidth)
plot([x0 x0],[y0 y0+2*Y],'r','linewidth',2*linewidth)
plot([x0 x0]+2*X,[y0 y0+2*Y],'r','linewidth',2*linewidth)
hold off
text(x0+X/2,y0+Y*3/2,'I_{90}','FontSize',I_size,'HorizontalAlignment','center')
text(x0+X/2,y0+Y/2,'I_{135}','FontSize',I_size,'HorizontalAlignment','center')
text(x0+X*3/2,y0+Y/2,'I_{0}','FontSize',I_size,'HorizontalAlignment','center')
text(x0+X*3/2,y0+Y*3/2,'I_{45}','FontSize',I_size,'HorizontalAlignment','center')
% text(x0+X,y0-Y*0.25,'''Super-pixel''','FontSize',labelsize,'HorizontalAlignment','center')

% Create second component of plot: sparse arrays

num_pixel = num_pixel/2;

tilenum = 3;
nexttile(tilenum)
plot(NaN,NaN);box on;grid off;ax=gca;ax.XTick=[];ax.YTick=[];
xlim([0 num_pixel]);ylim([0 num_pixel])

for n = 1:length(pos_starts)

    x0 = pos_starts(n);

    for m = 1:length(pos_starts)

        y0 = pos_starts(m);

        if mod(n,2) && ~mod(m,2) % upper left of superpixel

            plot_vert_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing)

        end

    end

end

% Plot pixel-by-pixel blocks
plot_horiz_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)
plot_vert_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)

text(num_pixel*X*1.1,num_pixel*Y*1.1,'Sparse Arrays','FontSize',titlesize,'HorizontalAlignment','center')

tilenum = 4;
nexttile(tilenum)
plot(NaN,NaN);box on;grid off;ax=gca;ax.XTick=[];ax.YTick=[];
xlim([0 num_pixel]);ylim([0 num_pixel])

for n = 1:length(pos_starts)

    x0 = pos_starts(n);

    for m = 1:length(pos_starts)

        y0 = pos_starts(m);

        if ~mod(n,2) && ~mod(m,2)

            plot_antidiag_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing) % upper right of superpixel

        end

    end

end

% Plot pixel-by-pixel blocks
plot_horiz_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)
plot_vert_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)

tilenum = 9;
nexttile(tilenum)
plot(NaN,NaN);box on;grid off;ax=gca;ax.XTick=[];ax.YTick=[];
xlim([0 num_pixel]);ylim([0 num_pixel])

for n = 1:length(pos_starts)

    x0 = pos_starts(n);

    for m = 1:length(pos_starts)

        y0 = pos_starts(m);

        if mod(n,2) && mod(m,2)

            plot_diag_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing) % lower left of superpixel

        end


    end

end

% Plot pixel-by-pixel blocks
plot_horiz_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)
plot_vert_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)

tilenum = 10;
nexttile(tilenum)
plot(NaN,NaN);box on;grid off;ax=gca;ax.XTick=[];ax.YTick=[];
xlim([0 num_pixel]);ylim([0 num_pixel])

for n = 1:length(pos_starts)

    x0 = pos_starts(n);

    for m = 1:length(pos_starts)

        y0 = pos_starts(m);

        if ~mod(n,2) && mod(m,2) % lower right of superpixel

            plot_horiz_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing)

        end

    end

end

% Plot pixel-by-pixel blocks
plot_horiz_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)
plot_vert_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)

% Create third component of plot: interpolated arrays

tilenum = 5;
nexttile(tilenum)
plot(NaN,NaN);box on;grid off;ax=gca;ax.XTick=[];ax.YTick=[];
xlim([0 num_pixel]);ylim([0 num_pixel])

for n = 1:length(pos_starts)

    x0 = pos_starts(n);

    for m = 1:length(pos_starts)

        y0 = pos_starts(m);

        plot_vert_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing)

    end

end

% Plot pixel-by-pixel blocks
plot_horiz_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)
plot_vert_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)

text(num_pixel*X*1.1,num_pixel*Y*1.1,'Interpolated Arrays','FontSize',titlesize,'HorizontalAlignment','center')

tilenum = 6;
nexttile(tilenum)
plot(NaN,NaN);box on;grid off;ax=gca;ax.XTick=[];ax.YTick=[];
xlim([0 num_pixel]);ylim([0 num_pixel])

for n = 1:length(pos_starts)

    x0 = pos_starts(n);

    for m = 1:length(pos_starts)

        y0 = pos_starts(m);

        plot_antidiag_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing) % upper right of superpixel

    end

end

% Plot pixel-by-pixel blocks
plot_horiz_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)
plot_vert_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)

tilenum = 11;
nexttile(tilenum)
plot(NaN,NaN);box on;grid off;ax=gca;ax.XTick=[];ax.YTick=[];
xlim([0 num_pixel]);ylim([0 num_pixel])

for n = 1:length(pos_starts)

    x0 = pos_starts(n);

    for m = 1:length(pos_starts)

        y0 = pos_starts(m);

        plot_diag_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing) % lower left of superpixel


    end

end

% Plot pixel-by-pixel blocks
plot_horiz_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)
plot_vert_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)

tilenum = 12;
nexttile(tilenum)
plot(NaN,NaN);box on;grid off;ax=gca;ax.XTick=[];ax.YTick=[];
xlim([0 num_pixel]);ylim([0 num_pixel])

for n = 1:length(pos_starts)

    x0 = pos_starts(n);

    for m = 1:length(pos_starts)

        y0 = pos_starts(m);

        plot_horiz_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing)

    end

end

% Plot pixel-by-pixel blocks
plot_horiz_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)
plot_vert_lines(fignum,tilenum,linewidth,0,0,num_pixel,num_pixel,X)

tlayout.TileSpacing = 'compact';

%% This is where we keep the individual plotting functions

    function plot_vert_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing)

        num_lines = floor(X/line_spacing)+1;

        x = x0:line_spacing:x0+X;

        figure(fignum)
        if ~isempty(tilenum)
            nexttile(tilenum)
        end
        hold on
        for i = 1:num_lines

            plot(x(i)*[1 1],y0+[0 Y],'k','linewidth',linewidth)

        end
        hold off

    end

    function plot_horiz_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing)

        num_lines = floor(Y/line_spacing)+1;

        y = y0:line_spacing:y0+Y;

        figure(fignum)
        if ~isempty(tilenum)
            nexttile(tilenum)
        end
        hold on
        for i = 1:num_lines

            plot(x0+[0 X],y(i)*[1 1],'k','linewidth',linewidth)

        end
        hold off

    end

    function plot_diag_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing)

        L = sqrt(X^2+Y^2);

        num_lines = floor(L/line_spacing);
        L_vals = 0:L/num_lines:L-L/num_lines;

        draw_ang = atand(Y/X);

        xvals = L_vals*cosd(draw_ang);

        figure(fignum)
        if ~isempty(tilenum)
            nexttile(tilenum)
        end
        hold on
        for i = 1:num_lines

            xstart = xvals(i);

            xplot_1 = 0;
            ystart = xstart*tand(draw_ang);
            yplot_1 = ystart + xstart*tand(draw_ang);

            xplot_2 = xstart + ystart*tand(draw_ang);
            yplot_2 = 0;

            if yplot_1 > Y

                yplot_1 = Y;
                xplot_1 = xstart - (Y-ystart)*tand(draw_ang);

                xplot_2 = X;
                yplot_2 = xplot_1;

            end

            plot([xplot_1 xplot_2]+x0,[yplot_1 yplot_2]+y0,'k','linewidth',linewidth)

        end
        hold off

    end

    function plot_antidiag_lines(fignum,tilenum,linewidth,x0,y0,X,Y,line_spacing)

        L = sqrt(X^2+Y^2);

        num_lines = floor(L/line_spacing);
        L_vals = 0:L/num_lines:L-L/num_lines;

        draw_ang = -atand(Y/X);

        xvals = L_vals*cosd(draw_ang);

        figure(fignum)
        if ~isempty(tilenum)
            nexttile(tilenum)
        end
        hold on
        for i = 1:num_lines

            xstart = xvals(i);

            xplot_1 = 0+X;
            ystart = xstart*tand(draw_ang)+Y;
            yplot_1 = ystart + xstart*tand(draw_ang);

            xplot_2 = xstart + ystart*tand(draw_ang)+X;
            yplot_2 = 0;

            if yplot_1 > Y

                yplot_1 = Y;
                xplot_1 = xstart - (Y-ystart)*tand(draw_ang)-X;

                xplot_2 = X;
                yplot_2 = xplot_1;

            end

            if yplot_1 < 0

                xplot_1 = xplot_1 - X;
                xplot_2 = xplot_2 - X;

                yplot_1 = yplot_1 + Y;
                yplot_2 = yplot_2 + Y;

            end

            plot([xplot_1 xplot_2]+x0,[yplot_1 yplot_2]+y0,'k','linewidth',linewidth)

        end
        hold off

    end

end
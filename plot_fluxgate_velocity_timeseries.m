%export flux gate velocity time series for figures
clearvars; close all;

%create a map
cd /Users/ellynenderlin/Research/miscellaneous/Greenland-ITSLIVE_1985-2018/
[A,R] = geotiffread('Greenland_ITSLIVEavg_1985-2018.tif');
vx = R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX;
vy = R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY;
v = A; clear A R;
map_fig = figure; set(map_fig,'position',[50 50 600 1200]); 
imagesc(vx,vy,v); axis xy equal; hold on;
vel_cmap = colormap(parula(10001)); vel_cmap(1,:) = [1 1 1]; colormap(gca,vel_cmap);
cbar = colorbar; set(gca,'clim',[0 500]); cbar.Label.String = 'speed (m/yr)';
set(gca,'ylim',[-3.35e6 -0.65e6],'xlim',[-6.5e5 8.5e5],...
    'ytick',[-3.25e6:0.5e6:-0.75e6],'yticklabel',[-3250:500:-750],...
    'xtick',[-6e5:2e5:8e5],'xticklabel',[-600:200:800],...
    'fontsize',18);
xlabel('Easting (km) ','fontsize',20); ylabel('Northing (km) ','fontsize',20); 

%load the Greenland_GIC_centerlines.mat that Kate has updated using
%interpolate_to_fluxgates6.m so that it has velocities at the flux gates
cd /users/ellynenderlin/Research/NASA_Greenland-Periph-Mapping/centerlines/ %need to change this path to where your file is located
load Greenland_GIC_centerlines.mat;
cd ../fluxgates %if your directory is set the same, this will jump up one directory from centerlines, then go into fluxgates and save the figures there

%default is a loop through all glaciers with flux gates but you can
%comment-out "for i = 1:length(term)" and the companion "end"  and simply
%specify "i=whatever" if you know the number glacier you want to create a figure for
for i = 1:length(term)
    
    if ~isempty(term(i).gateX) %check that there is a fluxgate line (we can also use this to check if you need to draw a fluxgate despite having a centerline for a glacier)
        
        speed_fig = figure; set(gcf,'position',[50 50 800 400]); %make a figure and set its position in pixel units [lower left X, lower left Y, width, height]
        plot(term(i).gateVdateavg,term(i).fluxVelavg,'xk','linewidth',2); %plot the points as Xs that are black (add a - before the x if you want to plot a line between points)
        grid on; set(gca,'xlim',[1985 2018],'xtick',[1985:5:2020],'fontsize',20); %if you want to set y-limits, you'd specific ylim like I did xlim here
        xlabel('Year','fontsize',20); ylabel('Speed (m/yr)','fontsize',20); %20 might not be big enough, maybe 24 or 28 is better
        leg = legend(['RGI ID: ',num2str(term(i).RGIref)]); %this adds a legend with the RGI reference, you could make it a title instead "title(['RGI ID: ',num2str(term(i).RGIref)],'fontsize',20);
        drawnow;
        
        %save the figure
        if i<10
            saveas(speed_fig,['00',num2str(i),'_fluxgate-speed-timeseries.png'],'png');
        elseif i >=10 && i<100
            saveas(speed_fig,['0',num2str(i),'_fluxgate-speed-timeseries.png'],'png');
        else
            saveas(speed_fig,[num2str(i),'_fluxgate-speed-timeseries.png'],'png');
        end
        
        close gcf; drawnow;
        
        %plot point on map
        figure(map_fig); plot(nanmean(term(i).X),nanmean(term(i).Y),'ok','markersize',10,'markerfacecolor','m'); hold on;
        
        %plot an X on top of circle if there are no associated velocity data
        if sum(sum(~isnan(term(i).fluxV))) == 0
            figure(map_fig); plot(nanmean(term(i).X),nanmean(term(i).Y),'xk','markersize',10,'linewidth',2); hold on;
        end
        drawnow;
    else
        %plot point on map
        figure(map_fig); plot(nanmean(term(i).X),nanmean(term(i).Y),'ok','markersize',10,'markerfacecolor','w'); hold on;
    end
    
end
%save the map figure
figure(map_fig);
saveas(map_fig,'GreenlandGIC_discharge-coverage-map.png','png');




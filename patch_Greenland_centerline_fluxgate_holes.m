%loop through all the elevation data for each study site to locate the grounding line
clearvars; close all; warning off; beep on

%make sure matlab knows where to look for supporting files
addpath('/users/ellynenderlin/mfiles','/users/ellynenderlin/mfiles/general');

%navigate to the appropriate parent directory (move files to your local computer to MASSIVELY speed up)
% root_dir = '/volumes/APicebergs/Greenland-GICs/';
% ArcticDEMdir = '/volumes/APicebergs/Greenland-GICs/';
root_dir = '/Users/ellynenderlin/Research/NASA_Greenland-Periph-Mapping/';
ArcticDEMdir = '/Users/ellynenderlin/Research/miscellaneous/ArcticDEM10m/';

%load the RGI shapefiles
cd_to_outlines = ['cd ',root_dir,'RGI-outlines/']; eval(cd_to_outlines);
S = shaperead('05_rgi50_GreenlandPeriphery_PS_BoxIDs_TW.shp');
for i = 1:length(S)
    BoxID(i) = S(i).BoxID;
end

%load the georeferencing info for the ArcticDEM (subset for each glacier in loop because it is HUGE)
cd_to_arcticDEM = ['cd ',ArcticDEMdir]; eval(cd_to_arcticDEM);
info = geotiffinfo('ArcticDEM_10m_Greenland_Mosaic.tif'); %pull spatial referencing info for the FULL dem
ArcticDEMx = single([info.SpatialRef.XWorldLimits(1)+0.5*info.SpatialRef.CellExtentInWorldX:info.SpatialRef.CellExtentInWorldX:info.SpatialRef.XWorldLimits(2)-0.5*info.SpatialRef.CellExtentInWorldX]);
ArcticDEMy = single([info.SpatialRef.YWorldLimits(2)-0.5*info.SpatialRef.CellExtentInWorldY:-info.SpatialRef.CellExtentInWorldY:info.SpatialRef.YWorldLimits(1)+0.5*info.SpatialRef.CellExtentInWorldY]);

%load the GIMP DEM to help with visualization
% cd_to_DEM = ['cd ',root_dir,'GIMP/']; eval(cd_to_DEM);
cd ../GIMP/
[I,R] = geotiffread('gimpdem_90m_v01.1.tif');
dem.x = single([R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX]);
dem.y = single([R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY]);
dem.z = single(I);
clear I R;
[demx,demy] = meshgrid(dem.x,dem.y);
%plot the DEM and the image mosaic
DEM_fig = figure; set(DEM_fig,'position',[25 50 800 600]);
sub1 = subplot(1,2,1);
imagesc(dem.x,dem.y,dem.z); axis xy equal; hold on;
dem_cmap = colormap(parula(10001)); dem_cmap(1,:) = [0 0 0];
colormap(gca,dem_cmap); cbar = colorbar;
set(gca,'xlim',[min(dem.x) max(dem.x)],'ylim',[min(dem.y) max(dem.y)]);
sub2 = subplot(1,2,2);
drawnow;

%load the velocity mosaic
cd_to_vels = ['cd ',root_dir,'Greenland-VelMosaic_1995-2015/']; eval(cd_to_vels);
[I,R] = geotiffread('greenland_vel_mosaic250_vx_v1');
V.x = single([R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX]);
V.y = single([R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY]);
V.vx = single(I); V.vx(V.vx==-2.0000e+09) = NaN; %velocity in x-direction in m/yr
clear I R;
[I,~] = geotiffread('greenland_vel_mosaic250_vy_v1');
V.vy = single(I); V.vy(V.vy==-2.0000e+09) = NaN; %velocity in y-direction in m/yr
clear I;
[VXgrid,VYgrid] = meshgrid(V.x,V.y);

%load the existing centerline & fluxgate data
cd_to_centerlines = ['cd ',root_dir,'centerlines/']; eval(cd_to_centerlines);
load Greenland_GIC_centerlines.mat;

%loop through the Landsat scenes & extract coordinates 
% cd_to_landsats = ['cd ',root_dir,'GIMP/landsat8_scenes/']; eval(cd_to_landsats);
cd_to_landsats = ['cd ',ArcticDEMdir]; eval(cd_to_landsats); cd ../GIMP/landsat8_scenes/
Landsats = dir('L*_PS.mat'); disp('Extracting corner coordinates from Landsat files... this will take several minutes');
for i = 1:length(Landsats)
    disp(['Landsat scene #',num2str(i),' of ',num2str(length(Landsats)),': ',Landsats(i).name]);
    load_landsat = ['load ',Landsats(i).name]; eval(load_landsat); %load the scene
    imx_lims(i,:) = [min(im.x) max(im.x)]; imy_lims(i,:) = [min(im.y) max(im.y)]; 
    clear im;
end

%loop through the term structure & identify what data are missing
for j = [266, 321, 423] %1:length(term)
    %trace the centerline & extract profiles if necessary
    if isempty(term(j).centerline)
        disp(['No centerline: ',num2str(j)]);
        %identify the landsat scene that has the terminus box closest to its center & load it
        term_distance = sqrt((nanmean(term(j).X(~isnan(term(j).X)))-nanmean(imx_lims,2)).^2 + (nanmean(term(j).Y(~isnan(term(j).Y)))-nanmean(imy_lims,2)).^2); 
        [mindist,minref] = min(term_distance);
        % cd_to_landsats = ['cd ',root_dir,'GIMP/landsat8_scenes/']; eval(cd_to_landsats);
        cd_to_landsats = ['cd ',ArcticDEMdir]; eval(cd_to_landsats); cd ../GIMP/landsat8_scenes/
        load_landsat = ['load ',Landsats(minref).name]; eval(load_landsat); %load the scene
        [imx,imy] = meshgrid(im.x,im.y);
        clear term_distance mindist minref;
        
        %plot the scene in the 2nd panel in the DEM figure
        figure(DEM_fig); subplot(sub2);
        imagesc(im.x,im.y,im.z_adjust); axis xy equal; hold on;
        plot(term(j).X(~isnan(term(j).X)),term(j).Y(~isnan(term(j).Y)),'-y','linewidth',2); hold on;
        colormap(sub2,'gray');
        set(sub2,'xlim',[min(im.x) max(im.x)],'ylim',[min(im.y) max(im.y)]);
        subplot(sub1);
        plot(term(j).X(~isnan(term(j).X)),term(j).Y(~isnan(term(j).Y)),'-m','linewidth',2); hold on;
        set(sub1,'clim',[0 max(max(dem.z(find(dem.y>=min(im.y) & dem.y<=max(im.y)),find(dem.x>=min(im.x) & dem.x<=max(im.x)))))]); %stretch colors in DEM
        set(sub1,'xlim',[min(im.x) max(im.x)],'ylim',[min(im.y) max(im.y)]);
        drawnow;
        
        %automatically zoom in to the glacier extent
        figure(DEM_fig); subplot(sub1);
        plot(S(term(j).RGIshaperef).X, S(term(j).RGIshaperef).Y, '-k','linewidth',2); hold on;
        set(gca,'xlim',[min([S(term(j).RGIshaperef).X term(j).X]) max([S(term(j).RGIshaperef).X term(j).X])],'ylim',[min([S(term(j).RGIshaperef).Y term(j).Y]) max([S(term(j).RGIshaperef).Y term(j).Y])]);
        set(gca,'clim',[0 max(max(dem.z(find(dem.y>=min(S(term(j).RGIshaperef).Y) & dem.y<=max(S(term(j).RGIshaperef).Y)),find(dem.x>=min(S(term(j).RGIshaperef).X) & dem.x<=max(S(term(j).RGIshaperef).X)))))]); %stretch colors in DEM
        figure(DEM_fig); subplot(sub2);
        plot(S(term(j).RGIshaperef).X,S(term(j).RGIshaperef).Y,'-r','linewidth',2); hold on;
        set(gca,'xlim',[min([S(term(j).RGIshaperef).X term(j).X]) max([S(term(j).RGIshaperef).X term(j).X])],'ylim',[min([S(term(j).RGIshaperef).Y term(j).Y]) max([S(term(j).RGIshaperef).Y term(j).Y])]);
        drawnow;
        
        %manually zoom in more
        beep;
        disp('in the image, click on the UL & LR corners of a region where you''ll draw the centerline');
        figure(DEM_fig); subplot(sub2); [a] = ginput(2); zoom_xlims = a(:,1); zoom_ylims = a(:,2); clear a;
        set(sub1,'xlim',sort(zoom_xlims),'ylim',sort(zoom_ylims));
        set(sub2,'xlim',sort(zoom_xlims),'ylim',sort(zoom_ylims));
        
        %trace the centerline
        disp('in the image, trace the centerline starting from the seaward side of the terminus box');
        figure(DEM_fig); subplot(sub2); [~,~,~,xi,yi] = improfile; %be patient, a + cursor will appear, double right click to terminate
        
        %create a regularly interpolated version of the centerline
        spacer = abs(info.SpatialRef.CellExtentInWorldX); %spacing increment (same as ArcticDEM)
        profx_points = []; profy_points = [];
        if max(xi)-min(xi) > max(yi) - min(yi)
            profy = fit(xi,yi,'smoothingspline'); %longer in x-direction so spline works better w/ x as the independent variable
            for k = 2:length(xi)
                center_orient = atand((yi(k)-yi(k-1))./(xi(k)-xi(k-1)));
                profx_points = [profx_points; [xi(k-1):sign(xi(k)-xi(k-1))*abs(spacer*sind(center_orient)):xi(k)]'];
                profy_points = [profy_points; profy(xi(k-1):sign(xi(k)-xi(k-1))*abs(spacer*sind(center_orient)):xi(k))];
                clear center_orient;
            end
            clear profy;
        else
            profx = fit(yi,xi,'smoothingspline'); %longer in y-direction so spline works better w/ y as the independent variable
            for k = 2:length(xi)
                center_orient = atand((yi(k)-yi(k-1))./(xi(k)-xi(k-1)));
                profy_points = [profy_points; [yi(k-1):sign(yi(k)-yi(k-1))*abs(spacer*cosd(center_orient)):yi(k)]'];
                profx_points = [profx_points; profx(yi(k-1):sign(yi(k)-yi(k-1))*abs(spacer*cosd(center_orient)):yi(k))];
                clear center_orient;
            end
            clear profx;
        end
        clear xi yi;
        
        %add the profile to a structure & export a shapefile
        term(j).RGIX = S(term(j).RGIshaperef).X; term(j).RGIY = S(term(j).RGIshaperef).Y;
        term(j).centerX = profx_points'; term(j).centerY = profy_points';
        %convert the centerline coordinates to along-profile distance from the origin
        term(j).centerline(1) = 0;
        for k = 2:length(term(j).centerX)
            term(j).centerline(k) = term(j).centerline(k-1)+sqrt((term(j).centerX(k)-term(j).centerX(k-1)).^2 + (term(j).centerY(k)-term(j).centerY(k-1)).^2);
        end
        cd_to_centerlines = ['cd ',root_dir,'centerlines/']; eval(cd_to_centerlines);
        s.Geometry = 'Polyline';
        s.BoundingBox = double([min(term(j).centerX) min(term(j).centerY); max(term(j).centerX) max(term(j).centerY)]);
        s.X = double(term(j).centerX); s.Y = double(term(j).centerY);
        s.RGIId = S(term(j).RGIshaperef).RGIId;
        s.BoxID = term(j).BoxID;
        if j<10
            shapewrite(s,['centerline_00',num2str(j),'.shp']);
        elseif j >=10 && j<100
            shapewrite(s,['centerline_0',num2str(j),'.shp']);
        else
            shapewrite(s,['centerline_',num2str(j),'.shp']);
        end
        save('Greenland_GIC_centerlines.mat','term','-v7.3');
        clear s profx* profy*;
        
        %extract the speed time series (m/yr) along the centerline
        cd_to_ITSLIVE = ['cd ',root_dir,'ITS_LIVE/']; eval(cd_to_ITSLIVE);
%         cd_to_landsats = ['cd ',ArcticDEMdir]; eval(cd_to_landsats); cd ../Greenland-ITSLIVE_1985-2018/
        speed_fig = figure; set(speed_fig,'position',[200 50 500 800]);
        ITSLIVEs = dir('GRE_*.nc'); speed_cmap = colormap(parula(length(ITSLIVEs)));
        for k=1:length(ITSLIVEs)
            %load the velocity coordinates
            v_x = ncread(ITSLIVEs(k).name,'x');
            v_y = ncread(ITSLIVEs(k).name,'y');
            
            %determine the spatial subset over which to extract data
            vxrefs = find(v_x >= min(term(j).centerX) & v_x <= max(term(j).centerX));
            vyrefs = find(v_y >= min(term(j).centerY) & v_y <= max(term(j).centerY));
            startloc = [min(vxrefs)-1 min(vyrefs)-1]; count = [range(vxrefs)+3 range(vyrefs)+3];
            
            %read in the subset
            vx = rot90(ncread(ITSLIVEs(k).name,'vx',startloc,count)); vx(vx==-32767) = NaN;
            vy = rot90(ncread(ITSLIVEs(k).name,'vy',startloc,count)); vy(vy==-32767) = NaN;
            verr = rot90(ncread(ITSLIVEs(k).name,'v_err',startloc,count));
            vdate = rot90(ncread(ITSLIVEs(k).name,'date',startloc,count));
            vdt = rot90(ncread(ITSLIVEs(k).name,'dt',startloc,count));
            [vyr,vmo,vday,vhh,vmm,vss] =  datevec(double(vdate)); vyr(vyr==0) = NaN; vmo(vmo==0) = NaN; vmo(isnan(vmo)) = 1;
            if mod(nanmean(vyr(~isnan(vyr))),4) ~=0; modays = [31 28 31 30 31 30 31 31 30 31 30 31]; else; modays = [31 29 31 30 31 30 31 31 30 31 30 31]; end
            cumdays = cumsum(modays); cumdays = [0 cumdays(1:11)];
            vdecidate = vyr+(cumdays(vmo)+vday)/sum(modays);
            vdecidt = vdt/sum(modays);
            
            %flip the y axis to agree with rotated grids
            v_y = flipud(v_y); vyrefs = find(v_y >= min(term(j).centerY) & v_y <= max(term(j).centerY));
            
            %interpolate the the centerline profile & add to structure
            [vxgrid,vygrid] = meshgrid(v_x(min(vxrefs)-1:max(vxrefs)+1),v_y(min(vyrefs)-1:max(vyrefs)+1));
            term(j).centerVdate(k,:) = interp2(vxgrid,vygrid,vdecidate,term(j).centerX,term(j).centerY,'nearest');
            term(j).centerVdt(k,:) = interp2(vxgrid,vygrid,vdecidt,term(j).centerX,term(j).centerY,'nearest');
            term(j).centerV(k,:) = interp2(vxgrid,vygrid,sqrt(vx.^2 + vy.^2),term(j).centerX,term(j).centerY);
            term(j).centerVdir(k,:) = interp2(vxgrid,vygrid,atand(vy./vx),term(j).centerX,term(j).centerY);
            term(j).centerVerr(k,:) = interp2(vxgrid,vygrid,verr,term(j).centerX,term(j).centerY);
            
            %plot
            figure(speed_fig);
            plot(term(j).centerline,term(j).centerV(k,:),'-','color',speed_cmap(k,:),'linewidth',1); hold on;
            drawnow;
            clear v_x v_y vx* vy* v*date v*dt;
        end
        figure(speed_fig);
        legend(num2str([1985:1:2018]'),'location','northeast','fontsize',12);
        set(gca,'fontsize',14); grid on;
        xlabel('Distance from terminus box edge (m)','fontsize',14); ylabel('Speed (m/yr)','fontsize',14);
        %save the speed profiles
        eval(cd_to_centerlines);
        if j<10
            saveas(speed_fig,['00',num2str(j),'_speed-profiles.png'],'png');
        elseif j >=10 && j<100
            saveas(speed_fig,['0',num2str(j),'_speed-profiles.png'],'png');
        else
            saveas(speed_fig,[num2str(j),'_speed-profiles.png'],'png');
        end
        
        
        %pull elevations for the profile from the 10 m-resolution ArcticDEM for the glacier
        cd_to_arcticDEM = ['cd ',ArcticDEMdir]; eval(cd_to_arcticDEM);
        dem_xsub = [find(ArcticDEMx<=min([S(term(j).RGIshaperef).X term(j).X]),1,'last'),find(ArcticDEMx>=max([S(term(j).RGIshaperef).X term(j).X]),1,'first')];
        dem_ysub = [find(ArcticDEMy>=max([S(term(j).RGIshaperef).Y term(j).Y]),1,'last'),find(ArcticDEMy<=min([S(term(j).RGIshaperef).Y term(j).Y]),1,'first')];
        I = imread('ArcticDEM_10m_Greenland_Mosaic.tif','PixelRegion',{dem_ysub,dem_xsub}); I(I==-9999) = NaN;
        ArcticDEMz = I; clear I;
        [ArcticDEM_xgrid,ArcticDEM_ygrid] = meshgrid(ArcticDEMx(min(dem_xsub):max(dem_xsub)),ArcticDEMy(min(dem_ysub):max(dem_ysub)));
        term(j).centerZdate(1,:) = 2008; term(j).centerZdate(2,:) = 2018;
        term(j).centerZ(1,:) = interp2(double(dem.x),double(dem.y),dem.z,double(term(j).centerX),double(term(j).centerY)); %GIMP elevation profile (from ~2008)
        term(j).centerZ(2,:) = interp2(double(ArcticDEMx(min(dem_xsub):max(dem_xsub))),double(ArcticDEMy(min(dem_ysub):max(dem_ysub))),ArcticDEMz,double(term(j).centerX),double(term(j).centerY)); %ArcticDEM profile (from ~2018)
        [lon,lat] = ps2wgs(term(j).centerX,term(j).centerY); geoid = geoidheight(lat,lon);
        %plot the elevation profiles & identify terminus positions
        prof_fig = figure; set(prof_fig,'position',[50 50 500 500]);
        figure(prof_fig); 
        if size(term(j).centerZ(1,:),2) ~= size(geoid,2)
            plot(term(j).centerline,term(j).centerZ(1,:)-geoid','-b','linewidth',1); hold on;
        else
            plot(term(j).centerline,term(j).centerZ(1,:)-geoid,'-b','linewidth',1); hold on;
        end
        legend('GIMP','location','northwest');
        disp('click on the base of the GIMP DEM terminus');
        beep;
        clear b prof_ref;
        figure(prof_fig); b = ginput(1); %get the x & y coordinates of the terminus point in the plot
        [~,prof_ref] = min(abs(term(j).centerline-b(1))); %find the reference for the closest neighboring centerline point
        term(j).centerZterm(1,:) = [term(j).centerX(prof_ref(1)) term(j).centerY(prof_ref(1))]; %x,y of centerline terminus position
        clear b prof_ref;
        figure(prof_fig); 
        if size(term(j).centerZ(2,:),2) ~= size(geoid,2)
            plot(term(j).centerline,term(j).centerZ(2,:)-geoid','-k','linewidth',2); hold on;
        else
            plot(term(j).centerline,term(j).centerZ(2,:)-geoid,'-k','linewidth',2); hold on;
        end
        legend('GIMP','ArcticDEM','location','northwest');
        disp('click on the base of the ArcticDEM terminus');
        figure(prof_fig); b = ginput(1); %get the x & y coordinates of the terminus point in the plot
        [~,prof_ref] = min(abs(term(j).centerline-b(1))); %find the reference for the closest neighboring centerline point
        term(j).centerZterm(2,:) = [term(j).centerX(prof_ref(1)) term(j).centerY(prof_ref(1))]; %x,y of centerline terminus position
        clear b prof_ref;
        figure(prof_fig); prof_zlim = get(gca,'ylim'); set(gca,'ylim',[0 max(prof_zlim)]); clear prof_zlim;
        clear lon lat;
        
        %plot the box-centerline intersection point on the DEM & image
        [box_interceptX,box_interceptY] = polyxpoly(term(j).centerX,term(j).centerY,term(j).X,term(j).Y); %find where the centerline intersects the terminus box
        intercept_dist = sqrt((term(j).centerX(1)-box_interceptX).^2 + (term(j).centerY(1)-box_interceptY).^2); %approximate distance from the first centerline point
        [~,maxref] = max(intercept_dist); clear intercept_dist; %use the most inland intersection in case the seaward edge of the terminus box also intersects the centerline
        [~,minref] = min(sqrt((term(j).centerX-box_interceptX(maxref)).^2 + (term(j).centerY-box_interceptY(maxref)).^2));
        figure(prof_fig); plot(term(j).centerline(minref),term(j).centerZ(2,minref)-geoid(minref),'xr','linewidth',2); hold on;
        term(j).centerTerm = term(j).centerline(minref); term(j).centerTermref = minref;
        eval(cd_to_centerlines);
        save('Greenland_GIC_centerlines.mat','term','-v7.3');
        
        %decide if there is a clear slope break indicating floating ice at the terminus
        prompt = 'Is there a flattening of slope below ~50m elevation suggesting floating ice INLAND of the box intercept (y/n)?';
        str = input(prompt,'s');
        if strmatch(str,'y')==1
            disp('in the elevation profile, click on the slope break');
            figure(prof_fig); b = ginput(1); %click on the point in the figure & output the centerline position to b
            [~,prof_ref] = min(abs(term(j).centerline-b(1))); %find the reference for the closest neighboring centerline point
            term(j).centerGL = [term(j).centerX(prof_ref) term(j).centerY(prof_ref) term(j).centerZ(2,prof_ref)-geoid(prof_ref)]; %x,y,z of centerline grounding line
            plot(term(j).centerline(prof_ref),term(j).centerZ(2,prof_ref)-geoid(prof_ref),'+r','linewidth',2); hold on;
            legend('GIMP','ArcticDEM','box intersection','grounding line','location','northwest');
            clear b prof_ref;
        else
            disp('automatically assigning the centerline flux gate to the inland edge of the terminus box');
            term(j).centerGL = [box_interceptX(maxref) box_interceptY(maxref) term(j).centerZ(2,minref)-geoid(minref)]; %x,y,z of centerline box intersection
            figure(prof_fig); legend('GIMP','ArcticDEM','box intersection','location','northwest')
        end
        figure(prof_fig); set(gca,'fontsize',14); grid on;
        xlabel('Distance from terminus box edge (m)','fontsize',14); ylabel('Elevation (m a.s.l.)','fontsize',14);
        %save the elevation profiles
        if j<10
            saveas(prof_fig,['00',num2str(j),'_elev-profiles.png'],'png');
        elseif j >=10 && j<100
            saveas(prof_fig,['0',num2str(j),'_elev-profiles.png'],'png');
        else
            saveas(prof_fig,[num2str(j),'_elev-profiles.png'],'png');
        end
        clear prof;
        close(prof_fig); 
    end
    
    
    %trace the fluxgate if necessary
    if isempty(term(j).gateZ)
        disp(['No flux gate: ',num2str(j)]);
        
        %identify the landsat scene that has the terminus box closest to its center & load it
        term_distance = sqrt((nanmean(term(j).X(~isnan(term(j).X)))-nanmean(imx_lims,2)).^2 + (nanmean(term(j).Y(~isnan(term(j).Y)))-nanmean(imy_lims,2)).^2); 
        [mindist,minref] = min(term_distance);
        % cd_to_landsats = ['cd ',root_dir,'GIMP/landsat8_scenes/']; eval(cd_to_landsats);
        cd_to_landsats = ['cd ',ArcticDEMdir]; eval(cd_to_landsats); cd ../GIMP/landsat8_scenes/
        load_landsat = ['load ',Landsats(minref).name]; eval(load_landsat); %load the scene
        [imx,imy] = meshgrid(im.x,im.y);
        clear term_distance mindist minref;
        
        %plot the scene in the 2nd panel in the DEM figure
        figure(DEM_fig); subplot(sub2);
        imagesc(im.x,im.y,im.z_adjust); axis xy equal; hold on;
        plot(term(j).X(~isnan(term(j).X)),term(j).Y(~isnan(term(j).Y)),'-y','linewidth',2); hold on;
        colormap(sub2,'gray');
        set(sub2,'xlim',[min(im.x) max(im.x)],'ylim',[min(im.y) max(im.y)]);
        subplot(sub1);
        plot(term(j).X(~isnan(term(j).X)),term(j).Y(~isnan(term(j).Y)),'-m','linewidth',2); hold on;
        set(sub1,'clim',[0 max(max(dem.z(find(dem.y>=min(im.y) & dem.y<=max(im.y)),find(dem.x>=min(im.x) & dem.x<=max(im.x)))))]); %stretch colors in DEM
        set(sub1,'xlim',[min(im.x) max(im.x)],'ylim',[min(im.y) max(im.y)]);
        drawnow;
        
        %automatically zoom in to the glacier extent
        figure(DEM_fig); subplot(sub1);
        plot(S(term(j).RGIshaperef).X, S(term(j).RGIshaperef).Y, '-k','linewidth',2); hold on;
        set(gca,'xlim',[min([S(term(j).RGIshaperef).X term(j).X]) max([S(term(j).RGIshaperef).X term(j).X])],'ylim',[min([S(term(j).RGIshaperef).Y term(j).Y]) max([S(term(j).RGIshaperef).Y term(j).Y])]);
        set(gca,'clim',[0 max(max(dem.z(find(dem.y>=min(S(term(j).RGIshaperef).Y) & dem.y<=max(S(term(j).RGIshaperef).Y)),find(dem.x>=min(S(term(j).RGIshaperef).X) & dem.x<=max(S(term(j).RGIshaperef).X)))))]); %stretch colors in DEM
        figure(DEM_fig); subplot(sub2);
        plot(S(term(j).RGIshaperef).X,S(term(j).RGIshaperef).Y,'-r','linewidth',2); hold on;
        set(gca,'xlim',[min([S(term(j).RGIshaperef).X term(j).X]) max([S(term(j).RGIshaperef).X term(j).X])],'ylim',[min([S(term(j).RGIshaperef).Y term(j).Y]) max([S(term(j).RGIshaperef).Y term(j).Y])]);
        drawnow;
        
        %crop the velocities to just the RGI outline
        disp('Prepping velocity vectors to aid mapping of the flux gate at the grounding line, perpendicular to flow...');
        vx_sub = V.vx(find(V.y>=min(im.y) & V.y<=max(im.y)),find(V.x>=min(im.x) & V.x<=max(im.x))); %crop vx to the landsat image extent
        vy_sub = V.vy(find(V.y>=min(im.y) & V.y<=max(im.y)),find(V.x>=min(im.x) & V.x<=max(im.x))); %crop vy to the landsat image extent
        x_sub = V.x(1,find(V.x>=min(im.x) & V.x<=max(im.x)));
        y_sub = V.y(1,find(V.y>=min(im.y) & V.y<=max(im.y)));
        [Xsubgrid,Ysubgrid] = meshgrid(x_sub,y_sub);
        in = inpolygon(Xsubgrid,Ysubgrid,[term(j).RGIX(1:end-1) term(j).RGIX(1)],[term(j).RGIY(1:end-1) term(j).RGIY(1)]); %mask velocities outside the glacier so velocity arrows scale approrpriately
        maskedvx = vx_sub; maskedvx(in==0) = NaN; %mask-out other regions by turning anything not in=1 into NaNs
        maskedvy = vy_sub; maskedvy(in==0) = NaN; %mask-out other regions by turning anything not in=1 into NaNs
        
        %plot the DEM with the velocities overlain
        GL_fig = figure; set(GL_fig,'position',[550 50 700 700]);
        cd_to_arcticDEM = ['cd ',ArcticDEMdir]; eval(cd_to_arcticDEM);
        dem_xsub = [find(ArcticDEMx<=min([S(term(j).RGIshaperef).X term(j).X]),1,'last'),find(ArcticDEMx>=max([S(term(j).RGIshaperef).X term(j).X]),1,'first')];
        dem_ysub = [find(ArcticDEMy>=max([S(term(j).RGIshaperef).Y term(j).Y]),1,'last'),find(ArcticDEMy<=min([S(term(j).RGIshaperef).Y term(j).Y]),1,'first')];
        I = imread('ArcticDEM_10m_Greenland_Mosaic.tif','PixelRegion',{dem_ysub,dem_xsub}); I(I==-9999) = NaN;
        ArcticDEMz = I; clear I;
        [ArcticDEM_xgrid,ArcticDEM_ygrid] = meshgrid(ArcticDEMx(min(dem_xsub):max(dem_xsub)),ArcticDEMy(min(dem_ysub):max(dem_ysub)));
        [lon,lat] = ps2wgs(term(j).centerX,term(j).centerY); geoid = geoidheight(lat,lon);
%         [ArcticDEM_xgrid,ArcticDEM_ygrid] = meshgrid(ArcticDEMx(min(dem_xsub):max(dem_xsub)),ArcticDEMy(min(dem_ysub):max(dem_ysub)));
        imagesc(double(ArcticDEMx(min(dem_xsub):max(dem_xsub))),double(ArcticDEMy(min(dem_ysub):max(dem_ysub))),ArcticDEMz-nanmedian(geoid)); axis xy equal; colormap gray; hold on;
        set(gca,'clim',[0 min([term(j).centerGL(3) 100])+20]); cbar = colorbar; cbar.Label.String = 'elevation (m a.s.l.)';
        set(gca,'xlim',[min([term(j).RGIX term(j).X]) max([term(j).RGIX term(j).X])],'ylim',[min([term(j).RGIY term(j).Y]) max([term(j).RGIY term(j).Y])],'fontsize',14);
        q = quiver(Xsubgrid,Ysubgrid,maskedvx,maskedvy,'-r'); hold on; %q.AutoScaleFactor = 50;
        plot(term(j).RGIX,term(j).RGIY,'-k','linewidth',2); %plot the RGI outline
        plot(term(j).X,term(j).Y,'-m','linewidth',2); %plot the terminus box
        plot(term(j).centerX,term(j).centerY,'--b','linewidth',2); %plot the centerline
        plot(term(j).centerGL(1),term(j).centerGL(2),'xc','linewidth',2,'markersize',12); %plot the centerline flux gate location
        drawnow;
        
        %manually zoom in & draw a flux gate spanning the
        %glacier width that is approx. perpendicular to flow
        disp('click on UL & LR corners of box bounding the area near the grounding line spanning the glacier width');
        figure(GL_fig); [b] = ginput(2);
        set(gca,'xlim',[min(b(:,1)) max(b(:,1))],'ylim',[min(b(:,2)) max(b(:,2))]); clear b;
        disp('draw a gate across the RGI outline, crossing the center grounding line, semi-perpendicular to flow');
        [~,~,~,xi,yi] = improfile;
        %interpolate profile to the velocity grid
        [vx,vy,~,~,~] = improfile(double(V.x),double(V.y),sqrt(V.vx.^2 + V.vy.^2),xi,yi,'nearest');
        vprof = interp2(Xsubgrid,Ysubgrid,sqrt(vx_sub.^2 + vy_sub.^2),vx',vy');
        [lon,lat] = ps2wgs(vx,vy); geoidprof = geoidheight(lat,lon); clear lat lon;
        if size(interp2(demx,demy,dem.z,vx',vy'),2) ~= size(geoidprof,2)
            Zprof(1,:) = interp2(demx,demy,dem.z,vx',vy')-geoidprof';
            Zprof(2,:) = interp2(ArcticDEM_xgrid,ArcticDEM_ygrid,ArcticDEMz,vx',vy')-geoidprof';
        else
            Zprof(1,:) = interp2(demx,demy,dem.z,vx',vy')-geoidprof;
            Zprof(2,:) = interp2(ArcticDEM_xgrid,ArcticDEM_ygrid,ArcticDEMz,vx',vy')-geoidprof;
        end
        term(j).gateX = vx'; term(j).gateY = vy'; term(j).gateV = vprof; term(j).gateZ = Zprof;
        clear xi yi vx vy vprof Zprof;
        plot(term(j).gateX,term(j).gateY,'+c'); hold on;
        set(gca,'xlim',[min([term(j).RGIX term(j).X]) max([term(j).RGIX term(j).X])],'ylim',[min([term(j).RGIY term(j).Y]) max([term(j).RGIY term(j).Y])],'fontsize',14);
        xlabel('Easting (m)','fontsize',14); ylabel('Northing (m)','fontsize',14);
        %save the map
        if j<10
            saveas(GL_fig,['00',num2str(j),'_fluxgate-map.png'],'png');
        elseif j >=10 && j<100
            saveas(GL_fig,['0',num2str(j),'_fluxgate-map.png'],'png');
        else
            saveas(GL_fig,[num2str(j),'_fluxgate-map.png'],'png');
        end
        
        %save the data to appropriate places
        save('Greenland_GIC_centerlines.mat','term','-v7.3');
        cd_to_fluxgates = ['cd ',root_dir,'fluxgates/']; eval(cd_to_fluxgates);
        s.Geometry = 'Polyline';
        s.BoundingBox = double([min(term(j).gateX) min(term(j).gateY); max(term(j).gateX) max(term(j).gateY)]);
        s.X = double(term(j).gateX); s.Y = double(term(j).gateY);
        s.RGIId = S(term(j).RGIshaperef).RGIId;
        s.BoxID = term(j).BoxID;
        if j<10
            shapewrite(s,['fluxgate_00',num2str(j),'.shp']);
        elseif j >=10 && j<100
            shapewrite(s,['fluxgate_0',num2str(j),'.shp']);
        else
            shapewrite(s,['fluxgate_',num2str(j),'.shp']);
        end
        clear s;
        
        %close site-specific figures, zoom out in overview figure, clear variables, and advance
        close(GL_fig);
        figure(DEM_fig); %subplot(sub1); imagesc(dem.x,dem.y,dem.z); axis xy equal; hold on; colormap(sub1,dem_cmap);
        set(sub1,'xlim',[min(im.x) max(im.x)],'ylim',[min(im.y) max(im.y)]);
        set(sub2,'xlim',[min(im.x) max(im.x)],'ylim',[min(im.y) max(im.y)]);
        drawnow;
        clear prof* s lon lat geoid minref maxref *intercept* *_sub *subgrid in maskedv* q RGI*lims geoid dem_*sub ArcitcDEMz ArcticDEM_*grid;
        disp('Done!');
        disp(' ');
        
    end
    
end





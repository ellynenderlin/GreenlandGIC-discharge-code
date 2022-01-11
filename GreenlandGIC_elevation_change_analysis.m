%%% Estimate Greenland GIC thickness change due to surface accumulation and
%%% melting between the flux gate and average terminus position, and over
%%% time at the flux gate (due to SMB and dynamics)

disp('This code combines several old codes so there may be some hiccups running it through!');

%loop through all the elevation data for each study site to locate the grounding line
clearvars; close all; warning off;
addpath('/users/ellynenderlin/mfiles','/users/ellynenderlin/mfiles/general');

%specify the root directory for project-specific files
root_dir = '/Users/ellynenderlin/Research/NASA_GreenlandPeriph-Mapping/'; %including trailing /
%specify the root directory for generic files
misc_dir = '/Users/ellynenderlin/Research/miscellaneous/'; %including trailing /

%load term structure
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
load Greenland_GIC_centerlines.mat;

%set up days of year for leap and non-leap years for date conversions
modays = [31 28 31 30 31 30 31 31 30 31 30 31]; cumdays = [0 cumsum(modays(1:end-1))];
leap_modays = [31 29 31 30 31 30 31 31 30 31 30 31]; leap_cumdays = [0 cumsum(leap_modays(1:end-1))];

%% identify up-glacier & down-glacier edges of the terminus boxes
%load the RGI shapefiles
cd_to_outlines = ['cd ',root_dir,'RGI-outlines/']; eval(cd_to_outlines);
S = shaperead('05_rgi50_GreenlandPeriphery_PS_BoxIDs_TW.shp');
for i = 1:length(S)
    RGIBoxID(i) = S(i).BoxID;
end

%load the centerline data structure
cd_to_centerlines = ['cd ',root_dir,'centerlines/']; eval(cd_to_centerlines);
% load Greenland_GIC_centerlines.mat;
    
%load the manual terminus delineations
cd_to_termini = ['cd ',root_dir,'terminus_traces/']; eval(cd_to_termini);
T = shaperead('delineations_2000.shp');
for i = 1:length(T)
        BoxID(i) = T(i).BoxID;
end

%find where terminus traces intersect box walls
for j = 1:length(term)
    termref = find(BoxID==term(j).BoxID); 

    %loop around the terminus box looking for intersections with the delineation
    siderefs = [];
    for k = 1:4
        [xi,yi] = polyxpoly(term(j).X(k:k+1),term(j).Y(k:k+1),T(termref).X,T(termref).Y);
        if ~isempty(xi)
            siderefs = [siderefs k]; xmid(k) = NaN; ymid(k) = NaN; 
        else
            xmid(k) = nanmean(term(j).X(k:k+1)); ymid(k) = nanmean(term(j).Y(k:k+1)); 
        end
        
        clear xi yi;
    end
    
    %identify the inland edge as the non-side closest to the RGI centroid
    if isempty(term(j).RGIX)
        term(j).RGIX = S(term(j).RGIshaperef).X; term(j).RGIY = S(term(j).RGIshaperef).Y; 
    end
    polyin = polyshape(term(j).RGIX,term(j).RGIY); [RGIX,RGIY] = centroid(polyin);
%     [~,minref] = min(sqrt((xmid-RGIX).^2 + (ymid-RGIY).^2)); 
    [sortedvals,sortrefs] = sort(sqrt((xmid-RGIX).^2 + (ymid-RGIY).^2));
    inlandX(j,:) = term(j).X(sortrefs(1):sortrefs(1)+1); inlandY(j,:) = term(j).Y(sortrefs(1):sortrefs(1)+1); 
    seawardX(j,:) = term(j).X(sortrefs(2):sortrefs(2)+1); seawardY(j,:) = term(j).Y(sortrefs(2):sortrefs(2)+1); 
    clear polyin;
    
    %buffer the inland & seaward lines to form boxes
    buffer = 10; %meters on either side of box edge
    edge_angle = atan2d((inlandY(j,2)-inlandY(j,1)),(inlandX(j,2)-inlandX(j,1)));
    %inland first
    if edge_angle == 0 || edge_angle == 180 || edge_angle == -180
        bufferX = [inlandX(j,1) inlandX(j,2) inlandX(j,2) inlandX(j,1)];
        bufferY = [inlandY(j,1)+buffer inlandY(j,2)+buffer inlandY(j,2)-buffer inlandY(j,1)-buffer];
    elseif edge_angle == 90 || edge_angle == -90
        bufferX = [inlandX(j,1)+buffer inlandX(j,2)+buffer inlandX(j,2)-buffer inlandX(j,1)-buffer];
        bufferY = [inlandY(j,1) inlandY(j,2) inlandY(j,2) inlandY(j,1)];
    else
        bufferX = [inlandX(j,1)+buffer*sind(edge_angle+90) inlandX(j,2)+buffer*sind(edge_angle+90) inlandX(j,2)-buffer*sind(edge_angle+90) inlandX(j,1)-buffer*sind(edge_angle+90)];
        bufferY = [inlandY(j,1)+buffer*cosd(edge_angle+90) inlandY(j,2)+buffer*cosd(edge_angle+90) inlandY(j,2)-buffer*cosd(edge_angle+90) inlandY(j,1)-buffer*cosd(edge_angle+90)];
    end
    inlandbuffer = polyshape(bufferX,bufferY); 
    termbox = polyshape(term(j).X(~isnan(term(j).X)),term(j).Y(~isnan(term(j).Y)));
    inlandbox = intersect(inlandbuffer,termbox);
    term(j).inlandX = inlandbox.Vertices(:,1); term(j).inlandY = inlandbox.Vertices(:,2); 
    clear bufferX bufferY edge_angle;
    %seaward second
    edge_angle = atan2d((seawardY(j,2)-seawardY(j,1)),(seawardX(j,2)-seawardX(j,1)));
    if edge_angle == 0 || edge_angle == 180 || edge_angle == -180
        bufferX = [seawardX(j,1) seawardX(j,2) seawardX(j,2) seawardX(j,1)];
        bufferY = [seawardY(j,1)+buffer seawardY(j,2)+buffer seawardY(j,2)-buffer seawardY(j,1)-buffer];
    elseif edge_angle == 90 || edge_angle == -90
        bufferX = [seawardX(j,1)+buffer seawardX(j,2)+buffer seawardX(j,2)-buffer seawardX(j,1)-buffer];
        bufferY = [seawardY(j,1) seawardY(j,2) seawardY(j,2) seawardY(j,1)];
    else
        bufferX = [seawardX(j,1)+buffer*sind(edge_angle+90) seawardX(j,2)+buffer*sind(edge_angle+90) seawardX(j,2)-buffer*sind(edge_angle+90) seawardX(j,1)-buffer*sind(edge_angle+90)];
        bufferY = [seawardY(j,1)+buffer*cosd(edge_angle+90) seawardY(j,2)+buffer*cosd(edge_angle+90) seawardY(j,2)-buffer*cosd(edge_angle+90) seawardY(j,1)-buffer*cosd(edge_angle+90)];
    end
    seawardbuffer = polyshape(bufferX,bufferY); 
    termbox = polyshape(term(j).X(~isnan(term(j).X)),term(j).Y(~isnan(term(j).Y)));
    seawardbox = intersect(seawardbuffer,termbox);
    term(j).seawardX = seawardbox.Vertices(:,1); term(j).seawardY = seawardbox.Vertices(:,2); 
    clear *buffer* *box edge_angle;
    
    %calculate the box width
    IDW(j) = term(j).BoxID;
    BoxW(j) = sqrt((term(j).X(sortrefs(1))-term(j).X(sortrefs(1)+1))^2 + (term(j).Y(sortrefs(1))-term(j).Y(sortrefs(1)+1))^2);
    disp(['BoxID=',num2str(term(j).BoxID),': sides=',num2str(siderefs),' & end=',num2str(sortrefs(1)),' with width=',num2str(BoxW(j)),'m']);
    
    %plot everything and zoom in on the terminus
    close all; figure; plot(term(j).RGIX,term(j).RGIY,'-k','linewidth',1.5); hold on;
    plot(term(j).X,term(j).Y,'-r','linewidth',1.5); hold on;
    plot([term(j).inlandX; term(j).inlandX(1)],[term(j).inlandY; term(j).inlandY(1)],'-b','linewidth',2); hold on;
    plot([term(j).seawardX; term(j).seawardX(1)],[term(j).seawardY; term(j).seawardY(1)],'-c','linewidth',2); hold on;
    set(gca,'xlim',[min(term(j).X) max(term(j).X)],'ylim',[min(term(j).Y) max(term(j).Y)]); 
    drawnow;
    
    clear termref xmid ymid siderefs RGIX RGIY minref;
end
cd_to_centerlines = ['cd ',root_dir,'centerlines/']; eval(cd_to_centerlines);
save('Greenland_GIC_centerlines.mat','term','-v7.3');

%export to a text files
sorted_widths = sortrows([IDW' BoxW'],1);
cd ..
csvwrite('Greenland_GIC_widths.csv',sorted_widths);
csvwrite('Greenland_GIC_inside-coords.csv',sortrows([IDW' inlandX inlandY]));
csvwrite('Greenland_GIC_outside-coords.csv',sortrows([IDW' seawardX seawardY]));

clear S IDW BoxW *box *buffer;

%% estimate thickness change between the flux gate and terminus

%load the RACMO SMB data (1958-2019; 62 years)
cd_to_RACMO = ['cd ',misc_dir,'RACMO_GreenlandPeriph/']; eval(cd_to_RACMO);
smb = ncread('SMB_rec.1958-2019.BN_RACMO2.3_ZGRN11.1km.YY.nc', 'SMB_rec'); smb(smb==-9999)=NaN;
lat = ncread('SMB_rec.1958-2019.BN_RACMO2.3_ZGRN11.1km.YY.nc', 'LAT');
lon = ncread('SMB_rec.1958-2019.BN_RACMO2.3_ZGRN11.1km.YY.nc', 'LON');
%convert from geographic coordniates to polar stereo
[x,y] = wgs2ps(lon,lat,'StandardParallel',70,'StandardMeridian',-45);

%loop through term structure & save BoxIDs to a vector
for i = 1:length(term)
    BoxID(i) = term(i).BoxID;
end

%load all the terminus position delineations
cd_to_termini = ['cd ',root_dir,'terminus_traces/']; eval(cd_to_termini);
Sold = shaperead('ManTermDelins_1985.shp');
for i = 1:length(Sold)
    term1ID(i) = Sold(i).BoxID;
end
S = shaperead('ManTermDelins_2000.shp');
for i = 1:length(S)
    term2ID(i) = S(i).BoxID;
end
Snew = shaperead('ManTermDelins_2015.shp');
for i = 1:length(Snew)
    term3ID(i) = Snew(i).BoxID;
end

%set-up a figure for plotting to check locations
figure;
imagesc(nanmean(smb,3)); hold on; %avg SMB over time
set(gca,'clim', [-8000 4000]);
cmap = [colormap(hot(8000)); colormap(cool(4000))]; colormap(cmap); cbar = colorbar;

%loop through all glaciers
for j = 1:length(BoxID)
    %find corresponding terminus from ~1985 (1980s-1990s)
    ID = find(term1ID == BoxID(j));
    %find center point of manual terminus delineation
    if ~isempty(ID)
        [t1ctrx,t1ctry] = polyxpoly(Sold(ID).X,Sold(ID).Y,term(j).centerX,term(j).centerY);
        %take avg if there are multiple intersections for a wavy terminus
        %or centerline somehow does not cross terminus trace
        if length(t1ctrx)>1 || isempty(t1ctrx) 
            t1ctrx = nanmean(t1ctrx); t1ctry = nanmean(t1ctry); 
        end 
    else
        t1ctrx = NaN; t1ctry = NaN; 
    end
    clear ID;
    
    %find corresponding terminus from 2000
    ID = find(term2ID == BoxID(j));
    %find center point of manual terminus delineation
    if ~isempty(ID)
        [t2ctrx,t2ctry] = polyxpoly(S(ID).X,S(ID).Y,term(j).centerX,term(j).centerY);
        %take avg if there are multiple intersections for a wavy terminus
        %or centerline somehow does not cross terminus trace
        if length(t2ctrx)>1 || isempty(t2ctrx)
            t2ctrx = nanmean(t2ctrx); t2ctry = nanmean(t2ctry);
        end
    else
        t2ctrx = NaN; t2ctry = NaN;
    end
    clear ID;
    
    %find corresponding terminus from 2015
    ID = find(term3ID == BoxID(j));
    %find center point of manual terminus delineation
    if ~isempty(ID)
        [t3ctrx,t3ctry] = polyxpoly(Snew(ID).X,Snew(ID).Y,term(j).centerX,term(j).centerY);
        if length(t3ctrx)>1 || isempty(t3ctrx)
            t3ctrx = nanmean(t3ctrx); t3ctry = nanmean(t3ctry);
        end
    else
        t3ctrx = NaN; t3ctry = NaN;
    end
    clear ID;
    
    %loop through fluxgate centroids
    [centX,centY] = polyxpoly(term(j).gateX,term(j).gateY,term(j).centerX,term(j).centerY);
    if isempty(centX)
        centX = nanmean(term(j).gateX);
    end
    if isempty(centY)
        centY = nanmean(term(j).gateY);
    end
    
    %calculate distance between fluxgate centroid and man delin centroid
    %distance to 1985 terminus
    dist1 = sqrt((centX - t1ctrx).^2 + (centY - t1ctry).^2);
    %distance to 2000 terminus
    dist2 = sqrt((centX - t2ctrx).^2 + (centY - t2ctry).^2);
    %distance to 2015 terminus
    dist3 = sqrt((centX - t3ctrx).^2 + (centY - t3ctry).^2);
    %grab average date for gate velocity
    date = term(j).gateVdateavg; vel = term(j).fluxVavg;
    
    %estimate travel times
    %find travel time (yrs) to 1985 term
    ttime1 = dist1/nanmean(vel(1:find(date<=1999,1,'last')));
    %find travel time (yrs) to 2000 term
    ttime2 = dist2/nanmean(vel(find(date>1999,1,'first'):find(date<=2005,1,'last')));
    %find travel time (yrs) to 2015 term
    ttime3 = dist3/nanmean(vel(find(date>2005,1,'first'):find(date<=2015,1,'last')));
    %check that travel times aren't outrageous... if they are then
    %recalculate with avg velocity over the full time period
    if ttime1 > max(date)-1999; ttime2 = max(date)-1999; end
    if ttime2 > max(date)-2005; ttime2 = max(date)-2005; end
    if ttime3 > max(date)-2015; ttime3 = max(date)-2015; end
    
    %save important measurements
    Dist(j,:) = [dist1 dist2 dist3];
    cent(j,:) = [centX centY];
    ttime(j,:) = [ttime1 ttime2 ttime3];
    
    %clear variables
    clear centX centY dist* ID t3* t2* ttime2 ttime3;
    
    
    %loop through centroids and find nearest RACMO value
    diff_map = sqrt((x-cent(j,1)).^2+(y-cent(j,2)).^2); %solve for the distance vector using the x&y distances
    diff_map(nanmean(smb,3)==0) = NaN;
    RACMO_ref(j) = find(diff_map==min(min(diff_map))); %find the minimum distance (reference for your grid cell is output)
    [RACMOy(j),RACMOx(j)] = ind2sub(size(nanmean(smb,3)),RACMO_ref(j)); %convert cell reference to an x- and y-cell index
    term(j).SMB = squeeze(smb(RACMOy(j),RACMOx(j),:));
    term(j).SMByrs = [1958:1:2019]';
    clear diff_map;
    
    %plot to make sure it worked
    plot(RACMOx(j),RACMOy(j),'ok','markerfacecolor','k','markeredgecolor','w','markersize',8); hold on;
%     plot(RACMOx(j),RACMOy(j),'xw','markersize',8); hold on;
    drawnow;
    
    %estimate thickness change due to accumulation or melt between the flux
    %gate & the terminus for each velocity observation
    for m = 1:length(date) %for each velocity date
        if ~isnan(date(m))
            date_closest = find(term(j).SMByrs<=date(m),1,'last'); %date that ice crosses fluxgate
            if date(m)<= 1999
                %find nearest RACMO date and add travel time
                %(date_farthest=time of discharge)
                date_farthest = find(term(j).SMByrs<=date(m)+ttime(j,1),1,'last');
                %sum SMB over the time between the fluxgate and the terminus
                %(SMBint = SMB integrated over time = mm water equiv)
                if date_farthest-date_closest >= 2
                    term(j).SMBint(m,1) = nansum([(ceil(date(m))-date(m))*(term(j).SMB(date_closest));
                        squeeze(term(j).SMB(date_closest+1:date_farthest-1));
                        ((date(m)+ttime(j,1))-floor((date(m)+ttime(j,1))))*term(j).SMB(date_farthest)]);
                else
                    term(j).SMBint(m,1) = nansum([(ceil(date(m))-date(m))*(term(j).SMB(date_closest));
                        ((date(m)+ttime(j,1))-floor((date(m)+ttime(j,1))))*term(j).SMB(date_farthest)]);
                end
            elseif date(m)>1999 && date(m)<=2005
                %find nearest RACMO date and add travel time
                %(date_farthest=time of discharge)
                date_farthest = find(term(j).SMByrs<=date(m)+ttime(j,2),1,'last');
                %sum SMB over the time between the fluxgate and the terminus
                %(SMBint = SMB integrated over time = mm water equiv)
                if date_farthest-date_closest >= 2
                    term(j).SMBint(m,1) = nansum([(ceil(date(m))-date(m))*(term(j).SMB(date_closest));
                        squeeze(term(j).SMB(date_closest+1:date_farthest-1));
                        ((date(m)+ttime(j,2))-floor((date(m)+ttime(j,2))))*term(j).SMB(date_farthest)]);
                else
                    term(j).SMBint(m,1) = nansum([(ceil(date(m))-date(m))*(term(j).SMB(date_closest));
                        ((date(m)+ttime(j,2))-floor((date(m)+ttime(j,2))))*term(j).SMB(date_farthest)]);
                end
            elseif date(m)>2005 && date(m)<=2015
                date_farthest = find(term(j).SMByrs<=date(m)+(ttime(j,2)+ttime(j,3))/2,1,'last');
                if date_farthest-date_closest >= 2
                    term(j).SMBint(m,1) = nansum([(ceil(date(m))-date(m))*(term(j).SMB(date_closest));
                        squeeze(term(j).SMB(date_closest+1:date_farthest-1));
                        ((date(m)+((ttime(j,2)+ttime(j,3))/2))-floor((date(m)+((ttime(j,2)+ttime(j,3))/2))))*term(j).SMB(date_farthest)]);
                else
                    term(j).SMBint(m,1) = nansum([(ceil(date(m))-date(m))*(term(j).SMB(date_closest));
                        ((date(m)+ttime(j,2))-floor((date(m)+ttime(j,2))))*term(j).SMB(date_farthest)]);
                end
            else
                date_farthest = find(term(j).SMByrs<=date(m)+ttime(j,3),1,'last');
                if date_farthest-date_closest >= 2
                    term(j).SMBint(m,1) = nansum([(ceil(date(m))-date(m))*(term(j).SMB(date_closest));
                        squeeze(term(j).SMB(date_closest+1:date_farthest-1));
                        ((date(m)+ttime(j,3))-floor((date(m)+ttime(j,3))))*term(j).SMB(date_farthest)]);
                else
                    term(j).SMBint(m,1) = nansum([(ceil(date(m))-date(m))*(term(j).SMB(date_closest));
                        ((date(m)+ttime(j,3))-floor((date(m)+ttime(j,3))))*term(j).SMB(date_farthest)]);
                end
            end
            clear date_closest date_farthest;
        else
            term(j).SMBint(m,1) = NaN;
        end
    end
    clear date
    
    %convert SMB between the flux gate & terminus to estimated thinning
    term(j).dH_SMB = (term(j).SMBint/1000)*(1000/917); %mm*(m/mm)*(density water/density ice)
    
    %estimate discharge change between flux gate and terminus due to
    %surface accumulation & ablation between the gate and terminus
    term(j).dD_SMB = nansum(repmat(term(j).dH_SMB,1,size(term(j).fluxH,2)).*term(j).fluxV.*repmat(term(j).fluxW,size(term(j).fluxD,1),1),2);
    term(j).dD_SMB(term(j).dD_SMB>term(j).fluxD) = term(j).fluxD(term(j).dD_SMB>term(j).fluxD); %surface melt cannot cause negative discharge
end

%save the structure
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
save('Greenland_GIC_centerlines.mat','term','-v7.3');
disp('RACMO SMB-derived thickness change added to structure');

clear S* ID dist* diff* date*;

%% extract elevation change at the flux gate
%navigate to the matlab structure with the flux gates
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
% load Greenland_GIC_centerlines.mat;
for i = 1:length(term)
    BoxID(i) = term(i).BoxID;
end

%navigate to and load the GIMP image mosaic for checking referencing of the
%various datasets
cd_to_GIMP = ['cd ',misc_dir,'GIMP/']; eval(cd_to_GIMP);
load GIMP_image_mosaic_150m.mat;
bigfig = figure; set(bigfig,'position',[50 50 1200 1200]);
imagesc(I.x,I.y,I.z); colormap gray; axis xy equal; hold on;

%navigate to the DEM directory
% cd_to_PGCDEMs = ['cd ',misc_dir,'PGC_DEMs/']; eval(cd_to_PGCDEMs);
cd /Users/katebollen/Documents/MS/data/DEMs/PGC_DEMs/
files = dir;
for i = 1:length(files)
    if files(i).isdir && isempty(strfind(files(i).name,'.'))
        termref = find(BoxID == str2num(files(i).name));
        
        %identify the DEMs & loop through them to extract elevations
        cd_to_dir = ['cd ',files(i).name]; eval(cd_to_dir);
        DEMs = dir('*dem*.tif'); 
        for j = 1:length(DEMs)
            %extract the date
            datestring = DEMs(j).name(6:13);
            if mod(str2num(datestring(1:4)),4) ~= 0
                elev_date(j) = str2num(datestring(1:4)) + ((cumdays(str2num(datestring(5:6)))+str2num(datestring(7:8)))./sum(modays));
            else
                elev_date(j) = str2num(datestring(1:4)) + ((leap_cumdays(str2num(datestring(5:6)))+str2num(datestring(7:8)))./sum(leap_modays));
            end
            
            %load the geotiff
            [A,R] = readgeoraster(DEMs(j).name); A(A==-9999) = NaN;
            
            %convert geospatial referencing info to vectors
            if strcmp(R.ColumnsStartFrom,'north')
                y = R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY;
            else
                y = R.YWorldLimits(1)+0.5*R.CellExtentInWorldY:R.CellExtentInWorldY:R.YWorldLimits(2)-0.5*R.CellExtentInWorldY;
            end
            if strcmp(R.RowsStartFrom,'west')
                x = R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX;
            else
                x = R.XWorldLimits(2)-0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(1)+0.5*R.CellExtentInWorldX;
            end

            %identify pixels in the DEM that are within the box edge 10m buffer
            [xgrid,ygrid] = meshgrid(x,y);
%             in = inpolygon(xgrid,ygrid,S(Sref).X,S(Sref).Y);
            inland = inpolygon(xgrid,ygrid,term(termref).inlandX,term(termref).inlandY);
            seaward = inpolygon(xgrid,ygrid,term(termref).seawardX,term(termref).seawardY);
         
            %compute the median elevation for the pixels in the box
            elev_medi(j) = nanmedian(A(inland)); elev_madi(j) = mad(A(inland),1);
            elev_meds(j) = nanmedian(A(seaward)); elev_mads(j) = mad(A(seaward),1);
            
            %check that the references for S, term, and the BoxID are right
            if j==1
                figure(bigfig); 
                plot(term(termref).X,term(termref).Y,'-','color',[255,255,191]/255,'linewidth',2); hold on;
                plot(term(termref).inlandX,term(termref).inlandY,'-','color',[44,123,182]/255,'linewidth',1.5); hold on;
                plot(term(termref).seawardX,term(termref).seawardY,'-','color',[253,174,97]/255,'linewidth',1.5); hold on;
                set(gca,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)]); drawnow;
            end
            
            %clear the variables
            clear datestring A R y x xgrid ygrid inland seaward;
        end
        
        %add the data to the term structure
        [sorted,sortrefs] = sort(elev_date);
        term(termref).inlandZmed = elev_medi(sortrefs); term(termref).inlandZmad = elev_madi(sortrefs); term(termref).inlandZyrs = sorted;
        term(termref).seawardZmed = elev_meds(sortrefs); term(termref).seawardZmad = elev_mads(sortrefs); term(termref).seawardZyrs = sorted;
        clear sorted sortrefs;
        
        %create elevation change time series and save
        elev_ts = figure; 
        plot(term(termref).inlandZyrs,term(termref).inlandZmed,'-','color',[44,123,182]/255,'linewidth',2); hold on;
        fill([term(termref).inlandZyrs fliplr(term(termref).inlandZyrs)],[term(termref).inlandZmed+term(termref).inlandZmad fliplr([term(termref).inlandZmed-term(termref).inlandZmad])],[44,123,182]/255,'facealpha',0.5,'edgecolor','none');
        if nanmedian(term(termref).inlandZmed) > nanmedian(term(termref).seawardZmed)
            term(termref).Zflag = 1; %first side identified as the inland box is correct
        else
            term(termref).Zflag = 2; %second side identified as the potential inland box is correct
        end
        grid on; set(gca,'fontsize',14);
        xlabel('Year ','fontsize',14); ylabel('Elevation (m)','fontsize',14); drawnow;
        cd ..
        saveas(elev_ts,['GIC_BoxID-',num2str(term(termref).BoxID),'_elev_timeseries.png'],'png');
        close(elev_ts);
        
        disp(['Elevation data extracted for DEM set #',files(i).name,' (term ref=',num2str(termref),')']);
        clear termref elev_* DEMs;
    end
end
disp('Elevation time series extracted for ALL DEMs!!!');

%save the structure
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
save('Greenland_GIC_centerlines.mat','term','-v7.3');


%% compare time series of velocity & thickness change at the flux gate
close all; clear all;

%load term structure
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
% load Greenland_GIC_centerlines.mat;

%set up figures
figure1=figure; set(gcf,'position',[100 100 800 500]);
figure2=figure; set(gcf,'position',[700 100 800 500]);

%loop through the data & compare speed & thickness changes over time
for i = 1:length(term)
    if term(i).MankoffFlag == 0
        if ~isempty(term(i).Zflag)
            %get rid of NaNs in elevation timeseries
            term(i).inlandZyrs(isnan(term(i).inlandZmed)) = []; term(i).inlandZmed(isnan(term(i).inlandZmed)) = [];
            term(i).seawardZyrs(isnan(term(i).seawardZmed)) = []; term(i).seawardZmed(isnan(term(i).seawardZmed)) = [];
            
            %choose seaward or inland edge for elevation timeseries: higher median
            %elevation (Zflag) = on-ice terminus box edge
            if term(i).Zflag == 1
                trueInland = term(i).inlandZmed;
                trueInlandyrs = term(i).inlandZyrs;
            else
                trueInland = term(i).seawardZmed;
                trueInlandyrs = term(i).seawardZyrs;
            end
            
            %skip glaciers with only one elevation value
            if length(trueInland(~isnan(trueInland))) > 1 && nansum(term(i).fluxVyrs(~isnan(term(i).fluxVyrs))) > 0
                
                %calculate median annual V across fluxgate and get difference
                %between each row (year)
                Vyrs = nanmedian(term(i).fluxVyrs,2); V = nanmedian(term(i).fluxV,2);
                Vyrs(isnan(V)) = []; V(isnan(V)) = [];
                %total chance in V across elevation change record
                [~,minref] = min(abs(trueInlandyrs(1)-Vyrs));
                totalVDelta = V(end)-V(minref);
                dtV = Vyrs(end)-Vyrs(minref);
                %change in V between years
                incVDelta = diff(V);
                dvrate = incVDelta./diff(Vyrs);
                
                %calculate change in H
                Hyrs = trueInlandyrs; H = trueInland;
                Hyrs(isnan(H)) = []; H(isnan(H)) = [];
                % get deltaH between first and last observation
                totalHDelta = H(end)-H(1);
                dtH = Hyrs(end)-Hyrs(1);
                % get deltaH between each observation (each DEM may be from the
                % same year or different years)
                incHDelta = diff(H); incHDelta(abs(incHDelta)>nanmedian(term(i).fluxH)) = NaN;
                dhrate = incHDelta./diff(Hyrs); dhrate(abs(dhrate)>nanmedian(term(i).fluxH)) = NaN;
                %         if max(dhrate)>100
                %             disp(['check ',num2str(i)]); keyboard
                %             for j = 1:length(dhrate)-1
                %                if abs(dhrate(j)) > 100 && abs(dhrate(j+1)) > 100
                %                   H(j+1) = NaN;
                %                end
                %             end
                %         end
                dhyrs = (Hyrs(1:end-1)+Hyrs(2:end))/2;
                
                %estimate rate of elevation change with changes in speed
                [~,minref] = min(abs(Hyrs(1)-Vyrs));
                for k = 1:size(incHDelta,2)
                    [~,earlyref] = min(abs(Hyrs(k)-Vyrs)); [~,lateref] = min(abs(Hyrs(k+1)-Vyrs));
                    
                    %if velocities are from different years, compute
                    %1) dH/dV between consecutive H observations
                    %2) dH/dV since the first H observation
                    if earlyref ~= lateref
                        dHdV(k) = incHDelta(k)./(V(lateref)-V(earlyref)); dt(k) = Hyrs(k+1)-Hyrs(k); dV(k) = V(lateref)-V(earlyref);
                        dHdVcum(k) = (H(k+1)-H(1))./(V(lateref)-V(minref)); dtcum(k) = Hyrs(k+1)-Hyrs(1); dVcum(k) = V(lateref)-V(minref); dHcum(k) = H(k+1)-H(1);
                    else
                        dHdV(k) = NaN; dt(k) = NaN; dV(k) = NaN;
                        dHdVcum(k) = (H(k+1)-H(1))./(V(lateref)-V(minref)); dtcum(k) = Hyrs(k+1)-Hyrs(1); dVcum(k) = V(lateref)-V(minref); dHcum(k) = H(k+1)-H(1);
                    end
                    clear earlyref lateref;
                end
                dtcum(isinf(dHdVcum)) = []; dVcum(isinf(dHdVcum)) = []; dHcum(isinf(dHdVcum)) = []; dHdVcum(isinf(dHdVcum)) = [];
                
                %         %plot thickness change over time
                %         figure(figure2);
                %         plot(dhyrs,dhrate,'xk','linewidth',2); hold on;
                
                %fractional change in H (tells us if H change is important)
                Hfrac = totalHDelta/nanmedian(term(i).fluxH);
                term(i).fluxHchange = Hfrac;
                figure(figure2); plot(nanmedian(term(i).fluxH),Hfrac,'xk','linewidth',2); hold on;
            end
            clear H* dt* dh* dH* V* dv*;
            clear glacier_* gate_* D Derr;
            clear trueInland trueInlandyrs totalVDelta incVDelta totalHDelta incHDelta minref;
            %     drawnow;
            
            
            %plot and fit change in velocity and RACMO to deltaH
            warning off;
            if term(i).Zflag == 1
                
                if length(term(i).inlandZmed) > 1
                    [fz, gofz] = fit(term(i).inlandZyrs', term(i).inlandZmed','poly1');
                    term(i).dZ_tsfit = fz; term(i).dZ_tsr2 = gofz.rsquare;
                    %compare elevation and velocity time series
                    if sum(~isnan(term(i).fluxVelavg)) > 0
                        for j = 1:length(term(i).inlandZyrs)
                            Vdiff = abs(term(i).inlandZyrs(j) - term(i).gateVdateavg);
                            Vref = find(Vdiff==min(Vdiff));
                            vel(j) = term(i).fluxVelavg(Vref);
                            clear Vdiff Vref;
                        end
                        [fz, gofz] = fit(vel', term(i).inlandZmed','poly1');
                        term(i).dZ_vfit = fz; term(i).dZ_vr2 = gofz.rsquare;
                    else
                        term(i).dZ_vfit = NaN; term(i).dZ_vr2 = NaN;
                    end
                    %compare elevation and SMB time series
                    for j = 1:length(term(i).inlandZyrs)
                        SMBref = find(term(i).SMByrs == floor(term(i).inlandZyrs(j)));
                        smb(j) = term(i).SMB(SMBref);
                        clear SMBref;
                    end
                    [fz, gofz] = fit(smb', term(i).inlandZmed','poly1');
                    term(i).dZ_smbfit = fz; term(i).dZ_smbr2 = gofz.rsquare;
                else
                    term(i).dZ_tsfit = NaN; term(i).dZ_tsr2 = NaN;
                    term(i).dZ_vfit = NaN; term(i).dZ_vr2 = NaN;
                    term(i).dZ_smbfit = NaN; term(i).dZ_smbr2 = NaN;
                end
            else
                if length(term(i).seawardZmed) > 1
                    [fz, gofz] = fit(term(i).seawardZyrs', term(i).seawardZmed','poly1');
                    term(i).dZ_tsfit = fz; term(i).dZ_tsr2 = gofz.rsquare;
                    %compare elevation and velocity time series
                    if sum(~isnan(term(i).fluxVelavg)) > 0
                        for j = 1:length(term(i).seawardZyrs)
                            Vdiff = abs(term(i).seawardZyrs(j) - term(i).gateVdateavg);
                            Vref = find(Vdiff==min(Vdiff));
                            vel(j) = term(i).fluxVelavg(Vref);
                            clear Vdiff Vref;
                        end
                        [fz, gofz] = fit(vel', term(i).seawardZmed','poly1');
                        term(i).dZ_vfit = fz; term(i).dZ_vr2 = gofz.rsquare;
                    else
                        term(i).dZ_vfit = NaN; term(i).dZ_vr2 = NaN;
                    end
                    %compare elevation and SMB time series
                    for j = 1:length(term(i).seawardZyrs)
                        SMBref = find(term(i).SMByrs == floor(term(i).seawardZyrs(j)));
                        smb(j) = term(i).SMB(SMBref);
                        clear SMBref;
                    end
                    [fz, gofz] = fit(smb', term(i).seawardZmed','poly1');
                    term(i).dZ_smbfit = fz; term(i).dZ_smbr2 = gofz.rsquare;
                else
                    term(i).dZ_tsfit = NaN; term(i).dZ_tsr2 = NaN;
                    term(i).dZ_vfit = NaN; term(i).dZ_vr2 = NaN;
                    term(i).dZ_smbfit = NaN; term(i).dZ_smbr2 = NaN;
                end
            end
            clear vel fz gofz;
            clear smb;
        end
    end
end

%save the structure
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
save('Greenland_GIC_centerlines.mat','term','-v7.3');


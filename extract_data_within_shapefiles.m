%extract Greenland peripheral glacier velocities from ITS_LIVE using the
%box shapefiles to determine which glaciers we should focus on
clear all; close all;

%navigate to the directory containing the velocities
cd /users/ellynenderlin/Research/NASA-GreenlandPeriph-Mapping/ITS_LIVE

%loop through the velocity data and concatenate to make a time series of
%velocity maps
disp('Compile ITS_LIVE velocities...');
vels = dir('GRE*.nc'); %locate the netCDF files containing the data
for i = 1:length(vels)
    disp(['extracting velocities for ',vels(i).name]);
    
    %extract year from file name
    velyr(i) = str2num(vels(i).name(11:14));
    
    %load the coordinates (only need to do once)
    if i == 1
        x = ncread(vels(i).name,'x'); y = ncread(vels(i).name,'y');
    end
    
    %extract velocities and concatenate (velocities in m/yr)
%     vx(:,:,i) = ncread(vels(i).name,'vx'); %x-component of velocity
%     vy(:,:,i) = ncread(vels(i).name,'vy'); %y-component of velocity
%     vx_err(:,:,i) = ncread(vels(i).name,'vx_err'); %uncertainty in x-component of velocity
%     vy_err(:,:,i) = ncread(vels(i).name,'vy_err'); %uncertainty in y-component of velocity
	v(:,:,i) = ncread(vels(i).name,'v'); %speed (velocity magnitude)
%     v_err(:,:,i) = ncread(vels(i).name,'v_err'); %uncertainty in speed
    
    %extract date for each pixel (days since 01 January, 0000)
    days(:,:,i) = ncread(vels(i).name,'date');
end
[x_grid,y_grid] = meshgrid(x,y); %create a coordinate grid


%loop through the files, pulling info from each box
disp('Loop through boxes and extract velocities...');
cd ../Boxes
boxes = dir('*.shp');
for i = 1:length(boxes)
    disp(['Box #',num2str(i),' of ',num2str(length(boxes)),':']);
    S = shaperead(boxes(i).name); %read the shapefile
    box_id(i) = i; box_x = S.X(~isnan(S.X)); box_y = S.Y(~isnan(S.Y)); %get rid on NaNs
    [in,on] = inpolygon(x_grid,y_grid,box_x,box_y); %identify pixels with centers that fall within the box or on its edge
    disp('extracting velocities...');
%     figure; set(gcf,'position',[50 50 1200 700]); cmap=colormap(parula(10001)); cmap(1,:)=[0 0 0]; %only uncomment if you want to plot maps
    for j = 1:size(v,3)
%         disp(num2str(velyr(j)));
        
        %format velocity data
        V = v(:,:,j); %pull out the velocities for the particular year
        V(V==-32767) = NaN; %replace -32767 (their non-data value) with NaN
        V(~in & ~on) = NaN; %replace everything outside the shapefile footprint with a NaN
        
        %uncomment 3 lines below if you want to plot
%         subplot(7,5,j);
%         imagesc(x,y,V); hold on;
%         set(gca,'xlim',[min(box_x) max(box_x)],'ylim',[min(box_y) max(box_y)]); drawnow; 

        %uncomment line below if you want to show the velocity range
%         disp(['V range = [',num2str(min(min(V))),' ',num2str(max(max(V))),']']);

        %calculate the fraction of pixels with centers in or on the box that
        %actually have data in them
        data(j) = sum(sum(~isnan(V)))./(sum(sum(in))+sum(sum(on)));
        medv(j) = nanmedian(V(~isnan(V)));
        clear V;
    end
    disp(['average fractional V coverage = ',num2str(nanmean(data))]);
    disp(['years with V coverage = ',num2str(sum(ceil(data)))]);
    
    %save the results
    x_v(i) = nanmean(S.X); y_v(i) = nanmean(S.Y);
    fraction_v(i) = nanmean(data);
    typical_v(i) = nanmean(medv);
    coverage_v(i) = sum(ceil(data));
    
    %clear the data for the particular box
    clear S box_x box_y in on data medv;
end

%read in the GIMP ice and ocean masks to have a good GPsraphic reference
%for the data
cd ..
[GIMPice,R] = geotiffread('GimpIceMask_90m_v1.1.tif');
GIMPx = R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX;
GIMPy = R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY;
[GIMPocean,~] = geotiffread('GimpOceanMask_90m_v1.1.tif');
GIMPland = int16(ones(size(GIMPice))); GIMPland = GIMPland-GIMPice-GIMPocean;
mapfig = figure; imagesc(GIMPx,GIMPy,GIMPland+2*GIMPice); axis xy equal; hold on; 
plot(x_v,y_v,'+k'); hold on


%divide the data into GPsraphic quadrants to ensure we sample all around
%the Greenland periphery (data coverage is probably best in SE in general)
NN_ind = find(y_v>=-1.00e6);
NW_ind = find(y_v<-1.00e6 & y_v>=-2.145e6 & x_v<=0);
SW_ind = find(y_v<-2.145e6 & x_v<=0);
SE_ind = find((y_v<-2.12e6 & x_v>=0) | (y_v>=-2.12e6 & y_v<-2.02e6 & x_v>=7.8e5));
NE_ind = [1:1:length(y_v)]; NE_ind([NN_ind NW_ind SW_ind SE_ind]) = [];

%identify the target glaciers in each region using the speed dataset
GPs(1).ind = NN_ind; GPs(2).ind = NW_ind; GPs(3).ind = SW_ind; GPs(4).ind = SE_ind; GPs(5).ind = NE_ind; %combine regional indices
for i = 1:5
    %pull out information for the region
    region_count = round(100*(length(GPs(i).ind)./length(typical_v))); %proportionate number of glaciers (out of 100)

    %sort the data according to the abundance of velocity observations
    [~,v_ind] = sort(coverage_v(GPs(i).ind),'descend');
    
    %isolate the best data
    region_refs = v_ind(1:region_count);
    region_boxes = box_id(GPs(i).ind(region_refs));
    region_x = x_v(GPs(i).ind(region_refs)); region_y = y_v(GPs(i).ind(region_refs));
    region_fracv = fraction_v(GPs(i).ind(region_refs));
    region_normv = typical_v(GPs(i).ind(region_refs));
    region_coverv = coverage_v(GPs(i).ind(region_refs));
    
    %plot histograms showing the range of speeds covered by the sample
    figure;
    histogram(typical_v(GPs(i).ind),[0:10:ceil(max(typical_v(GPs(i).ind)))+10]); hold on; %regional histogram
    if max(region_normv)<max(typical_v(GPs(i).ind)) %make sure the fastest-flowing glacier is included in the sample
        [maxv,maxref] = max(typical_v(GPs(i).ind(v_ind)));
        if coverage_v(GPs(i).ind(maxref)) > 0.5*size(v,3)
            region_refs = [region_refs v_ind(maxref)];
            region_boxes = [region_boxes box_id(GPs(i).ind(v_ind(maxref)))];
            region_x = [region_x x_v(GPs(i).ind(v_ind(maxref)))]; region_y = [region_y y_v(GPs(i).ind(v_ind(maxref)))];
            region_fracv = [region_fracv fraction_v(GPs(i).ind(v_ind(maxref)))];
            region_normv = [region_normv typical_v(GPs(i).ind(v_ind(maxref)))];
            region_coverv = [region_coverv coverage_v(GPs(i).ind(v_ind(maxref)))];
        end
    end
    histogram(region_normv,[0:10:ceil(max(typical_v(GPs(i).ind)))+10]); hold on; %sample histogram
    
    %compile the subsetted data
    if i == 1; GPs(i).name = 'NN'; elseif i == 2; GPs(i).name = 'NW'; elseif i == 3; GPs(i).name = 'SW'; elseif i == 4; GPs(i).name = 'SE'; else GPs(i).name = 'NE'; end
    GPs(i).sorted_refs = region_refs;
    GPs(i).boxIDs = region_boxes;
    GPs(i).X = region_x; GPs(i).Y = region_y;
    GPs(i).fractionV = region_fracv;
    GPs(i).typicalV = region_normv;
    GPs(i).numberV = region_coverv;

    
    %plot locations on the map
    figure(mapfig); plot(GPs(i).X,GPs(i).Y,'xm'); hold on;
    
    %write a text file with box numbers to focus on
    writematrix(region_boxes',['GreenlandPeriph_',GPs(i).name,'-focus-areas.txt'],'delimiter','tab');
        
    %clear regional data
    clear region_* max* v_ind;
    disp(['...identified the glaciers w/ the best velocity coverage for the ',GPs(i).name]);
end
disp('Done! Check output text files containing box IDs for each region.');




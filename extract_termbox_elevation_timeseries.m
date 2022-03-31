%extract elevation time series from interior edge of terminus boxes
clearvars; close all;

% %navigate to the box edge shapefile
% cd /Users/ellynenderlin/Research/NASA_GreenlandPeriph-Mapping/terminus_boxes/box-edges/
% S = shaperead('BoxBacks_buffered.shp'); 
% S = shaperead('insideEdge_buffered.shp');

%navigate to the matlab structure with the flux gates
cd /Users/ellynenderlin/Research/NASA_GreenlandPeriph-Mapping/
load Greenland_GIC_centerlines.mat;
for i = 1:length(term)
    BoxID(i) = term(i).BoxID;
end

%navigate to and load the GIMP image mosaic for checking referencing of the
%various datasets
cd /Users/ellynenderlin/Research/miscellaneous/GIMP
load GIMP_image_mosaic_150m.mat;
bigfig = figure; set(bigfig,'position',[50 50 1200 1200]);
imagesc(I.x,I.y,I.z); colormap gray; axis xy equal; hold on;

%set up days of year for leap and non-leap years for date conversions
modays = [31 28 31 30 31 30 31 31 30 31 30 31]; cumdays = [0 cumsum(modays(1:end-1))];
leap_modays = [31 29 31 30 31 30 31 31 30 31 30 31]; leap_cumdays = [0 cumsum(leap_modays(1:end-1))];

%navigate to the DEM directory
cd /Users/ellynenderlin/Research/NASA_GreenlandPeriph-Mapping/example_DEMs/
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
            [A,R] = geotiffread(DEMs(j).name); A(A==-9999) = NaN;
            
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
            
%             %identify the shapefile that corresponds to the box
%             for k = 1:length(S)
%                 if nanmean(S(k).X)>=min(x) && nanmean(S(k).X)<=max(x) && nanmean(S(k).Y)>=min(y) && nanmean(S(k).Y)<=max(y)
%                     Sref = k; %disp(num2str(k));
%                 end
%             end
            
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
%                 plot(S(Sref).X,S(Sref).Y,'-c'); hold on;
                set(gca,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)]); drawnow;
            end
            
            %clear the variables
            clear datestring A R y x xgrid ygrid inland seaward;
        end
        
        %add the data to the term structure
%         term(termref).bufferX = S(Sref).X; term(termref).bufferY = S(Sref).Y;
        [sorted,sortrefs] = sort(elev_date);
        term(termref).inlandZmed = elev_medi(sortrefs); term(termref).inlandZmad = elev_madi(sortrefs); term(termref).inlandZyrs = sorted;
        term(termref).seawardZmed = elev_meds(sortrefs); term(termref).seawardZmad = elev_mads(sortrefs); term(termref).seawardZyrs = sorted;
        clear sorted sortrefs;
        
        %create elevation change time series and save
        elev_ts = figure; 
        plot(term(termref).inlandZyrs,term(termref).inlandZmed,'-','color',[44,123,182]/255,'linewidth',2); hold on;
        fill([term(termref).inlandZyrs fliplr(term(termref).inlandZyrs)],[term(termref).inlandZmed+term(termref).inlandZmad fliplr([term(termref).inlandZmed-term(termref).inlandZmad])],[44,123,182]/255,'facealpha',0.5,'edgecolor','none');
        plot(term(termref).seawardZyrs,term(termref).seawardZmed,'-','color',[253,174,97]/255,'linewidth',2); hold on;
        fill([term(termref).seawardZyrs fliplr(term(termref).seawardZyrs)],[term(termref).seawardZmed+term(termref).seawardZmad fliplr([term(termref).seawardZmed-term(termref).seawardZmad])],[253,174,97]/255,'facealpha',0.5,'edgecolor','none');
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
        
%         disp(['Elevation data extracted for ',files(i).name,' (buffer ref=',num2str(Sref),' & term ref=',num2str(termref),')']);
        disp(['Elevation data extracted for DEM set #',files(i).name,' (term ref=',num2str(termref),')']);
        clear termref elev_* DEMs;
    end
end
disp('Elevation time series extracted for ALL DEMs!!!');
cd /Users/ellynenderlin/Research/NASA_GreenlandPeriph-Mapping/centerlines/
save('Greenland_GIC_centerlines.mat','term','-v7.3');


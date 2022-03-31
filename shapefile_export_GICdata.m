%export Greenland GIC terminus boxes that were automatically adjusted when
%1985 terminus positions were added so that they are each saved as a unique
%shapefile that can be edited in QGIS/ArcGIS

clearvars; close all;
cd /Users/ellynenderlin/Research/NASA_GreenlandPeriph-Mapping
load Greenland_GIC_centerlines.mat
cd terminus_boxes/

%land-terminating glacier box IDs!
land_BoxID = [484 485 222 140 142 163 103 608];

%load an old shapefile as a reference
% cd old
% S = shaperead('1.shp');
% cd ..

%loop through, checking for the modification flag
for i = 1:length(term)
    s.Geometry = 'Polygon';
    s.BoundingBox = [min(term(i).BoxX) min(term(i).BoxY); max(term(i).BoxX) max(term(i).BoxY)];
    s.X = term(i).BoxX; s.Y = term(i).BoxY;
    s.BoxID = term(i).BoxID;
    if term(i).BoxFlag == 1
        cd updated
        if term(i).BoxID < 10
            shapewrite(s,['Greenland-GIC_00',num2str(term(i).BoxID),'.shp']);
            copyfile('/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj',['./','Greenland-GIC_00',num2str(term(i).BoxID),'.prj']);
        elseif term(i).BoxID >= 10 && term(i).BoxID < 100
            shapewrite(s,['Greenland-GIC_0',num2str(term(i).BoxID),'.shp']);
            copyfile('/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj',['./','Greenland-GIC_0',num2str(term(i).BoxID),'.prj']);
        else
            shapewrite(s,['Greenland-GIC_',num2str(term(i).BoxID),'.shp']);
            copyfile('/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj',['./','Greenland-GIC_',num2str(term(i).BoxID),'.prj']);
        end
        
        cd ..
    else
        if term(i).BoxID < 10
            shapewrite(s,['Greenland-GIC_00',num2str(term(i).BoxID),'.shp']);
            copyfile('/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj',['./','Greenland-GIC_00',num2str(term(i).BoxID),'.prj']);
        elseif term(i).BoxID >= 10 && term(i).BoxID < 100
            shapewrite(s,['Greenland-GIC_0',num2str(term(i).BoxID),'.shp']);
            copyfile('/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj',['./','Greenland-GIC_0',num2str(term(i).BoxID),'.prj']);
        else
            shapewrite(s,['Greenland-GIC_',num2str(term(i).BoxID),'.shp']);
            copyfile('/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj',['./','Greenland-GIC_',num2str(term(i).BoxID),'.prj']);
        end
    end
    clear s;
    
end
disp('All Greenland GIC terminus boxes written to shapefiles');


%% export centerlines to shapefiles
cd /Users/ellynenderlin/Research/NASA_GreenlandPeriph-Mapping/centerlines
for i = 1:length(term)
    if term(i).MankoffFlag ~= 1 && isempty(find(land_BoxID == term(i).BoxID)) %added this filtering after initial centerline export
        
        %centerlines
        s.Geometry = 'Polyline';
        s.BoundingBox = [min(term(i).centerX) min(term(i).centerY); max(term(i).centerX) max(term(i).centerY)];
        s.X = term(i).centerX; s.Y = term(i).centerY;
        s.BoxID = term(i).BoxID;
        if term(i).BoxID < 10
            shapewrite(s,['Greenland-GICcenterline_00',num2str(term(i).BoxID),'.shp']);
            copyfile('/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj',['./','Greenland-GICcenterline_00',num2str(term(i).BoxID),'.prj']);
        elseif term(i).BoxID >= 10 && term(i).BoxID < 100
            shapewrite(s,['Greenland-GICcenterline_0',num2str(term(i).BoxID),'.shp']);
            copyfile('/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj',['./','Greenland-GICcenterline_0',num2str(term(i).BoxID),'.prj']);
        else
            shapewrite(s,['Greenland-GICcenterline_',num2str(term(i).BoxID),'.shp']);
            copyfile('/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj',['./','Greenland-GICcenterline_',num2str(term(i).BoxID),'.prj']);
        end
        clear s;
        
        
        %flux gates
        cd ../fluxgates
        s.Geometry = 'Polyline';
        s.BoundingBox = [min(term(i).gateX) min(term(i).gateY); max(term(i).gateX) max(term(i).gateY)];
        s.X = term(i).gateX; s.Y = term(i).gateY;
        s.BoxID = term(i).BoxID;
        if term(i).BoxID < 10
            shapewrite(s,['Greenland-GICfluxgate_00',num2str(term(i).BoxID),'.shp']);
            copyfile('/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj',['./','Greenland-GICfluxgate_00',num2str(term(i).BoxID),'.prj']);
        elseif term(i).BoxID >= 10 && term(i).BoxID < 100
            shapewrite(s,['Greenland-GICfluxgate_0',num2str(term(i).BoxID),'.shp']);
            copyfile('/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj',['./','Greenland-GICfluxgate_0',num2str(term(i).BoxID),'.prj']);
        else
            shapewrite(s,['Greenland-GICfluxgate_',num2str(term(i).BoxID),'.shp']);
            copyfile('/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj',['./','Greenland-GICfluxgate_',num2str(term(i).BoxID),'.prj']);
        end
        clear s;

        cd ../centerlines
    end
end
disp('All Greenland GIC centerlines & flux gates written to shapefiles');

disp('Note that all centerlines were originally exported but only used flux gates');



%% export terminus changes to csv
cd ..
for i = 1:length(term)
    if term(i).MankoffFlag ~= 1 && isempty(find(land_BoxID == term(i).BoxID))
        %calculate terminus change between observations
        for j = 1:2 %three times = two terminus change calculations
            if sum(isnan(term(i).center_termref(j:j+1))) < 1
                %change in length
                deltaL(i,j) = (term(i).centerline(term(i).center_termref(j+1)) - term(i).centerline(term(i).center_termref(j))); %meters retreat
                deltaLrate(i,j) = deltaL(i,j)./((term(i).year(j+1)+term(i).doy(j+1)/365)-(term(i).year(j)+term(i).doy(j)/365)); %meters retreat/yr ON AVERAGE
            else
                deltaL(i,j) = NaN; deltaLrate(i,j) = NaN;
            end
            
        end
        %calculate total terminus change over the time period
        deltaL_total(i,1) = nansum(deltaL(i,:));
        %calculate the average rate of change for the full time period
        if ~isnan(term(i).center_termref(1)) && ~isnan(term(i).center_termref(3))
            deltaLrate_total(i,1) = deltaL_total(i,1)./((term(i).year(3)+term(i).doy(3)/365)-(term(i).year(1)+term(i).doy(1)/365)); %meters retreat/yr ON AVERAGE
        else
            deltaLrate_total(i,1) = NaN;
        end
        
        %average terminus position
        Lx(i,1) = nanmean(term(i).X); Ly(i,1) = nanmean(term(i).Y);
        
        %decimal years for observations
        Ldates(i,:) = term(i).year+term(i).doy/365;
    else
        deltaL(i,1) = NaN; deltaL(i,2) = NaN; deltaLrate(i,1) = NaN; deltaLrate(i,2) = NaN;
        deltaL_total(i) = NaN; deltaLrate_total(i) = NaN;
    end
end

%combine average X & Y position, terminus observation decimal dates, & the
%terminus change and rates of change info
deltaL_combined = [Lx Ly Ldates deltaL deltaLrate deltaL_total deltaLrate_total];
csvwrite('Greenland-GIC_termchange.csv',deltaL_combined);



%% write RGI outlines to one big shapefile
clear s S; S = struct([]);
for i = 1:length(term)
    if term(i).MankoffFlag ~= 1 && isempty(find(land_BoxID == term(i).BoxID))
        if isempty(S)
            S(1).Geometry = 'Polygon';
            S(1).BoundingBox = [min(term(i).RGIX) min(term(i).RGIY); max(term(i).RGIX) max(term(i).RGIY)];
            S(1).X = term(i).RGIX; S(1).Y = term(i).RGIY;
            S(1).BoxID = term(i).BoxID;
        else
            Sref = length(S);
            S(Sref+1).Geometry = 'Polygon';
            S(Sref+1).BoundingBox = [min(term(i).RGIX) min(term(i).RGIY); max(term(i).RGIX) max(term(i).RGIY)];
            S(Sref+1).X = term(i).RGIX; S(Sref+1).Y = term(i).RGIY;
            S(Sref+1).BoxID = term(i).BoxID;
            
        end
    end
end
disp('All Greenland GIC RGI outlines written to a shapefile');
shapewrite(S,['Greenland-GIC_RGIoutlines.shp']);
copyfile('/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj',['./','Greenland-GIC_RGIoutlines.prj']);


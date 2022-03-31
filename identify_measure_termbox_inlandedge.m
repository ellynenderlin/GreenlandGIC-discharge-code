clearvars; close all; warning off;

%make sure matlab knows where to look for supporting files
addpath('/users/ellynenderlin/mfiles','/users/ellynenderlin/mfiles/general');

%navigate to the appropriate parent directory (move files to your local computer to MASSIVELY speed up)
root_dir = '/Users/ellynenderlin/Research/NASA_Greenland-Periph-Mapping/';

%load the RGI shapefiles
cd_to_outlines = ['cd ',root_dir,'RGI-outlines/']; eval(cd_to_outlines);
S = shaperead('05_rgi50_GreenlandPeriphery_PS_BoxIDs_TW.shp');
for i = 1:length(S)
        RGIBoxID(i) = S(i).BoxID;
end

%load the centerline data structure
cd_to_centerlines = ['cd ',root_dir,'centerlines/']; eval(cd_to_centerlines);
load Greenland_GIC_centerlines.mat;
    
%load the manual terminus delineations
cd_to_termini = ['cd ',root_dir,'terminus_traces/']; eval(cd_to_termini);
T = shaperead('delineations_2000.shp');
for i = 1:length(T)
        BoxID(i) = T(i).BoxID;
end

%find where terminus traces intersect box walls
for j = 1:length(term)
    termref = find(BoxID==term(j).BoxID); 
    %term(j).RGIshaperef = find(RGIBoxID <= term(j).RGIref,1,'last');
%     term(j).RGIX = S(term(j).RGIshaperef).X; term(j).RGIY = S(term(j).RGIshaperef).Y; 
    
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

    
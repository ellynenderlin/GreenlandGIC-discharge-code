%%% Use this script to pair manually-delineated terminus boxes with RGI
%%% outlines & manual terminus traces, draw centerlines, create centerline
%%% profiles & centerline terminus time series, & draw fluxgates.
%%%  

disp('This code combines several old codes so there may be some hiccups running it through!');

%% initialize (RUN THIS EVERY TIME)
%loop through all the elevation data for each study site to locate the grounding line
clearvars; close all; warning off; beep on
addpath('/users/ellynenderlin/mfiles','/users/ellynenderlin/mfiles/general');

%specify the root directory for project-specific files
root_dir = '/Users/ellynenderlin/Research/NASA_GreenlandPeriph-Mapping/'; %including trailing /
%specify the root directory for generic files
misc_dir = '/Users/ellynenderlin/Research/miscellaneous/'; %including trailing /

%load the RGI shapefiles
cd_to_outlines = ['cd ',root_dir,'/discharge/outlines/RGI/']; eval(cd_to_outlines);
RGI = shaperead('05_rgi50_GreenlandPeriphery_PS_BoxIDs_TW.shp');
for i = 1:length(RGI)
    RGIBoxID(i) = RGI(i).BoxID;
end

%make a directory for the centerline files if one doesn't already exist
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
status = mkdir('centerlines'); status = mkdir('fluxgates'); status = mkdir('terminus_traces');
cd_to_centerlines = ['cd ',root_dir,'centerlines/']; 

%set-up day of year variables
modays = [31 28 31 30 31 30 31 31 30 31 30 31]; cumdays = [0 cumsum(modays(1:end-1))];
leap_modays = [31 29 31 30 31 30 31 31 30 31 30 31]; leap_cumdays = [0 cumsum(leap_modays(1:end-1))];

%% check that the correct RGI outlines were assigned to each terminus box

%load the terminus boxes
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
if isfile('Greenland_GIC_centerlines.mat') == 0 %if the file containing the structure with centerline info doesn't exist, create it
    cd_to_boxes = ['cd ',root_dir,'terminus_boxes/']; eval(cd_to_boxes);
    term_files = dir('*.shp');
    for i = 1:length(term_files)
        if isempty(strfind(term_files(i).name, '._')) %ignore the ghost files on the external drive
        termshape=shaperead(term_files(i).name);
        term(i).BoxX = termshape.X; term(i).BoxY = termshape.Y; term(i).BoxID = termshape.BoxID;
        clear termshape;
        
        %identify the RGI outline corresponding to the terminus box
        IDdiff = RGIBoxID-term(i).BoxID; IDdiff(IDdiff>0) = NaN; [~,maxind] = max(IDdiff); %find the RGI outline (could have the same BoxID or a lesser value depending if there are multiple outlets)
        term(i).RGIref = maxind; clear IDdiff maxind;
        term(i).RGIshaperef = find(RGIBoxID <= term(i).BoxID,1,'last');
        term(i).centerline = [];
        term(i).centerGL = [];
        end
    end
    eval(cd_to_root);
    save('Greenland_GIC_centerlines.mat','term','-v7.3')
else
    eval(cd_to_root);
    load Greenland_GIC_centerlines.mat;
end

%load the velocity mosaic
cd_to_vels = ['cd ',misc_dir,'Greenland-VelMosaic_1995-2015/']; eval(cd_to_vels);
[I,R] = readgeoraster('greenland_vel_mosaic250_vx_v1.tif');
V.x = single([R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX]);
V.y = single([R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY]);
V.vx = single(I); V.vx(V.vx==-2.0000e+09) = NaN; %velocity in x-direction in m/yr
clear I R;
[I,~] = readgeoraster('greenland_vel_mosaic250_vy_v1.tif');
V.vy = single(I); V.vy(V.vy==-2.0000e+09) = NaN; %velocity in y-direction in m/yr
clear I;

%check that the correct RGI polygon was identified
disp('checking that the correct RGI polygon was identified for each terminus box...');
fig = figure; set(fig,'position',[550 50 800 800]);
imagesc(V.x,V.y,sqrt(V.vx.^2 + V.vy.^2)); axis xy equal; colormap gray; set(gca,'clim',[0 500]); hold on;
badRGI = [];
for j = 1:length(term)
%     disp(['for ',num2str(j),' of ',num2str(length(term))]);
    poly1 = polyshape(RGI(term(j).RGIshaperef).X, RGI(term(j).RGIshaperef).Y);
    poly2 = polyshape(term(j).BoxX(~isnan(term(j).BoxX)),term(j).BoxY(~isnan(term(j).BoxY)));
    polyout = intersect(poly1,poly2);
    if polyout.NumRegions == 0; %there are no overlapping regions so the incorrect RGI outline was selected!
        clear poly1 poly2 polyout; disp('wrong RGI polygon!');
        term(j).RGIref = []; term(j).RGIshaperef = [];
%         tic; 
        k = 1;
        while k <= length(RGI)
            if min(sqrt((nanmean(term(j).BoxX(~isnan(term(j).BoxX)))-RGI(k).X).^2 + (nanmean(term(j).BoxY(~isnan(term(j).BoxY)))-RGI(k).Y).^2)) < 10000
                poly1 = polyshape(RGI(k).X, RGI(k).Y);
                poly2 = polyshape(term(j).BoxX(~isnan(term(j).BoxX)),term(j).BoxY(~isnan(term(j).BoxY)));
                polyout = intersect(poly1,poly2);
                if polyout.NumRegions > 0
                    disp(['overlap for BoxID',num2str(term(j).BoxID),' w/ RGI ID ',num2str(RGI(k).BoxID)]);
                    term(j).RGIref = RGI(k).BoxID;
                    term(j).RGIshaperef = k;
                    clear poly1 poly2 polyout;
                    break
                end
            end
            k=k+1;
        end
%         toc
    end
    
    %plot to double-check
    if ~isempty(term(j).RGIshaperef)
    figure(fig); 
    plot(term(j).BoxX(~isnan(term(j).BoxX)),term(j).BoxY(~isnan(term(j).BoxY)),'--c','linewidth',2); hold on;
%     set(gca,'xlim',[min(term(j).BoxX(~isnan(term(j).BoxX))) max(term(j).BoxX(~isnan(term(j).BoxX)))],'ylim',[min(term(j).BoxY(~isnan(term(j).BoxY))) max(term(j).BoxY(~isnan(term(j).BoxY)))]); drawnow;
    plot(RGI(term(j).RGIshaperef).X, RGI(term(j).RGIshaperef).Y,'-r','linewidth',2); hold on;
    set(gca,'xlim',[min(RGI(term(j).RGIshaperef).X) max(RGI(term(j).RGIshaperef).X)],'ylim',[min(RGI(term(j).RGIshaperef).Y) max(RGI(term(j).RGIshaperef).Y)]);
    drawnow;
    else
       disp(['NO OVERLAPPING RGI OUTLINE FOR ',num2str(term(j).BoxID)]); 
       badRGI = [badRGI; j];
    end
end
%now assign the nearest neighboring RGI boxes to those that don't have any
%overlap & modify the boxes as necessary
for j = 1:length(badRGI)
    for k = 1:length(RGI)
        RGI_Box_dist(k) = min(sqrt((nanmean(term(badRGI(j)).X(~isnan(term(badRGI(j)).X)))-RGI(k).X).^2 + (nanmean(term(badRGI(j)).Y(~isnan(term(badRGI(j)).Y)))-RGI(k).Y).^2));
    end
    [minval,minref] = min(RGI_Box_dist);
    term(badRGI(j)).RGIref = RGI(minref).BoxID;
    term(badRGI(j)).RGIshaperef = minref;
    clear minval minref;
    
    %plot
    figure(fig); 
    plot(term(badRGI(j)).X(~isnan(term(badRGI(j)).X)),term(badRGI(j)).Y(~isnan(term(badRGI(j)).Y)),'--c','linewidth',2); hold on;
    plot(RGI(term(badRGI(j)).RGIshaperef).X, RGI(term(badRGI(j)).RGIshaperef).Y,'-r','linewidth',2); hold on;
    set(gca,'xlim',[min([term(badRGI(j)).X(~isnan(term(badRGI(j)).X)) RGI(term(badRGI(j)).RGIshaperef).X]) max([term(badRGI(j)).X(~isnan(term(badRGI(j)).X)) RGI(term(badRGI(j)).RGIshaperef).X])],...
        'ylim',[min([term(badRGI(j)).Y(~isnan(term(badRGI(j)).Y)) RGI(term(badRGI(j)).RGIshaperef).Y]) max([term(badRGI(j)).Y(~isnan(term(badRGI(j)).Y)) RGI(term(badRGI(j)).RGIshaperef).Y])]);
    drawnow;
end
%add RGI polygons to the structure
for j = 1:length(term)
    term(j).RGIX = RGI(term(j).RGIshaperef).X; term(j).RGIY = RGI(term(j).RGIshaperef).Y; 
end
disp('Done assigning RGI outlines to boxes'); clear badRGI; close(fig);
eval(cd_to_root);
save('Greenland_GIC_centerlines.mat','term','-v7.3');
clearvars; close all;

%% add manual terminus delineations as desired
disp('Add manual terminus delineations... need to modify extensively for different dates');
disp('currently assumes data exist for 2000 & 2015 and added date (range) is 1985-1990');
add_datedir = '1985-1990'; add_datestamp = '1985';

%load the term structure
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
load Greenland_GIC_centerlines.mat;

%specify whether you are starting from scratch or restarting the addition
%of termini traces for a particular year
prompt = 'Starting a new year of terminus delineations (y/n)?';
str = input(prompt,'s');

%load the existing terminus traces as a reference
cd_to_termini = ['cd ',root_dir,'terminus_traces/']; eval(cd_to_termini);
S = shaperead('ManTermDelins_2000.shp'); 
for j = 1:length(S)
    BoxID(j) = S(j).BoxID; %create list of BoxIDs for 2000 termini
    if strmatch(str,'y')==1
        s(j).X = S(j).X; s(j).Y = S(j).Y; %create dummy structure for new termini traces
        s(j).Geometry = 'Polyline'; s(j).BoundingBox = double([min(S(j).X) max(S(j).X); min(S(j).Y) max(S(j).Y)]);
        s(j).BoxID = 0; s(j).RGIId = 0;
        s(j).year = 1985; s(j).DOY = 0;
        s(j).SourceID = 'image name';
    end
end
Snew = shaperead('ManTermDelins_2015.shp'); 
for j = 1:length(Snew)
    BoxIDnew(j) = Snew(j).BoxID; %create list of BoxIDs for 2015 termini
end

%if starting from scratch, create flag vectors flags to identify if traces from other dates and/or terminus box should be redone 
if strmatch(str,'y')==1
    flag2000 = zeros(size(S)); flag2015 = zeros(size(S)); flagbox = zeros(size(S)); %1=redo, 0=OK
else
    %load the flags and the traces for this date so far
    s = shaperead(['ManTermDelins_',add_datestamp,'.shp']); 
    load GreenlandGIC_bad-termini_boxes.mat;
    cd ../2000; load GreenlandGIC_bad-termini_2000.mat;
    cd ../2015; load GreenlandGIC_bad-termini_2015.mat;
end
clear str;
    
%locate the Landsat image folders
cd_to_datedir = ['cd ../',add_datedir]; eval(cd_to_datedir); %use if starting from scratch
% cd ../added %use if only filling in holes w/ extra Landsat images
L5s = dir('LT05*'); L4s = dir('LT04*'); L7s = dir('LE07*'); 

%loop through the Landsat 5 folders and delineate termini contained in each image
%(check if the terminus delineation has already been created each time an image is opened so efforts are not duplicated)
disp('Looping through Landsat 5 images...');
for i = 1:length(L5s)
    disp(['Landsat scene #',num2str(i),' of ',num2str(length(L5s)),': ',L5s(i).name]);
    cd_to_dir = ['cd ',L5s(i).name]; eval(cd_to_dir);
    bands = dir('*B*PS.tif');
    [I,R] = readgeoraster(bands(3).name);
    im.x = R.XWorldLimits(1):R.SampleSpacingInWorldX:R.XWorldLimits(2);
    im.y = R.YWorldLimits(2):-R.SampleSpacingInWorldY:R.YWorldLimits(1);
    im.z = double(I); clear I R;
    im_scaled = (im.z-min(min(im.z)))./(max(max(im.z))-min(min(im.z)));
    im.z_adjust = imadjust(im_scaled); clear im_scaled;
    figure; set(gcf,'position',[50 50 800 800])
    imagesc(im.x,im.y,im.z_adjust); colormap gray; axis xy equal; hold on;
    [imx,imy] = meshgrid(im.x,im.y);
    drawnow;
    
    %loop through the glaciers
    for j = 1:length(term)
        
        %find if the terminus falls within the image
        if nanmean(term(j).BoxX) > min(im.x) && nanmean(term(j).BoxX) < max(im.x) && nanmean(term(j).BoxY) > min(im.y) && nanmean(term(j).BoxY) < max(im.y)
            
            %identify image pixels in the terminus box
            in_termbox = inpolygon(imx,imy,term(j).BoxX(~isnan(term(j).BoxX)),term(j).BoxY(~isnan(term(j).BoxY)));
            im_cropped = in_termbox.*im.z;
            plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'-k','linewidth',2); hold on;
            
            %check that the terminus overlies an area with data AND that
            %you haven't already drawn the terminus (Landsat scenes
            %overlap, so you otherwise could be asked to draw the
            %it several times)
            if s(j).BoxID == 0 && sum(sum(im_cropped)) > 0
                disp(['terminus ',num2str(j),' for BoxID ',num2str(term(j).BoxID),' in image bounds & undelineated']);
                
                %only draw the terminus if ALL of the image pixels in the terminus box have data
                im_data = zeros(size(im.z)); im_data(in_termbox==1 & im.z~=0) = 1;
%                 if sum(sum(im_data)) == sum(sum(in_termbox)) 
                    disp('...use image to draw terminus');
                    
                    %zoom in & overlay 2000 and 2015 termini as references
                    plot(S(BoxID==term(j).BoxID).X,S(BoxID==term(j).BoxID).Y,'-m','linewidth',2); hold on;
                    plot(Snew(BoxIDnew==term(j).BoxID).X,Snew(BoxIDnew==term(j).BoxID).Y,'-c','linewidth',2); hold on;
                    set(gca,'xlim',[min(term(j).BoxX)-2000 max(term(j).BoxX)+2000],'ylim',[min(term(j).BoxY)-2000 max(term(j).BoxY)+2000]);
                    drawnow;
                    
                    %ask the user if the image should be used for delineation (no if cloud-covered)
                    prompt = 'Is the terminus cloud-covered (y/n)?';
                    str = input(prompt,'s');
                    if strmatch(str,'n')==1
                        %fill in flags
                        prompt = 'Does the 2000 (MAGENTA) terminus need to be retraced (y/n)?';
                        str = input(prompt,'s');
                        if strmatch(str,'y')==1
                            flag2000(j) = 1;
                        end
                        prompt = 'Does the 2015 (CYAN) terminus need to be retraced (y/n)?';
                        str = input(prompt,'s');
                        if strmatch(str,'y')==1
                            flag2015(j) = 1;
                        end
                        prompt = 'Does the terminus box need to be redone (y/n)?';
                        str = input(prompt,'s');
                        if strmatch(str,'y')==1
                            flagbox(j) = 1;
                        end
                        
                        %delineate terminus
                        [~,~,~,xi,yi] = improfile; %be patient, a + cursor will appear, double right click to terminate
                        plot(xi,yi,'-y','linewidth',2); hold on; drawnow;
                        s(j).X = double(xi); s(j).Y = double(yi); %add to shapefile structure
                        s(j).Geometry = 'Polyline';
                        s(j).BoundingBox = double([min(xi) max(xi); min(yi) max(yi)]);
                        s(j).BoxID = term(j).BoxID; s(j).RGIId = term(j).RGIref;
                        s(j).year = str2num(L5s(i).name(18:21));
                        if mod(str2num(L5s(i).name(18:21)),4) == 0
                            s(j).DOY = leap_cumdays(str2num(L5s(i).name(22:23)))+str2num(L5s(i).name(24:25));
                        else
                            s(j).DOY = cumdays(str2num(L5s(i).name(22:23)))+str2num(L5s(i).name(24:25));
                        end
                        s(j).SourceID = L5s(i).name;
                        clear xi yi;
                    end
                    
                    %go back to full image shown
                    set(gca,'xlim',[min(im.x) max(im.x)],'ylim',[min(im.y) max(im.y)]);
                    
%                 else
%                     disp('... moving on b/c not all pixels in terminus box have data');
%                     
%                 end
                clear term_bbox;
            end
            clear in_* im_cropped im_data;
        end
    end
    clear bands im imx imy; 
    
    %check that no geometries got changed (no explanation for why this happens)
    for j = 1:length(s)
        s(j).Geometry = 'Polyline';
    end
    
    %save the shapefile & bad terminus trace flags
    cd ../../terminus_traces
    shapewrite(s,['ManTermDelins_',add_datestamp,'.shp']);
    save('GreenlandGIC_bad-termini_boxes.mat','flagbox');
    cd ../2000
    save('GreenlandGIC_bad-termini_2000.mat','flag2000');
    cd ../2015
    save('GreenlandGIC_bad-termini_2015.mat','flag2015');
    cd_to_datedir = ['cd ../',add_datedir]; eval(cd_to_datedir); 
%     cd ../added
    disp('moving on to next image...'); close all; drawnow;
    disp(' ');
end

%supplement the Landsat 5 data with Landsat 4
disp('Looping through Landsat 4 images...');
for i = 1:length(L4s)
    disp(['Landsat scene #',num2str(i),' of ',num2str(length(L4s)),': ',L4s(i).name]);
    cd_to_dir = ['cd ',L4s(i).name]; eval(cd_to_dir);
    bands = dir('*B*PS.tif');
    [I,R] = readgeoraster(bands(3).name);
    im.x = R.XWorldLimits(1):R.SampleSpacingInWorldX:R.XWorldLimits(2);
    im.y = R.YWorldLimits(2):-R.SampleSpacingInWorldY:R.YWorldLimits(1);
    im.z = double(I); clear I R;
    im_scaled = (im.z-min(min(im.z)))./(max(max(im.z))-min(min(im.z)));
    im.z_adjust = imadjust(im_scaled); clear im_scaled;
    figure; set(gcf,'position',[50 50 800 800])
    imagesc(im.x,im.y,im.z_adjust); colormap gray; axis xy equal; hold on;
    [imx,imy] = meshgrid(im.x,im.y);
    drawnow;
    
    %loop through the glaciers
    for j = 1:length(term)
        
        %find if the terminus falls within the image
        if nanmean(term(j).BoxX) > min(im.x) && nanmean(term(j).BoxX) < max(im.x) && nanmean(term(j).BoxY) > min(im.y) && nanmean(term(j).BoxY) < max(im.y)
            
            %identify image pixels in the terminus box
            in_termbox = inpolygon(imx,imy,term(j).BoxX(~isnan(term(j).BoxX)),term(j).BoxY(~isnan(term(j).BoxY)));
            im_cropped = in_termbox.*im.z;
            plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'-k','linewidth',2); hold on;
            
            %check that the terminus overlies an area with data AND that
            %you haven't already drawn the terminus (Landsat scenes
            %overlap, so you otherwise could be asked to draw the
            %it several times)
            if sum(sum(im_cropped)) > 0 && s(j).BoxID == 0
                disp(['terminus ',num2str(j),' for BoxID ',num2str(term(j).BoxID),' in image bounds & undelineated']);
                
                %only draw the terminus if ALL of the image pixels in the terminus box have data
                im_data = zeros(size(im.z)); im_data(in_termbox==1 & im.z~=0) = 1;
%                 if sum(sum(im_data)) == sum(sum(in_termbox)) 
                    disp('...use image to draw terminus');
                    
                    %zoom in & overlay 2000 and 2015 termini as references
                    plot(S(BoxID==term(j).BoxID).X,S(BoxID==term(j).BoxID).Y,'-m','linewidth',2); hold on;
                    plot(Snew(BoxIDnew==term(j).BoxID).X,Snew(BoxIDnew==term(j).BoxID).Y,'-c','linewidth',2); hold on;
                    set(gca,'xlim',[min(term(j).BoxX)-2000 max(term(j).BoxX)+2000],'ylim',[min(term(j).BoxY)-2000 max(term(j).BoxY)+2000]);
                    drawnow;
                    
                    %ask the user if the image should be used for delineation (no if cloud-covered)
                    prompt = 'Is the terminus cloud-covered (y/n)?';
                    str = input(prompt,'s');
                    if strmatch(str,'n')==1
                        %fill in flags
                        prompt = 'Does the 2000 (MAGENTA) terminus need to be retraced (y/n)?';
                        str = input(prompt,'s');
                        if strmatch(str,'y')==1
                            flag2000(j) = 1;
                        else
                            flag2000(j) = 0;
                        end
                        prompt = 'Does the 2015 (CYAN) terminus need to be retraced (y/n)?';
                        str = input(prompt,'s');
                        if strmatch(str,'y')==1
                            flag2015(j) = 1;
                        else
                            flag2015(j) = 0;
                        end
                        prompt = 'Does the terminus box need to be redone (y/n)?';
                        str = input(prompt,'s');
                        if strmatch(str,'y')==1
                            flagbox(j) = 1;
                        end
                        
                        %delineate terminus
                        [~,~,~,xi,yi] = improfile; %be patient, a + cursor will appear, double right click to terminate
                        plot(xi,yi,'-y','linewidth',2); hold on; drawnow;
                        s(j).X = double(xi); s(j).Y = double(yi); %add to shapefile structure
                        s(j).Geometry = 'Polyline';
                        s(j).BoundingBox = double([min(xi) max(xi); min(yi) max(yi)]);
                        s(j).BoxID = term(j).BoxID; s(j).RGIId = term(j).RGIref;
                        s(j).year = str2num(L4s(i).name(18:21));
                        if mod(str2num(L4s(i).name(18:21)),4) == 0
                            s(j).DOY = leap_cumdays(str2num(L4s(i).name(22:23)))+str2num(L4s(i).name(24:25));
                        else
                            s(j).DOY = cumdays(str2num(L4s(i).name(22:23)))+str2num(L4s(i).name(24:25));
                        end
                        s(j).SourceID = L4s(i).name;
                        clear xi yi;
                    end
                    
                    %go back to full image shown
                    set(gca,'xlim',[min(im.x) max(im.x)],'ylim',[min(im.y) max(im.y)]);
                    
%                 else
%                     disp('... moving on b/c not all pixels in terminus box have data');
%                     
%                 end
                clear term_bbox;
            end
            clear in_* im_cropped im_data;
        end
    end
    clear bands im imx imy;  
    
    %check that no geometries got changed (no explanation for why this happens)
    for j = 1:length(s)
        s(j).Geometry = 'Polyline';
    end
    
    %save the shapefile & bad terminus trace flags
    cd ../../terminus_traces
    shapewrite(s,['ManTermDelins_',add_datestamp,'.shp']);
    save('GreenlandGIC_bad-termini_boxes.mat','flagbox');
    cd ../2000
    save('GreenlandGIC_bad-termini_2000.mat','flag2000');
    cd ../2015
    save('GreenlandGIC_bad-termini_2015.mat','flag2015');
    cd_to_datedir = ['cd ../',add_datedir]; eval(cd_to_datedir); 
%     cd ../added
    disp('moving on to next image...'); close all; drawnow;
    disp(' ');
end

%fill in as many in the farth north  as possible with data from 1999
disp('Looping through Landsat 7 images...');
for i = 1:length(L7s)
    disp(['Landsat scene #',num2str(i),' of ',num2str(length(L7s)),': ',L7s(i).name]);
    cd_to_dir = ['cd ',L7s(i).name]; eval(cd_to_dir);
    bands = dir('*B*PS.tif');
    [I,R] = readgeoraster(bands(9).name);
    im.x = R.XWorldLimits(1):R.SampleSpacingInWorldX:R.XWorldLimits(2);
    im.y = R.YWorldLimits(2):-R.SampleSpacingInWorldY:R.YWorldLimits(1);
    im.z = double(I); clear I R;
    im_scaled = (im.z-min(min(im.z)))./(max(max(im.z))-min(min(im.z)));
    im.z_adjust = imadjust(im_scaled); clear im_scaled;
    figure; set(gcf,'position',[50 50 800 800])
    imagesc(im.x,im.y,im.z_adjust); colormap gray; axis xy equal; hold on;
    [imx,imy] = meshgrid(im.x,im.y);
    drawnow;
    
    %loop through the glaciers
    for j = 1:length(term)
        
        %find if the terminus falls within the image
        if nanmean(term(j).BoxX) > min(im.x) && nanmean(term(j).BoxX) < max(im.x) && nanmean(term(j).BoxY) > min(im.y) && nanmean(term(j).BoxY) < max(im.y)
            
            %identify image pixels in the terminus box
            in_termbox = inpolygon(imx,imy,term(j).BoxX(~isnan(term(j).BoxX)),term(j).BoxY(~isnan(term(j).BoxY)));
            im_cropped = in_termbox.*im.z;
            plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'-k','linewidth',2); hold on;
            
            %check that the terminus overlies an area with data AND that
            %you haven't already drawn the terminus (Landsat scenes
            %overlap, so you otherwise could be asked to draw the
            %it several times)
            if sum(sum(im_cropped)) > 0 && s(j).BoxID == 0
                disp(['terminus ',num2str(j),' for BoxID ',num2str(term(j).BoxID),' in image bounds & undelineated']);
                
                %only draw the terminus if ALL of the image pixels in the terminus box have data
                im_data = zeros(size(im.z)); im_data(in_termbox==1 & im.z~=0) = 1;
%                 if sum(sum(im_data)) == sum(sum(in_termbox)) 
                    disp('...use image to draw terminus');
                    
                    %zoom in & overlay 2000 and 2015 termini as references
                    plot(S(BoxID==term(j).BoxID).X,S(BoxID==term(j).BoxID).Y,'-m','linewidth',2); hold on;
                    plot(Snew(BoxIDnew==term(j).BoxID).X,Snew(BoxIDnew==term(j).BoxID).Y,'-c','linewidth',2); hold on;
                    set(gca,'xlim',[min(term(j).BoxX)-2000 max(term(j).BoxX)+2000],'ylim',[min(term(j).BoxY)-2000 max(term(j).BoxY)+2000]);
                    drawnow;
                    
                    %ask the user if the image should be used for delineation (no if cloud-covered)
                    prompt = 'Is the terminus cloud-covered (y/n)?';
                    str = input(prompt,'s');
                    if strmatch(str,'n')==1
                        %fill in flags
                        prompt = 'Does the 2000 (MAGENTA) terminus need to be retraced (y/n)?';
                        str = input(prompt,'s');
                        if strmatch(str,'y')==1
                            flag2000(j) = 1;
                        else
                            flag2000(j) = 0;
                        end
                        prompt = 'Does the 2015 (CYAN) terminus need to be retraced (y/n)?';
                        str = input(prompt,'s');
                        if strmatch(str,'y')==1
                            flag2015(j) = 1;
                        else
                            flag2015(j) = 0;
                        end
                        prompt = 'Does the terminus box need to be redone (y/n)?';
                        str = input(prompt,'s');
                        if strmatch(str,'y')==1
                            flagbox(j) = 1;
                        end
                        
                        %delineate terminus
                        [~,~,~,xi,yi] = improfile; %be patient, a + cursor will appear, double right click to terminate
                        plot(xi,yi,'-y','linewidth',2); hold on; drawnow;
                        s(j).X = double(xi); s(j).Y = double(yi); %add to shapefile structure
                        s(j).Geometry = 'Polyline';
                        s(j).BoundingBox = double([min(xi) max(xi); min(yi) max(yi)]);
                        s(j).BoxID = term(j).BoxID; s(j).RGIId = term(j).RGIref;
                        s(j).year = str2num(L7s(i).name(18:21));
                        if mod(str2num(L7s(i).name(18:21)),4) == 0
                            s(j).DOY = leap_cumdays(str2num(L7s(i).name(22:23)))+str2num(L7s(i).name(24:25));
                        else
                            s(j).DOY = cumdays(str2num(L7s(i).name(22:23)))+str2num(L7s(i).name(24:25));
                        end
                        s(j).SourceID = L7s(i).name;
                        clear xi yi;
                    end
                    
                    %go back to full image shown
                    set(gca,'xlim',[min(im.x) max(im.x)],'ylim',[min(im.y) max(im.y)]);
                    
%                 else
%                     disp('... moving on b/c not all pixels in terminus box have data');
%                     
%                 end
                clear term_bbox;
            end
            clear in_* im_cropped im_data;
        end
    end
    clear bands im imx imy;  
    
    %check that no geometries got changed (no explanation for why this happens)
    for j = 1:length(s)
        s(j).Geometry = 'Polyline';
    end
    
    %save the shapefile & bad terminus trace flags
    cd ../../terminus_traces
    shapewrite(s,['ManTermDelins_',add_datestamp,'.shp']);
    save('GreenlandGIC_bad-termini_boxes.mat','flagbox');
    cd ../2000
    save('GreenlandGIC_bad-termini_2000.mat','flag2000');
    cd ../2015
    save('GreenlandGIC_bad-termini_2015.mat','flag2015');
    cd_to_datedir = ['cd ../',add_datedir]; eval(cd_to_datedir); 
%     cd ../added
    disp('moving on to next image...'); close all; drawnow;
    disp(' ');
end
disp('Done adding termini using existing imagery!');

%create a matrix of Landsat path, Landsat row, center x, center y 
%identify imagery needed for more delineations as the path-row combos for
%imagery with centers closest to the terminus boxes
cd_to_datedir = ['cd ../',add_datedir]; eval(cd_to_datedir); 
L5s = dir('LT05*'); 
Lxcenter = []; Lycenter = []; Lpath = []; Lrow = [];
for i = 1:length(L5s)
    cd_to_dir = ['cd ',L5s(i).name]; eval(cd_to_dir);
    bands = dir('*B*PS.tif');
    [~,R] = readgeoraster(bands(3).name);
    Lxcenter = [Lxcenter; nanmean([R.XWorldLimits(1) R.XWorldLimits(2)])];
    Lycenter = [Lycenter; nanmean([R.YWorldLimits(1) R.YWorldLimits(2)])];
    Lpath = [Lpath; str2num(L5s(i).name(11:13))];
    Lrow = [Lrow; str2num(L5s(i).name(14:16))];
    cd ..
end
L4s = dir('LT04*'); 
for i = 1:length(L4s)
    cd_to_dir = ['cd ',L4s(i).name]; eval(cd_to_dir);
    bands = dir('*B*PS.tif');
    [~,R] = readgeoraster(bands(3).name);
    Lxcenter = [Lxcenter; nanmean([R.XWorldLimits(1) R.XWorldLimits(2)])];
    Lycenter = [Lycenter; nanmean([R.YWorldLimits(1) R.YWorldLimits(2)])];
    Lpath = [Lpath; str2num(L4s(i).name(11:13))];
    Lrow = [Lrow; str2num(L4s(i).name(14:16))];
    cd ..
end
L7s = dir('LE07*'); 
for i = 1:length(L7s)
    cd_to_dir = ['cd ',L7s(i).name]; eval(cd_to_dir);
    bands = dir('*B*PS.tif');
    [~,R] = readgeoraster(bands(9).name);
    Lxcenter = [Lxcenter; nanmean([R.XWorldLimits(1) R.XWorldLimits(2)])];
    Lycenter = [Lycenter; nanmean([R.YWorldLimits(1) R.YWorldLimits(2)])];
    Lpath = [Lpath; str2num(L7s(i).name(11:13))];
    Lrow = [Lrow; str2num(L7s(i).name(14:16))];
    cd ..
end


%determine images needed to finish added terminus positions
disp('Need to get more imagery to delineate the new time period''s termini...');
sref = []; Lref = []; Lpr = [];
for j = 1:length(s)
    if s(j).BoxID == 0
        imdist = sqrt((nanmean(s(j).X)-Lxcenter).^2 + (nanmean(s(j).Y)-Lycenter).^2);
        sref = [sref; j];
        Lref = [Lref; find(imdist == min(imdist))];
        Lpr = [Lpr; Lpath(find(imdist == min(imdist))) Lrow(find(imdist == min(imdist)))];
        disp([num2str(j)]);
    end
end
disp('... download the following path-rows:');
disp(unique(Lpr,'rows'));


%% use the flags & the path-row saved for each flagged glacier to download images needed to fix bad termini
add_datedir = '1985-1990'; add_datestamp = '1985';

%load the data
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
load Greenland_GIC_centerlines.mat;
cd terminus_traces
load GreenlandGIC_bad-termini_boxes.mat;

%create a matrix of Landsat path, Landsat row, center x, center y 
%identify imagery needed for more delineations as the path-row combos for
%imagery with centers closest to the terminus boxes
cd_to_datedir = ['cd ../',add_datedir]; eval(cd_to_datedir); 
L5s = dir('LT05*'); 
Lxcenter = []; Lycenter = []; Lpath = []; Lrow = [];
for i = 1:length(L5s)
    cd_to_dir = ['cd ',L5s(i).name]; eval(cd_to_dir);
    bands = dir('*B*PS.tif');
    [~,R] = readgeoraster(bands(3).name);
    Lxcenter = [Lxcenter; nanmean([R.XWorldLimits(1) R.XWorldLimits(2)])];
    Lycenter = [Lycenter; nanmean([R.YWorldLimits(1) R.YWorldLimits(2)])];
    Lpath = [Lpath; str2num(L5s(i).name(11:13))];
    Lrow = [Lrow; str2num(L5s(i).name(14:16))];
    cd ..
end
L4s = dir('LT04*'); 
for i = 1:length(L4s)
    cd_to_dir = ['cd ',L4s(i).name]; eval(cd_to_dir);
    bands = dir('*B*PS.tif');
    [~,R] = readgeoraster(bands(3).name);
    Lxcenter = [Lxcenter; nanmean([R.XWorldLimits(1) R.XWorldLimits(2)])];
    Lycenter = [Lycenter; nanmean([R.YWorldLimits(1) R.YWorldLimits(2)])];
    Lpath = [Lpath; str2num(L4s(i).name(11:13))];
    Lrow = [Lrow; str2num(L4s(i).name(14:16))];
    cd ..
end
L7s = dir('LE07*'); 
for i = 1:length(L7s)
    cd_to_dir = ['cd ',L7s(i).name]; eval(cd_to_dir);
    bands = dir('*B*PS.tif');
    [~,R] = readgeoraster(bands(9).name);
    Lxcenter = [Lxcenter; nanmean([R.XWorldLimits(1) R.XWorldLimits(2)])];
    Lycenter = [Lycenter; nanmean([R.YWorldLimits(1) R.YWorldLimits(2)])];
    Lpath = [Lpath; str2num(L7s(i).name(11:13))];
    Lrow = [Lrow; str2num(L7s(i).name(14:16))];
    cd ..
end

%identify imagery needed to fix bad delineations
disp('Identify Landsat images needed for years with prior delineations...');
cd ../2000
load GreenlandGIC_bad-termini_2000.mat;
disp('Need 2000 Landsat 7 images for path-row...');
for j = 1:length(flag2000)
    if flag2000(j) == 1
        disp([s(j).SourceID(11:13),'-',s(j).SourceID(14:16)]);
    end
end
cd ../2015
load GreenlandGIC_bad-termini_2015.mat;
disp('Need 2015 Landsat 8 images for path-row...');
for j = 1:length(flag2015)
    if flag2015(j) == 1
        disp([s(j).SourceID(11:13),'-',s(j).SourceID(14:16)]);
    end
end

disp('Run next section to remap bad termini from other dates after downloading needed imagery');

%% Fix terminus traces identified as potentially bad when the new traces were added
add_datedir = '1985-1990'; add_datestamp = '1985';

%load the data in the term structure
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
load Greenland_GIC_centerlines.mat;

%load the existing terminus traces
cd_to_termini = ['cd ',root_dir,'terminus_traces/']; eval(cd_to_termini);
S = shaperead('ManTermDelins_2000.shp'); 
for j = 1:length(S)
    BoxID(j) = S(j).BoxID; %create list of BoxIDs for 2000 termini
end
Snew = shaperead('ManTermDelins_2015.shp'); 
for j = 1:length(Snew)
    BoxIDnew(j) = Snew(j).BoxID; %create list of BoxIDs for 2015 termini
end
s = shaperead(['ManTermDelins_',add_datestamp,'.shp']); 
for j = 1:length(s)
    BoxIDold(j) = s(j).BoxID; %create list of BoxIDs for 2000 termini
end
load GreenlandGIC_bad-termini_boxes.mat;
cd ../2000; load GreenlandGIC_bad-termini_2000.mat;
cd ../2015; load GreenlandGIC_bad-termini_2015.mat;

% %check that terminus traces are assigned the correct BoxID
% disp('Checking that termini are assigned to the correct BoxIDs');
% for j = 1:length(term)
%     %2000
%     termref = find(BoxID==term(j).BoxID);
%     in = inpolygon(S(termref).X,S(termref).Y,term(j).BoxX,term(j).BoxY);
%     if isempty(in)
%         S(termref).BoxID = []; flag2000(j) = 1;
%         disp(['wrong 2000 terminus trace for ',num2str(j)]);
%         %find the correct terminus position
%         for k = 1:length(S)
%             in = inpolygon(S(k).X,S(k).Y,term(j).BoxX,term(j).BoxY);
%             if ~isempty(in)
%                S(k).BoxID = term(j).BoxID; BoxID(k) = term(j).BoxID;
%             end
%             clear in;  
%         end
%     end
%     clear termref;
%     
%     %2015
%     termref = find(BoxIDnew==term(j).BoxID);
%     in = inpolygon(Snew(termref).X,Snew(termref).Y,term(j).BoxX,term(j).BoxY);
%     if isempty(in)
%         Snew(termref).BoxID = []; flag2015(j) = 1;
%         disp(['wrong 2015 terminus trace for ',num2str(j)]);
%         %find the correct terminus position
%         for k = 1:length(Snew)
%             in = inpolygon(Snew(k).X,Snew(k).Y,term(j).BoxX,term(j).BoxY);
%             if ~isempty(in)
%                Snew(k).BoxID = term(j).BoxID; BoxIDnew(k) = term(j).BoxID;
%             end
%             clear in;  
%         end
%     end
%     clear termref;
%     
% end
% cd ../2000; save('GreenlandGIC_bad-termini_2000.mat','flag2000');
% cd ../2015; save('GreenlandGIC_bad-termini_2015.mat','flag2015');
% %now find the appropriate term reference for mislabeled traces
% %2000
% for j = 1:length(S)
%     if isempty(S(j).BoxID)
%         for k = 1:length(term)
%             in = inpolygon(S(j).X,S(j).Y,term(k).BoxX,term(k).BoxY);
%             if ~isempty(in)
%                S(j).BoxID = term(k).BoxID; BoxID(j) = term(k).BoxID; 
%             end
%             clear in;
%         end
%     end
% end
% %2015
% for j = 1:length(Snew)
%     if isempty(Snew(j).BoxID)
%         for k = 1:length(term)
%             in = inpolygon(Snew(j).X,Snew(j).Y,term(k).BoxX,term(k).BoxY);
%             if ~isempty(in)
%                Snew(j).BoxID = term(k).BoxID; BoxIDnew(j) = term(k).BoxID; 
%             end
%             clear in;
%         end
%     end
% end
% cd ../terminus_traces
% shapewrite(S,'ManTermDelins_2000.shp');
% shapewrite(Snew,'ManTermDelins_2015.shp');


%loop through the Landsat 7 images from 2000
cd ../2000; 
% cd ../added;
L7s = dir('LE07*'); 
disp('Need to fix 2000 traces for:');
disp(num2str(find(flag2000==1)));
disp('Looping through Landsat 7 images from 2000...');
for i = 1:length(L7s)
    disp(['Landsat scene #',num2str(i),' of ',num2str(length(L7s)),': ',L7s(i).name]);
    cd_to_dir = ['cd ',L7s(i).name]; eval(cd_to_dir);
    bands = dir('*B*PS.tif');
    [I,R] = readgeoraster(bands(9).name);
    im.x = R.XWorldLimits(1):R.SampleSpacingInWorldX:R.XWorldLimits(2);
    im.y = R.YWorldLimits(2):-R.SampleSpacingInWorldY:R.YWorldLimits(1);
    im.z = double(I); clear I R;
    im_scaled = (im.z-min(min(im.z)))./(max(max(im.z))-min(min(im.z)));
    im.z_adjust = imadjust(im_scaled); clear im_scaled;
    figure; set(gcf,'position',[50 50 800 800])
    imagesc(im.x,im.y,im.z_adjust); colormap gray; axis xy equal; hold on;
    [imx,imy] = meshgrid(im.x,im.y);
    drawnow;
    
    %loop through the glaciers
    for j = 1:length(term)
        
        %find if the terminus falls within the image
        if nanmean(term(j).BoxX) > min(im.x) && nanmean(term(j).BoxX) < max(im.x) && nanmean(term(j).BoxY) > min(im.y) && nanmean(term(j).BoxY) < max(im.y)
            
            %identify image pixels in the terminus box
            in_termbox = inpolygon(imx,imy,term(j).BoxX(~isnan(term(j).BoxX)),term(j).BoxY(~isnan(term(j).BoxY)));
            im_cropped = in_termbox.*im.z;
            plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'-k','linewidth',2); hold on;
            
            %check that the terminus overlies an area with data AND that
            %you need to update the terminus trace
            termref = find(BoxID == term(j).BoxID);
            if flag2000(j) == 1 && sum(sum(im_cropped)) > 0
                disp(['terminus ',num2str(j),' for BoxID ',num2str(term(j).BoxID),' in image bounds & undelineated']);
                
                %only draw the terminus if ALL of the image pixels in the terminus box have data
                im_data = zeros(size(im.z)); im_data(in_termbox==1 & im.z~=0) = 1;
%                 if sum(sum(im_data)) == sum(sum(in_termbox)) 
                    disp('...use image to draw terminus');
                    
                    %zoom in & overlay pre-2000s and 2015 termini as references
                    plot(term(j).centerX,term(j).centerY,'--g','linewidth',1); hold on;
                    plot(s(BoxIDold==term(j).BoxID).X,s(BoxIDold==term(j).BoxID).Y,'-m','linewidth',2); hold on;
                    plot(Snew(BoxIDnew==term(j).BoxID).X,Snew(BoxIDnew==term(j).BoxID).Y,'-c','linewidth',2); hold on;
                    plot(S(termref).X,S(termref).Y,'--y','linewidth',1); hold on;
                    set(gca,'xlim',[min(term(j).BoxX)-2000 max(term(j).BoxX)+2000],'ylim',[min(term(j).BoxY)-2000 max(term(j).BoxY)+2000]);
                    drawnow;
                    
                    %ask the user if the image should be used for delineation (no if cloud-covered)
                    prompt = 'Is the terminus cloud-covered (y/n)?';
                    str = input(prompt,'s');
                    if strmatch(str,'n')==1
                        %fill in flags
                        prompt = 'Does the terminus box need to be redone (y/n)?';
                        str = input(prompt,'s');
                        if strmatch(str,'y')==1
                            flagbox(j) = 1;
                        end
                        
                        %delineate terminus
                        [~,~,~,xi,yi] = improfile; %be patient, a + cursor will appear, double right click to terminate
                        plot(xi,yi,'-y','linewidth',2); hold on; drawnow;
                        S(termref).X = double(xi); S(termref).Y = double(yi); %add to shapefile structure
                        S(termref).Geometry = 'Polyline';
                        S(termref).BoundingBox = double([min(xi) max(xi); min(yi) max(yi)]);
                        S(termref).BoxID = term(j).BoxID; 
                        S(termref).year = str2num(L7s(i).name(18:21));
                        if mod(str2num(L7s(i).name(18:21)),4) == 0
                            S(termref).DOY = leap_cumdays(str2num(L7s(i).name(22:23)))+str2num(L7s(i).name(24:25));
                        else
                            S(termref).DOY = cumdays(str2num(L7s(i).name(22:23)))+str2num(L7s(i).name(24:25));
                        end
                        S(termref).SourceID = L7s(i).name;
                        clear xi yi;
                        flag2000(j) = 0;
                    end
                    
                    %go back to full image shown
                    set(gca,'xlim',[min(im.x) max(im.x)],'ylim',[min(im.y) max(im.y)]);
                    
%                 else
%                     disp('... moving on b/c not all pixels in terminus box have data');
%                     
%                 end
                clear term_bbox;
            end
            clear in_* im_cropped im_data;
        end
    end
    clear bands im imx imy; 
    
    %check that no geometries mysteriously switched from a polyline
    for j = 1:length(S)
        S(j).Geometry = 'Polyline';
    end
    
    %save the shapefile & bad terminus trace flags
    cd ../../terminus_traces
    save('GreenlandGIC_bad-termini_boxes.mat','flagbox');
    shapewrite(S,'ManTermDelins_2000.shp');
    cd ../2000
    save('GreenlandGIC_bad-termini_2000.mat','flag2000');
%     cd ../added
    disp('moving on to next image...'); close all; drawnow;
    disp(' ');
end
disp('still need to redo traces in 2000 version of Landsat scenes (old names):');
badref = find(flag2000==1);
if ~isempty(badref)
    for j = 1:length(badref)
        disp(Snew(find(BoxIDnew == BoxIDold(badref(j)))).SourceID);
    end
else
    disp('none!')
end
clear badref;

%loop through the Landsat 8 images from 2015
cd ../2015; 
% cd ../added;
L8s = dir('LC08*'); 
disp('Need to fix 2015 traces for:');
disp(num2str(find(flag2015==1)));
disp('Looping through Landsat 8 images from 2015...');
for i = 1:length(L8s)
    disp(['Landsat scene #',num2str(i),' of ',num2str(length(L8s)),': ',L8s(i).name]);
    cd_to_dir = ['cd ',L8s(i).name]; eval(cd_to_dir);
    bands = dir('*B*PS.tif');
    [I,R] = readgeoraster(bands(10).name);
    im.x = R.XWorldLimits(1):R.SampleSpacingInWorldX:R.XWorldLimits(2);
    im.y = R.YWorldLimits(2):-R.SampleSpacingInWorldY:R.YWorldLimits(1);
    im.z = double(I); clear I R;
    im_scaled = (im.z-min(min(im.z)))./(max(max(im.z))-min(min(im.z)));
    im.z_adjust = imadjust(im_scaled); clear im_scaled;
    figure; set(gcf,'position',[50 50 800 800])
    imagesc(im.x,im.y,im.z_adjust); colormap gray; axis xy equal; hold on;
    [imx,imy] = meshgrid(im.x,im.y);
    drawnow;
    
    %loop through the glaciers
    for j = 1:length(term)
        
        %find if the terminus falls within the image
        if nanmean(term(j).BoxX) > min(im.x) && nanmean(term(j).BoxX) < max(im.x) && nanmean(term(j).BoxY) > min(im.y) && nanmean(term(j).BoxY) < max(im.y)
            
            %identify image pixels in the terminus box
            in_termbox = inpolygon(imx,imy,term(j).BoxX(~isnan(term(j).BoxX)),term(j).BoxY(~isnan(term(j).BoxY)));
            im_cropped = in_termbox.*im.z;
            plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'-k','linewidth',2); hold on;
            
            %check that the terminus overlies an area with data AND that
            %you need to update the terminus trace
            termref = find(BoxIDnew == term(j).BoxID);
            if flag2015(j) == 1 && sum(sum(im_cropped)) > 0
                disp(['terminus ',num2str(j),' for BoxID ',num2str(term(j).BoxID),' in image bounds & undelineated']);
                
                %only draw the terminus if ALL of the image pixels in the terminus box have data
                im_data = zeros(size(im.z)); im_data(in_termbox==1 & im.z~=0) = 1;
%                 if sum(sum(im_data)) == sum(sum(in_termbox)) 
                    disp('...use image to draw terminus');
                    
                    %zoom in & overlay pre-2000s and 2000 termini as references
                    plot(term(j).centerX,term(j).centerY,'--g','linewidth',1); hold on;
                    plot(s(BoxIDold==term(j).BoxID).X,s(BoxIDold==term(j).BoxID).Y,'-m','linewidth',2); hold on;
                    plot(S(BoxID==term(j).BoxID).X,S(BoxID==term(j).BoxID).Y,'-y','linewidth',2); hold on;
                    plot(Snew(termref).X,Snew(termref).Y,'--c','linewidth',1); hold on;
                    set(gca,'xlim',[min(term(j).BoxX)-2000 max(term(j).BoxX)+2000],'ylim',[min(term(j).BoxY)-2000 max(term(j).BoxY)+2000]);
                    drawnow;
                    
                    %ask the user if the image should be used for delineation (no if cloud-covered)
                    prompt = 'Is the terminus cloud-covered (y/n)?';
                    str = input(prompt,'s');
                    if strmatch(str,'n')==1
                        %fill in flags
                        prompt = 'Does the terminus box need to be redone (y/n)?';
                        str = input(prompt,'s');
                        if strmatch(str,'y')==1
                            flagbox(j) = 1;
                        end
                        
                        %delineate terminus
                        [~,~,~,xi,yi] = improfile; %be patient, a + cursor will appear, double right click to terminate
                        plot(xi,yi,'-c','linewidth',2); hold on; drawnow;
                        Snew(termref).X = double(xi); Snew(termref).Y = double(yi); %add to shapefile structure
                        Snew(termref).Geometry = 'Polyline';
                        Snew(termref).BoundingBox = double([min(xi) max(xi); min(yi) max(yi)]);
                        Snew(termref).BoxID = term(j).BoxID; 
                        Snew(termref).year = str2num(bands(10).name(18:21));
                        if mod(str2num(bands(10).name(18:21)),4) == 0
                            Snew(termref).DOY = leap_cumdays(str2num(bands(10).name(22:23)))+str2num(bands(10).name(24:25));
                        else
                            Snew(termref).DOY = cumdays(str2num(bands(10).name(22:23)))+str2num(bands(10).name(24:25));
                        end
                        Snew(termref).SourceID = bands(10).name(1:end-9);
                        clear xi yi;
                        flag2015(j) = 0;
                    end
                    
                    %go back to full image shown
                    set(gca,'xlim',[min(im.x) max(im.x)],'ylim',[min(im.y) max(im.y)]);
                    
%                 else
%                     disp('... moving on b/c not all pixels in terminus box have data');
%                     
%                 end
                clear term_bbox;
            end
            clear in_* im_cropped im_data;
        end
    end
    clear bands im imx imy; 
    
    %check that no geometries mysteriously switched from polylines
    for j = 1:length(Snew)
        Snew(j).Geometry = 'Polyline';
    end
    
    %save the shapefile & bad terminus trace flags
    cd ../../terminus_traces
    save('GreenlandGIC_bad-termini_boxes.mat','flagbox');
    shapewrite(Snew,'ManTermDelins_2015.shp');
    cd ../2015
    save('GreenlandGIC_bad-termini_2015.mat','flag2015');
%     cd ../added;
    disp('moving on to next image...'); close all; drawnow;
    disp(' ');
end
disp('still need to redo traces in Landsat scenes (old names):');
badref = find(flag2015==1);
if ~isempty(badref)
    for j = 1:length(badref)
        disp(Snew(find(BoxIDnew == BoxIDold(badref(j)))).SourceID);
    end
else
    disp('none!')
end
clear badref;

%save the updated bad box flag
cd ../terminus_traces
save('GreenlandGIC_bad-termini_boxes.mat','flagbox');


%% Draw the centerlines
% clearvars; close all;
add_datedir = '1985-1990'; add_datestamp = '1985';

%load the data in the term structure
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
load Greenland_GIC_centerlines.mat;

%load the velocity mosaic
cd_to_vels = ['cd ',misc_dir,'Greenland-VelMosaic_1995-2015/']; eval(cd_to_vels);
[I,R] = readgeoraster('greenland_vel_mosaic250_vx_v1.tif');
V.x = single([R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX]);
V.y = single([R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY]);
V.vx = single(I); V.vx(V.vx==-2.0000e+09) = NaN; %velocity in x-direction in m/yr
clear I R;
[I,~] = readgeoraster('greenland_vel_mosaic250_vy_v1.tif');
V.vy = single(I); V.vy(V.vy==-2.0000e+09) = NaN; %velocity in y-direction in m/yr
clear I;

%specify whether you are starting from scratch or restarting the addition
%of termini traces for a particular year
prompt = 'Starting a new year of terminus delineations (y/n)?';
str = input(prompt,'s');
%if running for the first time, set up a flag used to determine if the
%centerline has been drawn & the terminus box checked
if strmatch(str,'y')==1
    checked = zeros(length(term));
end

%load the terminus traces
cd terminus_traces
S = shaperead('ManTermDelins_2000.shp'); 
for j = 1:length(S)
    BoxID(j) = S(j).BoxID; %create list of BoxIDs for 2000 termini
end
Snew = shaperead('ManTermDelins_2015.shp'); 
for j = 1:length(Snew)
    BoxIDnew(j) = Snew(j).BoxID; %create list of BoxIDs for 2015 termini
end
s = shaperead(['ManTermDelins_',add_datestamp,'.shp']); 
for j = 1:length(s)
    BoxIDold(j) = s(j).BoxID; %create list of BoxIDs for 2000 termini
end
% load GreenlandGIC_bad-termini_boxes.mat;

%use the time period with the most Landsat images
cd_to_datedir = ['cd ../',add_datedir]; eval(cd_to_datedir);
L5s = dir('LT05*'); L4s = dir('LT04*'); L7s = dir('LE07*'); 

%loop through the Landsat 5 images
for i = 1:length(L5s) %DONE!!!
    disp(['Landsat scene #',num2str(i),' of ',num2str(length(L5s)),': ',L5s(i).name]);
    cd_to_dir = ['cd ',L5s(i).name]; eval(cd_to_dir);
    bands = dir('*B*PS.tif');
    [I,R] = readgeoraster(bands(3).name);
    im.x = R.XWorldLimits(1):R.SampleSpacingInWorldX:R.XWorldLimits(2);
    im.y = R.YWorldLimits(2):-R.SampleSpacingInWorldY:R.YWorldLimits(1);
    im.z = double(I); clear I R;
    im_scaled = (im.z-min(min(im.z)))./(max(max(im.z))-min(min(im.z)));
    im.z_adjust = imadjust(im_scaled,[0 0.9],[0 1],0.7); clear im_scaled;
    oldfig = figure; set(gcf,'position',[50 50 800 800])
    imagesc(im.x,im.y,im.z_adjust); colormap gray; axis xy equal; hold on;
    newfig = figure; set(gcf,'position',[850 50 800 800])
    imagesc(im.x,im.y,im.z_adjust); colormap gray; axis xy equal; hold on;
    [imx,imy] = meshgrid(im.x,im.y);
    drawnow;
    
    %crop velocities to the RGI outline & overlay as a quiver plot on image
    disp('Prepping velocity vectors to aid mapping of the flux gate at the grounding line, perpendicular to flow...');
    vx_sub = V.vx(find(V.y>=min(im.y) & V.y<=max(im.y)),find(V.x>=min(im.x) & V.x<=max(im.x))); %crop vx to the landsat image extent
    vy_sub = V.vy(find(V.y>=min(im.y) & V.y<=max(im.y)),find(V.x>=min(im.x) & V.x<=max(im.x))); %crop vy to the landsat image extent
    x_sub = V.x(1,find(V.x>=min(im.x) & V.x<=max(im.x)));
    y_sub = V.y(1,find(V.y>=min(im.y) & V.y<=max(im.y)));
    [Xsubgrid,Ysubgrid] = meshgrid(x_sub,y_sub);
    in = inpolygon(Xsubgrid,Ysubgrid,[term(j).RGIX(1:end-1) term(j).RGIX(1)],[term(j).RGIY(1:end-1) term(j).RGIY(1)]); %mask velocities outside the glacier so velocity arrows scale approrpriately
    maskedvx = vx_sub; maskedvx(in==0) = NaN; %mask-out other regions by turning anything not in=1 into NaNs
    maskedvy = vy_sub; maskedvy(in==0) = NaN; %mask-out other regions by turning anything not in=1 into NaNs
    clear in;
    figure(oldfig);
    q = quiver(Xsubgrid,Ysubgrid,maskedvx,maskedvy,'-k'); hold on; %q.AutoScaleFactor = 50;
    
    %loop through the glaciers
    for j = 1:length(term)
        %check that the terminus box is within the image
        if nanmean(term(j).BoxX) > min(im.x) && nanmean(term(j).BoxX) < max(im.x) && nanmean(term(j).BoxY) > min(im.y) && nanmean(term(j).BoxY) < max(im.y)
            %identify image pixels in the terminus box
            in_termbox = inpolygon(imx,imy,term(j).BoxX(~isnan(term(j).BoxX)),term(j).BoxY(~isnan(term(j).BoxY)));
            im_cropped = in_termbox.*im.z;
            %             %identify pixels in the RGI outline
            %             in_RGI = inpolygon(imx,imy,term(j).RGIX(~isnan(term(j).RGIX)),term(j).RGIY(~isnan(term(j).RGIY)));
            %             im_cropped = in_RGI.*im.z;
            
            %check that the terminus overlies an area with data
            if sum(sum(im_cropped)) > 0
                %if the RGI outline is ENORMOUS, skip identifying pixels in the RGI polygon (takes too long)
                RGI_bbox = (max(term(j).RGIX)-min(term(j).RGIX))*(max(term(j).RGIY)-min(RGI(term(j).RGIY)));
                if RGI_bbox <= 1e+09
                    in_RGI = inpolygon(imx,imy,[term(j).RGIX(1:end-1) term(j).RGIX(1)],[term(j).RGIY(1:end-1) term(j).RGIY(1)]);
                    im_data = zeros(size(im.z_adjust)); im_data(in_RGI==1 & im.z_adjust~=0) = 1;
                end
                
                %add the centerline and fix the terminus box as necessary
                if (RGI_bbox > 1e+09 || sum(sum(im_data)) == sum(sum(in_RGI)))  && checked(j) == 0
                    figure(oldfig);
                    plot(term(j).RGIX, term(j).RGIY, '-r','linewidth',2); hold on;
                    plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'-k','linewidth',2); hold on; plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'--w','linewidth',2); hold on;
                    if ~isempty(find(BoxIDold==term(j).BoxID))
                        plot(s(BoxIDold==term(j).BoxID).X,s(BoxIDold==term(j).BoxID).Y,'-m','linewidth',2); hold on;
                    end
                    plot(S(BoxID==term(j).BoxID).X,S(BoxID==term(j).BoxID).Y,'-y','linewidth',2); hold on;
                    plot(Snew(BoxIDnew==term(j).BoxID).X,Snew(BoxIDnew==term(j).BoxID).Y,'-c','linewidth',2); hold on;
                    xmins = [min(term(j).RGIX) min(term(j).BoxX) min(s(BoxIDold==term(j).BoxID).X) min(S(BoxID==term(j).BoxID).X) min(Snew(BoxIDnew==term(j).BoxID).X)];
                    xmaxs = [max(term(j).RGIX) max(term(j).BoxX) max(s(BoxIDold==term(j).BoxID).X) max(S(BoxID==term(j).BoxID).X) max(Snew(BoxIDnew==term(j).BoxID).X)];
                    ymins = [min(term(j).RGIY) min(term(j).BoxY) min(s(BoxIDold==term(j).BoxID).Y) min(S(BoxID==term(j).BoxID).Y) min(Snew(BoxIDnew==term(j).BoxID).Y)];
                    ymaxs = [max(term(j).RGIY) max(term(j).BoxY) max(s(BoxIDold==term(j).BoxID).Y) max(S(BoxID==term(j).BoxID).Y) max(Snew(BoxIDnew==term(j).BoxID).Y)];
                    set(gca,'xlim',[nanmin(xmins) nanmax(xmaxs)],'ylim',[nanmin(ymins) nanmax(ymaxs)]);
                    drawnow;
                    
                    %check that the image is useable
                    prompt = 'Is the terminus cloud-covered (y/n)?';
                    str = input(prompt,'s');
                    if strmatch(str,'n')==1
                        checked(j) = 1;
                        disp('...use image to map centerline & draw the flux gate');
                        
                        %manually zoom in more
                        beep;
                        disp('in the image, click on the UL & LR corners of a region where you''ll draw the centerline');
                        figure(oldfig); [a] = ginput(2); zoom_xlims = a(:,1); zoom_ylims = a(:,2); clear a;
                        set(gca,'xlim',sort(zoom_xlims),'ylim',sort(zoom_ylims));
                        
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
                        s.RGIId = term(j).RGIref;
                        s.BoxID = term(j).BoxID;
                        if j<10
                            shapewrite(s,['centerline_00',num2str(j),'.shp']);
                        elseif j >=10 && j<100
                            shapewrite(s,['centerline_0',num2str(j),'.shp']);
                        else
                            shapewrite(s,['centerline_',num2str(j),'.shp']);
                        end
                        eval(cd_to_root);
                        save('Greenland_GIC_centerlines.mat','term','checked','-v7.3');
                        clear s profx* profy*;
                        
                        disp('Extracting centerline terminus positions...');
                        %1985 terminus
                        if ~isempty(find(BoxIDold==term(j).BoxID))
                            termref = find(BoxIDold==term(j).BoxID);
                            [xi,yi,ii] = polyxpoly(term(j).centerX,term(j).centerY,s(termref).X,s(termref).Y);
                            term(j).X(1) = xi(1); term(j).Y(1) = yi(1);
                            term(j).center_termref(1) = ii(1);
                            term(j).year(1) = s(termref).year; term(j).doy(1) = s(termref).DOY;
                            clear termref xi yi ii;
                        else
                            term(j).X(1) = NaN; term(j).Y(1) = NaN;
                            term(j).center_termref(1) = NaN;
                            term(j).year(1) = NaN; term(j).doy(1) = NaN;
                        end
                        %2000 terminus
                        termref = find(BoxID==term(j).BoxID);
                        [xi,yi,ii] = polyxpoly(term(j).centerX,term(j).centerY,S(termref).X,S(termref).Y);
                        term(j).X(2) = xi(1); term(j).Y(2) = yi(1);
                        term(j).center_termref(2) = ii(1);
                        term(j).year(2) = S(termref).year; term(j).doy(2) = S(termref).DOY;
                        clear termref xi yi ii;
                        %2015 terminus
                        termref = find(BoxIDnew==term(j).BoxID);
                        [xi,yi,ii] = polyxpoly(term(j).centerX,term(j).centerY,Snew(termref).X,Snew(termref).Y);
                        term(j).X(3) = xi(1); term(j).Y(3) = yi(1);
                        term(j).center_termref(3) = ii(1);
                        term(j).year(3) = Snew(termref).year; term(j).doy(3) = Snew(termref).DOY;
                        clear termref xi yi ii;
                        %replace 1985 terminus with NaN if = 2000 terminus
                        if sqrt((term(j).X(1)-term(j).X(2))^2 + (term(j).Y(1)-term(j).Y(2))^2) == 0
                            term(j).X(1) = NaN; term(j).Y(1) = NaN; term(j).center_termref(1) = NaN;
                            term(j).year(1) = NaN; term(j).doy(1) = NaN;
                        end
                        
                        
                        %automatically generate a new terminus box if the old one doesn't contain all terminus traces
                        prompt = 'Does the terminus box need to be re-made (y/n)?';
                        str = input(prompt,'s');
                        if strmatch(str,'y')==1
                            %plot a new terminus box using the centerline & terminus traces as guides
                            termref_bounds = [min(term(j).center_termref) max(term(j).center_termref)];
                            width_decrease = 120; %m (subtracted from min width)
                            %                             length_increase = 2000; %m (half added to each end)
                            box_orient = atan2d(-(term(j).centerY(termref_bounds(2))-term(j).centerY(termref_bounds(1))),-(term(j).centerX(termref_bounds(2))-term(j).centerX(termref_bounds(1)))); %angle of centerline between inland- and seaward-most terminus
                            %compute the box width as double the minimum distance between the centerline and terminus trace end
                            if ~isempty(find(BoxIDold==term(j).BoxID))
                                box_width = 2*min([sqrt((s(BoxIDold==term(j).BoxID).X(1)-term(j).X(1)).^2+(s(BoxIDold==term(j).BoxID).Y(1)-term(j).Y(1)).^2);...
                                    sqrt((S(BoxID==term(j).BoxID).X(1)-term(j).X(2)).^2+(S(BoxID==term(j).BoxID).Y(1)-term(j).Y(2)).^2);...
                                    sqrt((Snew(BoxIDnew==term(j).BoxID).X(1)-term(j).X(3)).^2+(Snew(BoxIDnew==term(j).BoxID).Y(1)-term(j).Y(3)).^2);...
                                    sqrt((term(j).X(1)-s(BoxIDold==term(j).BoxID).X(end-1)).^2+(term(j).Y(1)-s(BoxIDold==term(j).BoxID).Y(end-1)).^2);...
                                    sqrt((term(j).X(2)-S(BoxID==term(j).BoxID).X(end-1)).^2+(term(j).Y(2)-S(BoxID==term(j).BoxID).Y(end-1)).^2);...
                                    sqrt((term(j).X(3)-Snew(BoxIDnew==term(j).BoxID).X(end-1)).^2+(term(j).Y(3)-Snew(BoxIDnew==term(j).BoxID).Y(end-1)).^2)]);
                            else
                                box_width = 2*min([sqrt((S(BoxID==term(j).BoxID).X(1)-term(j).X(2)).^2+(S(BoxID==term(j).BoxID).Y(1)-term(j).Y(2)).^2);...
                                    sqrt((Snew(BoxIDnew==term(j).BoxID).X(1)-term(j).X(3)).^2+(Snew(BoxIDnew==term(j).BoxID).Y(1)-term(j).Y(3)).^2);...
                                    sqrt((term(j).X(2)-S(BoxID==term(j).BoxID).X(end-1)).^2+(term(j).Y(2)-S(BoxID==term(j).BoxID).Y(end-1)).^2);...
                                    sqrt((term(j).X(3)-Snew(BoxIDnew==term(j).BoxID).X(end-1)).^2+(term(j).Y(3)-Snew(BoxIDnew==term(j).BoxID).Y(end-1)).^2)]);
                            end
                            if box_width > 620
                                box_width = box_width-width_decrease; %crop box slightly if it is decently large
                            else
                                box_width = box_width-30;
                            end
                            length_increase = 2*box_width;
                            box_length = sqrt((term(j).centerX(termref_bounds(2))-term(j).centerX(termref_bounds(1))).^2+(term(j).centerY(termref_bounds(2))-term(j).centerY(termref_bounds(1))).^2)+length_increase;
                            boxXo = term(j).centerX(termref_bounds(2)) - 0.5*length_increase*cosd(box_orient) + 0.5*box_width*cosd(90-box_orient); %inland right x coordinate
                            boxYo = term(j).centerY(termref_bounds(2)) - 0.5*length_increase*sind(box_orient) - 0.5*box_width*sind(90-box_orient); %inland right y coordinate
                            polyin = polyshape([boxXo boxXo+box_length boxXo+box_length boxXo],...
                                [boxYo boxYo boxYo+box_width boxYo+box_width]);
                            term_box = rotate(polyin,box_orient,[boxXo,boxYo]);
                            term(j).BoxX = []; term(j).BoxY = []; %clear out old box
                            term(j).BoxX = [term_box.Vertices(:,1)' term_box.Vertices(1,1)];
                            term(j).BoxY = [term_box.Vertices(:,2)' term_box.Vertices(1,2)];
                            plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'-k','linewidth',2); hold on;
                            plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'--w','linewidth',2); hold on;
                            
                            
                            %NOTE: THE TERMINUS BOX MAY NOT BE IDEAL FOR THE
                            %'BOX METHOD' BECAUSE SOME GLACIERS CURVE A LOT &
                            %NOT ALL TERMINUS TRACES MAY INTERSECT THE BOX ON BOTH SIDES
                            %BUT IT SHOULD HOPEFULLY BE GOOD ENOUGH FOR CREATING
                            %THE TERMINUS CHANGE TIME SERIES USING WTMM
                            term(j).BoxFlag = 1; %retraced has a flag of 1
                            clear termref_bounds box_* box*o polyin term_box;
                        else
                            term(j).BoxFlag = 0; %original has a flag of 0
                        end
                        cd_to_termini = ['cd ',root_dir,'terminus_traces/']; eval(cd_to_termini);
                        figure(newfig);
                        plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'-k','linewidth',2); hold on;
                        plot(term(j).centerX,term(j).centerY,'--g','linewidth',1); hold on;
                        if ~isempty(find(BoxIDold==term(j).BoxID))
                            plot(s(BoxIDold==term(j).BoxID).X,s(BoxIDold==term(j).BoxID).Y,'-m','linewidth',2); hold on;
                        end
                        plot(S(BoxID==term(j).BoxID).X,S(BoxID==term(j).BoxID).Y,'-y','linewidth',2); hold on;
                        plot(Snew(BoxIDnew==term(j).BoxID).X,Snew(BoxIDnew==term(j).BoxID).Y,'-c','linewidth',2); hold on;
                        set(gca,'xlim',[min(xmins) max(xmaxs)],'ylim',[min(ymins) max(ymaxs)]);
                        drawnow;
                        saveas(newfig,['BoxID',num2str(term(j).BoxID),'terminus-map.png'],'png');
                        eval(cd_to_root);
                        save('Greenland_GIC_centerlines.mat','term','checked','-v7.3');
                        
                    end
                end
                clear RGI_bbox;  
            end
            clear in_* im_cropped im_data;
        end
    end
    clear im *_sub *subgrid masked*;
    close all;
    cd_to_datedir = ['cd ',root_dir,add_datedir]; eval(cd_to_datedir);
end

%loop through the Landsat 4 images
for i = 1:length(L4s) %DONE!!!
    disp(['Landsat scene #',num2str(i),' of ',num2str(length(L4s)),': ',L4s(i).name]);
    cd_to_dir = ['cd ',L4s(i).name]; eval(cd_to_dir);
    bands = dir('*B*PS.tif');
    [I,R] = readgeoraster(bands(3).name);
    im.x = R.XWorldLimits(1):R.SampleSpacingInWorldX:R.XWorldLimits(2);
    im.y = R.YWorldLimits(2):-R.SampleSpacingInWorldY:R.YWorldLimits(1);
    im.z = double(I); clear I R;
    im_scaled = (im.z-min(min(im.z)))./(max(max(im.z))-min(min(im.z)));
    im.z_adjust = imadjust(im_scaled,[0 0.9],[0 1],0.7); clear im_scaled;
    oldfig = figure; set(gcf,'position',[50 50 800 800])
    imagesc(im.x,im.y,im.z_adjust); colormap gray; axis xy equal; hold on;
    newfig = figure; set(gcf,'position',[850 50 800 800])
    imagesc(im.x,im.y,im.z_adjust); colormap gray; axis xy equal; hold on;
    [imx,imy] = meshgrid(im.x,im.y);
    drawnow;
    
    %crop velocities to the RGI outline & overlay as a quiver plot on image
    disp('Prepping velocity vectors to aid mapping of the flux gate at the grounding line, perpendicular to flow...');
    vx_sub = V.vx(find(V.y>=min(im.y) & V.y<=max(im.y)),find(V.x>=min(im.x) & V.x<=max(im.x))); %crop vx to the landsat image extent
    vy_sub = V.vy(find(V.y>=min(im.y) & V.y<=max(im.y)),find(V.x>=min(im.x) & V.x<=max(im.x))); %crop vy to the landsat image extent
    x_sub = V.x(1,find(V.x>=min(im.x) & V.x<=max(im.x)));
    y_sub = V.y(1,find(V.y>=min(im.y) & V.y<=max(im.y)));
    [Xsubgrid,Ysubgrid] = meshgrid(x_sub,y_sub);
    in = inpolygon(Xsubgrid,Ysubgrid,[term(j).RGIX(1:end-1) term(j).RGIX(1)],[term(j).RGIY(1:end-1) term(j).RGIY(1)]); %mask velocities outside the glacier so velocity arrows scale approrpriately
    maskedvx = vx_sub; maskedvx(in==0) = NaN; %mask-out other regions by turning anything not in=1 into NaNs
    maskedvy = vy_sub; maskedvy(in==0) = NaN; %mask-out other regions by turning anything not in=1 into NaNs
    clear in;
    figure(oldfig);
    q = quiver(Xsubgrid,Ysubgrid,maskedvx,maskedvy,'-k'); hold on; %q.AutoScaleFactor = 50;
    
    %loop through the glaciers
    for j = 1:length(term)
        %check that the terminus box is within the image
        if nanmean(term(j).BoxX) > min(im.x) && nanmean(term(j).BoxX) < max(im.x) && nanmean(term(j).BoxY) > min(im.y) && nanmean(term(j).BoxY) < max(im.y)
            %identify image pixels in the terminus box
            in_termbox = inpolygon(imx,imy,term(j).BoxX(~isnan(term(j).BoxX)),term(j).BoxY(~isnan(term(j).BoxY)));
            im_cropped = in_termbox.*im.z;
            %             %identify pixels in the RGI outline
            %             in_RGI = inpolygon(imx,imy,term(j).RGIX(~isnan(term(j).RGIX)),term(j).RGIY(~isnan(term(j).RGIY)));
            %             im_cropped = in_RGI.*im.z;
            
            %check that the terminus overlies an area with data
            if sum(sum(im_cropped)) > 0
                %if the RGI outline is ENORMOUS, skip identifying pixels in the RGI polygon (takes too long)
                RGI_bbox = (max(term(j).RGIX)-min(term(j).RGIX))*(max(term(j).RGIY)-min(RGI(term(j).RGIY)));
                if RGI_bbox <= 1e+09
                    in_RGI = inpolygon(imx,imy,[term(j).RGIX(1:end-1) term(j).RGIX(1)],[term(j).RGIY(1:end-1) term(j).RGIY(1)]);
                    im_data = zeros(size(im.z_adjust)); im_data(in_RGI==1 & im.z_adjust~=0) = 1;
                end
                
                %add the centerline and fix the terminus box as necessary
                if (RGI_bbox > 1e+09 || sum(sum(im_data)) == sum(sum(in_RGI)))  && checked(j) == 0
                    figure(oldfig);
                    plot(term(j).RGIX, term(j).RGIY, '-r','linewidth',2); hold on;
                    plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'-k','linewidth',2); hold on; plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'--w','linewidth',2); hold on;
                    if ~isempty(find(BoxIDold==term(j).BoxID))
                        plot(s(BoxIDold==term(j).BoxID).X,s(BoxIDold==term(j).BoxID).Y,'-m','linewidth',2); hold on;
                    end
                    plot(S(BoxID==term(j).BoxID).X,S(BoxID==term(j).BoxID).Y,'-y','linewidth',2); hold on;
                    plot(Snew(BoxIDnew==term(j).BoxID).X,Snew(BoxIDnew==term(j).BoxID).Y,'-c','linewidth',2); hold on;
                    xmins = [min(term(j).RGIX) min(term(j).BoxX) min(s(BoxIDold==term(j).BoxID).X) min(S(BoxID==term(j).BoxID).X) min(Snew(BoxIDnew==term(j).BoxID).X)];
                    xmaxs = [max(term(j).RGIX) max(term(j).BoxX) max(s(BoxIDold==term(j).BoxID).X) max(S(BoxID==term(j).BoxID).X) max(Snew(BoxIDnew==term(j).BoxID).X)];
                    ymins = [min(term(j).RGIY) min(term(j).BoxY) min(s(BoxIDold==term(j).BoxID).Y) min(S(BoxID==term(j).BoxID).Y) min(Snew(BoxIDnew==term(j).BoxID).Y)];
                    ymaxs = [max(term(j).RGIY) max(term(j).BoxY) max(s(BoxIDold==term(j).BoxID).Y) max(S(BoxID==term(j).BoxID).Y) max(Snew(BoxIDnew==term(j).BoxID).Y)];
                    set(gca,'xlim',[nanmin(xmins) nanmax(xmaxs)],'ylim',[nanmin(ymins) nanmax(ymaxs)]);
                    drawnow;
                    
                    %check that the image is useable
                    prompt = 'Is the terminus cloud-covered (y/n)?';
                    str = input(prompt,'s');
                    if strmatch(str,'n')==1
                        checked(j) = 1;
                        disp('...use image to map centerline & draw the flux gate');
                        
                        %manually zoom in more
                        beep;
                        disp('in the image, click on the UL & LR corners of a region where you''ll draw the centerline');
                        figure(oldfig); [a] = ginput(2); zoom_xlims = a(:,1); zoom_ylims = a(:,2); clear a;
                        set(gca,'xlim',sort(zoom_xlims),'ylim',sort(zoom_ylims));
                        
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
                        s.RGIId = term(j).RGIref;
                        s.BoxID = term(j).BoxID;
                        if j<10
                            shapewrite(s,['centerline_00',num2str(j),'.shp']);
                        elseif j >=10 && j<100
                            shapewrite(s,['centerline_0',num2str(j),'.shp']);
                        else
                            shapewrite(s,['centerline_',num2str(j),'.shp']);
                        end
                        eval(cd_to_root);
                        save('Greenland_GIC_centerlines.mat','term','checked','-v7.3');
                        clear s profx* profy*;
                        
                        disp('Extracting centerline terminus positions...');
                        %1985 terminus
                        if ~isempty(find(BoxIDold==term(j).BoxID))
                            termref = find(BoxIDold==term(j).BoxID);
                            [xi,yi,ii] = polyxpoly(term(j).centerX,term(j).centerY,s(termref).X,s(termref).Y);
                            term(j).X(1) = xi(1); term(j).Y(1) = yi(1);
                            term(j).center_termref(1) = ii(1);
                            term(j).year(1) = s(termref).year; term(j).doy(1) = s(termref).DOY;
                            clear termref xi yi ii;
                        else
                            term(j).X(1) = NaN; term(j).Y(1) = NaN;
                            term(j).center_termref(1) = NaN;
                            term(j).year(1) = NaN; term(j).doy(1) = NaN;
                        end
                        %2000 terminus
                        termref = find(BoxID==term(j).BoxID);
                        [xi,yi,ii] = polyxpoly(term(j).centerX,term(j).centerY,S(termref).X,S(termref).Y);
                        term(j).X(2) = xi(1); term(j).Y(2) = yi(1);
                        term(j).center_termref(2) = ii(1);
                        term(j).year(2) = S(termref).year; term(j).doy(2) = S(termref).DOY;
                        clear termref xi yi ii;
                        %2015 terminus
                        termref = find(BoxIDnew==term(j).BoxID);
                        [xi,yi,ii] = polyxpoly(term(j).centerX,term(j).centerY,Snew(termref).X,Snew(termref).Y);
                        term(j).X(3) = xi(1); term(j).Y(3) = yi(1);
                        term(j).center_termref(3) = ii(1);
                        term(j).year(3) = Snew(termref).year; term(j).doy(3) = Snew(termref).DOY;
                        clear termref xi yi ii;
                        %replace 1985 terminus with NaN if = 2000 terminus
                        if sqrt((term(j).X(1)-term(j).X(2))^2 + (term(j).Y(1)-term(j).Y(2))^2) == 0
                            term(j).X(1) = NaN; term(j).Y(1) = NaN; term(j).center_termref(1) = NaN;
                            term(j).year(1) = NaN; term(j).doy(1) = NaN;
                        end
                        
                        
                        %automatically generate a new terminus box if the old one doesn't contain all terminus traces
                        prompt = 'Does the terminus box need to be re-made (y/n)?';
                        str = input(prompt,'s');
                        if strmatch(str,'y')==1
                            %plot a new terminus box using the centerline & terminus traces as guides
                            termref_bounds = [min(term(j).center_termref) max(term(j).center_termref)];
                            width_decrease = 120; %m (subtracted from min width)
                            %                             length_increase = 2000; %m (half added to each end)
                            box_orient = atan2d(-(term(j).centerY(termref_bounds(2))-term(j).centerY(termref_bounds(1))),-(term(j).centerX(termref_bounds(2))-term(j).centerX(termref_bounds(1)))); %angle of centerline between inland- and seaward-most terminus
                            %compute the box width as double the minimum distance between the centerline and terminus trace end
                            if ~isempty(find(BoxIDold==term(j).BoxID))
                                box_width = 2*min([sqrt((s(BoxIDold==term(j).BoxID).X(1)-term(j).X(1)).^2+(s(BoxIDold==term(j).BoxID).Y(1)-term(j).Y(1)).^2);...
                                    sqrt((S(BoxID==term(j).BoxID).X(1)-term(j).X(2)).^2+(S(BoxID==term(j).BoxID).Y(1)-term(j).Y(2)).^2);...
                                    sqrt((Snew(BoxIDnew==term(j).BoxID).X(1)-term(j).X(3)).^2+(Snew(BoxIDnew==term(j).BoxID).Y(1)-term(j).Y(3)).^2);...
                                    sqrt((term(j).X(1)-s(BoxIDold==term(j).BoxID).X(end-1)).^2+(term(j).Y(1)-s(BoxIDold==term(j).BoxID).Y(end-1)).^2);...
                                    sqrt((term(j).X(2)-S(BoxID==term(j).BoxID).X(end-1)).^2+(term(j).Y(2)-S(BoxID==term(j).BoxID).Y(end-1)).^2);...
                                    sqrt((term(j).X(3)-Snew(BoxIDnew==term(j).BoxID).X(end-1)).^2+(term(j).Y(3)-Snew(BoxIDnew==term(j).BoxID).Y(end-1)).^2)]);
                            else
                                box_width = 2*min([sqrt((S(BoxID==term(j).BoxID).X(1)-term(j).X(2)).^2+(S(BoxID==term(j).BoxID).Y(1)-term(j).Y(2)).^2);...
                                    sqrt((Snew(BoxIDnew==term(j).BoxID).X(1)-term(j).X(3)).^2+(Snew(BoxIDnew==term(j).BoxID).Y(1)-term(j).Y(3)).^2);...
                                    sqrt((term(j).X(2)-S(BoxID==term(j).BoxID).X(end-1)).^2+(term(j).Y(2)-S(BoxID==term(j).BoxID).Y(end-1)).^2);...
                                    sqrt((term(j).X(3)-Snew(BoxIDnew==term(j).BoxID).X(end-1)).^2+(term(j).Y(3)-Snew(BoxIDnew==term(j).BoxID).Y(end-1)).^2)]);
                            end
                            if box_width > 620
                                box_width = box_width-width_decrease; %crop box slightly if it is decently large
                            else
                                box_width = box_width-30;
                            end
                            length_increase = 2*box_width;
                            box_length = sqrt((term(j).centerX(termref_bounds(2))-term(j).centerX(termref_bounds(1))).^2+(term(j).centerY(termref_bounds(2))-term(j).centerY(termref_bounds(1))).^2)+length_increase;
                            boxXo = term(j).centerX(termref_bounds(2)) - 0.5*length_increase*cosd(box_orient) + 0.5*box_width*cosd(90-box_orient); %inland right x coordinate
                            boxYo = term(j).centerY(termref_bounds(2)) - 0.5*length_increase*sind(box_orient) - 0.5*box_width*sind(90-box_orient); %inland right y coordinate
                            polyin = polyshape([boxXo boxXo+box_length boxXo+box_length boxXo],...
                                [boxYo boxYo boxYo+box_width boxYo+box_width]);
                            term_box = rotate(polyin,box_orient,[boxXo,boxYo]);
                            term(j).BoxX = []; term(j).BoxY = []; %clear out old box
                            term(j).BoxX = [term_box.Vertices(:,1)' term_box.Vertices(1,1)];
                            term(j).BoxY = [term_box.Vertices(:,2)' term_box.Vertices(1,2)];
                            plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'-k','linewidth',2); hold on;
                            plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'--w','linewidth',2); hold on;
                            
                            
                            %NOTE: THE TERMINUS BOX MAY NOT BE IDEAL FOR THE
                            %'BOX METHOD' BECAUSE SOME GLACIERS CURVE A LOT &
                            %NOT ALL TERMINUS TRACES MAY INTERSECT THE BOX ON BOTH SIDES
                            %BUT IT SHOULD HOPEFULLY BE GOOD ENOUGH FOR CREATING
                            %THE TERMINUS CHANGE TIME SERIES USING WTMM
                            term(j).BoxFlag = 1; %retraced has a flag of 1
                            clear termref_bounds box_* box*o polyin term_box;
                        else
                            term(j).BoxFlag = 0; %original has a flag of 0
                        end
                        cd_to_termini = ['cd ',root_dir,'terminus_traces/']; eval(cd_to_termini);
                        figure(newfig);
                        plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'-k','linewidth',2); hold on;
                        plot(term(j).centerX,term(j).centerY,'--g','linewidth',1); hold on;
                        if ~isempty(find(BoxIDold==term(j).BoxID))
                            plot(s(BoxIDold==term(j).BoxID).X,s(BoxIDold==term(j).BoxID).Y,'-m','linewidth',2); hold on;
                        end
                        plot(S(BoxID==term(j).BoxID).X,S(BoxID==term(j).BoxID).Y,'-y','linewidth',2); hold on;
                        plot(Snew(BoxIDnew==term(j).BoxID).X,Snew(BoxIDnew==term(j).BoxID).Y,'-c','linewidth',2); hold on;
                        set(gca,'xlim',[min(xmins) max(xmaxs)],'ylim',[min(ymins) max(ymaxs)]);
                        drawnow;
                        saveas(newfig,['BoxID',num2str(term(j).BoxID),'terminus-map.png'],'png');
                        eval(cd_to_root);
                        save('Greenland_GIC_centerlines.mat','term','checked','-v7.3');
                        
                    end
                end
                clear RGI_bbox;  
            end
            clear in_* im_cropped im_data;
        end
    end
    clear im *_sub *subgrid masked*;
    close all;
    cd_to_datedir = ['cd ',root_dir,add_datedir]; eval(cd_to_datedir);
end

%loop through the Landsat 7 images
for i = 1:length(L7s) %DONE!!!
    disp(['Landsat scene #',num2str(i),' of ',num2str(length(L7s)),': ',L7s(i).name]);
    cd_to_dir = ['cd ',L7s(i).name]; eval(cd_to_dir);
    bands = dir('*B*PS.tif');
    [I,R] = readgeoraster(bands(9).name);
    im.x = R.XWorldLimits(1):R.SampleSpacingInWorldX:R.XWorldLimits(2);
    im.y = R.YWorldLimits(2):-R.SampleSpacingInWorldY:R.YWorldLimits(1);
    im.z = double(I); clear I R;
    im_scaled = (im.z-min(min(im.z)))./(max(max(im.z))-min(min(im.z)));
    im.z_adjust = imadjust(im_scaled,[0 0.9],[0 1],0.7); clear im_scaled;
    oldfig = figure; set(gcf,'position',[50 50 800 800])
    imagesc(im.x,im.y,im.z_adjust); colormap gray; axis xy equal; hold on;
    newfig = figure; set(gcf,'position',[850 50 800 800])
    imagesc(im.x,im.y,im.z_adjust); colormap gray; axis xy equal; hold on;
    [imx,imy] = meshgrid(im.x,im.y);
    drawnow;
    
    %crop velocities to the RGI outline & overlay as a quiver plot on image
    disp('Prepping velocity vectors to aid mapping of the flux gate at the grounding line, perpendicular to flow...');
    vx_sub = V.vx(find(V.y>=min(im.y) & V.y<=max(im.y)),find(V.x>=min(im.x) & V.x<=max(im.x))); %crop vx to the landsat image extent
    vy_sub = V.vy(find(V.y>=min(im.y) & V.y<=max(im.y)),find(V.x>=min(im.x) & V.x<=max(im.x))); %crop vy to the landsat image extent
    x_sub = V.x(1,find(V.x>=min(im.x) & V.x<=max(im.x)));
    y_sub = V.y(1,find(V.y>=min(im.y) & V.y<=max(im.y)));
    [Xsubgrid,Ysubgrid] = meshgrid(x_sub,y_sub);
    in = inpolygon(Xsubgrid,Ysubgrid,[term(j).RGIX(1:end-1) term(j).RGIX(1)],[term(j).RGIY(1:end-1) term(j).RGIY(1)]); %mask velocities outside the glacier so velocity arrows scale approrpriately
    maskedvx = vx_sub; maskedvx(in==0) = NaN; %mask-out other regions by turning anything not in=1 into NaNs
    maskedvy = vy_sub; maskedvy(in==0) = NaN; %mask-out other regions by turning anything not in=1 into NaNs
    clear in;
    figure(oldfig);
    q = quiver(Xsubgrid,Ysubgrid,maskedvx,maskedvy,'-k'); hold on; %q.AutoScaleFactor = 50;
    
    %loop through the glaciers
    for j = 1:length(term)
        %check that the terminus box is within the image
        if nanmean(term(j).BoxX) > min(im.x) && nanmean(term(j).BoxX) < max(im.x) && nanmean(term(j).BoxY) > min(im.y) && nanmean(term(j).BoxY) < max(im.y)
            %identify image pixels in the terminus box
            in_termbox = inpolygon(imx,imy,term(j).BoxX(~isnan(term(j).BoxX)),term(j).BoxY(~isnan(term(j).BoxY)));
            im_cropped = in_termbox.*im.z;
            %             %identify pixels in the RGI outline
            %             in_RGI = inpolygon(imx,imy,term(j).RGIX(~isnan(term(j).RGIX)),term(j).RGIY(~isnan(term(j).RGIY)));
            %             im_cropped = in_RGI.*im.z;
            
            %check that the terminus overlies an area with data
            if sum(sum(im_cropped)) > 0
                %if the RGI outline is ENORMOUS, skip identifying pixels in the RGI polygon (takes too long)
                RGI_bbox = (max(term(j).RGIX)-min(term(j).RGIX))*(max(term(j).RGIY)-min(RGI(term(j).RGIY)));
                if RGI_bbox <= 1e+09
                    in_RGI = inpolygon(imx,imy,[term(j).RGIX(1:end-1) term(j).RGIX(1)],[term(j).RGIY(1:end-1) term(j).RGIY(1)]);
                    im_data = zeros(size(im.z_adjust)); im_data(in_RGI==1 & im.z_adjust~=0) = 1;
                end
                
                %add the centerline and fix the terminus box as necessary
                if (RGI_bbox > 1e+09 || sum(sum(im_data)) == sum(sum(in_RGI)))  && checked(j) == 0
                    figure(oldfig);
                    plot(term(j).RGIX, term(j).RGIY, '-r','linewidth',2); hold on;
                    plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'-k','linewidth',2); hold on; plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'--w','linewidth',2); hold on;
                    if ~isempty(find(BoxIDold==term(j).BoxID))
                     plot(s(BoxIDold==term(j).BoxID).X,s(BoxIDold==term(j).BoxID).Y,'-m','linewidth',2); hold on;
                    end
                    plot(S(BoxID==term(j).BoxID).X,S(BoxID==term(j).BoxID).Y,'-y','linewidth',2); hold on;
                    plot(Snew(BoxIDnew==term(j).BoxID).X,Snew(BoxIDnew==term(j).BoxID).Y,'-c','linewidth',2); hold on;
                    xmins = [min(term(j).RGIX) min(term(j).BoxX) min(s(BoxIDold==term(j).BoxID).X) min(S(BoxID==term(j).BoxID).X) min(Snew(BoxIDnew==term(j).BoxID).X)];
                    xmaxs = [max(term(j).RGIX) max(term(j).BoxX) max(s(BoxIDold==term(j).BoxID).X) max(S(BoxID==term(j).BoxID).X) max(Snew(BoxIDnew==term(j).BoxID).X)];
                    ymins = [min(term(j).RGIY) min(term(j).BoxY) min(s(BoxIDold==term(j).BoxID).Y) min(S(BoxID==term(j).BoxID).Y) min(Snew(BoxIDnew==term(j).BoxID).Y)];
                    ymaxs = [max(term(j).RGIY) max(term(j).BoxY) max(s(BoxIDold==term(j).BoxID).Y) max(S(BoxID==term(j).BoxID).Y) max(Snew(BoxIDnew==term(j).BoxID).Y)];
                    set(gca,'xlim',[nanmin(xmins) nanmax(xmaxs)],'ylim',[nanmin(ymins) nanmax(ymaxs)]);
                    drawnow;
                    
                    %check that the image is useable
                    prompt = 'Is the terminus cloud-covered (y/n)?';
                    str = input(prompt,'s');
                    if strmatch(str,'n')==1
                        checked(j) = 1;
                        disp('...use image to map centerline & draw the flux gate');
                        
                        %manually zoom in more
                        beep;
                        disp('in the image, click on the UL & LR corners of a region where you''ll draw the centerline');
                        figure(oldfig); [a] = ginput(2); zoom_xlims = a(:,1); zoom_ylims = a(:,2); clear a;
                        set(gca,'xlim',sort(zoom_xlims),'ylim',sort(zoom_ylims));
                        
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
                        s.RGIId = term(j).RGIref;
                        s.BoxID = term(j).BoxID;
                        if j<10
                            shapewrite(s,['centerline_00',num2str(j),'.shp']);
                        elseif j >=10 && j<100
                            shapewrite(s,['centerline_0',num2str(j),'.shp']);
                        else
                            shapewrite(s,['centerline_',num2str(j),'.shp']);
                        end
                        eval(cd_to_root);
                        save('Greenland_GIC_centerlines.mat','term','checked','-v7.3');
                        clear s profx* profy*;
                        
                        disp('Extracting centerline terminus positions...');
                        %1985 terminus
                        if ~isempty(find(BoxIDold==term(j).BoxID))
                            termref = find(BoxIDold==term(j).BoxID);
                            [xi,yi,ii] = polyxpoly(term(j).centerX,term(j).centerY,s(termref).X,s(termref).Y);
                            term(j).X(1) = xi(1); term(j).Y(1) = yi(1);
                            term(j).center_termref(1) = ii(1);
                            term(j).year(1) = s(termref).year; term(j).doy(1) = s(termref).DOY;
                            clear termref xi yi ii;
                        else
                            term(j).X(1) = NaN; term(j).Y(1) = NaN;
                            term(j).center_termref(1) = NaN;
                            term(j).year(1) = NaN; term(j).doy(1) = NaN;
                        end
                        %2000 terminus
                        termref = find(BoxID==term(j).BoxID);
                        [xi,yi,ii] = polyxpoly(term(j).centerX,term(j).centerY,S(termref).X,S(termref).Y);
                        term(j).X(2) = xi(1); term(j).Y(2) = yi(1);
                        term(j).center_termref(2) = ii(1);
                        term(j).year(2) = S(termref).year; term(j).doy(2) = S(termref).DOY;
                        clear termref xi yi ii;
                        %2015 terminus
                        termref = find(BoxIDnew==term(j).BoxID);
                        [xi,yi,ii] = polyxpoly(term(j).centerX,term(j).centerY,Snew(termref).X,Snew(termref).Y);
                        term(j).X(3) = xi(1); term(j).Y(3) = yi(1);
                        term(j).center_termref(3) = ii(1);
                        term(j).year(3) = Snew(termref).year; term(j).doy(3) = Snew(termref).DOY;
                        clear termref xi yi ii;
                        %replace 1985 terminus with NaN if = 2000 terminus
                        if sqrt((term(j).X(1)-term(j).X(2))^2 + (term(j).Y(1)-term(j).Y(2))^2) == 0
                            term(j).X(1) = NaN; term(j).Y(1) = NaN; term(j).center_termref(1) = NaN;
                            term(j).year(1) = NaN; term(j).doy(1) = NaN;
                        end
                        
                        
                        %automatically generate a new terminus box if the old one doesn't contain all terminus traces
                        prompt = 'Does the terminus box need to be re-made (y/n)?';
                        str = input(prompt,'s');
                        if strmatch(str,'y')==1
                            %plot a new terminus box using the centerline & terminus traces as guides
                            termref_bounds = [min(term(j).center_termref) max(term(j).center_termref)];
                            width_decrease = 120; %m (subtracted from min width)
                            %                             length_increase = 2000; %m (half added to each end)
                            box_orient = atan2d(-(term(j).centerY(termref_bounds(2))-term(j).centerY(termref_bounds(1))),-(term(j).centerX(termref_bounds(2))-term(j).centerX(termref_bounds(1)))); %angle of centerline between inland- and seaward-most terminus
                            %compute the box width as double the minimum distance between the centerline and terminus trace end
                            if ~isempty(find(BoxIDold==term(j).BoxID))
                                box_width = 2*min([sqrt((s(BoxIDold==term(j).BoxID).X(1)-term(j).X(1)).^2+(s(BoxIDold==term(j).BoxID).Y(1)-term(j).Y(1)).^2);...
                                    sqrt((S(BoxID==term(j).BoxID).X(1)-term(j).X(2)).^2+(S(BoxID==term(j).BoxID).Y(1)-term(j).Y(2)).^2);...
                                    sqrt((Snew(BoxIDnew==term(j).BoxID).X(1)-term(j).X(3)).^2+(Snew(BoxIDnew==term(j).BoxID).Y(1)-term(j).Y(3)).^2);...
                                    sqrt((term(j).X(1)-s(BoxIDold==term(j).BoxID).X(end-1)).^2+(term(j).Y(1)-s(BoxIDold==term(j).BoxID).Y(end-1)).^2);...
                                    sqrt((term(j).X(2)-S(BoxID==term(j).BoxID).X(end-1)).^2+(term(j).Y(2)-S(BoxID==term(j).BoxID).Y(end-1)).^2);...
                                    sqrt((term(j).X(3)-Snew(BoxIDnew==term(j).BoxID).X(end-1)).^2+(term(j).Y(3)-Snew(BoxIDnew==term(j).BoxID).Y(end-1)).^2)]);
                            else
                                box_width = 2*min([sqrt((S(BoxID==term(j).BoxID).X(1)-term(j).X(2)).^2+(S(BoxID==term(j).BoxID).Y(1)-term(j).Y(2)).^2);...
                                    sqrt((Snew(BoxIDnew==term(j).BoxID).X(1)-term(j).X(3)).^2+(Snew(BoxIDnew==term(j).BoxID).Y(1)-term(j).Y(3)).^2);...
                                    sqrt((term(j).X(2)-S(BoxID==term(j).BoxID).X(end-1)).^2+(term(j).Y(2)-S(BoxID==term(j).BoxID).Y(end-1)).^2);...
                                    sqrt((term(j).X(3)-Snew(BoxIDnew==term(j).BoxID).X(end-1)).^2+(term(j).Y(3)-Snew(BoxIDnew==term(j).BoxID).Y(end-1)).^2)]);
                            end
                            if box_width > 620
                                box_width = box_width-width_decrease; %crop box slightly if it is decently large
                            else
                                box_width = box_width-30;
                            end
                            length_increase = 2*box_width;
                            box_length = sqrt((term(j).centerX(termref_bounds(2))-term(j).centerX(termref_bounds(1))).^2+(term(j).centerY(termref_bounds(2))-term(j).centerY(termref_bounds(1))).^2)+length_increase;
                            boxXo = term(j).centerX(termref_bounds(2)) - 0.5*length_increase*cosd(box_orient) + 0.5*box_width*cosd(90-box_orient); %inland right x coordinate
                            boxYo = term(j).centerY(termref_bounds(2)) - 0.5*length_increase*sind(box_orient) - 0.5*box_width*sind(90-box_orient); %inland right y coordinate
                            polyin = polyshape([boxXo boxXo+box_length boxXo+box_length boxXo],...
                                [boxYo boxYo boxYo+box_width boxYo+box_width]);
                            term_box = rotate(polyin,box_orient,[boxXo,boxYo]);
                            term(j).BoxX = []; term(j).BoxY = []; %clear out old box
                            term(j).BoxX = [term_box.Vertices(:,1)' term_box.Vertices(1,1)];
                            term(j).BoxY = [term_box.Vertices(:,2)' term_box.Vertices(1,2)];
                            plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'-k','linewidth',2); hold on;
                            plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'--w','linewidth',2); hold on;
                            
                            
                            %NOTE: THE TERMINUS BOX MAY NOT BE IDEAL FOR THE
                            %'BOX METHOD' BECAUSE SOME GLACIERS CURVE A LOT &
                            %NOT ALL TERMINUS TRACES MAY INTERSECT THE BOX ON BOTH SIDES
                            %BUT IT SHOULD HOPEFULLY BE GOOD ENOUGH FOR CREATING
                            %THE TERMINUS CHANGE TIME SERIES USING WTMM
                            term(j).BoxFlag = 1; %retraced has a flag of 1
                            clear termref_bounds box_* box*o polyin term_box;
                        else
                            term(j).BoxFlag = 0; %original has a flag of 0
                        end
                        cd_to_termini = ['cd ',root_dir,'terminus_traces/']; eval(cd_to_termini);
                        figure(newfig);
                        plot([term(j).BoxX term(j).BoxX(1)],[term(j).BoxY term(j).BoxY(1)],'-k','linewidth',2); hold on;
                        plot(term(j).centerX,term(j).centerY,'--g','linewidth',1); hold on;
                        if ~isempty(find(BoxIDold==term(j).BoxID))
                            plot(s(BoxIDold==term(j).BoxID).X,s(BoxIDold==term(j).BoxID).Y,'-m','linewidth',2); hold on;
                        end
                        plot(S(BoxID==term(j).BoxID).X,S(BoxID==term(j).BoxID).Y,'-y','linewidth',2); hold on;
                        plot(Snew(BoxIDnew==term(j).BoxID).X,Snew(BoxIDnew==term(j).BoxID).Y,'-c','linewidth',2); hold on;
                        set(gca,'xlim',[min(xmins) max(xmaxs)],'ylim',[min(ymins) max(ymaxs)]);
                        drawnow;
                        saveas(newfig,['BoxID',num2str(term(j).BoxID),'terminus-map.png'],'png');
                        eval(cd_to_root);
                        save('Greenland_GIC_centerlines.mat','term','checked','-v7.3');
                        
                    end
                end
                clear RGI_bbox;  
            end
            clear in_* im_cropped im_data;
        end
    end
    clear im *_sub *subgrid masked*;
    close all;
    cd_to_datedir = ['cd ',root_dir,add_datedir]; eval(cd_to_datedir);
end
disp('Centerlines added to term structure for all glaciers');

%use edit_centerlines_termboxes.m to edit centerlines & term boxes if necessary
disp('If you noticed bad centerlines or terminus boxes, use edit_centerlines_termboxes.m to check all & edit');

%% extract centerline profiles & identify the grounding line position along the centerline
close all;
disp('loading data for centerline profile extraction');
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
load Greenland_GIC_centerlines.mat;

%load the georeferencing info for the ArcticDEM (subset for each glacier in loop because it is HUGE)
ArcticDEMdir = [misc_dir,'ArcticDEM10m/'];
cd_to_arcticDEM = ['cd ',ArcticDEMdir]; eval(cd_to_arcticDEM);
info = geotiffinfo('ArcticDEM_10m_Greenland_Mosaic.tif'); %pull spatial referencing info for the FULL dem
ArcticDEMx = single([info.SpatialRef.XWorldLimits(1)+0.5*info.SpatialRef.CellExtentInWorldX:info.SpatialRef.CellExtentInWorldX:info.SpatialRef.XWorldLimits(2)-0.5*info.SpatialRef.CellExtentInWorldX]);
ArcticDEMy = single([info.SpatialRef.YWorldLimits(2)-0.5*info.SpatialRef.CellExtentInWorldY:-info.SpatialRef.CellExtentInWorldY:info.SpatialRef.YWorldLimits(1)+0.5*info.SpatialRef.CellExtentInWorldY]);

%load the GIMP DEM to help with visualization
cd_to_DEM = ['cd ',misc_dir,'GIMP/']; eval(cd_to_DEM);
[I,R] = readgeoraster('gimpdem_90m_v01.1.tif');
dem.x = single([R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX]);
dem.y = single([R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY]);
dem.z = single(I);
clear I R;
disp('...term structure & DEM data loaded');

disp('Extracting centerline data for...');
counter = 0;
for j = 1:length(term)
    counter = counter+1; %count iterations for saving purposes
    disp([num2str(j),' (BoxID=',num2str(term(j).BoxID),' & RGIref=',num2str(term(j).RGIref),')']);
    
    %check that the centerline distance vector is correct
    if length(term(j).centerline) ~= length(term(j).centerX)
        term(j).centerline = [];
        term(j).centerline(1) = 0;
        for k = 2:length(term(j).centerX)
            term(j).centerline(k) = term(j).centerline(k-1)+sqrt((term(j).centerX(k)-term(j).centerX(k-1)).^2 + (term(j).centerY(k)-term(j).centerY(k-1)).^2);
        end
    end
    
    %velocity profile time series
    cd_to_vels = ['cd ',misc_dir,'Greenland-ITSLIVE_1985-2018/']; eval(cd_to_vels);
    speed_fig = figure; set(speed_fig,'position',[200 50 800 400]);
    ITSLIVEs = dir('GRE_*.nc'); speed_cmap = colormap(parula(length(ITSLIVEs)));
    for k=1:length(ITSLIVEs)
        %load the velocity coordinates
        v_x = ncread(ITSLIVEs(k).name,'x');
        v_y = ncread(ITSLIVEs(k).name,'y');
        
        %determine the spatial subset over which to extract data
        vxrefs = find(v_x >= min(term(j).centerX) & v_x <= max(term(j).centerX));
        if isempty(vxrefs) && unique(sign(gradient(v_x))) == 1
            vxrefs = [find(v_x >= min(term(j).centerX),1,'first') find(v_x <= max(term(j).centerX),1,'last')];
        elseif isempty(vxrefs) && unique(sign(gradient(v_x))) == -1
            vxrefs = [find(v_x >= min(term(j).centerX),1,'last') find(v_x <= max(term(j).centerX),1,'first')];
        end
        vyrefs = find(v_y >= min(term(j).centerY) & v_y <= max(term(j).centerY));
        if isempty(vyrefs) && unique(sign(gradient(v_y))) == -1
            vyrefs = [find(v_y >= min(term(j).centerY),1,'last') find(v_y <= max(term(j).centerY),1,'first')];
        elseif isempty(vyrefs) && unique(sign(gradient(v_y))) == 1
            vyrefs = [find(v_y >= min(term(j).centerY),1,'first') find(v_y <= max(term(j).centerY),1,'last')];
        end
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
        if isempty(vyrefs) && unique(sign(gradient(v_y))) == -1
            vyrefs = [find(v_y >= min(term(j).centerY),1,'last') find(v_y <= max(term(j).centerY),1,'first')];
        elseif isempty(vyrefs) && unique(sign(gradient(v_y))) == 1
            vyrefs = [find(v_y >= min(term(j).centerY),1,'first') find(v_y <= max(term(j).centerY),1,'last')];
        end
        
        %interpolate the the centerline profile & add to structure
        [vxgrid,vygrid] = meshgrid(v_x(min(vxrefs)-1:max(vxrefs)+1),v_y(min(vyrefs)-1:max(vyrefs)+1));
        term(j).centerVdate(k,:) = interp2(vxgrid,vygrid,vdecidate,term(j).centerX,term(j).centerY,'nearest');
        term(j).centerVdt(k,:) = interp2(vxgrid,vygrid,vdecidt,term(j).centerX,term(j).centerY,'nearest');
        term(j).centerV(k,:) = interp2(vxgrid,vygrid,sqrt(vx.^2 + vy.^2),term(j).centerX,term(j).centerY);
        term(j).centerVdir(k,:) = interp2(vxgrid,vygrid,atand(vy./vx),term(j).centerX,term(j).centerY);
        term(j).centerVerr(k,:) = interp2(vxgrid,vygrid,verr,term(j).centerX,term(j).centerY);
        
        %plot
        figure(speed_fig);
        plot(term(j).centerline-term(j).centerline(min(term(j).center_termref)),term(j).centerV(k,:),'-','color',speed_cmap(k,:),'linewidth',1); hold on;
        clear v_x v_y vx* vy* v*date v*dt;
    end
    figure(speed_fig);
    legend(num2str([1985:1:2018]'),'location','eastoutside','NumColumns',2,'fontsize',12);
    set(gca,'fontsize',14,'xlim',[-200 max(term(j).centerline-term(j).centerline(min(term(j).center_termref)))]); grid on;
    xlabel('Distance from terminus (m)','fontsize',14); ylabel('Speed (m/yr)','fontsize',14);
    drawnow;
    %save the speed profiles
    cd_to_centerlines = ['cd ',root_dir,'centerlines/']; eval(cd_to_centerlines);
    if term(j).BoxID<10
        saveas(speed_fig,['00',num2str(term(j).BoxID),'_speed-profiles.png'],'png');
    elseif term(j).BoxID>=10 && term(j).BoxID<100
        saveas(speed_fig,['0',num2str(term(j).BoxID),'_speed-profiles.png'],'png');
    else
        saveas(speed_fig,[num2str(term(j).BoxID),'_speed-profiles.png'],'png');
    end
    
    
    %pull elevations for the profile from the 10 m-resolution ArcticDEM for the glacier
    cd_to_arcticDEM = ['cd ',ArcticDEMdir]; eval(cd_to_arcticDEM);
    dem_xsub = [find(ArcticDEMx<=min([term(j).RGIX term(j).BoxX]),1,'last'),find(ArcticDEMx>=max([term(j).RGIX term(j).BoxX]),1,'first')];
    dem_ysub = [find(ArcticDEMy>=max([term(j).RGIY term(j).BoxY]),1,'last'),find(ArcticDEMy<=min([term(j).RGIY term(j).BoxY]),1,'first')];
    I = imread('ArcticDEM_10m_Greenland_Mosaic.tif','PixelRegion',{dem_ysub,dem_xsub}); I(I==-9999) = NaN;
    ArcticDEMz = I; clear I;
    [ArcticDEM_xgrid,ArcticDEM_ygrid] = meshgrid(ArcticDEMx(min(dem_xsub):max(dem_xsub)),ArcticDEMy(min(dem_ysub):max(dem_ysub)));
    term(j).centerZdate(1,:) = 2008; term(j).centerZdate(2,:) = 2018;
    term(j).centerZ(1,:) = interp2(double(dem.x),double(dem.y),dem.z,double(term(j).centerX),double(term(j).centerY)); %GIMP elevation profile (from ~2008)
    term(j).centerZ(2,:) = interp2(double(ArcticDEMx(min(dem_xsub):max(dem_xsub))),double(ArcticDEMy(min(dem_ysub):max(dem_ysub))),ArcticDEMz,double(term(j).centerX),double(term(j).centerY)); %ArcticDEM profile (from ~2018)
    [lon,lat] = ps2wgs(term(j).centerX,term(j).centerY); geoid = geoidheight(lat,lon);
    %plot the elevation profiles & identify terminus positions
    prof_fig = figure; set(prof_fig,'position',[50 50 400 400]);
    figure(prof_fig);
    %plot the GIMP profile
    if size(term(j).centerZ(1,:),2) ~= size(geoid,2)
        plot(term(j).centerline-term(j).centerline(min(term(j).center_termref)),term(j).centerZ(1,:)-geoid','-b','linewidth',1); hold on;
    else
        plot(term(j).centerline-term(j).centerline(min(term(j).center_termref)),term(j).centerZ(1,:)-geoid,'-b','linewidth',1); hold on;
    end
    %plot the ArcticDEM profile
    if size(term(j).centerZ(2,:),2) ~= size(geoid,2)
        plot(term(j).centerline-term(j).centerline(min(term(j).center_termref)),term(j).centerZ(2,:)-geoid','-k','linewidth',2); hold on;
    else
        plot(term(j).centerline-term(j).centerline(min(term(j).center_termref)),term(j).centerZ(2,:)-geoid,'-k','linewidth',2); hold on;
    end
    legend('GIMP','ArcticDEM','location','northwest');
    clear lon lat;
    
    %plot the box-centerline intersection point on the DEM & image
    [box_interceptX,box_interceptY] = polyxpoly(term(j).centerX,term(j).centerY,term(j).BoxX,term(j).BoxY); %find where the centerline intersects the terminus box
    intercept_dist = sqrt((term(j).centerX(1)-box_interceptX).^2 + (term(j).centerY(1)-box_interceptY).^2); %approximate distance from the first centerline point
    [~,maxref] = max(intercept_dist); clear intercept_dist; %use the most inland intersection in case the seaward edge of the terminus box also intersects the centerline
    [~,minref] = min(sqrt((term(j).centerX-box_interceptX(maxref)).^2 + (term(j).centerY-box_interceptY(maxref)).^2));
    figure(prof_fig); plot(term(j).centerline(minref)-term(j).centerline(min(term(j).center_termref)),term(j).centerZ(2,minref)-geoid(minref),'xr','linewidth',2); hold on;
    legend('GIMP','ArcticDEM','box intersection','location','northwest');
    
    
    %decide if there is a clear slope break indicating floating ice at the terminus
%     prompt = 'Is there a flattening of slope below ~50m elevation suggesting floating ice INLAND of the box intercept (y/n)?';
%     str = input(prompt,'s');
%     if strmatch(str,'y')==1
%         disp('in the elevation profile, click on the slope break');
%         figure(prof_fig); b = ginput(1); %click on the point in the figure & output the centerline position to b
%         [~,prof_ref] = min(abs(term(j).centerline-b(1))); %find the reference for the closest neighboring centerline point
%         term(j).centerGL = [term(j).centerX(prof_ref) term(j).centerY(prof_ref) term(j).centerZ(2,prof_ref)-geoid(prof_ref)]; %x,y,z of centerline grounding line
%         plot(term(j).centerline(prof_ref),term(j).centerZ(2,prof_ref)-geoid(prof_ref),'+r','linewidth',2); hold on;
%         legend('GIMP','ArcticDEM','box intersection','grounding line','location','northwest');
%         clear b prof_ref;
%     else
%         disp('automatically assigning the centerline flux gate to the inland edge of the terminus box');
%         term(j).centerGL = [box_interceptX(maxref) box_interceptY(maxref) term(j).centerZ(2,minref)-geoid(minref)]; %x,y,z of centerline box intersection
%         figure(prof_fig); legend('GIMP','ArcticDEM','box intersection','location','northwest')
%     end
    figure(prof_fig); set(gca,'fontsize',14,'xlim',...
        [-200 max(term(j).centerline-term(j).centerline(min(term(j).center_termref)))],...
        'ylim',[0 max(max(term(j).centerZ-geoid'))]); grid on;
    xlabel('Distance from terminus (m)','fontsize',14); ylabel('Elevation (m a.s.l.)','fontsize',14);
    %save the elevation profiles
    eval(cd_to_centerlines);
    if term(j).BoxID<10
        saveas(prof_fig,['00',num2str(term(j).BoxID),'_elev-profiles.png'],'png');
    elseif term(j).BoxID>=10 && term(j).BoxID<100
        saveas(prof_fig,['0',num2str(term(j).BoxID),'_elev-profiles.png'],'png');
    else
        saveas(prof_fig,[num2str(term(j).BoxID),'_elev-profiles.png'],'png');
    end
    clear prof;
    %save the grounding line to the structure
    eval(cd_to_root);
    if counter == 10;
        counter = 0;
        save('Greenland_GIC_centerlines.mat','term','-v7.3');
        disp(['term structure resaved after ',num2str(j),' loops']);
    end
    close all;
end

%% draw flux gates
clearvars; close all;

disp('Draw flux gates roughly perpendicular to flow and spanning the glacier');
eval(cd_to_root);
load Greenland_GIC_centerlines.mat;

%load the velocity mosaic
cd_to_vels = ['cd ',misc_dir,'Greenland-VelMosaic_1995-2015/']; eval(cd_to_vels);
[I,R] = readgeoraster('greenland_vel_mosaic250_vx_v1.tif');
V.x = single([R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX]);
V.y = single([R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY]);
V.vx = single(I); V.vx(V.vx==-2.0000e+09) = NaN; %velocity in x-direction in m/yr
clear I R;
[I,~] = readgeoraster('greenland_vel_mosaic250_vy_v1.tif');
V.vy = single(I); V.vy(V.vy==-2.0000e+09) = NaN; %velocity in y-direction in m/yr
clear I;
[VXgrid,VYgrid] = meshgrid(V.x,V.y);

%load the georeferencing info for the ArcticDEM (subset for each glacier in loop because it is HUGE)
ArcticDEMdir = [misc_dir,'ArcticDEM10m/'];
cd_to_arcticDEM = ['cd ',ArcticDEMdir]; eval(cd_to_arcticDEM);
info = geotiffinfo('ArcticDEM_10m_Greenland_Mosaic.tif'); %pull spatial referencing info for the FULL dem
ArcticDEMx = single([info.SpatialRef.XWorldLimits(1)+0.5*info.SpatialRef.CellExtentInWorldX:info.SpatialRef.CellExtentInWorldX:info.SpatialRef.XWorldLimits(2)-0.5*info.SpatialRef.CellExtentInWorldX]);
ArcticDEMy = single([info.SpatialRef.YWorldLimits(2)-0.5*info.SpatialRef.CellExtentInWorldY:-info.SpatialRef.CellExtentInWorldY:info.SpatialRef.YWorldLimits(1)+0.5*info.SpatialRef.CellExtentInWorldY]);

%load the GIMP DEM to help with visualization
cd_to_DEM = ['cd ',misc_dir,'GIMP/']; eval(cd_to_DEM);
[I,R] = readgeoraster('gimpdem_90m_v01.1.tif');
dem.x = single([R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX]);
dem.y = single([R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY]);
dem.z = single(I);
clear I R;
[demx,demy] = meshgrid(dem.x,dem.y);

%loop through the glaciers
for j = 1:length(term)
    
    %pull elevations for the profile from the 10 m-resolution ArcticDEM for the glacier
    cd_to_arcticDEM = ['cd ',ArcticDEMdir]; eval(cd_to_arcticDEM);
    dem_xsub = [find(ArcticDEMx<=min([term(j).RGIX term(j).BoxX]),1,'last'),find(ArcticDEMx>=max([term(j).RGIX term(j).BoxX]),1,'first')];
    dem_ysub = [find(ArcticDEMy>=max([term(j).RGIY term(j).BoxY]),1,'last'),find(ArcticDEMy<=min([term(j).RGIY term(j).BoxY]),1,'first')];
    I = imread('ArcticDEM_10m_Greenland_Mosaic.tif','PixelRegion',{dem_ysub,dem_xsub}); I(I==-9999) = NaN;
    ArcticDEMz = I; clear I;
    [ArcticDEM_xgrid,ArcticDEM_ygrid] = meshgrid(ArcticDEMx(min(dem_xsub):max(dem_xsub)),ArcticDEMy(min(dem_ysub):max(dem_ysub)));
    
    %crop velocities to the RGI outline & overlay as a quiver plot on image
    disp('Prepping velocity vectors to aid mapping of the flux gate at the grounding line, perpendicular to flow...');
    vx_sub = V.vx(find(V.y>=min(im.y) & V.y<=max(im.y)),find(V.x>=min(im.x) & V.x<=max(im.x))); %crop vx to the landsat image extent
    vy_sub = V.vy(find(V.y>=min(im.y) & V.y<=max(im.y)),find(V.x>=min(im.x) & V.x<=max(im.x))); %crop vy to the landsat image extent
    x_sub = V.x(1,find(V.x>=min(im.x) & V.x<=max(im.x)));
    y_sub = V.y(1,find(V.y>=min(im.y) & V.y<=max(im.y)));
    [Xsubgrid,Ysubgrid] = meshgrid(x_sub,y_sub);
    in = inpolygon(Xsubgrid,Ysubgrid,[term(j).RGIX(1:end-1) term(j).RGIX(1)],[term(j).RGIY(1:end-1) term(j).RGIY(1)]); %mask velocities outside the glacier so velocity arrows scale approrpriately
    maskedvx = vx_sub; maskedvx(in==0) = NaN; %mask-out other regions by turning anything not in=1 into NaNs
    maskedvy = vy_sub; maskedvy(in==0) = NaN; %mask-out other regions by turning anything not in=1 into NaNs
    clear in;
    
    %plot the DEM with the velocities overlain
    GL_fig = figure; set(GL_fig,'position',[550 50 700 700]);
    imagesc(double(ArcticDEMx(min(dem_xsub):max(dem_xsub))),double(ArcticDEMy(min(dem_ysub):max(dem_ysub))),ArcticDEMz-nanmedian(geoid)); axis xy equal; colormap gray; hold on;
    set(gca,'clim',[0 term(j).centerGL(3)+20]); cbar = colorbar; cbar.Label.String = 'elevation (m a.s.l.)';
    set(gca,'xlim',[min([term(j).RGIX term(j).BoxX]) max([term(j).RGIX term(j).BoxX])],'ylim',[min([term(j).RGIY term(j).BoxY]) max([term(j).RGIY term(j).BoxY])],'fontsize',14);
    q = quiver(Xsubgrid,Ysubgrid,maskedvx,maskedvy,'-r'); hold on; %q.AutoScaleFactor = 50;
    plot(term(j).RGIX,term(j).RGIY,'-k','linewidth',2); %plot the RGI outline
    plot(term(j).BoxX,term(j).BoxY,'-m','linewidth',2); %plot the terminus box
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
    Zprof(1,:) = interp2(demx,demy,dem.z,vx',vy')-geoidprof';
    Zprof(2,:) = interp2(ArcticDEM_xgrid,ArcticDEM_ygrid,ArcticDEMz,vx',vy')-geoidprof';
    term(j).gateX = vx'; term(j).gateY = vy'; term(j).gateV = vprof; term(j).gateZ = Zprof;
    clear xi yi vx vy vprof Zprof;
    plot(term(j).gateX,term(j).gateY,'+c'); hold on;
    set(gca,'xlim',[min([term(j).RGIX term(j).BoxX]) max([term(j).RGIX term(j).BoxX])],'ylim',[min([term(j).RGIY term(j).BoxY]) max([term(j).RGIY term(j).BoxY])],'fontsize',14);
    xlabel('Easting (m)','fontsize',14); ylabel('Northing (m)','fontsize',14);
    %save the map
    eval(cd_to_centerlines);
    if j<10
        saveas(GL_fig,['00',num2str(j),'_fluxgate-map.png'],'png');
    elseif j >=10 && j<100
        saveas(GL_fig,['0',num2str(j),'_fluxgate-map.png'],'png');
    else
        saveas(GL_fig,[num2str(j),'_fluxgate-map.png'],'png');
    end
    
    %save the data to appropriate places
    eval(cd_to_root);
    save('Greenland_GIC_centerlines.mat','term','-v7.3');
    cd_to_fluxgates = ['cd ',root_dir,'fluxgates/']; eval(cd_to_fluxgates);
    s.Geometry = 'Polyline';
    s.BoundingBox = double([min(term(j).gateX) min(term(j).gateY); max(term(j).gateX) max(term(j).gateY)]);
    s.X = double(term(j).gateX); s.Y = double(term(j).gateY);
    s.RGIId = term(j).RGIref;
    s.BoxID = term(j).BoxID;
    if j<10
        shapewrite(s,['fluxgate_00',num2str(j),'.shp']);
    elseif j >=10 && j<100
        shapewrite(s,['fluxgate_0',num2str(j),'.shp']);
    else
        shapewrite(s,['fluxgate_',num2str(j),'.shp']);
    end
    clear s;
    
    
end


%% extract flux gate speed time series
disp('Matlab REALLY slows down after running through ~100 loops so you may want to force quit & restart');
land_BoxID = [484 485 222 140 142 163 103 608]; %jth index = [428 429 138 47 49 72 6 566]; 
Mankoff_BoxID = [8 101 115 135 156 162 162 174 174 176 238 240 247 248 259 259 259 254 ...
    307 308 310 312 318 322 322 349 349 349 355 356 357 356 357 ...
    363 364 365 366 363 364 365 366 369 424 429 430 431 433 ...
    435 443 447 448 459 460 463 466 467 470 469 471 578 580];

cd_to_vels = ['cd ',misc_dir,'Greenland-ITSLIVE_1985-2018/']; eval(cd_to_vels);
ITSLIVEs = dir('GRE_*.nc'); %velocities in m/yr
disp('Extracting velocities across flux gate...');
counter = 0; disp('Looping through all glaciers...');
for i = 1:length(term)
    disp(['glacier ',num2str(i),' of ',num2str(length(term))]);
    %flag if it is in the Mankoff GrIS timeseries
    if ~isempty(find(Mankoff_BoxID==term(i).BoxID))
        term(i).MankoffFlag = 1;
    else
        term(i).MankoffFlag = 0;
    end
    
    %skip glaciers that have velocities already
    if isempty(term(i).fluxV)
        disp('...creating gate timeseries');
        disp(' ');
        
        %create a densely-spaced (1 m increment) vector for the flux gate
        dl = 1; %spacing increment
        gateX = []; gateY = []; %set-up gate structures
        for j = 1:length(term(i).gateX)-1
            flow_angle = atand((term(i).gateY(j+1)-term(i).gateY(j))./(term(i).gateX(j+1)-term(i).gateX(j)));
            dx = dl.*cosd(flow_angle); dy = dl.*sind(flow_angle);
            if dx == 0
                if sign(term(i).gateY(j+1)-term(i).gateY(j)) < 0
                    gateX = [gateX repmat(term(i).gateX(j),1,length(term(i).gateY(j):-abs(dy):term(i).gateY(j+1)))];
                else
                    gateX = [gateX repmat(term(i).gateX(j),1,length(term(i).gateY(j):abs(dy):term(i).gateY(j+1)))];
                end
            else
                if sign(term(i).gateX(j+1)-term(i).gateX(j)) < 0
                    gateX = [gateX term(i).gateX(j):-abs(dx):term(i).gateX(j+1)];
                else
                    gateX = [gateX term(i).gateX(j):abs(dx):term(i).gateX(j+1)];
                end
            end
            if dy == 0
                if sign(term(i).gateX(j+1)-term(i).gateX(j)) < 0
                    gateY = [gateY repmat(term(i).gateY(j),1,length(term(i).gateX(j):-abs(dx):term(i).gateX(j+1)))];
                else
                    gateY = [gateY repmat(term(i).gateY(j),1,length(term(i).gateX(j):abs(dx):term(i).gateX(j+1)))];
                end
            else
                if sign(term(i).gateY(j+1)-term(i).gateY(j)) < 0
                    gateY = [gateY term(i).gateY(j):-abs(dy):term(i).gateY(j+1)];
                else
                    gateY = [gateY term(i).gateY(j):abs(dy):term(i).gateY(j+1)];
                end
            end
            clear dx dy flow_angle;
        end
        
        %extract velocity cross section time series
        for k=1:length(ITSLIVEs)
            %load the velocity coordinates
            v_x = ncread(ITSLIVEs(k).name,'x');
            v_y = ncread(ITSLIVEs(k).name,'y');
            
            %determine the spatial subset over which to extract data
            vxrefs = find(v_x >= min(gateX) & v_x <= max(gateX));
            if isempty(vxrefs); vxrefs = find(v_x >= min(gateX),1,'first'); end
            vyrefs = find(v_y >= min(gateY) & v_y <= max(gateY));
            if isempty(vyrefs); vyrefs = find(v_y >= min(gateY),1,'last'); end
            if min(vyrefs) >1
                startloc = [min(vxrefs)-1 min(vyrefs)-1]; count = [range(vxrefs)+3 range(vyrefs)+3];
            else
                startloc = [min(vxrefs)-1 min(vyrefs)]; count = [range(vxrefs)+3 range(vyrefs)+2];
            end
            
            %read in the subset
            vx = rot90(ncread(ITSLIVEs(k).name,'vx',startloc,count)); vx(vx==-32767) = NaN;
            vy = rot90(ncread(ITSLIVEs(k).name,'vy',startloc,count)); vy(vy==-32767) = NaN;
            v = sqrt(vx.^2 + vy.^2);
            verr = rot90(ncread(ITSLIVEs(k).name,'v_err',startloc,count));
            vdate = rot90(ncread(ITSLIVEs(k).name,'date',startloc,count));
            vdt = rot90(ncread(ITSLIVEs(k).name,'dt',startloc,count));
            [vyr,vmo,vday,vhh,vmm,vss] =  datevec(double(vdate)); vyr(vyr==0) = NaN; vmo(vmo==0) = NaN; vmo(isnan(vmo)) = 1;
            if mod(nanmean(vyr(~isnan(vyr))),4) ~=0; modays = [31 28 31 30 31 30 31 31 30 31 30 31]; else; modays = [31 29 31 30 31 30 31 31 30 31 30 31]; end
            cumdays = cumsum(modays); cumdays = [0 cumdays(1:11)];
            vdecidate = vyr+(cumdays(vmo)+vday)/sum(modays);
            vdecidt = vdt/sum(modays);
            
            %flip the y axis to agree with rotated grids
            clear v_x; clear v_y;
            v_x = ncread(ITSLIVEs(k).name,'x',startloc(1),count(1));
            v_y = ncread(ITSLIVEs(k).name,'y',startloc(2),count(2));
            v_y = flipud(v_y); %vyrefs = find(v_y >= min(gateY) & v_y <= max(gateY));
            
            %extract data from the flux gate
            [VXgrid,VYgrid] = meshgrid(v_x,v_y);
            for j = 1:length(gateX)
                IDX(j) = find(v_x>=gateX(j),1,'first'); IDY(j) = find(v_y>=gateY(j),1,'first');
            end
            unique_coords = unique([IDX' IDY'],'rows');
            %loop through and pull data from each unique coordinate pair
            for j = 1:size(unique_coords,1)
                ID = find(IDX==unique_coords(j,1) & IDY==unique_coords(j,2));
                if max(ID) < length(gateX)
                    gate_X(j) = nanmean(gateX([ID(1:end) ID(end)+1])); gate_Y(j) = nanmean(gateY([ID(1:end) ID(end)+1]));
                    gate_avgX(j) = nanmean(gate_X); gate_avgY(j) = nanmean(gate_Y);
                    gate_width(j) = sqrt((gateX(ID(end)+1)-gateX(ID(1))).^2 + (gateY(ID(end)+1)-gateY(ID(1))).^2);
                    gate_angle(j) = atand((gateY(ID(end))-gateY(ID(1)))./(gateX(ID(end))-gateX(ID(1)))); %velocity perpendicular to the flux gate
                    gate_vel(k,j) = abs(v(unique_coords(j,2),unique_coords(j,1)).*cosd(gate_angle(j)));
                    gate_velerr(k,j) = verr(unique_coords(j,2),unique_coords(j,1)).*cosd(gate_angle(j));
                    gate_vdate(k,j) = vdecidate(unique_coords(j,2),unique_coords(j,1));
                    
                else
                    gate_X(j) = nanmean(gateX(ID(1:end))); gate_Y(j) = nanmean(gateY(ID(1:end)));
                    gate_width(j) = sqrt((gateX(ID(end))-gateX(ID(1))).^2 + (gateY(ID(end))-gateY(ID(1))).^2);
                    gate_angle(j) = atand((gateY(ID(end))-gateY(ID(1)))./(gateX(ID(end))-gateX(ID(1))));
                    gate_vel(k,j) = abs(v(unique_coords(j,2),unique_coords(j,1)).*cosd(gate_angle(j))); %velocity perpendicular to the flux gate
                    gate_velerr(k,j) = verr(unique_coords(j,2),unique_coords(j,1)).*cosd(gate_angle(j));
                    gate_vdate(k,j) = vdecidate(unique_coords(j,2),unique_coords(j,1));
                    
                end
                clear ID
                
            end
            %         gate_vdateavg(k) = nanmean(gate_vdate(k,:));
            clear v_x v_y v*refs V*grid vx vy v vdir v*date v*dt IDX IDY unique_coords startloc count;
            clear IDX IDY;
        end
        term(i).gateXavg = nanmean(gateX); term(i).gateYavg = nanmean(gateY);
        
        %add to term structure
        term(i).fluxX = gate_X; term(i).fluxY = gate_Y;
        term(i).fluxW = gate_width;
        term(i).fluxAng = gate_angle;
        term(i).fluxV = gate_vel; term(i).fluxVerr = gate_velerr; term(i).fluxVyrs = gate_vdate;
        term(i).fluxVavg = nanmean(gate_vel,2); term(i).gateVdateavg = nanmean(gate_vdate,2);
        clear gate_*;
        %name changes: term(i).fluxVerr = term(i).fluxVerr & term(i).fluxVavg = term(i).fluxVavg
        
        %save after every 50 glaciers
        counter = counter+1;
        if counter == 50
            disp('resaving');
            eval(cd_to_root);
            save('Greenland_GIC_centerlines.mat','term','-v7.3');
            counter = 0;
            eval(cd_to_vels);
        end
        
        %clear the old variables
        clear gate* v* V* 
        
    end
    
end
eval(cd_to_root);
save('Greenland_GIC_centerlines.mat','term','-v7.3');


%% assign to region
%1=W,2=SE,3=E,4=NE,5=N
eval(cd_to_root);
load Greenland_GIC_centerlines.mat;
%loop through & assign region flags
for j=1:length(term)
    if nanmean(term(j).BoxX) < 0
        if nanmean(term(j).BoxY) < -1e6 %west
            term(j).regionFlag = 1;
        else %north
            term(j).regionFlag = 5;
        end
    else
        if nanmean(term(j).BoxY) > -0.855e6 %north
            term(j).regionFlag = 5;
        elseif nanmean(term(j).BoxY) <= -0.855e6 && nanmean(term(j).BoxY) > -1.85e6 %northeast
            term(j).regionFlag = 4;
        elseif nanmean(term(j).BoxY) <= -1.85e6 && nanmean(term(j).BoxY) > -2.5e6 %central east
            term(j).regionFlag = 3;
        else %southeast
            term(j).regionFlag = 2;
        end
    end
end
save('Greenland_GIC_centerlines.mat','term','-v7.3');


%% continue analysis with other files

disp('Use GreenlandGIC_thickness_extrapolation_curves.m to create an empirical scaling relationship to estimate thickness & discharge');
disp('Use GreenlandGIC_discharge_timeseries.m to compute discharge. Combines code from the following:');
disp('interpolate_to_fluxgates.m, extract_termbox_elevation_timeseries.m, extract_racmo.m, deltaHvsdeltaV.m');


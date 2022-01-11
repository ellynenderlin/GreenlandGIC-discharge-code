%extract geometry & speed data for the SE GICs that retreated vs those that
%remained stable
clearvars; close all;

%specify the root directory for project-specific files
root_dir = '/Users/ellynenderlin/Research/NASA_GreenlandPeriph-Mapping/'; %including trailing /
%specify the root directory for generic files
misc_dir = '/Users/ellynenderlin/Research/miscellaneous/'; %including trailing /
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);

%load the Greenland GIC data structure
load Greenland_GIC_centerlines.mat;
for i = 1:length(term)
    BoxID(i) = term(i).BoxID;
end

%load the spatial reference info for the ArcticDEM
ArcticDEMdir = [misc_dir,'ArcticDEM10m/'];
cd_to_arcticDEM = ['cd ',ArcticDEMdir]; eval(cd_to_arcticDEM);
info = geotiffinfo('ArcticDEM_10m_Greenland_Mosaic.tif'); %pull spatial referencing info for the FULL dem
ArcticDEMx = single([info.SpatialRef.XWorldLimits(1)+0.5*info.SpatialRef.CellExtentInWorldX:info.SpatialRef.CellExtentInWorldX:info.SpatialRef.XWorldLimits(2)-0.5*info.SpatialRef.CellExtentInWorldX]);
ArcticDEMy = single([info.SpatialRef.YWorldLimits(2)-0.5*info.SpatialRef.CellExtentInWorldY:-info.SpatialRef.CellExtentInWorldY:info.SpatialRef.YWorldLimits(1)+0.5*info.SpatialRef.CellExtentInWorldY]);
%load the GIMP DEM
cd_to_DEM = ['cd ',misc_dir,'GIMP/']; eval(cd_to_DEM);
[I,R] = geotiffread('gimpdem_90m_v01.1.tif');
dem.x = single([R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX]);
dem.y = single([R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY]);
dem.z = single(I);
clear I R;

%identify retreaters vs stable glaciers
retreat_IDs = [70  71  76  77  87  93  95  99 120 121 122 124 127 129 154 158 159 166,...
 172 178 183 189 191 198 207 220 237 243 244 246 248 253 258 267 268 270,...
 274 275 277 282 294 296 308 311 320 325 327 332 352 355 357 360 369 380,...
 386 402 404 405 410 422 437 438];
stable_IDs = [339 252 108 324 331 383 239 336 323 142 254 419  96 238 411 372 354 321,...
 330 338 329 128 409 240 356 309  94 247 406 358 403 407 329 174 315 345,...
 415 400 249 177 326 353 250 408 421 333 389 173 307 199 433 426 364 375,...
 365 344 168 388 399 362 420 190 444 218 128 319 393 382 368 367 264 312,...
 303 385 390 394 337 441 286 430 245 132 395 391 431 278 373 377 228  90,...
 374 348 314 288];

%extract box width and average speed 2013-2015, average speed 2016, average speed
%2017-2018, & surface slope inland of the 2015 terminus position
for i = 1:length(retreat_IDs)
    termref = find(BoxID==retreat_IDs(i));
    
    %check that the velocity & elevation profiles exist and re-extract if they don't
    if isempty(term(termref).centerZ)
        disp(['Missing profiles for BoxID=',num2str(retreat_IDs(i))]);
        
        %velocity profile time series
        cd_to_vels = ['cd ',misc_dir,'Greenland-ITSLIVE_1985-2018/']; eval(cd_to_vels);
        ITSLIVEs = dir('GRE_*.nc'); speed_cmap = colormap(parula(length(ITSLIVEs)));
        for k=1:length(ITSLIVEs)
            %load the velocity coordinates
            v_x = ncread(ITSLIVEs(k).name,'x');
            v_y = ncread(ITSLIVEs(k).name,'y');
            
            %determine the spatial subset over which to extract data
            vxrefs = find(v_x >= min(term(termref).centerX) & v_x <= max(term(termref).centerX));
            if isempty(vxrefs) && unique(sign(gradient(v_x))) == 1
                vxrefs = [find(v_x >= min(term(termref).centerX),1,'first') find(v_x <= max(term(termref).centerX),1,'last')];
            elseif isempty(vxrefs) && unique(sign(gradient(v_x))) == -1
                vxrefs = [find(v_x >= min(term(termref).centerX),1,'last') find(v_x <= max(term(termref).centerX),1,'first')];
            end
            vyrefs = find(v_y >= min(term(termref).centerY) & v_y <= max(term(termref).centerY));
            if isempty(vyrefs) && unique(sign(gradient(v_y))) == -1
                vyrefs = [find(v_y >= min(term(termref).centerY),1,'last') find(v_y <= max(term(termref).centerY),1,'first')];
            elseif isempty(vyrefs) && unique(sign(gradient(v_y))) == 1
                vyrefs = [find(v_y >= min(term(termref).centerY),1,'first') find(v_y <= max(term(termref).centerY),1,'last')];
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
            v_y = flipud(v_y); vyrefs = find(v_y >= min(term(termref).centerY) & v_y <= max(term(termref).centerY));
            if isempty(vyrefs) && unique(sign(gradient(v_y))) == -1
                vyrefs = [find(v_y >= min(term(termref).centerY),1,'last') find(v_y <= max(term(termref).centerY),1,'first')];
            elseif isempty(vyrefs) && unique(sign(gradient(v_y))) == 1
                vyrefs = [find(v_y >= min(term(termref).centerY),1,'first') find(v_y <= max(term(termref).centerY),1,'last')];
            end
            
            %interpolate the the centerline profile & add to structure
            [vxgrid,vygrid] = meshgrid(v_x(min(vxrefs)-1:max(vxrefs)+1),v_y(min(vyrefs)-1:max(vyrefs)+1));
            term(termref).centerVdate(k,:) = interp2(vxgrid,vygrid,vdecidate,term(termref).centerX,term(termref).centerY,'nearest');
            term(termref).centerVdt(k,:) = interp2(vxgrid,vygrid,vdecidt,term(termref).centerX,term(termref).centerY,'nearest');
            term(termref).centerV(k,:) = interp2(vxgrid,vygrid,sqrt(vx.^2 + vy.^2),term(termref).centerX,term(termref).centerY);
            term(termref).centerVdir(k,:) = interp2(vxgrid,vygrid,atand(vy./vx),term(termref).centerX,term(termref).centerY);
            term(termref).centerVerr(k,:) = interp2(vxgrid,vygrid,verr,term(termref).centerX,term(termref).centerY);
            
            clear v_x v_y vx* vy* v*date v*dt;
        end
        disp('speed profiles extracted');
        
        %pull elevations for the profile from the 10 m-resolution ArcticDEM for the glacier
        cd_to_arcticDEM = ['cd ',ArcticDEMdir]; eval(cd_to_arcticDEM);
        dem_xsub = [find(ArcticDEMx<=min([term(termref).RGIX term(termref).BoxX]),1,'last'),find(ArcticDEMx>=max([term(termref).RGIX term(termref).BoxX]),1,'first')];
        dem_ysub = [find(ArcticDEMy>=max([term(termref).RGIY term(termref).BoxY]),1,'last'),find(ArcticDEMy<=min([term(termref).RGIY term(termref).BoxY]),1,'first')];
        I = imread('ArcticDEM_10m_Greenland_Mosaic.tif','PixelRegion',{dem_ysub,dem_xsub}); I(I==-9999) = NaN;
        ArcticDEMz = I; clear I;
        [ArcticDEM_xgrid,ArcticDEM_ygrid] = meshgrid(ArcticDEMx(min(dem_xsub):max(dem_xsub)),ArcticDEMy(min(dem_ysub):max(dem_ysub)));
        term(termref).centerZdate(1,:) = 2008; term(termref).centerZdate(2,:) = 2018;
        term(termref).centerZ(1,:) = interp2(double(dem.x),double(dem.y),dem.z,double(term(termref).centerX),double(term(termref).centerY)); %GIMP elevation profile (from ~2008)
        term(termref).centerZ(2,:) = interp2(double(ArcticDEMx(min(dem_xsub):max(dem_xsub))),double(ArcticDEMy(min(dem_ysub):max(dem_ysub))),ArcticDEMz,double(term(termref).centerX),double(term(termref).centerY)); %ArcticDEM profile (from ~2018)
        disp('elevation profiles extracted');
        
    end
    
    %save info to vectors
    termdate_retreat(i) = term(termref).year(3)+term(termref).doy(3)/365;
    width_retreat(i) = sqrt((term(termref).inlandX(2) - term(termref).inlandX(1)).^2 + (term(termref).inlandY(2) - term(termref).inlandY(1)).^2);
    slope05_retreat(i) = (term(termref).centerZ(2,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+500,1,'last')) - term(termref).centerZ(2,term(termref).center_termref(3)))./(term(termref).centerline(1,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+500,1,'last')) - term(termref).centerline(1,term(termref).center_termref(3)));
    disp(['Slope over 500 m = ',num2str(slope05_retreat(i)),'m']);
    disp(['elevation change = ',num2str((term(termref).centerZ(2,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+500,1,'last')) - term(termref).centerZ(2,term(termref).center_termref(3)))),'m']);
    disp(['distance = ',num2str(term(termref).centerline(1,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+500,1,'last')) - term(termref).centerline(1,term(termref).center_termref(3))),'m']);
    slope1_retreat(i) = (term(termref).centerZ(2,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+1000,1,'last')) - term(termref).centerZ(2,term(termref).center_termref(3)))./(term(termref).centerline(1,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+1000,1,'last')) - term(termref).centerline(1,term(termref).center_termref(3)));
    slope2_retreat(i) = (term(termref).centerZ(2,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+2000,1,'last')) - term(termref).centerZ(2,term(termref).center_termref(3)))./(term(termref).centerline(1,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+2000,1,'last')) - term(termref).centerline(1,term(termref).center_termref(3)));
    speedearly_retreat(i) = nanmean(term(termref).centerV(29:31,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+1000,1,'last')));
    speed2016_retreat(i) = term(termref).centerV(32,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+1000,1,'last'));
    speedlate_retreat(i) = nanmean(term(termref).centerV(33:34,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+1000,1,'last')));
    
    %move on
    cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
    clear termref;
end
for i = 1:length(stable_IDs)
    termref = find(BoxID==stable_IDs(i));
    
    %check that the velocity & elevation profiles exist and re-extract if they don't
    if isempty(term(termref).centerZ)
        disp(['Missing profiles for BoxID=',num2str(retreat_IDs(i))]);
        
        %velocity profile time series
        cd_to_vels = ['cd ',misc_dir,'Greenland-ITSLIVE_1985-2018/']; eval(cd_to_vels);
        ITSLIVEs = dir('GRE_*.nc'); speed_cmap = colormap(parula(length(ITSLIVEs)));
        for k=1:length(ITSLIVEs)
            %load the velocity coordinates
            v_x = ncread(ITSLIVEs(k).name,'x');
            v_y = ncread(ITSLIVEs(k).name,'y');
            
            %determine the spatial subset over which to extract data
            vxrefs = find(v_x >= min(term(termref).centerX) & v_x <= max(term(termref).centerX));
            if isempty(vxrefs) && unique(sign(gradient(v_x))) == 1
                vxrefs = [find(v_x >= min(term(termref).centerX),1,'first') find(v_x <= max(term(termref).centerX),1,'last')];
            elseif isempty(vxrefs) && unique(sign(gradient(v_x))) == -1
                vxrefs = [find(v_x >= min(term(termref).centerX),1,'last') find(v_x <= max(term(termref).centerX),1,'first')];
            end
            vyrefs = find(v_y >= min(term(termref).centerY) & v_y <= max(term(termref).centerY));
            if isempty(vyrefs) && unique(sign(gradient(v_y))) == -1
                vyrefs = [find(v_y >= min(term(termref).centerY),1,'last') find(v_y <= max(term(termref).centerY),1,'first')];
            elseif isempty(vyrefs) && unique(sign(gradient(v_y))) == 1
                vyrefs = [find(v_y >= min(term(termref).centerY),1,'first') find(v_y <= max(term(termref).centerY),1,'last')];
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
            v_y = flipud(v_y); vyrefs = find(v_y >= min(term(termref).centerY) & v_y <= max(term(termref).centerY));
            if isempty(vyrefs) && unique(sign(gradient(v_y))) == -1
                vyrefs = [find(v_y >= min(term(termref).centerY),1,'last') find(v_y <= max(term(termref).centerY),1,'first')];
            elseif isempty(vyrefs) && unique(sign(gradient(v_y))) == 1
                vyrefs = [find(v_y >= min(term(termref).centerY),1,'first') find(v_y <= max(term(termref).centerY),1,'last')];
            end
            
            %interpolate the the centerline profile & add to structure
            [vxgrid,vygrid] = meshgrid(v_x(min(vxrefs)-1:max(vxrefs)+1),v_y(min(vyrefs)-1:max(vyrefs)+1));
            term(termref).centerVdate(k,:) = interp2(vxgrid,vygrid,vdecidate,term(termref).centerX,term(termref).centerY,'nearest');
            term(termref).centerVdt(k,:) = interp2(vxgrid,vygrid,vdecidt,term(termref).centerX,term(termref).centerY,'nearest');
            term(termref).centerV(k,:) = interp2(vxgrid,vygrid,sqrt(vx.^2 + vy.^2),term(termref).centerX,term(termref).centerY);
            term(termref).centerVdir(k,:) = interp2(vxgrid,vygrid,atand(vy./vx),term(termref).centerX,term(termref).centerY);
            term(termref).centerVerr(k,:) = interp2(vxgrid,vygrid,verr,term(termref).centerX,term(termref).centerY);
            
            clear v_x v_y vx* vy* v*date v*dt;
        end
        disp('speed profiles extracted');
        
        %pull elevations for the profile from the 10 m-resolution ArcticDEM for the glacier
        cd_to_arcticDEM = ['cd ',ArcticDEMdir]; eval(cd_to_arcticDEM);
        dem_xsub = [find(ArcticDEMx<=min([term(termref).RGIX term(termref).BoxX]),1,'last'),find(ArcticDEMx>=max([term(termref).RGIX term(termref).BoxX]),1,'first')];
        dem_ysub = [find(ArcticDEMy>=max([term(termref).RGIY term(termref).BoxY]),1,'last'),find(ArcticDEMy<=min([term(termref).RGIY term(termref).BoxY]),1,'first')];
        I = imread('ArcticDEM_10m_Greenland_Mosaic.tif','PixelRegion',{dem_ysub,dem_xsub}); I(I==-9999) = NaN;
        ArcticDEMz = I; clear I;
        [ArcticDEM_xgrid,ArcticDEM_ygrid] = meshgrid(ArcticDEMx(min(dem_xsub):max(dem_xsub)),ArcticDEMy(min(dem_ysub):max(dem_ysub)));
        term(termref).centerZdate(1,:) = 2008; term(termref).centerZdate(2,:) = 2018;
        term(termref).centerZ(1,:) = interp2(double(dem.x),double(dem.y),dem.z,double(term(termref).centerX),double(term(termref).centerY)); %GIMP elevation profile (from ~2008)
        term(termref).centerZ(2,:) = interp2(double(ArcticDEMx(min(dem_xsub):max(dem_xsub))),double(ArcticDEMy(min(dem_ysub):max(dem_ysub))),ArcticDEMz,double(term(termref).centerX),double(term(termref).centerY)); %ArcticDEM profile (from ~2018)
        disp('elevation profiles extracted');
        
    end
    
    %save info to vectors
    termdate_stable(i) = term(termref).year(3)+term(termref).doy(3)/365;
    width_stable(i) = sqrt((term(termref).inlandX(2) - term(termref).inlandX(1)).^2 + (term(termref).inlandY(2) - term(termref).inlandY(1)).^2);
    slope05_stable(i) = (term(termref).centerZ(2,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+500,1,'last')) - term(termref).centerZ(2,term(termref).center_termref(3)))./(term(termref).centerline(1,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+500,1,'last')) - term(termref).centerline(1,term(termref).center_termref(3)));
    disp(['Slope over 500 m = ',num2str(slope05_stable(i)),'m']);
    disp(['elevation change = ',num2str((term(termref).centerZ(2,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+500,1,'last')) - term(termref).centerZ(2,term(termref).center_termref(3)))),'m']);
    disp(['distance = ',num2str(term(termref).centerline(1,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+500,1,'last')) - term(termref).centerline(1,term(termref).center_termref(3))),'m']);
    slope1_stable(i) = (term(termref).centerZ(2,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+1000,1,'last')) - term(termref).centerZ(2,term(termref).center_termref(3)))./(term(termref).centerline(1,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+1000,1,'last')) - term(termref).centerline(1,term(termref).center_termref(3)));
    slope2_stable(i) = (term(termref).centerZ(2,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+2000,1,'last')) - term(termref).centerZ(2,term(termref).center_termref(3)))./(term(termref).centerline(1,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+2000,1,'last')) - term(termref).centerline(1,term(termref).center_termref(3)));
    speedearly_stable(i) = nanmean(term(termref).centerV(29:31,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+1000,1,'last')));
    speed2016_stable(i) = term(termref).centerV(32,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+1000,1,'last'));
    speedlate_stable(i) = nanmean(term(termref).centerV(33:34,find(term(termref).centerline <= term(termref).centerline(term(termref).center_termref(3))+1000,1,'last')));
    
    %move on
    cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
    clear termref;
end
%resave the data structure if necessary
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
% save('Greenland_GIC_centerlines.mat','term','-v7.3');

%plot histograms to spot check if there are any obvious differences
figure; set(gcf,'position',[50 50 800 1200]);
subplot(3,1,1);
h(1) = histogram(width_stable,[0:500:4500],'FaceColor',[0 0 1],'EdgeColor','k'); hold on;
h(2) = histogram(width_retreat,[0:500:4500],'FaceColor',[1 0.5 0],'EdgeColor','k'); hold on;
subplot(3,1,2);
histogram(slope05_stable,[0:0.05:0.6],'FaceColor',[0 0 1],'EdgeColor','k'); hold on;
histogram(slope1_stable,[0:0.05:0.6],'FaceColor',[0 0 1],'EdgeColor','w'); hold on;
histogram(slope05_retreat,[0:0.05:0.6],'FaceColor',[1 0.5 0],'EdgeColor','k'); hold on;
histogram(slope1_retreat,[0:0.05:0.6],'FaceColor',[1 0.5 0],'EdgeColor','w'); hold on;
subplot(3,1,3);
histogram(speedearly_stable,[0:250:3000],'FaceColor',[0 0 1],'EdgeColor','k'); hold on;
histogram(speed2016_stable,[0:250:3000],'FaceColor',[0 1 1],'EdgeColor','k'); hold on;
histogram(speedlate_stable,[0:250:3000],'FaceColor',[0 1 0],'EdgeColor','k'); hold on;
histogram(speedearly_retreat,[0:250:3000],'FaceColor',[1 0.5 0],'EdgeColor','k'); hold on;
histogram(speed2016_retreat,[0:250:3000],'FaceColor',[1 0 0],'EdgeColor','k'); hold on;
histogram(speedlate_retreat,[0:250:3000],'FaceColor',[1 1 0],'EdgeColor','k'); hold on;

%hone in on the slope differences
figure; set(gcf,'position',[450 50 800 800]);
subplot(1,2,1);
histogram(slope05_stable,[0:0.05:0.6],'FaceColor',[0.5 0 1],'EdgeColor','k'); hold on;
histogram(slope1_stable,[0:0.05:0.6],'FaceColor',[0 0 1],'EdgeColor','w'); hold on;
subplot(1,2,2);
histogram(slope05_retreat,[0:0.05:0.6],'FaceColor',[1 0 0],'EdgeColor','k'); hold on;
histogram(slope1_retreat,[0:0.05:0.6],'FaceColor',[1 0.5 0],'EdgeColor','w'); hold on;


%export as csvs
retreat_info = [retreat_IDs' width_retreat' termdate_retreat' slope05_retreat' slope1_retreat' speedearly_retreat' speed2016_retreat' speedlate_retreat'];
csvwrite('SEGreenland-GIC_2016retreater-geoms.csv',retreat_info);
stable_info = [stable_IDs' width_stable' termdate_stable' slope05_stable' slope1_stable' speedearly_stable' speed2016_stable' speedlate_stable'];
csvwrite('SEGreenland-GIC_2016stable-geoms.csv',stable_info);

%% extract hypsometric curves for retreated and stable glaciers
close all;
retreat_hyp = figure;
for i = 1:length(retreat_IDs)
%     if i > 1;
%     close hfig; drawnow;
%     end
    termref = find(BoxID==retreat_IDs(i));
    
    %pull elevations for the profile from the 10 m-resolution ArcticDEM for the glacier
    cd_to_arcticDEM = ['cd ',ArcticDEMdir]; eval(cd_to_arcticDEM);
    dem_xsub = [find(ArcticDEMx<=min([term(termref).RGIX term(termref).BoxX]),1,'last'),find(ArcticDEMx>=max([term(termref).RGIX term(termref).BoxX]),1,'first')];
    dem_ysub = [find(ArcticDEMy>=max([term(termref).RGIY term(termref).BoxY]),1,'last'),find(ArcticDEMy<=min([term(termref).RGIY term(termref).BoxY]),1,'first')];
    I = imread('ArcticDEM_10m_Greenland_Mosaic.tif','PixelRegion',{dem_ysub,dem_xsub}); I(I==-9999) = NaN;
    ArcticDEMz = I; clear I;
    [ArcticDEM_xgrid,ArcticDEM_ygrid] = meshgrid(ArcticDEMx(min(dem_xsub):max(dem_xsub)),ArcticDEMy(min(dem_ysub):max(dem_ysub)));
    
    %use the RGI outline to mask the DEM
    RGIpolyrefs = find(isnan(term(termref).RGIX)==1);
    elevs(i).retreat = []; elevs(i).retreatmap = NaN(size(ArcticDEMz));
    for j = 1:length(RGIpolyrefs)
       if j == 1
           [in,~] = inpolygon(ArcticDEM_xgrid,ArcticDEM_ygrid,term(termref).RGIX(1:RGIpolyrefs(j)),term(termref).RGIY(1:RGIpolyrefs(j)));
           elevs(i).retreat = ArcticDEMz(in);
           elevs(i).retreatmap(in==1) = ArcticDEMz(in==1);
           clear in;
       else
           [in,~] = inpolygon(ArcticDEM_xgrid,ArcticDEM_ygrid,term(termref).RGIX(RGIpolyrefs(j-1):RGIpolyrefs(j)),term(termref).RGIY(RGIpolyrefs(j-1):RGIpolyrefs(j)));
           elevs(i).retreat = [elevs(i).retreat; ArcticDEMz(in)];
           elevs(i).retreatmap(in==1) = ArcticDEMz(in==1);
           clear in;
       end
    end
    
    %create output hypsometric curves
    hfig = figure;
    h = histogram(elevs(i).retreatmap,'Normalization','cdf','Orientation','horizontal');
    elevs(i).retreat_histZ = (h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
    elevs(i).retreat_histCDF = h.Values;
    figure(retreat_hyp); plot(elevs(i).retreat_histCDF./max(elevs(i).retreat_histCDF),elevs(i).retreat_histZ,'-','linewidth',2); hold on;
    clear h dem_*sub ArcticDEMz ArcticDEM_*grid;
end







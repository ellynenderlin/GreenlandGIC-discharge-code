%ESTIMATE FLUX ACROSS GATE
clearvars; close all; drawnow;

%load file containing flux gate data
% cd /Users/katebollen/Documents/MS/data/DEMs/centerlines/
cd /users/ellynenderlin/Research/NASA_Greenland-Periph-Mapping/centerlines/
load Greenland_GIC_centerlines.mat;
figure; set(gcf,'position',[50 50 1200 800]); cmap = colormap(hsv(length(term))); drawnow;
sub1 = subplot(1,3,1); sub2 = subplot(1,3,2); sub3 = subplot(1,3,3);
% cd ../../HxUxW
cd ..
% load H_vs_UdivW_functions.mat; 
load D_vs_UW_functions.mat;
ci = confint(f.fit); 

%loop through the glaciers & extract speeds from flux gates
% cd ../velocities/ITS_LIVE/NetCDF/
cd /users/ellynenderlin/Research/miscellaneous/Greenland-ITSLIVE_1985-2018/
disp('Extracting velocities across flux gate...');
Dtotal = [];
Dtotal_err = [];
for i = 1:length(term)
%     disp(['... glacier ',num2str(i),' of ',num2str(length(term))]);
    if ~isempty(term(i).gateX)
%     disp(['gate drawn for BoxID ',term(i).BoxID]);
    disp(['... glacier ',num2str(i),' of ',num2str(length(term))]);
    
    %densely interpolate the flux gate
    dl = 1; %spacing increment
    gateX = []; gateY = []; %set-up gate structures
    for j = 1%:length(term(i).gateX)-1
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
    
    
    %load the velocity data
    ITSLIVEs = dir('GRE_*.nc'); %velocities in m/yr
    for k=1:length(ITSLIVEs)
        %load the velocity coordinates
        v_x = ncread(ITSLIVEs(k).name,'x');
        v_y = ncread(ITSLIVEs(k).name,'y');
        
        %determine the spatial subset over which to extract data
        vxrefs = find(v_x >= min(gateX) & v_x <= max(gateX));
        if isempty(vxrefs); vxrefs = find(v_x >= min(gateX),1,'first'); end
        vyrefs = find(v_y >= min(gateY) & v_y <= max(gateY));
        if isempty(vyrefs); vyrefs = find(v_y >= min(gateY),1,'first'); end
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
    
    
    %estimate glacier thickness using gate_vel./gate_width
    gate_D = NaN(size(gate_vel)); gate_Derr = NaN(size(gate_vel)); 
    glacier_width = sqrt((term(i).gateX(1) - term(i).gateX(end)).^2 + (term(i).gateY(1)-term(i).gateY(end)).^2);
    gate_D(gate_vel*glacier_width<=max(f.xlims)) = f.fit.p1.*((gate_vel(gate_vel*glacier_width<=max(f.xlims))*glacier_width).^2) + f.fit.p2.*(gate_vel(gate_vel*glacier_width<=max(f.xlims))*glacier_width) + f.fit.p3;
    gate_Derr(gate_vel*glacier_width<=max(f.xlims)) = abs((ci(1,1).*((gate_vel(gate_vel*glacier_width<=max(f.xlims))*glacier_width).^2) + ci(1,2).*(gate_vel(gate_vel*glacier_width<=max(f.xlims))*glacier_width))-gate_D(gate_vel*glacier_width<=max(f.xlims)));
    gate_D(gate_vel*glacier_width>max(f.xlims)) = f.ext.p1.*gate_vel(gate_vel*glacier_width>max(f.xlims))*glacier_width + f.ext.p2;
    gate_Derr(gate_vel*glacier_width>max(f.xlims)) = abs((f.ext.ciu_p1.*gate_vel(gate_vel*glacier_width>max(f.xlims))*glacier_width + f.ext.ciu_p2)-gate_D(gate_vel*glacier_width>max(f.xlims)));
    gate_mean_vel = nanmean(gate_vel,2);
    gate_D(gate_D<=0.1e8) = 0.1e8;
    
    %estimate discharge for bins & sum to get glacier discharge
    for k = 1:size(gate_vel,1)
        D_yr(k) = 1985+k-1;
    end
    V_yr = nanmean(gate_vdate,2);
    D = nansum(gate_D,2); Derr = nansum(gate_Derr,2);
%     if nanmedian(D) < 0; disp(num2str(i)); end

    %plot
    subplot(sub1);
    plot(V_yr(~isnan(D)),D(~isnan(D)).*917./(1000*10^9),'-','linewidth',1.5,'color',cmap(i,:)); hold on; grid on; drawnow;
    title ('Annual Discharge','fontSize', 24); set(gca,'xlim',[1985 2018]);
    xlabel('Year'); ylabel('Discharge (Gt/yr)','fontSize',24); set(gca,'FontSize',20)
    subplot(sub2);
    plot(nanmedian(max(gate_vel(~isnan(D),:),[],2,'includenan')),nanmedian(D).*917./(1000*10^9),'x','linewidth',1.5,'color',cmap(i,:)); hold on; grid on; drawnow;
    title ('Discharge vs Velocity','fontSize',24);
    xlabel('Velocity (m/yr)'); ylabel('Median discharge (Gt/yr)','fontSize', 24); set(gca,'FontSize',20)
    subplot(sub3);
    plot(sqrt((term(i).gateX(1) - term(i).gateX(end)).^2 + (term(i).gateY(1)-term(i).gateY(end)).^2),nanmedian(D).*917./(1000*10^9),'x','linewidth',1.5,'color',cmap(i,:)); hold on; grid on; drawnow;
    title ('Discharge vs Fluxgate width (m)','fontSize',24);
    xlabel('Gate width (m)'); ylabel('Median discharge (Gt/yr)','fontSize', 24); set(gca,'FontSize',20)

    
    %add to term structure
    term(i).fluxX = gate_X; term(i).fluxY = gate_Y; 
    term(i).fluxW = gate_width; 
    term(i).fluxAng = gate_angle; 
    term(i).fluxV = gate_vel; term(i).Verr = gate_velerr; term(i).fluxVyrs = gate_vdate;
%     term(i).fluxH = gate_H; term(i).Herr = gate_Herr;
    term(i).fluxD = D; term(i).fluxDerr = Derr; term(i).fluxDyrs = D_yr';
    term(i).fluxVelavg = gate_mean_vel;
    term(i).gateavgX = nanmean(gateX); term(i).gateavgY = nanmean(gateY);
    term(i).gateVdateavg = V_yr;
    
    %add to the entire regional D
    Dtotal = [Dtotal; D'];
    Dtotal_err = [Dtotal_err; Derr'];
    
    %clear the old variables
    clear D D_yr gate* v* V* Derr Dbinned;
    end
end
set(sub1,'xlim',[1985 2018]);
% cd /Users/katebollen/Documents/MS/data/DEMs/centerlines/
cd /users/ellynenderlin/Research/NASA_Greenland-Periph-Mapping/centerlines/
save('Greenland_GIC_centerlines.mat','term','-v7.3');
cd ..
saveas(gcf,'Greenland-GIC-glacier-discharge.png','png'); saveas(gcf,'Greenland-GIC-glacier-discharge.eps','epsc');

%time series of total discharge
figure; set(gcf,'position',[1300 50 400 800]);
suba = subplot(2,1,1);
plot(term(1).fluxDyrs,nansum(Dtotal).*917./(1000*10^9),'.-k','linewidth',2); hold on;
errorbar(term(1).fluxDyrs,(nansum(Dtotal).*917./(1000*10^9)),sqrt(nansum(Dtotal_err.^2)).*917./(1000*10^9),sqrt(nansum(Dtotal_err.^2)).*917./(1000*10^9));
title ('Annual Greenland GIC Discharge','fontSize', 24);
xlabel('Year'); ylabel('Discharge (Gt/yr)','fontSize',24); set(gca,'FontSize',20,'xlim',[1985 2018]); grid on;
subb = subplot(2,1,2);
plot(term(1).fluxDyrs,100*nansum(Dtotal~=0)./size(Dtotal,1),'.-k','linewidth',2); hold on;
title ('Discharge Data Coverage','fontSize', 24);
xlabel('Year'); ylabel('Percent glaciers','fontSize',24); set(gca,'FontSize',20,'xlim',[1985 2018]); grid on;
saveas(gcf,'Greenland-GIC-discharge.png','png'); saveas(gcf,'Greenland-GIC-discharge.eps','epsc');

%plot a map of locations where we have discharge estimates
cd /Users/ellynenderlin/Research/miscellaneous/Greenland-ITSLIVE_1985-2018/
[A,R] = geotiffread('Greenland_ITSLIVEavg_1985-2018.tif');
vx = R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX;
vy = R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY;
v = A; clear A R;
figure; set(gcf,'position',[800 50 800 1200]);
imagesc(vx,vy,v); axis xy equal; hold on;
vel_cmap = colormap(parula(10001)); vel_cmap(1,:) = [1 1 1]; colormap(gca,vel_cmap);
cbar = colorbar; set(gca,'clim',[0 500]);
for i = 1:length(term)
    if ~isempty(term(i).gateX)
        plot(nanmean(term(i).X),nanmean(term(i).Y),'ok','markerfacecolor','m','markersize',10); hold on;
    else
        plot(nanmean(term(i).X),nanmean(term(i).Y),'ok','markerfacecolor','w','markersize',10); hold on;
    end
end

%EXTRACT SPEEDS PERPENDICULAR TO FLUX GATES
clear all; close all; drawnow;

%load file containing flux gate data
cd /Users/ellynenderlin/Research/NASA-Greenland-Periph-Mapping/centerlines/
load Greenland_GIC_centerlines.mat;
figure; set(gcf,'position',[50 50 1200 800]); cmap = colormap(hsv(length(term))); drawnow;
sub1 = subplot(1,3,1); sub2 = subplot(1,3,2); sub3 = subplot(1,3,3);
cd ..
load H_vs_UdivW_functions.mat; ci = confint(f.fit); 

%loop through the glaciers & extract speeds from flux gates
cd ../miscellaneous/Greenland-ITSLIVE_1985-2018/
disp('Extracting velocities across flux gate...');
for i = 1:length(term)
%     disp(['... glacier ',num2str(i),' of ',num2str(length(term))]);
    if ~isempty(term(i).gateX)
%     disp(['gate drawn for BoxID ',term(i).BoxID]);
    disp(['... glacier ',num2str(i),' of ',num2str(length(term))]);
    
    %densely interpolate the flux gate
    dl = 1; %spacing increment
    gateX = []; gateY = []; %set-up gate structures
    for j = 1:length(term(i).gateX)-1
        flow_angle = atand((term(i).gateY(j+1)-term(i).gateY(j))./(term(i).gateX(j+1)-term(i).gateX(j)));
        dx = dl.*cosd(flow_angle); dy = dl.*sind(flow_angle); 
        if dx == 0
            if sign(term(i).gateY(j+1)-term(i).gateY(j)) < 0
                gateX = [gateX repmat(gateX(end),1,length(term(i).gateY(j):-abs(dy):term(i).gateY(j+1)))];
            else
                gateX = [gateX repmat(gateX(end),1,length(term(i).gateY(j):abs(dy):term(i).gateY(j+1)))];
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
                gateY = [gateY repmat(gateY(end),1,length(term(i).gateX(j):-abs(dx):term(i).gateX(j+1)))];
            else
                gateY = [gateY repmat(gateY(end),1,length(term(i).gateX(j):abs(dx):term(i).gateX(j+1)))];
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
        vyrefs = find(v_y >= min(gateY) & v_y <= max(gateY));
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
        v_y = flipud(v_y); vyrefs = find(v_y >= min(gateY) & v_y <= max(gateY));
        
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
                gate_width(j) = sqrt((gateX(ID(end)+1)-gateX(ID(1))).^2 + (gateY(ID(end)+1)-gateY(ID(1))).^2);
                gate_angle(j) = atand((gateY(ID(end))-gateY(ID(1)))./(gateX(ID(end))-gateX(ID(1)))); %velocity perpendicular to the flux gate
                gate_vel(k,j) = v(unique_coords(j,2),unique_coords(j,1)).*cosd(gate_angle(j));
                gate_velerr(k,j) = verr(unique_coords(j,2),unique_coords(j,1)).*cosd(gate_angle(j));
            else
                gate_X(j) = nanmean(gateX(ID(1:end))); gate_Y(j) = nanmean(gateY(ID(1:end)));
                gate_width(j) = sqrt((gateX(ID(end))-gateX(ID(1))).^2 + (gateY(ID(end))-gateY(ID(1))).^2);
                gate_angle(j) = atand((gateY(ID(end))-gateY(ID(1)))./(gateX(ID(end))-gateX(ID(1))));
                gate_vel(k,j) = v(unique_coords(j,2),unique_coords(j,1)).*cosd(gate_angle(j)); %velocity perpendicular to the flux gate
                gate_velerr(k,j) = verr(unique_coords(j,2),unique_coords(j,1)).*cosd(gate_angle(j));
            end
            clear ID
            
        end
        clear v_x v_y v*refs V*grid vx vy v vdir v*date v*dt IDX IDY unique_coords startloc count;
        clear IDX IDY;
    end
    
    %estimate glacier thickness using gate_vel./gate_width
    gate_H = NaN(size(gate_vel)); gate_Herr = NaN(size(gate_vel)); 
    glacier_width = sqrt((term(i).gateX(1) - term(i).gateX(end)).^2 + (term(i).gateY(1)-term(i).gateY(end)).^2);
    UdivW = nanmean((gate_vel(1:10,:)./glacier_width),1);
    gate_H(UdivW<=max(f.xlims)) = f.fit.p1.*UdivW(UdivW<=max(f.xlims)).^2 + f.fit.p2.*UdivW(UdivW<=max(f.xlims)) + f.fit.p3;
    gate_Herr(UdivW<=max(f.xlims)) = abs((ci(1,1).*UdivW(UdivW<=max(f.xlims)).^2 + ci(1,2).*UdivW(UdivW<=max(f.xlims)))-gate_H(UdivW<=max(f.xlims)));
    gate_H(UdivW>max(f.xlims)) = f.ext.p1.*UdivW(UdivW>max(f.xlims)) + f.ext.p2;
    gate_Herr(UdivW>max(f.xlims)) = abs((f.ext.ciu_p1.*UdivW(UdivW>max(f.xlims)) + f.ext.ciu_p2)-gate_H(UdivW>max(f.xlims)));
    if ~isempty(gate_H(gate_H<0))
        disp('Negative thicknesses: ');
        disp(['W = ',num2str(gate_width)]);
        disp(['U = ',num2str(gate_vel(gate_H<0)')]);
        disp(['H = ',num2str(gate_H(gate_H<0)')])
        gate_H(gate_H<0) = 0;
    end
    
    %estimate discharge for bins & sum to get glacier discharge
    for k = 1:size(gate_H,1)
        if sum(isnan(gate_H(k,:))) == 0
            D(k) = sum(gate_H(k,:).*gate_vel(k,:).*gate_width); %m^3/yr
            if D(k) > 10^9; disp(['BoxID = ',term(i).BoxID,' & # = ',num2str(i)]); end
        else
            D(k) = NaN;
        end
        D_yr(k) = 1985+k-1;
    end
    %filter outliers (likely bad velocities)
%     Dmed = nanmedian(D); Dmad = mad(D(~isnan(D)),1);
%     D(D > Dmed + 2*1.4826*Dmed) = NaN; D(D < Dmed - 2*1.4826*Dmed) = NaN; 
%     clear Dmed Dmad;
    %plot
    subplot(sub1);
    plot(D_yr,D,'x','linewidth',1.5,'color',cmap(i,:)); hold on; grid on; drawnow; 
    subplot(sub2);
    plot(nanmean(max(gate_vel,[],2,'includenan')),nanmean(D),'x','linewidth',1.5,'color',cmap(i,:)); hold on; grid on; drawnow; 
    subplot(sub3);
    plot(sqrt((term(i).gateX(1) - term(i).gateX(end)).^2 + (term(i).gateY(1)-term(i).gateY(end)).^2),nanmean(D),'x','linewidth',1.5,'color',cmap(i,:)); hold on; grid on; drawnow; 

    
    %add to term structure
    term(i).fluxX = gate_X; term(i).fluxY = gate_Y; 
    term(i).fluxW = gate_width; 
    term(i).fluxAng = gate_angle; 
    term(i).fluxV = gate_vel; 
    term(i).fluxH = gate_H; 
    term(i).fluxD = D; term(i).fluxDyrs = D_yr; 
    
    clear D D_yr gate* v* V*;
    end
end
cd /Users/ellynenderlin/Research/NASA-Greenland-Periph-Mapping/centerlines/
save('Greenland_GIC_centerlines.mat','term','-v7.3');


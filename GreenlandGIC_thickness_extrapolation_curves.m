%create empirical function to estimate thickness
clearvars; close all;
addpath('/users/ellynenderlin/mfiles/general/','/users/ellynenderlin/mfiles/general/polyfitZero-1.3');

%load the data
root_dir = '/Users/ellynenderlin/Research/NASA_GreenlandPeriph-Mapping/';
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
load Greenland_GIC_centerlines.mat;
avg_thick = csvread('avg_thick_vel_2.csv');
cd ../miscellaneous/Greenland-ITSLIVE_1985-2018/
[A,R] = readgeoraster('Greenland_ITSLIVEavg_1985-2018.tif');
V.x = R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX;
V.y = R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY;
V.z = A; clear A R;
[vxgrid,vygrid] = meshgrid(V.x,V.y);
for i = 1:length(avg_thick(:,6))
    speed(i) = interp2(vxgrid,vygrid,V.z,avg_thick(i,1),avg_thick(i,2),'nearest');
end
avg_thick(:,6) = speed';
eval(cd_to_root);
%column key: 3=BoxID 4=thickness (m) 6=velocity (m/yr) 8=width (m) 9=u/w 10=uw 11=huw

%specify the colormap (5 regions)
region_cmap = [253,174,97; 215,25,28; 255,255,191; 171,217,233; 44,123,182]./255;

%load Ken Mankoff's data for the glaciers IDed as overlapping
cd /users/ellynenderlin/Research/NASA_GreenlandPeriph-Mapping/Mankoff
MD = readtable('D_sub.csv'); %Gt/yr (trimmed to only include glaciers in both datasets)
MH = readtable('th_sub.csv'); %m (trimmed to only include glaciers in both datasets)
MDgateID = table2array(MH(:,1));
Mg = shaperead('GrIS_gates.shp');
for i = 1:length(Mg)
   MgateID(i) = Mg(i).gate;
end


%box indices
land_BoxID = [484 485 222 140 142 163 103 608];
BoxIDs = unique(avg_thick(:,3));
cross_BoxIDs = [3 93 410 564];
Mankoff_BoxID = [8 101 115 135 156 162 162 174 174 176 238 240 247 248 259 259 259 254 ...
    307 308 310 312 318 322 322 349 349 349 355 356 357 356 357 ...
    363 364 365 366 363 364 365 366 369 424 429 430 431 433 ...
    435 443 447 448 459 460 463 466 467 470 469 471 578 580];
Mankoff_gates = [49 356 353 345 334 333 332 325 324 323 296 295 289 285 284 283 282 278 ...
    226 226 225 218 215 214 213 208 207 206 205 204 204 203 203 ...
    202 202 202 202 201 201 201 201 200 186 185 183 183 183 ...
    181 170 166 164 146 149 148 151 152 156 158 160 159 1 3];
GrIS_GIC_pair = [1 2 3 4 5 6 6 7 7 8 9 10 11 12 13 13 13 14 15 15 16 17 18 19 19 ...
    20 20 20 21 22 22 22 22 23 23 23 23 23 23 23 23 24 25 26 27 27 27 28 29 ...
    30 31 32 33 34 35 36 37 38 39 39 40 41];
GrIS_GIC_badpair = [23 41];

%extract average speed & thickness for each gate segment as thick
%end-members for curve fitting
unique_gates = unique(Mankoff_gates);
for i = 1:length(unique_gates)
    gaterefs = find(MgateID == unique_gates(i));
    
    %concatenate coordinates
    gateX = []; gateY = [];
    for j = 1:length(gaterefs)
        gateX = [gateX Mg(gaterefs(j)).x(~isnan(Mg(gaterefs(j)).x))];
        gateY = [gateY Mg(gaterefs(j)).y(~isnan(Mg(gaterefs(j)).y))];
    end
    
    %extract velocities
    M(i).gate = unique_gates(i);
    M(i).X = gateX; M(i).Y = gateY;
    M(i).speed = interp2(vxgrid,vygrid,V.z,gateX,gateY);
    clear gaterefs;
    
    %extract thicknesses
    gaterefs = find(MDgateID == unique_gates(i));
    for k = 1:length(gaterefs)
        M(i).discharge(1,k) = nanmean(table2array(MD(gaterefs(k),2:end-3)),2);
        M(i).thickness(1,k) = nanmean(table2array(MH(gaterefs(k),2:end-3)),2);
    end
    clear gaterefs;
    
    %create vectors of u,w,&H
    Mu(i) = nanmean(M(i).speed);
    Mw(i) = sqrt((M(i).X(1)-M(i).X(end)).^2 + (M(i).Y(1) - M(i).Y(end)).^2);
    Mh(i) = nanmean(M(i).thickness);
    Md(i) = nansum(M(i).discharge).*(10^9*1000)./917;
end

%filter-out ice sheet outlets as determined by Mankoff
for i = 1:length(BoxIDs)
    if ~isempty(find(Mankoff_BoxID == BoxIDs(i)))
        badref(find(avg_thick(:,3)==BoxIDs(i))) = 1;
    else
        badref(find(avg_thick(:,3)==BoxIDs(i))) = 0;
    end
end
avg_thick(badref==1,:) = [];

%identify the reference in the term structure & extract surface elevations
for j = 1:length(term); box(j) = term(j).BoxID; end %term BoxID vector
hvec = [];
for i = 1:length(BoxIDs)
    termref(i) = find(box == BoxIDs(i));
    if term(termref(i)).Zflag == 1
        hvec = [hvec; repmat(nanmedian(term(termref(i)).inlandZmed),length(find(avg_thick(:,3)==BoxIDs(i))),1)];
    else
        hvec = [hvec; repmat(nanmedian(term(termref(i)).seawardZmed),length(find(avg_thick(:,3)==BoxIDs(i))),1)];
    end
end

%extract widths for ALL glaciers, not just the best ones, & add to csv
for i = 1:length(BoxIDs)
    width = sqrt((term(termref(i)).gateX(1)-term(termref(i)).gateX(end)).^2 + (term(termref(i)).gateY(1)-term(termref(i)).gateY(end)).^2);
    avg_thick(avg_thick(:,3)==BoxIDs(i),8) = repmat(width,size(avg_thick(avg_thick(:,3)==BoxIDs(i),8)));
    clear width;
end
avg_thick(:,10) = avg_thick(:,6).*avg_thick(:,8);
avg_thick(:,9) = avg_thick(:,6)./avg_thick(:,8);
% avg_thick(:,11) = avg_thick(:,6).*avg_thick(:,8).*avg_thick(:,4);
avg_thick(:,11) = 240.*avg_thick(:,6).*avg_thick(:,4);
dlmwrite('avg_thick_vel_2.csv',avg_thick);

%filter questionable data
avg_thick(avg_thick(:,6)<15,:) = NaN;
avg_thick(avg_thick(:,4)<10,:) = NaN;
avg_thick(isnan(avg_thick(:,3)),:) = [];

% %extract data
% discharge = [avg_thick(:,3) avg_thick(:, 11)]; %D (yaxis)
% uxw = [avg_thick(:,3) avg_thick(:,10)]; %UW (xaxis)

%hone in on the glaciers with cross-flow data
n = [];
for j = 1:length(cross_BoxIDs)
    n = [n; find(avg_thick(:,3)==cross_BoxIDs(j))];
end
% uxw_lim = uxw(n,1:2); discharge_lim = discharge(n,1:2);

%flag thickness data according to region for color-coding 1=W,2=SE,3=E,4=NE,5=N
for i = 1:length(avg_thick(:,3))
    ref = find(box == avg_thick(i,3));
    regionref(i,1) = term(ref).regionFlag;
    region_colors(i,:) = region_cmap(regionref(i),:);
    clear ref;
end

%create glacier-averaged speed & thickness data
for i = 1:length(BoxIDs)
    GICu(i) = nanmean(avg_thick(find(avg_thick(:,3)==BoxIDs(i)),6));
    GICw(i) = nanmean(avg_thick(find(avg_thick(:,3)==BoxIDs(i)),8));
    GICh(i) = nanmean(avg_thick(find(avg_thick(:,3)==BoxIDs(i)),4));
    GICd(i) = nansum(avg_thick(find(avg_thick(:,3)==BoxIDs(i)),4).*avg_thick(find(avg_thick(:,3)==BoxIDs(i)),6).*240);
    GICregion_colors(i,:) = nanmean(region_colors(find(avg_thick(:,3)==BoxIDs(i)),:));
end

%assign regions to the GrIS data
for j=1:length(M)
    if nanmean(M(j).X) < 0
        if nanmean(M(j).Y) < -1e6 %west
            M(j).regionFlag = 1;
        else %north
            M(j).regionFlag = 5;
        end
    else
        if nanmean(M(j).Y) > -0.855e6 %north
            M(j).regionFlag = 5;
        elseif nanmean(M(j).Y) <= -0.855e6 && nanmean(M(j).Y) > -1.85e6 %northeast
            M(j).regionFlag = 4;
        elseif nanmean(M(j).Y) <= -1.85e6 && nanmean(M(j).Y) > -2.5e6 %central east
            M(j).regionFlag = 3;
        else %southeast
            M(j).regionFlag = 2;
        end
    end
end
%colormap GrIS data
for j = 1:length(M)
    GrISregion_colors(j,:) = region_cmap(M(j).regionFlag,:);
end

%re-export csv as a table with headers
Htable = array2table([avg_thick(:,1:4) avg_thick(:,6) avg_thick(:,8:11)]);
Htable.Properties.VariableNames = [{'X'},{'Y'},...
    {'BoxID'},{'Thickness (m)'},{'Speed (m/yr)'},{'Width (m)'},...
    {'Speed/Width'},{'SpeedxWidth'},{'SpeedxWidthxThickness'}];
Htable.Properties.VariableUnits = [{'m (Greenland Polar Stereo)'},{'m (Greenland Polar Stereo)'},...
    {'unitless'},{'m'},{'m yr^-1'},{'m'},...
    {'yr^-1'},{'m^2 yr^-1'},{'m^3 yr^-1'}];
writetable(Htable,'GreenlandGIC_UWH_table.csv');

%create a figure of all explored relationships
overviewfig = figure; set(gcf,'position',[50 50 1200 600]);
subplot(2,3,1);
scatter(avg_thick(:,6),avg_thick(:,4),36,region_colors,'filled','markeredgecolor','k','linewidth',1); hold on; %velocity vs thickness
plot(avg_thick(n,6),avg_thick(n,4),'+k','linewidth',1); %velocity vs thickness
scatter(GICu,GICh,36,GICregion_colors,'s','filled','markeredgecolor','k','linewidth',1.5); %average velocity vs thickness (all GICs with data)
scatter(Mu,Mh,36,GrISregion_colors,'d','filled','markeredgecolor','k','linewidth',1.5); %velocity vs thickness (all overlapping GrISs)
set(gca,'fontsize',14); grid on;
xlabel('speed (m yr^{-1})'); ylabel('thickness (m)');
subplot(2,3,2);
scatter(avg_thick(:,8),avg_thick(:,4),36,region_colors,'filled','markeredgecolor','k','linewidth',1); hold on;
plot(avg_thick(n,8),avg_thick(n,4),'+k','linewidth',1); %width vs thickness
scatter(GICw,GICh,36,GICregion_colors,'s','filled','markeredgecolor','k','linewidth',1.5); %width vs thickness (all GICs with data)
scatter(Mw,Mh,36,GrISregion_colors,'d','filled','markeredgecolor','k','linewidth',1.5); %width vs thickness (all overlapping GrISs)
set(gca,'fontsize',14); grid on;
xlabel('width (m)'); ylabel('thickness (m)');
subplot(2,3,3);
scatter(avg_thick(:,6).*avg_thick(:,8),avg_thick(:,4),36,region_colors,'filled','markeredgecolor','k','linewidth',1); hold on;
plot(avg_thick(n,10),avg_thick(n,4),'+k','linewidth',1); %velocity x width vs thickness
scatter(GICu.*GICw,GICh,36,GICregion_colors,'s','filled','markeredgecolor','k','linewidth',1.5); %average velocity X width vs thickness (all GICs with data)
scatter(Mu.*Mw,Mh,36,GrISregion_colors,'d','filled','markeredgecolor','k','linewidth',1.5); %velocity X width vs thickness (all overlapping GrISs)
set(gca,'fontsize',14); grid on;
xlabel('speed x width (m^2 yr^{-1})'); ylabel('thickness (m)');
subplot(2,3,4);
scatter(avg_thick(:,6)./avg_thick(:,8),avg_thick(:,4),36,region_colors,'filled','markeredgecolor','k','linewidth',1); hold on;
plot(avg_thick(n,9),avg_thick(n,4),'+k','linewidth',1); %velocity/width vs thickness
scatter(GICu./GICw,GICh,36,GICregion_colors,'s','filled','markeredgecolor','k','linewidth',1.5); %average velocity/width vs thickness (all GICs with data)
scatter(Mu./Mw,Mh,36,GrISregion_colors,'d','filled','markeredgecolor','k','linewidth',1.5); %velocity/width vs thickness (all overlapping GrISs)
set(gca,'fontsize',14); grid on;
xlabel('speed/width (yr^{-1})'); ylabel('thickness (m)');
subplot(2,3,5);
scatter(avg_thick(:,8)./avg_thick(:,6),avg_thick(:,4),36,region_colors,'filled','markeredgecolor','k','linewidth',1); hold on;
plot(avg_thick(n,8)./avg_thick(n,6),avg_thick(n,4),'+k','linewidth',1); %width/velocity vs thickness
scatter(GICw./GICu,GICh,36,GICregion_colors,'s','filled','markeredgecolor','k','linewidth',1.5); %width/average velocity vs thickness (all GICs with data)
scatter(Mw./Mu,Mh,36,GrISregion_colors,'d','filled','markeredgecolor','k','linewidth',1.5); %width/velocity vs thickness (all overlapping GrISs)
set(gca,'fontsize',14); grid on;
xlabel('width/speed (yr)'); ylabel('thickness (m)');
% subplot(2,3,6);
% % scatter((avg_thick(:,6).*avg_thick(:,8))./hvec,avg_thick(:,4),36,region_colors,'filled'); hold on;
% % plot(avg_thick(n,10)./hvec(n),avg_thick(n,4),'+k'); %(width x speed)/surface elevation vs thickness
% % title('(velocity x width)/elevation v. thickness');
% scatter(avg_thick(:,10),avg_thick(:,11),36,region_colors,'filled','markeredgecolor','k','linewidth',1); hold on;
% plot(avg_thick(n,10),avg_thick(n,11),'+k','linewidth',1); %(width x speed)/surface elevation vs thickness
% plot(GICu.*GICw,GICd,'+m','linewidth',1); %average velocity vs thickness (all GICs with data)
% plot(Mu.*Mw,Md,'xm','linewidth',1); %velocity vs thickness (all overlapping GrISs)
% set(gca,'fontsize',14); grid on;
% xlabel('speed x width (m^2 yr^{-1})'); ylabel('discharge (m^3 yr^{-1})');
saveas(gcf,'GreenlandGIC_test-empirical-UWH-relationship_subplots.png','png');
saveas(gcf,'GreenlandGIC_test-empirical-UWH-relationship_subplots.eps','epsc');

%plot distributions of the same combinations of velocities & widths for the
%entire dataset
U_distribution = []; W_distribution = []; WdivU_distribution = []; 
UxW_distribution = []; UdivW_distribution = [];
for i = 1:length(term)
    if isempty(find(land_BoxID == term(i).BoxID)) && isempty(find(Mankoff_BoxID == term(i).BoxID))
        glacier_width = sqrt((term(i).gateX(1) - term(i).gateX(end)).^2 + (term(i).gateY(1)-term(i).gateY(end)).^2);
        gate_width(i) = glacier_width;
        gate_vel = nanmedian(term(i).fluxV,1);
        U_distribution = [U_distribution; gate_vel'];
        W_distribution = [W_distribution; glacier_width];
        WdivU_distribution = [WdivU_distribution; glacier_width./gate_vel'];
        UxW_distribution = [UxW_distribution; glacier_width.*gate_vel'];
        UdivW_distribution = [UdivW_distribution; gate_vel'./glacier_width];
        clear glacier_width gate_vel;
    end
end

%plot histograms of velocity & width combinations
histogramfig = figure; set(gcf,'position',[50 50 1200 600]);
subplot(2,3,1);
histogram(U_distribution); xlabel('speed (m yr^{-1})');
subplot(2,3,2);
histogram(W_distribution); xlabel('width (m)');
subplot(2,3,3);
histogram(UxW_distribution); xlabel('speed x width (m^2 yr^{-1})');
subplot(2,3,4);
histogram(UdivW_distribution); xlabel('speed/width (yr^{-1})'); 
subplot(2,3,5);
histogram(WdivU_distribution); xlabel('width/speed (yr)'); 

%fit the data
totalu = [GICu Mu]; totalh = [GICh Mh];
ft = fittype( 'a*log10(x)+b*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.StartPoint = [0.769234158568257 0.73825160466599];
[UHfit, UHgof] = fit(totalu(~isnan(totalu))',totalh(~isnan(totalu))',ft,opts);
[GIC_UHfit, GIC_UHgof] = fit(GICu(~isnan(GICu))',GICh(~isnan(GICu))',ft,opts);
% ci = confint(UHfit,0.95); GIC_ci = confint(GIC_UHfit,0.95); 

%instead of using Matlab's built-in confidence interval solver, which
%appears to do a poor job with the custom equation, fit curves to the
%different data groupings as a way to define uncertainty
% thinrefs = find((avg_thick(:,8).*avg_thick(:,6)>=5*10^5 & regionref~=5) | regionref == 2);
% thickrefs = find(avg_thick(:,8).*avg_thick(:,6)<5*10^5 & regionref ~= 2);
% [GIC_UHfitmax, ~] = fit(avg_thick(thickrefs,6),avg_thick(thickrefs,4),ft,opts);
% [GIC_UHfitmin, ~] = fit(avg_thick(thinrefs,6),avg_thick(thinrefs,4),ft,opts);
thinrefs = find(totalh<=100); 
thickrefs = find(totalh>=200);
[GIC_UHfitmax, ~] = fit(totalu(thickrefs)',totalh(thickrefs)',ft,opts);
[GIC_UHfitmin, ~] = fit(totalu(thinrefs)',totalh(thinrefs)',ft,opts);

%compare estimated & measured thickness
eval(cd_to_root);
H_est = UHfit.a*log10(totalu)+UHfit.b*totalu; H_est(totalu==0) = 0;
HGIC_est = UHfit.a*log10(avg_thick(:,6))+UHfit.b*avg_thick(:,6); HGIC_est(avg_thick(:,6)==0) = 0;
figure; 
hp(1) = histogram(totalh-H_est); hold on;
hp(2) = histogram(avg_thick(:,4)-HGIC_est); hold on;
legend(hp,'GrIS & GIC (glacier averages)','GIC (all points)');
xlabel('Thickness residuals (m)'); ylabel('Count');
saveas(gcf,'GreenlandGIC_H-estimation_histograms.png','png'); saveas(gcf,'GreenlandGIC_H-estimation_histograms.eps','epsc');

%Plot a figure 
eval(cd_to_root);
Hfigure = figure; hold on; grid on; set(gcf,'Position',[100 0 800 800]); clear pl;
suba = subplot(2,1,1);
fill([1:1:15 sort(totalu(~isnan(totalu))) fliplr(sort(totalu(~isnan(totalu)))) 15:-1:1],...
    [GIC_UHfitmax.a*log10([1:1:15 sort(totalu(~isnan(totalu)))])+GIC_UHfitmax.b*[1:1:15 sort(totalu(~isnan(totalu)))] GIC_UHfitmin.a*log10([fliplr(sort(totalu(~isnan(totalu)))) 15:-1:1])+GIC_UHfitmin.b*[fliplr(sort(totalu(~isnan(totalu)))) 15:-1:1]],...
    [0.5 0.5 0.5]); hold on;
% fill([sort(GICu(~isnan(GICu))) fliplr(sort(GICu(~isnan(GICu))))],...
%     [GIC_ci(1,1)*log10(sort(GICu(~isnan(GICu))))+GIC_ci(1,2)*sort(GICu(~isnan(GICu))) GIC_ci(2,1)*log10(fliplr(sort(GICu(~isnan(GICu)))))+GIC_ci(2,2)*fliplr(sort(GICu(~isnan(GICu))))],...
%     [0.5 0.5 0.5]); hold on;
pl(1) = plot([1:1:15 sort(totalu(~isnan(totalu)))],UHfit.a*log10([1:1:15 sort(totalu(~isnan(totalu)))])+UHfit.b*[1:1:15 sort(totalu(~isnan(totalu)))],'-k','linewidth',2); hold on;
for k = 1:length(region_cmap)
pl(k+1) = scatter(GICu(1),GICh(1),40,region_cmap(k,:),'o','filled','linewidth',1,'markeredgecolor','k'); hold on;
end
scatter(avg_thick(:,6),avg_thick(:,4),40,region_colors,'+','linewidth',1); hold on;
pl(length(pl)+1) = scatter(GICu(1),GICh(1),36,'k','+','linewidth',1); hold on;
pl(length(pl)+1) = scatter(GICu(1),GICh(1),36,'w','s','filled','linewidth',1,'markeredgecolor','k'); hold on;
pl(length(pl)+1) = scatter(Mu(1),Mh(1),36,'w','d','filled','linewidth',1,'markeredgecolor','k'); hold on;
scatter(GICu,GICh,48,GICregion_colors,'s','filled','linewidth',1,'markeredgecolor','k'); hold on;
scatter(Mu,Mh,48,GrISregion_colors,'d','filled','linewidth',1,'markeredgecolor','k'); hold on;
set(gca,'FontSize',20,'ylim',[0 520],'xlim',[0 700],'xticklabel',[]);
ylabel('thickness (m)','fontsize',20);
legend(pl,...
    ['H=',num2str(round(UHfit.a,2)),'log10(U)+',num2str(round(UHfit.b,2)),'U'],...
    'west','southeast','central east','northeast','north','GIC_{observations}','GIC_{averages}','GrIS_{averages}');
legend('Location','southeast', 'FontSize',14);
subb = subplot(2,1,2);
hist_U = histogram(U_distribution,[0:10:600],'FaceColor',[0 0 0],'FaceAlpha',1);
set(subb,'fontsize',20,'xlim',[0 600],'ylim',[0 4000]);
set(gca,'YScale','log','ytick',[1 10 100 1000],'yticklabel',[1 10 100 1000]);
xlabel('speed (m yr^{-1})','fontsize',20); ylabel('count','fontsize',20); 
set(subb,'position',[0.1300 0.1100 0.7750 0.1800]);
subplot(subb); text(min(get(gca,'xlim'))+0.02*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.50*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'b)','FontSize',20);
set(suba,'position',[0.1300 0.3200 0.7750 0.6300]);
subplot(suba); text(min(get(gca,'xlim'))+0.02*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.975*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'a)','FontSize',20);
saveas(gcf,'GreenlandGIC_H-estimation_scatterplot.png','png'); saveas(gcf,'GreenlandGIC_H-estimation_scatterplot.eps','epsc');
 
%Save the function
% UW_lims = [min(uxw_lim(:,2)) max(uw)];
% regional_UWfits = []; 
% for i = 1:5
%     if ~isempty(region(i).UWfit)
%         regional_UWfits = [regional_UWfits region(i).UWfit.p1]; 
%     end
% end
eval(cd_to_root);
UHfit_min = GIC_UHfitmin; UHfit_max = GIC_UHfitmax; 
save('H_vs_U_logarithmic-functions.mat','UHfit','UHfit_min','UHfit_max','-v7.3');

%OLD
% save('D_vs_UW_linear-functions.mat','UWfit','regional_UWfits','-v7.3');
% save('H_vs_WdivU_linear-functions.mat','H_fit','H_ci','H_sigma','-v7.3');

 
%% fit trendlines (OLD... delete when settled on D estimation approach)
% u = avg_thick(:,6); uw = avg_thick(:,10); %bin speed (m/yr) and bin speed x glacier width (m^2/yr)
% H = avg_thick(:,4);
% uwH = 240.*avg_thick(:,6).*avg_thick(:,4); %bin width x bin speed x bin thickness (m^3/yr)
% 
% %fit a polynomial to the data: use the polyfitB0 function to force the x&y intercepts to (0,0)
% sigma = 2; %convert from 1 sigma (68% confidence interval) to 95% confidence interval
% % [p,S,mu] = polyfitB0(uw,uwH,2,0); %2 = degree, 0 = y-intercept
% % UWfit = p;
% opts = fitoptions('Method','LinearLeastSquares');
% opts.Lower = [-Inf 0]; opts.Upper = [Inf 0];
% opts.Robust = 'Bisquare'; %weight squared differences based on distance from line
% ft = fittype('poly1');
% [UWfit,UWgof] = fit(uw,uwH,ft,opts); 
% if sigma == 2
%     ci = confint(UWfit,0.95);
% else
%     ci = confint(UWfit,0.68);
% end
% ci(isnan(ci)) = 0;
% %fit to each region
% for i = 1:5
%     if length(find(regionref==i)) > 2
%         [region(i).UWfit,region(i).UWgof] = fit(uw(regionref==i),uwH(regionref==i),ft,opts); 
%     end
% end
% 
% 
% %set the function outside the data domain as a linear polynomial to avoid
% %erroneous negative relationships between H and U/W at very large values
% % [yval,yerr] = polyval(p,[max(uw(uw<max(uw))) max(uw)],S,mu);
% % UWext.p1 = (yval(2)-yval(1))./(max(uw)-max(uw(uw<max(uw))));
% % UWext.p2 = yval(2)-UWext.p1*max(uw);
% % UWext.cil_p1 = ((yval(2)-sigma*yerr(2))-(yval(1)-sigma*yerr(1)))./(max(uw)-max(uw(uw<max(uw))));
% % UWext.cil_p2 = (yval(2)-sigma*yerr(2))-UWext.cil_p1*max(uw);
% % UWext.ciu_p1 = ((yval(2)+sigma*yerr(2))-(yval(1)+sigma*yerr(1)))./(max(uw)-max(uw(uw<max(uw))));
% % UWext.ciu_p2 = (yval(2)+sigma*yerr(2))-UWext.ciu_p1*max(uw);
% %evaluate function over a fine interval for plotting
% clear UWH;
% xsub = sort(uw);
% UWH = UWfit.p1.*xsub+UWfit.p2;
% ci2_lower = ci(1,1).*xsub+ci(1,2); ci2_upper = ci(2,1).*xsub+ci(2,2);
% % xsub = [0:max(uw)/1000:max(uw)];
% % [UWH,UWHerr] = polyval(UWfit,xsub,S,mu);
% % ci2_lower = UWH-sigma*UWHerr; ci2_upper = UWH+sigma*UWHerr; 
% 
% 
% %Plot a figure 
% Dfigure = figure; hold on; grid on; set(gcf,'Position',[0 0 800 800]);
% suba = subplot(2,1,1);
% for k = 1:length(region_cmap)
% pl(k) = scatter(avg_thick(1,10),avg_thick(1,11),36,region_cmap(k,:),'filled','linewidth',1,'markeredgecolor','k'); hold on;
% end
% scatter(avg_thick(:,10),avg_thick(:,11),36,region_colors,'filled','markeredgecolor','k','linewidth',1); hold on;
% pl(length(region_cmap)+1) = plot(xsub,UWH,'-k','linewidth',2); hold on; 
% % plot(xsub,ci2_lower,'--k','linewidth',2); hold on; plot(xsub,ci2_upper,'--k','linewidth',2); hold on; 
% for i = 1:5
%     if ~isempty(region(i).UWfit)
%         plot(xsub,region(i).UWfit.p1.*xsub+region(i).UWfit.p2,'-k','linewidth',2); hold on;
%         plot(xsub,region(i).UWfit.p1.*xsub+region(i).UWfit.p2,'-','color',region_cmap(i,:),'linewidth',1); hold on;
%     end
% end
% %plot the extrapolated equations to check they make sense
% % plot([max(uw) 4e5],UWext.p1.*[max(uw) 4e5]+UWext.p2,'-k','linewidth',2); hold on;
% % plot([max(uw) 4e5],UWext.cil_p1.*[max(uw) 4e5]+UWext.cil_p2,'--k','linewidth',2); hold on;
% % plot([max(uw) 4e5],UWext.ciu_p1.*[max(uw) 4e5]+UWext.ciu_p2,'--k','linewidth',2); hold on;
% % col = parula(20);
% % for i=1:length(uxw_lim(:,1))
% %     if uxw_lim(i,1)==003
% %         p2 = plot(uxw_lim(i,2), discharge_lim(i,2), 's','color',col(4,:),'markersize',12, 'LineWidth', 3);
% %         %set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% %     elseif uxw_lim(i,1)==093
% %         p3 = plot(uxw_lim(i,2), discharge_lim(i,2), 's','color',col(8,:),'markersize',12, 'LineWidth', 3);
% %         %set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% %     elseif uxw_lim(i,1)==410
% %         p4 = plot(uxw_lim(i,2), discharge_lim(i,2), 's','color',col(13,:),'markersize',12, 'LineWidth', 3);
% %         %set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% %     elseif uxw_lim(i,1)==564
% %         p5 = plot(uxw_lim(i,2), discharge_lim(i,2), 's','color',col(18,:),'markersize',12, 'LineWidth', 3);
% %         %set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% %     end 
% % end 
% set(gca,'FontSize',20,'ylim',[0 ceil(max(uwH))],'xlim',[0 max(ceil(uw))]); 
% ylims = get(gca,'ylim'); yticks = get(gca,'ytick'); xticks = get(gca,'xtick');
% set(gca,'ylim',[0 max(ylims)],'ytick',yticks,'yticklabel',yticks./10^6,...
%     'xtick',xticks,'xticklabel',[]);
% ylabel('discharge_{bin} (10^6 m^3 yr^{-1})','fontsize',20);
% legend(pl,'west','southeast','central east','northeast','north',['D=',num2str(UWfit.p1),'UW']);
% % xlabel('speed_{bin} x width_{glacier}','fontsize',20); ylabel('thickness_{bin} x speed_{bin} x width_{glacier}','fontsize',20);
% % if p(2) == 0
% %     legend([pl(1) p2 p3 p4 p5],...
% %         ['y=',num2str(p(1)/10^6),'x10^6 x^2'],...
% %         'BoxID 003','BoxID 093','BoxID 410','BoxID 564');
% % else
% %     legend([pl(1) p2 p3 p4 p5],...
% %         ['y=',num2str(p(1)/10^6),'x10^6 x^2 + ',num2str(p(2)/10^6),' x'],...
% %         'BoxID 003','BoxID 093','BoxID 410','BoxID 564');
% % end
% legend('Location','northeast', 'FontSize',14);
% subb = subplot(2,1,2);
% histUW = histogram(UxW_distribution,'facecolor',[0.5 0.5 0.5]); 
% set(gca,'FontSize',20,'xlim',[0 max(ceil(uw))],'ylim',[0 max(ceil(histUW.Values/100)*100)],...
%     'xtick',xticks,'xticklabel',xticks/10^6); 
% xlabel('speed_{bin} x width_{glacier} (10^6 m^2 yr^{-1})','fontsize',20); ylabel('count ','fontsize',20);
% set(subb,'position',[0.1300 0.1100 0.7750 0.1800]);
% subplot(subb); text(min(get(gca,'xlim'))+0.95*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.10*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'b)','FontSize',20);
% set(suba,'position',[0.1300 0.3200 0.7750 0.6300]);
% subplot(suba); text(min(get(gca,'xlim'))+0.95*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.04*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'a)','FontSize',20);
% saveas(gcf,'GreenlandGIC_D-estimation_scatterplot.png','png'); saveas(gcf,'GreenlandGIC_D-estimation_scatterplot.eps','epsc');
% 
% %try a simpler H estimation approach
% WdivU = avg_thick(:,8)./avg_thick(:,6); H = avg_thick(:,4);
% % ft = fittype('poly1');
% % opts = fitoptions('Method','LinearLeastSquares');
% % opts.Lower = [-Inf 0]; opts.Upper = [Inf 0];
% % [WdivUfit, WdivUgof] = fit(WdivU,H,ft,opts); ci = confint(WdivUfit,0.95);
% H_fit = nanmedian(H./WdivU); 
% H_ci = [(nanmedian(H./WdivU)-1.4826*mad(H./WdivU,1)); (nanmedian(H./WdivU)+1.4826*mad(H./WdivU,1))];
% 
% %compare estimated & measured thickness
% % Dest = polyval(UWfit,avg_thick(:,10),S,mu); H_UWest = Dest./avg_thick(:,10);
% Dest = UWfit.p1.*uw + UWfit.p2; H_UWest = Dest./(240*u);
% H_UWest(avg_thick(:,10)>max(uw)) = (UWext.p1.*avg_thick((avg_thick(:,10)>max(uw)),10) + UWext.p2)./avg_thick((avg_thick(:,10)>max(uw)),10);
% H_WdivUest = nanmedian(H./WdivU)*WdivU; H_sigma = nanstd(H-H_WdivUest);
% figure; 
% hp(1) = histogram(H-H_UWest); hold on;
% hp(2) = histogram(H-H_WdivUest); hold on;
% legend(hp,'linear UW v. D fit','linear W/U v. H fit');
% xlabel('Thickness residuals (m)'); ylabel('Count');
% saveas(gcf,'GreenlandGIC_H-estimation_histograms.png','png'); saveas(gcf,'GreenlandGIC_H-estimation_histograms.eps','epsc');
% 
% %Plot a figure 
% Hfigure = figure; hold on; grid on; set(gcf,'Position',[100 0 800 800]); clear pl;
% suba = subplot(2,1,1);
% fill([sort(WdivU); flipud(sort(WdivU))],[(nanmedian(H./WdivU)+1.4826*mad(H./WdivU,1))*sort(WdivU); (nanmedian(H./WdivU)-1.4826*mad(H./WdivU,1))*flipud(sort(WdivU))],[0.5 0.5 0.5]);
% pl(1) = plot(sort(WdivU),nanmedian(H./WdivU)*sort(WdivU),'-k','linewidth',2); hold on;
% for k = 1:length(region_colors)
% pl(k+1) = scatter(avg_thick(1,8)./avg_thick(1,6),avg_thick(1,4),36,region_colors(k,:),'s','filled','linewidth',1,'markeredgecolor','k'); hold on;
% end
% scatter(avg_thick(:,8)./avg_thick(:,6),avg_thick(:,4),36,region_colors,'s','filled','linewidth',1,'markeredgecolor','k'); hold on;
% set(gca,'FontSize',14,'ylim',[0 400]);
% ylabel('thickness_{bin}','fontsize',20);
% legend(pl,...
%     ['H=',num2str(nanmedian(H./WdivU)),' (W/U)'],'west','southeast','central east','northeast','north');
% legend('Location','northwest', 'FontSize',14);
% subb = subplot(2,1,2);
% hist_WdivU = histogram(WdivU_distribution,[0:50:max(gate_width)],'FaceColor',[0.5 0.5 0.5]);
% set(gca,'fontsize',20,'xlim',[0 1000]);
% xlabel('width_{glacier}/speed_{bin}','fontsize',20); ylabel('count','fontsize',20); 
% set(subb,'position',[0.1300 0.1100 0.7750 0.1800]);
% subplot(subb); text(min(get(gca,'xlim'))+0.9*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.10*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'b)','FontSize',20);
% set(suba,'position',[0.1300 0.3200 0.7750 0.6300]);
% subplot(suba); text(min(get(gca,'xlim'))+0.9*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.04*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'a)','FontSize',20);
% % saveas(gcf,'GreenlandGIC_H-estimation_scatterplot.png','png'); saveas(gcf,'GreenlandGIC_H-estimation_scatterplot.eps','epsc');
%  
% %Save the function
% % clear f;
% % f.xlims = [min(uxw_lim(uxw_lim(:,1)==003,2)) max(uw)];
% % f.fit = f2; f.ext = f2ext; f.int = f2min;
% % save('D_vs_UW_parab-functions.mat','f','-v7.3');
% UW_lims = [min(uxw_lim(:,2)) max(uw)];
% regional_UWfits = []; 
% for i = 1:5
%     if ~isempty(region(i).UWfit)
%         regional_UWfits = [regional_UWfits region(i).UWfit.p1]; 
%     end
% end
% eval(cd_to_root);
% save('D_vs_UW_linear-functions.mat','UWfit','regional_UWfits','-v7.3');
% save('H_vs_WdivU_linear-functions.mat','H_fit','H_ci','H_sigma','-v7.3');

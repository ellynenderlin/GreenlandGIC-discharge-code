%%% Use this script to extract discharge time series using pre-drawn
%%% centerlines & flux gates.
%%%  

%% initialize (run every time)
clearvars; close all; warning off;
disp('Initialized discharge analysis code');

%add the path to supporting Matlab functions
addpath('/users/ellynenderlin/Research/miscellaneous/general-code');

%specify the root directory for project-specific files
root_dir = '/Users/ellynenderlin/Research/NASA_GreenlandPeriph-Mapping/'; %including trailing /
cd(root_dir);
%specify the root directory for generic supporting files
misc_dir = '/Users/ellynenderlin/Research/miscellaneous/'; %including trailing /

%grab the velocity years from the velocity map filenames
vel_mosaic = dir([misc_dir,'Greenland-ITSLIVE*']);
D_yrs = str2num(vel_mosaic.name(end-8:end-5)):1:str2num(vel_mosaic.name(end-3:end));

%set-up day of year variables
modays = [31 28 31 30 31 30 31 31 30 31 30 31]; cumdays = [0 cumsum(modays(1:end-1))];
leap_modays = [31 29 31 30 31 30 31 31 30 31 30 31]; leap_cumdays = [0 cumsum(leap_modays(1:end-1))];

%load term data structure containing flux gates & empirical scaling equation
load([root_dir,'Greenland_GIC_centerlines.mat']);

%flag Mankoff glaciers & those identified as land-terminating
land_BoxID = [484 485 222 140 142 163 103 608]; %jth index = [428 429 138 47 49 72 6 566]; 
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
GrIS_GIC_badpair = [23 41]; %Mankoff glaciers have flow-following gates

%% estimate discharge

%determine what you want to do
answer = questdlg('Do you want to calculate discharge or only update figures?',...
    'Discharge Handling','1) Calculate Discharge','2) Only update figures','1) Calculate Discharge');
switch answer
    case '1) Calculate Discharge'
        %output from GreenlandGIC_thickness_extrapolation_curves.m
        load([root_dir,'H_vs_U_logarithmic-functions.mat']);
        processing_flag = 0;
        disp('Loaded empirical thickness scaling functions');
    case '2) Only update figures'
        processing_flag = 1;
        disp('Working with existing discharge estimates');
end

Dtotal = []; Dtotal_err = []; Dx = []; Dy = []; DBoxID = []; Dupper = []; Dlower = [];
% counter = 0;
for i = 1:length(term)
    disp(['glacier ',num2str(i),' of ',num2str(length(term))]);
    BoxID(i) = term(i).BoxID;
    
    %calculate discharge only if specified
    if processing_flag == 0
        %create thickness cross sections at the flux gate
        glacier_width = sqrt((term(i).gateX(1) - term(i).gateX(end)).^2 + (term(i).gateY(1)-term(i).gateY(end)).^2);
        gate_vel = nanmedian(term(i).fluxV,1); gate_width = term(i).fluxW;
        gate_H = NaN(size(gate_vel)); gate_Hmax = NaN(size(gate_vel)); gate_Hmin = NaN(size(gate_vel));
        
        %apply logarithmic U vs H trendline
        gate_H(1,:) = UHfit.a.*log10(gate_vel(1,:))+UHfit.b.*gate_vel(1,:);
        gate_Hmin(1,:) = UHfit_min.a.*log10(gate_vel(1,:))+UHfit_min.b.*gate_vel(1,:);
        gate_Hmax(1,:) = UHfit_max.a.*log10(gate_vel(1,:))+UHfit_max.b.*gate_vel(1,:);
        
        %add H estimates to structure
        gate_H(gate_H<0) = 0; gate_Hmin(gate_Hmin<0) = 0; gate_Hmax(gate_Hmax<0) = 0;
        term(i).fluxH = gate_H; term(i).fluxHmax = gate_Hmax; term(i).fluxHmin = gate_Hmin;
        
        %construct discharge time series
        gate_D = term(i).fluxV.*repmat(term(i).fluxH,size(term(i).fluxV,1),1).*repmat(gate_width,size(term(i).fluxV,1),1);
        D = nansum(gate_D,2);
        %propagate errors
        gate_Verr = term(i).fluxVerr;
        gate_Werr = zeros(size(gate_width)); gate_Werr(1) = gate_Werr(1) + 15; gate_Werr(end) = gate_Werr(end) + 15; term(i).fluxWerr = gate_Werr;
        gate_Hposerr = repmat((term(i).fluxHmax-term(i).fluxH),size(term(i).fluxV,1),1); gate_Hnegerr = repmat((term(i).fluxHmin-term(i).fluxH),size(term(i).fluxV,1),1);
        gate_VHposcov = cov(gate_Hposerr,gate_Verr); gate_VHnegcov = cov(gate_Hnegerr,gate_Verr);
        gate_Dposerr = abs(gate_D).*sqrt((gate_Verr./term(i).fluxV).^2 + gate_Hposerr./term(i).fluxH.^2 +...
            2*repmat(gate_VHposcov(1,2),size(term(i).fluxV))./(term(i).fluxV.*repmat(term(i).fluxH,size(term(i).fluxV,1),1)) +...
            repmat(term(i).fluxWerr./term(i).fluxW,size(term(i).fluxV,1),1).^2);
        gate_Dnegerr = abs(gate_D).*sqrt((gate_Verr./term(i).fluxV).^2 + gate_Hnegerr./term(i).fluxH.^2 +...
            2*repmat(gate_VHnegcov(1,2),size(term(i).fluxV))./(term(i).fluxV.*repmat(term(i).fluxH,size(term(i).fluxV,1),1)) +...
            repmat(term(i).fluxWerr./term(i).fluxW,size(term(i).fluxV,1),1).^2);
        Dposerr = sqrt(nansum(gate_Dposerr.^2,2)); Dnegerr = sqrt(nansum(gate_Dnegerr.^2,2));
        clear glacier_width gate_*;
        term(i).fluxD = D; term(i).fluxDmax = term(i).fluxD + Dposerr; term(i).fluxDmin = term(i).fluxD - Dnegerr;
        term(i).fluxDerr = nanmean([Dposerr Dnegerr],2); term(i).fluxDyrs = D_yrs';
        term(i).fluxVavg(term(i).fluxVavg <= 0) = NaN;
        
        %clear the old variables
        clear D Derr gate*
    end
    
    %add to the entire GIC timeseries if it is confirmed as
    %marine-terminating and not in the Mankoff GrIS timeseries
    if isempty(find(land_BoxID == term(i).BoxID)) && isempty(find(Mankoff_BoxID == term(i).BoxID))
        Dtotal = [Dtotal; term(i).fluxD']; Dtotal_err = [Dtotal_err; term(i).fluxDerr'];
        Dx = [Dx; term(i).gateXavg]; Dy = [Dy; term(i).gateYavg];
        DBoxID = [DBoxID; term(i).BoxID];
        Dupper = [Dupper; term(i).fluxDmax'];
        Dlower = [Dlower; term(i).fluxDmin'];
    end
    
end
%save discharge estimates to term structure
if processing_flag == 0
    save([root_dir,'Greenland_GIC_centerlines.mat'],'term','-v7.3'); disp('date resaved');
end

% create discharge timeseries
% D_ann_avg = nanmean(Dtotal,1); D_median = nanmedian(Dtotal,1);
D_pre = nanmean(nansum(Dtotal(:,1:14),1));%mean of annual average discharge from 1985-1998, before mass anomaly
D_post = nanmean(nansum(Dtotal(:,15:end),1));
% D_anom = D_ann_avg(:,15:end)-D_pre; %subtract mean discharge pre-anomaly from annual average for remainder of record (1999-2018)
disp(['Greenland GIC discharge: ']);
disp(['Means: 1985-1998 = ',num2str(D_pre.*917./(1000*10^9)),' & 1999-2018 = ',num2str(D_post.*917./(1000*10^9)),' Gt/yr']);
disp(['1985-1998 median +/- MAD = ',num2str(nanmedian(nansum(Dtotal(:,1:14),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(Dtotal(:,1:14),1),1).*917./(1000*10^9)),' Gt/yr']);
disp(['1999-2018 median +/- MAD = ',num2str(nanmedian(nansum(Dtotal(:,15:end),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(Dtotal(:,15:end),1),1).*917./(1000*10^9)),' Gt/yr']);

% define RGB matrices for colors on discharge plot
cmap = [0 0 0; 123,50,148; 194,165,207; 166,219,160; 0,136,55]./255;

%filter to only include glaciers with complete timeseries
Dfiltered = NaN(length(term),length(term(1).fluxD)); Dfiltered_err = NaN(length(term),length(term(1).fluxD));
Dfiltered_upper = NaN(length(term),length(term(1).fluxD)); Dfiltered_lower = NaN(length(term),length(term(1).fluxD));
for i = 1:length(term)
    if sum(isnan(term(i).fluxVavg)) == 0 && isempty(find(land_BoxID == term(i).BoxID)) && isempty(find(Mankoff_BoxID == term(i).BoxID))
        Dfiltered(i,:) = term(i).fluxD'; Dfiltered_err(i,:) = term(i).fluxDerr';
        Dfiltered_upper(i,:) = term(i).fluxDmax'; Dfiltered_lower(i,:) = term(i).fluxDmin';
%     elseif term(i).MankoffFlag == 1
%         Dfiltered(i,:) = zeros(1,34);
    end
end
disp(['glaciers with complete time series = ', num2str(sum(~isnan(Dfiltered(:,1))))]);
disp(['Means: 1985-1998 = ',num2str(nanmean(nansum(Dfiltered(:,1:14),1)).*917./(1000*10^9)),' & 1999-2018 = ',num2str(nanmean(nansum(Dfiltered(:,15:end),1)).*917./(1000*10^9)),' Gt/yr']);
disp(['1985-1998 median +/ MAD = ',num2str(nanmedian(nansum(Dfiltered(:,1:14),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(Dfiltered(:,1:14),1),1).*917./(1000*10^9)),' Gt/yr']);
disp(['1999-2018 median +/- MAD = ',num2str(nanmedian(nansum(Dfiltered(:,15:end),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(Dfiltered(:,15:end),1),1).*917./(1000*10^9)),' Gt/yr']);

%fill in data gaps from 1985-1998 w/ 1985-1998 means & from 1999-2018
%with 1999-2018 means
Dfilled = Dtotal; Dfilled_err = Dtotal_err; Dfilled_upper = Dupper; Dfilled_lower = Dlower;
for i = 1:size(Dtotal)
    glacierD = Dtotal(i,:); glacierD(glacierD<0) = 0; glacierDerr = Dtotal_err(i,:); 
    glacierDupper = Dupper(i,:); glacierDlower = Dlower(i,:); 
    %before 1999
    Dearly_fill = nanmean(glacierD(glacierD(1:14)~=0)); Dearly_fillerr = sqrt(sum(glacierDerr(glacierD(1:14)~=0).^2)./(length(glacierD(glacierD(1:14)~=0)).^2));
    Dearly_fillupper = nanmean(glacierDupper(glacierD(1:14)~=0)); Dearly_filllower = nanmean(glacierDlower(glacierD(1:14)~=0)); 
    Dfilled(i,glacierD(1:14)==0) = Dearly_fill; Dfilled_err(i,glacierD(1:14)==0) = Dearly_fillerr;
    Dfilled_upper(i,glacierD(1:14)==0) = Dearly_fillupper; Dfilled_lower(i,glacierD(1:14)==0) = Dearly_filllower;
    
    %1999 to present
    Dlate_fill = nanmean(glacierD(glacierD(15:end)~=0)); Dlate_fillerr = sqrt(sum(glacierDerr(glacierD(15:end)~=0).^2)./(length(glacierD(glacierD(15:end)~=0)).^2));
    Dlate_fillupper = nanmean(glacierDupper(glacierD(15:end)~=0)); Dlate_filllower = nanmean(glacierDlower(glacierD(15:end)~=0)); 
    Dfilled(i,14+find(glacierD(15:end)==0)) = Dlate_fill; Dfilled_err(i,14+find(glacierD(15:end)==0)) = Dlate_fillerr;
    Dfilled_upper(i,14+find(glacierD(15:end)==0)) = Dlate_fillupper; Dfilled_lower(i,14+find(glacierD(15:end)==0)) = Dlate_filllower;
    clear glacierD Dearly* Dlate*;
    
    %add to the structure
    term_ref = find(DBoxID(i)-BoxID == 0);
    term(term_ref).fluxD_filled = Dfilled(i,:)'; term(term_ref).fluxDerr_filled = Dfilled_err(i,:)';
    term(term_ref).fluxDmax_filled = Dfilled_upper(i,:)'; term(term_ref).fluxDmin_filled = Dfilled_lower(i,:)';
    
    %add a flag to identify if the time series is complete w/ gap-filling
    if sum(isnan(Dfilled(i,:))) > 1; complete_flag(i) = 0; else complete_flag(i) = 1; end
end
Dfilled_err = sqrt(nansum((Dfilled-Dfilled_lower).^2));
disp('gap-filled time series...')
disp(['Means: 1985-1998 = ',num2str(nanmean(nansum(Dfilled(:,1:14),1)).*917./(1000*10^9)),' & 1999-2018 = ',num2str(nanmean(nansum(Dfilled(:,15:end),1)).*917./(1000*10^9)),' Gt/yr']);
disp(['1985-1998 median +/ MAD = ',num2str(nanmedian(nansum(Dfilled(:,1:14),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(Dfilled(:,1:14),1),1).*917./(1000*10^9)),' Gt/yr']);
disp(['1999-2018 median +/- MAD = ',num2str(nanmedian(nansum(Dfilled(:,15:end),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(Dfilled(:,15:end),1),1).*917./(1000*10^9)),' Gt/yr']);

%estimate average discharge for 1985-1998 for glaciers w/o data by assuming
%that the relative change is the same for those glaciers as the ones with
%data over the full time period
Dlate_incomplete = Dfilled(complete_flag == 0,15:end);
Dlate_complete = Dfilled(complete_flag == 1,15:end); Dearly_complete = Dfilled(complete_flag == 1,1:14); 
Dearly_incomplete = nanmean(Dlate_incomplete,2)./(1+nanmedian((nanmean(Dlate_complete,2)-nanmean(Dearly_complete,2))./nanmean(Dearly_complete,2)));
Dearly_incompleteMIN = nanmean(Dlate_incomplete,2)./(1+(nanmedian((nanmean(Dlate_complete,2)-nanmean(Dearly_complete,2))./nanmean(Dearly_complete,2))+mad((nanmean(Dlate_complete,2)-nanmean(Dearly_complete,2))./nanmean(Dearly_complete,2),1)));
Dearly_incompleteMAX = nanmean(Dlate_incomplete,2)./(1+(nanmedian((nanmean(Dlate_complete,2)-nanmean(Dearly_complete,2))./nanmean(Dearly_complete,2))-mad((nanmean(Dlate_complete,2)-nanmean(Dearly_complete,2))./nanmean(Dearly_complete,2),1)));
disp('gap-filled time series & uniform relative change...')
disp(['Means: 1985-1998 = ',num2str(nanmean(nansum(Dfilled(:,1:14),1)+nansum(Dearly_incomplete)).*917./(1000*10^9)),' & 1999-2018 = ',num2str(nanmean(nansum(Dfilled(:,15:end),1)).*917./(1000*10^9)),' Gt/yr']);
disp(['1985-1998 median +/ MAD = ',num2str(nanmedian(nansum(Dfilled(:,1:14),1)+nansum(Dearly_incomplete)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(Dfilled(:,1:14),1)+nansum(Dearly_incomplete),1).*917./(1000*10^9)),' Gt/yr']);
disp(['1999-2018 median +/- MAD = ',num2str(nanmedian(nansum(Dfilled(:,15:end),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(Dfilled(:,15:end),1),1).*917./(1000*10^9)),' Gt/yr']);


%plot glaciers with average annual D > 0.05 Gt/yr
Dbig = NaN(size(Dtotal)); Dbig_err = NaN(size(Dtotal_err)); Dbig_upper = NaN(size(Dtotal_err)); Dbig_lower = NaN(size(Dtotal_err));
for i = 1:size(Dtotal)
    if nanmean(Dtotal(i,:)).*917./(1000*10^9) > 0.05 && isempty(find(land_BoxID == term(i).BoxID)) && isempty(find(Mankoff_BoxID == term(i).BoxID))
        Dbig(i,:) = Dtotal(i,:); Dbig_err(i,:) = Dtotal_err(i,:);
        Dbig_upper(i,:) = Dupper(i,:); Dbig_lower(i,:) = Dlower(i,:);
    end
end
disp(['glaciers D>0.05 Gt/yr= ', num2str(sum(~isnan(Dbig(:,1))))]);
disp(['Means: 1985-1998 = ',num2str(nanmean(nansum(Dbig(:,1:14),1)).*917./(1000*10^9)),' & 1999-2018 = ',num2str(nanmean(nansum(Dbig(:,15:end),1)).*917./(1000*10^9)),' Gt/yr']);
disp(['1985-1998 median +/ MAD = ',num2str(nanmedian(nansum(Dbig(:,1:14),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(Dbig(:,1:14),1),1).*917./(1000*10^9)),' Gt/yr']);
disp(['1999-2018 median +/- MAD = ',num2str(nanmedian(nansum(Dbig(:,15:end),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(Dbig(:,15:end),1),1).*917./(1000*10^9)),' Gt/yr']);


%PLOT TIMESERIES
%total discharge timeseries
figure; set(gcf,'position',[100 100 800 800]);
suba = subplot(2,1,1);
pl(1) = plot(D_yrs,nansum(Dtotal).*917./(1000*10^9),'color',cmap(1,:),'linewidth',1.5); hold on;
errorbar(D_yrs,(nansum(Dtotal).*917./(1000*10^9)),sqrt(nansum((Dtotal-Dlower).^2)).*917./(1000*10^9),sqrt(nansum((Dupper-Dtotal).^2)).*917./(1000*10^9),'color',cmap(1,:),'linewidth',1.5);
%gaps filled with mean from 1985-1998 or 1999-2018
pl(2) = plot(D_yrs,nansum(Dfilled).*917./(1000*10^9),'color',cmap(2,:),'linewidth',1.5); hold on;
errorbar(D_yrs,(nansum(Dfilled).*917./(1000*10^9)),Dfilled_err.*917./(1000*10^9),Dfilled_err.*917./(1000*10^9),'color',cmap(2,:),'linewidth',1.5);
%gaps filled AND early discharge back-calculated assuming same relative change
pl(3) = plot(D_yrs,[(nansum(Dfilled(:,1:14))+nansum(Dearly_incomplete)) nansum(Dfilled(:,15:end))].*917./(1000*10^9),'color',cmap(3,:),'linewidth',1.5); hold on;
errorbar(D_yrs,[(nansum(Dfilled(:,1:14))+nansum(Dearly_incomplete)) nansum(Dfilled(:,15:end))].*917./(1000*10^9),...
    [sqrt(nansum((Dfiltered(:,1:14) - Dfiltered_lower(:,1:14)).^2)+nansum(max(abs([Dearly_incompleteMIN-Dearly_incomplete Dearly_incompleteMAX-Dearly_incomplete]),[],2)).^2) Dfilled_err(15:end)].*917./(1000*10^9),...
    [sqrt(nansum((Dfiltered(:,1:14) - Dfiltered_lower(:,1:14)).^2)+nansum(max(abs([Dearly_incompleteMIN-Dearly_incomplete Dearly_incompleteMAX-Dearly_incomplete]),[],2)).^2) Dfilled_err(15:end)].*917./(1000*10^9),'color',cmap(3,:),'linewidth',1.5);
%only glaciers with velocities for every year
pl(4) = plot(D_yrs,nansum(Dfiltered).*917./(1000*10^9),'color',cmap(4,:),'linewidth',1.5); hold on;
errorbar(D_yrs,(nansum(Dfiltered).*917./(1000*10^9)),sqrt(nansum((Dfiltered - Dfiltered_lower).^2)).*917./(1000*10^9),sqrt(nansum((Dfiltered_upper-Dfiltered).^2)).*917./(1000*10^9),'color',cmap(4,:),'linewidth',1.5);
%large glaciers
pl(5) = plot(D_yrs,nansum(Dbig).*917./(1000*10^9),'color',cmap(5,:),'linewidth',1.5); hold on;
errorbar(D_yrs,(nansum(Dbig).*917./(1000*10^9)),sqrt(nansum((Dbig-Dbig_lower).^2)).*917./(1000*10^9),sqrt(nansum((Dbig_upper-Dbig).^2)).*917./(1000*10^9),'color',cmap(5,:),'linewidth',1.5);
legend(pl,['all (n=',num2str(length(Dtotal)),')'],...
    ['filled time gaps (n=',num2str(round(nanmean(sum(~isnan(Dfilled(:,1:14)))))),')'],...
    ['filled time gaps + uniform change (n=',num2str(round(nanmean(sum(~isnan(Dfilled(:,15:end)))))),')'],...
    ['complete timeseries (n=',num2str(sum(~isnan(Dfiltered(:,end)))),')'],...
    ['D>0.05 Gt/yr (n=',num2str(sum(~isnan(Dbig(:,1)))),')'], 'Location', 'southeast','fontsize',16);
% title('Annual GICs Discharge 1985-2018 (Gt/yr)');
ylabel('Discharge (Gt/yr)'); set(gca,'FontSize',20,'xlim',[min(D_yrs) max(D_yrs)],...
    'xtick',[min(D_yrs):3:max(D_yrs)],'xticklabel',[],'ylim',[0 6.5]); grid on;


%add subplot showing total width of glaciers used to calculate discharge
Wtotal = []; coveragetotal = [];
for i = 1:length(term)
    if isempty(find(land_BoxID == term(i).BoxID)) && isempty(find(Mankoff_BoxID == term(i).BoxID))
        %concatenate width if there are velocities for that year
        tempwidth = NaN(size(term(i).fluxVavg));
        tempwidth(~isnan(term(i).fluxVavg)) = sum(term(i).fluxW);
        Wtotal = [Wtotal; tempwidth'];
        
        %concatenate count if there are velocities for that year
        tempcoverage = NaN(size(term(i).fluxVavg));
        tempcoverage(~isnan(term(i).fluxVavg) & term(i).fluxD ~= 0) = 1;
        coveragetotal = [coveragetotal; tempcoverage'];
        clear tempwidth tempcoverage;
    end
end
subb = subplot(2,1,2); ax= gca;
yyaxis left; pf(1) = plot(D_yrs,nansum(Wtotal,1)./1000, '-','color',cmap(1,:),'linewidth',2);
ylabel('Cumulative width (km)','FontSize', 20); ax.YAxis(1).Color = 'k'; ax.YLim = [500 1000];
yyaxis right; pf(2) = plot(D_yrs,100*nansum(coveragetotal,1)./size(coveragetotal,1), '--','color',cmap(1,:),'linewidth',2);
ylabel('Glacier coverage (%)','FontSize', 20); ax.YAxis(2).Color = 'k'; ax.YLim = [50 100];
set(subb,'xlim',[min(D_yrs) max(D_yrs)],'xtick',[min(D_yrs):3:max(D_yrs)],'FontSize', 20); grid on; xlabel('Year'); 
legf = legend(pf,'width','percentage', 'Location', 'southeast','fontsize',16);
set(subb,'position',[0.1300 0.1100 0.7750 0.1800]);
subplot(subb); text(min(D_yrs)+0.5,min(get(gca,'ylim'))+0.90*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'b','FontSize',20);
set(suba,'position',[0.1300 0.3200 0.7750 0.6300]);
subplot(suba); text(min(D_yrs)+0.5,min(get(gca,'ylim'))+0.95*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'a','FontSize',20);

%save the time series
saveas(gcf,[root_dir,'figures/','JOG-21-0146.Figure4.png'],'png'); 
saveas(gcf,[root_dir,'figures/','JOG-21-0146.Figure4.eps'],'epsc');

% %export data to a table
% Dtotal_err(Dtotal==0) = NaN; Dtotal(Dtotal==0) = NaN;
% Dtable = array2table([DBoxID,Dx,Dy,Dtotal]);
% for i = 1:length(D_yrs); Dyrs_cell(i) = {num2str(D_yrs(i))}; end
% Dtable.Properties.VariableNames = [{'BoxID'},{'X'},{'Y'},Dyrs_cell];
% Dtable.Properties.VariableUnits = [{'unitless'},{'m'},{'m'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},...
%     {'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},...
%     {'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},...
%     {'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'}];
% writetable(Dtable,'GreenlandGIC_discharge_timeseries.csv');
% %errors to table
% Dtable = array2table([DBoxID,Dx,Dy,Dtotal_err]);
% for i = 1:length(D_yrs); Dyrs_cell(i) = {num2str(D_yrs(i))}; end
% Dtable.Properties.VariableNames = [{'BoxID'},{'X'},{'Y'},Dyrs_cell];
% Dtable.Properties.VariableUnits = [{'unitless'},{'m'},{'m'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},...
%     {'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},...
%     {'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},...
%     {'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'},{'m^3 yr^-1'}];
% writetable(Dtable,'GreenlandGIC_discharge-uncertainty_timeseries.csv');


%display more info on the change in glacier data availability before & after 1999
fprintf('Percent discharge by largest glaciers before 1999 = %3.1f - %3.1f percent\n',nanmedian(nansum(Dbig(:,1:14))./(nansum(Dfilled(:,1:14))+nansum(Dearly_incomplete)))*100,nanmedian(nansum(Dbig(:,1:14))./nansum(Dfilled(:,1:14)))*100);
fprintf('Percent discharge by largest glaciers from 1999-2018 = %3.1fpercent\n',nanmedian(nansum(Dbig(:,15:end))./nansum(Dfilled(:,15:end)))*100);
fprintf('Percent width of glaciers w/o pre-1999 data = %3.1f percent\n',(nansum(nanmean(Wtotal(complete_flag==0,:),2))./nansum(nanmean(Wtotal,2)))*100);
fprintf('Percent discharge by glaciers w/o pre-1999 data = %3.1f percent\n',nanmean(nansum(Dlate_incomplete,1)./nansum(Dfilled(:,15:end),1))*100);
fprintf('Median percent change in discharge by glaciers w/ complete gap-filled records = %3.1f percent\n',nanmedian((nanmean(Dlate_complete,2)-nanmean(Dearly_complete,2))./nanmean(Dearly_complete,2))*100);
fprintf('MAD percent change in discharge by glaciers w/ complete gap-filled records = %3.1f percent\n',mad((nanmean(Dlate_complete,2)-nanmean(Dearly_complete,2))./nanmean(Dearly_complete,2),1)*100);
fprintf('Median discharge from 1985-1998 with filled gaps = %4.2f Gt/yr\n',nanmedian(nansum(Dfilled(:,1:14))).*917./(1000*10^9))
fprintf('Median discharge from 1999-2018 with filled gaps = %4.2f Gt/yr\n',nanmedian(nansum(Dfilled(:,15:end))).*917./(1000*10^9))
fprintf('Median discharge change from 1985-1998 to 1999-2018 with filled gaps = %4.2f Gt/yr\n',(nanmedian(nansum(Dfilled(:,15:end)))-nanmedian(nansum(Dfilled(:,1:14)))).*917./(1000*10^9))
fprintf('Median discharge from 1985-1998 with filled gaps AND assumed uniform change= %4.2f Gt/yr\n',(nanmedian(nansum(Dfilled(:,1:14)))+nansum(Dearly_incomplete)).*917./(1000*10^9))
fprintf('Median discharge change from 1985-1998 to 1999-2018 with filled gaps AND assumed uniform change = %4.2f Gt/yr\n',(nanmedian(nansum(Dfilled(:,15:end)))-(nanmedian(nansum(Dfilled(:,1:14)))+nansum(Dearly_incomplete))).*917./(1000*10^9))

%% plot regional discharge time series & velocity boxplots
close all;
region_cmap = [253,174,97; 215,25,28; 255,255,191; 171,217,233; 44,123,182]./255;

%set up dummy matrices for time series concatenation
west_normVall = []; west_normVcomplete = []; west_D = []; west_Derr = []; west_Dupper = []; west_Dlower = []; west_term = [];
southeast_normVall = []; southeast_normVcomplete = []; southeast_D = []; southeast_Derr = []; southeast_Dupper = []; southeast_Dlower = []; southeast_term = [];
centraleast_normVall = []; centraleast_normVcomplete = []; centraleast_D = []; centraleast_Derr = []; centraleast_Dupper = []; centraleast_Dlower = []; centraleast_term = [];
northeast_normVall = []; northeast_normVcomplete = []; northeast_D = []; northeast_Derr = []; northeast_Dupper = []; northeast_Dlower = [];northeast_term = [];
north_normVall = []; north_normVcomplete = []; north_D = []; north_Derr = []; north_Dupper = []; north_Dlower = []; north_term = [];

%loop through & parse normalized velocity & discharge time series according
%to region, with separate matrices for all glaciers and only those with
%complete velocity time series
for i = 1:length(term)
    if isempty(find(land_BoxID == term(i).BoxID)) && isempty(find(Mankoff_BoxID == term(i).BoxID))
        %calculate min & max velocities
        maxvel = nanmax(term(i).fluxVavg); minvel = nanmin(term(i).fluxVavg);
        
        %calculate the number of NaNs in the velocity time series
        nancount = sum(isnan(term(i).fluxVavg));
        
        %replace zeros in discharge time series with NaNs
        term(i).fluxD_filled(isnan(term(i).fluxVavg)) = NaN;
        term(i).fluxDerr_filled(isnan(term(i).fluxVavg)) = NaN;
        term(i).fluxDmax_filled(isnan(term(i).fluxVavg)) = NaN;
        term(i).fluxDmin_filled(isnan(term(i).fluxVavg)) = NaN;
        
        if term(i).regionFlag == 1; %west
            west_normVall = [west_normVall; ((term(i).fluxVavg'-minvel)/(maxvel-minvel))];
            if nancount == 0
                west_normVcomplete = [west_normVcomplete; ((term(i).fluxVavg'-minvel)/(maxvel-minvel))];
            end
            west_D = [west_D; term(i).fluxD_filled']; west_Derr = [west_Derr; term(i).fluxDerr_filled'];
            west_Dupper = [west_Dupper; term(i).fluxDmax_filled']; west_Dlower = [west_Dlower; term(i).fluxDmin_filled'];
            for j = 1:length(term(i).X)-1
                if ~isnan(term(i).center_termref(j)) && ~isnan(term(i).center_termref(j+1))
                    termchange(j) = -(term(i).centerline(term(i).center_termref(j+1)) - term(i).centerline(term(i).center_termref(j)));
                else
                    termchange(j) = NaN;
                end
            end
            west_term = [west_term; termchange];
        elseif term(i).regionFlag == 2; %southeast
            southeast_normVall = [southeast_normVall; ((term(i).fluxVavg'-minvel)/(maxvel-minvel))];
            if nancount == 0
                southeast_normVcomplete = [southeast_normVcomplete; ((term(i).fluxVavg'-minvel)/(maxvel-minvel))];
            end
            southeast_D = [southeast_D; term(i).fluxD_filled']; southeast_Derr = [southeast_Derr; term(i).fluxDerr_filled'];
            southeast_Dupper = [southeast_Dupper; term(i).fluxDmax_filled']; southeast_Dlower = [southeast_Dlower; term(i).fluxDmin_filled'];
            for j = 1:length(term(i).X)-1
                if ~isnan(term(i).center_termref(j)) && ~isnan(term(i).center_termref(j+1))
                    termchange(j) = -(term(i).centerline(term(i).center_termref(j+1)) - term(i).centerline(term(i).center_termref(j)));
                else
                    termchange(j) = NaN;
                end
            end
            southeast_term = [southeast_term; termchange];
        elseif term(i).regionFlag == 3; %central east
            centraleast_normVall = [centraleast_normVall; ((term(i).fluxVavg'-minvel)/(maxvel-minvel))];
            if nancount == 0
                centraleast_normVcomplete = [centraleast_normVcomplete; ((term(i).fluxVavg'-minvel)/(maxvel-minvel))];
            end
            centraleast_D = [centraleast_D; term(i).fluxD_filled']; centraleast_Derr = [centraleast_Derr; term(i).fluxDerr_filled'];
            centraleast_Dupper = [centraleast_Dupper; term(i).fluxDmax_filled']; centraleast_Dlower = [centraleast_Dlower; term(i).fluxDmin_filled'];
            for j = 1:length(term(i).X)-1
                if ~isnan(term(i).center_termref(j)) && ~isnan(term(i).center_termref(j+1))
                    termchange(j) = -(term(i).centerline(term(i).center_termref(j+1)) - term(i).centerline(term(i).center_termref(j)));
                else
                    termchange(j) = NaN;
                end
            end
            centraleast_term = [centraleast_term; termchange];
        elseif term(i).regionFlag == 4; %northeast
            northeast_normVall = [northeast_normVall; ((term(i).fluxVavg'-minvel)/(maxvel-minvel))];
            if nancount == 0
                northeast_normVcomplete = [northeast_normVcomplete; ((term(i).fluxVavg'-minvel)/(maxvel-minvel))];
            end
            northeast_D = [northeast_D; term(i).fluxD_filled']; northeast_Derr = [northeast_Derr; term(i).fluxDerr_filled'];
            northeast_Dupper = [northeast_Dupper; term(i).fluxDmax_filled']; northeast_Dlower = [northeast_Dlower; term(i).fluxDmin_filled'];
            for j = 1:length(term(i).X)-1
                if ~isnan(term(i).center_termref(j)) && ~isnan(term(i).center_termref(j+1))
                    termchange(j) = -(term(i).centerline(term(i).center_termref(j+1)) - term(i).centerline(term(i).center_termref(j)));
                else
                    termchange(j) = NaN;
                end
            end
            northeast_term = [northeast_term; termchange];
        else
            north_normVall = [north_normVall; ((term(i).fluxVavg'-minvel)/(maxvel-minvel))];
            if nancount == 0
                north_normVcomplete = [north_normVcomplete; ((term(i).fluxVavg'-minvel)/(maxvel-minvel))];
            end
            north_D = [north_D; term(i).fluxD_filled']; north_Derr = [north_Derr; term(i).fluxDerr_filled'];
            north_Dupper = [north_Dupper; term(i).fluxDmax_filled']; north_Dlower = [north_Dlower; term(i).fluxDmin_filled'];
            for j = 1:length(term(i).X)-1
                if ~isnan(term(i).center_termref(j)) && ~isnan(term(i).center_termref(j+1))
                    termchange(j) = -(term(i).centerline(term(i).center_termref(j+1)) - term(i).centerline(term(i).center_termref(j)));
                else
                    termchange(j) = NaN;
                end
            end
            north_term = [north_term; termchange];
        end
        
        clear maxvel minvel nancount termchange;
    end
end
disp('Regional discharges:');
disp('west');
disp(['1985-1998 median +/ MAD = ',num2str(nanmedian(nansum(west_D(:,1:14),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(west_D(:,1:14),1),1).*917./(1000*10^9)),' Gt/yr']);
disp(['1999-2018 median +/- MAD = ',num2str(nanmedian(nansum(west_D(:,15:end),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(west_D(:,15:end),1),1).*917./(1000*10^9)),' Gt/yr']);
disp(['percentage of glaciers = ',num2str(round(100*size(west_normVall,1)./size(Dtotal,1),1)),'%']);
disp(['percentage discharge: 1985-1998 = ',num2str(round(100*nanmedian(nansum(west_D(:,1:14),1))./nanmedian(nansum(Dtotal(:,1:14),1)),1)),'% & 1999-218 = ',...
    num2str(round(100*nanmedian(nansum(west_D(:,15:end),1))./nanmedian(nansum(Dtotal(:,15:end),1)),1)),'%']);
disp('southeast');
disp(['1985-1998 median +/ MAD = ',num2str(nanmedian(nansum(southeast_D(:,1:14),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(southeast_D(:,1:14),1),1).*917./(1000*10^9)),' Gt/yr']);
disp(['1999-2018 median +/- MAD = ',num2str(nanmedian(nansum(southeast_D(:,15:end),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(southeast_D(:,15:end),1),1).*917./(1000*10^9)),' Gt/yr']);
disp(['percentage of glaciers = ',num2str(round(100*size(southeast_normVall,1)./size(Dtotal,1),1)),'%']);
disp(['percentage discharge: 1985-1998 = ',num2str(round(100*nanmedian(nansum(southeast_D(:,1:14),1))./nanmedian(nansum(Dtotal(:,1:14),1)),1)),'% & 1999-218 = ',...
    num2str(round(100*nanmedian(nansum(southeast_D(:,15:end),1))./nanmedian(nansum(Dtotal(:,15:end),1)),1)),'%']);
disp('central east');
disp(['1985-1998 median +/ MAD = ',num2str(nanmedian(nansum(centraleast_D(:,1:14),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(centraleast_D(:,1:14),1),1).*917./(1000*10^9)),' Gt/yr']);
disp(['1999-2018 median +/- MAD = ',num2str(nanmedian(nansum(centraleast_D(:,15:end),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(centraleast_D(:,15:end),1),1).*917./(1000*10^9)),' Gt/yr']);
disp(['percentage of glaciers = ',num2str(round(100*size(centraleast_normVall,1)./size(Dtotal,1),1)),'%']);
disp(['percentage discharge: 1985-1998 = ',num2str(round(100*nanmedian(nansum(centraleast_D(:,1:14),1))./nanmedian(nansum(Dtotal(:,1:14),1)),1)),'% & 1999-218 = ',...
    num2str(round(100*nanmedian(nansum(centraleast_D(:,15:end),1))./nanmedian(nansum(Dtotal(:,15:end),1)),1)),'%']);
disp('northeast');
disp(['1985-1998 median +/ MAD = ',num2str(nanmedian(nansum(northeast_D(:,1:14),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(northeast_D(:,1:14),1),1).*917./(1000*10^9)),' Gt/yr']);
disp(['1999-2018 median +/- MAD = ',num2str(nanmedian(nansum(northeast_D(:,15:end),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(northeast_D(:,15:end),1),1).*917./(1000*10^9)),' Gt/yr']);
disp(['percentage of glaciers = ',num2str(round(100*size(northeast_normVall,1)./size(Dtotal,1),1)),'%']);
disp(['percentage discharge: 1985-1998 = ',num2str(round(100*nanmedian(nansum(northeast_D(:,1:14),1))./nanmedian(nansum(Dtotal(:,1:14),1)),1)),'% & 1999-218 = ',...
    num2str(round(100*nanmedian(nansum(northeast_D(:,15:end),1))./nanmedian(nansum(Dtotal(:,15:end),1)),1)),'%']);
disp('north');
disp(['1985-1998 median +/ MAD = ',num2str(nanmedian(nansum(north_D(:,1:14),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(north_D(:,1:14),1),1).*917./(1000*10^9)),' Gt/yr']);
disp(['1999-2018 median +/- MAD = ',num2str(nanmedian(nansum(north_D(:,15:end),1)).*917./(1000*10^9)),' +/- ',num2str(mad(nansum(north_D(:,15:end),1),1).*917./(1000*10^9)),' Gt/yr']);
disp(['percentage of glaciers = ',num2str(round(100*size(north_normVall,1)./size(Dtotal,1),1)),'%']);
disp(['percentage discharge: 1985-1998 = ',num2str(round(100*nanmedian(nansum(north_D(:,1:14),1))./nanmedian(nansum(Dtotal(:,1:14),1)),1)),'% & 1999-218 = ',...
    num2str(round(100*nanmedian(nansum(north_D(:,15:end),1))./nanmedian(nansum(Dtotal(:,15:end),1)),1)),'%']);

%create regional subplots
allV_figure = figure; set(gcf,'position',[50 50 800 1200]);
set(allV_figure, 'DefaultTextFontSize', 20);
%add boxplots for normalized velocity for each region
westVsub = subplot(5,3,3);
boxplot(west_normVall, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors',region_cmap(1,:), 'Widths', 0.6, 'LabelOrientation', 'inline','Symbol','','Whisker',0); hold on;
boxplot(west_normVall, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors','k', 'Widths', 0.25, 'LabelOrientation', 'inline','BoxStyle','outline'); hold on;
set(gca,'XTickLabel',{' '}); pos = get(gca,'position'); set(gca,'position',[0.65 pos(2) 0.32 1.25*pos(4)]);
set(gca,'fontsize',20); pos1 = get(gca,'position');
annotation('textbox',[0.65 0.7775 0.3 0.05],'string',['e '],'fontsize',16,'edgecolor','none')
annotation('textbox',[0.77 0.7775 0.3 0.05],'string',['(WE=',num2str(size(west_normVall,1)),')'],'fontsize',16,'edgecolor','none')
southeastVsub = subplot(5,3,6);
boxplot(southeast_normVall, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors',region_cmap(2,:), 'Widths', 0.6, 'LabelOrientation', 'inline','Symbol','','Whisker',0); hold on;
boxplot(southeast_normVall, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors','k', 'Widths', 0.25, 'LabelOrientation', 'inline','BoxStyle','outline'); hold on;
set(gca,'XTickLabel',{' '}); pos = get(gca,'position');  set(gca,'position',[0.65 pos(2)-0.005 0.32 1.25*pos(4)]);
set(gca,'fontsize',20);
annotation('textbox',[0.65 0.605 0.3 0.05],'string',['f '],'fontsize',16,'edgecolor','none')
annotation('textbox',[0.77 0.605 0.3 0.05],'string',['(SE=',num2str(size(southeast_normVall,1)),')'],'fontsize',16,'edgecolor','none')
centraleastVsub = subplot(5,3,9);
boxplot(centraleast_normVall, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors',region_cmap(3,:), 'Widths', 0.6, 'LabelOrientation', 'inline','Symbol','','Whisker',0); hold on;
boxplot(centraleast_normVall, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors','k', 'Widths', 0.25, 'LabelOrientation', 'inline','BoxStyle','outline'); hold on;
set(gca,'XTickLabel',{' '}); pos = get(gca,'position');  set(gca,'position',[0.65 pos(2)-0.01 0.32 1.25*pos(4)]);
set(gca,'fontsize',20);
annotation('textbox',[0.65 0.4275 0.3 0.05],'string',['g '],'fontsize',16,'edgecolor','none')
annotation('textbox',[0.77 0.4275 0.3 0.05],'string',['(CE=',num2str(size(centraleast_normVall,1)),')'],'fontsize',16,'edgecolor','none')
northeastVsub = subplot(5,3,12);
boxplot(northeast_normVall, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors',region_cmap(4,:), 'Widths', 0.6, 'LabelOrientation', 'inline','Symbol','','Whisker',0); hold on;
boxplot(northeast_normVall, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors','k', 'Widths', 0.25, 'LabelOrientation', 'inline','BoxStyle','outline'); hold on;
set(gca,'XTickLabel',{' '}); pos = get(gca,'position');  set(gca,'position',[0.65 pos(2)-0.015 0.32 1.25*pos(4)]);
set(gca,'fontsize',20); pos4 = get(gca,'position');
annotation('textbox',[0.65 0.25 0.3 0.05],'string',['h '],'fontsize',16,'edgecolor','none')
annotation('textbox',[0.77 0.25 0.3 0.05],'string',['(NE=',num2str(size(northeast_normVall,1)),')'],'fontsize',16,'edgecolor','none')
northVsub = subplot(5,3,15);
boxplot(north_normVall, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors',region_cmap(5,:), 'Widths', 0.6, 'LabelOrientation', 'inline','Symbol','','Whisker',0); hold on;
boxplot(north_normVall, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors','k', 'Widths', 0.25, 'LabelOrientation', 'inline','BoxStyle','outline',...
    'Labels',{num2str(min(D_yrs)),'','',num2str(min(D_yrs)+3),'','',num2str(min(D_yrs)+6),'','',num2str(min(D_yrs)+9),...
    '','',num2str(min(D_yrs)+12),'','',num2str(min(D_yrs)+15),'','',num2str(min(D_yrs)+18),...
    '','',num2str(min(D_yrs)+21),'','',num2str(min(D_yrs)+24),'','',num2str(min(D_yrs)+27),'','',...
    num2str(min(D_yrs)+30),'','',num2str(min(D_yrs)+33)}); hold on;
% set(gca,'XTickLabel',{' '})
pos = get(gca,'position'); set(gca,'position',[0.65 pos(2)-0.02 0.32 1.25*pos(4)]);
set(gca,'fontsize',20); pos5 = get(gca,'position');
xlabel('Year','fontsize',20); ylbl = ylabel('Normalized speed','fontsize',20);
set(ylbl,'position',[-4.4141 3.25 -1]);
annotation('textbox',[0.65 0.065 0.3 0.05],'string',['i '],'fontsize',16,'edgecolor','none')
annotation('textbox',[0.77 0.065 0.3 0.05],'string',['(NO=',num2str(size(north_normVall,1)),')'],'fontsize',16,'edgecolor','none')
% boxplot(west_normVcomplete, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors', [0.5 0.5 0.5], 'Widths', 0.4,'Symbol','x','LabelOrientation', 'inline'); hold on;
%add scatterplots for discharge for each region
Dsub = subplot(5,3,[1:2 4:5 7:8]);
errorbar(D_yrs,nansum(west_D).*917./(1000*10^9),sqrt(nansum((west_D-west_Dlower).^2)).*917./(1000*10^9),sqrt(nansum((west_Dupper-west_D).^2)).*917./(1000*10^9),'color',region_cmap(1,:),'linewidth',1.5); hold on;
vp(1) = plot(D_yrs,nansum(west_D).*917./(1000*10^9),'-','color',region_cmap(1,:),'linewidth',2); hold on;
errorbar(D_yrs,nansum(southeast_D).*917./(1000*10^9),sqrt(nansum((southeast_D-southeast_Dlower).^2)).*917./(1000*10^9),sqrt(nansum((southeast_Dupper-southeast_D).^2)).*917./(1000*10^9),'color',region_cmap(2,:),'linewidth',1.5);
vp(2) = plot(D_yrs,nansum(southeast_D).*917./(1000*10^9),'-','color',region_cmap(2,:),'linewidth',2); hold on;
errorbar(D_yrs,nansum(centraleast_D).*917./(1000*10^9),sqrt(nansum((centraleast_D-centraleast_Dlower).^2)).*917./(1000*10^9),sqrt(nansum((centraleast_Dupper-centraleast_D).^2)).*917./(1000*10^9),'color','k','linewidth',1.5);
errorbar(D_yrs,nansum(centraleast_D).*917./(1000*10^9),sqrt(nansum((centraleast_D-centraleast_Dlower).^2)).*917./(1000*10^9),sqrt(nansum((centraleast_Dupper-centraleast_D).^2)).*917./(1000*10^9),'color',region_cmap(3,:),'linewidth',1);
plot(D_yrs,nansum(centraleast_D).*917./(1000*10^9),'-','color','k','linewidth',2); hold on;
vp(3) = plot(D_yrs,nansum(centraleast_D).*917./(1000*10^9),'-','color',region_cmap(3,:),'linewidth',1.5); hold on;
errorbar(D_yrs,nansum(northeast_D).*917./(1000*10^9),sqrt(nansum((northeast_D-northeast_Dlower).^2)).*917./(1000*10^9),sqrt(nansum((northeast_Dupper-northeast_D).^2)).*917./(1000*10^9),'color',region_cmap(4,:),'linewidth',1.5);
vp(4) = plot(D_yrs,nansum(northeast_D).*917./(1000*10^9),'-','color',region_cmap(4,:),'linewidth',2); hold on;
errorbar(D_yrs,nansum(north_D).*917./(1000*10^9),sqrt(nansum((north_D-north_Dlower).^2)).*917./(1000*10^9),sqrt(nansum((north_Dupper-north_D).^2)).*917./(1000*10^9),'color',region_cmap(5,:),'linewidth',1.5);
vp(5) = plot(D_yrs,nansum(north_D).*917./(1000*10^9),'-','color',region_cmap(5,:),'linewidth',2); hold on;
set(gca,'ylim',[0 4],'xlim',[min(D_yrs)-0.5 max(D_yrs)+0.5],'xtick',[min(D_yrs):3:max(D_yrs)],...
    'xticklabel',[]); grid on;
pos = get(gca,'position'); set(gca,'position',[0.09 pos4(2)+pos4(4)+0.04 0.47 (pos1(2)+pos1(4))-(pos4(2)+pos4(4)+0.04)]); 
set(gca,'fontsize',20); ylabel('Discharge (Gt yr^{-1})','fontsize',20);
text(min(D_yrs)+0.5,3.85,'a ','fontsize',16);
vleg = legend(vp,'west (WE)','southeast (SE)','central east (CE)','northeast (NE)','north (NO)'); set(vleg,'location','northeast');
Dsub_perc = subplot(5,3,[10:11]);
plot(D_yrs,100*sum(~isnan(west_normVall))./size(west_normVall,1),'-','color',region_cmap(1,:),'linewidth',2); hold on;
plot(D_yrs,100*sum(~isnan(southeast_normVall))./size(southeast_normVall,1),'-','color',region_cmap(2,:),'linewidth',2); hold on;
plot(D_yrs,100*sum(~isnan(centraleast_normVall))./size(centraleast_normVall,1),'-','color','k','linewidth',2); hold on;
plot(D_yrs,100*sum(~isnan(centraleast_normVall))./size(centraleast_normVall,1),'-','color',region_cmap(3,:),'linewidth',1.5); hold on;
plot(D_yrs,100*sum(~isnan(northeast_normVall))./size(northeast_normVall,1),'-','color',region_cmap(4,:),'linewidth',2); hold on;
plot(D_yrs,100*sum(~isnan(north_normVall))./size(north_normVall,1),'-','color',region_cmap(5,:),'linewidth',2); hold on;
set(gca,'ylim',[0 110],'ytick',[0:25:100],'xlim',[min(D_yrs)-0.5 max(D_yrs)+0.5],'xtick',[min(D_yrs):3:max(D_yrs)],...
    'xticklabel',[min(D_yrs):3:max(D_yrs)],'xticklabelrotation',90); grid on;
pos = get(gca,'position'); set(gca,'position',[0.09 pos4(2)+0.045 0.47 pos4(4)-0.01]); 
set(gca,'fontsize',20); xlabel('Year','fontsize',20); ylabel('Glacier coverage (%)','fontsize',20);
text(min(D_yrs)+0.5,95,'b ','fontsize',16);
termsub1 = subplot(5,3,13);
he(1) = histogram(southeast_term(:,1),[-2500:100:500]); he(1).FaceColor = region_cmap(2,:); he(1).FaceAlpha = 1; he(1).EdgeColor = 'k'; hold on;
he(2) = histogram(west_term(:,1),[-2500:100:500]); he(2).FaceColor = region_cmap(1,:); he(2).EdgeColor = 'k'; hold on;
he(3) = histogram(centraleast_term(:,1),[-2500:100:500]); he(3).FaceColor = region_cmap(3,:); he(3).FaceAlpha = 0.75; he(3).EdgeColor = 'k'; hold on;
he(4) = histogram(north_term(:,1),[-2500:100:500]); he(4).FaceColor = region_cmap(5,:); he(4).EdgeColor = 'k'; hold on;
he(5) = histogram(northeast_term(:,1),[-2500:100:500]); he(5).FaceColor = region_cmap(4,:); he(5).FaceAlpha = 1; he(5).EdgeColor = 'k'; hold on;
set(gca,'ylim',[0 75],'xtick',[-3000:1000:1000],...
    'xticklabel',[-3:1:1],'xticklabelrotation',90); set(gca,'fontsize',20); 
ylabel('Glacier count','fontsize',20);
termpos = get(gca,'position'); 
text(-2475,65,'c ','fontsize',16);
termsub2 = subplot(5,3,14);
hl(1) = histogram(southeast_term(:,2),[-2500:100:500]); hl(1).FaceColor = region_cmap(2,:); hl(1).FaceAlpha = 1; hl(1).EdgeColor = 'k'; hold on;
hl(2) = histogram(west_term(:,2),[-2500:100:500]); hl(2).FaceColor = region_cmap(1,:); hl(2).EdgeColor = 'k'; hold on;
hl(3) = histogram(centraleast_term(:,2),[-2500:100:500]); hl(3).FaceColor = region_cmap(3,:); h1(3).FaceAlpha = 0.75; hl(3).EdgeColor = 'k'; hold on;
hl(4) = histogram(north_term(:,2),[-2500:100:500]); hl(4).FaceColor = region_cmap(5,:); hl(4).EdgeColor = 'k'; hold on;
hl(5) = histogram(northeast_term(:,2),[-2500:100:500]); hl(5).FaceColor = region_cmap(4,:); hl(5).FaceAlpha = 1; hl(5).EdgeColor = 'k'; hold on;
set(gca,'ylim',[0 75],'yticklabel',[],'xtick',[-3000:1000:1000],...
    'xticklabel',[-3:1:1],'xticklabelrotation',90); set(gca,'fontsize',20); 
xlbl = xlabel('Terminus change (km)','fontsize',20); set(xlbl,'position',[-2750 -12 -1]);
text(-2475,65,'d ','fontsize',16);
%shift to align with year span in plots above
set(termsub2,'position',[0.30 pos5(2) 0.21 pos5(4)-0.01]);
set(termsub1,'position',[0.09 pos5(2) 0.21 pos5(4)-0.01]);
%add time period labels
annotation('line',[0.09 0.09],[0.2 0.255]); annotation('line',[0.30 0.30],[0.2 0.255]); annotation('line',[0.51 0.51],[0.2 0.255]);
annotation('textbox',[0.13 0.205 0.3 0.05],'string',['1985-2000'],'fontsize',20,'edgecolor','none')
annotation('doublearrow',[.09 .30],[.255 .255])
annotation('textbox',[0.35 0.205 0.3 0.05],'string',['2000-2015'],'fontsize',20,'edgecolor','none')
annotation('doublearrow',[.30 0.51],[.255 .255])

%save the regional plots
saveas(gcf,[root_dir,'figures/','JOG-21-0146.Figure5.png'],'png'); 
saveas(gcf,[root_dir,'figures/','JOG-21-0146.Figure5.eps'],'epsc');


%create a separate figure showing the median regional speed changes for
%only glaciers with complete velocity timeseries and all glaciers then the
%boxplots for only complete glaciers below to highlight patterns
completeV_figure = figure; set(gcf,'position',[50 50 400 800]);
set(completeV_figure, 'DefaultTextFontSize', 20);
%add boxplots for normalized velocity for each region
westVsub = subplot(4,1,1);
boxplot(west_normVcomplete, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors',region_cmap(1,:), 'Widths', 0.6, 'LabelOrientation', 'inline','Symbol','','Whisker',0); hold on;
boxplot(west_normVcomplete, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors','k', 'Widths', 0.25, 'LabelOrientation', 'inline','BoxStyle','outline'); hold on;
set(gca,'XTickLabel',{' '}); pos = get(gca,'position'); set(gca,'position',[0.15 pos(2) pos(3) 1.25*pos(4)]);
set(gca,'fontsize',20); pos1 = get(gca,'position');
annotation('textbox',[0.15 0.7575 0.7 0.05],'string',['a '],'fontsize',16,'edgecolor','none')
annotation('textbox',[0.45 0.7575 0.7 0.05],'string',['(WE=',num2str(size(west_normVcomplete,1)),')'],'fontsize',16,'edgecolor','none')
southeastVsub = subplot(4,1,2);
boxplot(southeast_normVcomplete, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors',region_cmap(2,:), 'Widths', 0.6, 'LabelOrientation', 'inline','Symbol','','Whisker',0); hold on;
boxplot(southeast_normVcomplete, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors','k', 'Widths', 0.25, 'LabelOrientation', 'inline','BoxStyle','outline'); hold on;
set(gca,'XTickLabel',{' '}); pos = get(gca,'position');  set(gca,'position',[0.15 pos(2) pos(3) 1.25*pos(4)]);
set(gca,'fontsize',20);
annotation('textbox',[0.15 0.5425 0.7 0.05],'string',['b '],'fontsize',16,'edgecolor','none')
annotation('textbox',[0.45 0.5425 0.7 0.05],'string',['(SE=',num2str(size(southeast_normVcomplete,1)),')'],'fontsize',16,'edgecolor','none')
centraleastVsub = subplot(4,1,3);
boxplot(centraleast_normVcomplete, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors',region_cmap(3,:), 'Widths', 0.6, 'LabelOrientation', 'inline','Symbol','','Whisker',0); hold on;
boxplot(centraleast_normVcomplete, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors','k', 'Widths', 0.25, 'LabelOrientation', 'inline','BoxStyle','outline'); hold on;
set(gca,'XTickLabel',{' '}); pos = get(gca,'position');  set(gca,'position',[0.15 pos(2) pos(3) 1.25*pos(4)]);
set(gca,'fontsize',20);
annotation('textbox',[0.15 0.325 0.7 0.05],'string',['c '],'fontsize',16,'edgecolor','none')
annotation('textbox',[0.45 0.325 0.7 0.05],'string',['(CE=',num2str(size(centraleast_normVcomplete,1)),')'],'fontsize',16,'edgecolor','none')
northeastVsub = subplot(4,1,4);
boxplot(northeast_normVcomplete, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors',region_cmap(4,:), 'Widths', 0.6, 'LabelOrientation', 'inline','Symbol','','Whisker',0); hold on;
boxplot(northeast_normVcomplete, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors','k', 'Widths', 0.25, 'LabelOrientation', 'inline','BoxStyle','outline',...
    'Labels',{num2str(min(D_yrs)),'','',num2str(min(D_yrs)+3),'','',num2str(min(D_yrs)+6),'','',num2str(min(D_yrs)+9),...
    '','',num2str(min(D_yrs)+12),'','',num2str(min(D_yrs)+15),'','',num2str(min(D_yrs)+18),...
    '','',num2str(min(D_yrs)+21),'','',num2str(min(D_yrs)+24),'','',num2str(min(D_yrs)+27),'','',...
    num2str(min(D_yrs)+30),'','',num2str(min(D_yrs)+33)}); hold on;
% set(gca,'XTickLabel',{' '})
pos = get(gca,'position');  set(gca,'position',[0.15 pos(2) pos(3) 1.25*pos(4)]);
set(gca,'fontsize',20); pos4 = get(gca,'position');
annotation('textbox',[0.15 0.0975 0.7 0.05],'string',['d '],'fontsize',16,'edgecolor','none')
annotation('textbox',[0.45 0.0975 0.7 0.05],'string',['(NE=',num2str(size(northeast_normVcomplete,1)),')'],'fontsize',16,'edgecolor','none')
% northVsub = subplot(5,1,5);
% boxplot(north_normVcomplete, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors',region_cmap(5,:), 'Widths', 0.6, 'LabelOrientation', 'inline','Symbol','','Whisker',0); hold on;
% boxplot(north_normVcomplete, D_yrs, 'PlotStyle', 'compact', 'MedianStyle', 'line', 'Colors','k', 'Widths', 0.25, 'LabelOrientation', 'inline','BoxStyle','outline',...
%     'Labels',{num2str(min(D_yrs)),'','',num2str(min(D_yrs)+3),'','',num2str(min(D_yrs)+6),'','',num2str(min(D_yrs)+9),...
%     '','',num2str(min(D_yrs)+12),'','',num2str(min(D_yrs)+15),'','',num2str(min(D_yrs)+18),...
%     '','',num2str(min(D_yrs)+21),'','',num2str(min(D_yrs)+24),'','',num2str(min(D_yrs)+27),'','',...
%     num2str(min(D_yrs)+30),'','',num2str(min(D_yrs)+33)}); hold on;
% % set(gca,'XTickLabel',{' '})
% pos = get(gca,'position'); set(gca,'position',[0.13 pos(2) pos(3) 1.25*pos(4)]);
% set(gca,'fontsize',20); pos5 = get(gca,'position');
xlabel('Year','fontsize',20); ylbl = ylabel('Normalized speed','fontsize',20);
set(ylbl,'position',[-3.1 2.5 -1]);
% annotation('textbox',[0.13 0.0875 0.7 0.05],'string',['e) north speeds (n=',num2str(size(north_normVcomplete,1)),')'],'fontsize',16,'edgecolor','none')
saveas(gcf,[root_dir,'figures/','JOG-21-0146.FigureS2.png'],'png'); 
saveas(gcf,[root_dir,'figures/','JOG-21-0146.FigureS2.eps'],'epsc');


%% create plots of terminus change vs speed change
close all;
term_figure = figure; set(gcf,'position',[50 50 800 400]);
set(term_figure, 'DefaultTextFontSize', 20);
sub_early = subplot(1,2,1); sub_late = subplot(1,2,2);
speedyrs_early = [1985:1:2000]; speedyrs_late = [2000:1:2015];
%plot rate of terminus change vs speed change
GIC_dterm = NaN(length(term),2); GIC_dtermdyr = NaN(length(term),2);
GIC_termyrs = NaN(length(term),3);
for i = 1:length(term)
    if isempty(find(land_BoxID == term(i).BoxID)) && isempty(find(Mankoff_BoxID == term(i).BoxID))
        GIC_termyrs(i,:) = term(i).year;
        for j = 1:length(term(i).X)-1
            if ~isnan(term(i).center_termref(j)) && ~isnan(term(i).center_termref(j+1))
                termchange(j) = -((term(i).centerline(term(i).center_termref(j+1)) - term(i).centerline(term(i).center_termref(j)))./(term(i).year(j+1)-term(i).year(j)));
                GIC_dterm(i,j) = termchange(j).*(term(i).year(j+1)-term(i).year(j)); GIC_dtermdyr(i,j) = termchange(j);
            else
                termchange(j) = NaN;
            end
        end
        speedchange_early = term(i).fluxVavg(16) - term(i).fluxVavg(1:15);
        speedchange_late = term(i).fluxVavg(31) - term(i).fluxVavg(16:30);
        
        %plot
        if sum(~isnan(speedchange_early)) > 0
            subplot(sub_early); plot(termchange(1),speedchange_early(find(~isnan(speedchange_early),1,'first'))./(speedyrs_early(end) - speedyrs_early(find(~isnan(speedchange_early),1,'first'))),'dk','markerfacecolor',region_cmap(term(i).regionFlag,:)); hold on;
        end
        if sum(~isnan(speedchange_late)) > 0
            subplot(sub_late); plot(termchange(2),speedchange_late(find(~isnan(speedchange_late),1,'first'))./(speedyrs_late(end) - speedyrs_late(find(~isnan(speedchange_late),1,'first'))),'dk','markerfacecolor',region_cmap(term(i).regionFlag,:)); hold on;
        end
        clear termchange speedchange*;

    end
    drawnow;
end
subplot(sub_early); set(gca,'fontsize',20,'xlim',[-120 48],...
    'ylim',[-36 12],'ytick',[-36:12:12]); grid on;
xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
text(min(xlims)+0.05*(max(xlims)-min(xlims)),min(ylims)+0.05*(max(ylims)-min(ylims)),'a) ~1985-2000');
xlabel('Terminus change rate (m/yr)','fontsize',20); ylabel('Speed change rate (m/yr)','fontsize',20);
subplot(sub_late); set(gca,'fontsize',20,'xlim',[-120 48],...
    'ylim',[-36 12],'ytick',[-36:12:12],'yticklabel',[]); grid on;
text(min(xlims)+0.05*(max(xlims)-min(xlims)),min(ylims)+0.05*(max(ylims)-min(ylims)),'b) 2000-2015');
set(sub_early,'position',[0.13 0.15 0.4 0.7875]); set(sub_late,'position',[0.55 0.15 0.4 0.7875]); 
xlabel('Terminus change rate (m/yr)','fontsize',20); %ylabel('Speed change rate (m/yr)','fontsize',20);
saveas(gcf,[root_dir,'figures/','GreenlandGIC-speed-vs-termchange_subplots.png'],'png'); 
saveas(gcf,[root_dir,'figures/','GreenlandGIC-speed-vs-termchange_subplots.eps'],'epsc');

%create histograms of terminus delineation years
termyr_figure = figure; set(gcf,'position',[450 50 800 400]);
h1 = histogram(GIC_termyrs(:,1),'FaceColor',[123,50,148]/255); hold on;
h2 = histogram(GIC_termyrs(:,2),'FaceColor','k'); hold on;
h3 = histogram(GIC_termyrs(:,3),'FaceColor',[0,136,55]/255); hold on;
set(gca,'fontsize',20); grid on;
legend([h1 h2 h3],'Landsat 5','Landsat 7','Landsat 8','location','northwest');
xlabel('Year ','fontsize',20); ylabel('Observation count','fontsize',20);
saveas(gcf,[root_dir,'figures/','JOG-21-0146.Figure2.png'],'png'); 
saveas(gcf,[root_dir,'figures/','JOG-21-0146.Figure2.eps'],'epsc');

%spit out stats
disp('Terminus change magnitude:');
disp(['median: 1985-2000 = ',num2str(nanmedian(GIC_dterm(:,1))),' +/- ',num2str(mad(GIC_dterm(:,1),1)),...
    ' m & 2000-2015 = ',num2str(nanmedian(GIC_dterm(:,2))),' +/- ',num2str(mad(GIC_dterm(:,2),1)),' m']);
disp('Terminus change rate:');
disp(['median: 1985-2000 = ',num2str(nanmedian(GIC_dtermdyr(:,1))),' +/- ',num2str(mad(GIC_dtermdyr(:,1),1)),...
    ' m/yr & 2000-2015 = ',num2str(nanmedian(GIC_dtermdyr(:,2))),' +/- ',num2str(mad(GIC_dtermdyr(:,2),1)),' m/yr']);

%% compare regional speed time series to summer & winter NAO indices
disp('Must run the initialization section & the section immediately before this one');

%load the NAO data
T = readtable('NAO-index_1950-2021.csv');
NAOyr = table2array(T(:,1)); NAOmo = table2array(T(:,2));
NAOindex = table2array(T(:,3)); clear T;

%create a vector of unique years
NAOyrs = unique(NAOyr);

%create vectors of mean annual NAO indices for winter (Oct-Apr) & summer (Jun-Sept)
winter_refs = find((NAOmo >=1 & NAOmo <=4) | (NAOmo >=10 & NAOmo <=12));
summer_refs = find(NAOmo >=5 & NAOmo <=9);
NAOyrs_winter = NAOyr(winter_refs); NAOindices_winter = NAOindex(winter_refs);
NAOyrs_summer = NAOyr(summer_refs); NAOindices_summer = NAOindex(summer_refs);
%winter
[yrs,yr_refs] = unique(NAOyrs_winter);
for j = 1:length(yrs)
    if j < length(yrs)
        NAOindex_winter(j) = nanmean(NAOindices_winter(yr_refs(j):yr_refs(j+1)-1));
    else
        NAOindex_winter(j) = nanmean(NAOindices_winter(yr_refs(j):end));
    end
end
clear yrs yr_refs;
%summer
[yrs,yr_refs] = unique(NAOyrs_summer);
for j = 1:length(yrs)
    if j < length(yrs)
        NAOindex_summer(j) = nanmean(NAOindices_summer(yr_refs(j):yr_refs(j+1)-1));
    else
        NAOindex_summer(j) = nanmean(NAOindices_summer(yr_refs(j):end));
    end
end
clear yrs yr_refs;

%calculate moving means (2yrs, 5yrs, 10yrs)
for j = 1:length(D_yrs)
    %winter
    NAOindex_2yrwinter(j) = nanmean(NAOindex_winter(find(NAOyrs == D_yrs(j)-1):find(NAOyrs == D_yrs(j))));
    NAOindex_5yrwinter(j) = nanmean(NAOindex_winter(find(NAOyrs == D_yrs(j)-4):find(NAOyrs == D_yrs(j))));
    NAOindex_10yrwinter(j) = nanmean(NAOindex_winter(find(NAOyrs == D_yrs(j)-9):find(NAOyrs == D_yrs(j))));
    
    %summer
    NAOindex_2yrsummer(j) = nanmean(NAOindex_summer(find(NAOyrs == D_yrs(j)-1):find(NAOyrs == D_yrs(j))));
    NAOindex_5yrsummer(j) = nanmean(NAOindex_summer(find(NAOyrs == D_yrs(j)-4):find(NAOyrs == D_yrs(j))));
    NAOindex_10yrsummer(j) = nanmean(NAOindex_summer(find(NAOyrs == D_yrs(j)-9):find(NAOyrs == D_yrs(j))));
    
end
%plot the NAO indices
figure; 
%winter
plot(D_yrs,NAOindex_2yrwinter,'-b','linewidth',2); hold on;
plot(D_yrs,NAOindex_5yrwinter,'--b','linewidth',2); hold on;
plot(D_yrs,NAOindex_10yrwinter,'-.b','linewidth',2); hold on;
%summer
plot(D_yrs,NAOindex_2yrsummer,'-r','linewidth',2); hold on;
plot(D_yrs,NAOindex_5yrsummer,'--r','linewidth',2); hold on;
plot(D_yrs,NAOindex_10yrsummer,'-.r','linewidth',2); hold on;
set(gca,'fontsize',20,'xlim',[min(D_yrs),2020]); grid on;
leg = legend('winter: 2yr','winter: 5yr','winter: 10yr','summer: 2yr','summer: 5yr','summer: 10yr','NumColumns',2);
xlabel('Year ','fontsize',20); ylabel('NAO index ','fontsize',20);

%determine correlation between regional normalized speeds & NAO indices
%WEST
westRwinter = corrcoef(nanmedian(west_normVcomplete)',NAOindex_5yrwinter');
westRsummer = corrcoef(nanmedian(west_normVcomplete)',NAOindex_5yrsummer');
%SOUTHEAST
southeastRwinter = corrcoef(nanmedian(southeast_normVcomplete)',NAOindex_5yrwinter');
southeastRsummer = corrcoef(nanmedian(southeast_normVcomplete)',NAOindex_5yrsummer');
%CENTRAL EAST
centraleastRwinter = corrcoef(nanmedian(centraleast_normVcomplete)',NAOindex_5yrwinter');
centraleastRsummer = corrcoef(nanmedian(centraleast_normVcomplete)',NAOindex_5yrsummer');
%NORTHEAST
northeastRwinter = corrcoef(nanmedian(northeast_normVcomplete)',NAOindex_5yrwinter');
northeastRsummer = corrcoef(nanmedian(northeast_normVcomplete)',NAOindex_5yrsummer');
% %NORTH
% northRwinter = corrcoef(nanmedian(north_normVcomplete)',NAOindex_5yrwinter');
% northRsummer = corrcoef(nanmedian(north_normVcomplete)',NAOindex_5yrsummer');

%create normalized speeds vs NAO indices & list correlation coefficients in legend
clear pi pv;
NAO_figure = figure; set(gcf,'position',[50 50 800 400]);
set(NAO_figure, 'DefaultTextFontSize', 20);
% subtimeseries = subplot(2,1,1); set(gca,'fontsize',20);
yyaxis right; ax1 = gca; ax1.YColor = 'k'; ax1.YDir = 'reverse'; ax1.YLabel.String = 'Winter NAO index'; hold on;
pi(1) = plot(D_yrs,NAOindex_10yrwinter,'color','none'); hold on;
pi(2) = plot(D_yrs,NAOindex_10yrwinter,'color','none'); hold on;
pi(3) = plot(D_yrs,NAOindex_winter(find(NAOyrs == min(D_yrs)):find(NAOyrs == max(D_yrs))),'-k','linewidth',2); hold on;
pi(4) = plot(D_yrs,NAOindex_5yrwinter,'--k','linewidth',1); hold on;
pi(5) = plot(D_yrs,NAOindex_10yrwinter,'color','none'); hold on;
pi(6) = plot(D_yrs,NAOindex_10yrwinter,'color','none'); hold on;
ax1.XLim = [min(D_yrs),max(D_yrs)];
legi = legend(pi,'annual NAO','5-year average NAO','10-year average NAO','fontsize',16,'location','northeast');
yyaxis left; ax2 = gca; ax2.YColor = 'k'; ax2.YLabel.String = 'Median normalized speed';
pv(1) = plot(D_yrs,nanmedian(west_normVcomplete),'-','color',region_cmap(1,:),'linewidth',2); hold on;
pv(2) = plot(D_yrs,nanmedian(southeast_normVcomplete),'-','color',region_cmap(2,:),'linewidth',2); hold on;
plot(D_yrs,nanmedian(centraleast_normVcomplete),'-','color','k','linewidth',2.5); hold on;
pv(3) = plot(D_yrs,nanmedian(centraleast_normVcomplete),'-','color',region_cmap(3,:),'linewidth',2); hold on;
pv(4) = plot(D_yrs,nanmedian(northeast_normVcomplete),'-','color',region_cmap(4,:),'linewidth',2); hold on;
ax2.XLim = [min(D_yrs),max(D_yrs)]; grid on; set(gca,'box','on','fontsize',20);
xlabel('Year ','fontsize',20); %text(min(D_yrs)+0.5,0.075,'a) ','fontsize',20);
figpos = get(gca,'position');
legv = legend([pi pv],' ',' ','annual','5-year',' ',' ',['west (r=',num2str(round(westRwinter(2,1),2)),')'],['southeast (r=',num2str(round(southeastRwinter(2,1),2)),')'],...
    ['central east (r=',num2str(round(centraleastRwinter(2,1),2)),')'],['northeast (r=',num2str(round(northeastRwinter(2,1),2)),')'],...
    'fontsize',16,'location','northoutside','NumColumns',5,'box','on');
set(gca,'position',[figpos(1) figpos(2) figpos(3) 0.9*figpos(4)]);
legpos = get(legv,'position'); set(legv,'position',[figpos(1) figpos(2)+0.925*figpos(4) figpos(3) legpos(4)]);
annotation('textbox','position',[figpos(1) figpos(2)+0.975*figpos(4) 0.15 0.05],...
    'string','NAO_{winter}','fontsize',16,'fontweight','bold','linestyle','none');
annotation('textbox','position',[figpos(1)+0.285*figpos(3) figpos(2)+0.975*figpos(4) 0.05 0.05],...
    'string','Speed','fontsize',16,'fontweight','bold','linestyle','none');
annotation('line',[figpos(1)+0.28*figpos(3) figpos(1)+0.28*figpos(3)],[figpos(2)+0.925*figpos(4) figpos(2)+0.925*figpos(4)+legpos(4)],'linewidth',1.25);
%save the figure
saveas(gcf,[root_dir,'figures/','JOG-21-0146.Figure7.png'],'png'); 
saveas(gcf,[root_dir,'figures/','JOG-21-0146.Figure7.eps'],'epsc');



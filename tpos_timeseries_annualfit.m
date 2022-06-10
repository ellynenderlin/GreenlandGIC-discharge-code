%% Code to fit a Fourier model to terminus position timeseries
% Jukes Liu
% Last modified 11-14-2020
% Dependencies: getTestTrain.m, testmodelerror.m, Financial Toolbox
clear all; close all;

%load the average ITS_LIVE annual velocity map as background
cd /Users/ellynenderlin/Research/miscellaneous/Greenland-ITSLIVE_1985-2018/
map_fig = figure; set(map_fig,'position',[50 50 600 1200]); 
[A,R] = readgeoraster('Greenland_ITSLIVEavg_1985-2018.tif');
vx = R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX;
vy = R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY;
v = A; clear A R;
sub_map = subplot(7,1,[2:7]);
imagesc(vx,vy,v); axis xy equal; hold on;
vel_cmap = colormap(parula(10001)); vel_cmap(1,:) = [1 1 1]; colormap(gca,vel_cmap);
cbar = colorbar; set(gca,'clim',[0 500]); cbar.Label.String = 'speed (m/yr)';
set(gca,'ylim',[-3.35e6 -0.65e6],'xlim',[-6.5e5 8.5e5],...
    'ytick',[-3.25e6:0.5e6:-0.75e6],'yticklabel',[-3250:500:-750],...
    'xtick',[-6e5:2e5:8e5],'xticklabel',[-600:200:800],...
    'fontsize',18);
xlabel('Easting (km) ','fontsize',20); ylabel('Northing (km) ','fontsize',20); 

%plot the glacier locations
cd /Users/ellynenderlin/Research/NASA_GreenlandPeriph-Mapping/
load Greenland_GIC_centerlines.mat
for i = 1:length(term)
    X(i) = nanmean(term(i).X); Y(i) = nanmean(term(i).Y); Box(i) = term(i).BoxID;
    plot(X(i),Y(i),'ok','markerfacecolor','w','markersize',10); hold on;
%     text(X(i),Y(i),num2str(term(i).BoxID));
end
Wrefs = [1:1:63]; %find(X<-1e5 & Y<-1.2e6);
Nrefs = [524:1:641]; %find(Y>=-1.2e6);
NErefs = [501:1:523];
CErefs = [311:1:500];
SErefs = [64:1:310];

% set path to folder with filtered terminus position timeseries
tspath = '/Users/ellynenderlin/Research/NASA_GreenlandPeriph-Mapping/terminus_timeseries/';
cd_to_dir = ['cd ',tspath]; eval(cd_to_dir);
flowline = '50'; ts_files = dir(['*flowline',flowline,'.csv']);

%loop through the data and plot according to region
Wlims = []; Nlims = []; NElims = []; CElims = []; SElims = []; 
Wslope = NaN(length(ts_files),7); Nslope = NaN(length(ts_files),7); NEslope = NaN(length(ts_files),7); CEslope = NaN(length(ts_files),7); SEslope = NaN(length(ts_files),7); 
Wintercept = NaN(length(ts_files),7); Nintercept = NaN(length(ts_files),7); NEintercept = NaN(length(ts_files),7); CEintercept = NaN(length(ts_files),7); SEintercept = NaN(length(ts_files),7); 
Waccel = NaN(length(ts_files),7); Naccel = NaN(length(ts_files),7); NEaccel = NaN(length(ts_files),7); CEaccel = NaN(length(ts_files),7); SEaccel = NaN(length(ts_files),7); 

%% create a video that loops through the time series plots

%uncomment to make a video instead of a fig
% term_timeseries = VideoWriter('GreenlandGIC-auto-termchange-video.avi'); term_timeseries.FrameRate=5;
% open(term_timeseries); 

%loop through terminus timeseries and plot
for i = 1:length(ts_files)
    
    % identify the timeseries to plot
%     BoxID = '040'; disp(['Box', BoxID]); % the BoxID of the glacier
%     flowline = '50'; disp(['flowline', flowline]); % flowline (25, 50, or 75)
    BoxID = ts_files(i).name(end-17:end-15);
    [~,boxref] = min(abs(Box-str2num(BoxID))); 
    figure(map_fig); subplot(sub_map); plot(X(boxref),Y(boxref),'ok','markerfacecolor','m','markersize',10);

    % 1) LOAD TIMESERIES DATA
    tposfile = ['Tpos_filt_timeseries_Box',BoxID,'_flowline',flowline,'.csv']; % find the terminus position timeseries file
    tpos_table = readtable(strcat(tspath, tposfile)); % read tpos in as a table
    tpos_dates = table2array(tpos_table(:,1)); % grab the tpos timepoints
    tpos = table2array(tpos_table(:,2)); % grab the terminus positions
    
    % 2) CONVERT DATES TO DAYS AFTER 2013-01-01 (for model fitting)
    day1 = datetime('2013-01-01','InputFormat','yyyy-MM-dd');  % grab first day of 2013
    t_tpos = zeros(1, numel(tpos_dates)); % create vector to store the runoff timepoints (days after 2013-01-01)
    for b=1:length(tpos_dates)
        t_tpos(b) = datenum(tpos_dates(b)) - datenum(day1); % grab number of days after day 1
    end
    
    % 3) FOURIER MODEL-FITTING
    % OPTION A: Fit 3 Fourier curves to 90% of the data, auto choose the model with the lowest error
    pTrain = 0.9; % set percent training data
    [tpos_train, tpos_test] = getTrainTest([t_tpos' tpos], pTrain); % split into test and training datasets (use getTestTrain.m)
    n_tpos = round((max(datenum(tpos_dates))-min(datenum(tpos_dates)))/365); % grab number of years of tpos data
    
    % try fourier fits with 3 different number of terms to tpos data
    fouriern=strcat('fourier',string(n_tpos)); % # terms = # years of data
    fouriern_1=strcat('fourier',string(n_tpos-1)); % # terms = # years of data - 1
    fouriernplus1=strcat('fourier',string(n_tpos+1)); % # terms = # years of data + 1
    
    % fit Fourier curves to the tpos training data:
    [fn, g] = fit(tpos_train(:,1),tpos_train(:,2),fouriern);
    [fn_1, g_1] = fit(tpos_train(:,1),tpos_train(:,2),fouriern_1);
    [fnplus1, gplus1] = fit(tpos_train(:,1),tpos_train(:,2), fouriernplus1);
    
    % calculate error using testing data (use testmodelerror.m)
    error_fn = testmodelerror(fn, tpos_test);
    error_fn_1 = testmodelerror(fn_1, tpos_test);
    error_fnplus1 = testmodelerror(fnplus1, tpos_test);
    
    % Find best Fourier model
    [M,I] = min([error_fn_1 error_fn error_fnplus1]); % find the Fourier model with the lowest error
    if I == 1; fourierfit=fn_1; end % # terms = # years of data - 1
    if I == 2; fourierfit=fn; end % # terms = # years of data
    if I == 3; fourierfit=fnplus1; end  % # terms = # years of data + 1
    
    % OPTION B: Fit Fourier curve with number of terms manually using 100% of data (use if Option A looks bad)
    %fourierfit=fit(t_tpos', tpos, 'fourier8'); % number indicates number of terms
    
    % 4) LINEAR MODEL FITTING: fit annual linear trendlines
%     ts_fig = figure; set(ts_fig,'position',[850 50 800 400]); figure(ts_fig);
    figure(map_fig); subplot(7,1,1);
    for k = 2014:2020
        early = datetime([num2str(k),'-01-01'], 'InputFormat', 'yyyy-MM-dd'); % grab day representing May 2016
        later = datetime([num2str(k+1),'-01-01'], 'InputFormat', 'yyyy-MM-dd'); % grab day representing start of 2017
        if length(t_tpos(tpos_dates > early & tpos_dates <= later)) > 3
        flin = fit(t_tpos(tpos_dates > early & tpos_dates <= later)', tpos(tpos_dates > early & tpos_dates <= later), 'poly1'); % fit linear trend to pre2016 data
        yhat = feval(flin, t_tpos(tpos_dates > early & tpos_dates <= later)'); % evaluate only at the proper timepoints
        plot(tpos_dates(tpos_dates > early & tpos_dates <= later), yhat-min(tpos), 'k-', 'Linewidth', 1,'color',[0.5 0.5 0.5]); hold on % plot pre2016 linear trend
        
        % 5) PLOT THE MODELS AND DATA
        if str2num(BoxID) >= min(Wrefs) && str2num(BoxID) <= max(Wrefs)
%             subplot(subW); 
            Wlims = [Wlims; max(yhat)-min(yhat)];
            Wslope(i,k-2013) = flin.p1; Wintercept(i,k-2013) = flin.p2;
            if (29+(k-2014+1)) <= length(term(boxref).fluxVavg)
                Waccel(i,k-2013) = (term(boxref).fluxVavg(29+(k-2014+1))-term(boxref).fluxVavg(29+(k-2014)));
            end
        elseif str2num(BoxID) >= min(Nrefs) && str2num(BoxID) <= max(Nrefs)
%             subplot(subN); 
            Nlims = [Nlims; max(yhat)-min(yhat)];
            Nslope(i,k-2013) = flin.p1; Nintercept(i,k-2013) = flin.p2;
            if (29+(k-2014+1)) <= length(term(boxref).fluxVavg)
                Naccel(i,k-2013) = (term(boxref).fluxVavg(29+(k-2014+1))-term(boxref).fluxVavg(29+(k-2014)));
            end
        elseif str2num(BoxID) >= min(NErefs) && str2num(BoxID) <= max(NErefs)
%             subplot(subNE); 
            NElims = [NElims; max(yhat)-min(yhat)];
            NEslope(i,k-2013) = flin.p1; NEintercept(i,k-2013) = flin.p2;
            if (29+(k-2014+1)) <= length(term(boxref).fluxVavg)
                NEaccel(i,k-2013) = (term(boxref).fluxVavg(29+(k-2014+1))-term(boxref).fluxVavg(29+(k-2014)));
            end
        elseif str2num(BoxID) >= min(CErefs) && str2num(BoxID) <= max(CErefs)
%             subplot(subCE); 
            CElims = [CElims; max(yhat)-min(yhat)];
            CEslope(i,k-2013) = flin.p1; CEintercept(i,k-2013) = flin.p2;
            if (29+(k-2014+1)) <= length(term(boxref).fluxVavg)
                CEaccel(i,k-2013) = (term(boxref).fluxVavg(29+(k-2014+1))-term(boxref).fluxVavg(29+(k-2014)));
            end
        else
%             subplot(subSE); 
            SElims = [SElims; max(yhat)-min(yhat)];
            SEslope(i,k-2013) = flin.p1; SEintercept(i,k-2013) = flin.p2;
            if (29+(k-2014+1)) <= length(term(boxref).fluxVavg)
                SEaccel(i,k-2013) = (term(boxref).fluxVavg(29+(k-2014+1))-term(boxref).fluxVavg(29+(k-2014)));
            end
        end
        
        end
    end
    
%     plot(tpos_dates(tpos_dates < may2016), yhat_pre-min([yhat_pre; yhat_post]), 'k-', 'Linewidth', 2); hold on % plot pre2016 linear trend
%     plot(tpos_dates(tpos_dates > firstday2017), yhat_post-min([yhat_pre; yhat_post]), 'k-', 'Linewidth', 2); hold on % plot post2016 linear trend
    plot(tpos_dates, tpos-min(tpos), 'xm', 'MarkerSize',10); hold on % plot the tpos data
%     title(['Box ',BoxID]);  
    set(gca, 'FontSize', 16); % set title and fonts
    tdaily = 1:1:t_tpos(end); datesdaily = day1:caldays(1):tpos_dates(end)-1; % evaluate the Fourier fit at daily resolution
    plot(datesdaily, feval(fourierfit, tdaily)-min(tpos), 'k-','linewidth',2); hold on;% plot the Fourier fit
    ylabel('Terminus position (m)');% set y-axis label
    xlabel('Date'); % set xlabel
    figpos = get(gca,'position'); set(gca,'position',[0.13 figpos(2) figpos(3) 1.15*figpos(4)]);
%     set(gca,'ylim',[0 (max(tpos)-min(tpos))+0.2*(min(tpos))],'fontsize',18);
    if max(tpos)-min(tpos) < 500
        set(gca,'ylim',[0 500],'fontsize',18);
    elseif max(tpos)-min(tpos) >= 500 && max(tpos)-min(tpos) < 1000
        set(gca,'ylim',[0 1000],'fontsize',18);
    else
        set(gca,'ylim',[0 2000],'fontsize',18);
    end
    grid on; drawnow;
    clear *tpos* fourier* *fn* M I;
%     saveas(ts_fig,[BoxID,'_terminus-timeseries.png'],'png'); close(ts_fig);
    M(i) = getframe(map_fig);
    writeVideo(term_timeseries,M(i));

    clear flin* yhat*;
end
figure(map_fig); saveas(map_fig,'GreenlandGIC_terminus-timeseries-coverage-map.png','png');

%uncomment line below if making a video
% close(term_timeseries);


%% box plots of annual slopes (retreat rates)
boxplot_ylim = [-2 2]; boxplot_ytick = [-2:1:2];
figure; set(gcf,'position',[850 50 600 1200]);
subplot(5,1,1);
boxplot(-Wslope,'Labels',{'2014', '2015','2016','2017','2018','2019','2020'}); set(gca,'ylim',boxplot_ylim,'ytick',[-2:1:2],'fontsize',18);
subplot(5,1,2);
boxplot(-Nslope,'Labels',{'2014', '2015','2016','2017','2018','2019','2020'}); set(gca,'ylim',boxplot_ylim,'ytick',[-2:1:2],'fontsize',18);
subplot(5,1,3);
boxplot(-NEslope,'Labels',{'2014', '2015','2016','2017','2018','2019','2020'}); set(gca,'ylim',boxplot_ylim,'ytick',[-2:1:2],'fontsize',18);
subplot(5,1,4);
boxplot(-CEslope,'Labels',{'2014', '2015','2016','2017','2018','2019','2020'}); set(gca,'ylim',boxplot_ylim,'ytick',[-2:1:2],'fontsize',18);
subplot(5,1,5);
boxplot(-SEslope,'Labels',{'2014', '2015','2016','2017','2018','2019','2020'}); set(gca,'ylim',boxplot_ylim,'ytick',[-2:1:2],'fontsize',18);
xlabel('Year','fontsize',20); 
ylbl = ylabel('Terminus retreat rate (m/d)','fontsize',20); set(ylbl,'position',[0.1071 15 -1]);
saveas(gcf,'GreenlandGIC_terminus-retreatrate_boxplots.png','png');

%% box plots of inter-annual speed change
boxplot_ylim = [-13 13]; boxplot_ytick = [-12:4:12];
figure; set(gcf,'position',[1250 50 600 1200]);
subplot(5,1,1);
boxplot(Waccel,'Labels',{'2014', '2015','2016','2017','2018','2019','2020'}); set(gca,'ylim',boxplot_ylim,'ytick',boxplot_ytick,'fontsize',18);
subplot(5,1,2);
boxplot(Naccel,'Labels',{'2014', '2015','2016','2017','2018','2019','2020'}); set(gca,'ylim',boxplot_ylim,'ytick',boxplot_ytick,'fontsize',18);
subplot(5,1,3);
boxplot(NEaccel,'Labels',{'2014', '2015','2016','2017','2018','2019','2020'}); set(gca,'ylim',boxplot_ylim,'ytick',boxplot_ytick,'fontsize',18);
subplot(5,1,4);
boxplot(CEaccel,'Labels',{'2014', '2015','2016','2017','2018','2019','2020'}); set(gca,'ylim',boxplot_ylim,'ytick',boxplot_ytick,'fontsize',18);
subplot(5,1,5);
boxplot(SEaccel,'Labels',{'2014', '2015','2016','2017','2018','2019','2020'}); set(gca,'ylim',boxplot_ylim,'ytick',boxplot_ytick,'fontsize',18);
xlabel('Year','fontsize',20); 
ylbl = ylabel('Acceleration (m/yr^2)','fontsize',20); set(ylbl,'position',[0.1071 49 -1]);
saveas(gcf,'GreenlandGIC_terminus-acceleration_boxplots.png','png');

%% plot regional median & IQR (or MAD) terminus change rate time series on left axis and velocity time series on the right axis
figure; set(gcf,'position',[1050 50 600 1200]);
subplot(5,1,1);
yyaxis left; 
plot([2014:1:2020],-nanmedian(Wslope),'-k','linewidth',2); hold on;
fill([2014:1:2020 2020:-1:2014],[-nanmedian(Wslope)+mad(Wslope,1) fliplr(-nanmedian(Wslope)-mad(Wslope,1))],'k','edgecolor','none','facealpha',0.35);
set(gca,'ycolor','k','fontsize',18,'ylim',[-1 2]);
text(2014.1,1.65,'a) West','fontsize',18);
yyaxis right; 
plot([2014:1:2020],nanmedian(Waccel),'-r','linewidth',2); hold on;
velMAD = mad(Waccel,1); velYR = [2014:1:2020];
fill([velYR(~isnan(velMAD)) fliplr(velYR(~isnan(velMAD)))],[nanmedian(Waccel(:,(~isnan(velMAD))))+mad(Waccel(:,(~isnan(velMAD))),1) fliplr(nanmedian(Waccel(:,(~isnan(velMAD))))-mad(Waccel(:,(~isnan(velMAD))),1))],'r','edgecolor','none','facealpha',0.35);
set(gca,'ycolor','r','fontsize',18,'ylim',[-5 7]);

subplot(5,1,2);
yyaxis left; 
plot([2014:1:2020],-nanmedian(Nslope),'-k','linewidth',2); hold on;
fill([2014:1:2020 2020:-1:2014],[-nanmedian(Nslope)+mad(Nslope,1) fliplr(-nanmedian(Nslope)-mad(Nslope,1))],'k','edgecolor','none','facealpha',0.35);
set(gca,'ycolor','k','fontsize',18,'ylim',[-1 2]);
text(2014.1,1.65,'b) North','fontsize',18);
yyaxis right; 
plot([2014:1:2020],nanmedian(Naccel),'-r','linewidth',2); hold on;
velMAD = mad(Naccel,1); velYR = [2014:1:2020];
fill([velYR(~isnan(velMAD)) fliplr(velYR(~isnan(velMAD)))],[nanmedian(Naccel(:,(~isnan(velMAD))))+mad(Naccel(:,(~isnan(velMAD))),1) fliplr(nanmedian(Naccel(:,(~isnan(velMAD))))-mad(Naccel(:,(~isnan(velMAD))),1))],'r','edgecolor','none','facealpha',0.35);
set(gca,'ycolor','r','fontsize',18,'ylim',[-5 7]);

subplot(5,1,3);
yyaxis left; 
plot([2014:1:2020],-nanmedian(NEslope),'-k','linewidth',2); hold on;
fill([2014:1:2020 2020:-1:2014],[-nanmedian(NEslope)+mad(NEslope,1) fliplr(-nanmedian(NEslope)-mad(NEslope,1))],'k','edgecolor','none','facealpha',0.35);
set(gca,'ycolor','k','fontsize',18,'ylim',[-1 2]);
text(2014.1,1.65,'c) Northeast','fontsize',18);
yyaxis right; 
plot([2014:1:2020],nanmedian(NEaccel),'-r','linewidth',2); hold on;
velMAD = mad(NEaccel,1); velYR = [2014:1:2020];
fill([velYR(~isnan(velMAD)) fliplr(velYR(~isnan(velMAD)))],[nanmedian(NEaccel(:,(~isnan(velMAD))))+mad(NEaccel(:,(~isnan(velMAD))),1) fliplr(nanmedian(NEaccel(:,(~isnan(velMAD))))-mad(NEaccel(:,(~isnan(velMAD))),1))],'r','edgecolor','none','facealpha',0.35);
set(gca,'ycolor','r','fontsize',18,'ylim',[-5 7]);

subplot(5,1,4);
yyaxis left; 
plot([2014:1:2020],-nanmedian(CEslope),'-k','linewidth',2); hold on;
fill([2014:1:2020 2020:-1:2014],[-nanmedian(CEslope)+mad(CEslope,1) fliplr(-nanmedian(CEslope)-mad(CEslope,1))],'k','edgecolor','none','facealpha',0.35);
set(gca,'ycolor','k','fontsize',18,'ylim',[-1 2]);
text(2014.1,1.65,'d) East','fontsize',18);
yyaxis right; 
plot([2014:1:2020],nanmedian(CEaccel),'-r','linewidth',2); hold on;
velMAD = mad(CEaccel,1); velYR = [2014:1:2020];
fill([velYR(~isnan(velMAD)) fliplr(velYR(~isnan(velMAD)))],[nanmedian(CEaccel(:,(~isnan(velMAD))))+mad(CEaccel(:,(~isnan(velMAD))),1) fliplr(nanmedian(CEaccel(:,(~isnan(velMAD))))-mad(CEaccel(:,(~isnan(velMAD))),1))],'r','edgecolor','none','facealpha',0.35);
set(gca,'ycolor','r','fontsize',18,'ylim',[-5 7]);

subplot(5,1,5);
yyaxis left; 
plot([2014:1:2020],-nanmedian(SEslope),'-k','linewidth',2); hold on;
fill([2014:1:2020 2020:-1:2014],[-nanmedian(SEslope)+mad(SEslope,1) fliplr(-nanmedian(SEslope)-mad(SEslope,1))],'k','edgecolor','none','facealpha',0.35);
set(gca,'ycolor','k','fontsize',18,'ylim',[-1 2]);
ylbl1 = ylabel('Terminus retreat rate (m/d)','fontsize',20); set(ylbl1,'position',[2013.7 8.5 -1]);
text(2014.1,1.65,'e) Southeast','fontsize',18);
yyaxis right; 
plot([2014:1:2020],nanmedian(SEaccel),'-r','linewidth',2); hold on;
velMAD = mad(SEaccel,1); velYR = [2014:1:2020];
fill([velYR(~isnan(velMAD)) fliplr(velYR(~isnan(velMAD)))],[nanmedian(SEaccel(:,(~isnan(velMAD))))+mad(SEaccel(:,(~isnan(velMAD))),1) fliplr(nanmedian(SEaccel(:,(~isnan(velMAD))))-mad(SEaccel(:,(~isnan(velMAD))),1))],'r','edgecolor','none','facealpha',0.35);
set(gca,'ycolor','r','fontsize',18,'ylim',[-5 7]);
xlabel('Year ','fontsize',20); ylbl2 = ylabel('Acceleration (m/yr^2)','fontsize',20); set(ylbl2,'position',[2020.3 34 -1]);
saveas(gcf,'GreenlandGIC_terminuschange-acceleration_annual-timeseries.png','png');

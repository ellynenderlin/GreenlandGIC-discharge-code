%% Code to fit a Fourier model to terminus position timeseries
% Jukes Liu
% Last modified 11-14-2020
% Dependencies: getTestTrain.m, testmodelerror.m, Financial Toolbox
clear all;

%load GIMP masks as background
cd /Users/ellynenderlin/Research/miscellaneous/GIMP
load GIMPmasks.mat
map_fig = figure; set(map_fig,'position',[50 50 800 1600]); map_cmap = [1 1 1; 0.5 0.5 0.5; 0 1 1];
imagesc(GIMP.x,GIMP.y,2*GIMP.ice+GIMP.land); axis xy equal; colormap(map_cmap);

%plot the glacier locations
cd /Users/ellynenderlin/Research/NASA_Greenland-Periph-Mapping/centerlines/
load Greenland_GIC_centerlines.mat
for i = 1:length(term)
    X(i) = nanmean(term(i).X); Y(i) = nanmean(term(i).Y); Box(i) = term(i).BoxID;
    text(X(i),Y(i),num2str(term(i).BoxID));
end
Wrefs = [1:1:63]; %find(X<-1e5 & Y<-1.2e6);
Nrefs = [524:1:641]; %find(Y>=-1.2e6);
NErefs = [501:1:523];
CErefs = [311:1:500];
SErefs = [64:1:310];

% set path to folder with filtered terminus position timeseries
tspath = '/Users/ellynenderlin/Research/NASA_Greenland-Periph-Mapping/PeripheralGICs_timeseries/';
cd_to_dir = ['cd ',tspath]; eval(cd_to_dir);
flowline = '50'; ts_files = dir(['*flowline',flowline,'.csv']);

%loop through the data and plot according to region
ts_fig = figure; set(ts_fig,'position',[850 50 800 1200]);
subW = subplot(5,1,1); subN = subplot(5,1,2); subNE = subplot(5,1,3); subCE = subplot(5,1,4); subSE = subplot(5,1,5); 
Wlims = []; Nlims = []; NElims = []; CElims = []; SElims = []; 
Wprefit = []; Nprefit = []; NEprefit = []; CEprefit = []; SEprefit = []; 
Wpostfit = []; Npostfit = []; NEpostfit = []; CEpostfit = []; SEpostfit = []; 
for i = 1:length(ts_files)
    
    % identify the timeseries to plot
%     BoxID = '040'; disp(['Box', BoxID]); % the BoxID of the glacier
%     flowline = '50'; disp(['flowline', flowline]); % flowline (25, 50, or 75)
    BoxID = ts_files(i).name(end-17:end-15);

    % 1) LOAD TIMESERIES DATA
    tposfile = ['Tpos_filt_timeseries_Box',BoxID,'_flowline',flowline,'.csv']; % find the terminus position timeseries file
    tpos_table = readtable(strcat(tspath, tposfile)); % read tpos in as a table
    tpos_dates = table2array(tpos_table(:,1)); % grab the tpos timepoints
    tpos = table2array(tpos_table(:,2)); % grab the terminus positions
    
    % 2) CONVERT DATES TO DAYS AFTER 2013-01-01 (for model fitting)
    day1 = datetime('2013-01-01','InputFormat','yyyy-MM-dd');  % grab first day of 2013
    t_tpos = zeros(1, numel(tpos_dates)); % create vector to store the runoff timepoints (days after 2013-01-01)
    for b=1:length(tpos_dates)
        t_tpos(b) = daysact(day1, tpos_dates(b)); % grab number of days after day 1
    end
    
    % 3) FOURIER MODEL-FITTING
    % OPTION A: Fit 3 Fourier curves to 90% of the data, auto choose the model with the lowest error
    pTrain = 0.9; % set percent training data
    [tpos_train, tpos_test] = getTrainTest([t_tpos' tpos], pTrain); % split into test and training datasets (use getTestTrain.m)
    n_tpos = round(daysact(min(tpos_dates), max(tpos_dates))/365); % grab number of years of tpos data
    
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
    
    % 4) LINEAR MODEL FITTING: fit linear trendline to pre2016 data and post2016 data
    may2016 = datetime('2016-05-01', 'InputFormat', 'yyyy-MM-dd'); % grab day representing May 2016
    firstday2017 = datetime('2017-01-01', 'InputFormat', 'yyyy-MM-dd'); % grab day representing start of 2017
    flin_pre = fit(t_tpos(tpos_dates < may2016)', tpos(tpos_dates < may2016), 'poly1'); % fit linear trend to pre2016 data
    flin_post = fit(t_tpos(tpos_dates > firstday2017)', tpos(tpos_dates > firstday2017), 'poly1'); % fit linear trend to post2016 data
    yhat_pre = feval(flin_pre, t_tpos(tpos_dates < may2016)'); % evaluate only at the proper timepoints
    yhat_post = feval(flin_post, t_tpos(tpos_dates > firstday2017)'); % evaluate only at the proper timepoints
    
    % 5) PLOT THE MODELS AND DATA
    if str2num(BoxID) >= min(Wrefs) && str2num(BoxID) <= max(Wrefs)
        subplot(subW); Wlims = [Wlims; max([yhat_pre; yhat_post])-min([yhat_pre; yhat_post])];
        Wprefit = [Wprefit; flin_pre.p1 flin_pre.p2]; Wpostfit = [Wpostfit; flin_post.p1 flin_post.p2];
    elseif str2num(BoxID) >= min(Nrefs) && str2num(BoxID) <= max(Nrefs)
        subplot(subN); Nlims = [Nlims; max([yhat_pre; yhat_post])-min([yhat_pre; yhat_post])];
        Nprefit = [Nprefit; flin_pre.p1 flin_pre.p2]; Npostfit = [Npostfit; flin_post.p1 flin_post.p2];
    elseif str2num(BoxID) >= min(NErefs) && str2num(BoxID) <= max(NErefs)
        subplot(subNE); NElims = [NElims; max([yhat_pre; yhat_post])-min([yhat_pre; yhat_post])];
        NEprefit = [NEprefit; flin_pre.p1 flin_pre.p2]; NEpostfit = [NEpostfit; flin_post.p1 flin_post.p2];
    elseif str2num(BoxID) >= min(CErefs) && str2num(BoxID) <= max(CErefs)
        subplot(subCE); CElims = [CElims; max([yhat_pre; yhat_post])-min([yhat_pre; yhat_post])];
        CEprefit = [CEprefit; flin_pre.p1 flin_pre.p2]; CEpostfit = [CEpostfit; flin_post.p1 flin_post.p2];
    else
        subplot(subSE); SElims = [SElims; max([yhat_pre; yhat_post])-min([yhat_pre; yhat_post])];
        SEprefit = [SEprefit; flin_pre.p1 flin_pre.p2]; SEpostfit = [SEpostfit; flin_post.p1 flin_post.p2];
    end
    plot(tpos_dates(tpos_dates < may2016), yhat_pre-min([yhat_pre; yhat_post]), 'k-', 'Linewidth', 2); hold on % plot pre2016 linear trend
    plot(tpos_dates(tpos_dates > firstday2017), yhat_post-min([yhat_pre; yhat_post]), 'k-', 'Linewidth', 2); hold on % plot post2016 linear trend
%     plot(tpos_dates, tpos, 'b.', 'MarkerSize',15); hold on % plot the tpos data
%     title(['Box',BoxID,' flowline',flowline]);  set(gca, 'FontSize', 16); % set title and fonts
%     tdaily = 1:1:t_tpos(end); datesdaily = day1:caldays(1):tpos_dates(end)-1; % evaluate the Fourier fit at daily resolution
%     plot(datesdaily, feval(fourierfit, tdaily)-min(tpos), 'b-'); hold on;% plot the Fourier fit
%     set(gca,'xlim',[2013.5 2020]);
    ylabel('Terminus position (m)');% set y-axis label 
    xlabel('Date'); % set xlabel
    grid on; drawnow;
    legend off
    % saveas(gcf, ['/Users/jukesliu/Documents/AUTO-TERMINUS/Figures/Tpos_timeseries_Box',BoxID,'.jpg'], 'jpg'); % save
    
    clear *tpos* fourier* *fn* M I flin* yhat*;
end
set(subW,'ylim',[0 max(max(Wlims))]);
set(subN,'ylim',[0 max(max(Nlims))]);
set(subNE,'ylim',[0 max(max(NElims))]);
set(subCE,'ylim',[0 max(max(CElims))]);
set(subSE,'ylim',[0 max(max(SElims))]);


%box plots of pre-2016 & post-2016 slopes


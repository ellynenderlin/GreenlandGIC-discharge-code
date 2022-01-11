% H vs U/W
clearvars; close all;


%column key:
%3=BoxID 4=thickness 6=velocity 8=width 9=u/w 10=uw 11=huw
% velocities are m/yr

% cd /Users/katebollen/Documents/MS/data/HxUxW
cd /users/ellynenderlin/Research/NASA_Greenland-Periph-Mapping/
avg_thick = csvread('avg_thick_vel_2.csv');
thick = [avg_thick(:,3) avg_thick(:, 11)]; %D (yaxis)
uxw = [avg_thick(:,3) avg_thick(:,10)]; %UW (xaxis) this is u/w, i just left it as uxw for ease of use from a previous section


%Find all indices for the 4 boxes
n = [find(avg_thick(:,3)==003); find(avg_thick(:,3)==093); find(avg_thick(:,3)==410); find(avg_thick(:,3)==564)];
uxw_lim = uxw(n,1:2); thick_lim = thick(n,1:2);
 
%Fit linear and quadratic polynomials to the data
% f = fitlm(uxw_lim(:,2),thick_lim(:,2),'poly1', 'Intercept', false);
options1 = fitoptions('poly1','Lower',[-Inf 0],'Upper',[Inf Inf]);
options2 = fitoptions('poly2','Lower',[-Inf -Inf 0],'Upper',[Inf Inf Inf]);
[f1,gof1] = fit(uxw_lim(:,2),thick_lim(:,2),'poly1',options1);
disp(['Linear polynomial R^2=',num2str(gof1.rsquare),' and RMSE=',num2str(gof1.rmse)]);
[f2,gof2] = fit(uxw_lim(:,2),thick_lim(:,2),'poly2',options2);
disp(['Quadratic polynomial R^2=',num2str(gof2.rsquare),' and RMSE=',num2str(gof2.rmse)]);
%create a sigmoidal function so that D doesn't jump much after the data extent
ft3 = fittype('a+a*erf((x-b)/c)',...
    'dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c'});
[f3,gof3] = fit(uxw_lim(:,2),thick_lim(:,2),ft3,'start',[2.5e7,2e5,1e5]);
disp(['Sigmoidal polynomial R^2=',num2str(gof3.rsquare),' and RMSE=',num2str(gof3.rmse)]);
cif1 = confint(f1,0.68);cif2 = confint(f2,0.68); cif3 = confint(f3,0.68);
cif1(isnan(cif1))=0; cif2(isnan(cif2))=0;
x = sort(uxw_lim(:, 2));
% int = min(uxw_lim(:,2)):max(uxw_lim(:,2));
% ci_upper = (uxw_lim(:,2)*cif2(1,1)+cif2(1,2))*int;
% ci_lower = (uxw_lim(:,2)*cif2(2,1)+cif2(2,2))*int;
poly1 = feval(f1,x);
poly2 = feval(f2,x);
poly3 = feval(f3,x);
ci1_lower = cif1(1,1).*x + cif1(1,2); ci1_upper = cif1(2,1).*x + cif1(2,2);
ci2_lower = cif2(1,1).*x.^2 + cif2(1,2).*x + cif2(1,3); ci2_upper = cif2(2,1).*x.^2 + cif2(2,2).*x + cif2(2,3);
%set the function outside the data domain as a linear polynomial to avoid
%erroneous negative relationships between H and U/W at very large values
f2ext.p1 = ((f2.p1.*max(uxw_lim(:,2)).^2 + f2.p2.*max(uxw_lim(:,2)) + f2.p3)-(f2.p1.*max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2)).^2 + f2.p2.*max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2)) + f2.p3))./(max(uxw_lim(:,2))-max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2)));
f2ext.p2 = (f2.p1.*max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2)).^2 + f2.p2.*max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2)) + f2.p3)-f2ext.p1*max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2));
%use quadratic polynomial confidence interval to define a linear polynomial
%uncertainty envelope for the extrapolated region of data
f2ext.cil_p1 = ((cif2(1,1).*max(uxw_lim(:,2)).^2 + cif2(1,2).*max(uxw_lim(:,2)) + cif2(1,3))-(cif2(1,1).*max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2)).^2 + cif2(1,2).*max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2)) + cif2(1,3)))./(max(uxw_lim(:,2))-max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2)));
f2ext.cil_p2 = (cif2(1,1).*max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2)).^2 + cif2(1,2).*max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2)) + cif2(1,3))-f2ext.cil_p1*max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2));
f2ext.ciu_p1 = ((cif2(2,1).*max(uxw_lim(:,2)).^2 + cif2(2,2).*max(uxw_lim(:,2)) + cif2(2,3))-(cif2(2,1).*max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2)).^2 + cif2(2,2).*max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2)) + cif2(2,3)))./(max(uxw_lim(:,2))-max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2)));
f2ext.ciu_p2 = (cif2(2,1).*max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2)).^2 + cif2(2,2).*max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2)) + cif2(2,3))-f2ext.ciu_p1*max(uxw_lim(uxw_lim(:,2)<max(uxw_lim(:,2)),2));
%define uncertainty envelop of sigmoidal function by fitting curves to the
%two largest datasets
D=thick_lim(:,2); UW=uxw_lim(:,2); refs = uxw_lim(:,1);
% [cifu3,~] = fit(UW(refs==564),D(refs==564),ft3,'start',[nanmean(D(refs==564)),f3.b,f3.c]); ci3_upper = feval(cifu3,x);
% [cifl3,~] = fit(UW(refs==3),D(refs==3),ft3,'start',[nanmean(D(refs==3)),f3.b,f3.c]); ci3_lower = feval(cifl3,x);
% ci3_upper = poly3+(poly3-ci3_lower);
lower_err = feval(f3,UW(refs==3))-D(refs==3); upper_err = D(refs==564)-feval(f3,UW(refs==564));
overlap_refs = find(UW(refs==564)>= min(UW(refs==3)) & UW(refs==564)<= max(UW(refs==3)));
UWsub = UW(refs==564);
uncert = fit([UW(refs==3); UWsub(overlap_refs)],[lower_err; upper_err(overlap_refs)],'poly1','Lower',[-Inf 0]);


%Plot a figure
figure; hold on; grid on;
xlabel('UW (m^2/yr)'), ylabel('D (m^3/yr)'), title('D vs UW');
%     ci_fnc = ([int fliplr(int); ci_upper fliplr(ci_lower)]);
%     fill_shade = fill(ci_fnc(1,:),ci_fnc(2,:),[0.3010, 0.7450, 0.9330], 'FaceAlpha', 0.4);
pl(1) = plot(x,poly1,'-k','linewidth',2); hold on; %fill([x; flipud(x)],[ci1_upper; flipud(ci1_lower)],'k', 'FaceAlpha', 0.4); hold on;
pl(2) = plot(x,poly2,'--b','linewidth',2); hold on; %fill([x; flipud(x)],[ci2_upper; flipud(ci2_lower)],'b', 'FaceAlpha', 0.4); hold on;
pl(3) = plot(x,poly3,'--m','linewidth',2); hold on; %fill([x; flipud(x)],[ci3_upper; flipud(ci3_lower)],'m', 'FaceAlpha', 0.4); hold on;
% %plot the extrapolated equations to check they make sense
% plot([max(uxw_lim(:,2)) 4e5],f2ext.p1.*[max(uxw_lim(:,2)) 4e5]+f2ext.p2,'-b','linewidth',2); hold on;
% plot([max(uxw_lim(:,2)) 4e5],f2ext.cil_p1.*[max(uxw_lim(:,2)) 4e5]+f2ext.cil_p2,'-b','linewidth',1); hold on;
% plot([max(uxw_lim(:,2)) 4e5],f2ext.ciu_p1.*[max(uxw_lim(:,2)) 4e5]+f2ext.ciu_p2,'-b','linewidth',1); hold on;
col = parula(20);
    

for i=1:length(uxw_lim(:,1))
    if uxw_lim(i,1)==003
        p2 = plot(uxw_lim(i,2), thick_lim(i,2), 's','color',col(4,:),'markersize',12, 'LineWidth', 3);
        %set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    elseif uxw_lim(i,1)==093
        p3 = plot(uxw_lim(i,2), thick_lim(i,2), 's','color',col(8,:),'markersize',12, 'LineWidth', 3);
        %set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    elseif uxw_lim(i,1)==410
        p4 = plot(uxw_lim(i,2), thick_lim(i,2), 's','color',col(13,:),'markersize',12, 'LineWidth', 3);
        %set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    elseif uxw_lim(i,1)==564
        p5 = plot(uxw_lim(i,2), thick_lim(i,2), 's','color',col(18,:),'markersize',12, 'LineWidth', 3);
        %set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end 
end 
 
set(gca,'FontSize',12,'FontName','Calibri');
%      p1 = plot(x,poly,'-r','linewidth',2);
legend([pl(1) pl(2) pl(3) p2 p3 p4 p5],['y=',num2str(round(f1.p1)),'x (R^2=',num2str(round(gof1.rsquare,2)),')'],...
    ['y=',num2str(round(f2.p1,3)),'x^2+',num2str(round(f2.p2,2)),'x (R^2=',num2str(round(gof2.rsquare,2)),')'],...
    ['y=',num2str(round(f3.a),2),'+',num2str(round(f3.a),2),'erf(x-',num2str(round(f3.b),2),')/',num2str(round(f3.c),2),' (R^2=',num2str(round(gof2.rsquare,2)),')'],...
    '003','093','410','564');
legend('Location','northwest');
set(gcf,'Units','centimeters','Position',[10 10 40 30]);
 
%Save the function
clear f;
f.xlims = [min(uxw_lim(:,2)) max(uxw_lim(:,2))];
f.fit = f2; f.ext = f2ext;
% f.fit = f3; f.uncert = uncert;
save('D_vs_UW_parab-functions.mat','f','-v7.3');


saveas(gcf,'HvsUW_agu.png')
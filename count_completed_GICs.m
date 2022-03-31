%calculate the number of glaciers for which flux gates have been drawn
clearvars; close all;
root_dir = '/Users/ellynenderlin/Research/NASA-Greenland-Periph-Mapping/centerlines/';
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
load Greenland_GIC_centerlines.mat;

%loop through the term structure & sum up the number of glaciers with gate data
counter = 0;
for i = 1:length(term)
    if ~isempty(term(i).gateX)
    counter = counter+1;
    end
end
disp([num2str(counter),' of ',num2str(length(term)),' glaciers completed']);
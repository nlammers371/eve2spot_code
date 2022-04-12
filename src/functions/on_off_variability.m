% Look at intrinsic-extrinsic variability in fraction on, time on, and time
% off
clear
close all

% make figure directory
figPath = '../fig/on_off_variability/';
mkdir(figPath)

% load master data set
load('../data/masterSet.mat')

%% subset for nuclei that are present for minimum window of time
time_bounds = [3 40]*60;
minDP = 10;
ap_time = 25;
nBoots = 100;
max_on_time = 60*50;

long_flags = false(1,length(masterSet));
on_flags = false(1,length(masterSet));
off_flags = false(1,length(masterSet));

for i = 1:length(masterSet)
    time_vec = masterSet(i).time;
    fluo_vec = masterSet(i).fluo;
    long_flags(i) = time_vec(1)<=time_bounds(1) & time_vec(end)>=time_bounds(2);
    on_flags(i) = isnan(fluo_vec(1)) && sum(~isnan(fluo_vec))>=minDP && min(time_vec(~isnan(fluo_vec)))<=max_on_time;
    off_flags(i) = isnan(fluo_vec(end)) & sum(~isnan(fluo_vec))>=minDP;
end    

%%%%%%%%%%%%%%%%%%%%
%% Fraction ever on

% get unique list of nuclei
nucleus_id_vec = [masterSet.ncID];
nucleus_id_index = unique(nucleus_id_vec(long_flags));

% get number of active loci for each nucleus
ever_on_vec = zeros(size(nucleus_id_index));
ap_on_vec = NaN(size(nucleus_id_index));
for n = 1:length(nucleus_id_index)
    nc_ids = find(nucleus_id_vec==nucleus_id_index(n));
    na = 0;
    for i = 1:length(nc_ids)
        na = na + (1*(sum(~isnan(masterSet(nc_ids(i)).fluo))>minDP));
    end
    ever_on_vec(n) = na;
    % add ap info
    ap_vec = masterSet(nc_ids(1)).APPosNucleus;
    time_vec = masterSet(nc_ids(1)).time;
    [~,mi] = min(abs(ap_time*60-time_vec));
    ap_on_vec(n) = ap_vec(mi);
end

% get fraction on as a function of ap position
ap_bins = linspace(0,100,101);
ap_on_array = NaN(nBoots,length(ap_bins)-1);
ap_both_array = NaN(nBoots,length(ap_bins)-1);
ap_one_array = NaN(nBoots,length(ap_bins)-1);
ap_none_array = NaN(nBoots,length(ap_bins)-1);

for n = 1:nBoots
    boot_indices = randsample(1:length(ever_on_vec),length(ever_on_vec),true);
    ever_on_boot = ever_on_vec(boot_indices);
    ap_on_boot = ap_on_vec(boot_indices);
    ap_ids = discretize(ap_on_boot,ap_bins);
    for a = 1:length(ap_bins)-1
        if sum(ap_ids==a)
            ap_on_array(n,a) = mean(ever_on_boot(ap_ids==a));
            ap_both_array(n,a) = mean(ever_on_boot(ap_ids==a)==2);
            ap_one_array(n,a) = mean(ever_on_boot(ap_ids==a)==1);
            ap_none_array(n,a) = mean(ever_on_boot(ap_ids==a)==0);
        end
    end
end    

% calculate predictions
ap_on_mean = mean(ap_on_array,1)/2;
ap_on_ste = std(ap_on_array,[],1)/2;
ap_on_ub = ap_on_mean + ap_on_ste;
ap_on_lb = ap_on_mean - ap_on_ste;

both_pd = ap_on_mean.^2;
one_pd = ap_on_mean.*(1-ap_on_mean)*2;
none_pd = (1-ap_on_mean).^2;

% calculate results
ap_both_mean = mean(ap_both_array,1);
ap_both_ste = std(ap_both_array,[],1);
ap_both_ub = ap_both_mean + 2*ap_both_ste;
ap_both_lb = ap_both_mean - 2*ap_both_ste;

ap_one_mean = mean(ap_one_array,1);
ap_one_ste = std(ap_one_array,[],1);
ap_one_ub = ap_one_mean + 2*ap_one_ste;
ap_one_lb = ap_one_mean - 2*ap_one_ste;

ap_none_mean = mean(ap_none_array,1);
ap_none_ste = std(ap_none_array,[],1);
ap_none_ub = ap_none_mean + 2*ap_none_ste;
ap_none_lb = ap_none_mean - 2*ap_none_ste;

%% Plot results
pt_filter = ~isnan(ap_none_mean)&~isnan(ap_on_mean)&~isnan(ap_one_mean)&~isnan(ap_both_mean);
close all
ap_axis = ap_bins(1:end-1) + diff(ap_bins);

ever_on_fig = figure;
hold on
cmap = brewermap([],'Set2');

% pon
% fill([ap_axis(pt_filter) fliplr(ap_axis(pt_filter))], [ap_on_lb(pt_filter) ...
%         fliplr(ap_on_ub(pt_filter))], cmap(8,:), 'FaceAlpha', 0.3, 'EdgeAlpha',0)
% plot(ap_axis(pt_filter),ap_on_mean(pt_filter),'Color',cmap(8,:),'LineWidth',1.5);

% both on
fill([ap_axis(pt_filter) fliplr(ap_axis(pt_filter))], [ap_both_lb(pt_filter) ...
        fliplr(ap_both_ub(pt_filter))], cmap(2,:), 'FaceAlpha', 0.3, 'EdgeAlpha',.3)
% plot(ap_axis(pt_filter),ap_both_mean(pt_filter),'Color',cmap(2,:),'LineWidth',1.5);
p1 = plot(ap_axis(pt_filter),both_pd(pt_filter),'-.','Color',cmap(2,:),'LineWidth',2);

% neither on
fill([ap_axis(pt_filter) fliplr(ap_axis(pt_filter))], [ap_none_lb(pt_filter) ...
        fliplr(ap_none_ub(pt_filter))], cmap(8,:), 'FaceAlpha', 0.3, 'EdgeAlpha',.3)
% plot(ap_axis(pt_filter),ap_none_mean(pt_filter),'Color',cmap(8,:),'LineWidth',1.5);
p2 = plot(ap_axis(pt_filter),none_pd(pt_filter),'-.','Color',cmap(8,:),'LineWidth',2);

% either on
fill([ap_axis(pt_filter) fliplr(ap_axis(pt_filter))], [ap_one_lb(pt_filter) ...
        fliplr(ap_one_ub(pt_filter))], cmap(3,:), 'FaceAlpha', 0.3, 'EdgeAlpha',.3)
% plot(ap_axis(pt_filter),ap_one_mean(pt_filter),'Color',cmap(3,:),'LineWidth',1.5);
p3 = plot(ap_axis(pt_filter),one_pd(pt_filter),'-.','Color',cmap(3,:),'LineWidth',2);

xlim([20 90])
ylim([0 1.3])
set(gca,'Fontsize',14)
xlabel('% embryo length')
ylabel('fraction')
grid on
legend([p1 p2 p3], 'both (p^2)','either (2(1-p))','neither ((1-p)^2))','Location','northwest');%'Orientation','horizontal')
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ever_on_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(ever_on_fig,[figPath 'fraction_ever_on_ap.png'])
saveas(ever_on_fig,[figPath 'fraction_ever_on_ap.pdf'])

%% On time
nucleus_id_vec_on = nucleus_id_vec(on_flags);
on_indices = find(on_flags);
nucleus_id_index_on = unique(nucleus_id_vec(on_flags));
on_time_vec1 = [];
on_time_vec2 = [];
ap_on_vec = [];
for i = 1:length(nucleus_id_index)
    if sum(nucleus_id_vec_on==nucleus_id_index(i)) == 2
        on_inds = randsample(on_indices(nucleus_id_vec_on==nucleus_id_index(i)),2,false);
        f_ind1 = find(~isnan(masterSet(on_inds(1)).fluo),1);
        on_time_vec1(end+1) = masterSet(on_inds(1)).time(f_ind1);
        f_ind2 = find(~isnan(masterSet(on_inds(2)).fluo),1);
        on_time_vec2(end+1) = masterSet(on_inds(2)).time(f_ind2);
        ap_on_vec(end+1) = [masterSet(on_inds(1)).APPosNucleus(f_ind1) + masterSet(on_inds(2)).APPosNucleus(f_ind2)]/2;
    end
end  

on_mdl = fitlm(on_time_vec1/60,on_time_vec2/60);
on1_axis = linspace(0,50);
on2_pd = predict(on_mdl,on1_axis');

on_time_fig = figure;
hold on
cmap = flipud(brewermap([],'set2'));
% colormap(cmap)
% both on
scatter(on_time_vec1/60,on_time_vec2/60,25,ap_on_vec,'s','MarkerFaceColor',cmap(4,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
plot(on1_axis, on2_pd, 'Color', 'k', 'LineWidth', 2)

% xlim([0 25])
% ylim([0 25])
set(gca,'Fontsize',14)
xlabel('spot 1')
ylabel('spot 2')
grid on
% legend([p1 p2 p3], 'both (p^2)','either (2(1-p))','neither ((1-p)^2))','Location','northwest');%'Orientation','horizontal')
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% h = colorbar;
% ylabel(h,'% embryo length')
on_time_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(on_time_fig,[figPath 'on_time_scatter.png'])
% saveas(on_time_fig,[figPath 'on_time_scatter.pdf'])

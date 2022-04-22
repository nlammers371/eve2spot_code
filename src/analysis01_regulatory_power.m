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
minOverlap = 5;
ap_time = 25;
nBoots = 100;
max_on_time = 60*50;
min_dur = 3;

long_flags = false(1,length(masterSet));
on_flags = false(1,length(masterSet));
off_flags = false(1,length(masterSet));

for i = 1:length(masterSet)
    time_vec = masterSet(i).time;
    fluo_vec = masterSet(i).fluo;
    qc_flag = masterSet(i).TraceQCFlag;
    long_flags(i) = time_vec(1)<=time_bounds(1) & time_vec(end)>=time_bounds(2);
    if length(fluo_vec) > min_dur && ~isnan(qc_flag)
        on_flags(i) = all(isnan(fluo_vec(1:min_dur))) && sum(~isnan(fluo_vec))>=minDP && min(time_vec(~isnan(fluo_vec)))<=max_on_time;% && qc_flag;
        off_flags(i) = all(isnan(fluo_vec(end-min_dur+1:end)));% && sum(~isnan(fluo_vec))>=minDP && qc_flag;
    end
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

both_pd_mean = mean((ap_on_array/2).^2);
both_pd_ste = std((ap_on_array/2).^2,[],1);
both_pd_ub = both_pd_mean + both_pd_ste;
both_pd_lb = both_pd_mean - both_pd_ste;

one_pd_mean = mean(ap_on_array/2.*(1-ap_on_array/2)*2,1);
one_pd_ste = std(ap_on_array/2.*(1-ap_on_array/2)*2,[],1);
one_pd_ub = one_pd_mean + one_pd_ste;
one_pd_lb = one_pd_mean - one_pd_ste;

none_pd_mean = mean((1-ap_on_array/2).^2,1);
none_pd_ste = std((1-ap_on_array/2).^2,[],1);
none_pd_ub = none_pd_mean + none_pd_ste;
none_pd_lb = none_pd_mean - none_pd_ste;

% calculate results
ap_both_mean = mean(ap_both_array,1);
ap_both_ste = std(ap_both_array,[],1);
ap_both_ub = ap_both_mean + ap_both_ste;
ap_both_lb = ap_both_mean - ap_both_ste;

ap_one_mean = mean(ap_one_array,1);
ap_one_ste = std(ap_one_array,[],1);
ap_one_ub = ap_one_mean + ap_one_ste;
ap_one_lb = ap_one_mean - ap_one_ste;

ap_none_mean = mean(ap_none_array,1);
ap_none_ste = std(ap_none_array,[],1);
ap_none_ub = ap_none_mean + ap_none_ste;
ap_none_lb = ap_none_mean - ap_none_ste;

% Plot results
pt_filter = ~isnan(ap_none_mean)&~isnan(ap_on_mean)&~isnan(ap_one_mean)&~isnan(ap_both_mean);
close all
ap_axis = ap_bins(1:end-1) + diff(ap_bins);
fillAlpha = 0.6;
ever_on_fig = figure;
hold on
cmap = brewermap([],'Set2');

% pon
% fill([ap_axis(pt_filter) fliplr(ap_axis(pt_filter))], [ap_on_lb(pt_filter) ...
%         fliplr(ap_on_ub(pt_filter))], cmap(8,:), 'FaceAlpha', 0.3, 'EdgeAlpha',0)


% both on
fill([ap_axis(pt_filter) fliplr(ap_axis(pt_filter))], [both_pd_lb(pt_filter) ...
        fliplr(both_pd_ub(pt_filter))], cmap(2,:), 'FaceAlpha', fillAlpha, 'EdgeAlpha',0)
% plot(ap_axis(pt_filter),ap_both_mean(pt_filter),'Color',cmap(2,:),'LineWidth',1.5);
errorbar(ap_axis(pt_filter),ap_both_mean(pt_filter),ap_both_ste(pt_filter),'.','Color',cmap(2,:),'CapSize',0);
p1 = scatter(ap_axis(pt_filter),ap_both_mean(pt_filter),'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k', 'MarkerEdgeAlpha',0.5);

% neither on
fill([ap_axis(pt_filter) fliplr(ap_axis(pt_filter))], [none_pd_lb(pt_filter) ...
        fliplr(none_pd_ub(pt_filter))], cmap(8,:), 'FaceAlpha', fillAlpha, 'EdgeAlpha',0)
% plot(ap_axis(pt_filter),ap_both_mean(pt_filter),'Color',cmap(2,:),'LineWidth',1.5);
errorbar(ap_axis(pt_filter),ap_none_mean(pt_filter),ap_none_ste(pt_filter),'.','Color',cmap(8,:),'CapSize',0);
p2 = scatter(ap_axis(pt_filter),ap_none_mean(pt_filter),'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k', 'MarkerEdgeAlpha',0.5);

% either on
fill([ap_axis(pt_filter) fliplr(ap_axis(pt_filter))], [one_pd_lb(pt_filter) ...
        fliplr(one_pd_ub(pt_filter))], cmap(3,:), 'FaceAlpha', fillAlpha, 'EdgeAlpha',0)
% plot(ap_axis(pt_filter),ap_both_mean(pt_filter),'Color',cmap(2,:),'LineWidth',1.5);
errorbar(ap_axis(pt_filter),ap_one_mean(pt_filter),ap_one_ste(pt_filter),'.','Color',cmap(3,:),'CapSize',0);
p3 = scatter(ap_axis(pt_filter),ap_one_mean(pt_filter),'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k', 'MarkerEdgeAlpha',0.5);

plot(ap_axis(pt_filter),ap_on_mean(pt_filter),'--k','LineWidth',2);

xlim([20 90])
ylim([0 1.1])
set(gca,'Fontsize',14)
xlabel('% embryo length')
ylabel('fraction')
% grid on
% legend([p1 p2 p3], 'both (p^2)','either (2(1-p))','neither ((1-p)^2))','Location','northwest');%'Orientation','horizontal')
set(gca,'Color',[228,221,209]/255) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ever_on_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(ever_on_fig,[figPath 'fraction_ever_on_ap.png'])
saveas(ever_on_fig,[figPath 'fraction_ever_on_ap.pdf'])

%% On time
markerSize = 30;
markerAlpha = 0.5;
edgeAlpha = 0.2;


close all
rng(325);

nucleus_id_vec_on = nucleus_id_vec(on_flags);%&long_flags);
on_indices = find(on_flags);%&long_flags);
nucleus_id_index_on = unique(nucleus_id_vec(on_flags));%&long_flags));
on_time_vec1 = [];
on_time_vec2 = [];
ap_on_vec = [];
for i = 1:length(nucleus_id_index)
    if sum(nucleus_id_vec_on==nucleus_id_index(i)) == 2
        % check for overlap
        on_inds = randsample(on_indices(nucleus_id_vec_on==nucleus_id_index(i)),2,false); 
        f_ind_vec1 = find(~isnan(masterSet(on_inds(1)).fluo));
        f_ind_vec2 = find(~isnan(masterSet(on_inds(2)).fluo));
        if sum(ismember(f_ind_vec1,f_ind_vec2))>=minOverlap                       
            on_time_vec1(end+1) = masterSet(on_inds(1)).time(f_ind_vec1(1))/60;           
            on_time_vec2(end+1) = masterSet(on_inds(2)).time(f_ind_vec2(1))/60;
%             ap_on_vec(end+1) = [masterSet(on_inds(1)).APPosNucleus(f_ind1) + masterSet(on_inds(2)).APPosNucleus(f_ind2)]/2;
        end
    end
end  

on_mdl = fitlm(on_time_vec1,on_time_vec2,'Intercept',false);
on1_axis = linspace(0,50);
% on2_pd = predict(on_mdl,on1_axis');
on2_pd_full = predict(on_mdl,on_time_vec1');

on_time_fig = figure('Position',[100 100 428 420]);
hold on
cmap = brewermap(9,'Set2');
% colormap(cmap)
% both on
scatter(on_time_vec1,on_time_vec2,markerSize,'o','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha,'MarkerEdgeAlpha',edgeAlpha)
plot(on1_axis, on1_axis, '-.','Color', 'k', 'LineWidth', 1.5)

% xlim([0 25])
% ylim([0 25])
set(gca,'Fontsize',14)
xlabel('spot 1 ON time (minutes)')
ylabel('spot 2 ON time (minutes)')
box on
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
saveas(on_time_fig,[figPath 'on_time_scatter.pdf'])


on_var_prior = 0.5*mean((on_time_vec1 - randsample(on_time_vec2,length(on_time_vec2),false)).^2);
% on_var_post = on_mdl.RMSE^2;
on_var_post = 0.5*mean((on_time_vec1-on_time_vec2).^2);
% on_var_post_alt = 0.5*mean((on2_pd_full'-on_time_vec2).^2)
info_gain_on = 0.5*log(on_var_prior/on_var_post) * log2(exp(1));


%% Off times
rng(983);

nucleus_id_vec_off = nucleus_id_vec(off_flags);%&long_flags);
off_indices = find(off_flags);%&long_flags);
nucleus_id_index_off = unique(nucleus_id_vec(off_flags));%&long_flags));
off_time_vec1 = [];
off_time_vec2 = [];
ap_off_vec = [];
nc_id_array = [];
for i = 1:length(nucleus_id_index)
    if sum(nucleus_id_vec_off==nucleus_id_index(i)) == 2
        off_inds = randsample(off_indices(nucleus_id_vec_off==nucleus_id_index(i)),2,false);
        f_ind_vec1 = find(~isnan(masterSet(off_inds(1)).fluo));
        f_ind_vec2 = find(~isnan(masterSet(off_inds(2)).fluo));
        if sum(ismember(f_ind_vec1,f_ind_vec2))>=minOverlap                        
            off_time_vec1(end+1) = masterSet(off_inds(1)).time(f_ind_vec1(end))/60;            
            off_time_vec2(end+1) = masterSet(off_inds(2)).time(f_ind_vec2(end))/60;
%             ap_off_vec(end+1) = nanmean(masterSet(off_inds(1)).APPosNucleus(1));
%             nc_id_array = vertcat(nc_id_array, off_inds);
        end
    end
end  

off_mdl = fitlm(off_time_vec1,off_time_vec2);
off1_axis = linspace(0,50);
off2_pd = predict(off_mdl,off1_axis');
off2_pd_full = predict(off_mdl,off_time_vec1');

off_time_fig = figure('Position',[100 100 428 420]);
hold on
cmap = brewermap(9,'Set2');
% colormap(cmap)
% both on
scatter(off_time_vec1,off_time_vec2,markerSize,'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha,'MarkerEdgeAlpha',edgeAlpha)
plot(off1_axis, off1_axis, '-.', 'Color', 'k', 'LineWidth', 2)

xlim([0 50])
ylim([0 50])
set(gca,'Fontsize',14)
xlabel('spot 1 off time (minutes)')
ylabel('spot 2 off time (minutes)')
box on
% legend([p1 p2 p3], 'both (p^2)','either (2(1-p))','neither ((1-p)^2))','Location','northwest');%'Orientation','horizontal')
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% h = colorbar;
% ylabel(h,'% embryo length')
off_time_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(off_time_fig,[figPath 'off_time_scatter.png'])
saveas(off_time_fig,[figPath 'off_time_scatter.pdf'])

% calculate initial variance and residual
off_var_prior = 0.5*mean((off_time_vec1-randsample(off_time_vec2,length(off_time_vec2),false)).^2);%var([off_time_vec1 off_time_vec2]);
% off_var_post = off_mdl.RMSE^2;
off_var_post = 0.5*mean((off_time_vec1-off_time_vec2).^2);
info_gain_off = 0.5*log(off_var_prior/off_var_post) * log2(exp(1))

%% average fluo


close all
nucleus_id_vec_fluo = nucleus_id_vec;%(long_flags&on_flags&off_flags);
fluo_indices = 1:length(nucleus_id_vec);%find(long_flags&on_flags&off_flags);
nucleus_id_index_fluo = unique(nucleus_id_vec);%(long_flags&on_flags&off_flags));
fluo_vec1 = [];
fluo_vec2 = [];
ap_vec_fluo = [];
rng(451);

for i = 1:length(nucleus_id_index)
    if sum(nucleus_id_vec_fluo==nucleus_id_index(i)) == 2
        fluo_inds = randsample(fluo_indices(nucleus_id_vec_fluo==nucleus_id_index(i)),2,false);
        f_ind_vec1 = find(~isnan(masterSet(fluo_inds(1)).fluo));
        f_ind_vec2 = find(~isnan(masterSet(fluo_inds(2)).fluo));
        
        % calculate indices where there is overlap
        overlap_times1 = masterSet(fluo_inds(1)).time(ismember(f_ind_vec1,f_ind_vec2));
        overlap_times2 = masterSet(fluo_inds(2)).time(ismember(f_ind_vec2,f_ind_vec1));
        if length(overlap_times1) >= 2*minOverlap
            % translate to interpolated time
            t_filter1 = find(masterSet(fluo_inds(1)).timeInterp >=overlap_times1(1) & masterSet(fluo_inds(1)).timeInterp <=overlap_times1(end));
            fluo_vec1(end+1) = nanmean(masterSet(fluo_inds(1)).fluoInterp(t_filter1)/1e5);

            t_filter2 = find(masterSet(fluo_inds(2)).timeInterp >=overlap_times2(1) & masterSet(fluo_inds(2)).timeInterp <=overlap_times2(end));
            fluo_vec2(end+1) = nanmean(masterSet(fluo_inds(2)).fluoInterp(t_filter2)/1e5);

%             ap_vec_fluo(end+1) = nanmean(masterSet(fluo_inds(2)).APPosParticle);
        end
                
    end
end  

fluo_vec1 = fluo_vec1/nanmean([fluo_vec1 fluo_vec2]);
fluo_vec2 = fluo_vec2/nanmean([fluo_vec1 fluo_vec2]);

fluo_mdl = fitlm(fluo_vec1,fluo_vec2);
flim = 4;
fluo_axis = linspace(0,flim);
fluo2_pd = predict(fluo_mdl,fluo_axis');
fluo2_pd_full = predict(fluo_mdl,fluo_vec1');

fluo_fig = figure('Position',[100 100 428 420]);
hold on
cmap = brewermap(9,'Set2');
% colormap(cmap)
% both on
scatter(fluo_vec1,fluo_vec2,markerSize,'o','MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k','MarkerFaceAlpha',markerAlpha,'MarkerEdgeAlpha',edgeAlpha)
plot(fluo_axis, fluo_axis,'-.','Color', 'k', 'LineWidth', 1.5)

xlim([0 flim])
ylim([0 flim])
set(gca,'Fontsize',14)
xlabel('spot 1 intensity (au)')
ylabel('spot 2 intensity (au)')
box on
% legend([p1 p2 p3], 'both (p^2)','either (2(1-p))','neither ((1-p)^2))','Location','northwest');%'Orientation','horizontal')
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% h = colorbar;
% ylabel(h,'% embryo length')
fluo_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(fluo_fig,[figPath 'fluo_scatter.png'])
saveas(fluo_fig,[figPath 'fluo_scatter.pdf'])

% calculate initial variance and residual
fluo_var_prior = 0.5*nanmean((fluo_vec1-randsample(fluo_vec2,length(fluo_vec2))).^2);%^nanvar(fluo_vec2);
fluo_var_post = 0.5*nanmean((fluo_vec1-fluo_vec2).^2);
% fluo_var_post = fluo_mdl.RMSE^2;
info_gain = 0.5*log2(fluo_var_prior/fluo_var_post) * log2(exp(1))
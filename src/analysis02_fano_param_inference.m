clear
close all

% make figure directory
figPath = '../fig/fluo_noise_analysis/';
mkdir(figPath)

% load master data set
load('../data/masterSet.mat')

%% generate master arrays to track sister spot activity across time and space
time_bounds = [3 40]*60;
minDP = 10;
minOverlap = 5;

tres = masterSet(1).tresInterp;
masterTimeVec = 0:tres:60*60;

nc_id_vec = [masterSet.ncID];
nc_id_index = unique(nc_id_vec);
analysis_flag_vec = false(size(nc_id_index));

for i = 1:length(nc_id_index)
    nc_indices = find(nc_id_vec==nc_id_index(i));
    if length(nc_indices) == 2
        time_vec = masterSet(nc_indices(1)).time;
        f1_vec = masterSet(nc_indices(1)).fluo;
        f2_vec = masterSet(nc_indices(2)).fluo;
        analysis_flag_vec(i) = time_vec(1)<=time_bounds(1) && time_vec(end)>=time_bounds(2) && ...
                                sum(~isnan(f1_vec) & ~isnan(f2_vec))>=minOverlap && ...
                                sum(~isnan(f1_vec)) >= minDP && sum(~isnan(f2_vec)) >= minDP;
    end
end  

nc_id_index2 = nc_id_index(analysis_flag_vec);
fluo_array = NaN(length(masterTimeVec),length(nc_id_index2),2);
ap_array = NaN(length(masterTimeVec),length(nc_id_index2),1);
time_array = NaN(length(masterTimeVec),length(nc_id_index2),1);
rng(312)
for n = 1:length(nc_id_index2)
    nc_indices = randsample(find(nc_id_vec==nc_id_index2(n)),2,false);
    % extract values
    t1_vec = masterSet(nc_indices(1)).timeInterp;
    t2_vec = masterSet(nc_indices(2)).timeInterp;    
    f1_vec = masterSet(nc_indices(1)).fluoInterp;
    f2_vec = masterSet(nc_indices(2)).fluoInterp;
    ap_vec = masterSet(nc_indices(1)).APPosNucleusInterp;
    % combine AP-related fields
    [t_vec_i,ia,~] = intersect(t1_vec,t2_vec);
    % add to arrays
    fluo_array(ismember(masterTimeVec,t1_vec),n,1) = f1_vec;
    fluo_array(ismember(masterTimeVec,t2_vec),n,2) = f2_vec;
    ap_array(ismember(masterTimeVec,t_vec_i),n,1) = ap_vec(ia);
    time_array(ismember(masterTimeVec,t1_vec),n,1) = masterTimeVec(ismember(masterTimeVec,t1_vec));
    time_array(ismember(masterTimeVec,t2_vec),n,2) = masterTimeVec(ismember(masterTimeVec,t2_vec));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Look at Fano Factor

% For this analysis, we will group traces by fluorescence and align them by
% their start, the initial object is to examine how the Fano Factor for
% each group evolves over time

max_fano_time = 25*60; % relative to start of trace
fano_time_vec = 0:tres:max_fano_time;
start_lag = 1;
% look for nuclei where traces are long enough and
fano_analysis_flag_vec = false(size(nc_id_index));
mean_fluo_vec = NaN(size(nc_id_index));

for i = 1:length(nc_id_index)
    nc_indices = find(nc_id_vec==nc_id_index(i));
    if length(nc_indices) == 2
        t1_vec = masterSet(nc_indices(1)).timeInterp(start_lag:end);
        t2_vec = masterSet(nc_indices(2)).timeInterp(start_lag:end);
        f1_vec = masterSet(nc_indices(1)).fluoInterp;
        f2_vec = masterSet(nc_indices(2)).fluoInterp;
        o12_vec = ismember(t1_vec,t2_vec);
        o21_vec = ismember(t2_vec,t1_vec);
        n_overlap = sum(o12_vec);
        qc1_flag = masterSet(nc_indices(1)).TraceQCFlag;
        qc2_flag = masterSet(nc_indices(2)).TraceQCFlag;
        
        mean_fluo_vec(i) = mean([f1_vec(o12_vec) f2_vec(o21_vec)]);
        
        fano_analysis_flag_vec(i) = n_overlap >= length(fano_time_vec) & qc1_flag & qc2_flag;
%         fano_analysis_flag_vec(i) = length(f1_vec) >= length(fano_time_vec) & length(f2_vec) >= length(fano_time_vec) ...
%                                     & qc1_flag & qc2_flag ...
%                                     & abs(t1_vec(1)-t2_vec(1)) <= 300;
    end
end  
%% Divide elligible traces into inference groups
nc_id_index_fano = nc_id_index(fano_analysis_flag_vec);
mean_fluo_vec_fano = mean_fluo_vec(fano_analysis_flag_vec);
f_norm_factor = 10^4;
% assign to groups
n_fluo_groups = 9;
q_vec = linspace(0,1,n_fluo_groups+1);
fluo_bins = quantile(mean_fluo_vec_fano,q_vec);
fluo_bin_vec = discretize(mean_fluo_vec_fano, fluo_bins);

% generate inference arrays
inference_struct = struct;

% rng(312)
for f = 1:n_fluo_groups
    fluo_filter = fluo_bin_vec==f;
    fluo_indices = find(fluo_filter);
%     delta_array = NaN(length(fano_time_vec),sum(fluo_filter));
    trace_array = NaN(length(fano_time_vec),sum(fluo_filter),2);
    mean_array = NaN(length(fano_time_vec),sum(fluo_filter));

    for n = 1:length(fluo_indices)
        nc_indices = randsample(find(nc_id_vec==nc_id_index_fano(fluo_indices(n))),2,false);
        % extract values
        t1_vec = masterSet(nc_indices(1)).timeInterp(start_lag:end);
        t2_vec = masterSet(nc_indices(2)).timeInterp(start_lag:end);    
        f1_vec = masterSet(nc_indices(1)).fluoInterp(start_lag:end)/f_norm_factor;
        f2_vec = masterSet(nc_indices(2)).fluoInterp(start_lag:end)/f_norm_factor;
        
        % combine fields
        [t_vec_i,i12,i21] = intersect(t1_vec,t2_vec);
        f1_int = f1_vec(i12);
        f1_int = cumsum(f1_int(1:length(fano_time_vec)))'/140*15;
        f2_int = f2_vec(i21);
        f2_int = cumsum(f2_int(1:length(fano_time_vec)))'/140*15;
        
        % add to arrays
%         mean_array(:,n) = mean([f1_int f2_int],2);
        trace_array(:,n,:) = cat(3,f1_int,f2_int);
%         delta_array(:,n) = (f1_int-f2_int).^2;
%         total_var_array(:,n) = 0.5*mean((f1_int-randsample(f2_int,length(f2_int),false)).^2,2);        
    end
    trace_array_norm = trace_array - nanmean(trace_array,2);
%     rs_vec = randsample(1:length(fluo_indices),length(fluo_indices),false);
%     total_var_vec = 0.5*nanmean((trace_array(:,:,1) - trace_array(:,rs_vec,2)).^2,2);
%     total_var_vec = var([trace_array_norm(:,:,1) trace_array_norm(:,:,2)],[],2);    
    delta_array = (trace_array_norm(:,:,1)-trace_array_norm(:,:,2)).^2;
    
    inference_struct(f).mean_array = mean(trace_array,3);
%     inference_struct(f).mean_vec = nanmean(mean_array,2);
    inference_struct(f).delta_array = delta_array;
%     inference_struct(f).intrinsic_var_vec = 0.5*nanmean(delta_array,2);
%     inference_struct(f).fano_array = inference_struct(f).intrinsic_var_vec./nanmean(mean_array,2);
%     inference_struct(f).fano_vec = nanmean(delta_array./mean_array,2);
%     inference_struct(f).total_var_vec = total_var_vec;
    inference_struct(f).trace_array = trace_array;
end

% get bootstrap estimates of the end-point fano factor
nBoots = 100;
fano_array = NaN(length(fano_time_vec),n_fluo_groups,nBoots);
int_var_array = NaN(length(fano_time_vec),n_fluo_groups,nBoots);
fluo_array = NaN(length(fano_time_vec),n_fluo_groups,nBoots);
for f = 1:n_fluo_groups
    % extract arrays    
    mean_array = inference_struct(f).mean_array;
    delta_array = inference_struct(f).delta_array;
    % do bootstrap sampling
    option_vec = 1:size(mean_array,2);
    for n = 1:nBoots
        boot_indices = randsample(option_vec,length(option_vec),true);
        fano_array(:,f,n) = 0.5*nanmean(delta_array(:,boot_indices),2)./nanmean(mean_array(:,boot_indices),2);
        int_var_array(:,f,n) = 0.5*nanmean(delta_array(:,boot_indices),2);
        fluo_array(:,f,n) = nanmean(mean_array(:,boot_indices),2);
    end
end

fluo_mean_array = nanmean(fluo_array,3)./(1:length(fano_time_vec))';%*60;
fluo_ste_array = nanstd(fluo_array,[],3)./(1:length(fano_time_vec))';%
int_mean_array = nanmean(int_var_array,3)./(1:length(fano_time_vec))';%
int_ste_array = nanstd(int_var_array,[],3)./(1:length(fano_time_vec))';%
fano_mean_array = nanmean(fano_array,3);
fano_ste_array = nanstd(fano_array,[],3);

% Make figures
close all
fano_time_fig = figure;
hold on
cmap = flipud(brewermap([],'Spectral'));
colormap(cmap)
skip_vec = 1:3:length(fano_time_vec);
for f = 1:n_fluo_groups
%     errorbar(fluo_mean_array(skip_vec,f),fano_mean_array(skip_vec,f),fano_ste_array(skip_vec,f),'Color','k');
    plot(fluo_mean_array(skip_vec,f),fano_mean_array(skip_vec,f),'-k')
    scatter(fluo_mean_array(skip_vec,f),fano_mean_array(skip_vec,f),[],fano_time_vec(skip_vec)/60,'filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.8);
end    

xlabel('average spot intensity (au)')
ylabel('Fano Factor')
h = colorbar;
ylabel(h,'time (minutes)')
set(gca,'Fontsize',14)
% xlim([0 1])
set(gca,'Color',[228,221,209]/255) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

fano_time_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(fano_time_fig,[figPath 'fano_factor_over_time.png'])
saveas(fano_time_fig,[figPath 'fano_factor_over_time.pdf'])

%% Fano vs fluo
fub = fluo_ste_array;
flb = fluo_ste_array;
iub = int_ste_array;
ilb = int_ste_array;

int_vs_fluo = figure;
cmap = brewermap([],'Set2');
hold on
errorbar(fluo_mean_array(end,:),int_mean_array(end,:),ilb(end,:),iub(end,:),flb(end,:),fub(end,:),'.','Color','k','Capsize',0)
scatter(fluo_mean_array(end,:),int_mean_array(end,:),75,'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k')

xlabel('average spot intensity (au)')
ylabel('intrinsic noise (\delta_{int}^2)')

set(gca,'Fontsize',14)
% xlim([0 20])
% ylim([0 140])
set(gca,'Color',[228,221,209]/255) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

int_vs_fluo.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(int_vs_fluo,[figPath 'int_vs_fluo.png'])
saveas(int_vs_fluo,[figPath 'int_vs_fluo.pdf'])

%% Now let's see if we can recapitulate this trend with the results uncovered in eLife paper
T = 140/60; % approximate elongation time
npt = 1e3;
% koff_factor = 3;
eff2_factor = 4; % this derives from observation in Lammers2020 that equivalent 2 state system tends to have ~3-4 fold higher off rate
% bound kon and r values by calculating values at endpoints
kon_min = 0;
kon_max = 1.247;
kon_vec = linspace(kon_min,kon_max,npt);
koff_min = 0.84;
koff_max = 0.56;
koff_vec_raw = linspace(koff_min,koff_max,npt);
koff_vec = koff_vec_raw*eff2_factor;
r_min = 0;
r_max = 13.0;
r_vec_raw = linspace(r_min,r_max,npt);
r_vec = r_vec_raw.*(kon_vec+koff_vec)./(kon_vec+koff_vec_raw);


% pd_fun_mr = @(params) mean_rate_fun(kon_vec, params(1)*koff_vec, params(2)*r_vec, T);
% pd_fun_in = @(params) intrinsic_noise_fun(kon_vec, params(1)*koff_vec, params(2)*r_vec, T);
mean_trend_pd = mean_rate_fun(kon_vec, koff_vec, r_vec, T);
mean_trend_pd_raw = mean_rate_fun(kon_vec, koff_vec_raw, r_vec_raw, T);
int_trend_pd_raw = intrinsic_noise_fun(kon_vec, koff_vec_raw, r_vec_raw, T);
int_trend_pd = intrinsic_noise_fun(kon_vec, koff_vec, r_vec, T);

close all

int_vs_fluo = figure;
cmap = brewermap([],'Set2');
hold on
p1 = plot(mean_trend_pd_raw,int_trend_pd_raw,'-','Color','k','LineWidth',2);
p2 = plot(mean_trend_pd,int_trend_pd,'-.','Color','k','LineWidth',2);
errorbar(fluo_mean_array(end,:),int_mean_array(end,:),ilb(end,:),iub(end,:),flb(end,:),fub(end,:),'.','Color','k','Capsize',0)
s = scatter(fluo_mean_array(end,:),int_mean_array(end,:),75,'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k');

xlabel('average spot intensity (au)')
ylabel('intrinsic noise (\delta_{int}^2 f)')

set(gca,'Fontsize',14)
xlim([0 2.5])
ylim([0 22])
% ylim([0 140])
legend([p1 p2 s],'prediction (raw)','prediction (adjusted)','2 spot data','Location','northwest')
set(gca,'Color',[228,221,209]/255) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

int_vs_fluo.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(int_vs_fluo,[figPath 'int_vs_fluo_pd.png'])
saveas(int_vs_fluo,[figPath 'int_vs_fluo_pd.pdf'])

% %% calculate total and intrinsic variance for qualifying nuclei as a function of space and time
% time_window = 1.5*60;
% ap_window = 2.5;
% min_nc = 20;
% 
% ap_axis = 22.5:1:92.5;
% time_axis = (12.5:1:47.5)*60;
% 
% total_variance_array = NaN(length(time_axis),length(ap_axis));
% intrinsic_variance_array = NaN(length(time_axis),length(ap_axis));
% mean_array = NaN(length(time_axis),length(ap_axis));
% 
% for a = 1:length(ap_axis)
%     apc = ap_axis(a);
%     ap_filter = ap_array>=apc-ap_window & ap_array<=apc+ap_window;
%     for t = 1:length(time_axis)       
%        tc = time_axis(t);
%        ref_n = sum(masterTimeVec >= tc-time_window & masterTimeVec <= tc+time_window);
%        time_filter = time_array >= tc-time_window & time_array <= tc+time_window;
%        
%        % combine filters
%        at_filter_full = ap_filter & time_filter(:,:,1) & time_filter(:,:,2);
%        at_filter = nansum(ap_filter & time_filter(:,:,1) & time_filter(:,:,2),1);
%        f_vec = nansum(fluo_array.*repmat(ap_filter,1,1,2),1);
%        f1_vec = f_vec(:,:,1);
%        f2_vec = f_vec(:,:,2);
%        f1_mean_vec = f1_vec(at_filter==ref_n)/ref_n;
%        f2_mean_vec = f2_vec(at_filter==ref_n)/ref_n;
%        
%        % add values to array
%        if length(f1_mean_vec) >= min_nc
%            total_variance_array(t,a) = 0.5*mean((f1_mean_vec-randsample(f2_mean_vec,length(f2_mean_vec),false)).^2);
%            intrinsic_variance_array(t,a) = 0.5*mean((f1_mean_vec-f2_mean_vec).^2);
%            mean_array(t,a) = mean([f1_mean_vec f2_mean_vec]);
%        end
%     end
% end    

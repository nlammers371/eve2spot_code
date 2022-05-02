clear
close all

% make figure directory
figPath = '../fig/fluo_information/';
mkdir(figPath)

% load master data set
load('../data/masterSet.mat')
load('../data/trace_fit_struct.mat')


%% calculate total and intrinsic variance for qualitfying nuclei as a function of space and time
ap_array = trace_fit_struct.ap_array;
time_array = trace_fit_struct.time_array;
r_array = trace_fit_struct.r_array ./ nanmean(trace_fit_struct.r_array(:));
fluo_array = trace_fit_struct.fluo_array;
masterTimeVec = 0:20:60*60;
ap_window = 2.1; % approximately 1 cells
min_nc = 20;

% generate transcription rate arrays with various levels of temporal
% averaging
t_window_size_vec = [1:2:41];%
dt_vec = (t_window_size_vec)*20/60;
r_mean_array = NaN(size(r_array,1),size(r_array,2),size(r_array,3),length(t_window_size_vec));
for t = 1:length(t_window_size_vec)
    r_mean_array(:,:,:,t) = movmean(r_array,t_window_size_vec(t),1,'omitnan');
end
% this sets resolution of grid 
ap_axis = 22.5:1:92.5;
time_axis = masterTimeVec;

total_variance_array = NaN(length(time_axis),length(ap_axis),length(t_window_size_vec));
intrinsic_variance_array = NaN(length(time_axis),length(ap_axis),length(t_window_size_vec));
mean_array = NaN(length(time_axis),length(ap_axis),length(t_window_size_vec));

for a = 1:length(ap_axis)
    apc = ap_axis(a);
    ap_filter = ap_array>=apc-ap_window & ap_array<=apc+ap_window;
    
    for w = 1:length(t_window_size_vec) 
       % combine filters
%        at_filter_full = ap_filter & time_filter(:,:,1) & time_filter(:,:,2);
%        at_filter = nansum(ap_filter & time_filter(:,:,1) & time_filter(:,:,2),1);       
       r_array_w = r_mean_array(:,:,:,w);
       r_array_w_1 = r_array_w(:,:,1);
       r_array_w_2 = r_array_w(:,:,2);
       nn_filter = ~isnan(all(r_array_w,3));
       full_filter = ap_filter & nn_filter;
       r_array_w_1(~full_filter) = NaN;
       r_array_w_2(~full_filter) = NaN;
       
       % calculate statisticst row-wise
       diff_array = (r_array_w_2-r_array_w_1).^2;
       int_vec = 0.5*nanmean(diff_array,2);
       tot_vec = 0.5*nanvar(r_array_w_1,[],2) + 0.5*nanvar(r_array_w_2,[],2);
       mean_vec = 0.5*nanmean(r_array_w_1,2) + 0.5*nanmean(r_array_w_2,2);
       n_nc_vec = sum(full_filter,2);
       
       % add to array when appropriate
       total_variance_array(n_nc_vec>=min_nc,a,w) = tot_vec(n_nc_vec>=min_nc);
       intrinsic_variance_array(n_nc_vec>=min_nc,a,w) = int_vec(n_nc_vec>=min_nc);
       mean_array(n_nc_vec>=min_nc,a,w) = mean_vec(n_nc_vec>=min_nc);
      
    end
end    
%% simple fits to estimate time constant for Information accumulation
% Assume this is constant accross data set
% I(t) = I0 (1-exp(-t/tau))
% d I(t) / dt = I0/tau exp(-t/tau)
% d I(0) / dt = I0/tau
i_array_raw = log2(exp(1)).*log2(total_variance_array./intrinsic_variance_array);
i_array_raw(i_array_raw<0) = 0;
avg_I_vs_T = permute(nanmean(nanmean(i_array_raw,1),2),[3,2,1]);

i_fun = @(x) x(1) * (1 - exp(-dt_vec./x(2)));
ob_fun = @(x) avg_I_vs_T' - i_fun(x);
x_fit = lsqnonlin(ob_fun,[1 1],[0 0],[Inf Inf]);

% make plot
avg_fig = figure;
cmap = brewermap([],'Set2');
hold on

plot(dt_vec,avg_I_vs_T,'Color','k','LineWidth',3)
plot(dt_vec,i_fun(x_fit),'Color', cmap(2,:) ,'LineWidth',3)

set(gca,'Fontsize',14)
ylabel('Information (bits)')
xlabel('averaging window size (minutes) ')
legend('data', 'exponential fit','Location','southeast')
% xlabel('minutes into nc14')
set(gca,'Color',[228,221,209]/255) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

avg_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(avg_fig,[figPath 'averaging_fig.png'])
saveas(avg_fig,[figPath 'averaging_fig.pdf'])
%%

% estimate information rate and total information content
x_lim = [5 65];
y_lim = [12*3 50*3];
% For this analysis, we will group traces by fluorescence and align them by
% their start, the initial object is to examine how the Fano Factor for
% each group evolves over time

[~, plot_ind] = min(abs(dt_vec-x_fit(2)));%16;
info_array = i_array_raw(:,:,plot_ind);

% filtWidth = 5;
% filtSigma = 3;
% imageFilter=fspecial('gaussian',filtWidth,filtSigma);
% info_array_sm = nanconv(info_array,imageFilter, 'nanout');

% obtain smoothed estimate of raw information levels
n_reps = 100;
filterSigma = 1;
noise_mean = nanmean(info_array(:));
noise_sigma = nanstd(info_array(:));
sm_array = NaN(size(info_array,1),size(info_array,2),n_reps);
n_nan = sum(isnan(info_array(:)));
for r = 1:n_reps
    rand_frames = normrnd(noise_mean,noise_sigma,1,n_nan);
    slice = NaN(size(info_array));
    slice(~isnan(info_array)) = info_array(~isnan(info_array));
    slice(isnan(info_array)) = rand_frames;
    sm_array(:,:,r) = imgaussfilt(slice,filterSigma);
end    
info_sm_mean = nanmean(sm_array,3);
info_sm_mean(isnan(info_array)) = NaN;
  
% now calculate I0 and I0/dt
decay_time = 100;
I0_array = info_sm_mean/(1-exp(-1)) * (1-exp(-decay_time/x_fit(2)));
I_rate_array = I0_array / x_fit(2);

close all
info_fig = figure;
hold on
cmap = flipud(brewermap([],'Spectral'));
colormap(cmap);
p = pcolor(flipud(I_rate_array));
p.EdgeAlpha = 0.05;

colorbar
        
xlabel('AP position (%embryo length)')
ylabel('minutes into nc14')
set(gca,'Color',[228,221,209]/255) 
h = colorbar;
ylabel(h,'dI/dt (bits per minute)')

set(gca,'Fontsize',14)

box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
caxis([0 1]);
ylim(y_lim)
xlim(x_lim)
% caxis([4 6])
xt = get(gca,'xtick');
set(gca,'xticklabels',ap_axis(xt));
yt = get(gca,'ytick');
tflip = fliplr(masterTimeVec);
set(gca,'yticklabels',round(tflip(yt)/60));
info_fig.InvertHardcopy = 'off';
set(gcf,'color','w');


saveas(info_fig,[figPath 'info_fig_rate.png'])
saveas(info_fig,[figPath 'info_fig_rate.pdf'])

p = pcolor(flipud(I0_array));
p.EdgeAlpha = 0.05;

ylabel(h,'I_0 (bits)')
caxis([0 4]);
saveas(info_fig,[figPath 'info_fig.png'])
saveas(info_fig,[figPath 'info_fig.pdf'])
%%

mf_array = mean_array(:,:,plot_ind);

mf_mean = nanmean(mf_array(:));
mf_sigma = nanstd(mf_array(:));
sm_array = NaN(size(mf_array,1),size(mf_array,2),n_reps);
n_nan = sum(isnan(mf_array(:)));
for r = 1:n_reps
    rand_frames = normrnd(mf_mean,mf_sigma,1,n_nan);
    slice = NaN(size(mf_array));
    slice(~isnan(mf_array)) = mf_array(~isnan(mf_array));
    slice(isnan(mf_array)) = rand_frames;
    sm_array(:,:,r) = imgaussfilt(slice,filterSigma);
end    
mf_sm_mean = nanmean(sm_array,3);
mf_sm_mean(isnan(info_array)) = NaN;

mean_rate_fig = figure;
cmap = flipud(brewermap([],'Spectral'));
colormap(cmap);
p = pcolor(flipud(mf_sm_mean*3));
p.EdgeAlpha = 0.05;

colorbar
        
xlabel('AP position (%embryo length)')
ylabel('minutes into nc14')
set(gca,'Color',[228,221,209]/255) 
h = colorbar;
ylabel(h,'transcription rate (au/min)')

set(gca,'Fontsize',14)

box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ylim(y_lim)
xlim(x_lim)

xt = get(gca,'xtick');
set(gca,'xticklabels',ap_axis(xt));
yt = get(gca,'ytick');
tflip = fliplr(masterTimeVec);
set(gca,'yticklabels',round(tflip(yt)/60));
mean_rate_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(mean_rate_fig,[figPath 'rate_fig.png'])
saveas(mean_rate_fig,[figPath 'rate_fig.pdf'])

%% Total variance
close all
total_var_array = total_variance_array(:,:,plot_ind);

tot_var_fig = figure;%('Position',[100 100 768 256]);
cmap = (brewermap([],'Reds'));
colormap(cmap);
p = pcolor(flipud(total_var_array*3));
p.EdgeAlpha = 0.05;

colorbar
        
xlabel('AP position (%embryo length)')
ylabel('minutes into nc14')
set(gca,'Color',[228,221,209]/255) 
h = colorbar;
ylabel(h,'total variance (au)')

set(gca,'Fontsize',14)

box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ylim(y_lim)
xlim(x_lim)
caxis([0 4])
xt = get(gca,'xtick');
set(gca,'xticklabels',ap_axis(xt));
yt = get(gca,'ytick');
tflip = fliplr(masterTimeVec);
set(gca,'yticklabels',round(tflip(yt)/60));
tot_var_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(tot_var_fig,[figPath 'tot_var_fig.png'])
saveas(tot_var_fig,[figPath 'tot_var_fig.pdf'])


int_var_array = intrinsic_variance_array(:,:,plot_ind);

int_var_fig = figure;%('Position',[100 100 768 256]);
cmap = flipud(brewermap([],'Spectral'));
colormap(cmap);
p = pcolor(flipud(int_var_array*3));
p.EdgeAlpha = 0.05;

colorbar
        
xlabel('AP position (%embryo length)')
ylabel('minutes into nc14')
set(gca,'Color',[228,221,209]/255) 
h = colorbar;
ylabel(h,'intrinsic variance (au)')

set(gca,'Fontsize',14)

box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
caxis([0 4])
ylim(y_lim)
xlim(x_lim)

xt = get(gca,'xtick');
set(gca,'xticklabels',ap_axis(xt));
yt = get(gca,'ytick');
tflip = fliplr(masterTimeVec);
set(gca,'yticklabels',round(tflip(yt)/60));
int_var_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

saveas(int_var_fig,[figPath 'int_var_fig.png'])
saveas(int_var_fig,[figPath 'int_var_fig.pdf'])

%% Make example decoding figure

plot_ind = 10;

trace_fig = figure;%('Position',[100 100 768 256]);
cmap_s2 = flipud(brewermap([],'Set2'));

plot(masterTimeVec/60,fluo_array(:,plot_ind,2,1)/1e4,'Color','k','LineWidth',2)

set(gca,'Fontsize',14)

ylabel('spot intensity (f)')
xlabel('minutes into nc14')
set(gca,'Color',[228,221,209]/255) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';


trace_fig.InvertHardcopy = 'off';
% set(gcf,'color','w');

yyaxis right

stairs(masterTimeVec/60,r_array(:,plot_ind,2,1),'Color',cmap_s2(4,:),'LineWidth',2)

% set(gca,'Fontsize',14)

ylabel('transcription rate (r)')
% xlabel('minutes into nc14')
set(gca,'Color',[228,221,209]/255) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cmap_s2(4,:);
ax.XAxis(1).Color = 'k';

ylim([0 10])
trace_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(trace_fig,[figPath 'trace_fig.png'])
saveas(trace_fig,[figPath 'trace_fig.pdf'])
%%
% 
% time_window = 1*60;
% ap_window = 2.1;
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
%        f_vec = nansum(r_array.*repmat(ap_filter,1,1,2),1);
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
% 
% %%
% info_array = log2(exp(1)*total_variance_array./intrinsic_variance_array);
% figure;
% imagesc(info_array)
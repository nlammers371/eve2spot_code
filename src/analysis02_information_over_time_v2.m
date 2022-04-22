clear
close all

% make figure directory
figPath = '../fig/fluo_noise_analysis/';
mkdir(figPath)

% load master data set
load('../data/masterSet.mat')
load('../data/trace_fit_struct.mat')


%% calculate total and intrinsic variance for qualitfying nuclei as a function of space and time
ap_array = trace_fit_struct.ap_array;
time_array = trace_fit_struct.time_array;
r_array = trace_fit_struct.r_array;
fluo_array = trace_fit_struct.fluo_array;
masterTimeVec = 0:20:60*60;
time_window = 1.5*60;
ap_window = 2.1;
min_nc = 20;

ap_axis = 22.5:1:92.5;
time_axis = (12.5:1:47.5)*60;

total_variance_array = NaN(length(time_axis),length(ap_axis));
intrinsic_variance_array = NaN(length(time_axis),length(ap_axis));
mean_array = NaN(length(time_axis),length(ap_axis));

for a = 1:length(ap_axis)
    apc = ap_axis(a);
    ap_filter = ap_array>=apc-ap_window & ap_array<=apc+ap_window;
    for t = 1:length(time_axis)       
       tc = time_axis(t);
       ref_n = sum(masterTimeVec >= tc-time_window & masterTimeVec <= tc+time_window);
       time_filter = time_array >= tc-time_window & time_array <= tc+time_window;
       
       % combine filters
       at_filter_full = ap_filter & time_filter(:,:,1) & time_filter(:,:,2);
       at_filter = nansum(ap_filter & time_filter(:,:,1) & time_filter(:,:,2),1);
       f_vec = nansum(r_array.*repmat(ap_filter,1,1,2),1);
       f1_vec = f_vec(:,:,1);
       f2_vec = f_vec(:,:,2);
       f1_mean_vec = f1_vec(at_filter==ref_n)/ref_n;
       f2_mean_vec = f2_vec(at_filter==ref_n)/ref_n;
       
       % add values to array
       if length(f1_mean_vec) >= min_nc
           total_variance_array(t,a) = 0.5*mean((f1_mean_vec-randsample(f2_mean_vec,length(f2_mean_vec),false)).^2);
           intrinsic_variance_array(t,a) = 0.5*mean((f1_mean_vec-f2_mean_vec).^2);
           mean_array(t,a) = mean([f1_mean_vec f2_mean_vec]);
       end
    end
end    

%% Look at Fano Factor

% For this analysis, we will group traces by fluorescence and align them by
% their start, the initial object is to examine how the Fano Factor for
% each group evolves over time

info_array = log2(exp(1)*total_variance_array./intrinsic_variance_array);
filtWidth = 3;
filtSigma = 1;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
info_array_sm = nanconv(info_array,imageFilter, 'nanout');

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
  
    
close all
figure;
imagesc(info_sm_mean)
colorbar
        
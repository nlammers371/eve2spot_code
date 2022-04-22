% Script to Simulate Process of analyzing 2 spot data
%----------------------Specify Parameters and Load Data-------------------%
clear 
close all
%System Characteristics
K = 2;
w = 10;
R = [-0.5 1; 0.5 -1]/60;
r_emission = [0 10/60];
pi0 = [0.5 0.5];
noise = r_emission(2)*5;
t_res = 15;
seq_length = 100;
n_traces = 200;
alpha = 2;
% simulate traces
simulation_set = struct;
for n = 1:n_traces
    gillespie = synthetic_rate_gillespie(seq_length, alpha, ...
                                K, w, R, t_res, r_emission, noise, pi0);
    simulation_set(n).fluo = gillespie.fluo_MS2;
    simulation_set(n).nucleus = ceil(n/2);
    simulation_set(n).t_res = t_res;
    simulation_set(n).r_true = r_emission;
    simulation_set(n).k_on_true = R(2,1);
    simulation_set(n).k_off_true = R(1,2);
end
%%% -----Calculate Cumulative Fluorescence & Examine Noise ----------------%


t_res = simulation_set(1).t_res;
%set starting time step
start = 1;
%Account for jitter in trace start, end times
for i = 1:length(simulation_set)
    f = simulation_set(i).fluo; 
    simulation_set(i).fluo_normed = [0 f(start:find(f,1,'last'))] ;
    simulation_set(i).cf_fluo_normed = cumsum(simulation_set(i).fluo_normed);
    simulation_set(i).fluo_shuffled = randsample(simulation_set(i).fluo_normed,length(simulation_set(i).fluo_normed),false);
    simulation_set(i).cf_fluo_shuffled = cumsum(simulation_set(i).fluo_shuffled);
    simulation_set(i).duration = length(simulation_set(i).fluo_normed);    
end
%Determine maximum span of eligibillity for each nucleus given durations of
%consituent traces
nuc_index = [simulation_set.nucleus];
nuc_set = unique(nuc_index);
global_minimum = min([simulation_set.duration]);
global_maximum = max([simulation_set.duration]);
%Store Cumulative fluorecence for each pair for plotting
scatter_pairs = zeros(length(nuc_set),2);
i_vec = 1:global_maximum;
for i = 1:length(nuc_set)
    n_ids = find(nuc_index==i);
    min_len = min([simulation_set(n_ids).duration]);
    n_grp = 1;    
    for j = n_ids
        simulation_set(j).max_ind = min_len;        
        cf = simulation_set(j).cf_fluo_normed;
        cf = cf(i_vec <= min_len);                
        scatter_pairs(i,n_grp) = simulation_set(j).cf_fluo_normed(global_minimum);
        n_grp = n_grp + 1;
    end
end

noise_fig = figure;
colormap('winter')
cm = colormap;
hold on
max_f = .9*nanmax(nanmax(scatter_pairs));
min_f = .9*nanmin(nanmin(scatter_pairs));
plot(linspace(min_f,max_f),linspace(min_f,max_f),'black')
scatter(scatter_pairs(:,1),scatter_pairs(:,2), 15,'MarkerFaceColor', cm(30,:),...
    'MarkerFaceAlpha',.4,'MarkerEdgeAlpha',0)
grid on
axis([min_f max_f min_f max_f])
title('Examining Intrinsic and Extrinsic Variability in Cumualtive Fluorescence')
xlabel('Locus 1 (AU)')
ylabel('Locus 2 (AU)')
grid on
% saveas(noise_fig, [figpath '/noise_scatter.png'],'png');

% cum_fig = figure;
% hold on 
% for i = 1:nc_set
%% -------------------Extract Intrinsic and Extrinsic Terms --------------%
% fluo_per_mRNA = simulation_set(1).fluo_per_mRNA;
A_array_cf_n = NaN(global_maximum,length(nuc_set));
B_array_cf_n = NaN(global_maximum,length(nuc_set));

A_array_cf_s = NaN(global_maximum,length(nuc_set));
B_array_cf_s = NaN(global_maximum,length(nuc_set));

A_array_ff_n = NaN(global_maximum,length(nuc_set));
B_array_ff_n = NaN(global_maximum,length(nuc_set));

A_array_ff_s = NaN(global_maximum,length(nuc_set));
B_array_ff_s = NaN(global_maximum,length(nuc_set));

for n = nuc_set
    ids = find(nuc_index==n);
    mi = simulation_set(ids(1)).max_ind;
    A_array_cf_n(1:mi,n) = simulation_set(ids(1)).cf_fluo_normed(1:mi);
    B_array_cf_n(1:mi,n) = simulation_set(ids(2)).cf_fluo_normed(1:mi);
    A_array_cf_s(1:mi,n) = simulation_set(ids(1)).cf_fluo_shuffled(1:mi);
    B_array_cf_s(1:mi,n) = simulation_set(ids(2)).cf_fluo_shuffled(1:mi);
    
    A_array_ff_n(1:mi,n) = simulation_set(ids(1)).fluo_normed(1:mi);
    B_array_ff_n(1:mi,n) = simulation_set(ids(2)).fluo_normed(1:mi);
    A_array_ff_s(1:mi,n) = simulation_set(ids(1)).fluo_shuffled(1:mi);
    B_array_ff_s(1:mi,n) = simulation_set(ids(2)).fluo_shuffled(1:mi);
end

% A_array_cf_n = A_array_cf_n ;% ./ repmat(adjust',1,size(A_array_cf,2));
% B_array_cf_n = B_array_cf_n ;%./ repmat(adjust',1,size(A_array_cf,2));
% A_array_ff_n = A_array_ff_n ;%./ repmat(adjust',1,size(A_array_cf,2));
% B_array_ff_n = B_array_ff_n ;%./ repmat(adjust',1,size(A_array_cf,2));

mean_A_cf_n = nanmean(A_array_cf_n,2);
mean_B_cf_n = nanmean(B_array_cf_n,2);
mean_A_ff_n = nanmean(A_array_ff_n,2);
mean_B_ff_n = nanmean(B_array_ff_n,2);

mean_A_cf_s = nanmean(A_array_cf_s,2);
mean_B_cf_s = nanmean(B_array_cf_s,2);
mean_A_ff_s = nanmean(A_array_ff_s,2);
mean_B_ff_s = nanmean(B_array_ff_s,2);

diff_A_cf_n = A_array_cf_n - repmat(mean_A_cf_n,1,length(nuc_set));
diff_B_cf_n = B_array_cf_n - repmat(mean_B_cf_n,1,length(nuc_set));
diff_A_ff_n = A_array_ff_n - repmat(mean_A_ff_n,1,length(nuc_set));
diff_B_ff_n = B_array_ff_n - repmat(mean_B_ff_n,1,length(nuc_set));

diff_A_cf_s = A_array_cf_s - repmat(mean_A_cf_s,1,length(nuc_set));
diff_B_cf_s = B_array_cf_s - repmat(mean_B_cf_s,1,length(nuc_set));
diff_A_ff_s = A_array_ff_s - repmat(mean_A_ff_s,1,length(nuc_set));
diff_B_ff_s = B_array_ff_s - repmat(mean_B_ff_s,1,length(nuc_set));

%%% --------------------Modeling Intrinsic Noise---------------------------%
time_vec = (0:(size(B_array_ff_n,1)-1))*t_res;
%Number of bootstrap samples
n_boots = 100;
%Size of bootstrap
boot_size = length(nuc_set);

sample_index = 1:length(nuc_set);

%Sim params
avg_rate = simulation_set(1).r_true(2) ;

if K == 2
    factor = 1;
else
    factor = 2;
end

%arrays to store results
fano_emp_cf_array_n = zeros(length(time_vec),n_boots);
var_emp_cf_array_n = zeros(length(time_vec),n_boots);
mean_emp_cf_array_n = zeros(length(time_vec),n_boots);

fano_emp_ff_array_n = zeros(length(time_vec),n_boots);
var_emp_ff_array_n = zeros(length(time_vec),n_boots);
mean_emp_ff_array_n = zeros(length(time_vec),n_boots);

fano_emp_cf_array_s = zeros(length(time_vec),n_boots);
var_emp_cf_array_s = zeros(length(time_vec),n_boots);
mean_emp_cf_array_s = zeros(length(time_vec),n_boots);

fano_emp_ff_array_s = zeros(length(time_vec),n_boots);
var_emp_ff_array_s = zeros(length(time_vec),n_boots);
mean_emp_ff_array_s = zeros(length(time_vec),n_boots);

for boot = 1:n_boots
    boot_sample = randsample(sample_index,boot_size,true);
    
    %Calculate intrinsic noise;
    int_vec_cf_n = .5*(nanmean((diff_A_cf_n(:,boot_sample)-diff_B_cf_n(:,boot_sample)).^2,2));
    mean_vec_cf_n = nanmean([A_array_cf_n(:,boot_sample) B_array_cf_n(:,boot_sample)],2);
    
    var_vec_ff_n = mean(.5*((diff_A_ff_n(:,boot_sample).^2+diff_B_ff_n(:,boot_sample).^2)));
    int_vec_ff_n = .5*(nanmean((diff_A_ff_n(:,boot_sample)-diff_B_ff_n(:,boot_sample)).^2,2));
    mean_vec_ff_n = nanmean([A_array_ff_n(:,boot_sample) B_array_ff_n(:,boot_sample)],2);
    
    int_vec_cf_s = .5*(nanmean((diff_A_cf_s(:,boot_sample)-diff_B_cf_s(:,boot_sample)).^2,2));
    mean_vec_cf_s = nanmean([A_array_cf_s(:,boot_sample) B_array_cf_s(:,boot_sample)],2);
    
    int_vec_ff_s = .5*(nanmean((diff_A_ff_s(:,boot_sample)-diff_B_ff_s(:,boot_sample)).^2,2));
    mean_vec_ff_s = nanmean([A_array_ff_s(:,boot_sample) B_array_ff_s(:,boot_sample)],2);    
    
    fano_emp_cf_array_n(:,boot) = int_vec_cf_n ./ mean_vec_cf_n;
    var_emp_cf_array_n(:,boot) = int_vec_cf_n;
    mean_emp_cf_array_n(:,boot) = mean_vec_cf_n;
    
    fano_emp_ff_array_n(:,boot) = int_vec_ff_n ./ mean_vec_ff_n;
    var_emp_ff_array_n(:,boot) = int_vec_ff_n;
    mean_emp_ff_array_n(:,boot) = mean_vec_ff_n;
    
    fano_emp_cf_array_s(:,boot) = int_vec_cf_s ./ mean_vec_cf_s;
    var_emp_cf_array_s(:,boot) = int_vec_cf_s;
    mean_emp_cf_array_s(:,boot) = mean_vec_cf_s;
    
    fano_emp_ff_array_s(:,boot) = int_vec_ff_s ./ mean_vec_ff_s;
    var_emp_ff_array_s(:,boot) = int_vec_ff_s;
    mean_emp_ff_array_s(:,boot) = mean_vec_ff_s;
end
%% Define Full Functions (with transients) 
%NL: This isn't done yet, but important to make sure we understand what's
%going on a the start of these traces
% var_cf_p = @(x) ((x(2)*x(4))/(x(2) + x(3)) + (2*x(3)*x(2)*x(4)^2)/(x(2)+x(3))^3)*x(1) ...
%     + (2*x(3)*(x(2) + 2*x(4)) + x(2)^2*(1-2*x(4)*x(1))...
%     + x(3)^2*(1+2*x(4)*x(1)))*exp(-(x(2)+x(3))*x(1))/(x(2) + x(3))^4 ...
%     -(x(2)^2 * x(4)^2 *exp(-2*(x(2) + x(3))*x(1)))/(x(2)+x(3))^4 ...
%     + (x(2)^2 * x(4)^2 - x(2)*(x(2)+x(3))^2*x(4) - 2*x(2)*x(3)*x(4)^2);
% 
% mean_cf = @(x) (x(2)*x(4))/(x(2) + x(3)^2)*(-1 + (x(2)+x(3))*x(1) + exp(-(x(2)+x(3))*x(1)));
%     
%% Compare Steady State Peredictions to simulations
mean_fano_cf_n_emp = mean(fano_emp_cf_array_n,2);
err_fano_cf_n_emp = std(fano_emp_cf_array_n')';

mean_fano_ff_n_emp = mean(fano_emp_ff_array_n,2);
err_fano_ff_n_emp = std(fano_emp_ff_array_n')';

k_on = simulation_set(1).k_on_true;
k_off = simulation_set(1).k_off_true;

t_res = simulation_set(1).t_res;

%Compare emprical Fano Factor to Steady State expectation
%For pure poisson initiation
% 
% var_predicted_p = zeros(1,length(time_vec));
% mean_predicted = zeros(1,length(time_vec));
% for t = 1:length(time_vec)
%     input = [time_vec(t),k_on,k_off,avg_rate];
%     var_predicted_p(t) = factor*var_cf_p(input);
%     mean_predicted(t) = factor*mean_cf(input);
% end
% fano_p = var_predicted_p ./ mean_predicted;      
%For continuous initiation
T = w*t_res;
mean_predicted = mean_rate_fun(k_on,k_off,avg_rate,T)*(1:length(time_vec));%
var_predicted_c = intrinsic_noise_fun(k_on,k_off,avg_rate,T)*(1:length(time_vec))*w;
fano_c = var_predicted_c ./ mean_predicted;

%Compare empirical Fano to Steady State Expectation
fano_fig = figure;
hold on
errorbar(time_vec,mean_fano_cf_n_emp,err_fano_cf_n_emp);
plot(time_vec,fano_c,'.-','Color',cm(30,:));
legend('Simualtion','Prediction','Location','southeast');
title('Comparing Simulated Fano Factor with Steady State Predictions')
grid on


%%
d_F = mean_predicted ./ time_vec * t_res;
d_var_c = var_predicted_c / w ./ time_vec * t_res;
d_fano_c = d_var_c ./ d_F;

d_fano_fig = figure;
hold on
errorbar(time_vec,mean_fano_ff_n_emp,err_fano_ff_n_emp);
plot(time_vec,mean_fano_ff_n_emp,'-','Color',cm(10,:));
plot(time_vec,d_fano_c,'.-','Color',cm(30,:));
% plot(time_vec(w:end),d_fano_c(w:end),'--','Color',cm(50,:));
legend('Simulation', 'SS Prediction (Continuous)','Location','northeast')
title('Comparing Simulated Fano Factor with Steady State Predictions')
grid on
clear
close all

% load master data set
load('../data/masterSet.mat')
load('../data/hmm_results_t_window50_t_inf25.mat')
%%

% extract key inference parameters     
nStates = 3;
nSteps = 7;
alpha = hmm_results.alpha;
Tres = 20;
eps = 1e-4;

% generate down-sampled vectors for inference
Tres_data = masterSet(1).tresInterp;

% generate master arrays to track sister spot activity across time and space
time_bounds = [3 40]*60;
minDP = 10;
minOverlap = 5;

masterTimeVec = 0:Tres:60*60;

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
    
    % perform ds interpolation
    t1_vec_new = masterTimeVec(find(masterTimeVec<=t1_vec(1),1,'last'):find(masterTimeVec>=t1_vec(end),1));
    f1_vec_new = interp1(t1_vec,f1_vec,t1_vec_new,'linear',true);
    t2_vec_new = masterTimeVec(find(masterTimeVec<=t2_vec(1),1,'last'):find(masterTimeVec>=t2_vec(end),1));
    f2_vec_new = interp1(t2_vec,f2_vec,t2_vec_new,'linear',true);
    
    [t_full,ia_orig] = intersect(t1_vec,t2_vec);
    [t_vec_i,ia,~] = intersect(t1_vec_new,t2_vec_new);    
    if length(t_vec_i) >= minOverlap
        ap_vec_new = interp1(t_full,ap_vec(ia_orig),t_vec_i,'linear',true);     

        % add to arrays
        fluo_array(ismember(masterTimeVec,t1_vec_new),n,1) = f1_vec_new;
        fluo_array(ismember(masterTimeVec,t2_vec_new),n,2) = f2_vec_new;
        ap_array(ismember(masterTimeVec,t_vec_i),n,1) = ap_vec_new;
        time_array(ismember(masterTimeVec,t1_vec_new),n,1) = masterTimeVec(ismember(masterTimeVec,t1_vec_new));
        time_array(ismember(masterTimeVec,t2_vec_new),n,2) = masterTimeVec(ismember(masterTimeVec,t2_vec_new));
    end
end

%%


% loop through inference groups
% wb = waitbar(0, 'Conducting single trace fits...');



% get model parameters      
A = reshape(hmm_results.A_mean,nStates,nStates);
A_log = log(A);
[V,D] = eig(A);
[~,mi] = max(diag(D));
pi0 = V(:,mi)/sum(V(:,mi));
v = hmm_results.initiation_mean;
sigma = hmm_results.noise_mean;
pi0_log = log(pi0); 

% obtain subset of valid traces                

fluo_values = {};     
index_tracker = [];
spot_tracker = [];
for i = 1:size(fluo_array,2)
    f1_vec = fluo_array(:,i,1);
    f2_vec = fluo_array(:,i,2);
    if ~all(isnan(f1_vec))
        start_i_1 = find(~isnan(f1_vec),1);
        stop_i_1 = find(~isnan(f1_vec),1,'last');
        fluo_values(end+1) = {f1_vec(start_i_1:stop_i_1)'};
        index_tracker(end+1) = i;
        spot_tracker(end+1) = 1;
        start_i_2 = find(~isnan(f2_vec),1);
        stop_i_2 = find(~isnan(f2_vec),1,'last');
        fluo_values(end+1) = {f2_vec(start_i_2:stop_i_2)'};
        index_tracker(end+1) = i;
        spot_tracker(end+1) = 2;
    end    
end    

% initialize pool   
nWorkersMax = 25;
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool(nWorkersMax);
elseif p.NumWorkers~= nWorkersMax
    delete(p);
    parpool(nWorkersMax);
end       
%                 disp('conducting viterbi trace fits...')

soft_fits = struct;

parfor f = 1:length(fluo_values)

    local_em_outputs = local_em_MS2_reduced_memory (fluo_values(f), ...
                    v', sigma, pi0_log, A_log, nStates, nSteps, alpha, 1, eps);

    soft_fits(f).state_prob_array = exp(local_em_outputs.soft_struct.p_z_log_soft{1});
    soft_fits(f).transition_prob_array = exp(local_em_outputs.soft_struct.p_zz_log_soft{1});  
    disp(num2str(100*f/length(fluo_values)));
end

%% use soft-decoded traces to generate array of instantaneous transcription rates
initiation_array = NaN(size(fluo_array));
for f = 1:length(fluo_values)                    
    fnames = fieldnames(soft_fits(f));
    ss_array = soft_fits(f).state_prob_array;
    r_vec = v'*ss_array;
    f_vec = fluo_array(:,index_tracker(f),spot_tracker(f));
    initiation_array(~isnan(f_vec),index_tracker(f),spot_tracker(f)) = r_vec;
end                

trace_fit_struct = struct;
trace_fit_struct.fluo_array = fluo_array;
trace_fit_struct.r_array = initiation_array;
trace_fit_struct.hmm_results = hmm_results;
trace_fit_struct.ap_array = ap_array;
trace_fit_struct.time_array = time_array;

% save
save('../data/trace_fit_struct.mat','trace_fit_struct')
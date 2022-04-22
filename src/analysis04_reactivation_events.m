clear
close all

% make figure directory
figPath = '../fig/spot_drift_analysis/';
mkdir(figPath)

% load master data set
load('../data/masterSet.mat')

%% attempt to use scatters to visualize spatiotemporal activity patterns
fluo_factor = 1e4;
fluo_vec = [masterSet.fluo]/fluo_factor;
fluo_vec(fluo_vec<=1) = 1;
f_max = prctile(fluo_vec,95);
fluo_vec(fluo_vec>f_max) = f_max;
ap_vec = [masterSet.APPosParticle];
time_vec = [masterSet.time]/60;
fluo_ft = ~isnan(fluo_vec)&~isnan(time_vec)&~isnan(ap_vec);

close all
act_fig = figure;
cmap = flipud(brewermap([],'Spectral'));
colormap(cmap);
scatter(ap_vec(fluo_ft),time_vec(fluo_ft),fluo_vec(fluo_ft),fluo_vec(fluo_ft),'filled','MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.05)

xlabel('% embryo length')
ylabel('minutes into nc14')
set(gca,'ydir','reverse')
set(gca,'Fontsize',14)
% xlim([0 1])
% set(gca,'Color',[228,221,209]/255) 
ylim([0 50])
xlim([20 90])
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';



%% Identify "interesting" cases
% OFF -> ON
% ON -> OFF -> ON

big_gap = ceil(15*60/16);
late_on_time = 25*60;
time_bounds = [3 40]*60;
minDP = 10;
minOverlap = 5;
ap_time = 25;
nBoots = 100;
max_on_time = 60*50;
min_dur = 3;

long_flags = false(1,length(masterSet));
on_off_on_flags = false(1,length(masterSet));
off_on_flags = false(1,length(masterSet));

on_off_time_cell = {};
on_off_ap_cell = {};
on_off_fluo_cell = {};

late_on_time_cell = {};
late_on_ap_cell = {};
late_on_fluo_cell = {};

for i = 1:length(masterSet)
    time_vec = masterSet(i).time;
    fluo_vec = masterSet(i).fluo;
    qc_flag = masterSet(i).TraceQCFlag;
    if isnan(qc_flag)
        qc_flag = false;
    end
    long_flags(i) = time_vec(1)<=time_bounds(1) & time_vec(end)>=time_bounds(2);
    if long_flags(i) && qc_flag
        % look for long runs of zeros book-ended by activity
        nn_indices = find(~isnan(fluo_vec));
        dd_vec = diff(nn_indices);
        [max_diff,mi] = max(dd_vec);
        if max_diff >= big_gap && mi>=minDP && (length(nn_indices)-mi)>=minDP
            on_off_on_flags(i) = true;
            on_off_time_cell(end+1) = {time_vec/60};
            on_off_ap_cell(end+1) = {masterSet(i).APPosNucleus};
            fluo_vec(isnan(fluo_vec)) = 0;
            on_off_fluo_cell(end+1) = {fluo_vec};
        end
        off_on_flags(i) = min(time_vec(~isnan(fluo_vec))>=late_on_time) && length(nn_indices) >= minDP;
        if off_on_flags(i)
            late_on_time_cell(end+1) = {time_vec/60};
            late_on_ap_cell(end+1) = {masterSet(i).APPosNucleus};
            fluo_vec(isnan(fluo_vec)) = 0;
            late_on_fluo_cell(end+1) = {fluo_vec};
        end
    end
end    

%% Make grid indicating fraction active vs time
ap_axis = 0:90;
ap_array = repmat(ap_axis,length(time_axis),1);
time_axis = 0:50;
time_array = repmat(time_axis',1,length(ap_axis));

ap_vec_nuc = round([masterSet.APPosNucleus]);
time_vec_rd = round([masterSet.time]/60);
fluo_vec = round([masterSet.fluo]);

frac_active_array = NaN(length(time_axis),length(ap_axis));
for a = 1:length(ap_axis)
    for t = 1:length(time_axis)
        at_filter = ap_vec_nuc==ap_axis(a) & time_vec_rd==time_axis(t);
        if sum(at_filter) >= 25
            frac_active_array(t,a) = mean(~isnan(fluo_vec(at_filter)));
        end
    end
end

%% Make plot 
close all
rng(254);
n_plot = 100;
line_alpha = 0.75;

ra_fig = figure;
cmap = brewermap([],'Greys');
cmap2 = brewermap([],'Set2');
colormap(cmap);
p = pcolor(imgaussfilt(frac_active_array,0.5));
set(p,'EdgeAlpha',0);
set(p,'FaceAlpha',0.5);
% scatter(ap_array(:),time_array(:),50,frac_active_array(:),'s','filled','MarkerEdgeAlpha',0,'MArkerFaceAlpha',1)
hold on
set(gca,'YDir','reverse')

for ind = randsample(1:length(on_off_fluo_cell),n_plot,false)     
    f_vec = on_off_fluo_cell{ind};
    ap1_vec = on_off_ap_cell{ind};
    if max(abs(diff(ap1_vec))) < 2
        ap2_vec = ap1_vec;
    %     ap1_vec(0==(f_vec)) = NaN;
        ap2_vec(0~=(f_vec)) = NaN;
        t1_vec = on_off_time_cell{ind};
        t2_vec = t1_vec;
    %     t1_vec(0==(f_vec)) = NaN;
        t2_vec(0~=(f_vec)) = NaN;
        plot(ap1_vec,t1_vec,'LineWidth',1,'Color',[cmap2(5,:) 1])
        plot(ap2_vec,t2_vec,'LineWidth',1,'Color',[cmap2(3,:) 1])
    end
end    
% plot(repelem(35,100),linspace(0,50),'k','LineWidth',1.5)
xlim([20 90])
ylim([0 50])

xlabel('% embryo length')
ylabel('minutes into nc14')
h = colorbar;
ylabel(h,'fraction of active spots')
set(gca,'Fontsize',14)
% xlim([0 1])
% set(gca,'Color',[228,221,209]/255) 
ylim([5 50])
xlim([20 90])
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
saveas(ra_fig,[figPath 'drift_reactivation.png'])
saveas(ra_fig,[figPath 'drift_reactivation.pdf'])

%%
close all


late_on_fig = figure;
cmap = brewermap([],'Greys');
cmap2 = brewermap([],'Set2');
colormap(cmap);
p = pcolor(imgaussfilt(frac_active_array,0.5));
set(p,'EdgeAlpha',0);
set(p,'FaceAlpha',0.5);
% scatter(ap_array(:),time_array(:),50,frac_active_array(:),'s','filled','MarkerEdgeAlpha',0,'MArkerFaceAlpha',1)
hold on
set(gca,'YDir','reverse')

for ind = randsample(1:length(late_on_fluo_cell),n_plot,false)    
    f_vec = late_on_fluo_cell{ind};
    ap1_vec = late_on_ap_cell{ind};
    if nanstd(ap1_vec) < 3
        ap2_vec = ap1_vec;
%         ap1_vec(0==(f_vec)) = NaN;
        ap2_vec(0~=(f_vec)) = NaN;
        t1_vec = late_on_time_cell{ind};
        t2_vec = t1_vec;
%         t1_vec(0==(f_vec)) = NaN;
        t2_vec(0~=(f_vec)) = NaN;
        plot(ap1_vec,t1_vec,'LineWidth',1,'Color',[cmap2(5,:) 1])
        plot(ap2_vec,t2_vec,'LineWidth',1,'Color',[cmap2(3,:) 1])
    end
end    
% plot(repelem(35,100),linspace(0,50),'k','LineWidth',1.5)

xlabel('% embryo length')
ylabel('minutes into nc14')
h = colorbar;
ylabel(h,'fraction of active spots')
set(gca,'Fontsize',14)
% xlim([0 1])
% set(gca,'Color',[228,221,209]/255) 
ylim([5 50])
xlim([20 90])
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
saveas(ra_fig,[figPath 'drift_late_on.png'])
saveas(ra_fig,[figPath 'drift_late_on.pdf'])
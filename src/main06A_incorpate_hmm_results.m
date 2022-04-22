% Script to conduct "soft" decoding of MS2 traces
clear
close all
warning('off','all') %Shut off Warnings
addpath(genpath('functions'))


singleTraceFitInfo = struct;

ProjectName = 'hbBAC-MS2-25C/NC13';
inferenceModel = 'w7_K3_p0_ap9_t1_f2D';

singleTraceFitInfo.projectNameCell = ProjectName;
singleTraceFitInfo.inferenceModel = inferenceModel;

% set project identifiers (only applicable if running this on savio)
 %hbBAC-MS2-20C'}; % {'2xDl-Ven_hbP2P-mCh'};


% set inference options
singleTraceFitInfo.skipViterbiFitFlag = 0;
singleTraceFitInfo.skipSSFitFlag = 0;
singleTraceFitInfo.makeLongFormSet = 0;
singleTraceFitInfo.bootstrap_flag = 0;
singleTraceFitInfo.n_bootstraps = 1;

% Get basic project info and determing file paths
if contains(ProjectName(end), '/') | contains(ProjectName(end), '\') 
    ProjectName = ProjectName(1:end-1);
end

if ~contains(ProjectName, '/') & ~contains(ProjectName, '\') 
    liveProject = LiveEnrichmentProject(singleTraceFitInfo.projectNameCell);
else
    liveProject = LiveEnrichmentProject(fileparts(singleTraceFitInfo.projectNameCell));
end
singleTraceFitInfo.resultsDir = '';
singleTraceFitInfo.resultsRoot = ''; 
% save
slashes = regexp(liveProject.dataPath,'/|\');
dataDir = liveProject.dataPath(1:slashes(end-1));
singleTraceFitDir = [dataDir 'singleTraceFitsDirectory' filesep];
mkdir(singleTraceFitDir)
singleTraceFitDir2 = [singleTraceFitDir, ProjectName, filesep ];
mkdir(singleTraceFitDir2)
if ~isempty(singleTraceFitInfo.inferenceModel)
    singleTraceFitDir2 = [singleTraceFitDir2 singleTraceFitInfo.inferenceModel filesep];
end

mkdir(singleTraceFitDir2)
save([singleTraceFitDir2 'singleTraceFitInfo.mat'],'singleTraceFitInfo')


% copy bash file to inference directory
copyfile('GM_run_singleTraceFits.sh',[singleTraceFitDir2 'run_singleTraceFits.sh'])
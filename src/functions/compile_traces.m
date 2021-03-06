function masterSet = compile_traces(varargin)

% OUTPUT
% spot_struct: compiled data set contain key nucleus and
%                 particle attributes, contains the following fields:
%
%     setID: ID number for this dataset (experiment), see set_key.mat to
%            match setID to experiment Prefix
%     sourcePath: full path to the experiment's DynamicsResults folder
%     APFlag: indicates whether or not there is AP position infor for this
%             experiment
%     threeDFlag: indicates whether or not there is 3D spot fitting data
%                 for this experiment
%     xDim: # of pixels in X
%     yDim: # of pixels in Y
%     zDim: # of steps in Z
%     zStep: size of steps in z, in um
%     pixelSize: size of xy pixels, in um
%     particleID:
%     xPosParticle:
%     yPosParticle:
%     zPosParticle:
%     APPosParticle:
%     spotFrames:
%     xPosParticle3D:
%     yPosParticle3D:
%     zPosParticle3D:
%     fluo3D:
%     fluo:
%     fluoOffset:
%     qcFlag:
%     N: ???
%     sparsity:
%     xPosNucleus: x position of center of nucleus
%     yPosNucleus: y position of center of nucleus
%     rawNCProtein:
%     frames:
%     nucleusID
%     ncID
%     ncStart
%     minDP
%     minTime
%     pctSparsity
%
%     time
%
%     sisterIndex
%     sisterParticleID

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Set defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all force
addpath(genpath('utilities'));
% Defaults
firstNC = 14;   % first nuclear cycle to pull data from
minDP = 14;     % particles with fewer than minDP points will be flagged
pctSparsity = 50;   %
twoSpotFlag = true;
minTime = 0*60; % take no fluorescence data prior to this point
tresInterpFloor = 15;
SpotChannelIndex = 1; % this does nothing at the moment, but can serve as a starting point if ever we wish to analyze two-spot-two-color data

%% %%%%%%%%%%%%%%%%%%%%%%% Process input parameters %%%%%%%%%%%%%%%%%%%%%%%

for i = 1:numel(varargin)
    if ischar(varargin{i})
        if strcmpi(varargin{i}, 'firstNC')
            if  i < numel(varargin) && isnumeric(varargin{i+1}) && (9 <= varargin{i+1} && varargin{i+1} <= 14)
                firstNC = varargin{i+1};
            else
                error(['The ''firstNC'' option must be followed by an integer from 9 to 14, inclusive. Default firstNC is ' num2str(firstNC) '.'])
            end
        elseif strcmpi(varargin{i}, 'NC')
            if  i < numel(varargin) && isnumeric(varargin{i+1}) && (9 <= varargin{i+1} && varargin{i+1} <= 14)
                lastNC = varargin{i+1};
                firstNC = varargin{i+1};
            else
                error(['The ''NC'' option must be followed by an integer from 9 to 14, inclusive. Default NC is ' num2str(firstNC) '.'])
            end
        elseif strcmpi(varargin{i}, 'projectName')
            error('Pipeline currently does not support alternate names for projects (otehr than what is on the data status tab)')
        elseif i < numel(varargin)
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end

%% %%%%%%%%%%%%%%%%% Extract relevant processed data %%%%%%%%%%%%%%%%%%%%%%

% get list of files 
exp_files = dir('../data/embryos/*eve2sp*');
numExperiments = length(exp_files);

% Generate master structure with info on all nuclei and traces in
% constituent sets
masterSet = [];

% Add data from each experiment to the master structure
h = waitbar(0,'Compiling data ...');

for i = 1:numExperiments
    
    waitbar((i-1)/numExperiments,h, ['Compiling data for dataset ' num2str(i) ' of ' num2str(numExperiments)])
    
    setID = i;    
    Prefix = exp_files(i).name;
    %%%%%%%% Read in processed data from main mRNADyanmics pipeline %%%%%%%   
    load([exp_files(i).folder filesep Prefix filesep Prefix '_lin.mat'],'schnitzcells')    
    
    processedData = load([exp_files(i).folder filesep Prefix filesep 'CompiledParticles.mat']);   % this structure contains all info compiled for this experiment
    compiledParticles = processedData.CompiledParticles;    % this is the cell array containing the CompiledParticles cell array
    if iscell(compiledParticles)
        compiledParticles = compiledParticles{SpotChannelIndex};               % this assumes there is only one channel with spots to analyze
    end
        
    % Pull trace, time, and frame variables
    timeRaw = processedData.ElapsedTime*60; %time vector, [sec]
    framesRaw = 1:length(timeRaw); % Frame list
    tracesRaw = processedData.AllTracesVector;  %array with a column for each trace
    if iscell(tracesRaw)        %check to see if traces are stored in cell array
        tracesRaw = tracesRaw{SpotChannelIndex};
    end
    
    % Load FrameInfo array
    load([exp_files(i).folder filesep Prefix filesep 'FrameInfo.mat'],'FrameInfo')
    nc14_frame = processedData.nc14;    
    firstTime = timeRaw(nc14_frame);    

    % initialize structure to store nucleus and particle info
    compiledSchnitzCells = struct;

    nucleusCounter = 1;    %counter to ensure no empty spaces in the compiledSchnitzCells
    %struct if some schnitzcells have no frames in the
    %desired nuclear cycle

    firstNCFrame = nc14_frame;
    lastNCFrame = length(FrameInfo);
        
    % filter reference vectors
    framesNC = framesRaw(firstNCFrame:lastNCFrame);
    tracesNC = tracesRaw(firstNCFrame:lastNCFrame,:);
    timeNC = timeRaw(firstNCFrame:lastNCFrame);
    timeNC = timeNC - firstTime; %normalize to start of first nc
        
    %%%%%%%%%%%%%%%%%%%% Compile nucleus schnitz info %%%%%%%%%%%%%%%%%%%%%
        
    for s = 1:length(schnitzcells)
        schnitzFrames = schnitzcells(s).frames;

        % Filter for only the frames in this nc
        ncFilter = ismember(schnitzFrames,framesNC);
        rawNucleusFrames = schnitzFrames(ncFilter);

        if length(rawNucleusFrames) >= 1     %only grab data from the desired nuclear cycle(s)

            % Add info that's the same for all schnitzcells in this
            % experiment
            compiledSchnitzCells(nucleusCounter).setID = setID;
            compiledSchnitzCells(nucleusCounter).nc = 14;

            % Initialize particle fields--will be set to particle values
            % for nuclei with matching particle
            compiledSchnitzCells = initializeParticleFields(...
                compiledSchnitzCells,nucleusCounter,sum(ncFilter),0,1);

            % Add QC-related flags
            compiledSchnitzCells(nucleusCounter).TraceQCFlag = NaN;
            compiledSchnitzCells(nucleusCounter).FrameQCFlags = NaN(1,length(ncFilter));
            compiledSchnitzCells(nucleusCounter).N = NaN;
            compiledSchnitzCells(nucleusCounter).sparsity = NaN;

            % Add core nucleus info
            x = double(schnitzcells(s).cenx);
            y = double(schnitzcells(s).ceny);
            compiledSchnitzCells(nucleusCounter).xPosNucleus = x(ncFilter);
            compiledSchnitzCells(nucleusCounter).yPosNucleus = y(ncFilter);    
            compiledSchnitzCells(nucleusCounter).APPosNucleus =  NaN(1,sum(ncFilter));
            

            % Add protein and nucleus info           
            compiledSchnitzCells(nucleusCounter).frames = rawNucleusFrames';
            compiledSchnitzCells(nucleusCounter).nucleusID = s;
            compiledSchnitzCells(nucleusCounter).ncID = eval([num2str(setID) '.' sprintf('%04d',s)]);
            compiledSchnitzCells(nucleusCounter).minDP = minDP;
            compiledSchnitzCells(nucleusCounter).minTime = minTime;
            compiledSchnitzCells(nucleusCounter).pctSparsity = pctSparsity;

            % Add time and set info
            compiledSchnitzCells(nucleusCounter).time = timeNC(ismember(framesNC,rawNucleusFrames));
            if numel(unique(compiledSchnitzCells(nucleusCounter).time))~=numel(compiledSchnitzCells(nucleusCounter).time)
                error(['Non-unique time values. Check FrameInfo for Prefix: ' currExperiment.Prefix])
            end
            
            % initialize nucleus qc flag
            compiledSchnitzCells(nucleusCounter).missingNucleusFrames = 0;

            % Add sister spot fields
            compiledSchnitzCells(nucleusCounter).sisterIndex = NaN;
            compiledSchnitzCells(nucleusCounter).sisterParticleID = NaN;

            % increment counter to ensure no empty spaces in the
            % compiledSchnitzCells struct if some schnitzcells have no
            % frames in the desired nuclear cycle
            nucleusCounter = nucleusCounter + 1;
        end
    end
        
    % Index vector to cross-ref w/ particles
    schnitzIndex = [compiledSchnitzCells.nucleusID];

    % Iterate through traces
    for j = 1:size(tracesNC,2)
        % Raw fluo trace
        rawTrace = tracesNC(:,j);

        % Get nucleus ID
        schnitz = compiledParticles(j).schnitz;

        % Find first and last expression frames, requiring that spots
        % cannot apear earlier than some fixed time (minTime)
        traceStart = find(timeNC>=minTime&~isnan(rawTrace'),1);
        traceStop = find(~isnan(rawTrace),1,'last');

        % skip particles not in specified nc range
        if isempty(traceStart)
            continue
        end

        %Create versions with all intervening frames present (missing frames
        %appear as NaNs)
        traceFull = rawTrace(traceStart:traceStop)';
        traceFramesFull = framesNC(traceStart:traceStop);

        % Perform qc tests
        nDP = sum(~isnan(traceFull));
        sparsity = prctile(diff(find(~isnan(traceFull))),pctSparsity);
        TraceQCFlag = nDP >= minDP && sparsity == 1;

        % trace-nucleus mapping may be many-to-1
        ncIndex = find(schnitzIndex==schnitz);
        if length(ncIndex) ~= 1
            continue
            warning(['Problem with particle-nucleus crossreference for Prefix: ' Prefix])
        end

        % Identifier variable
        particle = compiledParticles(j).OriginalParticle;
        particleID = eval([num2str(setID) '.' sprintf('%04d',particle)]);

        % check to see if a different particle has already been assigned
        % if so, create a new entry
        if ~isnan(compiledSchnitzCells(ncIndex).particleID) && twoSpotFlag
            compiledSchnitzCells(end+1) = compiledSchnitzCells(ncIndex);

            % update cross-reference variables
            compiledSchnitzCells(ncIndex).sisterParticleID = particleID;
            compiledSchnitzCells(end).sisterParticleID = compiledSchnitzCells(ncIndex).particleID;
            compiledSchnitzCells(ncIndex).sisterIndex = length(compiledSchnitzCells);
            compiledSchnitzCells(end).sisterIndex = ncIndex;

            % redefine index variables
            ncIndex = length(compiledSchnitzCells);

            % reset particle fields
            nEntries = length(compiledSchnitzCells(ncIndex).xPosParticle);
            compiledSchnitzCells(ncIndex).ParticleID = NaN;

            % Initialize particle fields--will be set to particle values
            % for nuclei with matching particle
            compiledSchnitzCells = initializeParticleFields(...
                                compiledSchnitzCells,ncIndex,nEntries,0,1);
        end

        % Record particle identifier
        compiledSchnitzCells(ncIndex).particleID = particleID;

        % Find overlap between nucleus and trace
        rawNucleusFrames = compiledSchnitzCells(ncIndex).frames;
        spotFilter = ismember(rawNucleusFrames,traceFramesFull);

        compiledSchnitzCells(ncIndex).spotFrames = spotFilter;
        if sum(spotFilter) < length(traceFramesFull)
            warning(['Inconsistent particle and nucleus frames for Prefix: ' currExperiment.Prefix ' .Flagging...'])
            compiledSchnitzCells(ncIndex).missingNucleusFrames = 1;
        end

        % record fluorescence info
        compiledSchnitzCells(ncIndex).fluo(spotFilter) = traceFull(ismember(traceFramesFull,rawNucleusFrames));

        % Find intersection btw full frame range and CP frames
        rawParticleFrames = compiledParticles(j).Frame;        

        % make filters
        ncSpotFilter1 = ismember(rawNucleusFrames,rawParticleFrames);
        ncSpotFilter2 = ismember(rawParticleFrames,rawNucleusFrames);

        % add offset info
        compiledSchnitzCells(ncIndex).fluoOffset(ncSpotFilter1) = compiledParticles(j).Off(ncSpotFilter2);
        
        % x, y, and z info
        compiledSchnitzCells(ncIndex).xPosParticle(ncSpotFilter1) = compiledParticles(j).xPos(ncSpotFilter2);
        compiledSchnitzCells(ncIndex).yPosParticle(ncSpotFilter1) = compiledParticles(j).yPos(ncSpotFilter2);
%         compiledSchnitzCells(ncIndex).zPosParticle(ncSpotFilter1) = compiledParticles(j).zPos(ncSpotFilter2);        
        compiledSchnitzCells(ncIndex).APPosParticle(ncSpotFilter1) = compiledParticles(j).APposParticle(ncSpotFilter2)*100;
            
        % add qc info
        compiledSchnitzCells(ncIndex).N = nDP;
        compiledSchnitzCells(ncIndex).sparsity = sparsity;
        compiledSchnitzCells(ncIndex).TraceQCFlag = TraceQCFlag;
        compiledSchnitzCells(ncIndex).FrameQCFlags(ncSpotFilter1) = TraceQCFlag;
    end        
    masterSet = [masterSet  compiledSchnitzCells];
             
    waitbar(i/numExperiments,h, ['Compiling data for dataset ' num2str(i) ' of ' num2str(numExperiments)])
end
close(h)

% add nucleus AP info
disp('adding AP info for nuclei...')
set_vec = [masterSet.setID];
for n = 1:numExperiments
    set_filter = set_vec == n;
    resultsFolder = [exp_files(n).folder filesep exp_files(n).name filesep];    
    masterSet(set_filter) = convertToFractionalEmbryoLength(resultsFolder,masterSet(set_filter));    
end    

disp('interpolating data...')

% generate interpolation fields
interpFields = {'fluo','xPosParticle','yPosParticle','APPosNucleus'};

% calculate interpolation time
% sample true experimental res times
dt_vec_full = [];
multi_step_traces = find([masterSet.N]>2);
for i = randsample(multi_step_traces,min([50,length(multi_step_traces)]),false)
    t_vec = masterSet(i).time;
    dt_vec_full = [dt_vec_full diff(t_vec)];
end

% calculate rounded dt
dt_raw = median(dt_vec_full);
dt_round = floor(dt_raw/5)*5;

tresInterp = max([tresInterpFloor,dt_round]);
interpGrid = 0:tresInterp:60*60;

% use 97.5th percentile of fluorescence values to set scale for flagging
% unusual jump events
fluo_scale = prctile([masterSet.fluo],97.5);
big_blip_thresh = fluo_scale;

for i = 1:length(masterSet)
    
    % First identify isolated "islands" at start or end of trace. These
    % tend to be artifactual
    fluoVec = masterSet(i).fluo;
    frameIndex = 1:length(masterSet(i).frames);
    %         spotFrames = spot_struct(i).spotFrames;
    
    % NL: these sizes generally work ok, but may need to change if time res
    % is >> or  << ~20 seconds or reporates with elongation times >> or <<
    % 2 min
    startIndex = [];
    stopIndex = [];
    
    vecNans = ~isnan(fluoVec);
    startIndexRaw = find(vecNans,1);
    stopIndexRaw = find(vecNans,1,'last');
    nFramesRaw = stopIndexRaw-startIndexRaw+1;
    
    % this is designed to get rid of blips at start and end of traces
    if length(fluoVec) > 1
        bSize = 8;
        kernelBig = ones(1,bSize);
        sSize = 3;
        kernelSmall = ones(1,sSize);
        
        vecConvBig = conv(kernelBig,vecNans);
        vecConvSmall = conv(kernelSmall,vecNans);
        
        % generate starting flags
        vecConvStartBig = vecConvBig(bSize:end);
        vecConvStartSmall = vecConvSmall(sSize:end);
        
        startIndex = find(vecConvStartSmall>=sSize | vecConvStartBig>=floor(bSize/2.1) & vecNans,1);
        
        % generate ending flags
        vecConvEndBig = vecConvBig(1:end-bSize+1);
        vecConvEndSmall = vecConvSmall(1:end-sSize+1);
        
        stopIndex = find(vecConvEndSmall>=sSize | vecConvEndBig>=floor(bSize/2.1) & vecNans,1,'last');
    end

    % update QC flags
    if ~isempty(startIndex) && ~isempty(stopIndex)
        masterSet(i).FrameQCFlags(frameIndex<startIndex & vecNans) = false;
        masterSet(i).FrameQCFlags(frameIndex>stopIndex & vecNans) = false;
    else
        masterSet(i).FrameQCFlags(vecNans) = false;
    end
    
    if ~isnan(masterSet(i).TraceQCFlag)
        masterSet(i).TraceQCFlag =  masterSet(i).TraceQCFlag && nansum(masterSet(i).FrameQCFlags)>=minDP;
    end
    
    masterSet(i).truncatedFlag = 0;
    if ~isempty(startIndex) && ~isempty(stopIndex)
        nFramesNew = stopIndex-startIndex+1;
        masterSet(i).nFramesTruncated = nFramesRaw-nFramesNew;
        masterSet(i).truncatedFlag = masterSet(i).nFramesTruncated>0;
    elseif ~isempty(startIndexRaw) && ~isempty(stopIndexRaw)
        masterSet(i).nFramesTruncated = nFramesRaw;
        masterSet(i).truncatedFlag = 1;
    end
    
    % generate time reference vefctors for interpolation
    timeVec = masterSet(i).time;
    timeVec = timeVec(startIndex:stopIndex);
    
    if length(timeVec)>1
        
        % add caps before and after
        first_ind = find(interpGrid>=timeVec(1),1);
        first_ind = max([1 first_ind-1]);
        last_ind = find(interpGrid<=timeVec(end),1,'last');
        last_ind = min([length(interpGrid) last_ind+1]);
        timeInterp = interpGrid(first_ind:last_ind);
                
        if isempty(timeInterp)
            timeInterp = NaN;
            for  j = 1:numel(interpFields)
                masterSet(i).([interpFields{j} 'Interp']) = NaN;
            end
        else
            for  j = 1:length(interpFields)
                vec = masterSet(i).(interpFields{j})(startIndex:stopIndex);
                
                % perform basic QC on fluorescence fields
                if contains(interpFields{j},'fluo')
                    
                    %Look for clusters of 6 or more NaNs
                    kernel = [1,1,1,1,1];
                    vecNans = isnan(vec);
                    vecConv = conv(kernel,vecNans);
                    vecConv = vecConv(3:end-2);
                    zIDs = find(vecConv==numel(kernel));
                    zIDs = unique([zIDs-1 zIDs zIDs+1]); % get set of z_ids
                    vec(zIDs) = 0; % set clusters to zeros
                    vec(vec<0) = 0; % deal with negative values
                    
                    % find single dp "blips". These will be replaced via interpolation
                    % interpolate remaining NaNs
                    fluo_dd = [0 diff(vec,2) 0];
                    vec(abs(fluo_dd)>big_blip_thresh) = NaN;
                    
                    % find and remove suspsiciously large rises
                    fluo_d = [0 diff(vec,1) 0];
                    vec(abs(fluo_d)>big_blip_thresh*0.75) = NaN;
                    
                end
                
                % interpolate remaining NaNs
                referenceTime = timeVec(~isnan(vec));
                referenceVec = vec(~isnan(vec));
                
                if length(referenceTime)>1
                    % Interpolate to standardize spacing
                    vec_i = interp1(referenceTime,referenceVec,timeInterp,'linear','extrap');
                    vec_i(vec_i<0) = 0;
                    masterSet(i).([interpFields{j} 'Interp']) = vec_i;
                elseif length(referenceTime)==1 && referenceTime >= timeInterp(1) && referenceTime <= timeInterp(end)
                    vecTo = NaN(size(timeInterp));
                    [~,mi] = min(abs(referenceTime-timeInterp));
                    vecTo(mi) = referenceVec;
                    masterSet(i).([interpFields{j} 'Interp']) = vecTo;
                else
                    masterSet(i).([interpFields{j} 'Interp']) = NaN(size(timeInterp));
                end
            end
        end
    else
        timeInterp = NaN;
        for  j = 1:numel(interpFields)
            masterSet(i).([interpFields{j} 'Interp']) = NaN;
        end
    end  
    masterSet(i).timeInterp = timeInterp;
    masterSet(i).tresInterp = tresInterp;    
end

% save
save('../data/masterSet.mat' ,'masterSet')

disp('done.')

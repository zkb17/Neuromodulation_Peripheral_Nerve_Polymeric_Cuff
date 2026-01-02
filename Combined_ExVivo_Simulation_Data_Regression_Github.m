% Expansion of exvivo_data_analysis_post_trial_9_triggered_data to automate
% the generation and saving of figures across multiple data sets

close all
clear all

% Capture the path of the scripts file that the main matlab script 
% saved, so you can access the other saved scripts or matlab arrays easily
% as long as they are saved in the same relative positions across any
% computer used
scriptDir = fileparts(mfilename('fullpath'));  % location of the current .m filedataFile

% Navigate to the folder containing the data from all trials. Should note that every file in this folder
% will be run through the system. I am setting a default folder as this
% adds time to navigate to the folder everytime to debug. It allows to
% change from the default if desired.

% Default file path
allTrialsDir = "C:\Users\Green_Group\Zack\ExVivo\ExVivoData\All_Trials_Reduced_Data";

% Display the current default path and ask the user if they want to change it
fprintf('Current default directory for full data:\n%s\n', allTrialsDir);
changePath = input('Do you want to select a new directory? (y/n): ', 's');

% If user wants to change the path
if strcmpi(changePath, 'y')
    allTrialsDir = uigetdir;
    
    % If user cancels the selection, keep the default
    if allTrialsDir ~= 0
        dataPath = allTrialsDir;
    else
        fprintf('No new directory selected. Using default.\n');
        dataPath = allTrialsDir;
    end
else
    dataPath = allTrialsDir;
end

% Show the final path being used
fprintf('Using directory:\n%s\n', dataPath);

% Load the names of each trial data set
cd (allTrialsDir)
allTrialsDirContent = dir(allTrialsDir);
allTrialsDirContent = allTrialsDirContent(~ismember({allTrialsDirContent.name}, {'.', '..'}));
% sort trials in the dir contents
% Initialize variables
trialNumbers = []; % Array to store numerical trial numbers
validTrials = []; % Indices of valid trials in allTrialsDirContent

% Extract trial numbers and filter for valid "Trial X" patterns
for i = 1:numel(allTrialsDirContent)
    % Extract current name
    currentName = allTrialsDirContent(i).name;
    % Check if it starts with "Trial" and extract the trial number
    if startsWith(currentName, 'Trial')
        trialNumber = sscanf(currentName, 'Trial_%d');
        if ~isempty(trialNumber) % Ensure a valid trial number
            trialNumbers(end+1) = trialNumber; % Store trial number
            validTrials(end+1) = i; % Store index of the valid trial
        end
    end
end

% Sort trial numbers and reorder allTrialsDirContent accordingly
[sortedTrialNumbers, sortIndices] = sort(trialNumbers);

% Filter for trials in the range 8 to 12
validRangeIndices = sortedTrialNumbers >= 8 & sortedTrialNumbers <= 12;
finalIndices = validTrials(sortIndices(validRangeIndices));

% Update allTrialsDirContent to reflect the desired order
allTrialsDirContent = allTrialsDirContent(finalIndices);

trialNum = 0;
% For each element in the Trial folder
for i = 1:numel(allTrialsDirContent)
    trialNum = trialNum + 1;
    % delimiter string used in naming convention was underscore
    trialNames{trialNum} = split(allTrialsDirContent(i).name,'_');
    % remove file extension from end of split string
    trialNames{trialNum}{end} = split(trialNames{trialNum}{end},'.');
    trialNames{trialNum}{end} = trialNames{trialNum}{end}{1};
    % only thing needed for trial names is the second cell
    trialNames{trialNum} = [trialNames{trialNum}{1},'_',trialNames{trialNum}{2}];
end

% Navigate to the folder containing the stimulator config files matched to
% the ex vivo data in the selected folder

% Default file path
configDir = "C:\Users\Green_Group\Zack\ExVivo\StimConfigFiles\Experiments 9-~ - CURRENT USE";

% Display the current default path and ask the user if they want to change it
fprintf('Current default directory for config trials:\n%s\n', configDir);
changePath = input('Do you want to select a new directory? (y/n): ', 's');

% If user wants to change the path
if strcmpi(changePath, 'y')
    configDir = uigetdir;
    
    % If user cancels the selection, keep the default
    if configDir ~= 0
        dataPath = configDir;
    else
        fprintf('No new directory selected. Using default.\n');
        dataPath = configDir;
    end
else
    dataPath = configDir;
end

% Show the final path being used
fprintf('Using directory:\n%s\n', dataPath);

% Repeat for the config folder
cd(configDir)
configDirContents = dir(configDir);
% Remove any folders
configNotDir = ~[configDirContents.isdir];
configDirContents = configDirContents(configNotDir);
configNum = 0;
% For each element in the folder
for i = 1:numel(configDirContents)
    configNum = configNum+1;
    % Delimiter string we used in the naming convention was an
    % underscore
    configNames{configNum} = split(configDirContents(i).name,'_');
    % The only thing that's needed is the config number, given by the
    % second cell
    configNames{configNum} = configNames{configNum}{2};
end

% Navigate to the folder containing the stimulator config files for Trial 8

% Default file path
configDir8 = "C:\Users\Green_Group\Zack\ExVivo\StimConfigFiles\Trial_8";

% Display the current default path and ask the user if they want to change it
fprintf('Current default directory for config Trial 8:\n%s\n', configDir8);
changePath = input('Do you want to select a new directory? (y/n): ', 's');

% If user wants to change the path
if strcmpi(changePath, 'y')
    configDir8 = uigetdir;
    
    % If user cancels the selection, keep the default
    if configDir8 ~= 0
        dataPath = configDir8;
    else
        fprintf('No new directory selected. Using default.\n');
        dataPath = configDir8;
    end
else
    dataPath = configDir8;
end

% Repeat for the config folder
cd(configDir8)
configDir8Contents = dir(configDir8);
% Remove any folders
configNotDir8 = ~[configDir8Contents.isdir];
configDir8Contents = configDir8Contents(configNotDir8);
configNum8 = 0;
% For each element in the folder
for i = 1:numel(configDir8Contents)
    configNum8 = configNum8+1;
    % Delimiter string we used in the naming convention was an
    % underscore
    configNames8{configNum8} = split(configDir8Contents(i).name,'_');
    % The only thing that's needed is the config number, given by the
    % second cell
    configNames8{configNum8} = configNames8{configNum8}{2};
end

% User selects a folder to save figures into

% Default file path
saveDir = "C:\Users\Green_Group\Zack\ExVivo\Regression_Results";

% Display the current default path and ask the user if they want to change it
fprintf('Current default directory for saving data:\n%s\n', saveDir);
changePath = input('Do you want to select a new directory? (y/n): ', 's');

% If user wants to change the path
if strcmpi(changePath, 'y')
    saveDir = uigetdir;
    
    % If user cancels the selection, keep the default
    if saveDir ~= 0
        dataPath = saveDir;
    else
        fprintf('No new directory selected. Using default.\n');
        dataPath = saveDir;
    end
else
    dataPath = saveDir;
end

% Show the final path being used
fprintf('Using directory:\n%s\n', dataPath);

%% Expansion of ASCENT_Trend_Visualiser_ExVivo_Matched.m to work on a set of
% multiple sims to automate the figure generation process

% User inputs where the sample directory is that contains the different
% nerve trials

% Default file path
allSamplesDir = "C:\Users\Green_Group\Zack\ASCENT_NonGit\ascent\samples";

% Display the current default path and ask the user if they want to change it
fprintf('Current default directory for the simulation samples:\n%s\n', allSamplesDir);
changePath = input('Do you want to select a new directory? (y/n): ', 's');

% If user wants to change the path
if strcmpi(changePath, 'y')
    allSamplesDir = uigetdir;
    
    % If user cancels the selection, keep the default
    if allSamplesDir ~= 0
        dataPath = allSamplesDir;
    else
        fprintf('No new directory selected. Using default.\n');
        dataPath = allSamplesDir;
    end
else
    dataPath = allSamplesDir;
end

% Show the final path being used
fprintf('Using directory:\n%s\n', dataPath);

all_Nerves = [8 9 10 11 12];
Nerves_In_Study = [8 9 10 11 12]; % default all: [8 9 10 11 12] 
all_Samples = [21 22 23 24 25];
Samples_In_Study = [21 22 23 24 25]; % default all: [21 22 23 24 25] corresponds to trials 8, 9, 10, 11, 12
Models_In_Study = [1]; % corresponds to CE cuff (as opposed to Pt cuff)
all_Sims = [300 303 305 306 307];
Sims_In_Study = [300 303 305 306 307]; % default all: [300 303 305 306 307], can change depending on stimulation parameter of interest

% User selects the folder to save the output figures to

% Default file path
saveFolder = "C:\Users\Green_Group\Zack\ASCENT_NonGit\OutputGraphs\Combined_Regression_Outputs";

% Display the current default path and ask the user if they want to change it
fprintf('Current default directory for the simulation samples:\n%s\n', saveDir);
changePath = input('Do you want to select a new directory? (y/n): ', 's');

% If user wants to change the path
if strcmpi(changePath, 'y')
    saveFolder = uigetdir;
    
    % If user cancels the selection, keep the default
    if saveFolder ~= 0
        dataPath = saveFolder;
    else
        fprintf('No new directory selected. Using default.\n');
        dataPath = saveFolder;
    end
else
    dataPath = saveFolder;
end

% Show the final path being used
fprintf('Using directory:\n%s\n', dataPath);

%% Create cells that save data across every nerve and dataset
maxAmpsForNormUserDet_allTrials = cell(trialNum,5);
% Initialize a matrix to store the skipped row count for each nerve and
% dataset, each skipped row represents invalid data based off the peak
% finding mechanism and the stimulation artefact
skippedRowsCount = zeros(5, 5); % 5 nerves × 5 datasets
percentSkippedPerConfig = zeros(5,5);
fascicleNames = {'Tibial','Peroneal','Sural'};
figNum = 1;

regressionData = table( ...
    categorical([]), categorical([]), categorical([]), categorical([]), ... % Categorical variables describing basics
    [], [], [], [], ... % Numeric variables (Amplitude, Pulse Width)
    [], ... % Charge injected in each phase
    categorical([]), categorical([]), ... % categorical active channel numbers
    [], [], ... % Numeric distances
    categorical([]), categorical([]), ... % More categorical variables (Number of Poles, Phase Symmetry)
    [], [], ... % Normalized Fascicle Activation for ex vivo and simulation
    [], [], ... % Fascicle Selectivity metrics for ex vivo and simulation 
    'VariableNames', {'AnimalID', 'Fascicle', 'nsimID', 'configurationID', ...
    'Amplitude_Phase1', 'Amplitude_Phase2', 'PulseWidth_Phase1', 'PulseWidth_Phase2', ...
    'Charge_Injection', ...
    'Source_Electrode', 'Sink_Electrode', ...
    'Distance_Source_Fascicle', ...
    'Distance_Sink_Fascicle', ...
    'Number_Poles', 'Phase_Symmetry', ...
    'Norm_Fascicle_Activation_ExVivo', 'Norm_Fascicle_Activation_Simulation', ...
    'Fascicle_Selectivity_ExVivo','Fascicle_Selectivity_Simulation'} );

%% Load Presaved Data if it exists, otherwise set exist variable to false

allAveragedDataFile  = fullfile(scriptDir, 'all_averagedChannelData.mat');
% Check if the file exists
if exist(allAveragedDataFile) == 2
    % File exists, load it
    load(allAveragedDataFile,'all_averagedChannelData');
    allAveragedDataFileExists = true;
else
    % File does not exist, set flag to false
    allAveragedDataFileExists = false;
end

pulseParamFile  = fullfile(scriptDir, 'all_pulseParamTrends.mat');
% Check if the file exists
if exist(pulseParamFile) == 2
    % File exists, load it
    load(pulseParamFile,'all_pulseParamTrends');
    pulseParamFileExists = true;
else
    % File does not exist, set flag to false
    pulseParamFileExists = false;
end

maxAmpsUserDetFile  = fullfile(scriptDir, 'maxAmpsForNormUserDet_allTrials.mat');
% Check if the file exists
if exist(maxAmpsUserDetFile) == 2
    % File exists, load it
    load(maxAmpsUserDetFile,'maxAmpsForNormUserDet_allTrials');
    maxAmpsForNormFileExists = true;
else
    % File does not exist, set flag to false
    maxAmpsForNormFileExists = false;
end

% These two are not calcuated within this script, so they will always be
% loaded

electrodeFascicleRelsFile = fullfile(scriptDir, 'allTrialsElectrodeFascicleRels.mat');
load(electrodeFascicleRelsFile,'allTrialsElectrodeFascicleRels');

electrodeTripolarFascicleRelsFile = fullfile(scriptDir, 'tripolarAverages.mat');
allTrialsElectrodeTripolarFascicleRels = load(electrodeTripolarFascicleRelsFile);
allTrialsElectrodeTripolarFascicleRels = allTrialsElectrodeTripolarFascicleRels.allTrialsElectrodeTripolarFascicleRels;

%% Loading seaborn deep colors for plotting uses
% Seaborn "deep" palette approximate RGB values (0-indexed):
deep_blue   = [0.2980, 0.4470, 0.6900]; % Index 0
deep_green  = [0.3333, 0.6588, 0.4078]; % Index 1
deep_red    = [0.7765, 0.2941, 0.3216]; % Index 2
deep_purple = [0.5059, 0.4470, 0.6980]; % Index 3
deep_brown  = [0.8000, 0.7255, 0.4549]; % Index 4
deep_teal   = [0.3922, 0.7098, 0.8039]; % Index 5

%% Initiating the loop for each nerve trial

for nerveIndex = 1:numel(allTrialsDirContent)

    if ~any(ismember(Nerves_In_Study, all_Nerves(nerveIndex)))
        continue;
    else
    % Navigate to the folder containing the sets of ex vivo data. These should
    % all be from the same trial. Should note that every file in this folder
    % will be run through the system
    dataDir = [allTrialsDirContent(1).folder,'\',allTrialsDirContent(nerveIndex).name];
    
    % Load the names of the datasets to allow for cross-matching between each
    % dataset and the relevant config file/save location
    cd(dataDir)
    dataDirContents = dir(dataDir);
    % Remove any folders
    dataNotDir = ~[dataDirContents.isdir];
    dataDirContents = dataDirContents(dataNotDir);

    datasetNum = 0;

    % For each element in the folder
    for i = 1:numel(dataDirContents)
        datasetNum = datasetNum+1;
        % Delimiter string we used in the naming convention was an
        % underscore
        datasetNames{datasetNum} = split(dataDirContents(i).name,'_');
        % Remove the file extension from the end of the split string
        datasetNames{datasetNum}{end} = split(datasetNames{datasetNum}{end},'.');
        datasetNames{datasetNum}{end} = datasetNames{datasetNum}{end}{1};
    end
    
    % Create a folder to save the generated figures into, if it doesn't already
    % exist
    cd(saveDir)
    if ~isfolder(trialNames{nerveIndex})
        mkdir(trialNames{nerveIndex})
    end

    %% Ensure the simulation sample folder has been set for this nerve index loop
    sampleDir = fullfile(allSamplesDir, num2str(all_Samples(nerveIndex)));
    cd(sampleDir)
    sampleText = fileread('sample.json');
    sampleData = jsondecode(sampleText);
    
    %%

    top_max_Cnaps = cell(6,3);
    
    %%
    
    % Everything below needs to be done in an overall loop, with the correct
    % nsim folder for each of the sims in the selected top level folder
    % iteratively found
    
    % cd(topSimDir)
    % topSimFolders = dir;
    % topSimFolderNum = numel(topSimFolders);
    all_top20maxSels = cell(5,5);
    
    
    % Loop the below code over all of the datasets in the selected folder
    for datasetIndex = 1:datasetNum
        
        if ~any(ismember(Sims_In_Study, all_Sims(datasetIndex)))
            continue;
        else

            cd(dataDir)
            rawData = load(dataDirContents(datasetIndex).name);
            rawDataName = dataDirContents(datasetIndex).name;
            disp(['Working on: ',rawDataName])
            
            if sortedTrialNumbers(nerveIndex) ~= 8
                cd(configDir)
                % Find the config file that matches the current dataset, using the
                % extracted sim no.
                for i = 1:configNum
                    if strcmp(configNames{i},datasetNames{datasetIndex}{2})
                        configInd = i;
                    end
                end
                configFile = table2array(readtable(configDirContents(configInd).name));
            else 
                cd(configDir8)
                % Find the config file that matches the current dataset, using the
                % extracted sim no.
                for i = 1:configNum8
                    if strcmp(configNames8{i},datasetNames{datasetIndex}{2})
                        configInd8 = i;
                    end
                end
                configFile = table2array(readtable(configDir8Contents(configInd8).name));
            end
    
            %% Obtain the simulation data here as well by setting the model and sim directories
            modelDir = fullfile(sampleDir, 'models', num2str(Models_In_Study(1))); % always using the CE cuff for this, so only need the same value
            simDir = fullfile(modelDir, 'sims', num2str(all_Sims(datasetIndex)));
    
            %%
        
            % Data is segmented into three channels, with time unit of seconds and
            % voltage unit of uV. The timestep between subsequent values is 60 us.
            % Channel to fascicle mapping is determined by the experiment. The
            % Channel number and fascicle ID vary per experiment but will
            % always be mapped back to the below arrangement before running the
            % code
            % it's:
            % CH1 - Tibial
            % CH2 - Peroneal
            % CH3 - Sural
            chanNum = 3;
            channelCount = 8; % number of electrodes in array
            
            % Get the names of the fields in the raw data struct
            rawDataFieldnames = fieldnames(rawData);
            rawDataFieldnameNum = size(rawDataFieldnames,1);
            % Split each field name at underscores
            for i = 1:rawDataFieldnameNum
                rawDataFieldnames{i,2} = split(rawDataFieldnames{i},"_");
                rawDataFieldnames{i,3} = size(rawDataFieldnames{i,2},1);
            end
            
            % Find which fields of rawData correspond to channels 1, 2 and 3
            rawDataChannelFields = zeros(chanNum,1);
            for i = 1:chanNum
                % For channel N, the corresponding field will end with _ChN
                endString = ['Ch' num2str(i)];
                for j = 1:rawDataFieldnameNum
                    if strcmp(rawDataFieldnames{j,2}{end},endString)
                        rawDataChannelFields(i) = j;
                    end
                end
            end
            
            % Get the raw data separated out into the three channels values, and their
            % corresponding time vectors (these are offset from each other, presumably
            % due to an interleaving sampling method)
            rawChannelData = cell(chanNum,2);
            channelSamplePeriods = zeros(chanNum,1);
            for i = 1:chanNum
                % Get the channel data from the raw data struct
                channelData = rawData.(rawDataFieldnames{rawDataChannelFields(i),1});
                rawChannelData{i,1} = channelData.times;
                rawChannelData{i,2} = channelData.values;
                % Get the sampling rate of each channel from inverting its "interval"
                % property from the rawData structure
                channelSamplePeriods(i) = channelData.interval;
            end
            
            % From the chosen stimulator config file directory and the name of the data
            % file, load the correct stimulator config file
            
            % The sim number can be found from the second row of the split 
            % rawData struct fieldnames
            simNum = rawDataFieldnames{1,2}{2};
            
            % The stimulator config file will be called
            % "Sim_X_Stimulator_Config.dat"
            configFileName = ['Sim_',simNum,'_Stimulator_Config.dat'];
            
            % Load this file into an array
            if sortedTrialNumbers(nerveIndex) ~= 8
                cd(configDir)
            else
                cd(configDir8)
            end
            configFile = table2array(readtable(configFileName));
            
            % Find the rows of the array where the first value is NaN
            configFileNanRows = false(size(configFile,1),1);
            for i = 1:size(configFile,1)
                if isnan(configFile(i,1))
                    configFileNanRows(i,1) = true;
                end
            end
            configFileNanRowInds = find(configFileNanRows);
    
            %% Loading cruical data for each Stim Configuration
            switch simNum
                case ['300']
                    nsimNum = 80;
                    Number_Poles = 'Bipolar';
                    Phase_Symmetry = 'Even';
                    active_srcs_index = [1,2,3,4,5,6,7,8];
                    active_srcs_cat = [1,2,3,4,5,6,7,8];
                    active_sinks_index = [2,3,4,5,6,7,8,1];
                    active_sinks_cat = [2,3,4,5,6,7,8,1];
                case ['303']
                    nsimNum = 96;
                    Number_Poles = 'Bipolar';
                    Phase_Symmetry = 'Uneven';
                    active_srcs_index = [1,2,3,4,5,6,7,8];
                    active_srcs_cat = [1,2,3,4,5,6,7,8];
                    active_sinks_index = [2,3,4,5,6,7,8,1];
                    active_sinks_cat = [2,3,4,5,6,7,8,1];
                case ['305']
                    nsimNum = 96;
                    Number_Poles = 'Tripolar';
                    Phase_Symmetry = 'Uneven';
                    active_srcs_index = [1,2,3,4,5,6,7,8];
                    active_srcs_cat = [1,2,3,4,5,6,7,8];
                    active_sinks_index = [1,2,3,4,5,6,7,8];
                    active_sinks_cat = [[8,2],[1,3],[2,4],[3,5],[4,6],[5,7],[6,8],[7,1]];
                case ['306']
                    nsimNum = 96;
                    Number_Poles = 'Tripolar';
                    Phase_Symmetry = 'Uneven';
                    active_srcs_index = [1,2,3,4,5,6,7,8];
                    active_srcs_cat = [[8,2],[1,3],[2,4],[3,5],[4,6],[5,7],[6,8],[7,1]];
                    active_sinks_index = [1,2,3,4,5,6,7,8];
                    active_sinks_cat = [1,2,3,4,5,6,7,8];
                case ['307']
                    nsimNum = 128;
                    Number_Poles = 'Bipolar';
                    Phase_Symmetry = 'Uneven';
                    active_srcs_index = [1,2,3,4,5,6,7,8];
                    active_srcs_cat = [1,2,3,4,5,6,7,8];
                    active_sinks_index = [2,3,4,5,6,7,8,1];
                    active_sinks_cat = [2,3,4,5,6,7,8,1];
            end
            
            if pulseParamFileExists
                pulseParamTrends = all_pulseParamTrends{nerveIndex}{datasetIndex};
            else
                %% 
                % Because not all experiments have the same fascicle into the same
                % channel (electrode), you need to rearrange the raw data into Ch1
                % = Tibial, Ch2 = Peroneal, Ch3 = Sural
                % this variable currently only includes Trials 9-12, Trial 8 was
                % triggered differently, so unsure if it will run with this code
                % The 1 aligns to the tibial, 2 aligns to peroneal, and 3 aligns to
                % sural
                %channelToFascicleMap = [1 2 3;3 1 2;1 2 3;1 2 3]
                switch sortedTrialNumbers(nerveIndex)
                    case 8 % nerve trial number 
                        channelToFascicleMap = [1 2 3];
                    case 9
                        channelToFascicleMap = [1 2 3];
                    case 10
                        channelToFascicleMap = [3 2 1];
                    case 11
                        channelToFascicleMap = [1 2 3];
                    case 12
                        channelToFascicleMap = [1 2 3];
                end
                % Create a new cell array and array to hold the rearranged data
                rearrangedChannelData = cell(chanNum,2);
                rearrangedSamplePeriods = zeros(chanNum,1);
                
                for i = 1:chanNum
                    rearrangedChannelData{ channelToFascicleMap(i), 1 } = rawChannelData{i, 1};
                    rearrangedChannelData{ channelToFascicleMap(i), 2 } = rawChannelData{i, 2};
                    rearrangedSamplePeriods( channelToFascicleMap(i) ) = channelSamplePeriods(i);
                end
                
                % Now overwrite the original once we've finished
                rawChannelData = rearrangedChannelData;
                channelSamplePeriods = rearrangedSamplePeriods;
                
                % Get trigger times to clip data
                % first get the name of trigger channel and ID data in raw data structure
                triggerEndString = 'Ch4';
                for i = 1:rawDataFieldnameNum
                    rawDataFieldnames{i,2} = split(rawDataFieldnames{i},"_");
                    rawDataFieldnames{i,3} = size(rawDataFieldnames{i,2},1);
                    if strcmp(rawDataFieldnames{i,2}{end},triggerEndString)
                        triggerTimesDataField = i;
                    end
                end
                % isolate trigger time data into cell
                triggerTimes = cell(1,1)
                triggerTimes{1,1} = rawData.(rawDataFieldnames{triggerTimesDataField,1}).times
                %
                
                % Check where gaps in the trigger data
                for i = 1:length(triggerTimes{:})-1
                    trigDiff(i) = round(triggerTimes{:}(i+1)-triggerTimes{:}(i),5);
                end
                size04 = sum(round(trigDiff,3) == 0.04)
                size08 = sum(round(trigDiff,3) == 0.08)
                not0408Ind = find(round(trigDiff,3)~=0.04 & round(trigDiff,3)~=0.08)
                not0408Val = trigDiff(not0408Ind)
                %%
        
                % The channel data sets can be found from the non-sequential
                % discontinuities in configFileNanRowInds. If there is a discontinuity
                % between index i and index j, indices i+1:j-1 contain the channel data
                channelConfigData = cell(channelCount,1);
                configFileNanRowIndDisc = find((configFileNanRowInds-circshift(configFileNanRowInds,1))>1);
                for i = 1:channelCount
                    channelConfigData{i} = configFile(configFileNanRowInds(configFileNanRowIndDisc(i)-1)+1:configFileNanRowInds(configFileNanRowIndDisc(i))-1,:);
                end
                
                
                channelPulseTrains = cell(channelCount,1);
                for i = 1:channelCount
                    channelPulseTrains{i} = zeros(sum(channelConfigData{i}(:,2)),1);
                    startInd = 1;
                    for j = 1:size(channelConfigData{i},1)
                        channelPulseTrains{i}(startInd:startInd+channelConfigData{i}(j,2)-1) = ones(channelConfigData{i}(j,2),1)*channelConfigData{i}(j,1);
                        startInd = startInd+channelConfigData{i}(j,2);
                    end
                end
                
                % The channel pulse trains can be downsampled by a factor of 20, which was
                % the minimum time step, without losing any time resolution
                channelPulseTrainsDownsampled = cell(channelCount,1);
                downsampleRate = 20;
                for i = 1:channelCount
                    channelPulseTrainsDownsampled{i} = downsample(channelPulseTrains{i},downsampleRate);
                end
                pulseTrainDownsampleTime = 0:downsampleRate/1e6:size(channelPulseTrains{1},1)/1e6-downsampleRate/1e6;
                
                %%
                
                % By rolling a delayed version of the pulse train over the base version, 
                % find discontinuities. 
                % BELOW METHOD ONLY WORKS FOR CATHODIC LEADING PULSES
                channelPulseDisc = cell(channelCount,3);
                channelPulseNum = cell(channelCount,1);
                channelPulseParams = cell(channelCount,1);
                for i = 1:channelCount
                    % Find the discontinuities in the channel pulse train, padding with a
                    % zero to capture the first falling edge if it occurs on the first index
                    paddedChannelPulseTrain = [0; channelPulseTrains{i}];
                    channelPulseDisc{i,1} = paddedChannelPulseTrain(2:end)-paddedChannelPulseTrain(1:end-1);
                    % Find the falling edges in the pulse train
                    channelPulseDisc{i,2} = find(channelPulseDisc{i,1}<0);
                    % Find the rising edges in the pulse train
                    channelPulseDisc{i,3} = find(channelPulseDisc{i,1}>0);
                    % For each channel, the number of pulses is half the number of falling
                    % edges
                    channelPulseNum{i} = size(channelPulseDisc{i,2},1)/2;
                    channelPulseParams{i} = zeros(channelPulseNum{i},5);
                    riseStartInd = 1;
                    for j = 1:channelPulseNum{i}
                        % Find the jth and (j+1)th falling and rising edges
                        fall1 = channelPulseDisc{i,2}(1+(j-1)*2);
                        fall2 = channelPulseDisc{i,2}(2+(j-1)*2);
                        rise1 = channelPulseDisc{i,3}(riseStartInd);
                        rise2 = channelPulseDisc{i,3}(riseStartInd+1);
                        % The pulse beings on the first falling edge
                        channelPulseParams{i}(j,1) = fall1;
                        % The injection pulse width is always the first rise time minus the
                        % first fall time
                        channelPulseParams{i}(j,2) = rise1-fall1;
                        % Find if there's an interphase gap by seeing if the second rise is
                        % between the first and second falls
                        if rise2 < fall2 && rise2 > fall1
                            % If there is, the interphase gap is given by the second rise
                            % time minus the first rise time
                            channelPulseParams{i}(j,3) = rise2-rise1;
                            % The return pulse width is then given by the second fall time
                            % minus the second rise time
                            channelPulseParams{i}(j,4) = fall2-rise2;
                            % Two rises have been used, so increment the rise start index
                            % by 2
                            riseStartInd = riseStartInd+2;
                        else
                            % Else the interphase gap is zero
                            channelPulseParams{i}(j,3) = 0;
                            % Then the return pulse width is given by the second fall time
                            % minus the first rise time
                            channelPulseParams{i}(j,4) = fall2-rise1;
                            % One rise has been used, so increment the rise start index
                            % by 1
                            riseStartInd = riseStartInd+1;
                        end
                        % The magnitude of the pulse is found by the value of the odd
                        % falling edges
                        channelPulseParams{i}(j,5) = channelPulseDisc{i,1}(channelPulseDisc{i,2}(1+(j-1)*2));
                    end
                end
        
                %% Get the pulse parameters per Stim Configuration
                
                % Find the order and number of trends within the pulse parameters
                pulseParamTrends = cell(channelCount,4);
                for i = 1:channelCount
                    for j = 1:4
                        pulseParamTrends{i,j} = cell(2,1);
                        
                        % Get unique values
                        uniqueValues = unique(channelPulseParams{i}(:,j+1),'stable');
                        
                        % If there's only one unique value, repeat it 'channelCount' times
                        if length(uniqueValues) == 1
                            pulseParamTrends{i,j}{1} = repmat(uniqueValues, nsimNum/channelCount, 1); 
                        elseif length(uniqueValues) == 2 & all(uniqueValues) == 0
                            pulseParamTrends{i,j}{1} = repmat(uniqueValues(1), nsimNum/channelCount, 1); 
                        else
                            pulseParamTrends{i,j}{1} = uniqueValues; % Use normally extracted unique values
                        end
                        
                        % Store the length (number of unique values)
                        pulseParamTrends{i,j}{2} = length(pulseParamTrends{i,j}{1});
                    end
                end
        
                all_pulseParamTrends{nerveIndex}{datasetIndex} = pulseParamTrends;
            end
                
                %%
    
            if allAveragedDataFileExists
                averagedChannelData = all_averagedChannelData{nerveIndex}{datasetIndex};
            else
                
                % Find the discontinuities in the time data, which represent gaps between
                % triggers
                rawChannelDisc = cell(size(rawChannelData,1),1);
                rawChannelTimeDiffs = cell(size(rawChannelData,1),1);
                % For each channel
                for i = 1:size(rawChannelData,1)
                    % Find the difference between subsequent raw channel times, minus the
                    % sampling period
                    timeDiffs = rawChannelData{i,1}(2:end) - rawChannelData{i,1}(1:end-1);
                    rawChannelTimeDiffs{i} = timeDiffs;
                    rawChannelDisc{i} = find(timeDiffs>1.1*channelSamplePeriods(i));
                end
                
                %%
    
                % Find points in the data where there is at least a 40 ms gap between
                % subsequent samples. These points will exist if the contact configuration
                % is not tripolar, due to the method used to generate the trigger wave (for
                % tripolar, it was created deterministically and so triggers for zero
                % amplitude stim. For bipolar, it was generated dynamically and so doesn't
                % trigger for zero amplitude stim).
                singleLossInds = cell(size(rawChannelTimeDiffs));
                for i = 1:size(rawChannelTimeDiffs,1)
                    singleLossInds{i} = find(rawChannelTimeDiffs{i}>0.04);
                end
    
                % If there are segments with at least a 40 ms gap between them
                if ~isempty(singleLossInds{1})
                    % Find the difference in index between subsequent gaps to see if there are
                    % breaks in the pattern
                    singleLossIndDiffs = cell(size(rawChannelTimeDiffs));
                    for i = 1:1:size(rawChannelTimeDiffs,1)
                        singleLossIndDiffs{i} = singleLossInds{i}(2:end)-singleLossInds{i}(1:end-1);
                    end
    
    
                    %%
    
                    % As an algorithm, one way to extract the missing
                    % segments would be to:
                    % 1. Set burst one's integrity to good
                    % 2. Find the location of the first zero amp gap in the sample times, which
                    % be the burst one gap
                    % 3. Find the location of the next burst's zero amp gap
                    % 4. If there are 16000 samples (32 500 sample segments) between the
                    % current and next burst gaps, set the integrity of the next burst to good
                    % 5. Iterate the above across all bursts in the data
                    burstIntegrities = cell(size(rawChannelData,1),1);
                    for i = 1:size(rawChannelData,1)
                        burstIntegrities{i} = true;
                        currentBurst = 1;
                        while currentBurst*32*500 < length(rawChannelData{i,1})
                            if singleLossIndDiffs{1}(currentBurst) == 16000
                                burstIntegrities{i} = [burstIntegrities{i}; true];
                            else
                                burstIntegrities{i} = [burstIntegrities{i}; false];
                            end
                            currentBurst = currentBurst + 1;
                        end
                    end
    
                    % If burstIntegrities{i}(j) is false, then on channel i's data the burst
                    % between zero gap j-1 and zero gap j is missing a segment
                    repeatKey = cell(size(burstIntegrities));
                    repeats = 5;
                    nsimNum = length(burstIntegrities{1})/repeats;
                    if mod(nsimNum,1) ~= 0
                        disp('Issue with data integrity - number of segment bursts not divisible by repeat number')
                    else
                        for i = 1:size(repeatKey,1)
                            nsim = 1;
                            repeatKey{i} = zeros(size(burstIntegrities{i}));
                            for j = 1:length(repeatKey{i})
                                if burstIntegrities{i}(j)
                                    repeatKey{i}(j) = nsim;
                                end
                                % Every five bursts, increment the nsimNum
                                if mod(j,5) == 0
                                    nsim = nsim + 1;
                                end
                            end
                        end
                    end
    
                    % Separate the raw data into sets of bursts that are repeats of each other,
                    % omitting bursts that only contain 31 triggered segments of data
                    cutChannelData = cell(size(rawChannelData));
                    for i = 1:size(rawChannelData,1)
                        cutChannelData{i,1} = cell(nsimNum,1);
                        cutChannelData{i,2} = cell(nsimNum,1);
                        currentBurst = 0;
                        currentSegment = 0;
                        for j = 1:nsimNum
                            cutChannelData{i,1}{j} = cell(5,1);
                            cutChannelData{i,2}{j} = cell(5,1);
                            for k = 1:5
                                currentBurst = currentBurst + 1;
                                if burstIntegrities{i}(currentBurst)
                                    cutChannelData{i,1}{j}{k} = cell(pulseParamTrends{i,4}{2},1);
                                    cutChannelData{i,2}{j}{k} = cell(pulseParamTrends{i,4}{2},1);
                                    for m = 1:pulseParamTrends{i,4}{2}
                                        if currentSegment == 0
                                            cutChannelData{i,1}{j}{k}{m} = rawChannelData{i,1}(1:rawChannelDisc{i}(1));
                                            cutChannelData{i,2}{j}{k}{m} = rawChannelData{i,2}(1:rawChannelDisc{i}(1));
                                        elseif currentSegment == length(rawChannelDisc{i})
                                            cutChannelData{i,1}{j}{k}{m} = rawChannelData{i,1}(rawChannelDisc{i}(currentSegment)+1:end);
                                            cutChannelData{i,2}{j}{k}{m} = rawChannelData{i,2}(rawChannelDisc{i}(currentSegment)+1:end);
                                        else
                                            cutChannelData{i,1}{j}{k}{m} = rawChannelData{i,1}(rawChannelDisc{i}(currentSegment)+1:rawChannelDisc{i}(currentSegment+1));
                                            cutChannelData{i,2}{j}{k}{m} = rawChannelData{i,2}(rawChannelDisc{i}(currentSegment)+1:rawChannelDisc{i}(currentSegment+1));
                                        end
                                        currentSegment = currentSegment + 1;
                                    end
                                else
                                    currentSegment = currentSegment + pulseParamTrends{i,4}{2} - 1;
                                end
                            end
                        end
                    end
    
                % Else if there are no segments with a 40 ms gap between them, the zero
                % amplitude stims were triggered and so need to be skipped in the data
                % cutting. The relative location of the zero amplitude stims can be found
                % from channelPulseParams, since the first column of each cell contains the
                % time stamps of each non-zero stim - the gaps in this can then show where
                % a stim was at zero amplitude. This location only needs to be found for
                % the first set of amplitudes, since the amplitude order is consistent
                else
                    baseStimTimes = channelPulseParams{1}(1:pulseParamTrends{1,4}{2},1);
                    baseStimTimeDiffs = baseStimTimes(2:end)-baseStimTimes(1:end-1);
                    zeroAmpLoc = false(pulseParamTrends{1,4}{2}+1,1);
                    zeroAmpLoc(find(baseStimTimeDiffs>40000)+1) = true;
    
                    repeats = 5;
                    cutChannelData = cell(size(rawChannelData));
                    for i = 1:size(rawChannelData,1)
                        nsimNum = (size(rawChannelDisc{i},1)+1)/(repeats*(pulseParamTrends{1,4}{2}+1));
                        currentSegment = 0;
                        cutChannelData{i,1} = cell(nsimNum,1);
                        cutChannelData{i,2} = cell(nsimNum,1);
                        for j = 1:nsimNum
                            cutChannelData{i,1}{j} = cell(repeats,1);
                            cutChannelData{i,2}{j} = cell(repeats,1);
                            for k = 1:repeats
                                cutChannelData{i,2}{j}{k} = cell(pulseParamTrends{i,4}{2},1);
                                ampInd = 1;
                                for m = 1:pulseParamTrends{i,4}{2}+1
                                    if m ~= find(zeroAmpLoc)
                                         if currentSegment == 0
                                             cutChannelData{i,1}{j}{k}{ampInd} = rawChannelData{i,1}(1:rawChannelDisc{i}(1));
                                             cutChannelData{i,2}{j}{k}{ampInd} = rawChannelData{i,2}(1:rawChannelDisc{i}(1));
                                         elseif currentSegment == length(rawChannelDisc{i})
                                             cutChannelData{i,1}{j}{k}{ampInd} = rawChannelData{i,1}(rawChannelDisc{i}(currentSegment)+1:end);
                                             cutChannelData{i,2}{j}{k}{ampInd} = rawChannelData{i,2}(rawChannelDisc{i}(currentSegment)+1:end);
                                         else
                                             cutChannelData{i,1}{j}{k}{ampInd} = rawChannelData{i,1}(rawChannelDisc{i}(currentSegment)+1:rawChannelDisc{i}(currentSegment+1));
                                             cutChannelData{i,2}{j}{k}{ampInd} = rawChannelData{i,2}(rawChannelDisc{i}(currentSegment)+1:rawChannelDisc{i}(currentSegment+1));
                                         end
                                         ampInd = ampInd + 1;
                                    end
                                    currentSegment = currentSegment + 1;
                                end
                            end
                        end
                    end
                end
        
                %%
    
                % Iterating over the third level cell of cutChannelData iterates over the
                % segments that are repeats of each other
                % Average segments that are repeats. There are guaranteed to be 5 cells in
                % the third level, but one may be empty (0x0 double)
                averagedChannelData = cell(size(cutChannelData));
                for i = 1:size(averagedChannelData,1)
                    averagedChannelData{i,1} = cell(size(cutChannelData{i,1}));
                    averagedChannelData{i,2} = cell(size(cutChannelData{i,2}));
                    for j = 1:size(averagedChannelData{i,1},1)
                        averagedChannelData{i,1}{j} = cell(pulseParamTrends{i,4}{2},1);
                        averagedChannelData{i,2}{j} = cell(pulseParamTrends{i,4}{2},1);
                        for k = 1:pulseParamTrends{i,4}{2}
                            repeats = 5;
                            for m = 1:5
                                tempData = cutChannelData{i,1}{j}{m};
                                if isempty(tempData)
                                    i
                                    j
                                    k
                                    m
                                    repeats = 4;
                                else
                                    tempData1 = tempData{k};
                                    tempData2 = cutChannelData{i,2}{j}{m}{k};
                                    if isempty(averagedChannelData{i,1}{j}{k})
                                        % This is here becuase each segment should be the
                                        % same length, eg: 500 samples, but some are longer
                                        % (eg: 517). No idea why, may be related to the
                                        % issue with dropped segments. Coding around it for
                                        % now.
                                        averagedChannelData{i,1}{j}{k} = tempData1(1:length(cutChannelData{1,1}{1}{1}{1}));
                                        averagedChannelData{i,2}{j}{k} = tempData2(1:length(cutChannelData{1,1}{1}{1}{1}));
                                    else
                                        averagedChannelData{i,1}{j}{k} = averagedChannelData{i,1}{j}{k} + tempData1(1:length(cutChannelData{1,1}{1}{1}{1}));
                                        averagedChannelData{i,2}{j}{k} = averagedChannelData{i,2}{j}{k} + tempData2(1:length(cutChannelData{1,1}{1}{1}{1}));
                                    end
                                end
                            end
                            averagedChannelData{i,1}{j}{k} = averagedChannelData{i,1}{j}{k}/repeats;
                            averagedChannelData{i,2}{j}{k} = averagedChannelData{i,2}{j}{k}/repeats;
                        end
                    end
                end
        
                all_averagedChannelData{nerveIndex}{datasetIndex}=averagedChannelData;
            end
    
            %% Create a single vector that holds the amplitudes of each cut data segment
            ampNum = 33; %total number of amplitudes including 0, corrected in code to use 32
            repeats = 5;
            segmentAmplitudes = zeros(ampNum-1,repeats*nsimNum*(ampNum-1));
            for i = 1:size(segmentAmplitudes,2)
                segmentAmplitudes(:,i) = pulseParamTrends{1,4}{1}; 
            end
    
            segmentAmplitudeArray = reshape(segmentAmplitudes,numel(segmentAmplitudes),1);
    
            %% Need to exclude all data points that are above 300 us Pulse Widths 
            % pulse width threshold to ensure that the stimulation artefact does not 
            % have an impact on the CMAP calcuation so I must indicate 
            % which nsims to exclude per dataset, must have the google 
            % spreadsheet Simulation Logbook to match the nsims to waveforms  
    
            switch simNum
                case ['300']
                    nSimsToExclude = [80 79 78 77 70 69 68 67 60 59 58 57 50 49 48 47 40 39 38 37 30 29 28 27 20 19 18 17 10 9 8 7];
                case ['303']
                    nSimsToExclude = [96 95 94 93 92 84 83 82 81 80 72 71 70 69 68 60 59 58 57 56 48 47 46 45 44 36 35 34 33 32 24 23 22 21 20 12 11 10 9 8];
                case ['305']
                    nSimsToExclude = [96 95 94 93 92 84 83 82 81 80 72 71 70 69 68 60 59 58 57 56 48 47 46 45 44 36 35 34 33 32 24 23 22 21 20 12 11 10 9 8];
                case ['306']
                    nSimsToExclude = [96 95 94 93 92 84 83 82 81 80 72 71 70 69 68 60 59 58 57 56 48 47 46 45 44 36 35 34 33 32 24 23 22 21 20 12 11 10 9 8];
                case ['307']
                    nSimsToExclude = [128 127 126 125 112 111 110 109 96 95 94 93 80 79 78 77 64 63 62 61 48 47 46 45 32 31 30 29 16 15 14 13];
            end
        
            nsimNumPostExclusion = nsimNum - numel(nSimsToExclude);
            
            %%
            
            % Zero mean and smooth the averaged data with a Gaussian window, removing the
            % first 100 samples and last 225 samples as they're guaranteed to only 
            % contain noise
            windowWidth = 5;
            window = gausswin(windowWidth);
            normWindow = window/sum(window);
            smoothedChannelData = averagedChannelData;
            % making different data start and end times for each dataset and
            % nerve index to better optimize the peak finding process, doing
            % this in a separate file to reduce the number of lines in code
            for i = 1:size(averagedChannelData,1)
                for j = 1:size(averagedChannelData{i,2},1)
                    %set start and end times in this nest bc it does not depend
                    %on amplitude, just channel number and waveform which is
                    %all contained within nsimNum
                    dataStartArray = calculateAllDataStartTimes(sortedTrialNumbers,nerveIndex,simNum,nsimNum,channelCount,j,nSimsToExclude);
                    dataEnd = dataStartArray(i,mod(j-1,(nsimNum/channelCount))+1) + 40;
                    for k = 1:size(averagedChannelData{i,2}{j},1)          
                        smoothedChannelData{i,2}{j}{k} = conv(averagedChannelData{i,2}{j}{k}(dataStartArray(i,mod(j-1,(nsimNum/channelCount))+1):dataEnd)-mean(averagedChannelData{i,2}{j}{k}(dataStartArray(i,mod(j-1,(nsimNum/channelCount))+1):dataEnd)),normWindow,'same');
                    end
                end
            end
    
            for i = 1:size(averagedChannelData,1)
                for j = 1:size(averagedChannelData{i,2},1)
                    for k = 1:size(averagedChannelData{i,2}{j},1)
                        if all_Nerves(nerveIndex) == 8 && (simNum == string('300') || simNum == string('303') || simNum == string('307'))
                            averagedChannelData{i,1}{j}{k,2} = -segmentAmplitudeArray(k+(j-1)*(ampNum-1)*5-(floor(j/(nsimNum/channelCount))));
                            averagedChannelData{i,2}{j}{k,2} = -segmentAmplitudeArray(k+(j-1)*(ampNum-1)*5-(floor(j/(nsimNum/channelCount))));
                            smoothedChannelData{i,1}{j}{k,2} = -segmentAmplitudeArray(k+(j-1)*(ampNum-1)*5-(floor(j/(nsimNum/channelCount))));
                            smoothedChannelData{i,2}{j}{k,2} = -segmentAmplitudeArray(k+(j-1)*(ampNum-1)*5-(floor(j/(nsimNum/channelCount))));
                        else
                            averagedChannelData{i,1}{j}{k,2} = -segmentAmplitudeArray(k+(j-1)*(ampNum-1)*5);
                            averagedChannelData{i,2}{j}{k,2} = -segmentAmplitudeArray(k+(j-1)*(ampNum-1)*5);
                            smoothedChannelData{i,1}{j}{k,2} = -segmentAmplitudeArray(k+(j-1)*(ampNum-1)*5);
                            smoothedChannelData{i,2}{j}{k,2} = -segmentAmplitudeArray(k+(j-1)*(ampNum-1)*5);
                        end
                    end
                end
            end
    
            % Reorder the averaged, zeroed, filtered channel data by injection current amplitude within
            % each set of segments from each nsim
            for i = 1:chanNum
                for j = 1:nsimNum
                    reorderedAveragedChannelData{i,1}{j} = sortrows(averagedChannelData{i,1}{j},2,'ascend');
                    reorderedAveragedChannelData{i,2}{j} = sortrows(averagedChannelData{i,2}{j},2,'ascend');
                    reorderedSmoothedChannelData{i,1}{j} = sortrows(smoothedChannelData{i,1}{j},2,'ascend');
                    reorderedSmoothedChannelData{i,2}{j} = sortrows(smoothedChannelData{i,2}{j},2,'ascend');
                end
            end

            reorderedStimCycle = reorderedSmoothedChannelData{1,2}{1,1}(:,2);
            all_reorderedSmoothedChannelDat{nerveIndex}{datasetIndex} = smoothedChannelData;
                    
            %%
            
            % From a priori knowledge of the number of contacts and arrangement of the
            % input wave, create a mapping from nsim segments to contact and waveform
            % parameter numbers
            contactNum = 8;
            waveformNum = nsimNum/contactNum;
            waveformNumPostExclusion = nsimNumPostExclusion/contactNum;
            nsimList = repelem(1:nsimNum, 33)';   % 1 1 … 1  (33×)  2 2 … 2 (33×)  ... nsimNum
            ampListbyNsim = repmat([0; sort(-pulseParamTrends{1,4}{1})] , nsimNum , 1);   % transpose, replicate, keep column shape
            segmentNsimParameterMapping = [nsimList,ampListbyNsim,zeros(size(nsimList,1),6)];
            for i = 1:size(nsimList,1)
                segmentNsimParameterMapping(i,3) = ceil(segmentNsimParameterMapping(i,1)/waveformNum);
                waveform = segmentNsimParameterMapping(i,1);
                while waveform > waveformNum
                    waveform = waveform-waveformNum;
                end
                segmentNsimParameterMapping(i,4) = waveform;
                % Find the value of the changing waveform parameter automatically
                if pulseParamTrends{1,1}{1,1}(1) ~= pulseParamTrends{1,1}{1,1}(2)
                    segmentNsimParameterMapping(i,5) = pulseParamTrends{1,1}{1,1}(waveform);
                elseif pulseParamTrends{2,1}{1,1}(1) ~= pulseParamTrends{1,2}{1,1}(2)
                    segmentNsimParameterMapping(i,5) = pulseParamTrends{1,2}{1,1}(waveform);
                else
                    segmentNsimParameterMapping(i,5) = pulseParamTrends{1,3}{1,1}(waveform);
                end
                % Find the injection period for later calculation of the cropping
                % window
                segmentNsimParameterMapping(i,6) = pulseParamTrends{1,1}{1,1}(waveform);
                segmentNsimParameterMapping(i,7) = pulseParamTrends{1,2}{1,1}(waveform);
                segmentNsimParameterMapping(i,8) = pulseParamTrends{1,3}{1,1}(waveform);
            end
        
            %% Calculating Cnaps
            channelSegmentCmaps = cell(size(averagedChannelData,1),1);
            channelMaxCmaps_allNSims = cell(size(averagedChannelData,1),1);
            channelMaxCmaps_1NSim = cell(size(averagedChannelData,1),nsimNum);
            % For each low-pass filtered segment, find the cmap
            for i = 1:size(averagedChannelData,1)
                channelMaxCmaps_allNSims{i} = [0 0 0];
                for j = 1:nsimNum
                    channelMaxCmaps_1NSim{i,j} = [0 0];
                    for k = 1:pulseParamTrends{i,4}{2}
                        % Find the 2 largest positive and negative peaks 
                        cmapMax = max(reorderedSmoothedChannelData{i,2}{j}{k});
                        cmapMin = min(reorderedSmoothedChannelData{i,2}{j}{k});
                        cmapPP = cmapMax-cmapMin;
                        if abs(cmapMax) < 0.15*abs(cmapMin) | abs(cmapMin) < 0.25*abs(cmapMax)
                            cmapPP = 0;
                        end
                        channelSegmentCmaps{i}{j}{k} = cmapPP;
                        channelSegmentCmapsPosPeak{i}{j}{k} = cmapMax;
                        channelSegmentCmapsNegPeak{i}{j}{k} = cmapMin;
                        if cmapPP > channelMaxCmaps_allNSims{i}(1)
                            channelMaxCmaps_allNSims{i} = [cmapPP j k];
                        end
                        if cmapPP > channelMaxCmaps_1NSim{i,j}(1)
                            channelMaxCmaps_1NSim{i,j} = [cmapPP k];
                        end
                    end
                end
            end
    
            % Create utility variables for plotting
    
            fascicleColours = [1,0,0;0,1,0;0,0,1];
        
            % Plotting the specific sim of all max Cnaps to determine which one should
            % be used as the maximum for normalization
    
            % If this has been run and the values have been saved, I will just
            % import them and not have to reselect them everytime, otherwise,
            % the code will run again
            
            if maxAmpsForNormFileExists
                maxAmpsForNormUserDet = maxAmpsForNormUserDet_allTrials{nerveIndex,datasetIndex};
            else % determine the maxmized values manually
                % unrolling the segment cmaps (as you already do)
                unrolledSegmentCmaps = cell(chanNum,1);
                for i = 1:chanNum
                    for j = 1:nsimNum
                        if j ~= nSimsToExclude
                            unrolledSegmentCmaps{i} = [unrolledSegmentCmaps{i} cell2mat(channelSegmentCmaps{i}{j})];
                        end
                    end
                end
        
                figure(figNum);
                title('All Vpp of Cnaps in Sim');
                len_stim = size(unrolledSegmentCmaps{1},2);
                hold on
                plot(1:len_stim, sort(unrolledSegmentCmaps{1}), 'r')
                plot(1:len_stim, sort(unrolledSegmentCmaps{2}), 'g')
                plot(1:len_stim, sort(unrolledSegmentCmaps{3}), 'b')
                xlabel('Configuration in Sim (Unique Combo of Electrode, Amp, and Waveform)')
                ylabel('Amplitude (uV)')
                hold off
        
                % We’ll store the maximum amplitude (as “floor(y)”) for each of the 3 lines.
                maxAmpsForNormUserDet = zeros(3,1);
                % Also store the config indices (x-coordinates).
                configInd = zeros(3,1);
        
                % Optional: give descriptive names for each of the 3 lines
        
                % Loop through each of the 3 lines
                for lineIdx = 1:3
                    offsetFactor = 0;
                    offset = false;
                    confirmed = false;
                    while ~confirmed
                        % Ask user to click on the main figure
                        figure(figNum);
        
                        try
        
                            if ~offset
                                disp(['Click on the figure to select ', fascicleNames{lineIdx}, ' maximum amplitude...']);
                                [xClick, yClick] = ginput(1);
        
                                % Convert xClick and yClick as you do:
                                % (floor for y, plus +1 for x to convert from 0-based if you want)
                                configInd(lineIdx) = floor(xClick) + 1;
                                maxAmpsForNormUserDet(lineIdx) = floor(yClick);
                            else
                                maxAmpsForNormUserDet(lineIdx) = maxAmpsForNormUserDet(lineIdx)+offsetFactor;
                            end
        
                            % Locate the config and amp closest to this Cmap
                            % Example: We clicked on the figure for channel "lineIdx"
                            bigCell = channelSegmentCmaps{lineIdx};  % 80x1 cell array, each cell has 32 values
        
                            % Initialize tracking
                            minDiff = Inf;
                            closestVal = NaN;
                            bestConfigIndex = NaN;
                            bestAmpIndex = NaN;
        
                            % Loop through all 80 * 32 values
                            for j = 1:nsimNum
                                for k = 1:32
                                    val = bigCell{j}{k};  % the actual numeric value
                                    thisDiff = abs(val - maxAmpsForNormUserDet(lineIdx));
                                    if thisDiff < minDiff
                                        minDiff = thisDiff;
                                        closestVal = val;
                                        bestConfigIndex = j;   % which "row" in the 80
                                        bestAmpIndex = k;    % which "subcell" in the 32
                                    end
                                end
                            end
        
                            % Now open a verification figure            
                            % Plot only the line we’re verifying
                            figure('Position',[560 200 550 350])  %To have figure pop on the left side of the screen instead of the middle, easier for visualization
                            clf; 
                            plot(reorderedSmoothedChannelData{1,2}{bestConfigIndex}{bestAmpIndex},'r', 'DisplayName','Tibial')
                            hold on
                            plot(reorderedSmoothedChannelData{2,2}{bestConfigIndex}{bestAmpIndex},'g', 'DisplayName','Peroneal')
                            plot(reorderedSmoothedChannelData{3,2}{bestConfigIndex}{bestAmpIndex},'b', 'DisplayName','Sural')
                            title('Smoothed')
                            yline(channelSegmentCmapsPosPeak{1}{bestConfigIndex}{bestAmpIndex},'Color',[fascicleColours(1,:) 0.5], 'HandleVisibility','off')
                            yline(channelSegmentCmapsNegPeak{1}{bestConfigIndex}{bestAmpIndex},'Color',[fascicleColours(1,:) 0.5], 'HandleVisibility','off')
                            yline(channelSegmentCmapsPosPeak{2}{bestConfigIndex}{bestAmpIndex},'Color',[fascicleColours(2,:) 0.5], 'HandleVisibility','off')
                            yline(channelSegmentCmapsNegPeak{2}{bestConfigIndex}{bestAmpIndex},'Color',[fascicleColours(2,:) 0.5], 'HandleVisibility','off')
                            yline(channelSegmentCmapsPosPeak{3}{bestConfigIndex}{bestAmpIndex},'Color',[fascicleColours(3,:) 0.5], 'HandleVisibility','off')
                            yline(channelSegmentCmapsNegPeak{3}{bestConfigIndex}{bestAmpIndex},'Color',[fascicleColours(3,:) 0.5], 'HandleVisibility','off')
        
                            title(['Verification: ', fascicleNames{lineIdx}]);
                            legend show; grid on;
                            hold off;
        
                            % Ask user to confirm
                            userInput = input('Confirm selection (Y) or redo (R)? ','s');
        
                            if strcmpi(userInput, 'Y')
                                confirmed = true;
                                disp(['Confirmed for ', fascicleNames{lineIdx}, '.']);
                            elseif strcmpi(userInput, 'R')
                                disp(['Redoing selection for ', fascicleNames{lineIdx}, '...']);
                                offset = false;
                            elseif strcmpi(userInput, 'U')
                                disp(['Incrementing selection for ', fascicleNames{lineIdx}, '...']);
                                offsetFactor = offsetFactor + 1000;
                                offset = true;
                            else
                                disp(['Decrementing selection for ', fascicleNames{lineIdx}, '...']);
                                offsetFactor = offsetFactor - 1000;
                                offset = true; 
                            end
        
                            % If user closed the figure, ginput() might error out, so handle that:
                        catch ME
                            disp(['Error occurred or figure was closed: ', ME.message]);
                            break;
                        end
                    end
        
                    % If the figure was closed or an error occurred, we can decide to either
                    % break the whole process or keep going. Here we break:
                    if ~ishandle(figNum)
                        disp('Main figure was closed; stopping.');
                        break;
                    end
                end
        
                % Move on to the next figure number, etc.
                figNum = figNum + 1;
        
                % save each max data to the allTrials cell matrix, delete after
                % this run
                maxAmpsForNormUserDet_allTrials{nerveIndex,datasetIndex}=maxAmpsForNormUserDet
            end
    
            % Normalise the channel segment cmaps by the channel max
            % cmaps determined by user
            channelSegNormCmaps_allNSims = channelSegmentCmaps;
            for i = 1:size(averagedChannelData,1)
                for j = 1:nsimNum
                    for k = 1:pulseParamTrends{i,4}{2}
                        if j ~= nSimsToExclude
                            channelSegNormCmaps_allNSims{i}{j}{k} = channelSegNormCmaps_allNSims{i}{j}{k}/(maxAmpsForNormUserDet(i));
                        else 
                            channelSegNormCmaps_allNSims{i}{j}{k} = [];
                        end
                    end
                end
            end
            % Normalise the cnaps using the user determined per-channel maximum recorded cnaps
        
            % replacing all values above 1 with 0 that indicate non-significant Vpp values
            % (errant maximums)
     
            for i = 1:chanNum
                for j = 1:nsimNum
                    for k = 1:pulseParamTrends{i,4}{2}
                        if channelSegNormCmaps_allNSims{i}{j}{k}>1
                           channelSegNormCmaps_allNSims{i}{j}{k}
                           channelSegNormCmaps_allNSims{1}{j}{k} = 0;
                           channelSegNormCmaps_allNSims{2}{j}{k} = 0;
                           channelSegNormCmaps_allNSims{3}{j}{k} = 0;
                        end
                    end
                end
            end
    
            all_channelSegNormCmaps{nerveIndex}{datasetIndex} = channelSegNormCmaps_allNSims;
    
            % unrolling the segment cmaps to detect top 10 maximums easily and
            % to replot them to see differences
            unrolledFilteredSegNormCmaps = cell(chanNum,1);
            maxCmapValues = cell(chanNum,1);
            maxCmapConfigs = cell(chanNum,1);
            for i = 1:chanNum
                for j = 1:nsimNum
                    unrolledFilteredSegNormCmaps{i} = [unrolledFilteredSegNormCmaps{i} cell2mat(channelSegNormCmaps_allNSims{i}{j})];
                end
                [maxCmapValues{i,1},maxCmapConfigs{i,1}] = maxk(unrolledFilteredSegNormCmaps{i}(:),10);
                top_max_Cnaps{datasetIndex,i} = [maxCmapValues(i,:),maxCmapConfigs(i,:)];
            end
        
            % (Optional) Replotting the used normalized Cnaps
            % figure(figNum)
            % figNum = figNum+1;
            % title('All Vpp of Cnaps in Sim')
            % len_stim = size(unrolledFilteredSegNormCmaps{1},2);
            % hold on
            % plot(1:len_stim,sort(unrolledFilteredSegNormCmaps{1}),'r')
            % plot(1:len_stim,sort(unrolledFilteredSegNormCmaps{2}),'g')
            % plot(1:len_stim,sort(unrolledFilteredSegNormCmaps{3}),'b')
            % xlabel('Configuration in Sim (Unique Combo of Electrode, Amp, and Waveform)')
            % ylabel('Amplitude (uV)')
            % hold off
        
            %% Calculate Selectivities
    
            % Calculate the selectivities for each nsim/amplitude, for each channel
            % being the target, with the same off-target activation penalty as in the
            % simulation selectivity calculations
            
            % For averaged selectivities, the activation of the off-target
            % fascicle averaged and then used to find (1-offtargetActivation)^k
            % For split selectivities, each off-target fascicle's activation is used to
            % find 0.5*(1-offtargetActivation)^k, and then both are subtracted from the
            % on-target activation
            channelSplitSelectivities_allNSims = channelSegNormCmaps_allNSims;
            channelAvgSelectivities_allNSims = channelSegNormCmaps_allNSims;
            for i = 1:3
                targetInd = i;
                offtargetInd = 1:3;
                offtargetInd = offtargetInd(offtargetInd~=targetInd);
                for j = 1:nsimNum
                    channelAvgSelectivities_allNSims{i}{j} = zeros(pulseParamTrends{i,4}{2},2);
                    channelSplitSelectivities_allNSims{i}{j} = zeros(pulseParamTrends{i,4}{2},2);
                    for k = 1:pulseParamTrends{i,4}{2}
                        if j ~= nSimsToExclude
                            targetActivation_allNSims = channelSegNormCmaps_allNSims{i}{j}{k};
                            offtarget1Activation_allNSims = channelSegNormCmaps_allNSims{offtargetInd(1)}{j}{k};
                            offtarget2Activation_allNSims = channelSegNormCmaps_allNSims{offtargetInd(2)}{j}{k};
                            offtargetAvgActivation_allNSims = mean([offtarget1Activation_allNSims,offtarget2Activation_allNSims]);
                            channelAvgSelectivities_allNSims{i}{j}(k,2) = targetActivation_allNSims + (1-offtargetAvgActivation_allNSims) - 1;
                            channelSplitSelectivities_allNSims{i}{j}(k,2) = targetActivation_allNSims + 0.5*(1-offtarget1Activation_allNSims) + 0.5*(1-offtarget2Activation_allNSims) - 1;
                        end
                    end
                end
            end
            
            stdForEachAmp = zeros(nsimNum,32);
            for j = 1:nsimNum
                for k = 1:pulseParamTrends{1,4}{2}
                    if j ~= nSimsToExclude
                        stdForEachAmp(j,k) = std([channelAvgSelectivities_allNSims{1}{j}(k,2),channelAvgSelectivities_allNSims{2}{j}(k,2),channelAvgSelectivities_allNSims{3}{j}(k,2)]);
                    end
                end
            end
            [maxStdAcrossAmps,ind_maxStd] = max(stdForEachAmp,[],2);
           
    
            %% Calculate highest selectivities per NSim and per channel, then save for all datasets
            maxSelectivityPerNSim = zeros(nsimNum,1);
            mostSelectiveChannel = zeros(nsimNum,1);
            mostSelectiveAmplitude = zeros(nsimNum,1);
            for i = 1:nsimNum
                [Channel1MaxSelectivity(i),amplitudeInd_Ch1MaxSel(i)] = max(channelAvgSelectivities_allNSims{1}{i}(:,2));
                [Channel2MaxSelectivity(i),amplitudeInd_Ch2MaxSel(i)] = max(channelAvgSelectivities_allNSims{2}{i}(:,2));
                [Channel3MaxSelectivity(i),amplitudeInd_Ch3MaxSel(i)] = max(channelAvgSelectivities_allNSims{3}{i}(:,2));
                [maxSelectivityPerNSim(i), mostSelectiveChannel(i)] = max([Channel1MaxSelectivity(i), Channel2MaxSelectivity(i), Channel3MaxSelectivity(i)]);
            end
            [Ch1Top10MaxSel, nSimInd_Ch1Top10MaxSel] = maxk(Channel1MaxSelectivity,10)
            all_Ch1Top10MaxSel{nerveIndex}{datasetIndex}=[Ch1Top10MaxSel; amplitudeInd_Ch1MaxSel(nSimInd_Ch1Top10MaxSel)+(ampNum-1)*nSimInd_Ch1Top10MaxSel].'
            [Ch2Top10MaxSel, nSimInd_Ch2Top10MaxSel] = maxk(Channel2MaxSelectivity,10)
            all_Ch2Top10MaxSel{nerveIndex}{datasetIndex}=[Ch2Top10MaxSel; amplitudeInd_Ch2MaxSel(nSimInd_Ch2Top10MaxSel)+(ampNum-1)*nSimInd_Ch2Top10MaxSel].'
            [Ch3Top10MaxSel, nSimInd_Ch3Top10MaxSel] = maxk(Channel3MaxSelectivity,10)
            all_Ch3Top10MaxSel{nerveIndex}{datasetIndex}=[Ch3Top10MaxSel; amplitudeInd_Ch3MaxSel(nSimInd_Ch3Top10MaxSel)+(ampNum-1)*nSimInd_Ch3Top10MaxSel].'
        
            all_channelAvgSelectivites{nerveIndex}{datasetIndex} = channelAvgSelectivities_allNSims;
    
            %% Now upload and generate selectivity values for the simulation data with the same nerve and dataset Indices
    
            % Sample, Model, and Sim directories already established for each
            % nerve and dataset
    
            % Navigate to the nsim root folder and get the nsim folder names
            nsimDir = fullfile(simDir, 'n_sims');
            cd(nsimDir)
            nsimFiles = dir;
    
            % Navigate to the fibreset folder, and get the files containing the fibre
            % positions
            positionDir = fullfile(simDir, 'fibersets\0');
            cd(positionDir)
            positionFiles = dir;
    
            % Navigate to the waveform folder and get the waveform files
            waveformDir = fullfile(simDir, 'waveforms');
            cd(waveformDir)
            waveformFiles = dir;
    
            % Reorganise each of the file sets into natural sorted order
            natsortfilesDir = fullfile(scriptDir, 'natsortfiles');
            cd(natsortfilesDir);
    
            nsimFiles = natsortfiles(nsimFiles);
            positionFiles = natsortfiles(positionFiles);
            waveformFiles = natsortfiles(waveformFiles);
    
            nsimNum = size(nsimFiles,1)-2;
            positionNum = size(positionFiles,1)-2;
            waveformNum = size(waveformFiles,1)-2;
    
            %%
    
            % Generate parametric arrays in order to allow for formulation of ellipse
            % equations and their plotting
    
            % Set the length of the spatial arrays
            plotPrecision = 1000;
    
            % Create the parameter t, the angle to the point on the ellipse in degrees
            t = linspace(0,360,plotPrecision);
    
            % Values are in um (or um^2) in the file
    
            % Find the area of the nerve. Sample doesn't seem to contain nerve shape
            % data, though it's possible this is only because my test nerves are all
            % circular. Would be good to test with an elliptical nerve
            nerveArea = sampleData.Morphology.Nerve.area;
            % Set a flag if only the nerve area is given, since I assume this means the
            % nerve is circular. Use this assumption to calculate the nerve radius
            circularNerve = false;
            nerveRadius = NaN;
            if numel(sampleData.Morphology.Nerve) == 1
                circularNerve = true;
                nerveRadius = sqrt(nerveArea/pi);
                nerveX = nerveRadius*cosd(t);
                nerveY = nerveRadius*sind(t);
            end
    
            % Get the number of fascicles in the sample
            fascicleNum = numel(sampleData.Morphology.Fascicles);
            % For all fascicles, get the relevant parameters for plotting. There are
            % both "inners" and "outers" in the data, which each have 5 parameter
            % values (plus area, but that is unimportant for plotting). So data
            % structure needed is fascicleNumx2x5 to contain all data
            fascicleParams = zeros(fascicleNum,2,5);
            % For all fascicles
            for i = 1:fascicleNum
                fascicleFields = fieldnames(sampleData.Morphology.Fascicles(i));
                % For both outer and inners
                for j = 1:2
                    paramFields = fieldnames(getfield(sampleData.Morphology.Fascicles(i),fascicleFields{j}));
                    % For the relevant parameters, which comprise fields 2 to 6 of
                    % paramFields
                    for k = 2:6
                        fascicleParams(i,j,k-1) = getfield(getfield(sampleData.Morphology.Fascicles(i),fascicleFields{j}),paramFields{k});
                    end
                end
            end
    
            %%
    
            % Create a model of the cuff from the model file
    
            % Create the path to the cuff directory. This goes up one directory
            % from sampleDir to the ascentDir and then goes to the cuff
            % directory in config
            [ascentDir, ~, ~] = fileparts(allSamplesDir);
            cuffDir = fullfile(ascentDir, 'config\system\cuffs');
    
            % From the model file, find the cuff used and the geometric parameters.
            cd(modelDir)
            modelText = fileread('model.json');
            modelData = jsondecode(modelText);
            cuff = modelData.cuff.preset;
            cuffShiftX = modelData.cuff.shift.x;
            cuffShiftY = modelData.cuff.shift.y;
            cuffRot = modelData.cuff.rotate.pos_ang + modelData.cuff.rotate.add_ang;
    
            % Load the cuff file using the path and name
            cd(cuffDir)
            cuffText = fileread(cuff);
            cuffData = jsondecode(cuffText);
    
            % By finding which instances in the cuff model contain the string
            % 'Contact', find the number of contacts
            contactNum = 0;
            for i = 1:numel(cuffData.instances)
                if strfind(cuffData.instances(i).type,'Contact') > 0
                    contactNum = contactNum + 1;
                end
            end
            contactInstances = zeros(contactNum,1);
            contInstPointer = 1;
            for i = 1:numel(cuffData.instances)
                if strfind(cuffData.instances(i).type,'Contact') > 0
                    contactInstances(contInstPointer) = i;
                    contInstPointer = contInstPointer + 1;
                end
            end
    
            % The following code only functions for rectangular contacts. Could be
            % expanded to other kinds in the future if the need arises
    
            % Find the number of rectangular contacts
            rectangularContactNum = 0;
            for i = 1:contactNum
                if cuffData.instances(contactInstances(i)).type == "RectangleContact_Primitive"
                    rectangularContactNum = rectangularContactNum + 1;
                end
            end
            contactFlag = false;
    
            % If all the contacts are rectangular
            if rectangularContactNum == contactNum
                % Set the flag for plotting
                contactFlag = true;
    
                % Get the number of parameters in the cuff file
                cuffParamNum = numel(cuffData.params);
                % Create a cell lookup table for the parameters.
                % It needs 4 columns to account for units being inconsistently added to
                % the values
                cuffParams = cell(cuffParamNum,4);
                % Populate cuffParams
                for i = 1:cuffParamNum
                    cuffParams{i,1} = cuffData.params(i).name;
                    expression = split(cuffData.params(i).expression);
                    for j = 1:numel(expression)
                        cuffParams{i,j+1} = expression{j};
                    end
                    cuffParams{i,4} = cuffData.params(i).description;
                end
    
                % For each contact, 5 variables completely define the contact in 2D:
                % cuff radius, contact angle, contact width, recess thickness and
                % contact thickness.
                % For the custom cuffs based on Pitt (PittMod series), all of these
                % except for angle are common across all contacts so can be directly
                % read from cuffParams.
    
                % Get the contact angles (ACW angle from x axis to centre of cuff)
                contactRots = zeros(contactNum,1);
                for i = 1:contactNum
                    contactRotString = split(cuffData.instances(contactInstances(i)).def.Rotation_angle);
                    contactRotString = contactRotString{1};
                    for j = 1:cuffParamNum
                        if strfind(cuffParams{j,1},contactRotString) > 0
                            contactRots(i) = str2double(cuffParams{j,2});
                        end
                    end
                end
    
                % Remove any duplicates from contactRots, as these represent
                % vertically stacked contacts that can't be shown in a 2D plot. Change
                % the number of contacts accordingly
                contactRots = unique(contactRots);
                contactNum = length(contactRots);
    
                % Get the other contact parameters from cuffParams
                % Need to check the unit here, since if it's mm then need to multiply
                % by 1000 to match other um values
                for i = 1:cuffParamNum
                    if strfind(cuffParams{i,1},cuffData.instances(contactInstances(1)).def.R_in) > 0
                        if strfind(cuffParams{i,3},'mm') > 0
                            cuffRadius = str2double(cuffParams{i,2})*1000;
                        elseif strfind(cuffParams{i,3},'um') > 0
                            cuffRadius = str2double(cuffParams{i,2});
                        end
                    elseif strfind(cuffParams{i,1},cuffData.instances(contactInstances(1)).def.Rect_w) > 0
                        if strfind(cuffParams{i,3},'mm') > 0
                            contactWidth = str2double(cuffParams{i,2})*1000;
                        elseif strfind(cuffParams{i,3},'um') > 0
                            contactWidth = str2double(cuffParams{i,2});
                        end
                    elseif strfind(cuffParams{i,1},cuffData.instances(contactInstances(1)).def.Rect_thk) > 0
                        if strfind(cuffParams{i,3},'mm') > 0
                            contactThk = str2double(cuffParams{i,2})*1000;
                        elseif strfind(cuffParams{i,3},'um') > 0
                            contactThk = str2double(cuffParams{i,2});
                        end
                    elseif strfind(cuffParams{i,1},cuffData.instances(contactInstances(1)).def.Rect_recess) > 0
                        if strfind(cuffParams{i,3},'mm') > 0
                            contactRecessThk = str2double(cuffParams{i,2})*1000;
                        elseif strfind(cuffParams{i,3},'um') > 0
                            contactRecessThk = str2double(cuffParams{i,2});
                        end
                    end
                end
    
                % Create plotting arrays for the cuff, assuming it to be a circle with
                % internal radius R_in
                cuffX = cuffRadius*cosd(t)+cuffShiftX;
                cuffY = cuffRadius*sind(t)+cuffShiftY;
    
                % From the contact and cuff parameters, find the three half angles
                % swept by the base contact (centred on positive x-axis), in degrees
                contactTheta1 = atand(contactWidth/(2*cuffRadius));
                contactTheta2 = atand(contactWidth/(2*(cuffRadius+contactRecessThk)));
                contactTheta3 = atand(contactWidth/(2*(cuffRadius+contactRecessThk+contactThk)));
    
                % Create parameters to represent the angle to the point on the base
                % contact, for each of the three layers
                t2 = linspace(-contactTheta1,contactTheta1,plotPrecision);
                t3 = linspace(-contactTheta2,contactTheta2,plotPrecision);
                t4 = linspace(-contactTheta3,contactTheta3,plotPrecision);
    
            end
    
            %%
    
            % Extract the waveform data for each waveform and load it all into a cell
            % array
            cd(waveformDir)
            waveformEndIndices = zeros(waveformNum,1);
            for i = 1:waveformNum
                waveData = importdata(waveformFiles(i+1).name);
                % Take off the first two values, since they correspond to information
                % about the wave rather than wave data
                waveData = waveData(3:end);
                waveforms{cast(str2double(waveformFiles(i+1).name(1:end-4)),'int32')+1} = waveData;
                tempArray = find(waveData);
                % Find the final non-zero value in the wave data
                waveformEndIndices(i) = tempArray(end);
            end
    
            % Edit the waveforms so that they are all the same length, but excess zero
            % values from the end of the array
    
            % Find the maximum waveform end index
            waveformMaxEndIndex = max(waveformEndIndices);
            % Add zeroes onto the end of each wave array so that it is an integer
            % number of 200 us chunks (with a 200 us buffer if the wave is already an
            % integer number).
            paddingEndIndex = 200-mod(waveformMaxEndIndex,200)+waveformMaxEndIndex;
            for i = 1:waveformNum
                tempArray = find(waveforms{i});
                waveforms{i} = [waveforms{i}(1:tempArray(end));zeros(paddingEndIndex-tempArray(end),1)];
            end
            % Create the waveTime array for plotting
            waveTime = -0.1:0.001:paddingEndIndex*0.001;
            waveTime = reshape(waveTime,numel(waveTime),1);
            % Create a padding array to clarify the wave plots
            frontPadding = zeros(101,1);
    
            %%
    
            % Extract the fibre position and diameter data and load into an array
    
            % Move to the position data directory
            cd(positionDir)
            % Get the number of fibres in the sim
            fibreNum = size(positionFiles,1)-4;
            % Initialise a matrix to store fibre data
            fibreData = zeros(fibreNum,5);
            % For all fibres, extract the x-y positional data
            for i = 1:fibreNum
                fid = fopen(positionFiles(i+1).name);
                positionDatastring = textscan(fid, '%s %s %s', 2, 'Delimiter', ',');
                positionDatastring = positionDatastring{1};
                positionDatatemp = strsplit(positionDatastring{2},' ');
                for j = 1:2
                    fibreData(i,j) = str2double(positionDatatemp{j});
                end
                fclose(fid);
            end
    
            % For all fibres, extract the diameter in um from diams.txt
            diamData = importdata('diams.txt');
            for i = 1:fibreNum
                fibreData(i,3) = diamData(i);
            end
    
            %%
    
            % Create models of the fascicle inners, and use them to assign fibres to
            % fascicles
    
            fascicleBounds = cell(fascicleNum,1);
            % Get the area of each fascicle from the sample data for reordering due
            % to Grill lab. Want to order by area from largest to smallest,
            % assuming this maps to tibial-peroneal-sural
            fascicleAreas = zeros(fascicleNum,1);
            for i = 1:fascicleNum
                fascicleAreas(i) = sampleData.Morphology.Fascicles(i).outer.area;
            end
            % Sort the areas in order of size, descending
            [~,fascicleSortIndex] = sort(fascicleAreas,'descend');
            sortedFascicleParams = zeros(size(fascicleParams));
            for i = 1:fascicleNum
                sortedFascicleParams(i,:,:) = fascicleParams(fascicleSortIndex(i),:,:);
            end
            fascicleParams = sortedFascicleParams;
            for i = 1:fascicleNum
                % Create the spatial arrays for each fascicle, centred on (0,0)
                % to allow for rotation. The semi-major radius a is half of the
                % major axis length stored in fascicleParams(x,y,3); same for
                % the semi-minor radius b [(x,y,4)]
                fascicleMajorRadius = fascicleParams(i,2,3)/2;
                fascicleMinorRadius = fascicleParams(i,2,4)/2;
                fascicleRot = fascicleParams(i,2,5);
    
                fascicleX = fascicleMajorRadius*cosd(t);
                fascicleY = fascicleMinorRadius*sind(t);
    
                % Create the fascicle rotation matrix and rotate the fascicle
                Rfascicle = [cosd(fascicleRot) sind(fascicleRot); -sind(fascicleRot) cosd(fascicleRot)];
    
                fascicleCoords = [fascicleX;fascicleY]'*Rfascicle;
                % Shift the fascicle centre, and break up into X and Y arrays
                % for easier plotting
                fascicleX = fascicleCoords(:,1)+fascicleParams(i,2,1);
                fascicleY = fascicleCoords(:,2)+fascicleParams(i,2,2);
    
                % Create a boundary and store it
                fascicleBounds{i} = polyshape(fascicleX,fascicleY);
            end
    
            % Assign fibres to fascicles based on the generated inner boundaries
            for i = 1:fibreNum
                for j = 1:fascicleNum
                    if isinterior(fascicleBounds{j},fibreData(i,1:2))
                        fibreData(i,4) = j;
                    end
                end
            end
    
            %%
    
            % For all of the nsims, extract the nsim parameter and threshold data and
            % load into a cell array
            nsimParams = cell(nsimNum,1);
            nsimThresholds = cell(nsimNum,1);
            nsimContactWeights = cell(nsimNum,1);
            nsimWaveIndices = cell(nsimNum,1);
    
            for i = 1:nsimNum
                nsimParamDir = fullfile(nsimDir, num2str(i-1));
                nsimOutputDir = fullfile(nsimParamDir, 'data\outputs');
                cd(nsimParamDir)
                nsimParamText = fileread([num2str(i-1) '.json']);
                nsimParams{i} = jsondecode(nsimParamText);
                cd(nsimOutputDir)
                nsimOutputFiles = dir;
                cd(natsortfilesDir)
                nsimOutputFiles = natsortfiles(nsimOutputFiles);
                cd(nsimOutputDir)
                tempArray = zeros(fibreNum,1);
                for j = 1:fibreNum
                    tempArray(j) = importdata(nsimOutputFiles(j+2).name);
                end
                nsimThresholds{i} = tempArray;
                activeSourceFields = fields(nsimParams{i}.active_srcs);
                nsimContactWeights{i} = getfield(nsimParams{i}.active_srcs,activeSourceFields{1});
                nsimWaveIndices{i} = nsimParams{i}.waveform.EXPLICIT.index-nsimParams{1}.waveform.EXPLICIT.index+1;
            end
    
            %%
    
            % Find the unique fibre diameter values and assign fibre indices to those
            % values
            uniqueFibreDiams = unique(fibreData(:,3));
            uniqueFibreDiamNum = size(uniqueFibreDiams,1);
            for i = 1:uniqueFibreDiamNum
                fibreDiamMatchIndices = find(fibreData(:,3)==uniqueFibreDiams(i));
                fibreData(fibreDiamMatchIndices,5) = i;
            end
    
            % Create a diameter-threshold map for each nsim
    
            nsimUniqueFibreDiamThresholds = cell(nsimNum,1);
            for i = 1:nsimNum
                tempArray = zeros(uniqueFibreDiamNum,1);
                for j = 1:uniqueFibreDiamNum
                    tempArray(j) = mean(nsimThresholds{i}(fibreData(:,5)==j));
                end
                nsimUniqueFibreDiamThresholds{i} = tempArray;
            end
    
            %%
    
            % Calculate the per-fascicle normalised activation curves for each nsim
            nsimActivations = cell(nsimNum,fascicleNum);
            nsimNormActivations = cell(nsimNum,fascicleNum);
            for i = 1:nsimNum
                for j = 1:fascicleNum
                    peakCurrent = 1.6;
                    stepCurrent = 0.05;
                    fascicleThresholds = nsimThresholds{i}(fibreData(:,4)==j);
                    maxThreshold = max(nsimThresholds{i});
                    if maxThreshold < 0
                        thresholdTestRange = flip(-peakCurrent:stepCurrent:0);
                    else
                        thresholdTestRange = 0:stepCurrent:peakCurrent;
                    end
                    nsimActivations{i,j} = false(length(find(fibreData(:,4)==j)),length(thresholdTestRange));
                    nsimNormActivations{i,j} = zeros(1,length(thresholdTestRange));
                    for k = 1:length(thresholdTestRange)
                        if i ~= nSimsToExclude
                            if maxThreshold < 0
                                nsimActivations{i,j}(:,k) = fascicleThresholds >= thresholdTestRange(k);
                            else
                                nsimActivations{i,j}(:,k) = fascicleThresholds <= thresholdTestRange(k);
                            end
                            nsimNormActivations{i,j}(k) = length(find(nsimActivations{i,j}(:,k)))/length(nsimActivations{i,j}(:,k));
                        end 
                    end
                end
            end
    
            %%
    
            % Using the normalised activations and 3rd power off target activation,
            % calculate the nsim per-fascicle selectivities
            nsimFascicleSelectivities = zeros(nsimNum,fascicleNum,length(thresholdTestRange));
            for i = 1:nsimNum
                for j = 1:fascicleNum
                    if i ~= nSimsToExclude
                        targetInd = j;
                        offTargetInds = find(1:fascicleNum~=j);%%
                        targetNormActivation = nsimNormActivations{i,targetInd};
                        offTarget1NormActivation = nsimNormActivations{i,offTargetInds(1)};
                        offTarget2NormActivation = nsimNormActivations{i,offTargetInds(2)};
                        nsimFascicleSelectivities(i,j,:) = targetNormActivation - mean([offTarget1NormActivation;offTarget2NormActivation]);
                    end
                end
            end
    
            [amplitude_MaxSIs_sim, ind_amp_MaxSIs_sim] = max(nsimFascicleSelectivities,[],3)
            [waveform_elec_MaxSIs_sim, ind_waveform_elec_MaxSIs_sim] = max(amplitude_MaxSIs_sim,[],1)
            waveform_MaxSIs_sim = ceil(ind_waveform_elec_MaxSIs_sim./contactNum)
            elec_MaxSIs_sim = ceil(ind_waveform_elec_MaxSIs_sim./(length(amplitude_MaxSIs_sim)/8))
            all_top20maxSels{datasetIndex,1} = maxk(max(nsimFascicleSelectivities,[],3),20,1)
            all_top20maxSels{datasetIndex,2} = elec_MaxSIs_sim
            all_top20maxSels{datasetIndex,3} = waveform_MaxSIs_sim
            all_top20maxSels{datasetIndex,4} = ind_amp_MaxSIs_sim(ind_waveform_elec_MaxSIs_sim)*50
    
            %% Calculating and Appending Data to Dataframe for Regression Analysis
            
            % Correct for bipolar vs tripolar amplitudes
            switch simNum
               case ['300']
                   source_distance_array = allTrialsElectrodeFascicleRels;
                   sink_distance_array = allTrialsElectrodeFascicleRels;
                   amp_multiplier = 1;
               case ['303']
                   source_distance_array = allTrialsElectrodeFascicleRels;
                   sink_distance_array = allTrialsElectrodeFascicleRels;
                   amp_multiplier = 1;
               case ['305']
                   source_distance_array = allTrialsElectrodeFascicleRels;
                   sink_distance_array = allTrialsElectrodeTripolarFascicleRels;
                   amp_multiplier = 2;
               case ['306']
                   source_distance_array = allTrialsElectrodeTripolarFascicleRels;
                   sink_distance_array = allTrialsElectrodeFascicleRels;
                   amp_multiplier = 2;
               case ['307']
                   source_distance_array = allTrialsElectrodeFascicleRels;
                   sink_distance_array = allTrialsElectrodeFascicleRels;
                   amp_multiplier = 1;
            end
        
            
            for j = 1:nsimNum
                for k = 1:pulseParamTrends{1,4}{2}
                    for i = 1:chanNum
                        % Extract normalized cmap values
                        normActivation_Tibial = channelSegNormCmaps_allNSims{1}{j}{k};
                        normActivation_Peroneal = channelSegNormCmaps_allNSims{2}{j}{k};
                        normActivation_Sural = channelSegNormCmaps_allNSims{3}{j}{k};
                
                        % Check if any cmap value is missing
                        if isempty(normActivation_Tibial) || isempty(normActivation_Peroneal) || isempty(normActivation_Sural)
                            skippedRowsCount(nerveIndex, datasetIndex) = skippedRowsCount(nerveIndex, datasetIndex) + 1;
                            continue; % Skip this iteration
                        end
        
                        newRow = table( ...
                            string(nerveIndex), string(fascicleNames{i}), string(j), string(simNum), ...
                            amp_multiplier*reorderedStimCycle{k}, amp_multiplier*reorderedStimCycle{k} * pulseParamTrends{1,1}{1,1}(mod(j,waveformNum)) / pulseParamTrends{1,3}{1,1}(mod(j,waveformNum)), ...
                            pulseParamTrends{1,1}{1,1}(mod(j,waveformNum)), pulseParamTrends{1,3}{1,1}(mod(j,waveformNum)), ...
                            reorderedStimCycle{k}*pulseParamTrends{1,1}{1,1}(mod(j,waveformNum)), ...
                            string(active_srcs_cat(ceil(j/waveformNum))), string(active_sinks_cat(ceil(j/waveformNum))), ...
                            source_distance_array{nerveIndex}{active_srcs_index(ceil(j/waveformNum)),i}, ...
                            sink_distance_array{nerveIndex}{active_sinks_index(ceil(j/waveformNum)),i}, ...
                            string(Number_Poles), string(Phase_Symmetry), ...
                            channelSegNormCmaps_allNSims{i}{j}{k}, nsimNormActivations{j,i}(k+1), ...
                            channelAvgSelectivities_allNSims{i}{j}(k,2),nsimFascicleSelectivities(j,i,k+1), ...
                            'VariableNames', regressionData.Properties.VariableNames);
        
                        % Append new row to regressionData
                        regressionData = [regressionData; newRow];
        
                    end
                end
            end
            percentSkippedPerConfig(nerveIndex, datasetIndex) = skippedRowsCount(nerveIndex,datasetIndex)/(3*nsimNum*pulseParamTrends{1,4}{2});
    
            %% Plotting Heatmaps for all 8 electrodes with all waveforms and amps - Ex Vivo
            savefolder = "C:\Users\Green_Group\Zack\ExVivo\Results\Results_Heatmaps_ExVivo_Sim";
            configFolder = 'C:\Users\Green_Group\Zack\ExVivo\NerveOrientations\Configs';
            nerveCuffFolder = 'C:\Users\Green_Group\Zack\ExVivo\NerveOrientations\CuffOnNerveExVivo';
            
            fig = figure(figNum);
            figNum = figNum + 1;
            set(gcf,'Color',[1 1 1]);
            try 
                fig.WindowState = 'maximized';
            catch
                set(fig,'Units','normalized','OuterPosition',[0 0 1 1]);
            end

            % create a 3x3 tiled layout with minimal gaps
            t = tiledlayout(3,3, ...
                 'TileSpacing','compact', ...
                 'Padding','compact');
            
            circPlotMap = [2,1,4,7,8,9,6,3];
            plotAmplitudes = linspace(0,ampNum-1,1+(ampNum-1)/4)+1;
            plotAmplitudeLabels = cell(1,length(plotAmplitudes));
            for k = 1:length(plotAmplitudes)
                plotAmplitudeLabels{k} = num2str(abs(thresholdTestRange(plotAmplitudes(k))));
            end
            
            
            % ---- Heatmap plotting settings: show 6 PW columns (50–300 us) with tick labels centered on pixels ----
            plotPulseWidths = [50 100 150 200 250 300];
            plotWaveformNum = numel(plotPulseWidths);

            % Map desired pulse-width labels to the closest waveform indices in the data.
            pwAll = pulseParamTrends{1,3}{1,1};
            if max(pwAll) < 1
                pwAll_us = pwAll * 1e6;   % seconds -> microseconds
            else
                pwAll_us = pwAll;         % already in microseconds
            end
            plotPWindices = zeros(1, plotWaveformNum);
            for pwIdx = 1:plotWaveformNum
                [~, w] = min(abs(pwAll_us - plotPulseWidths(pwIdx)));
                plotPWindices(pwIdx) = w;
            end


            % draw each of the 8 heatmaps
            simID = str2double(datasetNames{1,datasetIndex}{2});
            isBipolar = ismember(simID, [300 303 307]);
            isTripolar = ismember(simID, [305 306]);

            for elec = 1:channelCount
                ax = nexttile(circPlotMap(elec));
                hold(ax,'on');

                % Only plot the first 0–300 us pulse-width set, even if waveformNum is larger
                baseIdx = (elec-1)*waveformNum;
                for pwIdx = 1:plotWaveformNum
                    w = plotPWindices(pwIdx);           % waveform column in full data
                    i = baseIdx + w;
                                            for j = 1:ampNum-1
                        ch1Act = channelSegNormCmaps_allNSims{1}{i}{j};
                        ch2Act = channelSegNormCmaps_allNSims{2}{i}{j};
                        ch3Act = channelSegNormCmaps_allNSims{3}{i}{j};
                        if isempty(ch1Act) || isempty(ch2Act) || isempty(ch3Act)
                            continue;
                        end
                            rectangle('Position',[pwIdx-1, j, 1, 1], ...
                                      'EdgeColor',[0 0 0], ...
                                      'FaceColor',[ch1Act,ch2Act,ch3Act], ...
                                      'Parent',ax);
                        end
                end

                % ---- Electrode naming (wrap-around cuff numbering) ----
                if isBipolar
                    e1 = elec;
                    e2 = elec + 1; if e2 == 9, e2 = 1; end
                    title(ax, sprintf('Electrodes %d - %d', e1, e2), 'FontSize',16,'FontName','Times New Roman');
                elseif isTripolar
                    e1 = elec - 1; if e1 == 0, e1 = 8; end
                    e2 = elec;
                    e3 = elec + 1; if e3 == 9, e3 = 1; end
                    title(ax, sprintf('Electrodes %d - %d - %d', e1, e2, e3), 'FontSize',16,'FontName','Times New Roman');
                end

                % ---- Axes: FIXED x-range so config 300 doesn't "shrink" the heatmaps when nsimNum differs ----
                xlim(ax, [0, plotWaveformNum]);
                ylim(ax, [1, ampNum]);
                yticks(ax, plotAmplitudes);
                yticklabels(ax, plotAmplitudeLabels);

                xline(ax, 1:(plotWaveformNum-1), 'Color',[0.5 0.5 0.5]);
                yline(ax, 1:ampNum, 'Color',[0.5 0.5 0.5]);
                xline(ax, [0, plotWaveformNum], 'k');
                yline(ax, [1, ampNum], 'k');

                xticks(ax, (0:(plotWaveformNum-1)) + 0.5);
                xticklabels(ax, plotPulseWidths);
                xlabel(ax, 'Pulse Width (us)', 'FontSize',16,'FontName','Times New Roman');
                ylabel(ax, 'Current (mA)',     'FontSize',16,'FontName','Times New Roman');
            end
            % grab & delete center tile
            axCenter = nexttile(5);
            pos = axCenter.Position;    % [x y width height]
            delete(axCenter);
            
            % scale factors (0→1): how much of the tile to occupy
            hScale = 1.75;  % 90% of tile width total (45% each side)
            vScale = 1.75;  % 80% of tile height
            
            % computed dimensions
            tileW = pos(3);
            tileH = pos(4);
            imgW  = (tileW * hScale) / 2;
            imgH  = tileH * vScale;
            xOff  = (tileW * (1-hScale))/2;
            yOff  = (tileH * (1-vScale))/2;
            
            configFile = sprintf('Config%d.png',str2num(datasetNames{1,datasetIndex}{2}))
            configFileName = fullfile(configFolder,configFile)
            % left image axes
            leftPos  = [ pos(1) + 1.4*xOff, ...
                         pos(2) + yOff, ...
                         imgW, ...
                         imgH ];
            axL = axes('Position', leftPos);
            imshow(imread(configFileName), 'Parent', axL);
            axis(axL,'off');
            
            nerveFile = sprintf('Trial%d.png',all_Nerves(nerveIndex))
            nerveFileName = fullfile(nerveCuffFolder,nerveFile)
            % right image axes
            rightPos = [ pos(1) + 1.7*xOff + imgW, ...
                         pos(2) + yOff, ...
                         imgW, ...
                         imgH ];
            axR = axes('Position', rightPos);
            imshow(imread(nerveFileName), 'Parent', axR);
            axis(axR,'off');

            
            % save and close
            fname = sprintf('Heatmap_ExVivo_Nerve%d_Config%d.png', ...
                            all_Nerves(nerveIndex), str2num(datasetNames{1,datasetIndex}{2}));
            exportgraphics(fig, fullfile(savefolder,fname), ...
                           'BackgroundColor','white','Resolution',300);
            close(fig);
            clear axR
            clear axL
            clear rightPos
            clear leftPos
            clear t
            clear tile_W
            clear tile_H
            
            %% Plotting Heatmaps for all 8 electrodes with all waveforms and amps - Simulation
            savefolder = "C:\Users\Green_Group\Zack\ExVivo\Results\Results_Heatmaps_ExVivo_Sim";
            configFolder = 'C:\Users\Green_Group\Zack\ExVivo\NerveOrientations\Configs';
            nerveCuffFolder = "C:\Users\Green_Group\Zack\ExVivo\NerveOrientations\CuffOnNerveSim";
            fig = figure(figNum);
            figNum = figNum + 1;
            set(gcf,'Color',[1 1 1]);
            try 
                fig.WindowState = 'maximized';
            catch
                set(fig,'Units','normalized','OuterPosition',[0 0 1 1]);
            end

            % create a 3x3 tiled layout with minimal gaps
            t = tiledlayout(3,3, ...
                 'TileSpacing','compact', ...
                 'Padding','compact');
            
            circPlotMap = [2,1,4,7,8,9,6,3];
            plotAmplitudes = linspace(0,ampNum-1,1+(ampNum-1)/4)+1;
            plotAmplitudeLabels = cell(1,length(plotAmplitudes));
            for k = 1:length(plotAmplitudes)
                plotAmplitudeLabels{k} = num2str(abs(thresholdTestRange(plotAmplitudes(k))));
            end
            
            
            % ---- Heatmap plotting settings: show 6 PW columns (50–300 us) with tick labels centered on pixels ----
            plotPulseWidths = [50 100 150 200 250 300];
            plotWaveformNum = numel(plotPulseWidths);

            % Map desired pulse-width labels to the closest waveform indices in the data.
            % IMPORTANT: This avoids confusing CONFIG=300 with PW=300, and handles cases where the dataset
            % does/does not include a 0-us column or stores PWs in seconds.
            pwAll = pulseParamTrends{1,3}{1,1};
            if max(pwAll) < 1
                pwAll_us = pwAll * 1e6;   % seconds -> microseconds
            else
                pwAll_us = pwAll;         % already in microseconds
            end
            plotPWindices = zeros(1, plotWaveformNum);
            for pwIdx = 1:plotWaveformNum
                [~, w] = min(abs(pwAll_us - plotPulseWidths(pwIdx)));
                plotPWindices(pwIdx) = w;
            end


            % draw each of the 8 heatmaps
            simID = str2double(datasetNames{1,datasetIndex}{2});
            isBipolar = ismember(simID, [300 303 307]);
            isTripolar = ismember(simID, [305 306]);

            for elec = 1:channelCount
                ax = nexttile(circPlotMap(elec));
                hold(ax,'on');

                % Only plot the first 0–300 us pulse-width set, even if waveformNum is larger
                baseIdx = (elec-1)*waveformNum;
                for pwIdx = 1:plotWaveformNum
                    w = plotPWindices(pwIdx);           % waveform column in full data
                    i = baseIdx + w;
                                            for j = 1:ampNum-1
                        ch1Act = nsimNormActivations{i,1}(j+1);
                        ch2Act = nsimNormActivations{i,2}(j+1);
                        ch3Act = nsimNormActivations{i,3}(j+1);
                        if isempty(ch1Act) || isempty(ch2Act) || isempty(ch3Act)
                            continue;
                        end
                            rectangle('Position',[pwIdx-1, j, 1, 1], ...
                                      'EdgeColor',[0 0 0], ...
                                      'FaceColor',[ch1Act,ch2Act,ch3Act], ...
                                      'Parent',ax);
                        end
                end

                % ---- Electrode naming (wrap-around cuff numbering) ----
                if isBipolar
                    e1 = elec;
                    e2 = elec + 1; if e2 == 9, e2 = 1; end
                    title(ax, sprintf('Electrodes %d - %d', e1, e2), 'FontSize',16,'FontName','Times New Roman');
                elseif isTripolar
                    e1 = elec - 1; if e1 == 0, e1 = 8; end
                    e2 = elec;
                    e3 = elec + 1; if e3 == 9, e3 = 1; end
                    title(ax, sprintf('Electrodes %d - %d - %d', e1, e2, e3), 'FontSize',16,'FontName','Times New Roman');
                end

                % ---- Axes: FIXED x-range so config 300 doesn't "shrink" the heatmaps when nsimNum differs ----
                xlim(ax, [0, plotWaveformNum]);
                ylim(ax, [1, ampNum]);
                yticks(ax, plotAmplitudes);
                yticklabels(ax, plotAmplitudeLabels);

                xline(ax, 1:(plotWaveformNum-1), 'Color',[0.5 0.5 0.5]);
                yline(ax, 1:ampNum, 'Color',[0.5 0.5 0.5]);
                xline(ax, [0, plotWaveformNum], 'k');
                yline(ax, [1, ampNum], 'k');

                xticks(ax, (0:(plotWaveformNum-1)) + 0.5);
                xticklabels(ax, plotPulseWidths);
                xlabel(ax, 'Pulse Width (us)', 'FontSize',16,'FontName','Times New Roman');
                ylabel(ax, 'Current (mA)',     'FontSize',16,'FontName','Times New Roman');
            end

            % grab & delete center tile
            axCenter = nexttile(5);
            pos = axCenter.Position;    % [x y width height]
            delete(axCenter);
            
            % scale factors (0→1): how much of the tile to occupy
            hScale = 1.75;  % 90% of tile width total (45% each side)
            vScale = 1.75;  % 80% of tile height
            
            % computed dimensions
            tileW = pos(3);
            tileH = pos(4);
            imgW  = (tileW * hScale) / 2;
            imgH  = tileH * vScale;
            xOff  = (tileW * (1-hScale))/2;
            yOff  = (tileH * (1-vScale))/2;
            
            configFile = sprintf('Config%d.png',str2num(datasetNames{1,datasetIndex}{2}))
            configFileName = fullfile(configFolder,configFile)
            % left image axes
            leftPos  = [ pos(1) + 1.4*xOff, ...
                         pos(2) + yOff, ...
                         imgW, ...
                         imgH ];
            axL = axes('Position', leftPos);
            imshow(imread(configFileName), 'Parent', axL);
            axis(axL,'off');
            
            nerveFile = sprintf('Trial%d_Sim.png',all_Nerves(nerveIndex))
            nerveFileName = fullfile(nerveCuffFolder,nerveFile)
            % right image axes
            rightPos = [ pos(1) + 1.7*xOff + imgW, ...
                         pos(2) + yOff, ...
                         imgW, ...
                         imgH ];
            axR = axes('Position', rightPos);
            imshow(imread(nerveFileName), 'Parent', axR);
            axis(axR,'off');

            fname = sprintf('Heatmap_Simulation_Nerve%d_Config%d.png',all_Nerves(nerveIndex),str2num(datasetNames{1,datasetIndex}{2}));
            exportgraphics(fig,fullfile(savefolder,fname), 'BackgroundColor','white','Resolution',300);
            close(fig);
            clear axR
            clear axL
            clear rightPosnn
            clear leftPos
            clear t
            clear tile_W
            clear tile_H
            end
        end
    end
end

            %% Common part for all types of interactive plots

            close all
            %amplitudes ordered from -1600 to 0       Used to replace pulseParamTrends{1,4}{2}
            amps_increasing_order = [0; sort(-pulseParamTrends{1,4}{1})]
            amp_size = 32
            figForAnalysis = figure(figNum)
            hold on
            title(['Sim ',simNum,' Channels Combined - Selectivity Heatmap - All Configurations'])
            xlabel('Configuration No.')
            ylabel('Injection Current Amplitude [mA]')
            % yticks([0.5 pulseParamTrends{1,4}{2}/2-0.5 pulseParamTrends{1,4}{2}-0.5])
            yticks([0.5 amp_size/2-0.5 amp_size-0.5])
            % yticklabels({num2str(angleNormCnaps{1,1}),num2str(angleNormCnaps{1,1}(pulseParamTrends{1,4}{2}/2,1)),num2str(angleNormCnaps{1,1}(end,1))})
            %yticklabels({num2str(angleNormCnaps{1,1}),num2str(angleNormCnaps{1,1}(amp_size/2,1)),num2str(angleNormCnaps{1,1}(end,1))})

            %% Interactive plot showing activations mapped on top of each other

            for i = 1:nsimNum % over all nsims
                faceColour = [0 0 0];
                for j = 1:amp_size % over all amplitudes
                    if i ~= nSimsToExclude
                        ch1Act = channelSegNormCmaps_allNSims{1}{i}{j}; % norm cnap for 1 sim
                        ch2Act = channelSegNormCmaps_allNSims{2}{i}{j};
                        ch3Act = channelSegNormCmaps_allNSims{3}{i}{j};


                        currentSegment = j + amp_size*(i-1);
                        % if ismember(currentSegment,all_Ch1Top10MaxSel{datasetIndex,1}(:,2)) | ismember(currentSegment,all_Ch2Top10MaxSel{datasetIndex,1}(:,2)) | ismember(currentSegment,all_Ch3Top10MaxSel{datasetIndex,1}(:,2))
                        %     edgeColor = 'white';
                        if ismember(currentSegment,all_Ch1Top10MaxSel{nerveIndex}{1,datasetIndex}(:,2)) 
                            edgeColor = 'red';
                        elseif ismember(currentSegment,all_Ch2Top10MaxSel{nerveIndex}{1,datasetIndex}(:,2))
                            edgeColor = 'green';
                        elseif ismember(currentSegment,all_Ch3Top10MaxSel{nerveIndex}{1,datasetIndex}(:,2))
                            edgeColor = 'blue';
                        elseif ismember(currentSegment,top_max_Cnaps{datasetIndex,1}{2}(:)) | ismember(currentSegment,top_max_Cnaps{datasetIndex,2}{2}(:)) | ismember(currentSegment,top_max_Cnaps{datasetIndex,3}{2}(:))
                            edgeColor = 'black';
                        else
                            edgeColor = 'none';
                        end

                        rectangle('Position',[i-1 j-1 1 1],'EdgeColor',edgeColor,'FaceColor',[ch1Act, ch2Act, ch3Act]);
                    end
                end
            end

            hold off

            figForAnalysis = figure(figNum) 
            stimOrder = amps_increasing_order;
            k = figForAnalysis; %Some Figure
            while true
                figure(figNum)
                try
                    [x, y] = ginput(1);
                    if ~ishandle(figure(figNum))
                        disp('Figure closed.');
                        break; % Break the loop if the figure is closed
                    end
                    ampIndForFig = floor(y)+1;
                    configIndForFig = floor(x)+1;

                    % Define text and variables
                    Sel1 = channelAvgSelectivities_allNSims{1,1}{configIndForFig}(ampIndForFig,2);
                    Sel2 = channelAvgSelectivities_allNSims{2,1}{configIndForFig}(ampIndForFig,2);
                    Sel3 = channelAvgSelectivities_allNSims{3,1}{configIndForFig}(ampIndForFig,2);
                    Act1 = channelSegNormCmaps_allNSims{1}{configIndForFig}{ampIndForFig};
                    Act2 = channelSegNormCmaps_allNSims{2}{configIndForFig}{ampIndForFig};
                    Act3 = channelSegNormCmaps_allNSims{3}{configIndForFig}{ampIndForFig};

                    %figure()
                    figure('Position',[5 200 550 350])  %To have figure pop on the left side of the screen instead of the middle, easier for visualization
                    plot(reorderedSmoothedChannelData{1,2}{configIndForFig}{ampIndForFig},'color',deep_red, 'DisplayName','Tibial','LineWidth',1.5)
                    hold on
                    plot(reorderedSmoothedChannelData{2,2}{configIndForFig}{ampIndForFig},'color',deep_green, 'DisplayName','Peroneal','LineWidth',1.5)
                    plot(reorderedSmoothedChannelData{3,2}{configIndForFig}{ampIndForFig},'color',deep_blue, 'DisplayName','Sural','LineWidth',1.5)
                    title('Fascicle CNAPs','FontName','Times New Roman')
                    line([1 25], [channelSegmentCmapsPosPeak{1}{configIndForFig}{ampIndForFig} channelSegmentCmapsPosPeak{1}{configIndForFig}{ampIndForFig}], 'Color', deep_red,'HandleVisibility','off','LineStyle','--','LineWidth',0.5);
                    line([1 25], [channelSegmentCmapsNegPeak{1}{configIndForFig}{ampIndForFig} channelSegmentCmapsNegPeak{1}{configIndForFig}{ampIndForFig}], 'Color', deep_red,'HandleVisibility','off','LineStyle','--','LineWidth',0.5);
                    line([1 25], [channelSegmentCmapsPosPeak{2}{configIndForFig}{ampIndForFig} channelSegmentCmapsPosPeak{2}{configIndForFig}{ampIndForFig}], 'Color', deep_green,'HandleVisibility','off','LineStyle','--','LineWidth',0.5);
                    line([1 25], [channelSegmentCmapsNegPeak{2}{configIndForFig}{ampIndForFig} channelSegmentCmapsNegPeak{2}{configIndForFig}{ampIndForFig}], 'Color', deep_green,'HandleVisibility','off','LineStyle','--','LineWidth',0.5);
                    line([1 25], [channelSegmentCmapsPosPeak{3}{configIndForFig}{ampIndForFig} channelSegmentCmapsPosPeak{3}{configIndForFig}{ampIndForFig}], 'Color', deep_blue,'HandleVisibility','off','LineStyle','--','LineWidth',0.5);
                    line([1 25], [channelSegmentCmapsNegPeak{3}{configIndForFig}{ampIndForFig} channelSegmentCmapsNegPeak{3}{configIndForFig}{ampIndForFig}], 'Color', deep_blue,'HandleVisibility','off','LineStyle','--','LineWidth',0.5);

                    % Creating textboxes for normalized activations
                    annotation('textbox', [0.62, 0.72, 0.3, 0.05], ...
                        'String', sprintf('Normalized CNAPs:'), ...
                        'Color', 'k', ...        % entire textbox in red
                        'EdgeColor','none', ...       % remove box around it
                        'FontName','Times New Roman',...
                        'FontSize',12);

                    % Annotation #1 (Tibial)
                    annotation('textbox', [0.62, 0.67, 0.3, 0.05], ...
                        'String', sprintf('Tibial: %.3f', Act1), ...
                        'Color', deep_red, ...        % entire textbox in red
                        'EdgeColor','none', ...       % remove box around it
                        'FontName','Times New Roman',...
                        'FontSize',12);

                    % Annotation #2 (Peroneal)
                    annotation('textbox', [0.62, 0.62, 0.3, 0.05], ...
                        'String', sprintf('Peroneal: %.3f', Act2), ...
                        'Color', deep_green, ...
                        'EdgeColor','none',...
                        'FontName','Times New Roman',...
                        'FontSize',12);

                    % Annotation #3 (Sural)
                    annotation('textbox', [0.62, 0.57, 0.3, 0.05], ...
                        'String', sprintf('Sural: %.3f', Act3), ...
                        'Color', deep_blue, ...
                        'EdgeColor','none',...
                        'FontName','Times New Roman',...
                        'FontSize',12);
                    % Creating textboxes for selectivity indices
                    % Annotation #1 (Selectivity Index Label)
                    annotation('textbox', [0.62, 0.4, 0.3, 0.05], ...
                        'String', sprintf('Selectivity Indices:'), ...
                        'Color', 'k', ...        % entire textbox in red
                        'EdgeColor','none', ...       % remove box around it
                        'FontName','Times New Roman',...
                        'FontSize',12);
                    % Annotation #2 (Tibial)
                    annotation('textbox', [0.62, 0.35, 0.3, 0.05], ...
                        'String', sprintf('Tibial: %.3f', Sel1), ...
                        'Color', deep_red, ...        % entire textbox in red
                        'EdgeColor','none', ...       % remove box around it
                        'FontName','Times New Roman',...
                        'FontSize',12);

                    % Annotation #3 (Peroneal)
                    annotation('textbox', [0.62, 0.3, 0.3, 0.05], ...
                        'String', sprintf('Peroneal: %.3f', Sel2), ...
                        'Color', deep_green, ...
                        'EdgeColor','none',...
                        'FontName','Times New Roman',...
                        'FontSize',12);

                    % Annotation #4 (Sural)
                    annotation('textbox', [0.62, 0.25, 0.3, 0.05], ...
                        'String', sprintf('Sural: %.3f', Sel3), ...
                        'Color', deep_blue, ...
                        'EdgeColor','none',...
                        'FontName','Times New Roman',...
                        'FontSize',12);
                    xlim([0 41])
                    %ylim([-4 4])
                    ylim(max(abs(ylim)).*[-1.2 1.2])
                    xlabel('Time (us)','FontName','Times New Roman','FontSize',12)
                    ax = gca;
                    yt = ax.YTick;                   % Get current tick values (in µV)
                    ax.YTickLabel = arrayfun(@(v) num2str(v/10), yt, 'UniformOutput', false);
                    ylabel('Voltage (uV)','FontName','Times New Roman','FontSize', 12);
                    % legend('FontName','Times New Roman')
                    ax = gca;
                    set(ax, 'Box', 'off');  

                    figure('Position',[5 625 550 350])  %To have figure pop on the left side of the screen instead of the middle, easier for visualization
                    plot(reorderedAveragedChannelData{1,2}{configIndForFig}{ampIndForFig},'color', deep_red, 'DisplayName','Tibial','LineWidth',1.5)
                    hold on
                    plot(reorderedAveragedChannelData{2,2}{configIndForFig}{ampIndForFig},'color', deep_green, 'DisplayName','Peroneal','LineWidth',1.5)
                    plot(reorderedAveragedChannelData{3,2}{configIndForFig}{ampIndForFig},'color', deep_blue, 'DisplayName','Sural','LineWidth',1.5)
                    title('Raw Data of CNAPs')
                    xlim([100 300])
                    %ylim([-4 4])
                    ylim(max(abs(ylim)).*[-1.2 1.2])
                    xlabel('Time (us)','FontName','Times New Roman','FontSize',12)
                    ax = gca;
                    yt = ax.YTick;                   % Get current tick values (in µV)
                    ax.YTickLabel = arrayfun(@(v) num2str(v/10), yt, 'UniformOutput', false);
                    ylabel('Voltage (uV)','FontName','Times New Roman','FontSize',12);
                    legend('FontSize',14,'FontName','Times New Roman')
                    ax = gca;
                    set(ax, 'Box', 'off'); 
                catch 
                    disp('Error occurred. Possibly the figure was closed.');
                    break; % Break the loop if an error occurs (e.g., figure closed)
                end
            end
        end
        figNum = figNum + 1
        close all
        end
    end 
end

% save("C:\Users\Green_Group\Zack\ExVivo\Scripts\regressionData","regressionData")

%% Regression Analysis with all data from every ex vivo trial 
load("C:\Users\Green_Group\Zack\ExVivo\Scripts\regressionData.mat",'regressionData')
%% --- Model 1: Ex Vivo — Amplitude & Pulse Widths
lme_exvivo_amp_pw = fitlme(regressionData, ...
    ['Fascicle_Selectivity_ExVivo ~ Fascicle * (Amplitude_Phase1 + PulseWidth_Phase1 + ' ...
     'Amplitude_Phase2 + PulseWidth_Phase2 + Distance_Source_Fascicle + Distance_Sink_Fascicle + ' ...
     'Number_Poles + Phase_Symmetry) + (1|AnimalID)']);

coeffs_exvivo_amp_pw = lme_exvivo_amp_pw.Coefficients;

%% --- Model 2: Ex Vivo — Charge Injection
lme_exvivo_charge = fitlme(regressionData, ...
    ['Fascicle_Selectivity_ExVivo ~ Fascicle * (Charge_Injection + Distance_Source_Fascicle + Distance_Sink_Fascicle + ' ...
     'Number_Poles + Phase_Symmetry) + (1|AnimalID)']);

coeffs_exvivo_charge = lme_exvivo_charge.Coefficients;

%% --- Model 3: Simulation — Amplitude & Pulse Widths
lme_sim_amp_pw = fitlme(regressionData, ...
    ['Fascicle_Selectivity_Simulation ~ Fascicle * (Amplitude_Phase1 + PulseWidth_Phase1 + ' ...
     'Amplitude_Phase2 + PulseWidth_Phase2 + Distance_Source_Fascicle + Distance_Sink_Fascicle + ' ...
     'Number_Poles + Phase_Symmetry) + (1|AnimalID)']);

coeffs_sim_amp_pw = lme_sim_amp_pw.Coefficients;


%% --- Model 4: Simulation — Charge Injection
lme_sim_charge = fitlme(regressionData, ...
    ['Fascicle_Selectivity_Simulation ~ Fascicle * (Charge_Injection + Distance_Source_Fascicle + Distance_Sink_Fascicle + ' ...
     'Number_Poles + Phase_Symmetry) + (1|AnimalID)']);

coeffs_sim_charge = lme_sim_charge.Coefficients;

%% Plotting independent variable vs 
% Define Seaborn-style colors
color_tibial = [0.8941, 0.1020, 0.1098];    % Red
color_peroneal = [0.3020, 0.6863, 0.2902];  % Green
color_sural = [0.2157, 0.4941, 0.7216];     % Blue

% Map fascicles to colors
fascicleColors = containers.Map( ...
    {'Tibial', 'Peroneal', 'Sural'}, ...
    {color_tibial, color_peroneal, color_sural});

% Independent variables to plot
xVars = {
    'Amplitude_Phase1', ...
    'Amplitude_Phase2', ...
    'PulseWidth_Phase1', ...
    'PulseWidth_Phase2', ...
    'Charge_Injection', ...
    'Distance_Source_Fascicle', ...
    'Distance_Sink_Fascicle'};

% Choose selectivity type
selectivityVar = 'Fascicle_Selectivity_ExVivo';  % or 'Fascicle_Selectivity_Simulation'

% Plot layout
numPlots = length(xVars);
nCols = 3;
nRows = ceil(numPlots / nCols);

figure('Name', 'Selectivity vs. Predictors by Fascicle', ...
       'Color', 'w', 'Position', [100 100 1200 800]);

% Plot loop
for i = 1:numPlots
    subplot(nRows, nCols, i); hold on;
    xVar = xVars{i};

    fascicles = {'Tibial', 'Peroneal', 'Sural'};
    for f = 1:length(fascicles)
        fascicleName = fascicles{f};
        mask = regressionData.Fascicle == fascicleName;

        % Use scatter with alpha blending
        scatter(regressionData.(xVar)(mask), ...
                regressionData.(selectivityVar)(mask), ...
                10, ...                        % smaller point size
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', fascicleColors(fascicleName), ...
                'MarkerFaceAlpha', 0.35, ...
                'DisplayName', fascicleName);
        ylim([0 1])
    end

    xlabel(strrep(xVar, '_', '\_'));
    ylabel(strrep(selectivityVar, '_', '\_'));
    title([strrep(selectivityVar, '_', '\_') ' vs. ' strrep(xVar, '_', '\_')]);
    grid on;
    legend('Location', 'best');
end

sgtitle(['Scatter Plots of ' strrep(selectivityVar, '_', '\_') ' vs. Predictors by Fascicle']);



%%
        % lme = fitlme(regressionData, ...
        % ['Fascicle_Selectivity ~ Fascicle * (Amplitude_Phase1 + PulseWidth_Phase1 + ' ...
        %  'Amplitude_Phase2 + PulseWidth_Phase2 + ' ...
        %  'Distance_Source_Fascicle + ' ...
        %  'Distance_Sink_Fascicle + Number_Poles + Phase_Symmetry) + (1|AnimalID)']);
        % 
        % 
        % pw_amp_interaction_lme = fitlme(regressionData, ...
        % ['Fascicle_Selectivity ~ Fascicle * (ChargeInjection + Distance_Source_Fascicle + Distance_Sink_Fascicle + ' ...
        %  'Number_Poles + Phase_Symmetry) + (1|AnimalID)']);
        % 
        % 
        % corrMatrix = corr(table2array(regressionData(:, {'Amplitude_Phase1', 'PulseWidth_Phase1', ...
        %     'Amplitude_Phase2', 'PulseWidth_Phase2'})), 'Rows', 'complete');
        % 
        % disp(array2table(corrMatrix, 'VariableNames', {'Amp1', 'PW1', 'Amp2', 'PW2'}));

% save("C:\Users\Green_Group\Zack\ExVivo\ExVivoData\maxAmpsForNormUserDet_allTrials","maxAmpsForNormUserDet_allTrials")



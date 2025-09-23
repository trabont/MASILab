clear; clc;

% === USER SETTINGS ===
baseDir    = fullfile(pwd,'Processed','MS');        % where subject folders live
msDir = baseDir;
upload     = fullfile(pwd,'Processed','FullFit_MS');
saveRoot   = fullfile(pwd,'Processed','FullFit_MS');   % where outputs go

if ~exist(saveRoot,'dir')
    mkdir(saveRoot);
end

% IDs you want to skip:
excludeIDs = fullfile(pwd,'Processed','FullFit_Analysis','exclusions.m');  
% exclusions is generated after running file_check_all.m

% find all subject folders
d       = dir(upload);
isSubj  = [d.isdir] & ~ismember({d.name},{'.','..'});
subjList= {d(isSubj).name};

for ii = 1:numel(subjList)
    fileID = subjList{ii};

    % skip excluded IDs
    if any(strcmp(fileID, excludeIDs))
        fprintf('Skipping subject %s\n', fileID);
        continue;
    end

    id = str2double(fileID);
    subjDir = fullfile(baseDir, fileID);
    saveDir = fullfile(saveRoot, fileID);
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end

    fullFit(baseDir, saveDir, id);
end
fprintf('MS OVER');

% === USER SETTINGS ===
baseDir    = fullfile(pwd,'Processed','HC');        % where subject folders live
hcDir = baseDir;
upload     = fullfile(pwd,'Processed','FullFit_HC');
saveRoot   = fullfile(pwd,'Processed','FullFit_HC');   % where outputs go

if ~exist(saveRoot,'dir')
    mkdir(saveRoot);
end

% find all subject folders
d       = dir(upload);
isSubj  = [d.isdir] & ~ismember({d.name},{'.','..'});
subjList= {d(isSubj).name};

for ii = 1:numel(subjList)
    fileID = subjList{ii};

    % skip excluded IDs
    if any(strcmp(fileID, excludeIDs))
        fprintf('Skipping subject %s\n', fileID);
        continue;
    end

    id = str2double(fileID);
    subjDir = fullfile(baseDir, fileID);
    saveDir = fullfile(saveRoot, fileID);
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end

    fullFit(baseDir, saveDir, id);
end
saveDir = fullfile(pwd,'Processed','FullFit_Analysis');
fitType = 1;
combine(msDir,hcDir,excludeIDs,saveDir,fitType)

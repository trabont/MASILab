clear; clc;

% === USER SETTINGS ===
baseDir    = fullfile(pwd,'Processed','MS');        % where subject folders live
upload     = fullfile(pwd,'Processed','FullFit_MS');
saveRoot   = fullfile(pwd,'Processed','FullFit_MS');   % where outputs go

if ~exist(saveRoot,'dir')
    mkdir(saveRoot);
end

% IDs you want to skip:
excludeIDs = {'MaskOverlays'};  % <-- modify this list

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

% === USER SETTINGS ===
baseDir    = fullfile(pwd,'Processed','HC');        % where subject folders live
upload     = fullfile(pwd,'Processed','FullFit_HC');
saveRoot   = fullfile(pwd,'Processed','FullFit_HC');   % where outputs go

if ~exist(saveRoot,'dir')
    mkdir(saveRoot);
end

% IDs you want to skip:
excludeIDs = {'MaskOverlays'};  % <-- modify this list

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

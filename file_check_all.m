% file_check_all.m
% ------------------
% Requirements to Run:
%       (1) file_checking
%       (2) mask_overlay.m
%       (3) med_dice_centr.m
% Functionality:
%       (1) First Process MS then HC groups
%       (2) Check file dimensions of registered files (export .xlsx)
%       (3) Using FILE_CHECK.xlsx generate subject exclusion list (based on
%           the dimension "NA"
%       (4) Produce all ROI mask overlays for slices with lymph ROI present
%           using the exclusion list to skip subjects with bad dimensions
%       (5) Produce med_dice_centr excel sheets and figures
%       (6) Produce and save 'exclusions.mat'
% ------------------

clear; clc;

% === USER SETTINGS ===
baseDir    = fullfile(pwd,'Processed','MS');        % where subject folders live
upload     = fullfile(pwd,'Processed','MS');
saveRoot   = fullfile(pwd,'Processed','FullFit_MS');   % where outputs go

file_checking(baseDir, saveRoot);

% find all subject folders
d       = dir(upload);
isSubj  = [d.isdir] & ~ismember({d.name},{'.','..'});
subjList= {d(isSubj).name};

xlsx = fullfile(saveRoot,'FILE_CHECK.xlsx');
% Read both sheets if present
T = table();
sn = sheetnames(xlsx);
if any(strcmpi(sn,'DIMS_SCT')),  T = [T; readtable(xlsx,'Sheet','DIMS_SCT')];  end
if any(strcmpi(sn,'DIMS_ANTS')), T = [T; readtable(xlsx,'Sheet','DIMS_ANTS')]; end

% Treat dims as strings so "NA" is testable regardless of import type
X = string(T.X); Y = string(T.Y); Z = string(T.Z); TT = string(T.T);
missing = (X=="NA") | (Y=="NA") | (Z=="NA") | (TT=="NA");

% Unique IDs to exclude (as a cell array of char)
excludeIDs = cellstr(unique(string(T.ID(missing))))';

mask_overlay(baseDir,subjList,excludeIDs);

for ii = 1:numel(subjList)
    fileID = subjList{ii};

    % skip excluded IDs
    if any(strcmp(fileID, excludeIDs))
        fprintf('Skipping subject %s\n', fileID);
        continue;
    end

    id      = str2double(fileID);
    subjDir = fullfile(baseDir, fileID);
    saveDir = fullfile(saveRoot, fileID);
    if ~exist(saveDir,'dir'), mkdir(saveDir); end

    try
        med_dice_centr(baseDir, saveRoot, id, excludeIDs);
    catch ME
        % Log and keep going
        fprintf(2, 'Subject %s failed in med_dice_centr: %s\n', fileID, ME.message);
        continue;
    end
end


excludeIDs_MS = excludeIDs;
clear("excludeIDs");

% === USER SETTINGS ===
baseDir    = fullfile(pwd,'Processed','HC');        % where subject folders live
upload     = fullfile(pwd,'Processed','HC');
saveRoot   = fullfile(pwd,'Processed','FullFit_HC');   % where outputs go

file_checking(baseDir, saveRoot);

% find all subject folders
d       = dir(upload);
isSubj  = [d.isdir] & ~ismember({d.name},{'.','..'});
subjList= {d(isSubj).name};

xlsx = fullfile(saveRoot,'FILE_CHECK.xlsx');
% Read both sheets if present
T = table();
sn = sheetnames(xlsx);
if any(strcmpi(sn,'DIMS_SCT')),  T = [T; readtable(xlsx,'Sheet','DIMS_SCT')];  end
if any(strcmpi(sn,'DIMS_ANTS')), T = [T; readtable(xlsx,'Sheet','DIMS_ANTS')]; end

% Treat dims as strings so "NA" is testable regardless of import type
X = string(T.X); Y = string(T.Y); Z = string(T.Z); TT = string(T.T);
missing = (X=="NA") | (Y=="NA") | (Z=="NA") | (TT=="NA");

% Unique IDs to exclude (as a cell array of char)
excludeIDs = cellstr(unique(string(T.ID(missing))))';

mask_overlay(baseDir,subjList,excludeIDs);

for ii = 1:numel(subjList)
    fileID = subjList{ii};

    % skip excluded IDs
    if any(strcmp(fileID, excludeIDs))
        fprintf('Skipping subject %s\n', fileID);
        continue;
    end

    id      = str2double(fileID);
    subjDir = fullfile(baseDir, fileID);
    saveDir = fullfile(saveRoot, fileID);
    if ~exist(saveDir,'dir'), mkdir(saveDir); end

    try
        med_dice_centr(baseDir, saveRoot, id, excludeIDs);
    catch ME
        % Log and keep going
        fprintf(2, 'Subject %s failed in med_dice_centr: %s\n', fileID, ME.message);
        continue;
    end
end

exclusions = unique([excludeIDs_MS,excludeIDs]);

exSave = fullfile(pwd,'Processed','FullFit_Analysis');
mkdir(exSave);
outname = fullfile(pwd,'Processed','FullFit_Analysis', sprintf('exclusions.mat'));
save(outname, 'exclusions');

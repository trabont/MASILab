%% combineFullFitMaps_MS.m
% Stacks FullFit .mat files into all_<line>_<ROI>.mat and creates MS_all_tissue
% Then generates figures per lineshape and map with consistent color scaling.

clear; clc;

% === USER SETTINGS ===
fullfitRoot = fullfile(pwd,'Processed','FullFit_MS');   % where all subject sub‑folders live
lines    = {'SL','L','G'};
rois    = {'GM','WM','CSF','LYMPH'};
accRoIs = setdiff(rois,'LYMPH'); % for the combined‐tissue sum, drop LYMPH

mapNames = {'PSR_map','kba_map','T2a_map','T2b_map','R1obs_map'};

% === GET SUBJECT LIST ===
d = dir(fullfitRoot);
isDir = [d.isdir] & ~ismember({d.name},{'.','..'});
subjList = {d(isDir).name};
if isempty(subjList)
    error('No subject folders found in %s', fullfitRoot);
end

% === DETERMINE DIMENSIONS & PREALLOCATE ===
firstLine = lines{1};  
firstRoi = rois{1};
mfFirst = struct('folder',{},'name',{});
for s = 1:numel(subjList)
    subjFolder = fullfile(fullfitRoot, subjList{s});
    pattern = sprintf('%s_slice*_%s_%s.mat', subjList{s}, firstLine, firstRoi);
    dlist = dir(fullfile(subjFolder,pattern));
    for j = 1:numel(dlist)
        mfFirst(end+1) = struct('folder',dlist(j).folder,'name',dlist(j).name);
    end
end
nF = numel(mfFirst);
if nF==0, error('No files found for %s_%s.', firstLine, firstRoi); end
sample = load(fullfile(mfFirst(1).folder,mfFirst(1).name));
[nx,ny] = size(sample.(mapNames{1}));
nM = numel(mapNames);
MS_all_tissue = zeros(nx, ny, numel(lines), nF, nM);

% === STACK & ACCUMULATE MAPS ===
for li = 1:numel(lines)
    line = lines{li};
    fprintf('Processing lineshape %s...\n', line);
    for ri = 1:numel(accRoIs)
        roi = accRoIs{ri};
        % collect files
        mf = struct('folder',{},'name',{});
        for s = 1:numel(subjList)
            subjFolder = fullfile(fullfitRoot, subjList{s});
            pattern = sprintf('%s_slice*_%s_%s.mat', subjList{s}, line, roi);
            dlist = dir(fullfile(subjFolder,pattern));
            for j = 1:numel(dlist)
                mf(end+1) = struct('folder',dlist(j).folder,'name',dlist(j).name);
            end
        end
        nFiles = numel(mf);
        if nFiles==0, warning('No files for %s/%s.', line, roi); continue; end
        if nFiles~=nF, warning('Mismatch %s/%s: %d vs %d files.', line, roi, nFiles,nF); end
        % accumulate
        for k = 1:nFiles
            S = load(fullfile(mf(k).folder,mf(k).name));
            for m = 1:nM
                MS_all_tissue(:,:,li,k,m) = MS_all_tissue(:,:,li,k,m) + S.(mapNames{m});
            end
        end
        % save per-ROI if needed
        data = zeros(nx,ny,nFiles,nM);
        for k=1:nFiles
            S = load(fullfile(mf(k).folder,mf(k).name));
            for m=1:nM, data(:,:,k,m)=S.(mapNames{m}); end
        end
        save(fullfile(fullfitRoot,sprintf('MS_all_%s_%s.mat',line,roi)), 'data','mapNames','mf');
    end
end
save(fullfile(fullfitRoot,'MS_all_tissue.mat'),'MS_all_tissue','mapNames','lines','subjList');



% === PREPARE PLOTTING ===
% get file names from first ROI for labeling
fileNames = {mfFirst.name};
sliceTokens = regexp(fileNames, '_slice(\d+)_', 'tokens', 'once');
sliceNamesTokens = regexp(fileNames, '(\d+)_', 'tokens', 'once');
sliceNames = cellfun(@(c) str2double(c{1}), sliceNamesTokens);
sliceNumbers = cellfun(@(c) str2double(c{1}), sliceTokens);
[xDim,yDim,nLines,nSlices,~] = size(MS_all_tissue);

% Select & compute: PSR, kba, T2aR1 (T2a×R1obs), T2b×1e6
tblIdx    = [1,2,3,4];
plotNames = {'PSR','kba (s^(-1)','T2R1','T2b (micro-s)'};
plotNames2 = {'PSR','kba','T2R1','T2b'};
nPlots    = numel(tblIdx);
plotData  = zeros(xDim,yDim,nLines,nSlices,nPlots);
for p=1:nPlots
    mi = tblIdx(p);
    switch mi
        case 1, plotData(:,:,:,:,p) = MS_all_tissue(:,:,:,:,1);
        case 2, plotData(:,:,:,:,p) = MS_all_tissue(:,:,:,:,2);
        case 3, plotData(:,:,:,:,p) = MS_all_tissue(:,:,:,:,3) .* MS_all_tissue(:,:,:,:,5);
        case 4, plotData(:,:,:,:,p) = MS_all_tissue(:,:,:,:,4) * 1e6;
    end
end
% global limits per plot
% global limits per plot
globalLimits = zeros(nPlots,2);
globalLimits(1,2) = 1;
globalLimits(2,2) = 100;
globalLimits(3,2) = 0.1;
globalLimits(4,2) = 100;


% plot each lineshape & selected map
for li=1:nLines
    lineName = lines{li};
    for p=1:nPlots
        pname = plotNames{p};
        pname2 = plotNames2{p};
        fig = figure('Name',sprintf('%s - %s',lineName,pname),...
             'NumberTitle','off',...
             'Units','normalized',...
             'Position',[0.1 0.1 0.55 0.65],'Visible','off');
        t   = tiledlayout('flow','TileSpacing','compact','Padding','compact');
        t.Title.String = sprintf('Lineshape %s, %s',lineName,pname);
        for si=1:nSlices
            ax = nexttile;
            imagesc(plotData(:,:,li,si,p));
            axis image off;
            title(ax,sprintf('Subject %d, Slice %d', sliceNames(si), sliceNumbers(si)));
            caxis(ax,globalLimits(p,:));
            colormap jet;
            zoom(2);
        end
        cb = colorbar('eastoutside');
        cb.Layout.Tile = 'east';
        saveas(fig, fullfile(fullfitRoot, sprintf('MS_all_%s_%s_MAP.png', lineName, pname2)));
        close(fig);
    end
end

disp('Done plotting selected MS_all_tissue maps.');

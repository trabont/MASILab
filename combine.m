function combine(msDir,hcDir,exclusions,saveDir,fitType)
% combine.m
% Stacks FullFit .mat files into all_<line>_<ROI>.mat and creates HC_all_tissue
% Then generates figures per lineshape and map with consistent color scaling.
clc;

% === USER SETTINGS ===
matRoot = [msDir; hcDir];   % where all subject sub‑folders live
lines    = {'SL','L','G'};
rois    = {'GM','WM','CSF','LYMPH'};
accRoIs = setdiff(rois,'LYMPH'); % for the combined‐tissue sum, drop LYMPH

if fitType == 1 % FullFit
    mapNames = {'PSR_map','kba_map','T2a_map','T2b_map','R1obs_map','chi2_map','chip_map','resn_map'};
    plotNames = {'PSR','kba (s(-1)','T2R1','T2b (micro-s)'};
    plotNames2 = {'PSR','kba','T2R1','T2b'};
elseif fitType == 2 % SinglePT
    mapNames = {'PSR_map','R1obs_map','chi2_map','chip_map','resn_map'};
    plotNames = {'PSR','R1obs'};
    plotNames2 = {'PSR','R1obs'};
end

for i = 1:2
    d = dir(matRoot(i,:));
    isDir = [d.isdir] & ~ismember({d.name},{'.','..'});
    subjList = {d(isDir).name};
    if isempty(subjList)
        error('No subject folders found in %s', matRoot(i,:));
    end
    subjList = setdiff(subjList,exclusions);

    firstLine = lines{1};  
    firstRoi = rois{1};
    mfFirst = struct('folder',{},'name',{});
    for s = 1:numel(subjList)
        subjFolder = fullfile(matRoot(i,:), subjList{s});
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
    all_tissues = zeros(nx, ny, numel(lines), nF, nM);

    % === STACK & ACCUMULATE MAPS ===
    for li = 1:numel(lines)
        line = lines{li};
        fprintf('Processing lineshape %s...\n', line);
        for ri = 1:numel(accRoIs)
            roi = accRoIs{ri};
            % collect files
            mf = struct('folder',{},'name',{});
            for s = 1:numel(subjList)
                subjFolder = fullfile(matRoot(i,:), subjList{s});
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
                    all_tissues(:,:,li,k,m) = all_tissues(:,:,li,k,m) + S.(mapNames{m});
                end
            end
            % save per-ROI if needed
            data = zeros(nx,ny,nFiles,nM);
            for k=1:nFiles
                S = load(fullfile(mf(k).folder,mf(k).name));
                for m=1:nM, data(:,:,k,m)=S.(mapNames{m}); end
            end
    
            if i == 1
                save(fullfile(matRoot(i,:),sprintf('MS_all_%s_%s.mat',line,roi)), 'data','mapNames','mf');
            elseif i == 2
                save(fullfile(matRoot(i,:),sprintf('HC_all_%s_%s.mat',line,roi)), 'data','mapNames','mf');
            end
        end
    end
    
    if i == 1
        save(fullfile(saveDir,'MS_all_tissue.mat'),'all_tissues','mapNames','lines','subjList');
    elseif i == 2
        save(fullfile(saveDir,'HC_all_tissue.mat'),'all_tissues','mapNames','lines','subjList');
    end
    
    % === STACK & ACCUMULATE MAPS FOR LYMPH ===
    for li = 1:numel(lines)
        line = lines{li};
        fprintf('Processing lineshape %s...\n', line);
        roi = rois{4};
        % collect files
        mf = struct('folder',{},'name',{});
        for s = 1:numel(subjList)
            subjFolder = fullfile(matRoot(i,:), subjList{s});
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
                all_tissues(:,:,li,k,m) = all_tissues(:,:,li,k,m) + S.(mapNames{m});
            end
        end
        % save per-ROI if needed
        data = zeros(nx,ny,nFiles,nM);
        for k=1:nFiles
            S = load(fullfile(mf(k).folder,mf(k).name));
            for m=1:nM, data(:,:,k,m)=S.(mapNames{m}); end
        end
        if i == 1
            save(fullfile(matRoot(i,:),sprintf('MS_all_%s_%s.mat',line,roi)), 'data','mapNames','mf');
        elseif i == 2
            save(fullfile(matRoot(i,:),sprintf('HC_all_%s_%s.mat',line,roi)), 'data','mapNames','mf');
        end
    end
    
    % === PREPARE PLOTTING ===
    % get file names from first ROI for labeling
    fileNames = {mfFirst.name};
    sliceTokens = regexp(fileNames, '_slice(\d+)_', 'tokens', 'once');
    sliceNamesTokens = regexp(fileNames, '(\d+)_', 'tokens', 'once');
    sliceNames = cellfun(@(c) str2double(c{1}), sliceNamesTokens);
    sliceNumbers = cellfun(@(c) str2double(c{1}), sliceTokens);
    [xDim,yDim,nLines,nSlices,~] = size(all_tissues);
    
    % Select & compute: PSR, kba, T2aR1 (T2a×R1obs), T2b×1e6
    tblIdx    = 1:1:length(plotNames);
    nPlots    = numel(tblIdx);
    globalLimits = zeros(nPlots,2); % global limits per plot
    plotData  = zeros(xDim,yDim,nLines,nSlices,nPlots);
    for p=1:nPlots
        mi = plotNames2(p);
        switch mi
            case 'PSR'
              plotData(:,:,:,:,p) = all_tissues(:,:,:,:,1);
              globalLimits(1,2) = 1;
            case 'kba'
              plotData(:,:,:,:,p) = all_tissues(:,:,:,:,2);
              globalLimits(2,2) = 100;
            case 'T2R1'
              plotData(:,:,:,:,p) = all_tissues(:,:,:,:,3) .* all_tissues(:,:,:,:,5);
              globalLimits(3,2) = 0.1;
            case 'T2b_us'
              plotData(:,:,:,:,p) = all_tissues(:,:,:,:,4) * 1e6;
              globalLimits(4,2) = 100;
            case 'R1obs'
              plotData(:,:,:,:,p) = all_tissues(:,:,:,:,5);
              globalLimits(2,2) = 100;
            otherwise
              continue;
        end
    end 
    
    % plot each lineshape & selected map
    for li=1:nLines
        lineName = lines{li};
        for p=1:nPlots
            pname = plotNames{p};
            pname2 = plotNames2{p};
            fig = figure('Name',sprintf('%s - %s',lineName,pname),...
                 'NumberTitle','off',...
                 'Units','normalized',...
                 'Position',[0.1 0.1 0.55 0.65], 'Visible','off');
            t   = tiledlayout('flow','TileSpacing','compact','Padding','compact');
            t.Title.String = sprintf('Lineshape %s, %s',lineName,pname);
            for si=1:nSlices
                ax = nexttile;
                imagesc(plotData(:,:,li,si,p));
                axis image off;
                title(ax,sprintf('Subject %d, Slice %d', sliceNames(si), sliceNumbers(si)));
                clim(ax,globalLimits(p,:));
                colormap jet;
                zoom(2);
            end
            cb = colorbar('eastoutside');
            cb.Layout.Tile = 'east';
            if i == 1
                saveas(fig, fullfile(saveDir, sprintf('MS_all_%s_%s_MAP.png', lineName, pname2)));
                close(fig);
            elseif i == 2
                saveas(fig, fullfile(saveDir, sprintf('HC_all_%s_%s_MAP.png', lineName, pname2)));
                close(fig);
            end
        end
    end
    
    disp('Done combining and producing maps.');
end

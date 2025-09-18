%% histBoxWhisker.m
% Load outputs from combine.m and produce:
%  1) per-group (HC/MS) summaries for ALL-TISSUE (med/iqr/mad/mean/std/n)
%  2) per-group summaries BY TISSUE (GM/WM/CSF/LYMPH)
%  3) HC vs MS comparisons (median diff, ranksum p, Cliff's delta)
%  4) Optional hist+KDE overlays for ALL-TISSUE (HC vs MS) per lineshape & param
%
% Run AFTER you've created:
%   Processed/FullFit_MS/MS_all_tissue.mat
%   Processed/FullFit_HC/HC_all_tissue.mat
%   and the per-ROI stacks: HC_all_<line>_<ROI>.mat, MS_all_<line>_<ROI>.mat
%
% Outputs go to: Processed/FullFit_Summary/
%
% NOTE: This script drops zeros & non-finite values. Optional bounds &
% optional outlier trimming are configurable below.

clear; clc;

%% ===================== USER SETTINGS =====================
msRoot = fullfile(pwd,'Processed','FullFit_MS');
hcRoot = fullfile(pwd,'Processed','FullFit_HC');
outDir = fullfile(pwd,'Processed','FullFit_Summary');
if ~exist(outDir,'dir'), mkdir(outDir); end

lineshapes = {'SL','L','G'};           % must match combine.m
rois       = {'GM','WM','CSF','LYMPH'}; % per-ROI files exist for each

% Parameters we will analyze (derived from maps in *_all_tissue)
% Map order in combine.m: {'PSR_map','kba_map','T2a_map','T2b_map','R1obs_map'}
paramList = {'PSR','kba','T2R1','T2b_us'};  % T2R1 = T2a * R1obs; T2b_us = T2b * 1e6

% Optional bounds (display-like sanity window). Set useBounds=false to ignore.
useBounds = true;
bounds.PSR    = [0, 1];
bounds.kba    = [0, 50];
bounds.T2R1   = [0, 0.1];
bounds.T2b_us = [0, 100];  % microseconds

% Outlier trimming (applied AFTER bounds). Options: 'none' | 'tukey' | 'mad'
trimMethod = 'none';     % change to 'tukey' or 'mad' if wanted
madK       = 3;          % if trimMethod='mad', keep |x-med| <= madK*MAD

% Histogram bin widths: set useFixedBW=true to force exact widths.
useFixedBW = false;
fixedBW.PSR    = 0.006;   % similar to Alex's paper
fixedBW.kba    = 0.65;
fixedBW.T2R1   = 0.01;
fixedBW.T2b_us = 0.5;

% Plot settings
makePlots   = true;     % hist+KDE overlays for ALL-TISSUE only (HC vs MS)
maxNBins    = 100;      % cap for FD-derived bins
minNBins    = 8;        % floor for FD-derived bins

%% ===================== LOAD DATA =====================
MS_all = load(fullfile(msRoot,'MS_all_tissue.mat'));  % MS_all_tissue, mapNames, lines, subjList
HC_all = load(fullfile(hcRoot,'HC_all_tissue.mat'));

assert(isfield(MS_all,'MS_all_tissue'), 'Missing MS_all_tissue.mat contents.');
assert(isfield(HC_all,'HC_all_tissue'), 'Missing HC_all_tissue.mat contents.');

MS_all_tissue = MS_all.MS_all_tissue;  % [nx,ny,nLines,nFiles,nMaps]
HC_all_tissue = HC_all.HC_all_tissue;
mapNames      = MS_all.mapNames;

[nx,ny,nLines,nFiles,nMaps] = size(MS_all_tissue);
assert(nLines == numel(lineshapes), 'lineshape count mismatch.');

%% ===================== HELPERS =====================
% (local function signatures provided at end of file)

%% ===================== SUMMARIES: ALL-TISSUE =====================
% This section produces the histogram + KDE images of all data for all
% parameters. The tissues are all grouped together.
% In addition this section generates a csv file that compares all tissues
rows = {};
compRowsAll = {};
for li = 1:numel(lineshapes)
  line = lineshapes{li};
  for pi = 1:numel(paramList)
    pname = paramList{pi};

    x_MS = gather_all_tissue(MS_all_tissue, li, pname);
    x_HC = gather_all_tissue(HC_all_tissue, li, pname);

    x_MS = clean_and_trim(x_MS, pname, useBounds, bounds, trimMethod, madK);
    x_HC = clean_and_trim(x_HC, pname, useBounds, bounds, trimMethod, madK);

    S_MS = summarize_vec(x_MS);
    S_HC = summarize_vec(x_HC);

    [pval, dCliff, medDiff] = compare_groups(x_HC, x_MS);

    rows(end+1,:) = {line, pname, 'MS', S_MS.n, S_MS.median, S_MS.iqr(1), S_MS.iqr(2), S_MS.mad, S_MS.mean, S_MS.std};
    rows(end+1,:) = {line, pname, 'HC', S_HC.n, S_HC.median, S_HC.iqr(1), S_HC.iqr(2), S_HC.mad, S_HC.mean, S_HC.std};

    compRowsAll(end+1,:) = {line, pname, S_HC.n, S_MS.n, medDiff, pval, dCliff};

    if makePlots
      fig = figure('Visible','off','Units','normalized','Position',[0.1 0.1 0.5 0.6], 'Name',sprintf('%s-%s HCvsMS',line,pname));
      ax = axes(fig); hold(ax,'on'); box(ax,'on');

      % Determine edges
      edges = choose_edges(x_HC, x_MS, pname, useFixedBW, fixedBW, useBounds, bounds, minNBins, maxNBins);

      % Histograms (PDF)
      h1 = histogram(ax, x_HC, 'BinEdges', edges, 'Normalization','pdf');
      h2 = histogram(ax, x_MS, 'BinEdges', edges, 'Normalization','pdf');
      % Set colors: HC blue, MS red
      set(h1,'FaceAlpha',0.35,'EdgeAlpha',0.2,'FaceColor',[0 0 1]);
      set(h2,'FaceAlpha',0.35,'EdgeAlpha',0.2,'FaceColor',[1 0 0]);

      % KDE
      try
        [f1,xi1] = ksdensity(x_HC);
        [f2,xi2] = ksdensity(x_MS);
        plot(ax, xi1, f1, 'LineWidth',1.75, 'Color',[0 0 1]);
        plot(ax, xi2, f2, 'LineWidth',1.75, 'Color',[1 0 0]);
      catch
        % ksdensity may fail if too few points; ignore gracefully
      end

      title(ax, sprintf('All-Tissue  |  %s  |  %s', line, pname));
      xlabel(ax, pname); ylabel(ax,'PDF');
      lg = legend(ax, {'HC hist','MS hist','HC KDE','MS KDE'}, 'Location','best');

      % Save
      outPng = fullfile(outDir, sprintf('ALLTISSUE_%s_%s_HCvMS_histKDE.png', line, pname));
      saveas(fig, outPng);
      close(fig);
    end
  end
end

T_all = cell2table(rows, 'VariableNames',{'Lineshape','Param','Group','N','Median','IQR_lo','IQR_hi','MAD','Mean','Std'});
writetable(T_all, fullfile(outDir,'summary_all_tissue.csv'));

T_comp_all = cell2table(compRowsAll, 'VariableNames',{'Lineshape','Param','N_HC','N_MS','MedianDiff_MS_minus_HC','RankSum_p','CliffsDelta'});
writetable(T_comp_all, fullfile(outDir,'compare_all_tissue.csv'));

%% ===================== SUMMARIES: BY TISSUE =====================
% This section generates a csv file that compares things base on tissues
rows = {}; compRows = {};
for li = 1:numel(lineshapes)
  line = lineshapes{li};
  for ri = 1:numel(rois)
    roi = rois{ri};

    % load per-ROI stacks for both groups (if present)
    fHC = fullfile(hcRoot, sprintf('HC_all_%s_%s.mat', line, roi));
    fMS = fullfile(msRoot, sprintf('MS_all_%s_%s.mat', line, roi));

    if ~exist(fHC,'file') || ~exist(fMS,'file')
      warning('Missing per-ROI file(s) for %s/%s; skipping.', line, roi);
      continue;
    end
    DHC = load(fHC);  % fields: data [nx,ny,nFiles,nMaps], mapNames, mf
    DMS = load(fMS);

    for pi = 1:numel(paramList)
      pname = paramList{pi};
      x_HC = gather_roi(DHC.data, pname);
      x_MS = gather_roi(DMS.data, pname);

      x_HC = clean_and_trim(x_HC, pname, useBounds, bounds, trimMethod, madK);
      x_MS = clean_and_trim(x_MS, pname, useBounds, bounds, trimMethod, madK);

      S_HC = summarize_vec(x_HC);
      S_MS = summarize_vec(x_MS);
      [pval, dCliff, medDiff] = compare_groups(x_HC, x_MS);

      rows(end+1,:) = {line, roi, pname, 'HC', S_HC.n, S_HC.median, S_HC.iqr(1), S_HC.iqr(2), S_HC.mad, S_HC.mean, S_HC.std};
      rows(end+1,:) = {line, roi, pname, 'MS', S_MS.n, S_MS.median, S_MS.iqr(1), S_MS.iqr(2), S_MS.mad, S_MS.mean, S_MS.std};

      compRows(end+1,:) = {line, roi, pname, S_HC.n, S_MS.n, medDiff, pval, dCliff};
    end
  end
end

T_roi = cell2table(rows, 'VariableNames',{'Lineshape','ROI','Param','Group','N','Median','IQR_lo','IQR_hi','MAD','Mean','Std'});
writetable(T_roi, fullfile(outDir,'summary_by_tissue.csv'));

T_comp_roi = cell2table(compRows, 'VariableNames',{'Lineshape','ROI','Param','N_HC','N_MS','MedianDiff_MS_minus_HC','RankSum_p','CliffsDelta'});
writetable(T_comp_roi, fullfile(outDir,'compare_by_tissue.csv'));

%% ===================== RAW BOX+POINTS (visual only) =====================
% This section produces boxplot figures of the parameters 
rawBox_enable = true;                  % turn off if you don't need these plots
rawBox_params = {'PSR','kba','T2R1','T2b_us'};
raw_includeZeros = false;               % include zeros like your example
raw_useBoundsForYLim = true;           % y-limits follow bounds for consistency
raw_maxPoints = 30000;                 % subsample per group if larger (for speed)

if rawBox_enable
  for li=1:numel(lineshapes)
    line = lineshapes{li};
    for ri=1:numel(rois)
      roi = rois{ri};

      fHC = fullfile(hcRoot, sprintf('HC_all_%s_%s.mat', line, roi));
      fMS = fullfile(msRoot, sprintf('MS_all_%s_%s.mat', line, roi));
      if ~exist(fHC,'file') || ~exist(fMS,'file')
        continue;
      end
      DHC = load(fHC);  % data: [nx,ny,nFiles,nMaps]
      DMS = load(fMS);

      for pi=1:numel(rawBox_params)
        pname = rawBox_params{pi};

        x_HC = gather_roi(DHC.data, pname);
        x_MS = gather_roi(DMS.data, pname);

        % keep finite; optionally keep zeros
        x_HC = x_HC(isfinite(x_HC));
        x_MS = x_MS(isfinite(x_MS));
        if ~raw_includeZeros
          x_HC = x_HC(x_HC~=0);
          x_MS = x_MS(x_MS~=0);
        end

        % (Optional) subsample points for plotting speed/clarity
        if numel(x_HC) > raw_maxPoints
          idx = randperm(numel(x_HC), raw_maxPoints); x_HC = x_HC(idx);
        end
        if numel(x_MS) > raw_maxPoints
          idx = randperm(numel(x_MS), raw_maxPoints); x_MS = x_MS(idx);
        end

        % Figure
        fig = figure('Visible','off','Units','normalized','Position',[0.08 0.08 0.55 0.72], ...
                     'Name', sprintf('RAW %s — %s (%s)', roi, pname, line));
        ax = axes(fig); hold(ax,'on'); box(ax,'on');

        % Jittered points
        jitter = 0.18;
        n1 = numel(x_HC); n2 = numel(x_MS);
        x1 = 1 + (rand(n1,1)-0.5)*2*jitter;
        x2 = 2 + (rand(n2,1)-0.5)*2*jitter;
        s1 = scatter(ax, x1, x_HC, 60, 'filled', 'MarkerFaceAlpha',0.25, 'MarkerEdgeAlpha',0.15, 'MarkerFaceColor',[0 0 1]);
        s2 = scatter(ax, x2, x_MS, 60, 'filled', 'MarkerFaceAlpha',0.25, 'MarkerEdgeAlpha',0.15, 'MarkerFaceColor',[1 0 0]);

        % Boxplots with per-group colors (pad with NaNs to equal length)
        maxLen = max(n1, n2);
        M = NaN(maxLen, 2);
        M(1:n1,1) = x_HC; M(1:n2,2) = x_MS;
        boxplot(ax, M, 'Colors', [0 0 0; 0 0 0], 'Symbol','', 'Positions',[1 2], 'Widths',0.6);
        set(findobj(ax,'Tag','Box'),'LineWidth',1.25);

        % Axes & labels
        xlim(ax,[0.5 2.5]); set(ax,'XTick',[1 2],'XTickLabel',{'HC','MS'});
        ylabel(ax, pname);
        title(ax, sprintf('RAW %s — %s (%s)', roi, pname, line));
        if raw_useBoundsForYLim
          b = bounds.(pname); ylim(ax, b);
        end
        legend(ax, [s1 s2], {'HC','MS'}, 'Location','northoutside','Orientation','horizontal');

        % Save
        outP = fullfile(outDir, sprintf('RAW_%s_%s_%s_BoxPoints.png', line, roi, strrep(pname,'_us','')));
        saveas(fig, outP); close(fig);
      end
    end
  end
end

%% ===================== HCMS BOX+POINTS PANELS (32 figs) =====================
% Figures with HC vs MS side-by-side **with scatter points** and median labels.
hcmsPanels_enable   = true;
hcms_params         = {'PSR','kba','T2R1','T2b_us'};   % parameters
hcms_lineshapes     = {'SL','L','G'};                   % lineshapes
hcms_rois           = {'WM','GM','CSF','LYMPH'};        % tissues
hcms_includeZeros   = false;                            % drop zeros
hcms_useBoundsYLim  = false;                            % use bounds for y-limits
hcms_scatterSize    = 6;                                % point size
hcms_jitter         = 0.12;                             % half-width jitter

if hcmsPanels_enable
  make_boxpanels_hc_ms(hcRoot, msRoot, hcms_params, hcms_lineshapes, hcms_rois, ...
    hcms_includeZeros, hcms_useBoundsYLim, bounds, outDir, hcms_scatterSize, hcms_jitter);
end

%% ===================== HCMS HIST+KDE PANELS (32 figs) =====================
% Figures with HC vs MS side-by-side with KDE lineshape.
histPanels_enable     = true;
hist_params           = {'PSR','kba','T2R1','T2b_us'};   % parameters
hist_lineshapes       = {'SL','L','G'};                   % lineshapes
hist_rois             = {'WM','GM','CSF','LYMPH'};        % tissues
hist_includeZeros     = false;                            % drop zeros for density

if histPanels_enable
  make_hist_panels_hc_ms(hcRoot, msRoot, hist_params, hist_lineshapes, hist_rois, ...
    hist_includeZeros, useBounds, bounds, useFixedBW, fixedBW, minNBins, maxNBins, outDir);
end

%% ===================== DONE =====================
disp('Summaries written to Processed/FullFit_Summary/.');

%% ===================== LOCAL FUNCTIONS =====================
function make_hist_panels_hc_ms(hcRoot, msRoot, params, lineshapes, rois, includeZeros, useBounds, bounds, useFixedBW, fixedBW, minNBins, maxNBins, outDir)
% Build 32 histogram+KDE figures mirroring the box panels layout.
% A) Param + Lineshape: per (param,line) -> 4 tissues (4 subplots)
% B) Param only: per param -> 12 categories (3x4 grid of line×tissue)
% C) Param + Tissue: per (param,tissue) -> 3 lines (1x3 grid)

  blue = [0 0 1]; red = [1 0 0];

  % ---------- A) Param + Lineshape (4 tissues) ----------
  for pi=1:numel(params)
    pname = params{pi};
    for li=1:numel(lineshapes)
      line = lineshapes{li};
      hcC = cell(1,numel(rois)); msC = hcC;
      for ri=1:numel(rois)
        [hcC{ri}, msC{ri}] = collect_pair(hcRoot, msRoot, line, rois{ri}, pname, includeZeros);
      end
      if all(cellfun(@isempty,hcC)) && all(cellfun(@isempty,msC)), continue; end
      % unified edges across tissues for this (param,line)
      edges = edges_from_cells(hcC, msC, pname, useFixedBW, fixedBW, useBounds, bounds, minNBins, maxNBins);

      fig = figure('Visible','off','Units','normalized','Position',[0.05 0.07 0.72 0.74], ...
                   'Name', sprintf('%s × %s × All Tissues (HCvsMS)', pname, line));
      t = tiledlayout(fig, 2, 2, 'TileSpacing','compact','Padding','compact');
      title(t, sprintf('%s × %s × All Tissues', pname, line));
      for ri=1:numel(rois)
        ax = nexttile(t);
        draw_hist_kde(ax, hcC{ri}, msC{ri}, edges, blue, red);
        xlabel(ax, pretty_name(pname)); ylabel(ax,'PDF');
        title(ax, rois{ri});
        xlim(ax,[edges(1) edges(end)]);
      end
      saveas(fig, fullfile(outDir, sprintf('HIST_%sx%sxAllTissues_HCvMS.png', strip_us(pname), line)));
      close(fig);
    end
  end

  % ---------- B) Param only (12 cats = 3 lines × 4 tissues) ----------
  for pi=1:numel(params)
    pname = params{pi};
    hcC = cell(1,numel(lineshapes)*numel(rois)); msC = hcC; labels = msC; idx=1;
    for li=1:numel(lineshapes)
      for ri=1:numel(rois)
        line = lineshapes{li}; roi = rois{ri};
        [hcC{idx}, msC{idx}] = collect_pair(hcRoot, msRoot, line, roi, pname, includeZeros);
        labels{idx} = sprintf('%s-%s', line, roi); idx=idx+1;
      end
    end
    if all(cellfun(@isempty,hcC)) && all(cellfun(@isempty,msC)), continue; end
    edges = edges_from_cells(hcC, msC, pname, useFixedBW, fixedBW, useBounds, bounds, minNBins, maxNBins);

    fig = figure('Visible','off','Units','normalized','Position',[0.03 0.06 0.92 0.78], ...
                 'Name', sprintf('%s × All Lineshapes × All Tissues (HCvsMS)', pname));
    t = tiledlayout(fig, 3, 4, 'TileSpacing','compact','Padding','compact');
    title(t, sprintf('%s × All Lineshapes × All Tissues', pname));
    for j=1:numel(hcC)
      ax = nexttile(t);
      draw_hist_kde(ax, hcC{j}, msC{j}, edges, blue, red);
      xlabel(ax, pretty_name(pname)); ylabel(ax,'PDF');
      title(ax, labels{j});
      xlim(ax,[edges(1) edges(end)]);
    end
    saveas(fig, fullfile(outDir, sprintf('HIST_%sxAllLineshapesxAllTissues_HCvMS.png', strip_us(pname))));
    close(fig);
  end

  % ---------- C) Param + Tissue (3 lines) ----------
  for pi=1:numel(params)
    pname = params{pi};
    for ri=1:numel(rois)
      roi = rois(ri);
    end
  end
  for pi=1:numel(params)
    pname = params{pi};
    for ri=1:numel(rois)
      roi = rois{ri};
      hcC = cell(1,numel(lineshapes)); msC = hcC;
      for li=1:numel(lineshapes)
        line = lineshapes{li};
        [hcC{li}, msC{li}] = collect_pair(hcRoot, msRoot, line, roi, pname, includeZeros);
      end
      if all(cellfun(@isempty,hcC)) && all(cellfun(@isempty,msC)), continue; end
      edges = edges_from_cells(hcC, msC, pname, useFixedBW, fixedBW, useBounds, bounds, minNBins, maxNBins);

      fig = figure('Visible','off','Units','normalized','Position',[0.06 0.12 0.7 0.5], ...
                   'Name', sprintf('%s × All Lineshapes × %s (HCvsMS)', pname, roi));
      t = tiledlayout(fig, 1, 3, 'TileSpacing','compact','Padding','compact');
      title(t, sprintf('%s × All Lineshapes × %s', pname, roi));
      for li=1:numel(lineshapes)
        ax = nexttile(t);
        draw_hist_kde(ax, hcC{li}, msC{li}, edges, blue, red);
        xlabel(ax, pretty_name(pname)); ylabel(ax,'PDF');
        title(ax, lineshapes{li});
        xlim(ax,[edges(1) edges(end)]);
      end
      saveas(fig, fullfile(outDir, sprintf('HIST_%sxAllLineshapesx%s_HCvMS.png', strip_us(pname), roi)));
      close(fig);
    end
  end
end

function edges = edges_from_cells(hcC, msC, pname, useFixedBW, fixedBW, useBounds, bounds, minNBins, maxNBins)
% Build common edges across multiple categories using existing choose_edges logic.
  xHC = []; xMS = [];
  for j=1:numel(hcC)
    if ~isempty(hcC{j}), xHC = [xHC; hcC{j}(:)]; end
    if ~isempty(msC{j}), xMS = [xMS; msC{j}(:)]; end
  end
  if isempty(xHC) && isempty(xMS)
    edges = 0:0.1:1; return;
  end
  edges = choose_edges(xHC, xMS, pname, useFixedBW, fixedBW, useBounds, bounds, minNBins, maxNBins);
end

function draw_hist_kde(ax, xHC, xMS, edges, blue, red)
% Draw overlapped histograms + KDE for one category.
  if nargin<5, blue=[0 0 1]; red=[1 0 0]; end
  hold(ax,'on'); box(ax,'on');
  if ~isempty(xHC)
    histogram(ax, xHC, 'BinEdges', edges, 'Normalization','pdf', 'FaceAlpha',0.35, 'EdgeAlpha',0.15, 'FaceColor',blue);
    try
      [f,xi] = ksdensity(xHC); plot(ax, xi, f, 'LineWidth',1.5, 'Color',blue);
    catch, end
  end
  if ~isempty(xMS)
    histogram(ax, xMS, 'BinEdges', edges, 'Normalization','pdf', 'FaceAlpha',0.35, 'EdgeAlpha',0.15, 'FaceColor',red);
    try
      [f,xi] = ksdensity(xMS); plot(ax, xi, f, 'LineWidth',1.5, 'Color',red);
    catch, end
  end
end

function name = pretty_name(p)
  if strcmp(p,'T2b_us'), name = 'T2b (µs)'; else, name = p; end
end

function s = strip_us(p)
  s = strrep(p,'_us','');
end

function make_boxpanels_hc_ms(hcRoot, msRoot, params, lineshapes, rois, includeZeros, useBounds, bounds, outDir, sz, jitter)
% Build 32 figures with HC vs MS, scatter points, and median labels.
% A) Param + Lineshape → per param & line: 4 tissues ⇒ 8 boxes
% B) Param only → per param: 12 categories (line×tissue) ⇒ 24 boxes
% C) Param + Tissue → per param & tissue: 3 lines ⇒ 6 boxes

  % ---- A) Param + Lineshape ----
  for pi=1:numel(params)
    pname = params{pi};
    for li=1:numel(lineshapes)
      line = lineshapes{li};
      catLbl = rois; hcC = cell(1,numel(rois)); msC = cell(1,numel(rois));
      for ri=1:numel(rois)
        roi = rois{ri};
        [hcC{ri}, msC{ri}] = collect_pair(hcRoot, msRoot, line, roi, pname, includeZeros);
      end
      titleStr = sprintf('%s × %s × All Tissues', pname, line);
      outP = fullfile(outDir, sprintf('HCMS_%sx%sxAllTissues_BoxPoints.png', strip_us(pname), line));
      plot_hcms_panel(catLbl, hcC, msC, pname, titleStr, useBounds, bounds, outP, sz, jitter);
    end
  end

  % ---- B) Param only (All Lineshapes × All Tissues) ----
  for pi=1:numel(params)
    pname = params{pi};
    catLbl = cell(1,numel(lineshapes)*numel(rois));
    hcC = catLbl; msC = catLbl; idx=1;
    for li=1:numel(lineshapes)
      for ri=1:numel(rois)
        line = lineshapes{li}; roi = rois{ri};
        [hcC{idx}, msC{idx}] = collect_pair(hcRoot, msRoot, line, roi, pname, includeZeros);
        catLbl{idx} = sprintf('%s-%s', line, roi); idx=idx+1;
      end
    end
    titleStr = sprintf('%s × All Lineshapes × All Tissues', pname);
    outP = fullfile(outDir, sprintf('HCMS_%sxAllLineshapesxAllTissues_BoxPoints.png', strip_us(pname)));
    plot_hcms_panel(catLbl, hcC, msC, pname, titleStr, useBounds, bounds, outP, sz, jitter);
  end

  % ---- C) Param + Tissue (All Lineshapes × tissue) ----
  for pi=1:numel(params)
    pname = params{pi};
    for ri=1:numel(rois)
      roi = rois{ri}; catLbl = lineshapes; hcC = cell(1,numel(lineshapes)); msC = hcC;
      for li=1:numel(lineshapes)
        line = lineshapes{li};
        [hcC{li}, msC{li}] = collect_pair(hcRoot, msRoot, line, roi, pname, includeZeros);
      end
      titleStr = sprintf('%s × All Lineshapes × %s', pname, roi);
      outP = fullfile(outDir, sprintf('HCMS_%sxAllLineshapesx%s_BoxPoints.png', strip_us(pname), roi));
      plot_hcms_panel(catLbl, hcC, msC, pname, titleStr, useBounds, bounds, outP, sz, jitter);
    end
  end
end

function [vHC, vMS] = collect_pair(hcRoot, msRoot, line, roi, pname, includeZeros)
% Return vectors for HC and MS for a given line/roi/param.
  vHC = []; vMS = [];
  fHC = fullfile(hcRoot, sprintf('HC_all_%s_%s.mat', line, roi));
  if exist(fHC,'file')
    D = load(fHC);
    try
      vHC = gather_roi(D.data, pname);
    catch, vHC = []; end
  end
  fMS = fullfile(msRoot, sprintf('MS_all_%s_%s.mat', line, roi));
  if exist(fMS,'file')
    D = load(fMS);
    try
      vMS = gather_roi(D.data, pname);
    catch, vMS = []; end
  end
  vHC = vHC(isfinite(vHC)); vMS = vMS(isfinite(vMS));
  if ~includeZeros
    vHC = vHC(vHC~=0); vMS = vMS(vMS~=0);
  end
end

function plot_hcms_panel(catLabels, hcCells, msCells, pname, titleStr, useBounds, bounds, outPath, sz, jitter)
% Draw side-by-side HC (blue) and MS (red) boxes per category, with scatter + medians.
  nCat = numel(catLabels);
  % Prepare NaN-padded matrices
  kHC = 0; kMS = 0;
  for j=1:nCat
    kHC = max(kHC, numel(hcCells{j}));
    kMS = max(kMS, numel(msCells{j}));
  end
  Mhc = NaN(kHC, nCat); Mms = NaN(kMS, nCat);
  for j=1:nCat
    v1 = hcCells{j}; v2 = msCells{j};
    if ~isempty(v1), Mhc(1:numel(v1),j) = v1(:); end
    if ~isempty(v2), Mms(1:numel(v2),j) = v2(:); end
  end

  figW = 0.1 + 0.04*max(12, 2*nCat);  % scale width a bit for many boxes
  figW = min(figW, 0.95);
  fig = figure('Visible','off','Units','normalized','Position',[0.03 0.08 figW 0.7], ...
               'Name', titleStr);
  ax = axes(fig); hold(ax,'on'); box(ax,'on');

  posBase = 1:nCat;
  posHC = posBase - 0.17; posMS = posBase + 0.17;

  % Plot boxes separately to color them
  boxplot(ax, Mhc, 'Positions', posHC, 'Colors', [0 0 0], 'Symbol','', 'Widths', 0.28);
  boxplot(ax, Mms, 'Positions', posMS, 'Colors', [0 0 0], 'Symbol','', 'Widths', 0.28);
  set(findobj(ax,'Tag','Box'),'LineWidth',1.1);

  % Scatter overlays
  for j=1:nCat
    if ~all(isnan(Mhc(:,j)))
      y = Mhc(~isnan(Mhc(:,j)), j); x = posHC(j) + (rand(size(y))-0.5)*2*jitter;
      scatter(ax, x, y, sz, 'filled', 'MarkerFaceAlpha',0.25, 'MarkerEdgeAlpha',0.12, 'MarkerFaceColor',[0 0 1]);
    end
    if ~all(isnan(Mms(:,j)))
      y = Mms(~isnan(Mms(:,j)), j); x = posMS(j) + (rand(size(y))-0.5)*2*jitter;
      scatter(ax, x, y, sz, 'filled', 'MarkerFaceAlpha',0.25, 'MarkerEdgeAlpha',0.12, 'MarkerFaceColor',[1 0 0]);
    end
  end

  % Median labels
  for j=1:nCat
    medHC = median(Mhc(:,j),'omitnan');
    medMS = median(Mms(:,j),'omitnan');
    if isfinite(medHC)
      text(ax, posHC(j), medHC, sprintf('%.3g', medHC), 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'Color',[0 0 0], 'FontSize',9,'FontWeight','bold');
    end
    if isfinite(medMS)
      text(ax, posMS(j), medMS, sprintf('%.3g', medMS), 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'Color',[0 0 0], 'FontSize',9,'FontWeight','bold');
    end
  end

  xlim(ax, [0.5, nCat+0.5]);
  set(ax,'XTick',posBase,'XTickLabel',catLabels);
  xtickangle(ax, min(35, 10 + nCat));
  ylabel(ax, pretty_name(pname));
  title(ax, titleStr);
  if useBounds && isfield(bounds,pname), ylim(ax, bounds.(pname)); end
  legend(ax, {'HC','MS'}, 'Location','northoutside','Orientation','horizontal');

  saveas(fig, outPath);
  close(fig);
end

function make_box_panels(hcRoot, msRoot, params, lineshapes, rois, includeZeros, useBounds, bounds, outDir)
% Generates the 32 requested figures of simple boxplots combining HC+MS voxels.
% A) Parameter + Lineshape conditional: for each param & line -> 4 boxes (tissues)
% B) Parameter conditional: each param -> 12 boxes (line×tissue)
% C) Parameter + Tissue conditional: for each param & tissue -> 3 boxes (lineshapes)

  % ---------- A) Param + Lineshape (4 boxes) ----------
  for pi=1:numel(params)
    pname = params{pi};
    for li=1:numel(lineshapes)
      line = lineshapes{li};
      X = cell(1,numel(rois));
      for ri=1:numel(rois)
        roi = rois{ri};
        vec = [];
        fHC = fullfile(hcRoot, sprintf('HC_all_%s_%s.mat', line, roi));
        if exist(fHC,'file')
          D = load(fHC); v = gather_roi(D.data, pname); vec = [vec; v(:)];
        end
        fMS = fullfile(msRoot, sprintf('MS_all_%s_%s.mat', line, roi));
        if exist(fMS,'file')
          D = load(fMS); v = gather_roi(D.data, pname); vec = [vec; v(:)];
        end
        vec = vec(isfinite(vec)); if ~includeZeros, vec = vec(vec~=0); end
        X{ri} = vec;
      end
      if all(cellfun(@isempty,X)), continue; end
      k = max(cellfun(@numel,X)); M = NaN(k, numel(rois));
      for r=1:numel(rois), v=X{r}; M(1:numel(v),r)=v; end

      fig = figure('Visible','off','Units','normalized','Position',[0.1 0.1 0.6 0.6], ...
                   'Name', sprintf('%s x %s x All Tissues', pname, line));
      ax = axes(fig); box(ax,'on');
      boxplot(ax, M, 'Labels', rois, 'Symbol','');
      ylabel(ax, pretty_name(pname));
      title(ax, sprintf('%s × %s × All Tissues', pname, line));
      if useBounds && isfield(bounds,pname), ylim(ax, bounds.(pname)); end
      saveas(fig, fullfile(outDir, sprintf('BOX_%sx%sxAllTissues.png', strip_us(pname), line)));
      close(fig);
    end
  end

  % ---------- B) Param-only (12 boxes = lineshape×tissue) ----------
  for pi=1:numel(params)
    pname = params{pi};
    labels = cell(1, numel(lineshapes)*numel(rois));
    X = cell(size(labels)); idx = 1;
    for li=1:numel(lineshapes)
      line = lineshapes{li};
      for ri=1:numel(rois)
        roi = rois{ri};
        vec = [];
        fHC = fullfile(hcRoot, sprintf('HC_all_%s_%s.mat', line, roi));
        if exist(fHC,'file'), D = load(fHC); vec = [vec; gather_roi(D.data, pname)]; end
        fMS = fullfile(msRoot, sprintf('MS_all_%s_%s.mat', line, roi));
        if exist(fMS,'file'), D = load(fMS); vec = [vec; gather_roi(D.data, pname)]; end
        vec = vec(isfinite(vec)); if ~includeZeros, vec = vec(vec~=0); end
        X{idx} = vec; labels{idx} = sprintf('%s-%s', line, roi); idx = idx+1;
      end
    end
    if all(cellfun(@isempty,X)), continue; end
    k = max(cellfun(@numel,X)); M = NaN(k, numel(X));
    for j=1:numel(X), v = X{j}; M(1:numel(v),j)=v; end

    fig = figure('Visible','off','Units','normalized','Position',[0.08 0.08 0.85 0.6], ...
                 'Name', sprintf('%s x All Lineshapes x All Tissues', pname));
    ax = axes(fig); box(ax,'on');
    boxplot(ax, M, 'Labels', labels, 'LabelOrientation','inline', 'Symbol','');
    xtickangle(ax, 30);
    ylabel(ax, pretty_name(pname));
    title(ax, sprintf('%s × All Lineshapes × All Tissues', pname));
    if useBounds && isfield(bounds,pname), ylim(ax, bounds.(pname)); end
    saveas(fig, fullfile(outDir, sprintf('BOX_%sxAllLineshapesxAllTissues.png', strip_us(pname))));
    close(fig);
  end

  % ---------- C) Param + Tissue (3 boxes = SL,L,G) ----------
  for pi=1:numel(params)
    pname = params{pi};
    for ri=1:numel(rois)
      roi = rois{ri};
      X = cell(1,numel(lineshapes));
      for li=1:numel(lineshapes)
        line = lineshapes{li};
        vec = [];
        fHC = fullfile(hcRoot, sprintf('HC_all_%s_%s.mat', line, roi));
        if exist(fHC,'file'), D = load(fHC); vec = [vec; gather_roi(D.data, pname)]; end
        fMS = fullfile(msRoot, sprintf('MS_all_%s_%s.mat', line, roi));
        if exist(fMS,'file'), D = load(fMS); vec = [vec; gather_roi(D.data, pname)]; end
        vec = vec(isfinite(vec)); if ~includeZeros, vec = vec(vec~=0); end
        X{li} = vec;
      end
      if all(cellfun(@isempty,X)), continue; end
      k = max(cellfun(@numel,X)); M = NaN(k, numel(lineshapes));
      for j=1:numel(lineshapes), v=X{j}; M(1:numel(v),j)=v; end

      fig = figure('Visible','off','Units','normalized','Position',[0.1 0.1 0.55 0.6], ...
                   'Name', sprintf('%s x All Lineshapes x %s', pname, roi));
      ax = axes(fig); box(ax,'on');
      boxplot(ax, M, 'Labels', lineshapes, 'Symbol','');
      ylabel(ax, pretty_name(pname));
      title(ax, sprintf('%s × All Lineshapes × %s', pname, roi));
      if useBounds && isfield(bounds,pname), ylim(ax, bounds.(pname)); end
      saveas(fig, fullfile(outDir, sprintf('BOX_%sxAllLineshapesx%s.png', strip_us(pname), roi)));
      close(fig);
    end
  end
end

function simple_boxplots(msRoot,hcRoot,lineshapes,rois,params,bounds,useBounds,outDir)
% Minimal, robust box-and-whisker figures.
% Produces three sets:
%   1) BOXONLY_HCvMS_<line>_<ROI>_<param>.png
%   2) BOXONLY_HC_LINES_<ROI>_<param>.png (across lineshapes)
%   3) BOXONLY_MS_LINES_<ROI>_<param>.png (across lineshapes)

  % ---- (1) HC vs MS per lineshape/ROI/param ----
  for li=1:numel(lineshapes)
    line = lineshapes{li};
    for ri=1:numel(rois)
      roi = rois{ri};
      fHC = fullfile(hcRoot, sprintf('HC_all_%s_%s.mat', line, roi));
      fMS = fullfile(msRoot, sprintf('MS_all_%s_%s.mat', line, roi));
      if ~exist(fHC,'file') || ~exist(fMS,'file'), continue; end
      DHC = load(fHC); DMS = load(fMS);  % expect D*.data

      for pi=1:numel(params)
        pname = params{pi};
        try
          x_HC = gather_roi(DHC.data, pname); x_HC = x_HC(isfinite(x_HC));
          x_MS = gather_roi(DMS.data, pname); x_MS = x_MS(isfinite(x_MS));
        catch, continue; end

        % Build NaN-padded matrix for boxplot
        n1 = numel(x_HC); n2 = numel(x_MS); k = max(n1,n2);
        M = NaN(k,2); M(1:n1,1)=x_HC; M(1:n2,2)=x_MS;

        fig = figure('Visible','off','Units','normalized','Position',[0.1 0.1 0.42 0.6], ...
                     'Name', sprintf('BOX %s — %s (%s) HCvMS', roi, pname, line));
        ax = axes(fig); box(ax,'on');
        boxplot(ax, M, 'Labels',{'HC','MS'}, 'Symbol','');
        ylabel(ax, pname);
        title(ax, sprintf('%s — %s  (%s)', roi, pname, line));
        if useBounds && isfield(bounds,pname)
          ylim(ax, bounds.(pname));
        end
        saveas(fig, fullfile(outDir, sprintf('BOXONLY_HCvMS_%s_%s_%s.png', line, roi, strrep(pname,'_us',''))));
        close(fig);
      end
    end
  end

  % ---- (2,3) Across lineshapes (HC-only and MS-only) ----
  for ri=1:numel(rois)
    roi = rois{ri};
    for pi=1:numel(params)
      pname = params{pi};
      Xhc = cell(1,numel(lineshapes)); Xms = cell(1,numel(lineshapes));
      for lj=1:numel(lineshapes)
        fHC = fullfile(hcRoot, sprintf('HC_all_%s_%s.mat', lineshapes{lj}, roi));
        fMS = fullfile(msRoot, sprintf('MS_all_%s_%s.mat', lineshapes{lj}, roi));
        if exist(fHC,'file')
          try
            v = gather_roi(load(fHC).data, pname); v = v(isfinite(v)); Xhc{lj}=v;
          catch, Xhc{lj} = []; end
        end
        if exist(fMS,'file')
          try
            v = gather_roi(load(fMS).data, pname); v = v(isfinite(v)); Xms{lj}=v;
          catch, Xms{lj} = []; end
        end
      end
      % HC only
      if any(~cellfun(@isempty,Xhc))
        k = max(cellfun(@numel,Xhc)); Mhc = NaN(k,numel(lineshapes));
        for jj=1:numel(lineshapes), v=Xhc{jj}; Mhc(1:numel(v),jj)=v; end
        fig = figure('Visible','off','Units','normalized','Position',[0.1 0.1 0.5 0.6], ...
                     'Name', sprintf('BOX HC lines — %s %s', roi, pname));
        ax = axes(fig); box(ax,'on');
        boxplot(ax, Mhc, 'Labels', lineshapes, 'Symbol','');
        ylabel(ax, pname);
        title(ax, sprintf('%s — %s (HC across lineshapes)', roi, pname));
        if useBounds && isfield(bounds,pname), ylim(ax, bounds.(pname)); end
        saveas(fig, fullfile(outDir, sprintf('BOXONLY_HC_LINES_%s_%s.png', roi, strrep(pname,'_us',''))));
        close(fig);
      end
      % MS only
      if any(~cellfun(@isempty,Xms))
        k = max(cellfun(@numel,Xms)); Mms = NaN(k,numel(lineshapes));
        for jj=1:numel(lineshapes), v=Xms{jj}; Mms(1:numel(v),jj)=v; end
        fig = figure('Visible','off','Units','normalized','Position',[0.1 0.1 0.5 0.6], ...
                     'Name', sprintf('BOX MS lines — %s %s', roi, pname));
        ax = axes(fig); box(ax,'on');
        boxplot(ax, Mms, 'Labels', lineshapes, 'Symbol','');
        ylabel(ax, pname);
        title(ax, sprintf('%s — %s (MS across lineshapes)', roi, pname));
        if useBounds && isfield(bounds,pname), ylim(ax, bounds.(pname)); end
        saveas(fig, fullfile(outDir, sprintf('BOXONLY_MS_LINES_%s_%s.png', roi, strrep(pname,'_us',''))));
        close(fig);
      end
    end
  end
end

function v = gather_all_tissue(M, li, pname)
% M: [nx,ny,nLines,nFiles,nMaps]. Compute derived params and flatten.
  switch pname
    case 'PSR'
      v = M(:,:,li,: ,1);
    case 'kba'
      v = M(:,:,li,: ,2);
    case 'T2R1'
      v = M(:,:,li,: ,3) .* M(:,:,li,: ,5);
    case 'T2b_us'
      v = M(:,:,li,: ,4) * 1e6;
    otherwise
      error('Unknown param %s',pname);
  end
  v = v(:);
end

function v = gather_roi(D, pname)
% D: [nx,ny,nFiles,nMaps] from per-ROI stack. Compute derived params & flatten.
  switch pname
    case 'PSR'
      v = D(:,:,:,1);
    case 'kba'
      v = D(:,:,:,2);
    case 'T2R1'
      v = D(:,:,:,3) .* D(:,:,:,5);
    case 'T2b_us'
      v = D(:,:,:,4) * 1e6;
    otherwise
      error('Unknown param %s',pname);
  end
  v = v(:);
end

function x = clean_and_trim(x, pname, useBounds, bounds, trimMethod, madK)
  % drop non-finite and zeros
  x = x(isfinite(x));
  x = x(x ~= 0);
  % optional bounds
  if useBounds
    b = bounds.(pname);
    x = x(x >= b(1) & x <= b(2));
  end
  if isempty(x), return; end
  % optional trimming
  switch lower(trimMethod)
    case 'none'
      % do nothing
    case 'tukey'
      q = quantile(x,[0.25 0.75]); I = q(2)-q(1);
      lo = q(1) - 1.5*I; hi = q(2) + 1.5*I;
      x = x(x>=lo & x<=hi);
    case 'mad'
      medx = median(x);
      MAD  = median(abs(x - medx));
      if MAD > 0
        x = x(abs(x - medx) <= madK*MAD);
      end
    otherwise
      error('Unknown trimMethod %s', trimMethod);
  end
end

function S = summarize_vec(x)
  if isempty(x)
    S = struct('n',0,'median',NaN,'iqr',[NaN NaN],'mad',NaN,'mean',NaN,'std',NaN);
    return;
  end
  S.n      = numel(x);
  S.median = median(x);
  q        = quantile(x,[0.25 0.75]);
  S.iqr    = q;
  S.mad    = median(abs(x - S.median));
  S.mean   = mean(x);
  S.std    = std(x,0);
end

function [pval, dCliff, medDiff] = compare_groups(xHC, xMS)
  medDiff = median(xMS) - median(xHC);
  pval = NaN; dCliff = NaN;
  try
    [pval,~,stats] = ranksum(xHC, xMS);  % requires Statistics Toolbox
    n = numel(xHC); m = numel(xMS);
    if isfield(stats,'ranksum') && n>0 && m>0
      U = stats.ranksum - n*(n+1)/2;  % Mann-Whitney U for HC
      % Convert to Cliff's delta: delta = 2*U/(n*m) - 1, but U must be for MS over HC.
      % We computed U for HC; for symmetry, use U_MS = n*m - U.
      U_ms = n*m - U;
      dCliff = (2*U_ms/(n*m)) - 1;
    end
  catch
    % leave as NaN if ranksum unavailable
  end
end

function edges = choose_edges(x1, x2, pname, useFixedBW, fixedBW, useBounds, bounds, minNBins, maxNBins)
  x = [x1(:); x2(:)];
  if isempty(x)
    edges = linspace(0,1,16); return;
  end
  if useBounds
    b = bounds.(pname);
    lo = b(1); hi = b(2);
  else
    lo = min(x); hi = max(x);
    if lo==hi, lo = lo - eps; hi = hi + eps; end
  end
  if useFixedBW
    bw = fixedBW.(pname);
    loE = floor(lo/bw)*bw;
    hiE = ceil(hi/bw)*bw;
    edges = loE:bw:hiE;
    if numel(edges)<8
      edges = linspace(lo,hi,16);
    end
  else
    % Freedman–Diaconis with caps
    q = quantile(x,[0.25 0.75]); I = q(2)-q(1);
    if I==0
      edges = linspace(lo,hi,16);
    else
      bwFD = 2*I/(numel(x)^(1/3));
      if bwFD<=0 || ~isfinite(bwFD)
        edges = linspace(lo,hi,16);
      else
        nBins = max(min(round((hi-lo)/bwFD), maxNBins), minNBins);
        edges = linspace(lo,hi,nBins+1);
      end
    end
  end
end


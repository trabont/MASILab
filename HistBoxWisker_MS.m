%% HistBoxWisker_MS.m
% 1) HISTOGRAMS + KDE OVERLAY with Buffered Truncation & Data‐Driven Bins
% 2) COMBINED BOX & WHISKER PER PARAMETER
clear; close all; clc;

% === SELECT FOLDERS ===
SUBJ_combined = fullfile(pwd,'Processed','FullFit_MS');
save_dir    = fullfile(pwd,'Processed','FullFit_MS');

% === CONFIGURATION ===
lineshapes = {'SL','L','G'};
tissues    = {'WM','GM','CSF','LYMPH'};
params     = {'kba','PSR','T2R1','T2b'};
groups     = {'SUBJ'};
groupDirs  = {SUBJ_combined};
groupCols  = struct('SUBJ',[0 0 1]);

% Fit boundaries for buffered truncation
initBounds = struct( ...
  'PSR',  [1e-3, 0.6], ...
  'kba',  [1e-3, 30], ...
  'T2R1', [0,     0.1], ...
  'T2b',  [5,   100] ...
);

bufFrac = 0.025;  % 2.5%

% mapping from param → combined‐stack map index
mapNamesComb = {'PSR_map','kba_map','T2a_map','T2b_map','R1obs_map'};

%% 1) HISTOGRAMS + KDE
for ls = 1:numel(lineshapes)
  shape = lineshapes{ls};
  for ti = 1:numel(tissues)
    tissue = tissues{ti};
    for pi = 1:numel(params)
      param  = params{pi};
      bounds = initBounds.(param);

      % load & convert raw values per group from combined stacks
      raw.SUBJ = [];
      for g = 1
        grp   = groups{g};
        fn    = fullfile(groupDirs{g}, sprintf('all_%s_%s.mat',shape,tissue));
        S     = load(fn);
        D4    = S.data;           % X×Y×Nfiles×5
        switch param
          case 'kba'
            idx = find(strcmp(mapNamesComb,'kba_map'));
            tmp = D4(:,:,:,idx);
          case 'PSR'
            idx = find(strcmp(mapNamesComb,'PSR_map'));
            tmp = D4(:,:,:,idx);
          case 'T2R1'
            i1 = find(strcmp(mapNamesComb,'T2a_map'));
            i2 = find(strcmp(mapNamesComb,'R1obs_map'));
            tmp = squeeze(D4(:,:,:,i1)) .* squeeze(D4(:,:,:,i2));
          case 'T2b'
            idx = find(strcmp(mapNamesComb,'T2b_map'));
            tmp = D4(:,:,:,idx).*1e6;
        end
        tmp = tmp(:);
        raw.(grp) = tmp(tmp>0 & ~isnan(tmp));
      end

      rawAll = [raw.SUBJ];
      if isempty(rawAll), continue; end

      % buffered truncation
      rngVal  = bounds(2) - bounds(1);
      lowEdge = bounds(1) + bufFrac*rngVal;
      highEdge= bounds(2) - bufFrac*rngVal;
      Up.SUBJ   = raw.SUBJ(raw.SUBJ>=lowEdge & raw.SUBJ<=highEdge);
      UpAll   = [Up.SUBJ];
      if isempty(UpAll)
        Up.SUBJ = raw.SUBJ; UpAll = rawAll;
        fprintf('Buffer removed all UpData for %s-%s-%s; using raw data.\n', ...
                shape,tissue,param);
      end

      % FD bin width + fallback
      n    = numel(UpAll);
      IQRv = iqr(UpAll);
      bw   = 2*IQRv / n^(1/3);
      if bw <= 0
        bw = (max(UpAll)-min(UpAll))/50;
        fprintf('FD bin width <=0 for %s-%s-%s; using fallback bw.\n',shape,tissue,param);
      end
      edges = min(UpAll):bw:max(UpAll);

      % too few bins?
      numBins = numel(edges)-1;
      if numBins < 8
        fprintf('Only %d bins for %s-%s-%s; switching to 15 uniform bins.\n', ...
                numBins,shape,tissue,param);
        edges = linspace(min(UpAll), max(UpAll), 15);
        bw    = edges(2)-edges(1);
      end

      ctrs = edges(1:end-1) + bw/2;

      % percent per bin
      N_SUBJ   = histcounts(Up.SUBJ, edges);
      pct_SUBJ = N_SUBJ / sum(N_SUBJ) * 100;

      % KDE overlay
      xi  = linspace(edges(1), edges(end), 200);
      if ~isempty(Up.SUBJ), fSUBJ = ksdensity(Up.SUBJ, xi); else fSUBJ = zeros(size(xi)); end
      ySUBJ = fSUBJ * bw * 100;

      % plot
      fig = figure('Visible','off','Position',[100 100 800 600]);
      ax  = axes(fig); hold(ax,'on');
      bar(ax, ctrs, pct_SUBJ, 1.0, 'FaceColor',groupCols.SUBJ,'FaceAlpha',0.3,'EdgeColor','none');
      plot(ax, xi, ySUBJ, 'Color',groupCols.SUBJ,'LineWidth',1.5);
      title(ax, sprintf('%s — %s (%s)', tissue, param, shape), 'Interpreter','none');
      xlabel(ax, param); ylabel(ax,'Voxel % (sum=100)');
      xlim(ax, [edges(1), edges(end)]);
      legend(ax,{'SUBJ'}, 'Location','northeast');
      saveas(fig, fullfile(save_dir, sprintf('%s_%s_%s_Hist.png',shape,tissue,param)));
      close(fig);
    end
  end
end

%% 2) COMBINED BOX & WHISKER PER PARAMETER
for pi = 1:numel(params)
  param  = params{pi};
  bounds = initBounds.(param);
  allData = [];
  allCats = string([]);
  allGroup= string([]);

  for ls = 1:numel(lineshapes)
    shape = lineshapes{ls};
    for ti = 1:numel(tissues)
      tissue = tissues{ti};

      % load & truncate per group
      for g = 1
        grp = groups{g};
        fn  = fullfile(groupDirs{g}, sprintf('all_%s_%s.mat',shape,tissue));
        S   = load(fn);
        D4  = S.data;
        switch param
          case 'kba'
            idx = find(strcmp(mapNamesComb,'kba_map'));
            tmp = D4(:,:,:,idx);
          case 'PSR'
            idx = find(strcmp(mapNamesComb,'PSR_map'));
            tmp = D4(:,:,:,idx);
          case 'T2R1'
            i1  = find(strcmp(mapNamesComb,'T2a_map'));
            i2  = find(strcmp(mapNamesComb,'R1obs_map'));
            tmp = squeeze(D4(:,:,:,i1)) .* squeeze(D4(:,:,:,i2));
          case 'T2b'
            idx = find(strcmp(mapNamesComb,'T2b_map'));
            tmp = D4(:,:,:,idx)*1e6;
        end
        vals = tmp(:);
        vals = vals(vals>0 & ~isnan(vals));

        % buffered truncation
        rngV  = bounds(2)-bounds(1);
        lowE  = bounds(1)+bufFrac*rngV;
        highE = bounds(2)-bufFrac*rngV;
        Up    = vals(vals>=lowE & vals<=highE);
        if isempty(Up)
          Up = vals;
        end

        % accumulate
        nUp = numel(Up);
        allData  = [allData;  Up];
        allCats  = [allCats;  repmat(string(sprintf('%s_/%s',shape,tissue)), nUp,1)];
        allGroup = [allGroup; repmat(grp, nUp,1)];
      end
    end
  end

  % define category order
  catsOrder = { ...
    'SL_/WM','L_/WM','G_/WM', ...
    'SL_/GM','L_/GM','G_/GM', ...
    'SL_/CSF','L_/CSF','G_/CSF', ...
    'SL_/LYMPH','L_/LYMPH','G_/LYMPH' ...
  };
  C = categorical(allCats, catsOrder, 'Ordinal',true);

  % colors per group
  cols = zeros(numel(allGroup), 3);
  for i = 1:numel(allGroup)
    cols(i,:) = groupCols.(allGroup(i));
  end

  % plot
  fig = figure('Visible','off','Units','normalized','Position',[.2 .2 .6 .6]);
  ax  = axes(fig); hold(ax,'on');
  boxchart(ax, C, allData,'BoxWidth',0.6);
  xnum = double(C);
  xjit = xnum + (rand(size(xnum))-0.5)*0.2;
  scatter(ax, xjit, allData, 8, cols, 'filled','MarkerFaceAlpha',0.3);

  % medians
  for k = 1:numel(catsOrder)
    mask = C==catsOrder{k};
    if any(mask)
      mval = median(allData(mask));
      scatter(ax, k, mval, 20,'k','d','filled');
      text(ax, k, mval, sprintf('%.3f',mval), ...
           'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',16,'FontWeight','bold');
    end
  end

  title(ax, sprintf('Box & Whisker + Points: %s',param), 'Interpreter','none');
  xlabel(ax,'Lineshape_/Tissue'); ylabel(ax,param);
  xtickangle(ax,45);
  saveas(fig, fullfile(save_dir, sprintf('Box_%s.png',param)));
  close(fig);
end

disp('All histograms and boxplots complete.');

%% 3) OVERLAID HISTOGRAMS + KDE ACROSS TISSUES PER LINESHAPE & PARAMETER
tissueCols = struct( ...
  'WM',   [1 0 0], ...    % red
  'GM',   [0 0 1], ...    % blue
  'CSF',  [0 1 1], ...    % cyan
  'LYMPH',[1 0 1] ...     % magenta
);

for ls = 1:numel(lineshapes)
  shape = lineshapes{ls};
  for pi = 1:numel(params)
    param  = params{pi};
    bounds = initBounds.(param);

    % gather & truncate data for each tissue
    UpAll = [];
    for ti = 1:numel(tissues)
      tissue = tissues{ti};
      fn     = fullfile(SUBJ_combined, sprintf('all_%s_%s.mat',shape,tissue));
      S      = load(fn);
      D4     = S.data;
      switch param
        case 'kba'
          tmp = D4(:,:,:,strcmp(mapNamesComb,'kba_map'));
        case 'PSR'
          tmp = D4(:,:,:,strcmp(mapNamesComb,'PSR_map'));
        case 'T2R1'
          tmp = squeeze(D4(:,:,:,strcmp(mapNamesComb,'T2a_map'))) .* ...
                squeeze(D4(:,:,:,strcmp(mapNamesComb,'R1obs_map')));
        case 'T2b'
          tmp = D4(:,:,:,strcmp(mapNamesComb,'T2b_map')) * 1e6;
      end
      v = tmp(:);
      v = v(v>0 & ~isnan(v));
      rngV  = bounds(2)-bounds(1);
      lowE  = bounds(1)+bufFrac*rngV;
      highE = bounds(2)-bufFrac*rngV;
      Up.(tissue) = v(v>=lowE & v<=highE);
      if isempty(Up.(tissue)), Up.(tissue)=v; end
      UpAll = [UpAll; Up.(tissue)];
    end

    if isempty(UpAll)
      warning('No data for %s-%s, skipping.', shape, param);
      continue;
    end

    % compute shared bins
    n    = numel(UpAll);
    IQRv = iqr(UpAll);
    bw   = 2*IQRv/n^(1/3);
    if bw<=0, bw=(max(UpAll)-min(UpAll))/50; end
    edges = min(UpAll):bw:max(UpAll);
    if numel(edges)-1<8
      edges = linspace(min(UpAll),max(UpAll),15);
      bw    = edges(2)-edges(1);
    end
    ctrs = edges(1:end-1)+bw/2;

    % plot
    fig = figure('Visible','off','Position',[100 100 800 600]);
    ax  = axes(fig); hold(ax,'on');
    xi  = linspace(edges(1),edges(end),200);
    for ti = 1:numel(tissues)
      tissue = tissues{ti};
      data   = Up.(tissue);
      if isempty(data), continue; end

      % histogram bar
      N   = histcounts(data, edges);
      pct = N/sum(N)*100;
      bar(ax, ctrs, pct, 1.0, ...
          'FaceColor', tissueCols.(tissue), ...
          'FaceAlpha', 0.3, 'EdgeColor','none');

      % KDE curve
      f   = ksdensity(data, xi);
      y   = f * bw * 100;
      plot(ax, xi, y, 'LineWidth',1.5, ...
           'Color', tissueCols.(tissue));
    end

    title(ax, sprintf('%s — %s (All Tissues)', param, shape), 'Interpreter','none');
    xlabel(ax, param); ylabel(ax, 'Voxel % (sum=100)');
    xlim(ax, [edges(1) edges(end)]);
    legend(ax, tissues, 'Location','northeast');
    saveas(fig, fullfile(save_dir, sprintf('%s_%s_AllTissues_Hist_KDE.png', shape, param)));
    close(fig);
  end
end

disp('Section 3: Overlaid histograms + KDE complete.');

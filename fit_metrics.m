%% fit_metrics.m
% Builds HC medians (excluding zeros) directly from per-ROI stacks.
% Exports:
%   - Processed/FullFit_Summary/median_values.mat  (medians [Tissue x Param x Lineshape])
%   - Processed/FullFit_Summary/median_values.xlsx (ONE sheet: Tissue | Lineshape | Parameter | Median)

clear; clc;

% ---- I/O ----
hcRoot = fullfile(pwd,'Processed','FullFit_HC');
outDir = fullfile(pwd,'Processed','FullFit_Summary');
if ~exist(outDir,'dir'), mkdir(outDir); end

% ---- Labels ----
lineshapes = {'SL','L','G'};
rois       = {'GM','WM','CSF','LYMPH'};

% Parameter names (lowercase in Excel/struct)
params     = {'psr','kba','t2a','t2b','r1obs','t2r1'};

nR = numel(rois); nP = numel(params); nL = numel(lineshapes);
medians = NaN(nR, nP, nL);   % [Tissue x Param x Lineshape]

% Rows for Excel: Tissue | Lineshape | Parameter | Median
xls_T = strings(0,1);
xls_L = strings(0,1);
xls_P = strings(0,1);
xls_M = [];

for li = 1:nL
  line = lineshapes{li};
  for ri = 1:nR
    roi = rois{ri};
    fHC = fullfile(hcRoot, sprintf('HC_all_%s_%s.mat', line, roi));
    if ~exist(fHC,'file')
      warning('Missing %s', fHC);
      continue;
    end
    S = load(fHC);
    if ~isfield(S,'data')
      warning('No ''data'' field in %s', fHC);
      continue;
    end
    D = S.data;                       % [nx,ny,nFiles,nMaps]
    mapNames = []; if isfield(S,'mapNames'), mapNames = S.mapNames; end
    idx = default_indices();          % fallback order
    if ~isempty(mapNames), idx = resolve_indices(mapNames, idx); end

    for pi = 1:nP
      p = params{pi};
      v = get_param_vec(D, p, idx);   % flatten vector for this param
      v = v(isfinite(v) & v~=0);      % EXCLUDE zeros + non-finite
      med = NaN; if ~isempty(v), med = median(v); end

      medians(ri,pi,li) = med;

      xls_T(end+1,1) = roi;           
      xls_L(end+1,1) = line;
      xls_P(end+1,1) = p;
      xls_M(end+1,1) = med;
    end
  end
end

% ---- Save .mat ----
save(fullfile(outDir,'median_values.mat'), 'medians', 'rois', 'params', 'lineshapes');

% ---- Save ONE-sheet Excel (hard overwrite to avoid stale columns) ----
xlsx = fullfile(outDir,'median_values.xlsx');
if exist(xlsx,'file')
  try, delete(xlsx); catch, error('Close %s in Excel and rerun.', xlsx); end
end
T = table(xls_T, xls_L, xls_P, xls_M, ...
  'VariableNames', {'Tissue','Lineshape','Parameter','Median'});
writetable(T, xlsx);

fprintf('Wrote:\n  %s\n  %s\n', ...
  fullfile(outDir,'median_values.mat'), xlsx);

% ================= helpers =================
function idx = default_indices()
% Default order produced by combine.m:
% {'PSR_map','kba_map','T2a_map','T2b_map','R1obs_map'}
  idx.psr   = 1;
  idx.kba   = 2;
  idx.t2a   = 3;
  idx.t2b   = 4;
  idx.r1obs = 5;
end

function idx = resolve_indices(mapNames, idx)
% Tolerant name matching if map order differs
  for i = 1:numel(mapNames)
    m = lower(string(mapNames{i}));
    if contains(m,'psr'),      idx.psr   = i; end
    if contains(m,'kba'),      idx.kba   = i; end
    if contains(m,'t2a'),      idx.t2a   = i; end
    if contains(m,'t2b'),      idx.t2b   = i; end
    if contains(m,'r1obs')||strcmp(m,'r1')||contains(m,'r_1')
      idx.r1obs = i;
    end
  end
end

function v = get_param_vec(D, p, idx)
% D: [nx,ny,nFiles,nMaps]; p is one of {'psr','kba','t2a','t2b','r1obs','t2r1'}
  switch p
    case 'psr',    v = D(:,:,:,idx.psr);
    case 'kba',    v = D(:,:,:,idx.kba);
    case 't2a',    v = D(:,:,:,idx.t2a);
    case 't2b',    v = D(:,:,:,idx.t2b);
    case 'r1obs',  v = D(:,:,:,idx.r1obs);
    case 't2r1',   v = D(:,:,:,idx.t2a) .* D(:,:,:,idx.r1obs);
    otherwise, error('Unknown parameter %s', p);
  end
  v = v(:);
end

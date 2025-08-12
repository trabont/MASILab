% Dice + Centroid + Median ROI Excel Analysis (no voxelwise loops)
function med_dice_centr(baseDir, saveDir, fileID)
% Running Requirements:
%   - functions/Analysis_Yarnykh_Full_Fit.m on path
%   - baseDir/<ID>/ contains "*_registered_*.nii.gz"
%   - saveDir exists or will be created
%
% What this does:
%   • Loads MT, MFA, B0, B1 from both SCT(_1) and ANTS where relevant
%   • Loads masks: WM, GM, CSF (registered_*.nii.gz), LYMPH (registered_*_ANTS.nii.gz)
%   • For each slice that has LYMPH (>10 voxels), and for each ROI:
%       - Selects the correct registration source (LYMPH→ANTS, others→SCT)
%       - Builds median MT signals (normalize per voxel first, then take median per offset)
%       - Uses median B0, B1, and median Ernst (MFA) signals
%       - Runs Full Fit once per lineshape (SL, L, G)
%   • Writes one Excel file: <saveDir>/<ID>_medianROI.xlsx with columns:
%       ID, Slice, ROI, Reg, Lineshape, PSR, kba, T2a, T2b, R1obs
%
% Notes:
%   - If no LYMPH slices are found, it returns (keeps your previous behavior).
%   - Dice/Centroid not computed here (name preserved for continuity). Add later if needed.

addpath("functions/");

if ~isfolder(saveDir), mkdir(saveDir); end
idstr = char(string(fileID));

% --- filename helpers (use %s so string or numeric IDs both work) ---
fn_1 = @(key) fullfile(baseDir, sprintf('%s/%s_registered_%s_1.nii.gz',   idstr, idstr, key));   % SCT
fn_2 = @(key) fullfile(baseDir, sprintf('%s/%s_registered_%s_ANTS.nii.gz', idstr, idstr, key)); % ANTS
fn_3 = @(key) fullfile(baseDir, sprintf('%s/%s_registered_%s.nii.gz',     idstr, idstr, key));  % masks (WM/GM/CSF)

% --- Load data (SCT & ANTS) ---
data.MT_SCT    = double(niftiread(fn_1('MT')));    % X×Y×Z×16
data.B0_SCT    = double(niftiread(fn_1('B0')));    % X×Y×Z
data.B1_SCT    = double(niftiread(fn_1('B1')));    % X×Y×Z
data.MFA_SCT   = double(niftiread(fn_1('MFA')));   % X×Y×Z×MFA
data.T2AX_SCT  = double(niftiread(fn_1('T2AX')));   % X×Y×Z
data.mFFE      = double(niftiread(fn_1('mFFE')));

data.MT_ANTS   = double(niftiread(fn_2('MT')));    % X×Y×Z×16
data.B0_ANTS   = double(niftiread(fn_2('B0')));
data.B1_ANTS   = double(niftiread(fn_2('B1')));
data.MFA_ANTS  = double(niftiread(fn_2('MFA')));
data.T2AX_ANTS = double(niftiread(fn_2('T2AX')));   % X×Y×Z

% --- ROI masks ---
masks.LYMPH = logical(niftiread(fn_2('LYMPH'))); % ANTS
masks.WM    = logical(niftiread(fn_3('WM')));    % registered_*.nii.gz
masks.GM    = logical(niftiread(fn_3('GM')));
masks.CSF   = logical(niftiread(fn_3('CSF')));
roiNames = {'LYMPH','WM','GM','CSF'};

% --- Spinal Cord masks ---
segFiles = struct( ...
  'MT0',  fullfile(baseDir, sprintf('%s/%s_MT_dyn_0000_reg_sc.nii.gz', idstr, idstr)), ...
  'mFFE',  fullfile(baseDir, sprintf('%s/%s_registered_mFFE_1_sc.nii.gz', idstr, idstr)), ...
  'MT8',  fullfile(baseDir, sprintf('%s/%s_MT_dyn_0008_reg_sc.nii.gz', idstr, idstr)), ...
  'MFA0', fullfile(baseDir, sprintf('%s/%s_MFA_dyn_0000_in_MT_sc.nii.gz', idstr, idstr)), ...
  'T1',   fullfile(baseDir, sprintf('%s/%s_registered_T1_1_sc.nii.gz', idstr, idstr)), ...
  'T2',   fullfile(baseDir, sprintf('%s/%s_registered_T2AX_1_sc.nii.gz', idstr, idstr)) );

seg = struct();
names = fieldnames(segFiles);
for ii = 1:numel(names)
    f = segFiles.(names{ii});
    if exist(f,'file')
        seg.(names{ii}) = logical(niftiread(f));
    else
        seg.(names{ii}) = false( size(masks.LYMPH) ); % safe empty
        warning('Missing SC mask: %s', f);
    end
end

seg_Name = {'MT0','MT8','MFA0','T1','T2'};
seg_CT = numel(seg_Name);

% --- determine slices to process (based on LYMPH as you had) ---
sliceSum = squeeze(sum(sum(masks.LYMPH,1),2));
slices   = find(sliceSum > 10);
if isempty(slices)
  warning('Subject %s: no LYMPH mask → skipping', fileID);
  return;
end

% --- qMT acquisition parameters ---
rf_offset = [1000,1500,2000,2500,8000,16000,32000,100000];  % 8
deg       = [360, 820];                                     % 2 blocks
scan      = numel(rf_offset);                               % 8

BaseParms.T1flip   = 30;
BaseParms.MFA      = 6;       % uses first 6 MFA echoes
BaseParms.T1TR     = 50;
BaseParms.deltaMT  = rf_offset;
BaseParms.pwMT     = [20e-3, 20e-3];
BaseParms.MT_flip  = deg;
BaseParms.qMTflip  = 6;
BaseParms.TR       = 50;

lines   = {'SL','L','G'};
shapes  = {'super-lorentzian','lorentzian','gaussian'};

% --- results accumulator for Excel ---
acc = cell(0,10);  % ID,Slice,ROI,Reg,Lineshape,PSR,kba,T2a,T2b,R1obs
bcc = cell(0,5);  % ID,Slice,Seg,DICE,Centroid Diff
dcc = cell(0,5); %ID, x-pixels, y-pixels, z-dimension, t-dimension

% --- main loop over slices ---
for zz = slices(:)'
  % pre-slice views
  MT_SCT  = squeeze(data.MT_SCT(:,:,zz,:));   % X×Y×16
  MT_ANTS = squeeze(data.MT_ANTS(:,:,zz,:));
  B0_SCT  = data.B0_SCT(:,:,zz);  
  B1_SCT = data.B1_SCT(:,:,zz);
  B0_ANTS = data.B0_ANTS(:,:,zz); 
  B1_ANTS = data.B1_ANTS(:,:,zz);
  MFA_SCT = squeeze(data.MFA_SCT(:,:,zz, 1:BaseParms.MFA));   % X×Y×M
  MFA_ANTS= squeeze(data.MFA_ANTS(:,:,zz,1:BaseParms.MFA));

  for r = 1:numel(roiNames)
    roi = roiNames{r};
    mask2D = masks.(roi)(:,:,zz);

    if ~any(mask2D(:)), continue; end

    % pick source volume set: LYMPH → ANTS; others → SCT (matches your vz logic)
    if strcmp(roi,'LYMPH')
      MT   = MT_ANTS;   B0 = B0_ANTS;   B1 = B1_ANTS;   MF = MFA_ANTS;   regUsed = "ANTS";
    else
      MT   = MT_SCT;    B0 = B0_SCT;    B1 = B1_SCT;    MF = MFA_SCT;    regUsed = "SCT";
    end

    % --- gather voxel stacks within ROI ---
    % MT: [Nvox × 16], split into two blocks of 8 (per your acquisition)
    voxIdx = find(mask2D);
    [X,Y,~] = size(MT); %#ok<ASGLU>

    % extract MT per voxel (normalize per voxel within each block, then median per offset)
    MTv = reshape(MT, [], size(MT,3));         % [XY × 16]
    MTv = MTv(voxIdx, :);                      % [Nvox × 16]
    M1  = MTv(:, 1:scan);                      % [Nvox × 8]
    M2  = MTv(:, scan+1:end);                  % [Nvox × 8]

    % normalize per voxel by last point (guard zeros/NaNs)
    den1 = M1(:,end); den1(~isfinite(den1) | den1==0) = NaN;
    den2 = M2(:,end); den2(~isfinite(den2) | den2==0) = NaN;
    M1n  = M1 ./ den1;
    M2n  = M2 ./ den2;

    M1_med = median(M1n, 1, 'omitnan');
    M2_med = median(M2n, 1, 'omitnan');
    if all(~isfinite(M1_med)) || all(~isfinite(M2_med)), continue; end
    BaseParms.M = [M1_med;M2_med];

    B0_med = median(B0(voxIdx), 'omitnan');
    BaseParms.B0 = B0_med;
    B1_med = median(B1(voxIdx), 'omitnan');
    BaseParms.B1 = B1_med;

    mfa_use = min(BaseParms.MFA, size(MF,3));
    MFv = reshape(MF, [], size(MF,3));
    MFv = MFv(voxIdx, 1:mfa_use);
    Ernst_med = median(MFv, 1, 'omitnan');
    BaseParms.Ernst=Ernst_med;

    % run three lineshapes
    for L = 1:3
      BaseParms.lineshape = shapes{L};

      % Full Fit on ROI-median signals
      [PSR, kba, T2a, T2b, R1obs] = Analysis_Yarnykh_Full_Fit(BaseParms, shapes{L});

        corrB1  = BaseParms.B1/100;
        [B1MT,tMT] = philipsRFpulse_FA(BaseParms.MT_flip, BaseParms.pwMT, 'am_sg_100_100_0');
        B1eMT   = CWEqMTPulse(B1MT*corrB1, tMT, BaseParms.pwMT);
        thetaEX = ([BaseParms.qMTflip BaseParms.qMTflip]*pi/180).*corrB1;
        ts = 1e-3;
        TRs = [BaseParms.TR*1e-3, BaseParms.TR*1e-3];    % seconds
        delta = BaseParms.deltaMT + BaseParms.B0;        % Hz

        M0 = [1, PSR]';
        R1 = [R1obs, R1obs]';
        T2 = [T2a, T2b]';

        [~,Mzn,~] = yarnykh_pulseMT(M0,R1,T2,TRs,kba,BaseParms.pwMT,ts,thetaEX,delta,B1eMT,shapes{L});
        model_b1 = Mzn(1,:);  
        model_b2 = Mzn(2,:);

        fig = figure('Visible','off'); set(fig,'Position',[200 200 760 480]);
        semilogx(rf_offset, M1_med, 'o-'); hold on;
        semilogx(rf_offset, M2_med, 'o-');
        semilogx(rf_offset, model_b1, '-');
        semilogx(rf_offset, model_b2, '-');
        xlabel('RF offset (Hz)'); ylabel('Normalized signal');
        ttl = sprintf('%s | Slice %d | %s (%s) | %s', idstr, zz, roi, regUsed, lines{L});
        title(ttl);
        lg = {'360° data','820° data','360° Fit','820° Fit'};
        legend(lg,'Location','best');
        saveDirID = fullfile(saveDir,idstr);
        outpng = fullfile(saveDirID, sprintf('%s_slice%02d_%s_%s_FitZ.png', idstr, zz, roi, lines{L}));
        exportgraphics(fig, outpng); close(fig);

      % Append one row
      acc(end+1,:) = {fileID, zz, roi, char(regUsed), lines{L}, PSR, kba, T2a, T2b, R1obs}; %#ok<AGROW>
    end
  end % ROI

  for ii = 1:seg_CT
    % Compute dice comparison and centroid
    name = seg_Name{ii};
    A = seg.mFFE(:,:,zz);
    B = seg.(name)(:,:,zz);
    if ~isequal(size(A), size(B))
        warning('Size mismatch mFFE vs %s on slice %d (skipping)', name, zz);
        continue;
    end
    A = logical(A);
    B = logical(B);
    intersection = sum(A(:) & B(:));
    d = (2 * intersection) / (sum(A(:)) + sum(B(:))); % DICE
    
    [rowA, colA] = find(A~=0);
    if isempty(colA)
        cxA = NaN; cyA = NaN;
    else
        cxA = mean(colA);  % x = column index
        cyA = mean(rowA);  % y = row index
    end

    [rowB, colB] = find(B~=0);
    if isempty(colB)
        cxB = NaN; cyB = NaN;
    else
        cxB = mean(colB);  % x = column index
        cyB = mean(rowB);  % y = row index
    end

    if any(isnan([cxA,cyA,cxB,cyB]))
        cdist = NaN;
    else
        cdist = hypot(cxA-cxB, cyA-cyB); % pixels
    end
    bcc(end+1,:) = {fileID, zz, name, d, cdist};
  end

end % slice


% --- Write / Append to Excel by lineshape + QC ---
if ~isempty(acc)
  T = cell2table(acc, 'VariableNames', ...
      {'ID','Slice','ROI','Reg','Lineshape','PSR','kba','T2a','T2b','R1obs'});
  xlsx = fullfile(saveDir, 'FF_MED_DICE_CENTROID.xlsx');

  lsheets = {'SL','L','G'};
  for i = 1:numel(lsheets)
      Ti = T(strcmp(T.Lineshape, lsheets{i}), :);
      if isempty(Ti), continue; end
      if exist(xlsx,'file') && any(strcmpi(sheetnames(xlsx), lsheets{i}))
          % append rows (no headers)
          writetable(Ti, xlsx, 'Sheet', lsheets{i}, ...
              'WriteMode','append', 'WriteVariableNames', false);
      else
          % first time: write with headers
          writetable(Ti, xlsx, 'Sheet', lsheets{i}, 'FileType','spreadsheet');
      end
  end

  if exist('bcc','var') && ~isempty(bcc)
      QC = cell2table(bcc, 'VariableNames', {'ID','Slice','Seg','Dice','CentroidDiff_px'});
      if exist(xlsx,'file') && any(strcmpi(sheetnames(xlsx), 'DICE_CENTROID'))
          writetable(QC, xlsx, 'Sheet','DICE_CENTROID', ...
              'WriteMode','append','WriteVariableNames', false);
      else
          writetable(QC, xlsx, 'Sheet','DICE_CENTROID','FileType','spreadsheet');
      end
  end

  fprintf('Updated %s (sheets: %s + DICE_CENTROID)\n', xlsx, strjoin(lsheets, ', '));
else
  T = table();
end

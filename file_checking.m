function file_checking(baseDir,outDir)

% file_check.m
% Appends dimensions for MT, B0, B1, T2AX, mFFE, MFA to:
%   - Sheet "DIMS_SCT"  ( *_registered_*_1.nii.gz )
%   - Sheet "DIMS_ANTS" ( *_registered_*_ANTS.nii.gz )
% Columns: ID, Modality, X, Y, Z, T. Missing → "NA".

% ============= USER SETTINGS =============
%baseDir   = fullfile(pwd,'Processed','HC');         % subjects live here: baseDir/<ID>/
%outDir    = fullfile(pwd,'Processed','FullFit_HC'); % workbook lives/accumulates here
xlsx      = fullfile(outDir,'FILE_CHECK.xlsx');     % same workbook as fits             
% =========================================

if ~exist(outDir,'dir'), mkdir(outDir); end

% discover subjects
d = dir(baseDir);
isSubj = [d.isdir] & ~ismember({d.name},{'.','..'});
subjList = {d(isSubj).name};

% --- per-reg keys & patterns (order must match length) ---
KEYS_SCT = {'MT','B0','B1','T2AX','T1','mFFE','MFA','GM','WM','CSF'};
PAT_SCT  = { ...
  '%s_registered_MT_1.nii', ...
  '%s_registered_B0_1.nii', ...
  '%s_registered_B1_1.nii', ...
  '%s_registered_T2AX_1.nii', ...
  '%s_registered_T1_1.nii', ...
  '%s_registered_mFFE_1.nii', ...
  '%s_registered_MFA_1.nii', ...
  '%s_registered_GM.nii', ...
  '%s_registered_WM.nii', ...
  '%s_registered_CSF.nii' ...
};

% ANTS may not have mFFE or tissues in your pipeline; include them if you want NA rows
KEYS_ANTS = {'MT','B0','B1','T2AX','T1','MFA','LYMPH'};
PAT_ANTS  = { ...
  '%s_registered_MT_ANTS.nii', ...
  '%s_registered_B0_ANTS.nii', ...
  '%s_registered_B1_ANTS.nii', ...
  '%s_registered_T2AX_ANTS.nii', ...
  '%s_registered_MFA_ANTS.nii', ...
  '%s_registered_LYMPH_ANTS.nii', ...
  '%s_registered_T1_ANTS.nii' ...
};

dcc_SCT  = cell(0,6);   % {ID, Modality, X, Y, Z, T}
dcc_ANTS = cell(0,6);

for ii = 1:numel(subjList)
    idstr  = subjList{ii};
    subjDir= fullfile(baseDir, idstr);

    % ---- SCT set ----
    for k = 1:numel(KEYS_SCT)
        f = first_existing( ...
              fullfile(subjDir, sprintf([PAT_SCT{k} '.gz'], idstr)), ...
              fullfile(subjDir, sprintf( PAT_SCT{k}      , idstr)) );
        dims = get_nifti_dims(f);           % [X Y Z T] or NaNs
        dcc_SCT(end+1,:) = [{idstr, KEYS_SCT{k}}, dims_to_str(dims)];
    end

    % ---- ANTS set ----
    for k = 1:numel(KEYS_ANTS)
        f = first_existing( ...
              fullfile(subjDir, sprintf([PAT_ANTS{k} '.gz'], idstr)), ...
              fullfile(subjDir, sprintf( PAT_ANTS{k}      , idstr)) );
        dims = get_nifti_dims(f);
        dcc_ANTS(end+1,:) = [{idstr, KEYS_ANTS{k}}, dims_to_str(dims)];
    end
end

% ----- write/append to Excel -----
if ~isempty(dcc_SCT)
    T_SCT = cell2table(dcc_SCT, 'VariableNames', {'ID','Modality','X','Y','Z','T'});
    if exist(xlsx,'file') && any(strcmpi(sheetnames(xlsx),'DIMS_SCT'))
        writetable(T_SCT, xlsx, 'Sheet','DIMS_SCT', ...
            'WriteMode','append', 'WriteVariableNames', false);
    else
        writetable(T_SCT, xlsx, 'Sheet','DIMS_SCT', 'FileType','spreadsheet');
    end
end

if ~isempty(dcc_ANTS)
    T_ANTS = cell2table(dcc_ANTS, 'VariableNames', {'ID','Modality','X','Y','Z','T'});
    if exist(xlsx,'file') && any(strcmpi(sheetnames(xlsx),'DIMS_ANTS'))
        writetable(T_ANTS, xlsx, 'Sheet','DIMS_ANTS', ...
            'WriteMode','append', 'WriteVariableNames', false);
    else
        writetable(T_ANTS, xlsx, 'Sheet','DIMS_ANTS', 'FileType','spreadsheet');
    end
end

fprintf('Updated %s (sheets: DIMS_SCT, DIMS_ANTS)\n', xlsx);

% ================= helpers =================
function path = first_existing(varargin)
% returns first existing path; '' if none
path = '';
for i=1:nargin
    if exist(varargin{i},'file'), path = varargin{i}; return; end
end
end

function dims = get_nifti_dims(fpath)
% Returns [X Y Z T]; if missing/corrupt -> NaNs; T=1 for 3D.
if ~isempty(fpath) && exist(fpath,'file')
    try
        info = niftiinfo(fpath);
        sz = info.ImageSize;
        if numel(sz) < 4, sz(end+1:4) = 1; end
        dims = double(sz(1:4));
    catch
        dims = [NaN NaN NaN NaN];
    end
else
    dims = [NaN NaN NaN NaN];
end
end

function out = dims_to_str(dims)
% dims: [X Y Z T] doubles (may be NaN) -> 1x4 cellstr with "NA" for NaNs
s = string(dims);
s(isnan(dims)) = "NA";
out = cellstr(s(:).');
end

end

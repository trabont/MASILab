function fullFit(baseDir, saveDir, fileID)
% Running Requirements for fullFit.m:
%   baseDir: directory with "*_registered_*.nii.gz" files ('/MS' or '/HC')
%   saveDir: output directory ('FullFit_MS' or 'FullFit_HC')
%   fileID: subject ID (numeric or string)
% -------
% qMT Voxelwise Analysis
% Non-interactive full-fit for one subject.
% Jul 21 2025
% -------
% Program Details:
% (1) Retrieve data: MT, T1, T2AX, MFA, B0, and B1
% (2) Retieve masks: WM, GM, CSF, and Lymph Node
% (3) Perform Yarnykh Full Fit Analysis
% (4) Produce and save the maps: PSR, R1Obs, kba, etc.
% (5) Produce and save images: PSR overlay, Frequency Offset Figure
% -------

% Ensure saveDir exists
addpath("functions/")
  % ensure saveDir exists
  if ~isfolder(saveDir)
    mkdir(saveDir);
  end

  % helper to build filenames
  fn_1 = @(key) fullfile(baseDir, sprintf('%d/%d_registered_%s_1.nii.gz', fileID, fileID, key));
  fn_2 = @(key) fullfile(baseDir, sprintf('%d/%d_registered_%s_ANTS.nii.gz', fileID, fileID, key));
  fn_3 = @(key) fullfile(baseDir, sprintf('%d/%d_registered_%s.nii.gz', fileID, fileID, key));

  %--- Load data ---
  data.MT_SCT   = double(niftiread(fn_1('MT')));    % X×Y×Z×16
  data.T2AX_SCT = double(niftiread(fn_1('T2AX')));   % X×Y×Z
  data.B0_SCT   = double(niftiread(fn_1('B0')));
  data.B1_SCT   = double(niftiread(fn_1('B1')));
  data.MFA_SCT  = double(niftiread(fn_1('MFA')));  % X×Y×Z×dyn

  data.MT_ANTS   = double(niftiread(fn_2('MT')));    % X×Y×Z×16
  data.T2AX_ANTS = double(niftiread(fn_2('T2AX')));   % X×Y×Z
  data.B0_ANTS   = double(niftiread(fn_2('B0')));
  data.B1_ANTS   = double(niftiread(fn_2('B1')));
  data.MFA_ANTS  = double(niftiread(fn_2('MFA')));  % X×Y×Z×dyn

  %--- Load ROI masks ---
  masks.LYMPH = logical(niftiread(fn_2('LYMPH')));
  masks.WM    = logical(niftiread(fn_3('WM')));
  masks.GM    = logical(niftiread(fn_3('GM')));
  masks.CSF   = logical(niftiread(fn_3('CSF')));
  roiNames = fieldnames(masks);

  %--- Determine slices with LYMPH mask ---
  sliceSum = squeeze(sum(sum(masks.LYMPH,1),2));
  slices = find(sliceSum > 6);
  % Exclude if no Lymph Masks greater than 6 voxels
  if isempty(slices)
    warning('Subject %d: no LYMPH mask → skipping', fileID);
    return;
  end
  % Exclude slice if lymph mask is on first or last slice (i.e. slice 1 or 11)
  for ct = 1:length(slices)
      if slices(ct) < 2
          warning('Subject %d: MASK ON BAD SLICE', fileID);
          slices(ct) = 0;
      elseif slices(ct) > 10
          warning('Subject %d: MASK ON BAD SLICE', fileID);
          slices(ct) = 0;
      else
          slices(ct) = slices(ct);
      end
  end
  
  slices = slices(find(slices > 0));
  % Check to see if lymph mask slice/slices meet inclusion criterion (above)
  if isempty(slices)
    warning('Subject %d: no LYMPH mask on reasonable slice → skipping', fileID);
    % return;
  end

  %--- qMT parameters ---
  rf_offset = [1000,1500,2000,2500,8000,16000,32000,100000];
  deg       = [360, 820];
  scan = length(rf_offset);
  Parms.T1flip = 30;
  Parms.MFA     = 6;
  Parms.T1TR    = 50;
  Parms.deltaMT = rf_offset;
  Parms.pwMT    = [20e-3,20e-3];
  Parms.MT_flip = deg;
  Parms.qMTflip = 6;
  Parms.TR      = 50;

  lines  = {'SL','L','G'};
  shapes = {'super-lorentzian','lorentzian','gaussian'};
  total_PSR_map = zeros(256,256,3);

  %--- Loop over slices ---
  for zz = slices'
    % Imaging Data for WM, GM, and CSF (SCT registration)
    MT2D(:,:,:,1)  = squeeze(data.MT_SCT(:,:,zz,:));
    B02D(:,:,1)  = data.B0_SCT(:,:,zz);
    B12D(:,:,1)  = data.B1_SCT(:,:,zz);
    MFA2D(:,:,:,1) = squeeze(data.MFA_SCT(:,:,zz,:));

    % Imaging Data for Lymph Nodes (ANTs registration)
    MT2D(:,:,:,2)  = squeeze(data.MT_ANTS(:,:,zz,:));
    B02D(:,:,2)  = data.B0_ANTS(:,:,zz);
    B12D(:,:,2)  = data.B1_ANTS(:,:,zz);
    MFA2D(:,:,:,2) = squeeze(data.MFA_ANTS(:,:,zz,:));

    %--- Voxelwise fits for each lineshape & ROI ---
    for L = 1:3 % Loop through lineshapes
      Parms.lineshape = shapes{L};
      for R = 1:4 % Loop through tissues 
        if R == 1
            vz=2; % Lymph node imaging data
        else
            vz=1; % WM, GM, and CSF imaging data
        end

        mask2D = masks.(roiNames{R})(:,:,zz);
        [X,Y,~] = size(MT2D);
        % Generate empty array for tissue
        PSR_map   = zeros(X,Y);
        kba_map   = zeros(X,Y);
        T2a_map   = zeros(X,Y);
        T2b_map   = zeros(X,Y);
        R1obs_map = zeros(X,Y);

        fprintf('Subject %d slice %d → %s on %s\n', fileID, zz, shapes{L}, roiNames{R});
        for i = 1:X % loop through x coordinates
          for j = 1:Y % loop through y coordinates
            if mask2D(i,j) == 1 % only run if there is a voxel for the mask coordinate
              % MT Data:
              M1 = squeeze(MT2D(i,j,1:scan,vz));
              M2 = squeeze(MT2D(i,j,scan+1:end,vz));
              M1n = M1 ./ M1(end);
              M2n = M2 ./ M2(end);
              Parms.M = [M1n(:), M2n(:)]; % Normalized MT Data
              Parms.B1    = B12D(i,j,vz);
              Parms.B0    = B02D(i,j,vz);
              Parms.Ernst = squeeze(MFA2D(i,j,1:Parms.MFA,vz))';
              % Full-fit
              [PSR, kba, T2a, T2b, R1obs] = Analysis_Yarnykh_Full_Fit(Parms,shapes{L});
              % Add voxel values to maps at coordinate location
              PSR_map(i,j)   = PSR;
              kba_map(i,j)   = kba;
              T2a_map(i,j)   = T2a;
              T2b_map(i,j)   = T2b;
              R1obs_map(i,j) = R1obs;
            end % mask loop
          end % y-coord loop
        end % x-coord loop

        %--- Save matrices ---
        % Saving the maps for each slice x lineshape x tissue
        outname = fullfile(saveDir, sprintf('%d_slice%02d_%s_%s.mat', fileID, zz, lines{L}, roiNames{R}));
        save(outname, 'PSR_map', 'kba_map', 'T2a_map', 'T2b_map', 'R1obs_map');
        total_PSR_map(:,:,L) = total_PSR_map(:,:,L) + PSR_map; % Creating one PSR map with all tissues, seperated by lineshape
      end % tissue loop
    end % lineshape loop

    %--- MT dynamics figure ---
    % figure shows all contrast differences for each offset and power in a 2degx8offset grid
    fig1 = figure('Visible','off');
    set(fig1,'Position', [750, 800, 1455, 400]);
    tiledlayout(numel(deg), numel(rf_offset), 'TileSpacing','none','Padding','none');
    for rr = 1:numel(deg)
      for cc = 1:numel(rf_offset)
        nexttile;
        imagesc(MT2D(:,:, (rr-1)*numel(rf_offset) + cc));
        sgtitle(sprintf('MT Dynamics: File %d Slice %d', fileID, zz),'FontWeight','bold');
        axis image off; colormap gray;
        if cc==1, ylabel(sprintf('%d°', deg(rr)),'FontWeight','bold'); end
        if rr==1, title(sprintf('%d Hz', rf_offset(cc)),'FontWeight','bold'); end
      end
    end
    exportgraphics(fig1, fullfile(saveDir, sprintf('%d_slice%02d_MTdynamics.png', fileID, zz)));
    close(fig1);

    %--- PSR overlay (wrapped) ---
    base = mat2gray(MT2D(:,:,scan+1));
    cmap = jet(256);
    alpha = 0.6;
    PSR_LYM = total_PSR_map(:,:,1);              % using the combined PSR map (1 = super-lorent)
    vals = PSR_LYM(PSR_LYM>0);
    psr_norm = (PSR_LYM-min(vals))/(max(vals)-min(vals));
    rgb = ind2rgb(uint8(psr_norm*255),cmap);
    overlay = repmat(base,1,1,3);
    mask2  = PSR_LYM>0;
    for c=1:3
  	tmp   = overlay(:,:,c);
     	rgb_c = rgb(:,:,c);                       % extract channel c
     	tmp(mask2) = (1-alpha)*tmp(mask2) ...     % only index the masked pixels
              + alpha*rgb_c(mask2);
  	overlay(:,:,c) = tmp;
    end

    fig2 = figure('Visible','off');
    imshow(overlay);
    exportgraphics(fig2, fullfile(saveDir, sprintf('%d_slice%02d_PSR_overlay.png',fileID,zz)));
    close(fig2);

  end
end

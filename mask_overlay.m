function mask_overlay(rootDir)
% MASK_OVERLAY Loop through ID folders, find a lymph slice with >10 voxels,
% overlay GM/WM/CSF/LYMPH on MT/T2AX images with transparency, and save JPGs.
%
%   mask_overlay;                             % uses pwd
%   mask_overlay('/full/path/to/HC');         % specify root

if nargin<1 || isempty(rootDir)
    rootDir = pwd;
end

maskDir = fullfile(rootDir,'MaskOverlays');
if ~exist(maskDir,'dir')
    mkdir(maskDir);
end

listing = dir(rootDir);
for i = 1:numel(listing)
    name = listing(i).name;
    % only numeric folders, skip MaskOverlays
    if listing(i).isdir && ~ismember(name, {'.','..','MaskOverlays'}) ...
            && ~isempty(regexp(name,'^\d+$','once'))
        
        idDir = fullfile(rootDir,name);
        
        % locate files
        MT1_f   = dir(fullfile(idDir,'*_registered_MT_1.nii.gz'));
        T2AX1_f = dir(fullfile(idDir,'*_registered_T2AX_1.nii.gz'));
        MTa_f   = dir(fullfile(idDir,'*_registered_MT_ANTS.nii.gz'));
        T2a_f   = dir(fullfile(idDir,'*_registered_T2AX_ANTS.nii.gz'));
        GM_f    = dir(fullfile(idDir,'*_registered_GM.nii.gz'));
        WM_f    = dir(fullfile(idDir,'*_registered_WM.nii.gz'));
        CSF_f   = dir(fullfile(idDir,'*_registered_CSF.nii.gz'));
        LYM_f   = dir(fullfile(idDir,'*_registered_LYMPH_ANTS.nii.gz'));
        
        try
            MT_1      = niftiread(fullfile(idDir, MT1_f(1).name));
            T2AX_1    = niftiread(fullfile(idDir, T2AX1_f(1).name));
            MT_ANTS   = niftiread(fullfile(idDir, MTa_f(1).name));
            T2AX_ANTS = niftiread(fullfile(idDir, T2a_f(1).name));
            GM        = niftiread(fullfile(idDir, GM_f(1).name));
            WM        = niftiread(fullfile(idDir, WM_f(1).name));
            CSF       = niftiread(fullfile(idDir, CSF_f(1).name));
            LYMPH     = niftiread(fullfile(idDir, LYM_f(1).name));
        catch
            warning('Skipping %s: missing or unreadable file.', name);
            continue;
        end
        
        % find first slice in lymph mask with >10 voxels
        nz = size(LYMPH,3);
        slice = find(arrayfun(@(z) nnz(LYMPH(:,:,z)), 1:nz) > 8, 1);
        if isempty(slice)
            warning('No slice > 8 voxels in LYMPH for %s', name);
            continue;
        end
        
        % extract 2D slices
        MT1_slice   = squeeze(MT_1(:,:,slice,9));
        T2AX1_slice = squeeze(T2AX_1(:,:,slice));
        MTa_slice   = squeeze(MT_ANTS(:,:,slice,9));
        T2a_slice   = squeeze(T2AX_ANTS(:,:,slice));
        
        gm  = GM(:,:,slice) > 0;
        wm  = WM(:,:,slice) > 0;
        csf = CSF(:,:,slice) > 0;
        ly  = LYMPH(:,:,slice) > 0;
        
        [H,W] = size(gm);
        alphaVal = 0.75;
        redMask    = cat(3, ones(H,W), zeros(H,W), zeros(H,W));
        blueMask   = cat(3, zeros(H,W), zeros(H,W), ones(H,W));
        yellowMask = cat(3, ones(H,W), ones(H,W), zeros(H,W));
        greenMask  = cat(3, zeros(H,W), ones(H,W), zeros(H,W));
        
        hFig = figure('Visible','off','Position',[100 100 800 800]);
        
        % MT_1 + GM(red)/WM(blue)/CSF(yellow)
        subplot(2,2,1)
        imshow(MT1_slice,[],'InitialMagnification','fit'), hold on
        if any(gm(:)),    h=imshow(redMask);    set(h,'AlphaData',gm*alphaVal);    end
        if any(wm(:)),    h=imshow(blueMask);   set(h,'AlphaData',wm*alphaVal);    end
        if any(csf(:)),   h=imshow(yellowMask); set(h,'AlphaData',csf*alphaVal);   end
        title(sprintf('%s: MT_1', name),'Interpreter','none')
        hold off
        
        % T2AX_1 + same masks
        subplot(2,2,2)
        imshow(T2AX1_slice,[],'InitialMagnification','fit'), hold on
        if any(gm(:)),    h=imshow(redMask);    set(h,'AlphaData',gm*alphaVal);    end
        if any(wm(:)),    h=imshow(blueMask);   set(h,'AlphaData',wm*alphaVal);    end
        if any(csf(:)),   h=imshow(yellowMask); set(h,'AlphaData',csf*alphaVal);   end
        title('T2AX_1')
        hold off
        
        % MT_ANTS + LYMPH(green)
        subplot(2,2,3)
        imshow(MTa_slice,[],'InitialMagnification','fit'), hold on
        if any(ly(:)),    h=imshow(greenMask);  set(h,'AlphaData',ly*alphaVal);    end
        title('MT_ANTS')
        hold off
        
        % T2AX_ANTS + LYMPH(green)
        subplot(2,2,4)
        imshow(T2a_slice,[],'InitialMagnification','fit'), hold on
        if any(ly(:)),    h=imshow(greenMask);  set(h,'AlphaData',ly*alphaVal);    end
        title('T2AX_ANTS')
        hold off
        
        % save figure
        outname = fullfile(maskDir, [name '_maskOverlay.jpg']);
        saveas(hFig, outname);
        close(hFig);
    end
end
end

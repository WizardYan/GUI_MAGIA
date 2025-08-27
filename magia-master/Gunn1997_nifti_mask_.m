function parametric_images = Gunn1997_nifti_mask_(lb, ub, theta3_lb,theta3_ub,nBases,decaytime,refTAC,frametimes,filenamedyn,maskfilename,outputdir, atlasfilename, roi_ids, roi_names)
% NEW args at end (optional):
%   atlasfilename : 3D integer label atlas in PET space (NIfTI)
%   roi_ids       : vector of ROI label IDs to report (default: all nonzero labels)
%   roi_names     : cellstr of names (default: 'ROI_<id>')
%
% Returns (unchanged): parametric_images (voxelwise NIfTIs)
%
% Also writes ROI CSV to: <outputdir>/<PETbasename>_bfsrtm_roi_results.csv

frameMid = 0.5.*(frametimes(:,1)+frametimes(:,2));
frameDur = frametimes(:,2)-frametimes(:,1);
lambda = log(2)/decaytime;%C-11 decay constant: 20.4 min
decayRemove = 2.^-(frameMid/decaytime);

%% wzyan (keep your overrides)
decayRemove = ones(length(frametimes),1); %%%%% revised by Weizheng Yan;
lambda = 0; %%%% revised by Weizheng Yan;

%% 
dynhdr = spm_vol_nifti(filenamedyn);
dynVol = nifti(filenamedyn);

if dynhdr.private.dat.dim(4) ~= length(frameMid)
    error(['Number of frame times ' num2str(length(frameMid)) ...
        ' does not match number of dynamic frames in image '  ...
        num2str(dynhdr.private.dat.dim(4))]);
end

%masking
maskhdr = spm_vol_nifti(maskfilename);
maskVol = nifti(maskfilename);
MASK = zeros(maskhdr.dim);
MASK(:,:,:) = maskVol.dat(:,:,:);

%% calculate weights
tots = zeros(length(frameMid),1);
DYNframe = zeros(dynhdr.dim);
for f=1:length(frameMid) %calculate mean framewise
    helper_img = dynVol.dat(:,:,:,f);
    helper_img(isnan(helper_img)) = 0;
    DYNframe(:,:,:) = helper_img;
    tots(f) = decayRemove(f)*sum(DYNframe(:));
end
clear DYNframe
weights = (frameDur.^2)./ max(eps,tots);
weights(isinf(weights)) = 0; %if frame mean==0;
weights = weights/mean(weights); % normalize
W = diag(sqrt(weights));

% compute basis functions: Cr @ exp(-theta3*t)
nFrames = size(frametimes, 1);
B = zeros([nFrames nBases]); % the basis functions
bases_lb=log(theta3_lb)/log(10);
bases_ub=log(theta3_ub)/log(10);
theta3 = logspace(bases_lb, bases_ub, nBases);

% evaluate input function Cr
Cr =  decayRemove.*refTAC;
interp_refdata = pchip([0 ; frameMid],[0 ; Cr]);
inTol = 1e-4; %integration tolerance
trace=0; %quad option
convIntegrand = @(tau, funH, parm, theta3, t, varargin)  funH(parm, t - tau, varargin{:}) .* exp(-theta3 * tau);

for i = 1:nBases
    for j = 1:nFrames
        B(j,i) = quad(convIntegrand, 0, frameMid(j), inTol, trace, @ppval, interp_refdata, theta3(i), frameMid(j));
    end
end

% Precompute QR projectors (unchanged)
M = zeros([2*nBases nFrames]);
for i = 1:nBases
    A = [Cr B(:,i)];
    [Q, R] = qr(W * A);
    M(2*(i-1)+[1 2], :) = R \ (Q.');
end
clear A Q R;

Ct = zeros([nFrames nBases]);
theta = zeros([2 nBases]);

BP_img_3D = zeros(maskhdr.dim);
RI_img_3D = zeros(maskhdr.dim);
k2_img_3D = zeros(maskhdr.dim);
theta3_img_3D = zeros(maskhdr.dim);

DYNplane = zeros([dynVol.dat.dim(1) dynVol.dat.dim(2) dynVol.dat.dim(4)]);
fprintf(1,' \nStarting to fit bf-SRTM (voxelwise).\n');
tic

for pl = 1:dynVol.dat.dim(3)
    if (mod(pl,10)==0), disp(['fitting plane: ' int2str(pl)]); end
    DYNplane(:,:,:) = dynVol.dat(:,:,pl,:);
    [mask_iX,mask_iY] = find(MASK(:,:,pl));
    TACs = zeros(nFrames,length(mask_iX));
    nTACs = length(mask_iX);
    RI = zeros([1 nTACs]);
    k2 = zeros([1 nTACs]);
    BP = zeros([1 nTACs]);
    th3 = zeros([1 nTACs]);
    for T = 1:nTACs
        CPET    = DYNplane(mask_iX(T), mask_iY(T),:);
        TACs(:,T)=decayRemove.*CPET(:);
    end
    TACs(isnan(TACs)) = 0 ;   %replace NaNs with zeros
    clear CPET
    
    for j = 1:nTACs
        H = W * TACs(:, j);
        for i = 1:nBases
            theta(:,i) = M(2*(i-1)+[1 2], :) * H;
            Ct(:,i) = theta(1,i) * Cr + theta(2,i) * B(:,i);
        end
        RSS = sum(repmat(weights, [1 nBases]) .* ((repmat(TACs(:, j), [1 nBases]) - Ct) .^ 2), 1);
        [RI(j),k2(j),BP(j),th3(j)] = bound_srtm_parameters(RSS,theta,theta3,lambda,lb,ub);
    end
    
    for T=1:length(mask_iX)
        BP_img_3D(mask_iX(T),mask_iY(T),pl)   = BP(T);
        RI_img_3D(mask_iX(T),mask_iY(T),pl)   = RI(T);
        k2_img_3D(mask_iX(T),mask_iY(T),pl)   = k2(T);
        theta3_img_3D(mask_iX(T),mask_iY(T),pl)   = th3(T);
    end
    clear RI k2 BP RSS TACs th3
end

[~,filename] = fileparts(filenamedyn);

parametric_images = cell(3,1);
parametric_images{1} = fullfile(outputdir,[filename '_bfsrtm_BP.nii']);
parametric_images{2} = fullfile(outputdir,[filename '_bfsrtm_R1.nii']);
parametric_images{3} = fullfile(outputdir,[filename '_bfsrtm_k2.nii']);

VO = maskhdr;
VO.dt = [spm_type('int16') spm_platform('bigend')];
VO.pinfo = [Inf Inf Inf]';
VO.fname = parametric_images{1}; spm_write_vol(VO,BP_img_3D);
VO.fname = parametric_images{2}; spm_write_vol(VO,RI_img_3D);
VO.fname = parametric_images{3}; spm_write_vol(VO,k2_img_3D);
VO.fname = fullfile(outputdir,[filename '_bfsrtm_theta3.nii']);
spm_write_vol(VO,theta3_img_3D);


%% ================= ROI table from voxel maps (NEW) =================
% Build ROI-wise table by aggregating voxelwise BFSRTM maps (BP/R1/k2/theta3)
% Requires: atlasfilename (PET-space), optional roi_ids, roi_names
% Uses existing variables: MASK, BP_img_3D, RI_img_3D, k2_img_3D, theta3_img_3D, lambda,
%                          maskhdr/dynhdr (for dimension check), filenamedyn, outputdir

if nargin >= 12 && ~isempty(atlasfilename)
    fprintf('\nBuilding ROI table from voxel maps using atlas: %s\n', atlasfilename);

    % ---- Load and validate atlas ----
    athdr  = spm_vol_nifti(atlasfilename);
    atlVol = nifti(atlasfilename);

    % Dimension check (3D shape must match PET/mask)
    if ~isequal(athdr.dim(1:3), maskhdr.dim(1:3))
        error('Atlas and PET/mask dimensions differ; please resample the atlas to PET space first.');
    end

    ATLAS = double(atlVol.dat(:,:,:));
    ATLAS(isnan(ATLAS)) = 0;

    % ---- ROI selection / names ----
    if nargin < 13 || isempty(roi_ids)
        roi_ids = unique(ATLAS(:));
        roi_ids = roi_ids(roi_ids ~= 0); % exclude background
    else
        roi_ids = roi_ids(:);
    end

    if nargin < 14 || isempty(roi_names)
        roi_names = arrayfun(@(x) sprintf('ROI_%d', x), roi_ids, 'UniformOutput', false);
    end
    if numel(roi_names) ~= numel(roi_ids)
        error('roi_names must match roi_ids in length.');
    end

    % ---- Pull voxel maps ----
    BPv  = BP_img_3D;      % BPnd
    R1v  = RI_img_3D;      % R1
    k2v  = k2_img_3D;      % k2
    th3v = theta3_img_3D;  % theta3 = k2a + lambda

    % ---- Preallocate outputs ----
    nR = numel(roi_ids);
    Nvox   = zeros(nR,1);
    R1_out = nan(nR,1);
    k2_out = nan(nR,1);
    BP_out = nan(nR,1);
    k2a_out= nan(nR,1);

    % ---- Aggregate per ROI ----
    useMask = (MASK ~= 0);        % analysis mask
    for r = 1:nR
        rid = roi_ids(r);
        m = (ATLAS == rid) & useMask;  % ROI ∩ analysis mask
        Nvox(r) = nnz(m);
        if Nvox(r) == 0
            warning('ROI %s (ID=%d) has 0 voxels in mask; skipping stats.', roi_names{r}, rid);
            continue;
        end

        % Mean over voxels (omit NaNs)
        BP_out(r)  = mean(BPv(m),  'omitnan');          % BPnd
        R1_out(r)  = mean(R1v(m),  'omitnan');          % R1
        k2_out(r)  = mean(k2v(m),  'omitnan');          % k2
        k2a_out(r) = mean(th3v(m) - lambda, 'omitnan'); % k2a = theta3 - lambda
    end

    % ---- Assemble table (schema must match exactly) ----
    roi_table = table( ...
        roi_ids(:), string(roi_names(:)), Nvox(:), ...
        R1_out(:), k2_out(:), BP_out(:), k2a_out(:), ...
        'VariableNames', {'ROI_ID','ROI_Name','Nvox','R1','k2','BPnd','k2a'});

    % ---- Save CSV ----
    [~,filename] = fileparts(filenamedyn);
    csv_path = fullfile(outputdir, [filename '_bfsrtm_voxels_results.csv']);
    writetable(roi_table, csv_path);
    fprintf('ROI-from-voxel maps table saved: %s\n', csv_path);
end
%% ===================================================================


%% ================= ROI-wise BFS-SRTM (NEW) =================
if nargin >= 12 && ~isempty(atlasfilename)
    fprintf('\nStarting ROI-wise BFS-SRTM using atlas: %s\n', atlasfilename);
    athdr  = spm_vol_nifti(atlasfilename);
    atlVol = nifti(atlasfilename);
    if ~isequal(athdr.dim(1:3), dynhdr.dim(1:3))
        error('Atlas and PET dimensions differ; please resample to PET space first.');
    end

    atlasData = atlVol.dat(:,:,:);
    atlasData(isnan(atlasData)) = 0;

    % Choose ROI IDs
    if nargin < 13 || isempty(roi_ids)
        roi_ids = unique(atlasData(:));
        roi_ids = roi_ids(roi_ids~=0); % exclude background
    else
        roi_ids = roi_ids(:);
    end

    % Names
    if nargin < 14 || isempty(roi_names)
        roi_names = arrayfun(@(x) sprintf('ROI_%d', x), roi_ids, 'UniformOutput', false);
    end
    if numel(roi_names) ~= numel(roi_ids)
        error('roi_names must match roi_ids in length.');
    end

    % Prealloc
    nROI = numel(roi_ids);
    R1s = nan(nROI,1); k2s = nan(nROI,1); BPs = nan(nROI,1); TH3 = nan(nROI,1); Nvox = zeros(nROI,1);

    % Reuse W, Cr, B, M, theta3, weights, frameMid from above
    Ct_roi = zeros(nFrames, nBases);
    theta_roi = zeros(2, nBases);

    for r = 1:nROI
        rid = roi_ids(r);
        roiMask = (atlasData == rid);
        Nvox(r) = nnz(roiMask);
        if Nvox(r) == 0
            warning('ROI %s (ID=%d) has 0 voxels. Skipping.', roi_names{r}, rid);
            continue;
        end

        % Mean TAC over ROI per frame
        TAC = zeros(nFrames,1);
        for f=1:nFrames
            img = dynVol.dat(:,:,:,f);
            vals = img(roiMask);
            vals(isnan(vals)) = 0;
            TAC(f) = mean(vals);
        end
        TAC = decayRemove .* TAC;

        % BFS sweep (reuse projectors)
        H = W * TAC;
        for i=1:nBases
            theta_roi(:,i) = M(2*(i-1)+[1 2], :) * H;
            Ct_roi(:,i)    = theta_roi(1,i) * Cr + theta_roi(2,i) * B(:,i);
        end
        RSS = sum(repmat(weights, [1 nBases]) .* ((repmat(TAC, [1 nBases]) - Ct_roi) .^ 2), 1);
        [R1s(r),k2s(r),BPs(r),TH3(r)] = bound_srtm_parameters(RSS,theta_roi,theta3,lambda,lb,ub);
    end

    % Write CSV
    roi_table = table(roi_ids(:), string(roi_names(:)), Nvox(:), R1s(:), k2s(:), BPs(:), TH3(:), ...
        'VariableNames', {'ROI_ID','ROI_Name','Nvox','R1','k2','BPnd','k2a'});

    csv_path = fullfile(outputdir, [filename '_bfsrtm_roi_results.csv']);

    writetable(roi_table, csv_path);
    fprintf('ROI results saved: %s\n', csv_path);
end
%% ===========================================================

end

function [RI,k2,BP,optim_theta3] = bound_srtm_parameters(RSS,theta,theta3,lambda, lb, ub)
[~,I] = sort(RSS,'ascend');
sorted_theta(1,:) = theta(1,I);
sorted_theta(2,:) = theta(2,I);
sorted_theta3(1,:) = theta3(I);

RI = sorted_theta(1, 1);
k2 = sorted_theta(2, 1) + RI*(sorted_theta3(1) - lambda);
BP = k2 / (sorted_theta3(1) - lambda) - 1;
i = 1;
optim_theta3 = sorted_theta3(i);

while(RI <= lb(1) || RI > ub(1) ||  k2 <= lb(2) || k2 >= ub(2) || BP < lb(3) || BP > ub(3) )
    i = i + 1;
    optim_theta3 = sorted_theta3(i);
    RI = sorted_theta(1, i);
    k2 = sorted_theta(2, i) + RI*(sorted_theta3(i) - lambda);
    BP = k2 / (sorted_theta3(i) - lambda) - 1;
    if(i==length(RSS))
        RI = NaN; k2 = NaN; BP = NaN; break;
    end
end
end

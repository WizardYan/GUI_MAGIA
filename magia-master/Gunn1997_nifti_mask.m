function parametric_images = Gunn1997_nifti_mask(lb, ub, theta3_lb,theta3_ub,nBases,decaytime,refTAC,frametimes,filenamedyn,maskfilename,outputdir)

frameMid = 0.5.*(frametimes(:,1)+frametimes(:,2));
%frameMid = frametimes(:,2);
frameDur = frametimes(:,2)-frametimes(:,1);
lambda = log(2)/decaytime;%C-11 decay constant: 20.4 min
decayRemove = 2.^-(frameMid/decaytime);

%% wzyan
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



%% ---------- Auto-select theta3 range ----------wzyan
% fprintf('Estimating optimal theta3 range using global TAC...\n');
% nFrames = size(frametimes,1);
% global_TAC = zeros(nFrames,1);
% for f = 1:nFrames
%     frame_data = dynVol.dat(:,:,:,f);
%     frame_masked = frame_data(MASK > 0);
%     global_TAC(f) = mean(frame_masked(:));
% end
% global_TAC = decayRemove .* global_TAC;
% % BFS prefit
% wide_theta3 = logspace(-2, 1, 50);
% Cr = decayRemove .* refTAC;
% interp_refdata = pchip([0 ; frameMid],[0 ; Cr]);
% inTol = 1e-4; trace = 0;
% Btest = zeros(nFrames, numel(wide_theta3));
% for i = 1:numel(wide_theta3)
%     for j = 1:nFrames
%         Btest(j,i) = quad(@(tau) ppval(interp_refdata, frameMid(j)-tau).*exp(-wide_theta3(i)*tau), 0, frameMid(j), inTol, trace);
%     end
% end
% rss_all = zeros(numel(wide_theta3),1);
% for i = 1:numel(wide_theta3)
%     A = [Cr, Btest(:,i)];
%     x = A \ global_TAC;
%     fit_TAC = A * x;
%     rss_all(i) = sum((global_TAC - fit_TAC).^2);
% end
% [~, idx_min] = min(rss_all);
% theta3_best = wide_theta3(idx_min);
% theta3_lb = theta3_best * 0.5;
% theta3_ub = theta3_best * 5;
% fprintf('Auto-selected theta3 range: %.4f to %.4f\n', theta3_lb, theta3_ub);


% ---------- Visualization of Global BFS Fit ----------
% figure;
% semilogx(wide_theta3, rss_all, '-o', 'LineWidth', 2);
% xlabel('k2a (theta3)');
% ylabel('Residual Sum of Squares (RSS)');
% title('Global BFS Prefit: RSS vs. k2a');
% grid on;
% 
% % Mark best theta3
% hold on;
% plot(theta3_best, rss_all(idx_min), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% text(theta3_best, rss_all(idx_min), sprintf('  Best k2a = %.4f', theta3_best), ...
%     'VerticalAlignment', 'bottom', 'FontSize', 10);
% 
% % Plot best model fit vs. actual global TAC
% A_best = [Cr, Btest(:,idx_min)];
% x_best = A_best \ global_TAC;
% fit_TAC_best = A_best * x_best;
% 
% figure;
% plot(frameMid, global_TAC, 'ko-', 'DisplayName', 'Global TAC'); hold on;
% plot(frameMid, fit_TAC_best, 'r--', 'LineWidth', 2, 'DisplayName', 'Best BFS Fit');
% xlabel('Time (min)');
% ylabel('Activity');
% title('Best BFS Fit vs. Global TAC');
% legend;
% grid on;


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
weights = (frameDur.^2)./tots;
weights(isinf(weights)) = 0; %if frame mean==0;
%weights = ones(length(frameMid),1); %no weights
W = diag(sqrt(weights));

% compute basis functions: Cr @ exp(-theta3*t)
nFrames = size(frametimes, 1);
B = zeros([nFrames nBases]); % the basis functions
%theta3 = linspace(theta3_lb,theta3_ub,nBases);
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
    % use quad to evaluate the convolution integrals at mid frame time
    % Basis function = integral_0^t (Cr(t-tau) * exp(-theta3*tau)) dtau
    for j = 1:nFrames
        B(j,i) = quad(convIntegrand, 0, frameMid(j), inTol, trace, @ppval, interp_refdata, theta3(i), frameMid(j));
    end
end

M = zeros([2*nBases nFrames]);
for i = 1:nBases
    A = [Cr B(:,i)];
    [Q, R] = qr(W * A);
    M(2*(i-1)+[1 2], :) = R \ (Q.');
end

clear A, clear Q, clear R;

Ct = zeros([nFrames nBases]);
theta = zeros([2 nBases]);

BP_img_3D = zeros(maskhdr.dim);
RI_img_3D = zeros(maskhdr.dim);
k2_img_3D = zeros(maskhdr.dim);
theta3_img_3D = zeros(maskhdr.dim);

DYNplane = zeros([dynVol.dat.dim(1) dynVol.dat.dim(2) dynVol.dat.dim(4)]);
fprintf(1,' \nStarting to fit bf-SRTM.\n');

tic

for pl = 1:dynVol.dat.dim(3)
    
    if (mod(pl,10)==0), disp(['fitting plane: ' int2str(pl)]); end;
    DYNplane(:,:,:) = dynVol.dat(:,:,pl,:);
    % process all TACs from the corresponding plane in one go
    [mask_iX,mask_iY] = find(MASK(:,:,pl));
    TACs        = zeros(nFrames,length(mask_iX));
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
    
    for T=1:length(mask_iX),
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
VO.fname = parametric_images{1};
spm_write_vol(VO,BP_img_3D);
VO.fname = parametric_images{2};
spm_write_vol(VO,RI_img_3D);
VO.fname = parametric_images{3};
spm_write_vol(VO,k2_img_3D);
VO.fname = fullfile(outputdir,[filename '_bfsrtm_theta3.nii']);
spm_write_vol(VO,theta3_img_3D);

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
        RI = NaN;
        k2 = NaN;
        BP = NaN;
        break;
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN_human_NVC.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% author: Brad Rauscher (2025/Dec/16)
% 
% Main analysis script for NVC modeling in humans paper. Uses fMRI BOLD,
% EEG, and pupil data acquired by the Lewis lab at MIT.

%% load data

data_dir = '/Users/bcraus/Library/Mobile Documents/com~apple~CloudDocs/Documents/BostonU/Research/DevorLab/AnalysisCode/NVC_human/Zinong_data.mat';
save_dir = '/Users/bcraus/Library/Mobile Documents/com~apple~CloudDocs/Documents/BostonU/Research/DevorLab/AnalysisCode/NVC_human/Images';

subjects = whos('-file',data_dir);
subjects = {subjects.name};
sub_idx = 1;

data = f_load_Zinong_data(subjects{sub_idx});

save_dir = fullfile(save_dir,subjects{sub_idx});
[~, ~, ~] = mkdir(save_dir);

save_plots = 1;

data.pupil(data.pupil==0) = NaN;
data.pupil = fillmissing(data.pupil,'linear');

%% process EEG

ephys = struct;
[ephys.spg, ~, ephys.f] = f_morlet(data.EEG, data.settings.EEG_fs, [1.5, 100], 200);

%%
bands = struct;
[~, bands.delta] = min(abs(ephys.f' - [0.5, 4]));
[~, bands.theta] = min(abs(ephys.f' - [4, 8]));
[~, bands.alpha] = min(abs(ephys.f' - [8, 12]));
[~, bands.beta] = min(abs(ephys.f' - [12, 30]));
[~, bands.gamma] = min(abs(ephys.f' - [30,50]));

bands.labels = {'delta', 'theta', 'alpha', 'beta', 'gamma'};

for i = 1:numel(bands.labels)
    label = bands.labels{i};
    ephys.(label) = log10(squeeze(mean(ephys.spg(:,bands.(label)(1):bands.(label)(2),:),2)));
    % lb = min(ephys.(label), [], 'all');
    ephys.(label) = resample(ephys.(label),data.settings.MR_fs,data.settings.EEG_fs);
    % ephys.(label)(ephys.(label) < lb) = lb;
end

%% plot signals

t = (1:numel(data.pupil)) / data.settings.MR_fs;
t_EEG = (1:size(data.EEG,1)) / data.settings.EEG_fs;

f = figure;
tiledlayout(3,1);

ax = nexttile;
plot(t,mean(data.BOLD,2));
ylabel('BOLD');
box off;
title('Global Average');

ax(2) = nexttile;
plot(t,data.pupil);
ylabel('Pupil');
box off;

ax(3) = nexttile;
EEG_channel = 1;
f_morlet(data.EEG(:,EEG_channel),500,[1.5,45],300,plot=1);
title(data.EEG_labels{EEG_channel});

for i = 1:numel(ax); ax(i).XLim = [0, t(end)]; end
set(ax(1:2),XTickLabel=[]);

f_savePNG(fullfile(save_dir,'Hb_peaks.png'), flag=save_plots, dim=[4,5]);

%% analyze global ephys patterns

cmp = cmpvir;
cmp = cmp(round(linspace(1,240,5)),:);

f = figure;
tiledlayout(3,2);

for i = 1:5
    ax(i) = nexttile;
    imagesc(corrcoef(ephys.(bands.labels{i})));
    axis image off;
    clim([0, 1]);
    title(bands.labels{i},color=cmp(i,:));
    ylabel('channels');
    ax(i).YLabel.Visible = 'on';
end
c = colorbar;
c.Label.String = 'r';

ax(6) = nexttile;hold on;
idx = ~tril(ones(numel(data.EEG_labels),numel(data.EEG_labels)));
for i = 1:5
    tmp = corrcoef(ephys.(bands.labels{i}));
    bar(i,mean(tmp(idx),'all'),FaceColor=cmp(i,:));
end
ylabel('Average r');
ax(6).XAxis.Visible = 'off';

f_savePNG(fullfile(save_dir,'global_corr.png'), flag=save_plots, dim=[5,4]);

%% analyze local ephys patterns

corr = struct;
corr.EEG = NaN(5,5,numel(data.EEG_labels));

for i = 1:numel(data.EEG_labels)
    tmp = [];
    for c = 1:5
        tmp = [tmp, ephys.(bands.labels{c})(:,i)];
    end
    corr.EEG(:,:,i) = corrcoef(tmp);
end

figure;
imagesc(mean(corr.EEG,3));
caxis([-1, 1]);
colormap cmpbbr;
axis image off;
c = colorbar;
c.Label.String = 'r';
title('Average Band Power R');
ax = gca;
ax.YLabel.Visible = 'on';
ylabel('\gamma    \beta    \alpha    \theta    \delta');

f_savePNG(fullfile(save_dir,'EEG_local_corr.png'), flag=save_plots, dim=[3,4]);

%% correlate signals

% Calculate correlation between EEG bands and pupil

corr.EEG_pupil = NaN(numel(data.EEG_labels),5);
for i = 1:5
    corr.EEG_pupil(:,i) = f_corr(data.pupil,ephys.(bands.labels{i}),1);
end

figure;
imagesc(corr.EEG_pupil);
axis image off;
clim(0.5*[-1, 1]);
colormap cmpbbr;
c = colorbar;
c.Label.String = 'r';
title('\delta \theta  \alpha  \beta \gamma');
ylabel('EEG Ch.');
ax = gca;
ax.YLabel.Visible = 'on';

f_savePNG(fullfile(save_dir,'EEG_vs_pupil.png'), flag=save_plots, dim=[5,3]);

%% Calculate correlation between BOLD and pupil
corr.BOLD_pupil = f_corr(data.BOLD,data.pupil,1);

figure;
bar(corr.BOLD_pupil);
ylabel('r');
xlabel('MRI Ch.')
title('BOLD vs. pupil');
box off;
f_savePNG(fullfile(save_dir,'BOLD_vs_pupil.png'), flag=save_plots, dim=[3,4]);

%% Calculate correlation between BOLD and EEG bands
corr.BOLD_EEG = NaN(numel(data.MR_labels), numel(data.EEG_labels), 5);

tmpBOLD = data.BOLD .* ones(1,1,numel(data.EEG_labels));
for i = 1:5
    label = bands.labels{i};
    tmpEEG = ephys.(label) .* ones(1,1,numel(data.MR_labels));
    tmpEEG = permute(tmpEEG,[1,3,2]);
    corr.BOLD_EEG(:,:,i) = f_corr(tmpBOLD,tmpEEG,1);
end

figure;
tiledlayout(5,1);

for i = 1:5
    ax(i) = nexttile;
    imagesc(corr.BOLD_EEG(:,:,i));
    axis image off;
    clim(0.5*[-1, 1]);
    c = colorbar;
    colormap cmpbbr;
    c.Label.String = 'r';
    ylabel(bands.labels{i})
    ax(i).YLabel.Visible = 'on';
end

ax(1).Title.String = 'MR x EEG Ch.';

f_savePNG(fullfile(save_dir,'BOLD_vs_EEG.png'), flag=save_plots, dim=[5,3]);

%% estimate multi-IRF

fs = data.settings.MR_fs;
win = fs*[-10, 10];

meas_idx = [1,2,3,4,5,6];

perf = zeros(numel(data.MR_labels), numel(data.EEG_labels));
IRFs = zeros(sum(abs(win))+1, numel(meas_idx), numel(data.MR_labels), numel(data.EEG_labels));
params = zeros(3, numel(meas_idx), numel(data.MR_labels), numel(data.EEG_labels));

for MR_channel = 1:numel(data.MR_labels)
    disp(MR_channel);
    for EEG_channel = 1:numel(data.EEG_labels)
    
        X = [data.pupil,...
            ephys.delta(:,EEG_channel),...
            ephys.theta(:,EEG_channel),...
            ephys.alpha(:,EEG_channel),...
            ephys.beta(:,EEG_channel),...
            ephys.gamma(:,EEG_channel)];
        X = X(:,meas_idx);
        % X(:,2:end) = log10(X(:,2:end));
        X = X./std(X,0);

        Y = data.BOLD(:,MR_channel);
        Y = f_bpf(Y, [0, 0.5], fs);
        Y = Y./std(Y,0);

        [IRFs(:, :, MR_channel, EEG_channel), perf(MR_channel, EEG_channel), ~, params(:, :, MR_channel, EEG_channel)] = f_gamma_IRF(X, Y, win);
    end
end

mean(perf,'all')

%% plot 

plotLabels = {'Pupil', '\delta', '\theta', '\alpha', '\beta', '\gamma'};
plotLabels = plotLabels(meas_idx);

t = (win(1):win(2))/fs;
colors = cmpvir;
colors = colors(round(linspace(1,255,numel(meas_idx))),:);

figure;hold on;
tiledlayout(numel(meas_idx),1);
for i = 1:numel(meas_idx)
    ax = nexttile;f_plotLineError(t, mean(IRFs(:,i,:,:), [3,4]), std(IRFs(:,i,:,:), 0, [3,4]),color=colors(i,:));
    axis off;
    ax.YLabel.Visible = 'on';
    ylabel(plotLabels{i});
end
ax.XAxis.Visible = 'on';
xlabel('Time (s)');

f_savePNG(fullfile(save_dir,'multi_IRFs.png'), flag=save_plots, dim=[5,2.5]);

%% plot results

f = figure;
imagesc(perf);
axis image off;
c = colorbar;
clim([0, 1]);
c.Label.String = 'r';
colormap cmpvir;
ylabel('MR channel');
xlabel('EEG channel');
ax = gca;
ax.XLabel.Visible = 'on';
ax.YLabel.Visible = 'on';

f_savePNG(fullfile(save_dir,'multi_IRFs_r.png'), flag=save_plots, dim=[3,5]);

%% FC

win = data.settings.MR_fs * [10, 3];

FC = f_funConGram(data.BOLD, win);
ds_pupil = movmean(data.pupil, win(1), Endpoints='discard');
ds_pupil = ds_pupil(1:win(2):end);

high_pupil = ds_pupil > prctile(ds_pupil, 60);
low_pupil = ds_pupil < prctile(ds_pupil, 40);

FC_high = mean(FC(:,:,high_pupil), 3);
FC_low = mean(FC(:,:,low_pupil), 3);

figure;
tiledlayout(2,1);
ax = nexttile;imagesc(FC_low);
axis image off;
clim([0, 1]);
colormap cmpvir;
c = colorbar;
c.Label.String = 'r';
ax.YLabel.Visible = 'on';
ylabel('Low Pupil');

ax = nexttile;imagesc(FC_high);
axis image off;
clim([0, 1]);
colormap cmpvir;
c = colorbar;
c.Label.String = 'r';
ax.YLabel.Visible = 'on';
ylabel('High Pupil');

%%

f = figure;
tiledlayout(numel(meas_idx),3);
for i = 1:numel(meas_idx)
    ax = nexttile(3*(i-1)+1);
    tmp = squeeze(params(1,i,:,:));
    imagesc(tmp);
    axis image off;
    c = colorbar;
    clim(max(abs(prctile(tmp(:), [3, 97]))) * [-1, 1]);
    ax.YLabel.Visible = 'on';
    ylabel(plotLabels{i});
    if i == 1
        title('A');
    end
    colormap cmpbbr;
end

for i = 1:numel(meas_idx)
    ax = nexttile(3*(i-1)+2);
    tmp = squeeze(params(2,i,:,:));
    imagesc(tmp);
    axis image off;
    c = colorbar;
    clim([0, prctile(tmp(:), 97)]);
    if i == 1
        title('\tau');
    end
    set(ax,'Colormap',cmpinf);
end

for i = 1:numel(meas_idx)
    ax = nexttile(3*(i-1)+3);
    tmp = squeeze(params(3,i,:,:));
    imagesc(tmp / fs);
    axis image off;
    c = colorbar;
    clim(win / fs);
    if i == 1
        title('t_0');
    end
end

f_savePNG(fullfile(save_dir,'multi_IRFs_params.png'), flag=save_plots, dim=[5,5.5]);

%%

% figure;
% plot(f_gamma_IRF(data.pupil, data.BOLD(:, 1), win));

%% IRF functions
function [IRF, perf, pred, y_real, params_hat] = fir_direct_gamma(X, y, win)
num_sigs = size(X,2);

% remove mean
X = X - mean(X);
y = y - mean(y);
y_real = y;

% IRF length
L = sum(abs(win)) + 1;
dt = 1;  % sampling interval
t = (win(1):win(2)) * dt;

%% Objective function
% params = [A1, shape1, scale1, t01, ..., An, shapn, scalen, t0n]
obj_fun = @(params) gamma_multi_err_shift(params, X, y, t, L);

% initial guess
params0 = repmat([1, 3, 1, 0], 1, num_sigs);  % t0 = 0

% bounds
lb = repmat([0, 0.5, 0.1, win(1)], 1, num_sigs);  % allow negative t0
ub = repmat([Inf, 10, 10, win(2)], 1, num_sigs);

% options: silent
opts = optimoptions('lsqnonlin', 'Display','off', ...
    'MaxIterations',200,'MaxFunctionEvaluations',1000);

% fit
params_hat = lsqnonlin(obj_fun, params0, lb, ub, opts);

%% reconstruct IRFs
IRF = zeros(L, num_sigs);
for i = 1:num_sigs
    idx = 4*(i-1)+1;
    IRF(:,i) = gamma_IRF_shift(t, params_hat(idx), params_hat(idx+1), ...
                               params_hat(idx+2), params_hat(idx+3));
end

%% prediction
pred = zeros(size(y));
for i = 1:num_sigs
    h = IRF(:,i);
    Xi = convmtx(X(:,i), L);
    Xi = Xi(1:length(y), :);
    pred = pred + Xi * h(:);
end

perf = f_corr(pred, y_real, 1);

end

%% --- gamma IRF with shift ---
function h = gamma_IRF_shift(t, A, shape, scale, t0)
    eps0 = 1e-6;
    t_safe = max(t - t0, eps0);
    h_raw = t_safe.^(shape-1) .* exp(-t_safe / scale);
    h_raw(t < t0) = 0;  % zero for times before shift
    h = A * h_raw / (scale^shape * gamma(shape));
end

%% --- multi-predictor error for shifted gamma ---
function err = gamma_multi_err_shift(params, X, y, t, L)
num_sigs = size(X,2);
pred = zeros(size(y));
for i = 1:num_sigs
    idx = 4*(i-1)+1;
    A = params(idx);
    shape = params(idx+1);
    scale = params(idx+2);
    t0 = params(idx+3);
    h = gamma_IRF_shift(t, A, shape, scale, t0);
    Xi = convmtx(X(:,i), L);
    Xi = Xi(1:length(y), :);
    pred = pred + Xi * h(:);
end
err = y - pred;
end

function [IRF, perf, pred, y_real] = fir_direct(x, y, win)

num_sigs = size(x,2);

x = x - mean(x);
y = y - mean(y);

y_real = y;

N = size(x,1);
L = sum(abs(win))+1;

% ---------- build temporal basis ----------
nbasis = 8;   % <<<<<< adjust (6–10 typical)
B = make_basis(L, nbasis);   % L × nbasis

% ---------- build convolution matrix ----------
X = [];
for i = 1:num_sigs
    Xi = convmtx(x(:,i), L);     % (N+L-1) × L
    Xi = Xi * B;                % project to basis → (N+L-1) × nbasis
    X  = [X, Xi];
end

x_real = X;
x_real(end-win(2)+1:end,:) = [];
x_real(1:-1*win(1),:) = [];

X(N+1:end,:) = [];
X(1:L-1,:) = [];
y(end+win(1)+1:end) = [];
y(1:win(2)) = [];

% ---------- estimate basis coefficients ----------
lambda = 1e-2 * trace(X'*X) / size(X,2);
cb = (X'*X + lambda*eye(size(X,2))) \ (X'*y);

% ---------- prediction ----------
pred = x_real * cb;

% ---------- reconstruct full IRFs ----------
IRF = zeros(L, num_sigs);
idx = 0;
for i = 1:num_sigs
    ci = cb(idx + (1:nbasis));
    IRF(:,i) = B * ci;
    idx = idx + nbasis;
end

perf = f_corr(pred,y_real,1);

end

function B = make_basis(L, nbasis)

t = linspace(0,1,L)';
centers = linspace(0,1,nbasis);
width = centers(2) - centers(1);

B = zeros(L, nbasis);
for k = 1:nbasis
    x = (t - centers(k)) / width;
    B(:,k) = (cos(pi*min(max(x,-1),1)) + 1)/2;
end

B = B ./ sqrt(sum(B.^2,1));  % normalize
end

% function [IRF, perf, pred, y_real] = fir_direct(x, y, win)
% 
% x = x(:);
% y = y(:);
% 
% x = x - mean(x);
% y = y - mean(y);
% 
% y_real = y;
% 
% 
% N = numel(x);
% L = sum(abs(win))+1;
% 
% x = convmtx(x(:),L);
% x_real = x;
% x_real(end-win(2)+1:end,:) = [];
% x_real(1:-1*win(1),:) = [];
% 
% x(N+1:end,:) = [];
% x(1:L-1,:) = [];
% y(end+win(1)+1:end) = [];
% y(1:win(2)) = [];
% 
% % IRF = x \ y;
% lambda = 1e-2 * trace(x'*x) / size(x,2); %%
% IRF = (x'*x + lambda*eye(size(x,2))) \ (x'*y); %%
% 
% pred = x_real * IRF;
% 
% perf = f_corr(pred,y_real,1);
% 
% end

function IRF = fir_kernel(x,y,fs,l)

alpha = 1e-3;

N = min(length(x),length(y));
x = x(1:N);
y = y(1:N);

x = x-mean(x);
y = y-mean(y);

M = round(l*fs);
nfft = 2^nextpow2(2*N);

X = fft(x,nfft);
Y = fft(y,nfft);

Sxx = X.*conj(X);
Syx = Y.*conj(X);

rxx_full = ifft(Sxx, nfft, 'symmetric');
ryx_full = ifft(Syx,nfft,'symmetric');

rxx = rxx_full(1:M);
ryx = ryx_full(1:M);

Rxx = toeplitz(rxx);

lam = alpha * median(abs(rxx));
Rxx = Rxx + lam*eye(M);

IRF = Rxx \ ryx;

end

function [H_all, tlag] = fir_multi_kernel(Xmat, y, fs, l)
% Xmat : N x P matrix of predictors (each col is a signal)
% y    : N x 1 target signal
% fs   : sampling rate
% Lsec : FIR length (seconds)

alpha = 1e-3;
[N, P] = size(Xmat);
M = round(l * fs);
nfft = 2^nextpow2(2*N);

% Demean all
Xmat = Xmat - mean(Xmat, 1);
y    = y - mean(y);

% Preallocate block autocorr matrix R and crosscorr vector r
R = zeros(P*M, P*M);
r = zeros(P*M, 1);

% FFTs of predictors
Xfft = fft(Xmat, nfft);
Yfft = fft(y, nfft);

% Loop over predictors for crosscorr with y
for p = 1:P
    S_yxp = Yfft .* conj(Xfft(:, p));
    r_full = ifft(S_yxp, nfft, 'symmetric');
    r((p-1)*M + (1:M)) = r_full(1:M);
end

% Loop over predictor pairs for block Toeplitz R
for p = 1:P
    for q = 1:P
        S_xpxq = Xfft(:, p) .* conj(Xfft(:, q));
        rxx_full = ifft(S_xpxq, nfft, 'symmetric');
        rxx = rxx_full(1:M);
        blk = toeplitz(rxx);
        R((p-1)*M + (1:M), (q-1)*M + (1:M)) = blk;
    end
end

% Ridge regularization
lam = alpha * median(abs(R(:)));
R = R + lam * eye(P*M);

% Solve for all kernels at once
H_all = R \ r;    % length P*M vector
H_all = reshape(H_all, [M, P]);  % each column is one kernel

% Time lags
tlag = (0:M-1)' / fs;
end

function [perf,ex,IRF,pred] = perf_multi(X,Y,fs,l)

IRF = fir_multi_kernel(X,Y,fs,l);

pred = zeros(size(Y));
ex = zeros(size(IRF,2),1);

for p = 1:size(IRF,2)
    pred = pred+filter(IRF(:,p),1,X(:,p));
    ex(p) = f_corr(filter(IRF(:,p),1,X(:,p)),Y,1);
end

perf = corrcoef(Y,pred);
perf = perf(2,1);

end
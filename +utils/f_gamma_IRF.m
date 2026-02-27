function [IRF, r, pred, params] = f_gamma_IRF(X, Y, win)

N = size(X, 2);
T = size(X, 1);

X = X - mean(X);
Y = Y - mean(Y);

L = sum(abs(win)) + 1;
t = win(1):win(2);

t_rows = (1 - win(1)) : (T + L - 1 - win(2));
valid_rows = (win(2)) : (T + win(1));

% define minimization problem
obj_fun = @(params) gamma_multi_err_shift(params, X, Y, t, L, t_rows, valid_rows);

params0 = repmat([1, 3, 0], 1, N);  % t0 = 0

% bounds
lb = repmat([-2, 0.5, win(1)], 1, N);  % allow negative t0
ub = repmat([2, 10, win(2)], 1, N);

% options: silent
opts = optimoptions('lsqnonlin', 'Display','off', ...
    'MaxIterations',200,'MaxFunctionEvaluations',1000);

% fit
params_hat = lsqnonlin(obj_fun, params0, lb, ub, opts);

IRF = [];
for i = 1:N
    idx = 3 * (i - 1) + 1;
    IRF = [IRF; gamma_IRF_shift(t, params_hat(idx), params_hat(idx + 1), params_hat(idx + 2))];
end
IRF = IRF';

pred = apply_gamma_multi(params_hat, X, Y, t, L, t_rows);
params = reshape(params_hat,3,N);

r = f_corr(Y, pred, 1);

%% additional functions

function h = gamma_IRF_shift(t, A, tau, t0)
    t = t - t0;

    D = (t / tau).^3 .* exp(-t / tau);
    D(t < 0) = 0;
    
    h = A * D;
end

%% --- multi-predictor error for shifted gamma ---
function err = gamma_multi_err_shift(params, X, Y, t, L, t_rows, valid_rows)
N = size(X,2);

pred = zeros(numel(find(valid_rows)),1);

for i = 1:N
    idx = 3 * (i - 1) + 1;
    
    % gather params
    A = params(idx);
    tau = params(idx+1);
    t0 = params(idx+2);

    h = gamma_IRF_shift(t, A, tau, t0);

    Xi = convmtx(X(:,i), L);
    Xi = Xi(t_rows, :);
    Xi = Xi(valid_rows, :);
    
    pred = pred + Xi * h(:);
end

err = Y(valid_rows) - pred;
end

function pred = apply_gamma_multi(params, X, Y, t, L, t_rows)

N = size(X,2);
pred = zeros(size(Y));

for i = 1:N
    idx = 3 * (i - 1) + 1;
    
    % gather params
    A = params(idx);
    tau = params(idx+1);
    t0 = params(idx+2);

    h = gamma_IRF_shift(t, A, tau, t0);

    Xi = convmtx(X(:,i), L);
    Xi = Xi(t_rows, :);
    
    pred = pred + Xi * h(:);
end

end

end
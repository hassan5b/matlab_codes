%% ===== Four-approach comparison + fixed-T sweep + full-run metrics =====
% Approaches included in this run:
% 1) Original SD baseline, c = c_TI only, T = 0.1
% 2) Igarashi-style TV reciprocal HOCBF-ZOH baseline, T = 0.1
% 3) Corrected TV SD-HOCBF, fixed T = 0.020
% 4) Corrected TV SD-HOCBF, fixed T = 0.100
% 5) Corrected TV SD-HOCBF, fixed T = 0.120
% 6) Corrected TV SD-HOCBF, variable T
%
% Required external dependencies:
% - TIRA and your System_description case system_choice = 99
% - Optimization Toolbox: quadprog, fmincon
%
% Output figures:
% - Full trajectory comparison over the complete rollout
% - Barrier comparison over the complete rollout
% - One snapshot at t = 5.3 with color + line style only
% - Thin black disturbed/measured paths Z overlaid on the snapshot
%
% Output tables printed in the command window:
% - updates, average sampling interval, control effort, path length,
%   RMS trajectory deviation versus the fixed T = 0.020 TV SD-HOCBF reference
% - update reduction of variable-T versus fixed T = 0.020, 0.100, and 0.120

clear; clc; close all; rng(0)

%% ===================== Shared scenario =====================
T_main  = 0.1;      % original baselines and variable-T initial value
T_mpc   = 0.1;      % MPC discretization kept fixed, as in your previous script
t_snap  = 5.3;      % snapshot time only, not the metrics horizon

% Full-run horizon. Each controller runs until it reaches the goal or this
% maximum physical time is exhausted. Metrics are computed on the full result.
max_sim_time = 90;  % matches the old N = 900, T = 0.1 rollout length

fixed_T_sweep = [0.020, 0.100, 0.120];

obs.c_fun     = @(t) [32 + 25*cos(0.3*t); 25];
obs.cdot_fun  = @(t) [-25*0.3*sin(0.3*t); 0];
obs.cddot_fun = @(t) [-25*(0.3^2)*cos(0.3*t); 0];

obs.D_fun     = @(t) 5 + 2*sin(0.4*t);
obs.Ddot_fun  = @(t) 2*0.4*cos(0.4*t);
obs.Dddot_fun = @(t) -2*(0.4^2)*sin(0.4*t);

% Use one shared initial condition so all approaches are comparable.
x0   = [5;25;-pi/2;0.8];
goal = [45;21];

u_min = [-1;-2];
u_max = [ 1; 2];

gamma = 0.3;
epsM  = 0.5;
W     = eye(2);

ig.alpha_h = 1.0;
ig.K       = 0.5;
ig.C       = 0.5;
ig.L_B     = 1e-5;

varpar.beta      = 0.6;
varpar.Tmin      = 0.02;
varpar.Tmax_hw   = 0.12;
varpar.deltaRT   = 0.005;
varpar.g_release = 8.0;
varpar.grow_max  = 1.35;

export_figures = false;
export_tables  = false;

%% ===================== Shared randomness =====================
% All approaches receive the same measurement-noise and disturbance sequences
% by controller-update index. Allocate enough samples for the smallest fixed T.
T_min_for_budget = min([fixed_T_sweep, varpar.Tmin, T_main]);
N_budget = ceil(max_sim_time/T_min_for_budget) + 20;

meas_noise = zeros(4, N_budget+1);
dist_seq   = zeros(2, N_budget);

for k = 1:N_budget+1
    meas_noise(:,k) = [(2*rand-1)*epsM;
                       (2*rand-1)*epsM;
                       0;
                       0];
end

for k = 1:N_budget
    dist_seq(:,k) = gamma*(2*rand(2,1)-1);
end

%% ===================== Run all approaches =====================
global system_choice;
system_choice = 99;

res = repmat(make_result(zeros(4,1), zeros(2,0), zeros(4,1), 0, 0), 1, 6);

fprintf('\nRunning 1/6: original SD baseline, T = %.3f...\n', T_main);
N_main = ceil(max_sim_time/T_main);
tmp = run_original_sd_baseline(x0, goal, obs, T_main, N_main, ...
    u_min, u_max, gamma, epsM, W, meas_noise, dist_seq);
tmp.name = sprintf('Original SD, T=%.3f', T_main);
res(1) = tmp;

fprintf('\nRunning 2/6: Igarashi-style TV reciprocal baseline, T = %.3f...\n', T_main);
tmp = run_igarashi_baseline(x0, goal, obs, T_main, N_main, ...
    u_min, u_max, gamma, epsM, ig, meas_noise, dist_seq);
tmp.name = sprintf('Igarashi TV reciprocal, T=%.3f', T_main);
res(2) = tmp;

fixed_tv_indices = 3:5;
for j = 1:numel(fixed_T_sweep)
    Tj = fixed_T_sweep(j);
    idx_res = fixed_tv_indices(j);
    Nj = ceil(max_sim_time/Tj);
    fprintf('\nRunning %d/6: corrected TV fixed-T SD-HOCBF, T = %.3f...\n', idx_res, Tj);
    tmp = run_tv_fixed_baseline(x0, goal, obs, Tj, T_mpc, Nj, ...
        u_min, u_max, gamma, epsM, W, meas_noise, dist_seq);
    tmp.name = sprintf('TV SD-HOCBF fixed T=%.3f', Tj);
    res(idx_res) = tmp;
end

variable_idx = 6;
N_var = ceil(max_sim_time/varpar.Tmin) + 20;
fprintf('\nRunning 6/6: corrected TV variable-T SD-HOCBF...\n');
tmp = run_tv_variable_baseline(x0, goal, obs, T_main, T_mpc, N_var, ...
    u_min, u_max, gamma, epsM, W, varpar, meas_noise, dist_seq);
tmp.name = 'TV SD-HOCBF variable T';
res(variable_idx) = tmp;

% Do not truncate here: the tables below use each complete rollout.

%% ===================== Metrics and command-window report =====================
ref_idx = 3;  % fixed T = 0.020 is the high-rate TV SD-HOCBF reference
[summary_table, reduction_table] = build_comparison_tables(res, ref_idx, fixed_tv_indices, variable_idx, goal);

fprintf('\n==================== SUMMARY METRICS ====================\n');
disp(summary_table);

fprintf('\n===== VARIABLE-T UPDATE REDUCTION VS FIXED-T TV SD-HOCBF =====\n');
disp(reduction_table);

if export_tables
    writetable(summary_table, 'summary_metrics_four_approaches_plus_fixedT_sweep.csv');
    writetable(reduction_table, 'variableT_update_reduction_vs_fixedT.csv');
end

%% ===================== Styling =====================
colors = [0.000 0.447 0.741;   % blue
          0.850 0.325 0.098;   % orange
          0.466 0.674 0.188;   % green
          0.494 0.184 0.556;   % purple
          0.635 0.078 0.184;   % dark red
          0.301 0.745 0.933];  % light blue

line_styles = {'-', '--', ':', '-.', '-', '--'};

%% ===================== Full path comparison =====================
fig_full = figure('Color','w'); hold on; axis equal; grid on
xlabel('x [m]'); ylabel('y [m]');
title('Trajectory comparison over complete rollout');

for i = 1:numel(res)
    plot(res(i).X(1,:), res(i).X(2,:), ...
        'LineStyle', line_styles{i}, 'Color', colors(i,:), 'LineWidth', 2.0, ...
        'DisplayName', res(i).name);
end

plot(goal(1), goal(2), 'kx', 'MarkerSize', 10, 'LineWidth', 2.0, 'HandleVisibility','off');

th = linspace(0,2*pi,250);
max_t_plot = max(arrayfun(@(r) r.time(end), res));
for tt = linspace(0, max_t_plot, 5)
    c = obs.c_fun(tt);
    D = obs.D_fun(tt);
    plot(c(1)+D*cos(th), c(2)+D*sin(th), 'k--', 'LineWidth', 0.9, 'HandleVisibility','off');
end
legend('Location','eastoutside');
xlim([0,60]); ylim([0,50]);

if export_figures
    exportgraphics(fig_full, 'trajectory_comparison_with_fixedT_sweep.png', 'Resolution', 300);
end

%% ===================== Barrier comparison =====================
fig_h = figure('Color','w'); hold on; grid on
for i = 1:numel(res)
    plot(res(i).time, res(i).h, ...
        'LineStyle', line_styles{i}, 'Color', colors(i,:), 'LineWidth', 1.8, ...
        'DisplayName', res(i).name);
end
yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility','off');
xlabel('time [s]'); ylabel('h(x,t)');
title('Barrier comparison over complete rollout');
legend('Location','eastoutside');

if export_figures
    exportgraphics(fig_h, 'barrier_comparison_with_fixedT_sweep.png', 'Resolution', 300);
end

%% ===================== Single snapshot at t = 5.3 =====================
for i = 1:numel(res)
    snaps(i) = snapshot_at_time(res(i), t_snap, obs);
end

fig_snap = figure('Color','w'); hold on; axis equal; grid on
xlabel('x [m]'); ylabel('y [m]');
title(sprintf('Snapshot at t = %.1f s, color + line style', t_snap));

c = obs.c_fun(t_snap);
D = obs.D_fun(t_snap);
plot(c(1)+D*cos(th), c(2)+D*sin(th), 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');
plot(goal(1), goal(2), 'kx', 'MarkerSize', 10, 'LineWidth', 2.0, 'HandleVisibility','off');

snap_handles = gobjects(1, numel(res));
for i = 1:numel(res)
    % Thin black disturbed/measured path z for this approach.
    pz = snaps(i).z_path;
    if ~isempty(pz)
        plot(pz(1,:), pz(2,:), ...
            'LineStyle', line_styles{i}, ...
            'Color', [0 0 0], ...
            'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
    end

    % Main path: color + line style only. No markers, no robot shapes.
    p = snaps(i).path;
    snap_handles(i) = plot(p(1,:), p(2,:), ...
        'LineStyle', line_styles{i}, ...
        'Color', colors(i,:), ...
        'LineWidth', 2.2, ...
        'DisplayName', res(i).name);
end

legend(snap_handles, 'Location', 'eastoutside');
xlim([0,60]); ylim([0,50]);

if export_figures
    exportgraphics(fig_snap, 'snapshot_t53_fixedT_sweep_variableT_with_zpaths.png', 'Resolution', 300);
end

%% ========================================================================
%% Method 1: original robust SD baseline, c = c_TI only
%% ========================================================================
function res = run_original_sd_baseline(x0, goal, obs, T, N, ...
    u_min, u_max, gamma, epsM, W, meas_noise, dist_seq)

    global system_choice;
    system_choice = 99;

    X = zeros(4,N+1); X(:,1)=x0;
    U = zeros(2,N);
    Z = nan(4,N+1);
    h = nan(1,N+1);

    x_true = x0;
    opts_qp = optimoptions('quadprog','Display','off');

    for k = 1:N
        t = (k-1)*T;

        xhat = x_true + meas_noise(:,k);
        Z(:,k) = xhat;

        u_perf = Kperf_MPC(xhat, goal, u_min, u_max, T);
        [a_hat, c_hat] = fv_affine_at_TV_original(xhat, t, obs);

        Z_lo = xhat + [-epsM; -epsM; 0; 0];
        Z_hi = xhat + [ epsM;  epsM; 0; 0];
        p_lo = [u_min; -gamma; -gamma];
        p_hi = [u_max; +gamma; +gamma];

        [R_lo, R_hi] = TIRA([0, T], Z_lo, Z_hi, p_lo, p_hi);
        Mk = compute_margin_sup_TV_original(R_lo, R_hi, xhat, t, obs, u_min, u_max, gamma);

        A = -a_hat;
        b = c_hat - Mk;
        H = 2*W;
        f = -2*W*u_perf;
        Aqp = [A; eye(2); -eye(2)];
        bqp = [b; u_max; -u_min];

        [u,~,flag] = quadprog(H,f,Aqp,bqp,[],[],[],[],[],opts_qp);
        if isempty(u) || flag <= 0
            u = clip(u_perf, u_min, u_max);
        end
        U(:,k) = u;

        x_true = propagate_true_plant(x_true, u, dist_seq(:,k), T);
        X(:,k+1) = x_true;
        Z(:,k+1) = x_true + meas_noise(:,k+1);
        h(k+1) = h_tv(x_true, k*T, obs);

        if norm(X(1:2,k+1)-goal) < 0.5
            X = X(:,1:k+1); U = U(:,1:k); Z = Z(:,1:k+1); h = h(1:k+1);
            break;
        end
    end

    time = 0:T:T*(size(X,2)-1);
    for j = 1:size(X,2)
        h(j) = h_tv(X(:,j), time(j), obs);
    end

    res = make_result(X, U, Z, h, time);
    res.T_log = T*ones(1,size(U,2));
end

%% ========================================================================
%% Method 2: Igarashi-style TV reciprocal HOCBF-ZOH baseline
%% ========================================================================
function res = run_igarashi_baseline(x0, goal, obs, T, N, ...
    u_min, u_max, gamma, epsM, ig, meas_noise, dist_seq)

    X = zeros(4,N+1); X(:,1)=x0;
    Z = nan(4,N+1);
    U_total = zeros(2,N);
    U_delta = zeros(2,N);

    h_log    = nan(1,N+1);
    psi1_log = nan(1,N);
    B_log    = nan(1,N);
    I_log    = nan(1,N);
    J_log    = nan(1,N);

    x_true = x0;

    for k = 1:N
        t = (k-1)*T;
        xhat = x_true + meas_noise(:,k);
        Z(:,k) = xhat;

        u_h = Kperf_MPC(xhat, goal, u_min, u_max, T);
        [u_delta, u, info] = igarashi_tv_reciprocal_hocbf( ...
            xhat, t, u_h, obs, u_min, u_max, ig.alpha_h, ig.K, ig.C, ig.L_B);

        U_delta(:,k) = u_delta;
        U_total(:,k) = u;
        psi1_log(k) = info.psi1;
        B_log(k)    = info.B;
        I_log(k)    = info.I;
        J_log(k)    = info.J;

        x_true = propagate_true_plant(x_true, u, dist_seq(:,k), T);
        X(:,k+1) = x_true;
        Z(:,k+1) = x_true + meas_noise(:,k+1);
        h_log(k+1) = h_tv(X(:,k+1), k*T, obs);

        if norm(X(1:2,k+1)-goal) < 0.5
            X = X(:,1:k+1); Z = Z(:,1:k+1);
            U_total = U_total(:,1:k); U_delta = U_delta(:,1:k);
            h_log = h_log(1:k+1); psi1_log = psi1_log(1:k);
            B_log = B_log(1:k); I_log = I_log(1:k); J_log = J_log(1:k);
            break;
        end
    end

    time = 0:T:T*(size(X,2)-1);
    for j = 1:size(X,2)
        h_log(j) = h_tv(X(:,j), time(j), obs);
    end

    res = make_result(X, U_total, Z, h_log, time);
    res.U_delta = U_delta;
    res.psi1 = psi1_log;
    res.B = B_log;
    res.I = I_log;
    res.J = J_log;
    res.T_log = T*ones(1,size(U_total,2));
end

%% ========================================================================
%% Method 3: corrected TV SD-HOCBF, fixed T
%% ========================================================================
function res = run_tv_fixed_baseline(x0, goal, obs, T, T_mpc, N, ...
    u_min, u_max, gamma, epsM, W, meas_noise, dist_seq)

    global system_choice;
    system_choice = 99;

    X = zeros(4,N+1); X(:,1)=x0;
    U = zeros(2,N);
    Z = nan(4,N+1);
    h = nan(1,N+1);

    x_true = x0;
    opts_qp = optimoptions('quadprog','Display','off');

    for k = 1:N
        t = (k-1)*T;
        xhat = x_true + meas_noise(:,k);
        Z(:,k) = xhat;

        u_perf = Kperf_MPC(xhat, goal, u_min, u_max, T_mpc);
        [a_hat, c_hat] = fv_affine_at_TV_corrected(xhat, t, obs);
        lFt = tv_time_lipschitz(xhat, t, obs);

        Z_lo = xhat + [-epsM; -epsM; 0; 0];
        Z_hi = xhat + [ epsM;  epsM; 0; 0];
        p_lo = [u_min; -gamma; -gamma];
        p_hi = [u_max; +gamma; +gamma];

        [R_lo, R_hi] = TIRA([0, T], Z_lo, Z_hi, p_lo, p_hi);
        Mk = compute_margin_sup_TV_corrected(R_lo, R_hi, xhat, t, obs, u_min, u_max, gamma);

        A = -a_hat;
        b = c_hat - Mk - lFt*T;
        H = 2*W;
        f = -2*W*u_perf;
        Aqp = [A; eye(2); -eye(2)];
        bqp = [b; u_max; -u_min];

        [u,~,flag] = quadprog(H,f,Aqp,bqp,[],[],[],[],[],opts_qp);
        if isempty(u) || flag <= 0
            u = clip(u_perf, u_min, u_max);
        end
        U(:,k) = u;

        x_true = propagate_true_plant(x_true, u, dist_seq(:,k), T);
        X(:,k+1) = x_true;
        Z(:,k+1) = x_true + meas_noise(:,k+1);
        h(k+1) = h_tv(x_true, k*T, obs);

        if norm(X(1:2,k+1)-goal) < 0.5
            X = X(:,1:k+1); U = U(:,1:k); Z = Z(:,1:k+1); h = h(1:k+1);
            break;
        end
    end

    time = 0:T:T*(size(X,2)-1);
    for j = 1:size(X,2)
        h(j) = h_tv(X(:,j), time(j), obs);
    end

    res = make_result(X, U, Z, h, time);
    res.T_log = T*ones(1,size(U,2));
end

%% ========================================================================
%% Method 4: corrected TV SD-HOCBF, variable T
%% ========================================================================
function res = run_tv_variable_baseline(x0, goal, obs, T_init, T_mpc, N, ...
    u_min, u_max, gamma, epsM, W, varpar, meas_noise, dist_seq)

    global system_choice;
    system_choice = 99;

    T = T_init;
    X = zeros(4,N+1); X(:,1)=x0;
    U = zeros(2,N);
    Z = nan(4,N+1);
    h = nan(1,N+1);
    thist = nan(1,N+1); thist(1)=0;
    T_log = nan(1,N);
    Tnext_log = nan(1,N);

    x_true = x0;
    t = 0;
    opts_qp = optimoptions('quadprog','Display','off');

    for k = 1:N
        xhat = x_true + meas_noise(:,k);
        Z(:,k) = xhat;

        u_perf = Kperf_MPC(xhat, goal, u_min, u_max, T_mpc);
        [a_hat, c_hat] = fv_affine_at_TV_corrected(xhat, t, obs);
        lFt = tv_time_lipschitz(xhat, t, obs);

        Z_lo = xhat + [-epsM; -epsM; 0; 0];
        Z_hi = xhat + [ epsM;  epsM; 0; 0];
        p_lo = [u_min; -gamma; -gamma];
        p_hi = [u_max; +gamma; +gamma];

        [R_lo, R_hi] = TIRA([0, T], Z_lo, Z_hi, p_lo, p_hi);
        Mk = compute_margin_sup_TV_corrected(R_lo, R_hi, xhat, t, obs, u_min, u_max, gamma);

        A = -a_hat;
        b = c_hat - Mk - lFt*T;
        H = 2*W;
        f = -2*W*u_perf;
        Aqp = [A; eye(2); -eye(2)];
        bqp = [b; u_max; -u_min];

        t_qp = tic;
        [u,~,flag] = quadprog(H,f,Aqp,bqp,[],[],[],[],[],opts_qp);
        solve_time_k = toc(t_qp);

        if isempty(u) || flag <= 0
            u = clip(u_perf, u_min, u_max);
        end
        U(:,k) = u;

        g_now = a_hat*u + c_hat - Mk - lFt*T;
        if g_now > varpar.g_release
            T_next = min(varpar.Tmax_hw, max(varpar.Tmin, varpar.grow_max*T));
        else
            if feasible_at_TV_corrected(varpar.Tmin, xhat, u, a_hat, c_hat, lFt, t, obs, ...
                    u_min, u_max, gamma, epsM)
                T_lo = varpar.Tmin;
                T_hi = varpar.Tmax_hw;
                for it = 1:8
                    T_mid = 0.5*(T_lo + T_hi);
                    if feasible_at_TV_corrected(T_mid, xhat, u, a_hat, c_hat, lFt, t, obs, ...
                            u_min, u_max, gamma, epsM)
                        T_lo = T_mid;
                    else
                        T_hi = T_mid;
                    end
                end
                T_next = varpar.beta*T_lo;
            else
                T_next = varpar.Tmin;
            end
            T_next = min(max(T_next, varpar.Tmin), varpar.Tmax_hw);
        end

        T_next = max(T_next, solve_time_k + varpar.deltaRT);
        T_next = min(max(T_next, varpar.Tmin), varpar.Tmax_hw);

        x_true = propagate_true_plant(x_true, u, dist_seq(:,k), T);
        X(:,k+1) = x_true;
        Z(:,k+1) = x_true + meas_noise(:,k+1);

        T_log(k) = T;
        Tnext_log(k) = T_next;
        t = t + T;
        thist(k+1) = t;
        h(k+1) = h_tv(x_true, t, obs);

        T = T_next;

        if norm(X(1:2,k+1)-goal) < 0.5
            X = X(:,1:k+1); U = U(:,1:k); Z = Z(:,1:k+1);
            h = h(1:k+1); thist = thist(1:k+1);
            T_log = T_log(1:k); Tnext_log = Tnext_log(1:k);
            break;
        end
    end

    for j = 1:size(X,2)
        h(j) = h_tv(X(:,j), thist(j), obs);
    end

    res = make_result(X, U, Z, h, thist);
    res.T_log = T_log;
    res.Tnext_log = Tnext_log;
end

%% ========================================================================
%% Affine pieces and margins
%% ========================================================================
function [a, c] = fv_affine_at_TV_original(x, t, obs)
% Original baseline version: compute TV terms, but use only c_TI.
    [a, c_TI, ~] = fv_affine_parts(x, t, obs);
    c = c_TI;
end

function [a, c] = fv_affine_at_TV_corrected(x, t, obs)
% Corrected TV version: c = c_TI + explicit time-varying correction.
    [a, c_TI, c_TVcorr] = fv_affine_parts(x, t, obs);
    c = c_TI + c_TVcorr;
end

function [a, c_TI, c_TVcorr] = fv_affine_parts(x, t, obs)
    xx=x(1); yy=x(2); th=x(3); v=x(4);

    co   = obs.c_fun(t);      Dx  = co(1);      Dy  = co(2);
    cDot = obs.cdot_fun(t);   Dxdt = cDot(1);   Dydt = cDot(2);
    cDD  = obs.cddot_fun(t);  Dxd2 = cDD(1);    Dyd2 = cDD(2);
    D    = obs.D_fun(t);
    Dd   = obs.Ddot_fun(t);
    Ddd  = obs.Dddot_fun(t);

    dx = xx - Dx; dy = yy - Dy;
    phi = dx*cos(th) + dy*sin(th);

    h    = dx*dx + dy*dy - D*D;
    Lfh  = 2*v*phi;
    Lf2h = 2*v*v;

    a1 = 2*v*(-dx*sin(th) + dy*cos(th));
    a2 = 2*phi;
    a  = [a1, a2];

    h_t   = -2*(dx*Dxdt + dy*Dydt) - 2*D*Dd;
    Lf_ht = -2*v*(Dxdt*cos(th) + Dydt*sin(th));
    dtLfh = -2*v*(Dxdt*cos(th) + Dydt*sin(th));
    h_tt  = -2*(dx*Dxd2 + dy*Dyd2) + 2*(Dxdt*Dxdt + Dydt*Dydt) ...
            - 2*(Dd*Dd) - 2*D*Ddd;

    c_TI = Lf2h + 2*Lfh + 2*h;
    c_TVcorr = Lf_ht + dtLfh + h_tt + 2*h_t;
end

function Mk = compute_margin_sup_TV_original(R_lo, R_hi, xhat, t, obs, umin, umax, gamma)
    Mk = compute_margin_sup_generic(R_lo, R_hi, xhat, t, obs, umin, umax, gamma, false);
end

function Mk = compute_margin_sup_TV_corrected(R_lo, R_hi, xhat, t, obs, umin, umax, gamma)
    Mk = compute_margin_sup_generic(R_lo, R_hi, xhat, t, obs, umin, umax, gamma, true);
end

function Mk = compute_margin_sup_generic(R_lo, R_hi, xhat, t, obs, umin, umax, gamma, use_tv_correction)
    if use_tv_correction
        [a_hat, c_hat] = fv_affine_at_TV_corrected(xhat, t, obs);
    else
        [a_hat, c_hat] = fv_affine_at_TV_original(xhat, t, obs);
    end

    support_u = @(Da) max(Da(1)*umin(1), Da(1)*umax(1)) + ...
                      max(Da(2)*umin(2), Da(2)*umax(2));

    xs = linspace(0,1,4);     % 4^4 = 256 samples, including corners
    [X1,X2,X3,X4] = ndgrid(xs,xs,xs,xs);
    gridPts = [X1(:) X2(:) X3(:) X4(:)].';

    c_now = obs.c_fun(t);
    sup_all = -inf;

    for ii = 1:size(gridPts,2)
        r = gridPts(:,ii);
        x = R_lo + r.*(R_hi - R_lo);

        if use_tv_correction
            [a, c] = fv_affine_at_TV_corrected(x, t, obs);
        else
            [a, c] = fv_affine_at_TV_original(x, t, obs);
        end

        Dc = c - c_hat;
        Da = a - a_hat;
        sup_u = support_u(Da);

        v = x(4); th = x(3);
        dx = x(1) - c_now(1);
        dy = x(2) - c_now(2);
        c1 = 4*(v*cos(th) + dx);
        c2 = 4*(v*sin(th) + dy);
        sup_d = gamma*(abs(c1) + abs(c2)) + 4*gamma^2;

        sup_all = max(sup_all, Dc + sup_u + sup_d);
    end

    Mk = sup_all;
end

function Mk_u = compute_margin_sup_TV_given_u_corrected(R_lo, R_hi, xhat, t, obs, u, gamma)
    [a_hat, c_hat] = fv_affine_at_TV_corrected(xhat, t, obs);

    xs = linspace(0,1,3);
    [X1,X2,X3,X4] = ndgrid(xs,xs,xs,xs);
    gridPts = [X1(:) X2(:) X3(:) X4(:)].';

    c_now = obs.c_fun(t);
    sup_all = -inf;

    for ii = 1:size(gridPts,2)
        r = gridPts(:,ii);
        x = R_lo + r.*(R_hi - R_lo);

        [a, c] = fv_affine_at_TV_corrected(x, t, obs);
        Dc = c - c_hat;
        Du = (a - a_hat)*u;

        v = x(4); th = x(3);
        dx = x(1) - c_now(1);
        dy = x(2) - c_now(2);
        c1 = 4*(v*cos(th) + dx);
        c2 = 4*(v*sin(th) + dy);
        sup_d = gamma*(abs(c1) + abs(c2)) + 4*gamma^2;

        sup_all = max(sup_all, Dc + Du + sup_d);
    end

    Mk_u = sup_all;
end

function ok = feasible_at_TV_corrected(Tcand, xhat, u, a_hat, c_hat, lFt, t, obs, ...
                                      u_min, u_max, gamma, epsM)
    Z_lo = xhat + [-epsM; -epsM; 0; 0];
    Z_hi = xhat + [ epsM;  epsM; 0; 0];
    p_lo = [u_min; -gamma; -gamma];
    p_hi = [u_max; +gamma; +gamma];

    [R_lo, R_hi] = TIRA([0, Tcand], Z_lo, Z_hi, p_lo, p_hi);
    Mk_T = compute_margin_sup_TV_given_u_corrected(R_lo, R_hi, xhat, t, obs, u, gamma);
    g_T = a_hat*u + c_hat - Mk_T - lFt*Tcand;
    ok = (g_T >= 0);
end

function lFt = tv_time_lipschitz(xhat, t, obs)
    dt = 1e-3;
    [~, c0] = fv_affine_at_TV_corrected(xhat, t,    obs);
    [~, c1] = fv_affine_at_TV_corrected(xhat, t+dt, obs);
    lFt = abs(c1 - c0)/dt;
end

%% ========================================================================
%% Igarashi helpers
%% ========================================================================
function [u_delta, u_total, info] = igarashi_tv_reciprocal_hocbf( ...
    x, t, u_h, obs, u_min, u_max, alpha_h, K, C, L_B)

    [B, psi1, gradB, Bt, valid] = reciprocal_lifted_barrier_data(x, t, obs, alpha_h, L_B);

    f = drift_unicycle(x);
    G = input_matrix_unicycle();

    LfB = gradB.'*f;
    LgB = gradB.'*G;

    I = LfB + LgB*u_h + Bt;
    J = K*B + C;

    if valid && I > J && norm(LgB) > 1e-10
        u_delta_exact = -((I - J)/(norm(LgB)^2))*LgB.';
    else
        u_delta_exact = zeros(2,1);
    end

    H = eye(2);
    fqp = zeros(2,1);
    Aineq = LgB;
    bineq = J - I;
    lb = u_min - u_h;
    ub = u_max - u_h;

    opts = optimoptions('quadprog','Display','off');
    [u_delta_qp,~,flag] = quadprog(H, fqp, Aineq, bineq, [], [], lb, ub, [], opts);

    if valid && ~isempty(u_delta_qp) && flag > 0
        u_delta = u_delta_qp;
    else
        u_delta = u_delta_exact;
    end

    u_total = clip(u_h + u_delta, u_min, u_max);
    u_delta = u_total - u_h;

    info.B = B;
    info.psi1 = psi1;
    info.I = I;
    info.J = J;
    info.valid = valid;
end

function [B, psi1, gradB, Bt, valid] = reciprocal_lifted_barrier_data(x, t, obs, alpha_h, L_B)
    [psi1, gradPsi1, psi1_t] = psi1_data(x, t, obs, alpha_h);

    eps_psi = 1e-8;
    valid = psi1 > eps_psi;
    psi_safe = max(psi1, eps_psi);

    B = 1/psi_safe + L_B*(x.'*x);

    if valid
        gradB = -(1/psi1^2)*gradPsi1 + 2*L_B*x;
        Bt = -(1/psi1^2)*psi1_t;
    else
        gradB = zeros(4,1);
        Bt = 0;
    end
end

function [psi1, gradPsi1, psi1_t] = psi1_data(x, t, obs, alpha_h)
    px = x(1); py = x(2); th = x(3); v = x(4);

    c    = obs.c_fun(t);
    cdot = obs.cdot_fun(t);
    cdd  = obs.cddot_fun(t);
    D    = obs.D_fun(t);
    Dd   = obs.Ddot_fun(t);
    Ddd  = obs.Dddot_fun(t);

    r = [px - c(1); py - c(2)];
    q = [cos(th); sin(th)];
    qperp = [-sin(th); cos(th)];

    h = r.'*r - D^2;
    phi = r.'*q;
    Lfh = 2*v*phi;
    h_t = -2*(r.'*cdot) - 2*D*Dd;

    psi1 = Lfh + h_t + alpha_h*h;

    grad_Lfh = [2*v*cos(th);
                2*v*sin(th);
                2*v*(r.'*qperp);
                2*phi];

    grad_ht = [-2*cdot(1);
               -2*cdot(2);
                0;
                0];

    grad_h = [2*r(1);
              2*r(2);
              0;
              0];

    gradPsi1 = grad_Lfh + grad_ht + alpha_h*grad_h;

    Lfh_t = -2*v*(cdot.'*q);
    ht_t = 2*(cdot.'*cdot) - 2*(r.'*cdd) - 2*(Dd^2) - 2*D*Ddd;
    psi1_t = Lfh_t + ht_t + alpha_h*h_t;
end

%% ========================================================================
%% Dynamics, MPC, and result helpers
%% ========================================================================
function res = make_result(X, U, Z, h, time)
    res = struct();
    res.name = '';
    res.X = X;
    res.U = U;
    res.Z = Z;
    res.h = h;
    res.time = time;
    res.psi1 = [];
    res.B = [];
    res.I = [];
    res.J = [];
    res.U_delta = [];
    res.T_log = [];
    res.Tnext_log = [];
end

function h = h_tv(x, t, obs)
    c = obs.c_fun(t);
    D = obs.D_fun(t);
    r = x(1:2) - c;
    h = r.'*r - D^2;
end

function x_next = propagate_true_plant(x, u, d, T)
    G = input_matrix_unicycle();
    P = [1 0;
         0 1;
         0 0;
         0 0];

    dyn = @(xx) drift_unicycle(xx) + G*u + P*d;

    k1 = dyn(x);
    k2 = dyn(x + 0.5*T*k1);
    k3 = dyn(x + 0.5*T*k2);
    k4 = dyn(x + T*k3);

    x_next = x + (T/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function f = drift_unicycle(x)
    f = [x(4)*cos(x(3));
         x(4)*sin(x(3));
         0;
         0];
end

function G = input_matrix_unicycle()
    G = [0 0;
         0 0;
         1 0;
         0 1];
end

function xnext = rk4_unicycle(x, u, T)
    dyn = @(xx) drift_unicycle(xx) + input_matrix_unicycle()*u;
    k1 = dyn(x);
    k2 = dyn(x + 0.5*T*k1);
    k3 = dyn(x + 0.5*T*k2);
    k4 = dyn(x + T*k3);
    xnext = x + (T/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function u0 = Kperf_MPC(x0, goal, umin, umax, T)
    Np = 50;
    Qp = 5.0;
    Qf = 400;
    Rw = diag([0.1, 0.05]);

    U0 = zeros(2*Np,1);
    lb = repmat([umin(1); umin(2)], Np, 1);
    ub = repmat([umax(1); umax(2)], Np, 1);

    function J = cost(U)
        x = x0;
        J = 0;
        for kk = 1:Np
            u = U(2*kk-1:2*kk);
            x = rk4_unicycle(x, u, T);
            pos_err = x(1:2) - goal(:);
            J = J + Qp*(pos_err.'*pos_err) + u.'*Rw*u;
        end
        term_err = x(1:2) - goal(:);
        J = J + Qf*(term_err.'*term_err);
    end

    opts = optimoptions('fmincon','Display','off','MaxIterations',100,'Algorithm','sqp');
    Uopt = fmincon(@cost, U0, [], [], [], [], lb, ub, [], opts);
    u0 = Uopt(1:2);
end

function y = clip(x, lo, hi)
    y = min(max(x, lo), hi);
end

%% ========================================================================
%% Metric, snapshot, and plotting helpers
%% ========================================================================
function res_out = truncate_result_at_time(res, t_stop, obs)
% Truncate a result structure to t_stop. If t_stop falls inside a held-input
% interval, the final state and measurement are linearly interpolated and the
% last sampling interval is shortened for fair control-effort accounting.
    res_out = res;
    time = res.time(:).';

    if isempty(time) || time(end) <= t_stop + 1e-12
        return;
    end

    X = res.X;
    Z = res.Z;
    U = res.U;
    K = size(U,2);

    idx_pre = find(time <= t_stop, 1, 'last');
    if isempty(idx_pre)
        idx_pre = 1;
    end

    exact_hit = abs(time(idx_pre) - t_stop) <= 1e-12;

    if exact_hit
        new_time = time(1:idx_pre);
        new_X = X(:,1:idx_pre);
        new_Z = Z(:,1:idx_pre);
        new_K = max(0, idx_pre - 1);
        new_U = U(:,1:min(new_K,K));
    else
        x_stop = interp1(time.', X.', t_stop, 'linear').';
        z_stop = interp1(time.', Z.', t_stop, 'linear').';

        new_time = [time(1:idx_pre), t_stop];
        new_X = [X(:,1:idx_pre), x_stop];
        new_Z = [Z(:,1:idx_pre), z_stop];

        % The update that started at time(idx_pre) is active until t_stop.
        new_K = min(idx_pre, K);
        new_U = U(:,1:new_K);
    end

    new_h = zeros(1, numel(new_time));
    for j = 1:numel(new_time)
        new_h(j) = h_tv(new_X(:,j), new_time(j), obs);
    end

    res_out.X = new_X;
    res_out.Z = new_Z;
    res_out.U = new_U;
    res_out.h = new_h;
    res_out.time = new_time;

    if isfield(res_out, 'T_log') && ~isempty(res.T_log)
        Tlog = res.T_log(:).';
        Tlog = Tlog(1:min(new_K, numel(Tlog)));
        if ~exact_hit && ~isempty(Tlog)
            Tlog(end) = t_stop - time(idx_pre);
        end
        res_out.T_log = Tlog;
    end

    if isfield(res_out, 'Tnext_log') && ~isempty(res.Tnext_log)
        Tnext = res.Tnext_log(:).';
        res_out.Tnext_log = Tnext(1:min(new_K, numel(Tnext)));
    end

    if isfield(res_out, 'psi1') && ~isempty(res.psi1)
        res_out.psi1 = res.psi1(1:min(new_K, numel(res.psi1)));
    end
    if isfield(res_out, 'B') && ~isempty(res.B)
        res_out.B = res.B(1:min(new_K, numel(res.B)));
    end
    if isfield(res_out, 'I') && ~isempty(res.I)
        res_out.I = res.I(1:min(new_K, numel(res.I)));
    end
    if isfield(res_out, 'J') && ~isempty(res.J)
        res_out.J = res.J(1:min(new_K, numel(res.J)));
    end
    if isfield(res_out, 'U_delta') && ~isempty(res.U_delta)
        res_out.U_delta = res.U_delta(:,1:min(new_K, size(res.U_delta,2)));
    end
end

function [summary_table, reduction_table] = build_comparison_tables(res, ref_idx, fixed_tv_indices, variable_idx, goal)
    n = numel(res);

    Approach = cell(n,1);
    FinalTime_s = zeros(n,1);
    Updates = zeros(n,1);
    AvgSampling_s = zeros(n,1);
    ControlEffort_u2s = zeros(n,1);
    PathLength_m = zeros(n,1);
    RMSDevVsFixedT0020_m = zeros(n,1);
    MinH = zeros(n,1);
    FinalDistanceToGoal_m = zeros(n,1);

    ref = res(ref_idx);
    for i = 1:n
        Approach{i} = res(i).name;
        FinalTime_s(i) = res(i).time(end);
        Updates(i) = size(res(i).U,2);
        AvgSampling_s(i) = average_sampling_interval(res(i));
        ControlEffort_u2s(i) = control_effort_metric(res(i));
        PathLength_m(i) = path_length_metric(res(i));
        RMSDevVsFixedT0020_m(i) = rms_position_deviation(res(i), ref, 0.01);
        MinH(i) = min(res(i).h);
        FinalDistanceToGoal_m(i) = norm(res(i).X(1:2,end) - goal(:));
    end

    summary_table = table(Approach, FinalTime_s, Updates, AvgSampling_s, ...
        ControlEffort_u2s, PathLength_m, RMSDevVsFixedT0020_m, MinH, FinalDistanceToGoal_m);

    m = numel(fixed_tv_indices);
    FixedT_s = zeros(m,1);
    FixedFinalTime_s = zeros(m,1);
    FixedUpdates = zeros(m,1);
    VariableFinalTime_s = zeros(m,1);
    VariableAvgT_s = zeros(m,1);
    VariableUpdates = zeros(m,1);
    UpdateReductionPct = zeros(m,1);

    var_updates = size(res(variable_idx).U,2);
    var_avg_T = average_sampling_interval(res(variable_idx));
    var_final_time = res(variable_idx).time(end);

    for j = 1:m
        idx = fixed_tv_indices(j);
        FixedT_s(j) = average_sampling_interval(res(idx));
        FixedFinalTime_s(j) = res(idx).time(end);
        FixedUpdates(j) = size(res(idx).U,2);
        VariableFinalTime_s(j) = var_final_time;
        VariableAvgT_s(j) = var_avg_T;
        VariableUpdates(j) = var_updates;
        UpdateReductionPct(j) = 100*(1 - VariableUpdates(j)/FixedUpdates(j));
    end

    reduction_table = table(FixedT_s, FixedFinalTime_s, FixedUpdates, ...
        VariableFinalTime_s, VariableAvgT_s, VariableUpdates, UpdateReductionPct);
end

function avgT = average_sampling_interval(res)
    K = size(res.U,2);
    if K == 0
        avgT = NaN;
        return;
    end

    if isfield(res, 'T_log') && ~isempty(res.T_log)
        Tlog = res.T_log(:).';
        Tlog = Tlog(1:min(K,numel(Tlog)));
        Tlog = Tlog(~isnan(Tlog));
        avgT = mean(Tlog);
    else
        avgT = (res.time(end) - res.time(1))/K;
    end
end

function effort = control_effort_metric(res)
% Integral-like control effort: sum_k ||u_k||^2 T_k.
    K = size(res.U,2);
    if K == 0
        effort = 0;
        return;
    end

    if isfield(res, 'T_log') && ~isempty(res.T_log)
        Tlog = res.T_log(:).';
        Tlog = Tlog(1:min(K,numel(Tlog)));
    else
        Tlog = diff(res.time);
    end

    Kuse = min(K, numel(Tlog));
    effort = sum(sum(res.U(:,1:Kuse).^2, 1) .* Tlog(1:Kuse));
end

function L = path_length_metric(res)
    if size(res.X,2) < 2
        L = 0;
        return;
    end
    dX = diff(res.X(1:2,:), 1, 2);
    L = sum(sqrt(sum(dX.^2, 1)));
end

function rms_dev = rms_position_deviation(res, ref, dt_grid)
% RMS position deviation over the common time interval.
    t_end = min(res.time(end), ref.time(end));
    if t_end <= 0
        rms_dev = 0;
        return;
    end

    tq = 0:dt_grid:t_end;
    if tq(end) < t_end
        tq = [tq, t_end];
    end

    p = interp1(res.time(:), res.X(1:2,:).', tq(:), 'linear', 'extrap');
    pref = interp1(ref.time(:), ref.X(1:2,:).', tq(:), 'linear', 'extrap');

    e2 = sum((p - pref).^2, 2);
    rms_dev = sqrt(mean(e2));
end

function snap = snapshot_at_time(res, t_snap, obs)
    time = res.time(:).';
    X = res.X;
    Z = [];
    if isfield(res, 'Z') && ~isempty(res.Z)
        Z = res.Z;
    end

    if t_snap <= time(1)
        x_snap = X(:,1);
        if ~isempty(Z), z_snap = Z(:,1); else, z_snap = []; end
        t_eff = time(1);
        idx_pre = 1;
    elseif t_snap >= time(end)
        x_snap = X(:,end);
        if ~isempty(Z), z_snap = Z(:,end); else, z_snap = []; end
        t_eff = time(end);
        idx_pre = numel(time);
    else
        x_snap = interp1(time.', X.', t_snap, 'linear').';
        if ~isempty(Z)
            z_snap = interp1(time.', Z.', t_snap, 'linear').';
        else
            z_snap = [];
        end
        t_eff = t_snap;
        idx_pre = find(time <= t_snap, 1, 'last');
    end

    path = X(:,1:idx_pre);
    if norm(path(:,end) - x_snap) > 1e-12
        path = [path, x_snap];
    end

    z_path = [];
    if ~isempty(Z)
        z_path = Z(:,1:idx_pre);
        if ~isempty(z_snap)
            if any(isnan(z_path(:,end))) || norm(z_path(:,end) - z_snap) > 1e-12
                z_path = [z_path, z_snap];
            end
        end
    end

    snap.t = t_eff;
    snap.x = x_snap;
    snap.z = z_snap;
    snap.h = h_tv(x_snap, t_eff, obs);
    snap.path = path;
    snap.z_path = z_path;
end

%% ===== Compare ORIGINAL baseline codes + snapshot at t = 5.3 =====
% Keeps both methods exactly as originally used:
%
% 1) Robust SD-HOCBF code (original baseline version)
%    -> keeps c = c_TI exactly as in the uploaded script
%
% 2) Igarashi-style TV reciprocal HOCBF-ZOH baseline
%
% Output:
%   - full path comparison
%   - h(x,t) comparison
%   - snapshot at t = 5.3 for each method
%   - overlay snapshot at t = 5.3

clear; clc; close all; rng(0)

%% ===================== Shared scenario =====================
T = 0.1;
N = 900;
t_snap = 5.3;

obs.c_fun     = @(t) [32 + 25*cos(0.3*t); 25];
obs.cdot_fun  = @(t) [-25*0.3*sin(0.3*t); 0];
obs.cddot_fun = @(t) [-25*(0.3^2)*cos(0.3*t); 0];

obs.D_fun     = @(t) 5 + 2*sin(0.4*t);
obs.Ddot_fun  = @(t) 2*0.4*cos(0.4*t);
obs.Dddot_fun = @(t) -2*(0.4^2)*sin(0.4*t);

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

%% ===================== Shared randomness =====================
meas_noise = zeros(4, N+1);
dist_seq   = zeros(2, N);

for k = 1:N+1
    meas_noise(:,k) = [(2*rand-1)*epsM;
                       (2*rand-1)*epsM;
                       0;
                       0];
end

for k = 1:N
    dist_seq(:,k) = gamma*(2*rand(2,1)-1);
end

%% ===================== Run both original methods =====================
fprintf('\nRunning ORIGINAL robust SD-HOCBF baseline...\n');
res_sd = run_original_sd_baseline(x0, goal, obs, T, N, ...
    u_min, u_max, gamma, epsM, W, meas_noise, dist_seq);

fprintf('\nRunning ORIGINAL Igarashi-style baseline...\n');
res_ig = run_original_igarashi_baseline(x0, goal, obs, T, N, ...
    u_min, u_max, gamma, epsM, ig, meas_noise, dist_seq);

%% ===================== Summary =====================
fprintf('\n==================== SUMMARY ====================\n');
fprintf('Original SD baseline:\n');
fprintf('  min h(x,t)     = %.4f\n', min(res_sd.h));
fprintf('  final distance = %.4f\n', norm(res_sd.X(1:2,end)-goal));
fprintf('  updates        = %d\n', size(res_sd.U,2));

fprintf('\nOriginal Igarashi baseline:\n');
fprintf('  min h(x,t)     = %.4f\n', min(res_ig.h));
fprintf('  min psi1       = %.4f\n', min(res_ig.psi1));
fprintf('  final distance = %.4f\n', norm(res_ig.X(1:2,end)-goal));
fprintf('  updates        = %d\n', size(res_ig.U,2));

%% ===================== Full path comparison =====================
figure('Color','w'); hold on; axis equal; grid on
xlabel('x [m]'); ylabel('y [m]');
title('Full trajectory comparison');

plot(res_sd.X(1,:), res_sd.X(2,:), 'b-', 'LineWidth', 2.0);
plot(res_ig.X(1,:), res_ig.X(2,:), 'r-', 'LineWidth', 2.0);

plot(x0(1), x0(2), 'ko', 'MarkerSize', 7, 'LineWidth', 1.5);
plot(goal(1), goal(2), 'gx', 'MarkerSize', 10, 'LineWidth', 2.0);

th = linspace(0,2*pi,250);
snap_times = linspace(0, max(res_sd.time(end), res_ig.time(end)), 7);
for tt = snap_times
    c = obs.c_fun(tt);
    D = obs.D_fun(tt);
    plot(c(1)+D*cos(th), c(2)+D*sin(th), 'k--', 'LineWidth', 1.0);
end

legend('Original SD baseline', ...
       'Original Igarashi baseline', ...
       'Start', 'Goal', 'Obstacle snapshots', ...
       'Location', 'best');
xlim([0,60]); ylim([0,50]);

%% ===================== h(x,t) comparison =====================
figure('Color','w'); hold on; grid on
plot(res_sd.time, res_sd.h, 'b-', 'LineWidth', 1.8);
plot(res_ig.time, res_ig.h, 'r-', 'LineWidth', 1.8);
yline(0, 'k--', 'LineWidth', 1.2);
xlabel('time [s]');
ylabel('h(x,t)');
title('Barrier comparison');
legend('Original SD baseline', 'Original Igarashi baseline', 'h=0', ...
       'Location', 'best');

%% ===================== Snapshot at t = 5.3 =====================
snap_sd = get_snapshot(res_sd, t_snap, T, obs);
snap_ig = get_snapshot(res_ig, t_snap, T, obs);

% --- side-by-side snapshots ---
figure('Color','w');

subplot(1,2,1); hold on; axis equal; grid on
plot_snapshot(res_sd, snap_sd, obs, t_snap, goal, [0 0.45 0.95], ...
    'Original SD baseline at t = 5.3');
xlim([0,60]); ylim([0,50]);

subplot(1,2,2); hold on; axis equal; grid on
plot_snapshot(res_ig, snap_ig, obs, t_snap, goal, [0.85 0.2 0.2], ...
    'Original Igarashi baseline at t = 5.3');
xlim([0,60]); ylim([0,50]);

% --- overlay snapshot ---
figure('Color','w'); hold on; axis equal; grid on
xlabel('x [m]'); ylabel('y [m]');
title('Overlay snapshot at t = 5.3');

% obstacle at t_snap
c = obs.c_fun(t_snap);
D = obs.D_fun(t_snap);
plot(c(1)+D*cos(th), c(2)+D*sin(th), 'k--', 'LineWidth', 1.5);

% paths up to t_snap
plot(res_sd.X(1,1:snap_sd.idx), res_sd.X(2,1:snap_sd.idx), ...
    'b-', 'LineWidth', 2.0);
plot(res_ig.X(1,1:snap_ig.idx), res_ig.X(2,1:snap_ig.idx), ...
    'r-', 'LineWidth', 2.0);

% robot glyphs
[vx_sd, vy_sd] = triangle_pose(snap_sd.x(1), snap_sd.x(2), snap_sd.x(3), 0.9);
patch(vx_sd, vy_sd, [0 0.45 0.95], 'FaceAlpha', 0.45, 'EdgeColor', 'k');

[vx_ig, vy_ig] = triangle_pose(snap_ig.x(1), snap_ig.x(2), snap_ig.x(3), 0.9);
patch(vx_ig, vy_ig, [0.85 0.2 0.2], 'FaceAlpha', 0.45, 'EdgeColor', 'k');

plot(goal(1), goal(2), 'gx', 'MarkerSize', 10, 'LineWidth', 2.0);
legend('Obstacle at t=5.3', ...
       'Original SD baseline path', ...
       'Original Igarashi baseline path', ...
       'Original SD baseline robot', ...
       'Original Igarashi baseline robot', ...
       'Goal', ...
       'Location', 'best');
xlim([0,60]); ylim([0,50]);

%% ===================== Optional: export snapshot figures =====================
% exportgraphics(gcf, 'snapshot_overlay_t53.png', 'Resolution', 300);


%% ========================================================================
%% Method 1: ORIGINAL robust SD baseline (kept exactly as originally used)
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

        z = x_true + meas_noise(:,k);
        xhat = z;
        Z(:,k) = z;

        u_perf = Kperf_MPC(xhat, goal, u_min, u_max, T);

        [a_hat, c_hat] = fv_affine_at_TV_original(xhat, t, obs);

        Z_lo = xhat + [-epsM; -epsM; 0; 0];
        Z_hi = xhat + [ epsM;  epsM; 0; 0];

        p_lo = [u_min; -gamma; -gamma];
        p_hi = [u_max; +gamma; +gamma];

        [R_lo, R_hi] = TIRA([0, T], Z_lo, Z_hi, p_lo, p_hi);

        Mk = compute_margin_sup_TV_original(R_lo, R_hi, xhat, t, obs, u_min, u_max, gamma);

        A = -a_hat;
        b = (c_hat - Mk);

        H = 2*W;
        f = -2*W*u_perf;

        Aqp = [A; eye(2); -eye(2)];
        bqp = [b; u_max; -u_min];

        [u,~,flag] = quadprog(H,f,Aqp,bqp,[],[],[],[],[],opts_qp);
        if isempty(u) || flag<=0
            u = clip(u_perf, u_min, u_max);
        end
        U(:,k) = u;

        d = dist_seq(:,k);
        x_true = propagate_true_plant(x_true, u, d, T);

        X(:,k+1) = x_true;
        Z(:,k+1) = X(:,k+1) + meas_noise(:,k+1);

        h(k+1) = h_tv(x_true, k*T, obs);

        if norm(X(1:2,k+1)-goal) < 0.5
            X = X(:,1:k+1);
            U = U(:,1:k);
            Z = Z(:,1:k+1);
            h = h(1:k+1);
            break;
        end
    end

    time = 0:T:T*(size(X,2)-1);
    for j = 1:size(X,2)
        h(j) = h_tv(X(:,j), time(j), obs);
    end

    res.X = X;
    res.U = U;
    res.Z = Z;
    res.h = h;
    res.time = time;
end


function [a, c] = fv_affine_at_TV_original(x, t, obs)
% EXACTLY as in the original uploaded baseline script:
% c = c_TI   (TV correction computed but NOT used)

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
    c_TVcorr = Lf_ht + dtLfh + h_tt + 2*h_t; %#ok<NASGU>

    % KEEP ORIGINAL BASELINE EXACTLY:
    c = c_TI;
end


function Mk = compute_margin_sup_TV_original(R_lo, R_hi, xhat, t, obs, umin, umax, gamma)
    [a_hat, c_hat] = fv_affine_at_TV_original(xhat, t, obs);

    support_u = @(Da) max(Da(1)*umin(1), Da(1)*umax(1)) + ...
                      max(Da(2)*umin(2), Da(2)*umax(2));

    M = 256;
    xs = linspace(0,1,round(M^(1/4)));
    [X1,X2,X3,X4] = ndgrid(xs,xs,xs,xs);
    gridPts = [X1(:) X2(:) X3(:) X4(:)].';

    c_now = obs.c_fun(t);

    sup_all = -inf;
    for i = 1:size(gridPts,2)
        r = gridPts(:,i);
        x = R_lo + r.*(R_hi - R_lo);

        [a, c]  = fv_affine_at_TV_original(x, t, obs);
        Dc = c - c_hat;
        Da = a - a_hat;
        sup_u   = support_u(Da);

        v = x(4); th = x(3); dx = x(1)- c_now(1); dy = x(2)- c_now(2);
        c1 = 4*(v*cos(th) + dx);
        c2 = 4*(v*sin(th) + dy);
        sup_d = gamma*(abs(c1) + abs(c2)) + 4*gamma^2;

        sup_all = max(sup_all, Dc + sup_u + sup_d);
    end
    Mk = sup_all;
end


%% ========================================================================
%% Method 2: ORIGINAL Igarashi baseline
%% ========================================================================

function res = run_original_igarashi_baseline(x0, goal, obs, T, N, ...
    u_min, u_max, gamma, epsM, ig, meas_noise, dist_seq)

    X = zeros(4, N+1);
    X(:,1) = x0;

    Z = nan(4, N+1);
    U_total = zeros(2, N);
    U_delta = zeros(2, N);

    h_log    = nan(1, N+1);
    psi1_log = nan(1, N);
    B_log    = nan(1, N);
    I_log    = nan(1, N);
    J_log    = nan(1, N);

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

        d = dist_seq(:,k);
        x_true = propagate_true_plant(x_true, u, d, T);

        X(:,k+1) = x_true;
        Z(:,k+1) = x_true + meas_noise(:,k+1);

        h_log(k+1) = h_tv(X(:,k+1), k*T, obs);

        if norm(X(1:2,k+1)-goal) < 0.5
            X = X(:,1:k+1);
            Z = Z(:,1:k+1);
            U_total = U_total(:,1:k);
            U_delta = U_delta(:,1:k);
            h_log    = h_log(1:k+1);
            psi1_log = psi1_log(1:k);
            B_log    = B_log(1:k);
            I_log    = I_log(1:k);
            J_log    = J_log(1:k);
            break;
        end
    end

    time = 0:T:T*(size(X,2)-1);
    for j = 1:size(X,2)
        h_log(j) = h_tv(X(:,j), time(j), obs);
    end

    res.X = X;
    res.Z = Z;
    res.U = U_total;
    res.U_delta = U_delta;
    res.h = h_log;
    res.psi1 = psi1_log;
    res.B = B_log;
    res.I = I_log;
    res.J = J_log;
    res.time = time;
end


function [u_delta, u_total, info] = igarashi_tv_reciprocal_hocbf( ...
    x, t, u_h, obs, u_min, u_max, alpha_h, K, C, L_B)

    [B, psi1, gradB, Bt, valid] = reciprocal_lifted_barrier_data( ...
        x, t, obs, alpha_h, L_B);

    f = drift_unicycle(x);
    G = input_matrix_unicycle();

    LfB = gradB.' * f;
    LgB = gradB.' * G;

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

    opts = optimoptions('quadprog', 'Display', 'off');
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


function [B, psi1, gradB, Bt, valid] = reciprocal_lifted_barrier_data( ...
    x, t, obs, alpha_h, L_B)

    [psi1, gradPsi1, psi1_t] = psi1_data(x, t, obs, alpha_h);

    eps_psi = 1e-8;
    valid = psi1 > eps_psi;
    psi_safe = max(psi1, eps_psi);

    B = 1/psi_safe + L_B*(x.'*x);

    if valid
        gradB = -(1/psi1^2)*gradPsi1 + 2*L_B*x;
        Bt    = -(1/psi1^2)*psi1_t;
    else
        gradB = zeros(4,1);
        Bt = 0;
    end
end


function [psi1, gradPsi1, psi1_t] = psi1_data(x, t, obs, alpha_h)
    px = x(1);
    py = x(2);
    th = x(3);
    v  = x(4);

    c    = obs.c_fun(t);
    cdot = obs.cdot_fun(t);
    cdd  = obs.cddot_fun(t);

    D    = obs.D_fun(t);
    Dd   = obs.Ddot_fun(t);
    Ddd  = obs.Dddot_fun(t);

    rx = px - c(1);
    ry = py - c(2);

    q = [cos(th); sin(th)];
    qperp = [-sin(th); cos(th)];
    r = [rx; ry];

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

    grad_h = [2*rx;
              2*ry;
              0;
              0];

    gradPsi1 = grad_Lfh + grad_ht + alpha_h*grad_h;

    Lfh_t = -2*v*(cdot.'*q);

    ht_t = 2*(cdot.'*cdot) ...
           - 2*(r.'*cdd) ...
           - 2*(Dd^2) ...
           - 2*D*Ddd;

    psi1_t = Lfh_t + ht_t + alpha_h*h_t;
end


%% ========================================================================
%% Shared helpers
%% ========================================================================

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
    Uopt = fmincon(@cost, U0, [],[],[],[], lb, ub, [], opts);
    u0   = Uopt(1:2);
end


function snap = get_snapshot(res, t_snap, T, obs)
    idx = round(t_snap / T) + 1;
    idx = max(1, min(idx, size(res.X,2)));

    snap.idx = idx;
    snap.t   = res.time(idx);
    snap.x   = res.X(:,idx);
    snap.h   = res.h(idx);
    snap.c   = obs.c_fun(snap.t);
    snap.D   = obs.D_fun(snap.t);
end


function plot_snapshot(res, snap, obs, t_snap, goal, color_rgb, fig_title)
    th = linspace(0,2*pi,250);

    % Path up to snapshot
    plot(res.X(1,1:snap.idx), res.X(2,1:snap.idx), '-', ...
        'Color', color_rgb, 'LineWidth', 2.0);

    % Obstacle at snapshot
    plot(snap.c(1)+snap.D*cos(th), snap.c(2)+snap.D*sin(th), ...
        'k--', 'LineWidth', 1.5);

    % Robot pose at snapshot
    [vx,vy] = triangle_pose(snap.x(1), snap.x(2), snap.x(3), 0.9);
    patch(vx, vy, color_rgb, 'FaceAlpha', 0.4, 'EdgeColor', 'k');

    plot(goal(1), goal(2), 'gx', 'MarkerSize', 10, 'LineWidth', 2.0);

    xlabel('x [m]');
    ylabel('y [m]');
    title(sprintf('%s\n(t = %.1f, h = %.2f)', fig_title, t_snap, snap.h));
end


function [vx,vy] = triangle_pose(x, y, th, s)
    tri = s * [1 0 -0.6;
               0 0.3 -0.3];
    R = [cos(th) -sin(th);
         sin(th)  cos(th)];
    pts = R*tri + [x; y];
    vx = pts(1,:);
    vy = pts(2,:);
end


function y = clip(x, lo, hi)
    y = min(max(x, lo), hi);
end
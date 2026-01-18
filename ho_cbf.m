%% ===== TV-HOCBF (moving obstacle) + Variable Sampling Time (merged) =====
% States: x=[x;y;theta;v], Inputs: u=[u1;u2]=[theta_dot; v_dot]
% Safety: stay outside a time-varying circle: h(x,t) >= 0
%
% Baseline TV version: tv_cbf.txt
% Variable sampling logic: variable_t.txt
%
% NOTE: MPC discretization is kept FIXED at 0.1s (per user request).

clear; clc; close all; rng(0)

% ---------- initial sampling time & horizon ----------
T = 0.1;                         % initial sample time (will vary online)
T_mpc = 0.1;                     % MPC discretization (kept constant)
N = 900;

% ---------- moving obstacle (same as TV code) ----------
obs.c_fun     = @(t) [32 + 25*cos(0.3*t); 25];
obs.cdot_fun  = @(t) [-25*0.3*sin(0.3*t); 0];      % = [-7.5*sin(0.3*t); 0]
obs.cddot_fun = @(t) [-25*(0.3^2)*cos(0.3*t); 0];  % = [-2.25*cos(0.3*t); 0]

obs.D_fun     = @(t) 5 + 2*sin(0.4*t);
obs.Ddot_fun  = @(t) 2*0.4*cos(0.4*t);             % = 0.8*cos(0.4*t)
obs.Dddot_fun = @(t) -2*(0.4^2)*sin(0.4*t);        % = -0.32*sin(0.4*t)

% ---------- robot + controller setup ----------
x_true = [5;25;0;0];       % initial state
goal   = [45;21];

u_min = [-1;-2];  u_max = [1;2]; % box limits

gamma  = 0.3;                    % d_i in [-gamma,gamma]
epsM   = 0.5;                    % measurement error on x,y
k1 = 1; k2 = 1;                  % (kept from your codes)
W = eye(2);                      % objective weight

% ---------- variable sampling parameters (from variable_t) ----------
beta      = 0.6;     % safety factor on (approx) T_max
Tmin      = 0.02;    % lower bound [s]
Tmax_hw   = 0.12;    % upper bound [s]
deltaRT   = 0.005;   % guard above solver time [s]
g_release = 8.0;     % slack threshold to skip bisection
grow_max  = 1.35;    % per-step growth cap

% ---------- logs ----------
X = zeros(4,N+1); X(:,1)=x_true;
U = zeros(2,N);
Z = nan(4, N+1);                % measurements
T_log = nan(1,N);               % T used on step k
Tnext_log = nan(1,N);           % chosen T for next step
thist = nan(1,N+1); thist(1)=0; % time stamps for each X(:,j)

% ---------- TIRA system selector ----------
global system_choice;
system_choice = 99;

% ---------- visualization ----------
t = 0;
Viz = initScene(X(:,1), obs, T);
h_min = +inf;

for k = 1:N
    t_cycle = tic;

    % --- measurement (x,y noisy; theta,v perfect) ---
    z = x_true + [ (2*rand-1)*epsM; (2*rand-1)*epsM; 0; 0 ];
    xhat = z;
    Z(:,k) = z;

    % --- K_perf: MPC with FIXED discretization (requested) ---
    u_perf = Kperf_MPC(xhat, goal, u_min, u_max, T_mpc);

    % --- time-varying affine pieces a(x,t), c(x,t) ---
    [a_hat, c_hat] = fv_affine_at_TV(xhat, t, obs);

    % --- time-Lipschitz bound lFt (same as TV code) ---
    lFt = tv_time_lipschitz(xhat, t, obs);

    % ========== TIRA call for current step length T ==========
    Z_lo = xhat + [-epsM; -epsM; 0; 0];
    Z_hi = xhat + [ epsM;  epsM; 0; 0];

    % p = [u1;u2;d1;d2] ranges (kept as in your codes)
    p_lo = [u_min; -gamma; -gamma];
    p_hi = [u_max; +gamma; +gamma];

    [R_lo, R_hi] = TIRA([0, T], Z_lo, Z_hi, p_lo, p_hi);

    % --- Eq.(12) margin with moving obstacle ---
    Mk = compute_margin_sup_TV(R_lo, R_hi, xhat, t, obs, u_min, u_max, gamma);

    % ---------- QP solve (same structure as TV code) ----------
    % constraint: a_hat*u + c_hat - Mk - lFt*T >= 0
    A = -a_hat;
    b = (c_hat - Mk - lFt*T);

    H = 2*W;
    f = -2*W*u_perf;

    Aqp = [A;  eye(2); -eye(2)];
    bqp = [b;  u_max;   -u_min];

    opts = optimoptions('quadprog','Display','off');
    t_qp = tic;
    [u,~,flag] = quadprog(H,f,Aqp,bqp,[],[],[],[],[],opts);
    solve_time_k = toc(t_qp);

    if isempty(u) || flag<=0
        u = max(min(u_perf,u_max),u_min);
    end
    U(:,k) = u;

    % Slack for release rule (include TV correction term)
    g_now = a_hat*u + c_hat - Mk - lFt*T;

    % ---------- self-triggered update of T (for next step) ----------
    if g_now > g_release
        T_next = min(Tmax_hw, max(Tmin, grow_max*T));
    else
        if feasible_at_TV(Tmin, xhat, u, a_hat, c_hat, lFt, t, obs, u_min, u_max, gamma, epsM)
            T_lo = Tmin;
            T_hi = Tmax_hw;
            for it = 1:8
                T_mid = 0.5*(T_lo + T_hi);
                if feasible_at_TV(T_mid, xhat, u, a_hat, c_hat, lFt, t, obs, u_min, u_max, gamma, epsM)
                    T_lo = T_mid;
                else
                    T_hi = T_mid;
                end
            end
            T_next = beta * T_lo;
        else
            T_next = Tmin;
        end
        T_next = min(max(T_next, Tmin), Tmax_hw);
    end

    % Real-time guard (optional, consistent with variable_t)
    T_next = max(T_next, solve_time_k + deltaRT);
    T_next = min(max(T_next, Tmin), Tmax_hw);

    % ---------- propagate TRUE plant: RK4 + ZOH ----------
    d = gamma*(2*rand(2,1)-1);
    f_dyn = @(x)[ x(4)*cos(x(3)); x(4)*sin(x(3)); 0; 0 ];
    G = [0 0; 0 0; 1 0; 0 1];
    P = [1 0; 0 1; 0 0; 0 0];
    dyn = @(x) f_dyn(x) + G*u + P*d;

    k1r = dyn(x_true);
    k2r = dyn(x_true + 0.5*T*k1r);
    k3r = dyn(x_true + 0.5*T*k2r);
    k4r = dyn(x_true + T*k3r);
    x_true = x_true + (T/6)*(k1r + 2*k2r + 2*k3r + k4r);

    % ---------- log, update time, visualize ----------
    X(:,k+1) = x_true;
    Z(:,k+1) = X(:,k+1) + [ (2*rand-1)*epsM; (2*rand-1)*epsM; 0; 0 ];

    T_log(k) = T;
    Tnext_log(k) = T_next;

    t = t + T;
    thist(k+1) = t;

    c = obs.c_fun(t);
    D = obs.D_fun(t);
    h_curr = (X(1,k+1)-c(1))^2 + (X(2,k+1)-c(2))^2 - D^2;
    h_min  = min(h_min, h_curr);

    comp_ms = 1000*toc(t_cycle);
    updateScene(Viz, X(:,k+1), Z(:,k+1), R_lo, R_hi, obs, t, h_curr, h_min, comp_ms, T_log(k));

    % advance to next step length
    T = T_next;

    if norm(X(1:2,k+1)-goal) < 0.5
        X = X(:,1:k+1);
        U = U(:,1:k);
        Z = Z(:,1:k+1);
        thist = thist(1:k+1);
        T_log = T_log(1:k);
        Tnext_log = Tnext_log(1:k);
        break;
    end
end

% ---- plots (after loop) ----
idx = all(~isnan(Z(1:2,:)), 1);

figure; hold on; axis equal; grid on
th = linspace(0, 2*pi, 200);
co = obs.c_fun(t);
cx = co(1); cy = co(2);
D  = obs.D_fun(t);
plot(cx + D*cos(th), cy + D*sin(th), 'k--', 'LineWidth', 1.5);

plot(Z(1,idx), Z(2,idx), 'k-','LineWidth', 1.2)
plot(goal(1),goal(2),'gx','LineWidth',2)
xlabel x; ylabel y;
xlim([0,60]); ylim([0,50]);
title('TV SD-HOCBF (Eq. 12) with variable sampling')

% Sample period over time (variable sampling visualization)
figure; plot(thist(1:end-1), T_log, 'LineWidth', 1.5); grid on;
xlabel('time t [s]'); ylabel('T_k [s]');
title('Variable sample period');

% Barrier over time: evaluate h(x_j, t_j) at the actual time stamps
hvals = zeros(1, size(X,2));
for j = 1:size(X,2)
    cc = obs.c_fun(thist(j));
    Dj = obs.D_fun(thist(j));
    dx = X(1,j)-cc(1); dy = X(2,j)-cc(2);
    hvals(j) = dx.^2 + dy.^2 - Dj.^2;
end
figure; plot(thist, hvals, 'LineWidth', 1.8); yline(0,'k--'); grid on;
xlabel('time t [s]'); ylabel('h(x,t)');
title('Barrier over time (TV-HOCBF, variable T)');

%% ===== Helper functions (TV version + variable sampling helpers) =====

function [a, c] = fv_affine_at_TV(x, t, obs)
% Time-varying HOCBF (m=2) for unicycle + circular keep-out.
% a(x,t) = L_g L_f h(x,t) (1x2),
% c(x,t) = L_f^2 h + 2 L_f h + 2 h + [ L_f h_t + d/dt(L_f h) + h_tt + 2 h_t ]

    xx=x(1); yy=x(2); th=x(3); v=x(4);

    % obstacle and radius and derivatives
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

    % time derivative terms
    h_t   = -2*(dx*Dxdt + dy*Dydt) - 2*D*Dd;
    Lf_ht = -2*v*(Dxdt*cos(th) + Dydt*sin(th));
    dtLfh = -2*v*(Dxdt*cos(th) + Dydt*sin(th));
    h_tt  = -2*(dx*Dxd2 + dy*Dyd2) + 2*(Dxdt*Dxdt + Dydt*Dydt) ...
            - 2*(Dd*Dd) - 2*D*Ddd;

    c_TI = Lf2h + 2*Lfh + 2*h;
    c_TVcorr = Lf_ht + dtLfh + h_tt + 2*h_t;
    c = c_TI + c_TVcorr;
end

function Mk = compute_margin_sup_TV(R_lo, R_hi, xhat, t, obs, umin, umax, gamma)
% Eq.(12) margin for the TV version (same as your TV code)

    [a_hat, c_hat] = fv_affine_at_TV(xhat, t, obs);

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

        [a, c]  = fv_affine_at_TV(x, t, obs);
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

function Mk_u = compute_margin_sup_TV_given_u(R_lo, R_hi, xhat, t, obs, u, gamma)
% Eq.(12)-style margin but using the ACTUAL held input u (no sup over U).
% Used only for selecting T via feasibility/bisection.

    [a_hat, c_hat] = fv_affine_at_TV(xhat, t, obs);

    xs = linspace(0,1,3);               % 3^4 = 81 samples (faster for bisection)
    [X1,X2,X3,X4] = ndgrid(xs,xs,xs,xs);
    gridPts = [X1(:) X2(:) X3(:) X4(:)].';

    c_now = obs.c_fun(t);

    sup_all = -inf;
    for i = 1:size(gridPts,2)
        r = gridPts(:,i);
        x = R_lo + r.*(R_hi - R_lo);

        [a, c] = fv_affine_at_TV(x, t, obs);
        Dc = c - c_hat;
        Du = (a - a_hat) * u;

        v = x(4); th = x(3); dx = x(1)- c_now(1); dy = x(2)- c_now(2);
        c1 = 4*(v*cos(th) + dx);
        c2 = 4*(v*sin(th) + dy);
        sup_d = gamma*(abs(c1) + abs(c2)) + 4*gamma^2;

        sup_all = max(sup_all, Dc + Du + sup_d);
    end
    Mk_u = sup_all;
end

function ok = feasible_at_TV(Tcand, xhat, u, a_hat, c_hat, lFt, t, obs, ...
                             u_min, u_max, gamma, epsM)
% Check if holding the CURRENT u is safe up to Tcand (TV version).
% Uses a tube from TIRA and an actual-u margin, and includes the -lFt*Tcand term.

    Z_lo = xhat + [-epsM; -epsM; 0; 0];
    Z_hi = xhat + [ epsM;  epsM; 0; 0];

    p_lo = [u_min; -gamma; -gamma];
    p_hi = [u_max; +gamma; +gamma];

    [R_lo, R_hi] = TIRA([0, Tcand], Z_lo, Z_hi, p_lo, p_hi);

    Mk_T = compute_margin_sup_TV_given_u(R_lo, R_hi, xhat, t, obs, u, gamma);

    g_T = a_hat*u + c_hat - Mk_T - lFt*Tcand;
    ok = (g_T >= 0);
end

function u0 = Kperf_MPC(x0, goal, umin, umax, T)
% Same MPC structure as your TV code, but T is supplied by caller.
% In this merged script, the caller always uses T_mpc = 0.1.

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
        for tt = 1:Np
            u = U(2*tt-1:2*tt);
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

function xnext = rk4_unicycle(x, u, T)
    f  = @(x)[ x(4)*cos(x(3)); x(4)*sin(x(3)); 0; 0 ];
    G  = [0 0; 0 0; 1 0; 0 1];
    dyn = @(x) f(x) + G*u;
    k1 = dyn(x);
    k2 = dyn(x + 0.5*T*k1);
    k3 = dyn(x + 0.5*T*k2);
    k4 = dyn(x + T*k3);
    xnext = x + (T/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function lFt = tv_time_lipschitz(xhat, t, obs)
% Conservative time-Lipschitz bound on c(xhat,t)
    dt = 1e-3;
    [~, c0] = fv_affine_at_TV(xhat, t,     obs);
    [~, c1] = fv_affine_at_TV(xhat, t+dt,  obs);
    lFt = abs(c1 - c0)/dt;
end

function H = initScene(x0, obs, T)
    figure('Color','w'); hold on; axis equal; grid on
    xlim([0, 60]); ylim([0, 50]);
    xlabel('x [m]'); ylabel('y [m]');
    title('Unicycle - robust SD-HOCBF (Eq. 12) + variable T')

    th = linspace(0,2*pi,200);
    H.th = th;
    co = obs.c_fun(0); D = obs.D_fun(0);
    H.obs = plot(co(1)+D*cos(th), co(2)+D*sin(th),'k--','LineWidth',1.5);

    H.path_true = animatedline('LineWidth',1.6,'Color',[0 0 0]);
    H.path_meas = animatedline('LineWidth',1.0,'Color',[0.6 0.6 0.6]);

    [vx,vy] = triangle_pose(x0(1), x0(2), x0(3), 0.7);
    H.rob = patch('XData',vx,'YData',vy,'FaceColor',[0.2 0.6 1],...
                  'EdgeColor','k','LineWidth',1);

    H.tira = patch('XData',[],'YData',[],'FaceColor',[1 0.6 0],...
                   'FaceAlpha',0.2,'EdgeColor',[1 0.4 0],'LineWidth',1);

    H.txt = text(1,49,sprintf('t=%.1f  T=%.3f  h=%.2f  min h=%.2f  comp=%.0f ms',0,T,NaN,NaN,0),...
                 'FontName','Consolas','FontSize',10,'VerticalAlignment','top');

    drawnow;
end

function updateScene(H, x_true, z_meas, R_lo, R_hi, obs, t, h_curr, h_min, comp_ms, T)
    c = obs.c_fun(t); D = obs.D_fun(t);
    set(H.obs,'XData',c(1)+D*cos(H.th),'YData',c(2)+D*sin(H.th));

    addpoints(H.path_true, x_true(1), x_true(2));
    if ~isempty(z_meas)
        addpoints(H.path_meas, z_meas(1), z_meas(2));
    end

    [vx,vy] = triangle_pose(x_true(1), x_true(2), x_true(3), 0.7);
    set(H.rob,'XData',vx,'YData',vy);

    if ~isempty(R_lo) && ~isempty(R_hi)
        Xr = [R_lo(1) R_hi(1) R_hi(1) R_lo(1)];
        Yr = [R_lo(2) R_lo(2) R_hi(2) R_hi(2)];
        set(H.tira,'XData',Xr,'YData',Yr);
    end

    set(H.txt,'String',sprintf('t=%.1f  T=%.3f  h=%.2f  min h=%.2f  comp=%d ms', ...
        t, T, h_curr, h_min, round(comp_ms)));

    drawnow limitrate
end

function [vx,vy] = triangle_pose(x, y, th, s)
    tri = s * [ 1  0 -0.6; 0  0.3 -0.3 ];
    R = [cos(th) -sin(th); sin(th) cos(th)];
    pts = R*tri + [x; y];
    vx = pts(1,:); vy = pts(2,:);
end

%% ===== HOCBF (relative-degree-2) — Unicycle circle avoidance via QP =====
% States: x=[x;y;theta;v], Inputs: u=[u1;u2]=[theta_dot; v_dot]
% Safety: stay outside a circle: h(x) >= 0
clear; clc; close all; rng(0)

% ---------- fixed scenario ----------
T = 0.1;                         % sample time
N = 900;
obs.c = [32;25];
D = 5;       % obstacle center & radius (D)
x_true = [5;25;-pi/2;0.8];       % initial state
goal   = [45;21];

u_min = [-1;-2];  u_max = [1;2]; % box limits (Eq. 11)
gamma  = 0.3;                    % d_i in [-0.3,0.3]
epsM   = 0.5;                    % measurement error on x,y
k1 = 1; k2 = 1;                  % (same as paper)
W = eye(2);                      % objective weight
s_bar  = 2.45;           % bound on ||f(x)+g(x)u|| (pick a safe value)

params = struct();
params.epsM = epsM;
params.s_bar = s_bar;
beta     = 0.6;     % safety factor on T_max (0<beta<=1)
Tmin     = 0.02;    % lower bound [s]  (>= solver time + a little margin)
Tmax_hw  = 0.12;    % upper bound [s]  (hardware/OS cap)
deltaRT  = 0.005;   % guard over solver time [s]
g_release = 8.0;        % margin threshold to skip bisection
grow_max  = 1.35;       % cap per-step growth (35%)

X = zeros(4,N+1); X(:,1)=x_true; 
U = zeros(2,N);
Z = nan(4, N+1);   % to store measurements
T_log     = zeros(1,N);     % T_k each step
dT_log    = zeros(1,N);     % absolute change ΔT_k = T_k - T_{k-1}
dTpct_log = zeros(1,N);     % percent change 100*ΔT_k/T_{k-1}
t=0;
global system_choice;
system_choice = 99;

for k = 1:N
    % --- measurement (x,y noisy; theta,v perfect) ---
    z = x_true + [ (2*rand-1)*epsM; (2*rand-1)*epsM; 0; 0 ];
    xhat = z;
    Z(:,k) = z;
    % --- K_perf: minimal MPC (receding horizon, see function below) ---
    u_perf = Kperf_MPC(xhat, goal, u_min, u_max, T);

    % --- Eq. (12) pieces at xhat: a_hat*u + c_hat  ---
    [a_hat, c_hat] = fv_affine_at(xhat, obs.c, D, k1, k2);

    % ========== TIRA call: one-step tube R(xhat,T) ==========
    
    % measurement set Z around xhat
    Z_lo = xhat + [-epsM; -epsM; 0; 0];
    Z_hi = xhat + [ epsM;  epsM; 0; 0];

    % inputs p = [u1;u2;d1;d2] ranges (constant over [0,T])
    p_lo = [u_min; -gamma; -gamma];
    p_hi = [u_max; +gamma; +gamma];


    % reachability of the continuous-time system over [0,T]
    [R_lo, R_hi] = TIRA([0, T], Z_lo, Z_hi, p_lo, p_hi);   % uses your System_description case=99

    % --- compute the margin M_k = sup_{x∈R,u∈U,d∈P} margin ---
    % Build margin pieces for this step
    Mk_sup = compute_margin_sup(R_lo, R_hi, xhat, obs.c, D, u_min, u_max, gamma);
    
    Mk = Mk_sup;                 % Eq.(12) margin from TIRA for current T
    A  = -a_hat; 
    b  = (c_hat - Mk);
    H = 2*W;       f = -2*W*u_perf;
    Aqp = [A;  eye(2); -eye(2)];
    bqp = [b;  u_max;   -u_min];
    opts = optimoptions('quadprog','Display','off');
    
    tStart = tic;
    [u,~,flag] = quadprog(H,f,Aqp,bqp,[],[],[],[],[],opts);
    solve_time_k = toc(tStart);
    if isempty(u) || flag<=0, u = max(min(u_perf,u_max),u_min); end
    U(:,k) = u;
    g_now = a_hat*u + c_hat - Mk;

    % ---------- self-triggered update of T_k (for next step) ----------
    % We look for the largest T in [Tmin, Tmax_hw] such that
    %   a_hat*u + c_hat - Mk(T) - lFt*T >= 0
    % using the current u, xhat, and time t.
    if g_now > g_release
        T_next = min(Tmax_hw, max(Tmin, min(grow_max*T, Tmax_hw)));
        T_log(k)   = T_next;
    else
    
        % 1) check feasibility at Tmin
        if feasible_at_TI(Tmin, xhat, u, a_hat, c_hat, obs.c, D, ...
                          u_min, u_max, gamma, epsM)
            % 2) bisection on [Tmin, Tmax_hw]
            T_lo = Tmin;
            T_hi = Tmax_hw;
            for it = 1:8   % 8 iterations is usually plenty
                T_mid = 0.5*(T_lo + T_hi);
                if feasible_at_TI(T_mid, xhat, u, a_hat, c_hat, obs.c, D, ...
                                  u_min, u_max, gamma, epsM)
                    % still safe at T_mid -> try longer hold
                    T_lo = T_mid;
                else
                    % unsafe at T_mid -> shorten
                    T_hi = T_mid;
                end
            end
            T_next = beta * T_lo;   % safety factor
        else
            % cannot even guarantee Tmin with this u -> clamp to Tmin
            T_next = Tmin;
        end
    
        % clip and log
        T_next     = min(max(T_next, Tmin), Tmax_hw);
        T_log(k)   = T_next;
    end

    % ---------- propagate TRUE plant: RK4 + ZOH (Eq. 11) ----------
    d = gamma*(2*rand(2,1)-1);  % adversarial within bounds
    f = @(x)[ x(4)*cos(x(3)); x(4)*sin(x(3)); 0; 0 ];
    G = [0 0; 0 0; 1 0; 0 1];
    P = [1 0; 0 1; 0 0; 0 0];
    dyn = @(x) f(x) + G*u + P*d;

    k1r = dyn(x_true);
    k2r = dyn(x_true + 0.5*T*k1r);
    k3r = dyn(x_true + 0.5*T*k2r);
    k4r = dyn(x_true + T*k3r);
    x_true = x_true + (T/6)*(k1r + 2*k2r + 2*k3r + k4r);

    X(:,k+1) = x_true;
    Z(:,k+1) = X(:,k+1) + [ (2*rand-1)*epsM; (2*rand-1)*epsM; 0; 0 ];
    t = t + T;
    T = T_next;

    if norm(X(1:2,k+1)-goal) < 1 , X=X(:,1:k+1); U=U(:,1:k); break; end
end

% ---- plots (after loop) ----
idx = all(~isnan(Z(1:2,:)), 1);    % keep samples where x,y measurements exist

figure; hold on; axis equal; grid on
th = linspace(0,2*pi,256);
plot(obs.c(1)+D*cos(th), obs.c(2)+D*sin(th),'k--','LineWidth',1.5)
%plot(X(1,:),X(2,:),'-','LineWidth',2)

plot(Z(1,idx), Z(2,idx), 'k-','LineWidth', 1.2)        % line

plot(goal(1),goal(2),'gx','LineWidth',2)
xlabel x; ylabel y;
xlim([0,50]);
ylim([0,50]);

title('Eq.(12) with TIRA margin — unicycle')

figure; plot(1:N, T_log, 'LineWidth',1.5); grid on;
xlabel('step k'); ylabel('T_k [s]'); title('Sample period');
xlim([1 100]);                 % show only first 100 on x-axi

k_col     = (1:N).';
T_col     = T_log(:);                        % Nx1
dT_col    = [0; diff(T_col)];                % Nx1
dTpct_col = [0; 100 * (diff(T_col) ./ max(T_col(1:end-1), eps))];  % Nx1, guarded /0

log_mat = [k_col, T_col, dT_col, dTpct_col]; % N x 4
writematrix(log_mat, 'T_log.csv');
xlim([1 100]);                 % show only first 100 on x-axis

% Barrier check
dx = X(1,:)-obs.c(1); dy = X(2,:)-obs.c(2);
hvals = dx.^2 + dy.^2 - D^2;
figure; plot(hvals,'LineWidth',1.8); yline(0,'k--');
xlabel step; ylabel('h(x)'); title('Barrier over time (should stay ≥ 0)')

%% ===== Helper: projection onto halfspace ∩ box (no toolbox) =====
function [a, c] = fv_affine_at(x, cObs, D, k1, k2)
% a(x) = L_g L_f h(x) (1x2 row),  c(x) = L_f^2 h + 2 L_f h + 2 h  (Eq. 12)
    xx=x(1); yy=x(2); th=x(3); v=x(4);
    dx=xx-cObs(1); dy=yy-cObs(2);
    phi = dx*cos(th) + dy*sin(th);
    h   = dx*dx + dy*dy - D^2;
    Lfh = 2*v*phi;
    Lf2h= 2*v^2;
    a1 = 2*v*(-dx*sin(th) + dy*cos(th));
    a2 = 2*phi;
    a  = [a1, a2];
    c  = Lf2h + 2*Lfh + 2*h;
end
function Mk_sup= compute_margin_sup(R_lo, R_hi, xhat, cObs, D, umin, umax, gamma)
    % === your existing sup over the tube (keep this) ===
    [a_hat, c_hat] = fv_affine_at(xhat, cObs, D, 1, 1); %#ok<NASGU>
    sup_all = -inf;
    M = 256;
    xs = linspace(0,1,round(M^(1/4)));
    [X1,X2,X3,X4] = ndgrid(xs,xs,xs,xs);
    gridPts = [X1(:) X2(:) X3(:) X4(:)].';
    for i = 1:size(gridPts,2)
        r = gridPts(:,i);
        x = R_lo + r.*(R_hi - R_lo);
        [a, c] = fv_affine_at(x, cObs, D, 1, 1);
        Dc = c - c_hat;
        Da = a - a_hat;
        sup_u = max(Da(1)*umin(1), Da(1)*umax(1)) + ...
                max(Da(2)*umin(2), Da(2)*umax(2));
        v = x(4); th = x(3); dx = x(1)-cObs(1); dy = x(2)-cObs(2);
        c1 = 4*(v*cos(th) + dx);
        c2 = 4*(v*sin(th) + dy);
        sup_d = gamma*(abs(c1) + abs(c2)) + 4*gamma^2;
        sup_all = max(sup_all, Dc + sup_u + sup_d);
    end
    Mk_sup = sup_all;  % <- constant tube/disturbance margin (keep)

end

    

function u0 = Kperf_MPC(x0, goal, umin, umax, T)
    Np = 70;                                   % horizon
    Qp = 5.0;                                  % stage position weight
    Qf = 400;                                   % terminal position weight
    Rw = diag([0.1, 0.05]);                    % input weight

    U0 = zeros(2*Np,1);
    lb = repmat([umin(1); umin(2)], Np, 1);
    ub = repmat([umax(1); umax(2)], Np, 1);

    function J = cost(U)
        x = x0;
        J = 0;
        for t = 1:Np
            u = U(2*t-1:2*t);
            x = rk4_unicycle(x, u, T);         % roll dynamics one step
            pos_err = x(1:2) - goal(:);
            J = J + Qp*(pos_err.'*pos_err) + u.'*Rw*u;   % stage cost
        end
        % ---------- TERMINAL COST (add here) ----------
        term_err = x(1:2) - goal(:);           % x is x_{k+Np}
        J = J + Qf*(term_err.'*term_err);      % terminal position penalty
        % ---------------------------------------------
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

function ok = feasible_at_TI(Tcand, xhat, u, a_hat, c_hat, cObs, D, ...
                             u_min, u_max, gamma, epsM)
%FEASIBLE_AT_TI  Check if fixed u is safe for a candidate Tcand (time-invariant).
% Implements the same inequality as in the QP:
%       a_hat*u + c_hat - Mk(Tcand) >= 0,
% where Mk(Tcand) is the Eq.(12) margin computed with TIRA and compute_margin_sup.

    % measurement set Z around xhat
    Z_lo = xhat + [-epsM; -epsM; 0; 0];
    Z_hi = xhat + [ epsM;  epsM; 0; 0];

    % input/disturbance bounds p = [u1;u2;d1;d2]
    p_lo = [u_min; -gamma; -gamma];
    p_hi = [u_max; +gamma; +gamma];

    % reachable tube over [0, Tcand]
    [R_lo, R_hi] = TIRA([0, Tcand], Z_lo, Z_hi, p_lo, p_hi);

    % margin for this Tcand (time-invariant version)
    Mk_T = compute_margin_sup_given_u(R_lo, R_hi, xhat, cObs, D, u, gamma);
    % same inequality as in the QP
    g_T = a_hat*u + c_hat - Mk_T;

    ok = (g_T >= 0);
end
function Mk_u = compute_margin_sup_given_u(R_lo, R_hi, xhat, cObs, D, u, gamma)
% Sup over the tube using the ACTUAL held input u (no sup over the full box).
    [a_hat, c_hat] = fv_affine_at(xhat, cObs, D, 1, 1);
    sup_all = -inf;

    % Coarse but fast grid; increase to 3^4 or 4^4 if needed
    xs = linspace(0,1,3);                         % 3 points/axis → 81 samples
    [X1,X2,X3,X4] = ndgrid(xs,xs,xs,xs);
    gridPts = [X1(:) X2(:) X3(:) X4(:)].';

    for i = 1:size(gridPts,2)
        r  = gridPts(:,i);
        x  = R_lo + r.*(R_hi - R_lo);

        [a, c] = fv_affine_at(x, cObs, D, 1, 1);
        Dc     = c - c_hat;
        Du     = (a - a_hat) * u;                 % <-- actual input, not sup over U

        % disturbance support (same as before)
        v = x(4); th = x(3); dx = x(1)-cObs(1); dy = x(2)-cObs(2);
        c1 = 4*(v*cos(th) + dx);
        c2 = 4*(v*sin(th) + dy);
        sup_d = gamma*(abs(c1) + abs(c2)) + 4*gamma^2;

        sup_all = max(sup_all, Dc + Du + sup_d);
    end
    Mk_u = sup_all;
end

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
beta     = 0.005;     % safety factor on T_max (0<beta<=1)
Tmin     = 0.02;    % lower bound [s]  (>= solver time + a little margin)
Tmax_hw  = 0.12;    % upper bound [s]  (hardware/OS cap)
deltaRT  = 0.005;   % guard over solver time [s]

X = zeros(4,N+1); X(:,1)=x_true; 
U = zeros(2,N);
Z = nan(4, N+1);   % to store measurements
T_log     = zeros(1,N);     % T_k each step
dT_log    = zeros(1,N);     % absolute change ΔT_k = T_k - T_{k-1}
dTpct_log = zeros(1,N);     % percent change 100*ΔT_k/T_{k-1}

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
    [Mk, coeff] = compute_margin_sup(R_lo, R_hi, xhat, obs.c, D, u_min, u_max, gamma, params);
    
    %QP solver
    A = -a_hat;    b = (c_hat - Mk);
    H = 2*W;       f = -2*W*u_perf;
    Aqp = [A;  eye(2); -eye(2)];
    bqp = [b;  u_max;   -u_min];
    opts = optimoptions('quadprog','Display','off');
    
    tStart = tic;
    [u,~,flag] = quadprog(H,f,Aqp,bqp,[],[],[],[],[],opts);
    solve_time_k = toc(tStart);
    
    if isempty(u) || flag<=0, u = max(min(u_perf,u_max),u_min); end
    U(:,k) = u;
    % A0 = a_hat*u + c_hat - ||L_p psi_{m-1}||*gamma   (your existing pieces)
    A0 = a_hat*u + c_hat - coeff.LpPsi_m1*gamma;
    
    % Prefer tighter L(u_k) if you exposed Lipschitz parts, else fall back to Lu_max
    if isfield(coeff,'l_Lfpsi')
        L_u = coeff.l_Lfpsi + coeff.l_alpha_psi + coeff.l_Lgpsi * norm(u);
    else
        L_u = coeff.Lu_max;  % conservative but OK
    end
    
    eps_k  = coeff.eps;
    Theta  = coeff.Theta;
    
    num = A0 - L_u*eps_k;
    den = L_u * max(Theta, 1e-12);
    if num <= 0
        T = Tmin;
    else
        T = beta * (num/den);
        T = min(max(T, Tmin), Tmax_hw);
    end
    T = max(T, solve_time_k + deltaRT);   % causality guard
    if k == 1
        T_log(k)     = T;
        dT_log(k)    = 0;
        dTpct_log(k) = 0;
    else
        dT           = T - T_log(k-1);
        T_log(k)     = T;
        dT_log(k)    = dT;
        dTpct_log(k) = 100 * dT / T_log(k-1);
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


    if norm(X(1:2,k+1)-goal) < 0.5, X=X(:,1:k+1); U=U(:,1:k); break; end
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

figure; stem(2:N, dT_log(2:N), 'filled'); grid on;
xlabel('step k'); ylabel('\Delta T_k [s]'); title('Change in T');


k_col     = (1:N).';
T_col     = T_log(:);                        % Nx1
dT_col    = [0; diff(T_col)];                % Nx1
dTpct_col = [0; 100 * (diff(T_col) ./ max(T_col(1:end-1), eps))];  % Nx1, guarded /0

log_mat = [k_col, T_col, dT_col, dTpct_col]; % N x 4
writematrix(log_mat, 'T_log.csv');

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

%function Mk = compute_margin_sup(R_lo, R_hi, xhat, cObs, D, k1, k2, umin, umax, gamma)
% M_k = sup_{x∈R,u∈U,d∈P} [ fv(x,u)-fv(xhat,u) + f_d(x,d) ]
%     = sup_{x∈R} [ Δc(x) + sup_u Δa(x)u ]  + sup_{x∈R,d∈P} f_d(x,d)
    % precompute at xhat
 %   [a_hat, c_hat] = fv_affine_at(xhat, cObs, D, k1, k2);

    % sampling on the TIRA box (increase M for tighter bound)
  %  M = 400;
   % sup1 = -inf;  sup2 = -inf;

    %for i = 1:M
     %   r = rand(4,1);
      %  x = R_lo + r.*(R_hi - R_lo);

        % Δc and Δa
       % [a, c] = fv_affine_at(x, cObs, D, k1, k2);
        %Dc = c - c_hat;
        %Da = a - a_hat;

        % sup over u in the input box (componentwise)
        %uabs = [max(abs(umin(1)),abs(umax(1))), max(abs(umin(2)),abs(umax(2)))];
        %sup_u = sum(abs(Da).*uabs);

        %sup1 = max(sup1, Dc + sup_u);

        % sup over d in [-gamma,gamma]^2 for f_d(x,d)
        % f_d(x,d) = 4*v(d1 cosθ + d2 sinθ) + 2(d1^2 + d2^2) + 4(dx d1 + dy d2)
        %v=x(4); th=x(3); dx=x(1)-cObs(1); dy=x(2)-cObs(2);
        %c1 = 4*(v*cos(th) + dx);
        %c2 = 4*(v*sin(th) + dy);
        %sup_d = gamma*(abs(c1) + abs(c2)) + 4*gamma^2;

        %sup2 = max(sup2, sup_d);
    %end
    %Mk = sup1; %+ sup2;
%end
function [Mk, coeff] = compute_margin_sup(R_lo, R_hi, xhat, cObs, D, umin, umax, gamma, params)
% Paper-faithful Eq. (12):
% margin = sup_{x∈R, u∈U, d∈P} [ fv(x,u) - fv(xhat,u) + f_d(x,d) ]
% We over-approximate sup_{x∈R} by sampling the TIRA box [R_lo, R_hi]
% (implementation detail; the formula itself is exactly Eq. (12)).

    % Precompute fv pieces at xhat
    [a_hat, c_hat] = fv_affine_at(xhat, cObs, D, 1, 1); %#ok<NASGU> k1,k2 unused here

    % Support function of the input box U for the affine term Δa(x)·u
    uabs = [max(abs(umin(1)),abs(umax(1))), max(abs(umin(2)),abs(umax(2)))];

    % Deterministic grid samples over the tube (covers corners)
    M = 256;
    xs = linspace(0,1,round(M^(1/4)));
    [X1,X2,X3,X4] = ndgrid(xs,xs,xs,xs);
    gridPts = [X1(:) X2(:) X3(:) X4(:)].';

    sup_all = -inf;
    for i = 1:size(gridPts,2)
        r = gridPts(:,i);
        x = R_lo + r.*(R_hi - R_lo);

        % Δfv(x,u) = [Δc(x) + Δa(x)·u] with u ∈ U (box)
        [a, c] = fv_affine_at(x, cObs, D, 1, 1);
        Dc = c - c_hat;
        Da = a - a_hat;
        sup_u = max(Da(1)*umin(1), Da(1)*umax(1)) + ...
        max(Da(2)*umin(2), Da(2)*umax(2));

        % f_d(x,d) from Eq. (12), with exact box sup over d ∈ [-γ,γ]^2
        v = x(4); th = x(3); dx = x(1)-cObs(1); dy = x(2)-cObs(2);
        c1 = 4*(v*cos(th) + dx);
        c2 = 4*(v*sin(th) + dy);
        sup_d = gamma*(abs(c1) + abs(c2)) + 4*gamma^2;

        sup_all = max(sup_all, Dc + sup_u + sup_d);
    end
    Mk = sup_all;
    if nargin>=9 && isfield(params,'epsM')
        eps_k = sqrt(2) * params.epsM;   % from ±epsM on x,y ⇒ Euclidean radius
    else
        % fallback if you already keep a sensor box elsewhere; edit if needed
        eps_k = sqrt(2) * 1e-2;          % <--- replace with your value
    end
    
    % (b) Theta = s_bar + ||L_p h||_sup * gamma   over the reachable x,y box
    % h(x,y) = (x-cx)^2 + (y-cy)^2 - D^2  ⇒  ∇h = [2(x-cx), 2(y-cy), 0, 0]
    % ||L_p h|| uses only x,y (your P projects to x,y), and the max over a box
    % occurs at the farthest corner from the obstacle center c.
    corners = [R_lo(1) R_lo(2);
               R_lo(1) R_hi(2);
               R_hi(1) R_lo(2);
               R_hi(1) R_hi(2)];
    dxy     = corners - cObs(1:2).';
    rmax    = max( sqrt(sum(dxy.^2,2)) );    % farthest (x,y) from c
    Lp_h_sup = 2*rmax;                        % ||∇h|| at the farthest corner
    
    if nargin>=9 && isfield(params,'s_bar')
        s_bar = params.s_bar;
    else
        % fallback (conservative): edit to your plant; for unicycle, s_bar ≈ vmax + ||u||_max
        umax_norm = norm([max(abs(u_min(1)),abs(u_max(1))), ...
                          max(abs(u_min(2)),abs(u_max(2)))],2);
        s_bar = umax_norm;   % minimal placeholder; replace with your bound if you have one
    end
    
    Theta = s_bar + Lp_h_sup * gamma;
    
    % (c) L(u) term used in M_k(T) — the QP needs an affine constraint, so we
    % keep the **u-bounded** version (with ||u|| ≤ ||u||_max). This is what you
    % already use inside Mk; we just expose it as Lu_max.
    if nargin>=9 && all(isfield(params, {'l_Lfpsi','l_alpha_psi','l_Lgpsi'}))
        umax_norm = norm([max(abs(u_min(1)),abs(u_max(1))), ...
                          max(abs(u_min(2)),abs(u_max(2)))],2);
        Lu_max = params.l_Lfpsi + params.l_alpha_psi + params.l_Lgpsi * umax_norm;
    else
        % fallback: if you already computed a scalar L inside Mk, reuse it here
        % (rename your internal variable to Lu_max and assign it below).
        % Otherwise set a conservative placeholder and tighten later.
        Lu_max = 1.0;
    end
    
    % (d) Pack outputs for use in the T_k rule outside this function
    coeff = struct();
    coeff.eps    = eps_k;     % ε(z_k)
    coeff.Theta  = Theta;     % s_bar + ||L_p h||_sup * γ
    coeff.Lu_max = Lu_max;    % L with ||u|| replaced by ||u||_max
    
    % Optional: expose the individual Lipschitz pieces if you want to compute L(u_k) later
    if nargin>=9 && all(isfield(params, {'l_Lfpsi','l_alpha_psi','l_Lgpsi'}))
        coeff.l_Lfpsi     = params.l_Lfpsi;
        coeff.l_alpha_psi = params.l_alpha_psi;
        coeff.l_Lgpsi     = params.l_Lgpsi;
    end
    % ||L_p psi_1|| at xhat (rel.-deg 2 case)
    dx = xhat(1) - cObs(1);  dy = xhat(2) - cObs(2);
    th = xhat(3);         v  = xhat(4);
    gx = 2*v*cos(th) + 2*dx;
    gy = 2*v*sin(th) + 2*dy;
    coeff.LpPsi_m1 = hypot(gx, gy);
    % NEW: expose the exact constant and slope used by Mk(T)
    coeff.Mk_const = coeff.Lu_max * coeff.eps;         % = L_u * eps_k
    coeff.Mk_slope = coeff.Lu_max * coeff.Theta;       % = L_u * Theta


end

function u0 = Kperf_MPC(x0, goal, umin, umax, T)
    Np = 50;                                   % horizon
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


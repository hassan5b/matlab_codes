%% ===== HOCBF (relative-degree-2) — Unicycle circle avoidance via QP =====
% States: x=[x;y;theta;v], Inputs: u=[u1;u2]=[theta_dot; v_dot]
% Safety: stay outside a circle: h(x) >= 0
clear; clc; close all; rng(0)

% ---------- fixed scenario ----------
T = 0.1;                         % sample time
N = 200;
% Scale (goal x = 45 -> 8 m)
s = 8/45;

% Original (grid)
x0y0_grid = [5; 25];

obs_c_grid = [32; 25];
goal_grid  = [45; 21];

% Translate in grid so start is (0,0)
obs.c = s * (obs_c_grid - x0y0_grid);   % = s*[27; 0]  = [4.8; 0] m
D     = s * 5;                           % = 0.8888889 m

x_true = [0; 0; 0; 0];                   % start at origin after transform
goal   = s * (goal_grid - x0y0_grid);    % = s*[40; -4] = [7.1111111; -0.7111111] m

u_min = [-1.7825; -3.5];
u_max = [ 1.7825;  3.5];

gamma = s * 0;                         % = 0.0533333 m (only if gamma is x/y position disturbance)
epsM  = s * 0;                         % = 0.0888889 m (x/y measurement error)

k1 = 1; k2 = 1;                  % (same as paper)
W = eye(2);                      % objective weight

X = zeros(4,N+1); X(:,1)=x_true; 
U = zeros(2,N);
Z = nan(4, N+1);   % to store measurements
flag_log  = zeros(N,1);
h_log     = zeros(N,1);
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
    Mk = compute_margin_sup(R_lo, R_hi, xhat, obs.c, D, u_min, u_max, gamma);

    % ---------- Solve Eq. (12): QP ----------
    % constraint: a_hat*u + c_hat - M_k >= 0  ->  [-a_hat] u <= -(c_hat - M_k)
    A = -a_hat;    b = (c_hat - Mk);
    H = 2*W;       f = -2*W*u_perf;
    Aqp = [A;  eye(2); -eye(2)];
    bqp = [b;  u_max;   -u_min];
    opts = optimoptions('quadprog','Display','off');
    [u,~,flag] = quadprog(H,f,Aqp,bqp,[],[],[],[],[],opts);
    if isempty(u) || flag<=0, u = max(min(u_perf,u_max),u_min); end
    U(:,k) = u;
    flag_log(k) = flag;


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


% Barrier check
dx = X(1,:)-obs.c(1); dy = X(2,:)-obs.c(2);
hvals = dx.^2 + dy.^2 - D^2;
figure; plot(hvals,'LineWidth',1.8); yline(0,'k--');
xlabel step; ylabel('h(x)'); title('Barrier over time (should stay ≥ 0)')

figure; plot(flag_log,'LineWidth',1.2); grid on;
xlabel('k'); ylabel('quadprog exitflag'); title('QP exitflag');
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
function Mk = compute_margin_sup(R_lo, R_hi, xhat, cObs, D, umin, umax, gamma)
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




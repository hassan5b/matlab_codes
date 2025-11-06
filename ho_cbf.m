%% ===== HOCBF (relative-degree-2) — Unicycle circle avoidance via QP =====
% States: x=[x;y;theta;v], Inputs: u=[u1;u2]=[theta_dot; v_dot]
% Safety: stay outside a circle: h(x) >= 0
clear; clc; close all; rng(0)

% ---------- fixed scenario ----------
T = 0.1;                         % sample time
N = 900;
%obs.c_fun     = @(tt) [32;25];         % o(t)
%obs.cdot_fun  = @(tt) [0;0];           % \dot o(t)
%obs.cddot_fun = @(tt) [0;0];           % \ddot o(t)
%obs.D_fun     = @(tt) 5.0;             % D(t)
%obs.Ddot_fun  = @(tt) 0.0;             % \dot D(t)
%obs.Dddot_fun = @(tt) 0.0;             % \ddot D(t)

%moving obstacle
obs.c_fun     = @(t) [32 + 25*cos(0.3*t); 25];
obs.cdot_fun  = @(t) [-25*0.3*sin(0.3*t); 0];      % = [-7.5*sin(0.3*t); 0]
obs.cddot_fun = @(t) [-25*(0.3^2)*cos(0.3*t); 0];  % = [-2.25*cos(0.3*t); 0]

obs.D_fun     = @(t) 5 + 2*sin(0.4*t);
obs.Ddot_fun  = @(t) 2*0.4*cos(0.4*t);             % = 0.8*cos(0.4*t)
obs.Dddot_fun = @(t) -2*(0.4^2)*sin(0.4*t);        % = -0.32*sin(0.4*t)


x_true = [5;25;-pi/2;0.8];       % initial state
goal   = [45;21];

u_min = [-1;-2];  u_max = [1;2]; % box limits (Eq. 11)
gamma  = 0.3;                    % d_i in [-0.3,0.3]
epsM   = 0.5;                    % measurement error on x,y
k1 = 1; k2 = 1;                  % (same as paper)
W = eye(2);                      % objective weight

X = zeros(4,N+1); X(:,1)=x_true; 
U = zeros(2,N);
Z = nan(4, N+1);   % to store measurements

global system_choice;
system_choice = 99;
t=0;
Viz = initScene(X(:,1), obs, T);   % from the helper I gave
h_min = +inf;

for k = 1:N
    % --- measurement (x,y noisy; theta,v perfect) ---
    z = x_true + [ (2*rand-1)*epsM; (2*rand-1)*epsM; 0; 0 ];
    xhat = z;
    Z(:,k) = z;

    c_now = obs.c_fun(t);
    D_now = obs.D_fun(t);
    
    % --- K_perf (unchanged) ---
    u_perf = Kperf_MPC(xhat, goal, u_min, u_max, T);
    
    % --- TIME-VARYING affine pieces a(x,t), c(x,t) ---
    [a_hat, c_hat] = fv_affine_at_TV(xhat, t, obs);
    
    % --- time-Lipschitz bound (new) ---
    %lFt = tv_time_lipschitz(xhat, t, obs);


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
    Mk = compute_margin_sup_TV(R_lo, R_hi, xhat, t, obs, u_min, u_max, gamma);


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
    t = k*T;      
    c = obs.c_fun(t);           % 2x1 vector [cx; cy]
    D = obs.D_fun(t);
    
    % --- barrier at current step ---
    h_curr = (X(1,k+1)-c(1))^2 + (X(2,k+1)-c(2))^2 - D^2;
    h_min  = min(h_min, h_curr);
    
    % --- (optional) compute cycle time if you wrapped tic/toc around TIRA+QP
    % comp_ms = 1000*toc(t0);
    comp_ms = 0;                % set to 0 if you didn't time this step
    
    % --- update visualization (R_lo/R_hi may be [] if not available) ---
    updateScene(Viz, X(:,k+1), Z(:,k+1), R_lo, R_hi, obs, t, h_curr, h_min, comp_ms);
    

    if norm(X(1:2,k+1)-goal) < 0.5, X=X(:,1:k+1); U=U(:,1:k); break; end
end

% ---- plots (after loop) ----
idx = all(~isnan(Z(1:2,:)), 1);    % keep samples where x,y measurements exist

figure; hold on; axis equal; grid on
th = linspace(0, 2*pi, 200);
co = obs.c_fun(t);      % evaluate center at time t  -> 2x1
cx = co(1); cy = co(2);
D  = obs.D_fun(t);      % evaluate radius at time t  -> scalar
plot(cx + D*cos(th), cy + D*sin(th), 'k--', 'LineWidth', 1.5);
%plot(X(1,:),X(2,:),'-','LineWidth',2)

plot(Z(1,idx), Z(2,idx), 'k-','LineWidth', 1.2)        % line

plot(goal(1),goal(2),'gx','LineWidth',2)
xlabel x; ylabel y;
xlim([0,50]);
ylim([0,50]);

title('Eq.(12) with TIRA margin — unicycle')

t_series = 0; hvals = zeros(1, size(X,2));
for j = 1:size(X,2)
    cc = obs.c_fun(t_series);
    Dj = obs.D_fun(t_series);
    dx = X(1,j)-cc(1); dy = X(2,j)-cc(2);
    hvals(j) = dx.^2 + dy.^2 - Dj.^2;    % h(x_j, t_j)
    t_series = t_series + T;             % if you log T_k, sum those instead
end
figure; plot(hvals,'LineWidth',1.8); yline(0,'k--'); xlabel step; ylabel('h(x,t)');
title('Barrier over time (TV-HOCBF, sampled-data)')

% Barrier check
dx = X(1,:)-obs.c_fun(1); dy = X(2,:)-obs.c_fun(2);
hvals = dx.^2 + dy.^2 - D^2;
figure; plot(hvals,'LineWidth',1.8); yline(0,'k--');
xlabel step; ylabel('h(x)'); title('Barrier over time (should stay ≥ 0)')

%% ===== Helper: projection onto halfspace ∩ box (no toolbox) =====
function [a, c] = fv_affine_at_TV(x, t, obs)
% Time-varying HOCBF (m=2) for unicycle + circular keep-out.
% a(x,t) = L_g L_f h(x,t) (1x2),  c(x,t) = L_f^2 h + 2 L_f h + 2 h  +  [ L_f h_t + (d/dt)(L_f h) + h_tt + 2 h_t ]

    xx=x(1); yy=x(2); th=x(3); v=x(4);

    % obstacle & radius and their time derivatives at t
    co   = obs.c_fun(t);      Dx  = co(1);      Dy  = co(2);
    cDot = obs.cdot_fun(t);   Dxdt = cDot(1);   Dydt = cDot(2);
    cDD  = obs.cddot_fun(t);  Dxd2 = cDD(1);    Dyd2 = cDD(2);
    D    = obs.D_fun(t);
    Dd   = obs.Ddot_fun(t);
    Ddd  = obs.Dddot_fun(t);

    % geometry
    dx = xx - Dx; dy = yy - Dy;
    phi = dx*cos(th) + dy*sin(th);       % projection along heading

    % barrier and standard Lie derivatives
    h    = dx*dx + dy*dy - D*D;
    Lfh  = 2*v*phi;                       % ∇h·f
    Lf2h = 2*v*v;                         % f·∇(Lfh)

    % input channel (same as before)
    a1 = 2*v*(-dx*sin(th) + dy*cos(th));  % ∂/∂u1 of Lf h via θ̇
    a2 = 2*phi;                           % ∂/∂u2 via v̇
    a  = [a1, a2];

    % ---- time-derivative pieces (new) ----
    % h_t
    h_t   = -2*(dx*Dxdt + dy*Dydt) - 2*D*Dd;

    % (d/dt)(L_f h) and L_f(h_t)  (for this structure they coincide)
    Lf_ht = -2*v*(Dxdt*cos(th) + Dydt*sin(th));
    dtLfh = -2*v*(Dxdt*cos(th) + Dydt*sin(th));

    % h_tt
    h_tt  = -2*(dx*Dxd2 + dy*Dyd2) + 2*(Dxdt*Dxdt + Dydt*Dydt) ...
            - 2*(Dd*Dd) - 2*D*Ddd;

    % original (time-invariant) c  +  time-varying correction
    c_TI = Lf2h + 2*Lfh + 2*h;
    %c_TVcorr = Lf_ht + dtLfh + h_tt + 2*h_t;

    c = c_TI ;
end


function Mk = compute_margin_sup_TV(R_lo, R_hi, xhat, t, obs, umin, umax, gamma)
% Paper-faithful Eq. (12):
% margin = sup_{x∈R, u∈U, d∈P} [ fv(x,u) - fv(xhat,u) + f_d(x,d) ]
% We over-approximate sup_{x∈R} by sampling the TIRA box [R_lo, R_hi]
% (implementation detail; the formula itself is exactly Eq. (12)).

    % Precompute fv pieces at xhat
    [a_hat, c_hat] = fv_affine_at_TV(xhat, t, obs);

    % Support function of the input box U for the affine term Δa(x)·u
    support_u = @(Da) max(Da(1)*umin(1), Da(1)*umax(1)) + ...
                       max(Da(2)*umin(2), Da(2)*umax(2));


    % Deterministic grid samples over the tube (covers corners)
    M = 256;
    xs = linspace(0,1,round(M^(1/4)));
    [X1,X2,X3,X4] = ndgrid(xs,xs,xs,xs);
    gridPts = [X1(:) X2(:) X3(:) X4(:)].';

    c_now = obs.c_fun(t);

    sup_all = -inf;
    for i = 1:size(gridPts,2)
        r = gridPts(:,i);
        x = R_lo + r.*(R_hi - R_lo);

        % Δfv(x,u) = [Δc(x) + Δa(x)·u] with u ∈ U (box)
        [a, c]  = fv_affine_at_TV(x, t, obs);
        Dc = c - c_hat;
        Da = a - a_hat;
        sup_u   = support_u(Da);

        % f_d(x,d) from Eq. (12), with exact box sup over d ∈ [-γ,γ]^2
        v = x(4); th = x(3); dx = x(1)- c_now(1); dy = x(2)- c_now(2);
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

function lFt = tv_time_lipschitz(xhat, t, obs)
% Conservative time-Lipschitz bound on the nominal left-hand side at (xhat,t):
% l_Ft ≈ |c(xhat,t+dt) - c(xhat,t)| / dt
% (Here c includes the TV corrections; a(x,t) has no explicit t for this model.)
    dt = 1e-3;   % small step (tune if needed)
    [~, c0] = fv_affine_at_TV(xhat, t,     obs);
    [~, c1] = fv_affine_at_TV(xhat, t+dt,  obs);
    lFt = abs(c1 - c0)/dt;
end

function H = initScene(x0, obs, T)
    figure('Color','w'); hold on; axis equal; grid on
    xlim([0, 60]); ylim([0, 50]);
    xlabel('x [m]'); ylabel('y [m]');
    title('Unicycle — robust SD-HOCBF (Eq. 12)')

    % obstacle (will update every step)
    th = linspace(0,2*pi,200);
    H.th = th;
    co = obs.c_fun(0); D = obs.D_fun(0);
    H.obs = plot(co(1)+D*cos(th), co(2)+D*sin(th),'k--','LineWidth',1.5);

    % paths
    H.path_true = animatedline('LineWidth',1.6,'Color',[0 0 0]);
    H.path_meas = animatedline('LineWidth',1.0,'Color',[0.6 0.6 0.6]);

    % robot glyph
    [vx,vy] = triangle_pose(x0(1), x0(2), x0(3), 0.7);
    H.rob = patch('XData',vx,'YData',vy,'FaceColor',[0.2 0.6 1],...
                  'EdgeColor','k','LineWidth',1);

    % TIRA (x,y) rectangle (translucent)
    H.tira = patch('XData',[],'YData',[],'FaceColor',[1 0.6 0],...
                   'FaceAlpha',0.2,'EdgeColor',[1 0.4 0],'LineWidth',1);

    % HUD text
    H.txt = text(1,49,sprintf('t=%.1f  h=%.2f  min h=%.2f  comp=%.0f ms',0,NaN,NaN,0),...
                 'FontName','Consolas','FontSize',10,'VerticalAlignment','top');

    drawnow;
end

function updateScene(H, x_true, z_meas, R_lo, R_hi, obs, t, h_curr, h_min, comp_ms)
    % obstacle
    c = obs.c_fun(t); D = obs.D_fun(t);
    set(H.obs,'XData',c(1)+D*cos(H.th),'YData',c(2)+D*sin(H.th));

    % paths
    addpoints(H.path_true, x_true(1), x_true(2));
    if ~isempty(z_meas)
        addpoints(H.path_meas, z_meas(1), z_meas(2));
    end

    % robot glyph
    [vx,vy] = triangle_pose(x_true(1), x_true(2), x_true(3), 0.7);
    set(H.rob,'XData',vx,'YData',vy);

    % TIRA endpoint set projected on (x,y)
    if ~isempty(R_lo) && ~isempty(R_hi)
        Xr = [R_lo(1) R_hi(1) R_hi(1) R_lo(1)];
        Yr = [R_lo(2) R_lo(2) R_hi(2) R_hi(2)];
        set(H.tira,'XData',Xr,'YData',Yr);
    end

    % HUD
    set(H.txt,'String',sprintf('t=%.1f  h=%.2f  min h=%.2f  comp=%d ms', ...
        t, h_curr, h_min, round(comp_ms)));

    drawnow limitrate
end

function [vx,vy] = triangle_pose(x, y, th, s)
    tri = s * [ 1  0 -0.6; 0  0.3 -0.3 ];
    R = [cos(th) -sin(th); sin(th) cos(th)];
    pts = R*tri + [x; y];
    vx = pts(1,:); vy = pts(2,:);
end
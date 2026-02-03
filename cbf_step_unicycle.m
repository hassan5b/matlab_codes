function [u, debug, u_perf1, a_hat1, c_hat1, Mk1, A1, b1, R_lo1, R_hi1, Tnext1] = cbf_step_unicycle(xhat, goal, obs, params)
% One-step sampled-data robust HOCBF safety filter (TV obstacle + variable sampling).
%
% States: x = [x; y; theta; v]
% Inputs: u = [theta_dot; v_dot]
%
% This version matches your working QBot integration constraints:
%   - Returns only fixed-size numerics (no structs) for Simulink logging
%   - Uses TIRA to compute a tube + robust margin M_k
%   - Adds TV correction via -lFt*T term (finite-difference time Lipschitz)
%   - Selects a variable next sampling time T_next (bisection rule) like your MATLAB script
%
% Inputs
%   xhat   : 4x1 measured state at time t_k
%   goal   : 2x1 goal position
%   obs    : struct with numeric motion parameters (no function handles), fields:
%              .c0    [2x1] base center
%              .c_amp [2x1] cosine amplitude for center
%              .c_w   scalar frequency (cos)
%              .D0    scalar base radius
%              .D_amp scalar amplitude (sin)
%              .D_w   scalar frequency (sin)
%   params : struct with fields
%              .T        sample time used for THIS solve
%              .t        current time (for TV obstacle evaluation)
%              .T_mpc    MPC discretization (kept fixed)
%              .u_min/.u_max box bounds on u
%              .gamma    disturbance bound
%              .epsM     measurement noise bound on x,y
%              .W        QP weight matrix
%              .beta, .Tmin, .Tmax_hw, .deltaRT, .g_release, .grow_max, .bisectIts
%
% Outputs
%   u      : safe input
%   debug  : quadprog exitflag
%   u_perf1, a_hat1, c_hat1, Mk1, A1, b1, R_lo1, R_hi1: debug scalars/vectors
%   Tnext1 : selected next sampling time (continuous; scheduler quantizes)

    % You still need System_description case 99 set up for TIRA
    global system_choice;
    if isempty(system_choice)
        system_choice = 99;
    end

    % ---- unpack params ----
    T     = params.T;
    t     = params.t;

    u_min = params.u_min;
    u_max = params.u_max;

    gamma = params.gamma;
    epsM  = params.epsM;
    W     = params.W;

    % MPC discretization (fixed)
    if isfield(params,'T_mpc')
        T_mpc = params.T_mpc;
    else
        T_mpc = T;
    end

    % variable sampling parameters (defaults match updated_matlab_cbf)
    beta      = getfield_def(params,'beta',0.85);
    Tmin      = getfield_def(params,'Tmin',0.05);
    Tmax_hw   = getfield_def(params,'Tmax_hw',0.12);
    deltaRT   = getfield_def(params,'deltaRT',0.0);
    g_release = getfield_def(params,'g_release',1);
    grow_max  = getfield_def(params,'grow_max',1.2);
    bisectIts = getfield_def(params,'bisectIts',8);

    % -------------------------------------------------------------
    % 1) Nominal controller K_perf (MPC) at xhat (FIXED discretization)
    % -------------------------------------------------------------
    u_perf = Kperf_MPC(xhat, goal, u_min, u_max, T_mpc);

    % -------------------------------------------------------------
    % 2) TV affine CBF pieces at (xhat,t): a_hat * u + c_hat
    % -------------------------------------------------------------
    [a_hat, c_hat] = fv_affine_at_TV(xhat, t, obs);

    % -------------------------------------------------------------
    % 3) Time-Lipschitz bound (finite difference) for TV correction
    % -------------------------------------------------------------
    lFt = tv_time_lipschitz(xhat, t, obs);

    % -------------------------------------------------------------
    % 4) TIRA tube over [0, T] around measurement xhat
    % -------------------------------------------------------------
    Z_lo = xhat + [-epsM; -epsM; 0; 0];
    Z_hi = xhat + [ epsM;  epsM; 0; 0];

    p_lo = [u_min; -gamma; -gamma];
    p_hi = [u_max;  gamma;  gamma];

    [R_lo, R_hi] = TIRA([0, T], Z_lo, Z_hi, p_lo, p_hi);

    % -------------------------------------------------------------
    % 5) Sampled-data margin M_k over the tube (TV version)
    % -------------------------------------------------------------
    Mk = compute_margin_sup_TV(R_lo, R_hi, xhat, t, obs, u_min, u_max, gamma);

    % -------------------------------------------------------------
    % 6) Solve QP: min ||u - u_perf||_W^2
    %    s.t. a_hat*u + c_hat - Mk - lFt*T >= 0
    %         u_min <= u <= u_max
    % -------------------------------------------------------------
    A = -a_hat;
    b =  (c_hat - Mk - lFt*T);

    H   = 2*W;
    f   = -2*W*u_perf;
    Aqp = [A;  eye(2); -eye(2)];
    bqp = [b;  u_max;   -u_min];

    opts = optimoptions('quadprog','Display','off');

    t_qp = tic;
    [u,~,flag] = quadprog(H, f, Aqp, bqp, [], [], [], [], [], opts);
    solve_time_k = toc(t_qp);

    % Fallback to saturated nominal if QP fails
    if isempty(u) || flag <= 0
        u = max(min(u_perf, u_max), u_min);
    end

    % -------------------------------------------------------------
    % 7) Variable sampling: choose T_next using feasibility checks
    % -------------------------------------------------------------
    g_now = a_hat*u + c_hat - Mk - lFt*T;

    if g_now > g_release
        T_next = min(Tmax_hw, max(Tmin, grow_max*T));
    else
        % bisection on T_max s.t. holding current u stays feasible
        if feasible_at_TV(Tmin, xhat, u, a_hat, c_hat, lFt, t, obs, u_min, u_max, gamma, epsM)
            T_lo = Tmin;
            T_hi = Tmax_hw;
            for it = 1:bisectIts
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

    % Optional real-time guard (kept from the MATLAB script)
    if deltaRT > 0
        T_next = max(T_next, solve_time_k + deltaRT);
        T_next = min(max(T_next, Tmin), Tmax_hw);
    end

    % -------------------------------------------------------------
    % Outputs
    % -------------------------------------------------------------
    debug  = flag;
    u_perf1 = u_perf;
    a_hat1  = a_hat;
    c_hat1  = c_hat;
    Mk1     = Mk;
    A1      = A;
    b1      = b;
    R_lo1   = R_lo;
    R_hi1   = R_hi;
    Tnext1  = T_next;
end

% =================================================================
% Helper functions (TV + variable sampling)
% =================================================================

function [a, c] = fv_affine_at_TV(x, t, obs)
% Time-varying HOCBF (relative degree 2) for unicycle + circular keep-out.
% a(x,t) = L_g L_f h(x,t) (1x2)
% c(x,t) = L_f^2 h + 2 L_f h + 2 h + [ L_f h_t + d/dt(L_f h) + h_tt + 2 h_t ]

    xx=x(1); yy=x(2); th=x(3); v=x(4);

    [co, cDot, cDD, D, Dd, Ddd] = obs_at(t, obs);
    Dx   = co(1);     Dy   = co(2);
    Dxdt = cDot(1);   Dydt = cDot(2);
    Dxd2 = cDD(1);    Dyd2 = cDD(2);

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

    c_TI    = Lf2h + 2*Lfh + 2*h;
    c_TVcorr = Lf_ht + dtLfh + h_tt + 2*h_t;

    c = c_TI + c_TVcorr;
end

function Mk = compute_margin_sup_TV(R_lo, R_hi, xhat, t, obs, umin, umax, gamma)
% Eq.(12) margin for the TV version.

    [a_hat, c_hat] = fv_affine_at_TV(xhat, t, obs);

    support_u = @(Da) max(Da(1)*umin(1), Da(1)*umax(1)) + ...
                      max(Da(2)*umin(2), Da(2)*umax(2));

    M = 256;
    xs = linspace(0,1,round(M^(1/4)));
    [X1,X2,X3,X4] = ndgrid(xs,xs,xs,xs);
    gridPts = [X1(:) X2(:) X3(:) X4(:)].';

    co = obs_at(t, obs);

    sup_all = -inf;
    for i = 1:size(gridPts,2)
        r = gridPts(:,i);
        x = R_lo + r.*(R_hi - R_lo);

        [a, c]  = fv_affine_at_TV(x, t, obs);
        Dc = c - c_hat;
        Da = a - a_hat;
        sup_u   = support_u(Da);

        v = x(4); th = x(3); dx = x(1)- co(1); dy = x(2)- co(2);
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

    xs = linspace(0,1,3);               % 3^4 = 81 samples
    [X1,X2,X3,X4] = ndgrid(xs,xs,xs,xs);
    gridPts = [X1(:) X2(:) X3(:) X4(:)].';

    co = obs_at(t, obs);

    sup_all = -inf;
    for i = 1:size(gridPts,2)
        r = gridPts(:,i);
        x = R_lo + r.*(R_hi - R_lo);

        [a, c] = fv_affine_at_TV(x, t, obs);
        Dc = c - c_hat;
        Du = (a - a_hat) * u;

        v = x(4); th = x(3); dx = x(1)- co(1); dy = x(2)- co(2);
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

function lFt = tv_time_lipschitz(xhat, t, obs)
% Conservative time-Lipschitz bound on c(xhat,t) via finite differences.
    dt = 1e-3;
    [~, c0] = fv_affine_at_TV(xhat, t,    obs);
    [~, c1] = fv_affine_at_TV(xhat, t+dt, obs);
    lFt = abs(c1 - c0)/dt;
end

function [co, cDot, cDD, D, Dd, Ddd] = obs_at(t, obs)
% Option-B obstacle model, numeric-parameter form (no function handles).
%
% Center:
%   c(t) = c0 + c_amp .* cos(c_w*t)
% Radius:
%   D(t) = D0 + D_amp * sin(D_w*t)

    c0    = obs.c0;
    c_amp = obs.c_amp;
    c_w   = obs.c_w;
    D0    = obs.D0;
    D_amp = obs.D_amp;
    D_w   = obs.D_w;

    co = [c0(1) + c_amp(1)*cos(c_w*t);
          c0(2) + c_amp(2)*cos(c_w*t)];

    cDot = [-c_amp(1)*c_w*sin(c_w*t);
            -c_amp(2)*c_w*sin(c_w*t)];

    cDD  = [-c_amp(1)*(c_w^2)*cos(c_w*t);
            -c_amp(2)*(c_w^2)*cos(c_w*t)];

    D   = D0 + D_amp*sin(D_w*t);
    Dd  = D_amp*D_w*cos(D_w*t);
    Ddd = -D_amp*(D_w^2)*sin(D_w*t);
end

% =================================================================
% MPC + RK4 helpers (same structure as your working version)
% =================================================================

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
        for tt = 1:Np
            u = U(2*tt-1:2*tt);
            x = rk4_unicycle(x, u, T);
            pos_err = x(1:2) - goal(:);
            J = J + Qp*(pos_err.'*pos_err) + u.'*Rw*u;
        end
        term_err = x(1:2) - goal(:);
        J = J + Qf*(term_err.'*term_err);
    end

    opts = optimoptions('fmincon','Display','off',...
                        'MaxIterations',100,'Algorithm','sqp');
    Uopt = fmincon(@cost, U0, [], [], [], [], lb, ub, [], opts);
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

function v = getfield_def(S, field, default)
% Small helper to safely read optional params fields.
    if isstruct(S) && isfield(S, field)
        v = S.(field);
    else
        v = default;
    end
end

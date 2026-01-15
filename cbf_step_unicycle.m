function [u, debug,u_perf1,a_hat1,c_hat1,Mk1,A1,b1,R_lo1,R_hi1] = cbf_step_unicycle(xhat, goal, obs, params)
% One-step SD-HOCBF controller for the unicycle system
% ... (rest of the function body and subfunctions) ...
% One-step SD-HOCBF controller for the unicycle system (relative degree 2)
% States: x = [x; y; theta; v]
% Inputs: u = [u1; u2] = [theta_dot; v_dot]
%
% xhat   : 4x1 measured state at time t_k
% goal   : 2x1 desired position
% obs    : struct with fields
%            .c  [2x1] obstacle center
%            .D  scalar obstacle radius
% params : struct with fields
%            .T      sample time
%            .u_min  2x1 lower bounds on [u1; u2]
%            .u_max  2x1 upper bounds
%            .gamma  disturbance bound
%            .epsM   measurement noise bound on x,y
%            .W      2x2 weight matrix in QP
%
% u      : 2x1 safe input [u1; u2],
% debug  : struct with some extra info (optional)

    % You still need System_description case 99 set up for TIRA
    global system_choice;
    if isempty(system_choice)
        system_choice = 99;   % your unicycle case in System_description
    end

    T      = params.T;
    u_min  = params.u_min;
    u_max  = params.u_max;
    gamma  = params.gamma;
    epsM   = params.epsM;
    W      = params.W;

    % -------------------------------------------------------------
    % 1) Nominal controller K_perf (MPC) at xhat
    % -------------------------------------------------------------
    u_perf = Kperf_MPC(xhat, goal, u_min, u_max, T);

    % -------------------------------------------------------------
    % 2) Local affine CBF pieces at xhat: a_hat * u + c_hat
    % -------------------------------------------------------------
    [a_hat, c_hat] = fv_affine_at(xhat, obs.c, obs.D, 1, 1);  % k1=k2=1

    % -------------------------------------------------------------
    % 3) TIRA tube over [0, T] around the measurement xhat
    % -------------------------------------------------------------
    % measurement set Z around xhat (uncertainty on x,y only)
    Z_lo = xhat + [-epsM; -epsM; 0; 0];
    Z_hi = xhat + [ epsM;  epsM; 0; 0];

    % p = [u1; u2; d1; d2] ranges (constant over [0,T])
    p_lo = [u_min; -gamma; -gamma];
    p_hi = [u_max;  gamma;  gamma];

    % Compute continuous-time reachable tube R([0,T]) using your
    % System_description case 99 inside TIRA
    [R_lo, R_hi] = TIRA([0, T], Z_lo, Z_hi, p_lo, p_hi);

    % -------------------------------------------------------------
    % 4) Sampled-data margin M_k over the tube (Eq. 12)
    % -------------------------------------------------------------
    Mk = compute_margin_sup(R_lo, R_hi, xhat, obs.c, obs.D, u_min, u_max, gamma);
    
    % -------------------------------------------------------------
    % 5) Solve QP: min ||u - u_perf||_W^2
    %     s.t. a_hat*u + c_hat - M_k >= 0   (SD-HOCBF)
    %          u_min <= u <= u_max
    % -------------------------------------------------------------
    % inequality a_hat*u + c_hat - M_k >= 0
    %  ->  [-a_hat] u <= -(c_hat - M_k)
    A = -a_hat;
    b =  (c_hat - Mk);

    H   = 2*W;
    f   = -2*W*u_perf;
    Aqp = [A;  eye(2); -eye(2)];
    bqp = [b;  u_max;   -u_min];

    opts = optimoptions('quadprog','Display','off');
    [u,~,flag] = quadprog(H, f, Aqp, bqp, [], [], [], [], [], opts);

    % Fallback to saturated nominal if QP fails
    if isempty(u) || flag <= 0
        u = max(min(u_perf, u_max), u_min);
    end
    % Optional debug info
    if nargout > 1
        debug=flag;
        u_perf1 = u_perf;
        a_hat1  = a_hat;
        c_hat1  = c_hat;
        Mk1     = Mk;
        A1      = A;
        b1      = b;
        R_lo1   = R_lo;
        R_hi1   = R_hi;
    end
end

% =================================================================
% Subfunctions copied from your original script
% =================================================================

function [a, c] = fv_affine_at(x, cObs, D, k1, k2) %#ok<INUSD>
% a(x) = L_g L_f h(x) (1x2 row)
% c(x) = L_f^2 h(x) + 2 L_f h(x) + 2 h(x)
    xx = x(1); yy = x(2); th = x(3); v = x(4);
    dx = xx - cObs(1); dy = yy - cObs(2);
    phi = dx*cos(th) + dy*sin(th);
    h   = dx*dx + dy*dy - D^2;
    Lfh = 2*v*phi;
    Lf2h = 2*v^2;
    a1 = 2*v*(-dx*sin(th) + dy*cos(th));
    a2 = 2*phi;
    a  = [a1, a2];
    c  = Lf2h + 2*Lfh + 2*h;
end

function Mk = compute_margin_sup(R_lo, R_hi, xhat, cObs, D, umin, umax, gamma)
% Sample based over-approximation of the margin in Eq. (12).

    [a_hat, c_hat] = fv_affine_at(xhat, cObs, D, 1, 1);

    uabs = [max(abs(umin(1)),abs(umax(1))), ...
            max(abs(umin(2)),abs(umax(2)))];

    M = 128;
    xs = linspace(0,1,round(M^(1/4)));
    [X1,X2,X3,X4] = ndgrid(xs,xs,xs,xs);
    gridPts = [X1(:) X2(:) X3(:) X4(:)].';

    sup_all = -inf;
    for i = 1:size(gridPts,2)
        r = gridPts(:,i);
        x = R_lo + r.*(R_hi - R_lo);

        [a, c] = fv_affine_at(x, cObs, D, 1, 1);
        Dc = c - c_hat;
        Da = a - a_hat;

        % sup over u in the box
        sup_u = max(Da(1)*umin(1), Da(1)*umax(1)) + ...
                max(Da(2)*umin(2), Da(2)*umax(2));

        % disturbance term, exact sup over d in [-gamma,gamma]^2
        v = x(4); th = x(3); dx = x(1)-cObs(1); dy = x(2)-cObs(2);
        c1 = 4*(v*cos(th) + dx);
        c2 = 4*(v*sin(th) + dy);
        sup_d = gamma*(abs(c1) + abs(c2)) + 4*gamma^2;

        sup_all = max(sup_all, Dc + sup_u + sup_d);
    end
    Mk = sup_all;
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
        for t = 1:Np
            u = U(2*t-1:2*t);
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

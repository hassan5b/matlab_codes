function [v_cmd, w_cmd] = cbf_to_vw_wrapper(x,y,th,v, goal, obs_c, obs_D, T, u_min, u_max, gamma, epsM, W)

% Build obs and params as structs locally (no struct signals in Simulink)
obs = struct();
obs.c = obs_c;
obs.D = obs_D;

params = struct();
params.T     = T;
params.u_min = u_min;
params.u_max = u_max;
params.gamma = gamma;
params.epsM  = epsM;
params.W     = W;

% State for controller
xhat = [x; y; th; v];
dx = goal(1) - x;
dy = goal(2) - y;
if (dx*dx + dy*dy) < (0.2^2)   % 15 cm tolerance, tune 0.1–0.25
    v_cmd = 0;
    w_cmd = 0;
    return
end

% Call your existing controller (uses fmincon/quadprog/TIRA/etc.)
u = cbf_step_unicycle(xhat, goal, obs, params);   % u = [theta_dot; v_dot]

% Map to body-speed command
w_cmd = u(1);
v_cmd = v + u(2)*T;

% Clamp to QBot limits (these are [v; w] bus limits)
v_cmd = min(max(v_cmd, -0.35), 0.35);
w_cmd = min(max(w_cmd, -1.7825), 1.7825);

end

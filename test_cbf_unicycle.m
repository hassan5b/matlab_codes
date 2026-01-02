clear; clc; close all;

global system_choice;
system_choice = 99;        % your unicycle case in System_description

params.T     = 0.1;
params.u_min = [-1; -2];
params.u_max = [ 1;  2];
params.gamma = 0.3;
params.epsM  = 0.5;
params.W     = eye(2);

obs.c = [32; 25];
obs.D = 5;

x_true = [5;25;-pi/2;0.8];
goal   = [45;21];

N = 120;
X = zeros(4,N+1);
X(:,1) = x_true;

for k = 1:N
    xhat = x_true;     % later you can add measurement noise here

    [u, dbg] = cbf_step_unicycle(xhat, goal, obs, params);

    % propagate the plant one step (use your existing RK4 + disturbance)
    T = params.T;
    x_true = x_true + T * [ x_true(4)*cos(x_true(3));
                            x_true(4)*sin(x_true(3));
                            u(1);
                            u(2) ];
    dx = x_true(1) - obs.c(1);
    dy = x_true(2) - obs.c(2);
    h_log(k) = dx*dx + dy*dy - obs.D^2;

    X(:,k+1) = x_true;
end
figure; hold on; axis equal; grid on;
plot(X(1,:), X(2,:), 'LineWidth', 1.5);

% obstacle circle
th = linspace(0, 2*pi, 200);
plot(obs.c(1) + obs.D*cos(th), obs.c(2) + obs.D*sin(th), 'r', 'LineWidth', 2);

% start and goal
plot(X(1,1), X(2,1), 'go', 'MarkerSize', 8, 'LineWidth', 2);
plot(goal(1), goal(2), 'kx', 'MarkerSize', 10, 'LineWidth', 2);

legend('trajectory','obstacle','start','goal');
figure; grid on;
plot(h_log, 'LineWidth', 1.5);
yline(0,'r--');
title('h(x) over time'); xlabel('k'); ylabel('h');

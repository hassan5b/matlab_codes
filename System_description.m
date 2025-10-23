%% Descripion of the time-varying system with input
% either continuous-time dynamics: dx/dt = f(t,x,p)
% or discrete-time successor: x^+ = f(t,x,p)

% List of inputs
%   t: time
%   x: state
%   p: input value

% List of outputs
%   dx: 
%    continuous-time: vector field evaluated at time t, state x and input p
%    discrete-time: one-step successor from state x at time t with input p  

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 19th of February 2019

function dx = System_description(t,x,p)
% ...
global system_choice                 % <-- add this line

n_x = length(x);
dx   = NaN(n_x,1);

switch system_choice
    case 99
        th = x(3); v = x(4);
        u1 = p(1); u2 = p(2); d1 = p(3); d2 = p(4);
        dx = [ v*cos(th) + d1;
               v*sin(th) + d2;
               u1;
               u2 ];
    % (keep other cases as they were)
    otherwise
        error('System_description: unsupported system_choice=%g', system_choice);
end
end



%% User-provided system description
% For continuous-time system: dx is the time derivative dx/dt
% For discrete-time system: dx is the one-step successor x^+



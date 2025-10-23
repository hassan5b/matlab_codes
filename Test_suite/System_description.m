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
%   Neelay Junnarkar, <neelay.junnarkar -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 2nd of December 2021

function dx = System_description(t,x,p)

%% Default value as NaN vector (not a number)
n_x = length(x);
dx = NaN(n_x,1);

%% User-provided system description
% For continuous-time system: dx is the time derivative dx/dt
% For discrete-time system: dx is the one-step successor x^+
global system_choice
global u        % Some systems need an external control input

switch system_choice
    case 1
        %% Unicycle (continuous-time)
        % A simple example of a non-holonomic system
        % This system represents a point mass moving in the plane
        % The speed of the point and its angular velocity can be controlled
        
        % states:
        %     x(1), x(2): planar (x,y) coordinates of the unicycle
        %     x(3): unicycle facing angle
        
        % control:
        %     u(1): unicycle speed
        %     u(2): unicycle angular velocity
        
        % inputs:
        %     p: disturbance inputs to the three states
        
        dx = [u(1)*cos(x(3))+p(1); ...
              u(1)*sin(x(3))+p(2); ...
              u(2)+p(3)];

    case 4
        %% Traffic diverge 3 link (discrete-time)
        % This system is a simplified representation of the flow of traffic
        % through a three-way intersection, using the cell transmission
        % model. Here, x represents the traffic flow, in vehicles/period,
        % and the links are organized as follows:
        %           /-----(2)---
        %   (1)-----|
        %           \-----(3)---
        % Vehicles move from link 1 to links 2 and 3, making this a
        % "diverge" link. 
        % A more detailed description of this model is available in:
        % Samuel Coogan and Murat Arcak, "A Benchmark Problem in 
        % Transportation Networks". arxiv:1803.00367
        
        % states:
        %     x(i): traffic flow through the i-th link
        
        % inputs:
        %     p: traffic flow injected into the first link

        
        % Parameters
        v = 0.5;            % free-flow speed, in links/period
        w = 1/6;            % congestion-wave speed, in links/period
        c = 40;             % capacity (max downstream flow), in vehicles/period
        xbar = 320;         % max occupancy when jammed, in vehicles

        % Discrete-time successor
        dx = x + [p-min([c ; v*x(1) ; 2*w*(xbar-x(2)) ; 2*w*(xbar-x(3))]); ...
                  min([c/2 ; v*x(1)/2 ; w*(xbar-x(2)) ; w*(xbar-x(3))])-min(c , v*x(2)); ...
                  min([c/2 ; v*x(1)/2 ; w*(xbar-x(2)) ; w*(xbar-x(3))])-min(c , v*x(3))];
    case 6
        %% Time-varying system (continuous-time)
        % This is a simple linear periodic time-varying system, 
        % without a particular physical interpretation attached to it.
        % dx = A(t)*x+p
        
        % states:
        %     x(1), x(2): nonphysical
        
        % inputs: 
        %     p: disturbance inputs to the two states
        
        a = 0.8;
        A = [sin(t) cos(t) a;0 sin(t) cos(t);0 0 sin(t)];
        dx = A*x+p;
    case 7
        %% Bicycle (continuous-time)
        % Another cousin of the unicycle system, the Bicycle system models 
        % the planar motion of a vehicle with two wheels.
        % For a more detailed description of this model, see:
        % Majid Zamani, Giordano Pola, Manuel Mazo and Paulo Tabuada,
        % "Symbolic models for nonlinear control systems without stability 
        % assumptions". IEEE Transactions on Automatic Control, v. 57(7), 
        % pp. 1804-1809, 2012. DOI: 10.1109/TAC.2011.2176409
        
        % states:
        %     x(1), x(2): planar (x,y) coordinates of the bicycle
        %     x(3): heading of the bicycle.
        
        % control:
        %     u(1): bicycle forward velocity.
        %     u(2): bicycle angular velocity.
        
        % inputs:
        %     p: disturbance inputs to the three states

        a = atan(tan(u(2))/2);
        dx = [u(1)*cos(a+x(3))/cos(a);u(1)*sin(a+x(3))/cos(a);u(1)*tan(u(2))] + p;
    case 8
        %% Protein interaction (continuous-time)
        % This system is a simplified model of a positive feedback
        % interaction between protein and mRNA. In this case, the protein
        % encoded by a gene stimulates transcription through the nonlinear
        % term below. The salient features of this system are that it is 
        % monotone, and that for a and b values where a*b<0.5 the system 
        % exhibits a "bistable switch" behaviour, where trajectories will 
        % go to one of two stable attractors.
        
        % For a detailed discussion of protein kinematics models of this
        % kind, see the book:
        % Lee A. Segel, "Modeling dynamic phenomena in molecular and 
        % cellular biology". Cambridge University Press, 1984.
        
        % states:
        %    x(1): protein concentration
        %    x(2): mRNA concentration
        
        % inputs: none.
        
        a = 1;   % protein decay rate
        b = 0.3; % mRNA decay rate
        dx = [-a*x(1)+x(2);x(1)^2/(1+x(1)^2)-b*x(2)];
    case 11
        %% Pursuer-evader game with 2 Dubin's vehicles (continuous-time)
        % Stack of 2 unicycle models as in case 1, 
        % but where the inputs of the first unicycle are controlled
        % while those of the second one are disturbances

        % states:
        %     x(1), x(2): planar (x,y) coordinates of the unicycle 1
        %     x(3): unicycle 1 facing angle
        %     x(4), x(5): planar (x,y) coordinates of the unicycle 2
        %     x(6): unicycle 2 facing angle
        
        % control:
        %     u(1): unicycle 1 speed
        %     u(2): unicycle 1 angular velocity
        
        % inputs:
        %     p(1): unicycle 2 speed
        %     p(2): unicycle 2 angular velocity
        
        dx = [u(1)*cos(x(3)); ...
              u(1)*sin(x(3)); ...
              u(2); ...
              p(1)*cos(x(6)); ...
              p(1)*sin(x(6)); ...
              p(2)];
    case 12
        %% Simple oscillator model (continuous-time)
        % The trajectories of this system are circles centered on 0 (when
        % the input p is [0;0])
        
        % states:
        %     x(1), x(2): nonphysical
        
        % inputs: 
        %     p: disturbance inputs to the two states
        
        dx = [x(2);-x(1)]+p;
    case 14
        %% Simple unstable model (discrete-time)
        
        % states:
        %     x(1), x(2): nonphysical
        
        % inputs:
        %     p: disturbance inputs to the two states
        
        % Discrete-time successor
        a = 2;
        dx = a*x + p;
end

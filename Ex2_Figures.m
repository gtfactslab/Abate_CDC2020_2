
% Title: Computing Robustly Forward Invariant Sets for Mixed-Monotone
%        Systems
% To appear in 2020 IEEE 59th Conference on Decision and Control (CDC)
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 8/18/2020
% Description:  This script generates Figure 2.
%               An invariant rectangle is computed via analysis of 
%               equilibria in the embedding space.  The smallest robustly 
%               forward invariant set in X is also computed via exhaustive
%               simulation.

clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global W
W = [-2, 2]; % Didurbance Bound


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Hyperrectangular RFI Set Using MM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute an equilibrium in embedding space.
% This point defines a robustly forward invairant set for the original
% dynamics.
thing = fsolve(@E, zeros(4, 1));
XE = reshape(thing, 2, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Smallest Attractive RFI Set Via Exhaustive Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = .002;   % Time step for simulation
T  = 1;      % Simulation time-horizon

% Simulate all points on the boundary of XE forward in time using every
% possible disturbance.
XE_Boundary = makeRectangle(XE);
REACH = XE_Boundary;
T_size = size(0:dt:T, 2);
for t = 1:T_size
    t
    % get next FORWARD TIME reachable set
    holder = [];
    for i = 1:size(REACH, 2)
            x = REACH(:, i);
            for w = W(1):1:W(2)
                x_next = x + dt*dxdt(x, w);
                holder = [holder, x_next];
            end
    end
    k = boundary(holder(1, :)',holder(2, :)',0.8);
    REACH = holder(:, k);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
hold on; grid on;

% Plot RFI set from MM approach
patch(XE_Boundary(1, :), XE_Boundary(2, :), 'r', ...
            'LineWidth', 1.15, ...
            'FaceAlpha', .2, ...
            'HandleVisibility', 'off');
scatter(XE(1, :), XE(2, :), 'k', 'filled', ...
            'HandleVisibility', 'off');
% Plot RFI set from exhaustive simulation
patch(REACH(1, 1:8:end), REACH(2, 1:8:end), 'g', ...
            'LineWidth', 1.25, ...
            'FaceAlpha', .8, ...
            'HandleVisibility', 'off');

xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-1.5, 1.5, -2.25, 2.25]);
xticks([-1 0 1])
yticks([-2 -1 0 1 2])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = dxdt(x, w)
    out = [-x(1) - x(1)^3 - x(2) - w; ...
           -x(2) - x(2)^3 + x(1) + w^3]; 
end

function out = E(in)
    global W
    x(1) = in(1);
    x(2) = in(2);
    xh(1) = in(3);
    xh(2) = in(4);
    
    out(1:2, 1) = [-x(1) - x(1)^3 - xh(2) - W(2); ...
                   -x(2) - x(2)^3 + x(1) + W(1)^3]; 
    out(3:4, 1) = [-xh(1) - xh(1)^3 - x(2) - W(1); ...
                   -xh(2) - xh(2)^3 + xh(1) + W(2)^3];  
end

function out = makeRectangle(X0)
    d = [X0(1, 2) - X0(1, 1); X0(2, 2) - X0(2, 1)];
    [X0_x, X0_y] = meshgrid(X0(1, 1): d(1)/5 :X0(1, 2), ...
                            X0(2, 1): d(2)/5 :X0(2, 2));
    X_int = [X0_x(:), X0_y(:)];
    [k,av] = convhull(X_int);
    out = X_int(k, :)';
end

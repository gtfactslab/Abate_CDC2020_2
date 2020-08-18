
% Title: Computing Robustly Forward Invariant Sets for Mixed-Monotone
%        Systems
% To appear in 2020 IEEE 59th Conference on Decision and Control (CDC)
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 8/18/2020
% Description:  This script generates Figure 3.
%               Forward and backward-time reachable sets are approximated
%               via the MM property.

clc; clear all;

% Initial Set
X0 = [ -.25 ,.25; ....
       -.5 , 0 ];
% Disturbance Bound
W = [0, .25];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = .002;   % Timestep for simulation
T  = 1;      % Simulation time-horizon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X0_Boundary = makeRectangle(X0);
Phi0 = X0_Boundary;
Phi_size = size(Phi0, 2);

Phi = Phi0;
nPhi = Phi0;

holder = Phi;
nholder = Phi;

T_size = size(0:dt:T, 2);
xu = zeros(2, T_size + 1);  xu(:, 1) = X0(:, 1);
xo = zeros(2, T_size + 1);  xo(:, 1) = X0(:, 2);

nxu = zeros(2, T_size + 1);  nxu(:, 1) = X0(:, 1);
nxo = zeros(2, T_size + 1);  nxo(:, 1) = X0(:, 2);
 
for t = 1:T_size
    t
    % get next FORWARD TIME reachable set
    holder2 = [];
    for i = 1:size(holder, 2)
            x = holder(:, i);
            for w = W(1):.25:W(2)
                x_next = x + dt*dxdt(x, w);
                holder2 = [holder2, x_next];
            end
    end
    k = boundary(holder2(1, :)',holder2(2, :)',0.02);
    holder = holder2(:, k);
    Phi = [Phi, holder];
    % propegate corners with decomp function
    xu(:, t + 1) = xu(:, t) + dt*g(xu(:, t), xo(:, t), W(1), W(2));
    xo(:, t + 1) = xo(:, t) + dt*g(xo(:, t), xu(:, t), W(2), W(1));
    
    
    
    % get next BACKWARD TIME reachable set
    nholder2 = [];
    for i = 1:size(nholder, 2)
            nx = nholder(:, i);
            for w = W(1):.25:W(2)
                nx_next = nx + dt*ndxdt(nx, w);
                nholder2 = [nholder2, nx_next];
            end
    end
    k = boundary(nholder2(1, :)', nholder2(2, :)',0.02);
    nholder = nholder2(:, k);
    nPhi = [nPhi, nholder];
    % propegate corners with decomp function
    nxu(:, t + 1) = nxu(:, t) + dt*ng(nxu(:, t), nxo(:, t), W(1), W(2));
    nxo(:, t + 1) = nxo(:, t) + dt*ng(nxo(:, t), nxu(:, t), W(2), W(1));
    
end
x_T = holder;
xu_T = xu(:, t + 1);
xo_T = xo(:, t + 1);

nx_T = nholder;
nxu_T = nxu(:, t + 1);
nxo_T = nxo(:, t + 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
hold on; grid on;

% Plot Initial Set X0
patch(X0_Boundary(1, :), X0_Boundary(2, :), 'r', 'LineWidth', 1.25);

% Plot forward time reachable set and MM approximation
SF = makeRectangle([xu_T, xo_T]);
patch(SF(1, :), SF(2, :), 'g', 'FaceAlpha', .1, 'LineWidth', 1.25);
patch(x_T(1, :), x_T(2, :), 'g', 'FaceAlpha', .7, 'LineWidth', 1.25);

scatter([xu_T(1, 1), xo_T(1, 1)], [xu_T(2, 1), xo_T(2, 1)], 'k', 'filled');

% Plot backward time reachable set and MM approximation
NSF = makeRectangle([nxu_T, nxo_T]);
patch(NSF(1, :), NSF(2, :), 'b', 'FaceAlpha', .05, 'LineWidth', 1.25);
patch(nx_T(1, :), nx_T(2, :), 'b', 'FaceAlpha', .6, 'LineWidth', 1.25);
scatter(nxu_T(1, 1), nxu_T(2, 1), 'k', 'filled');
scatter(nxo_T(1, 1), nxo_T(2, 1), 'k', 'filled');

xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-1.5, 1, -2.25, 1.5]);
xticks([-1 0 1])
yticks([-2 -1 0 1])

Leg = legend();
set(Leg,'visible','off');

grid on;
ax.Layer = 'top';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = g(x, xh, w, wh)
    out = [g1(x, xh, w, wh); ...
           x(1) + 1]; 
end

function out = ng(x, xh, w, wh)
    out = [ng1(x, xh, w, wh); ...
           - xh(1) - 1]; 
end

function out = g1(x, xh, w, wh)
    if x(1, 1) >= 0  
        out = x(1, 1)*x(2, 1) + w;   
    else
        out = x(1, 1)*xh(2, 1) + w;
    end
end

function out = ng1(x, xh, w, wh)
    if x(1, 1) >= 0  
        out = -x(1, 1)*xh(2, 1) - wh;   
    else
        out = -x(1, 1)*x(2, 1) - wh;
    end
end

function out = makeRectangle(X0)
    d = [X0(1, 2) - X0(1, 1); X0(2, 2) - X0(2, 1)];
    [X0_x, X0_y] = meshgrid(X0(1, 1): d(1)/5 :X0(1, 2), ...
                            X0(2, 1): d(2)/5 :X0(2, 2));
    X_int = [X0_x(:), X0_y(:)];
    [k,av] = convhull(X_int);
    out = X_int(k, :)';
end

function out = dxdt(x, w)
    out = [ x(1, 1)*x(2, 1) + w;...
            x(1, 1) + 1];
end

function out = ndxdt(x, w)
    out = [ - x(1, 1)*x(2, 1) - w;...
            - x(1, 1) - 1];
end

% Title: Computing Robustly Forward Invariant Sets for Mixed-Monotone
%        Systems
% To appear in 2020 IEEE 59th Conference on Decision and Control (CDC)
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 8/18/2020
% Description:  This script generates Figure 4.
%               An gloally attractive region is computed via from an
%               equilibrium in the embeddign space, and a robustly forward
%               invariant region is computed from an equilibrium in the
%               embedding space of the backward time dynamics.  The
%               smallest attractive set is also computed via exhaustive
%               simulation


clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global W
% Bound on disturbance
W = [-3/4, 3/4; ...
     -3/4, 3/4];
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute forward invariant regions using MM approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yeq = fsolve(@E, [-5; -5; 5; 5]);       % get forward time eq
yeq_back = fsolve(@EG, [-1; -1; 1; 1]); % get backward time eq

% create invariant regions
yeq = reshape(yeq, 2, 2)    % Invariant and attractive region
[a, b] = meshgrid(yeq(1, :), flip(yeq(2, :)));
X_EQ = [a(:)'; b(1, 1), b(2, 1), b(2, 2), b(1, 2)];

yeq_back = reshape(yeq_back, 2, 2)  % Invariant region
[a, b] = meshgrid(yeq_back(1, :), flip(yeq_back(2, :)));
X_EQB = [a(:)'; b(1, 1), b(2, 1), b(2, 2), b(1, 2)];

if prod(yeq(:, 1) >= yeq(:, 2)) ||...
   prod(yeq_back(:, 1) >= yeq_back(:, 2))
    error('Equilibrium does not have correct ordering');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute smallest attractive region via exhaustive simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global dt 
dt = .005;  % Simulation timestep 
T = 2;      % Simulation time-horizon

% Get set of initial conditions on boundary of attractive region
disc = 20;
[X0i, X0j] = meshgrid(linspace(yeq(1, 1), yeq(1, 2), disc), ...
                      linspace(yeq(2, 1), yeq(2, 2), disc) );
X0 = [X0i(:), X0j(:)]';
holder = [];
for i = 1:size(X0, 2)
    if X0(1, i) <= yeq_back(1, 1) || ...
       X0(1, i) >= yeq_back(1, 2) || ...
       X0(2, i) <= yeq_back(2, 1) || ...
       X0(2, i) >= yeq_back(2, 2)
   
       holder = [holder, X0(:, i)];
    end
end
X0 = holder; % Set of initial conditions

% Simulate initial conditions forward in time to get starting set
TRAJ = X0(:);
for t = 0:dt:T
    w = [0; 0];
    dxdxt = 0;
    for i = 1:size(X0, 2)
        dxdt(2*i-1:2*i, 1) = F(TRAJ(2*i-1:2*i, end), w); %1:2, 3:4, etc
    end
    TRAJ(:, end+1) = TRAJ(:, end) + dt*dxdt;
end

% Simulate starting set forward to get smallest attractive set
REACH = reshape(TRAJ(:, end), 2, []);
disc = 3;
[W1_points, W2_points] = meshgrid(linspace(W(1, 1), W(1, 2), disc), ...
                                  linspace(W(2, 1), W(2, 2), disc));
W_points = [W1_points(:)'; ...
            W2_points(:)'];

holder2 = REACH;
dt = .02; % new simulation timestep 
for t = 0:dt:5
    t
    holder = [];
    for i = 1:size(holder2, 2)
        for j = 1:size(W_points, 2)
            w = W_points(:, j);
            dxdt = F(holder2(:, i), w);
            z = .001;
            holder = [holder, z*round((holder2(:, i) + dt*dxdt)/z)];
        end
    end
    % holder contains all the points reachable from holder2
    
    holder2 = unique([holder]', 'rows')';
    if isequal(holder2, [])
        break
    end
    [inner, outer] = boundary2(holder2);
    holder2 = setdiff([inner, outer]', REACH', 'rows')';
    REACH = [REACH, holder2];
end
[inner_REACH, outer_REACH] = boundary2(REACH);
k = boundary(inner_REACH(1, :)', inner_REACH(2, :)');
inner_REACH = inner_REACH(:, k);
k = boundary(outer_REACH(1, :)', outer_REACH(2, :)');
outer_REACH = outer_REACH(:, k);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
hold on; grid on;
ax = gca

patch(X_EQ(1, :), X_EQ(2, :), 'w', ...
                              'HandleVisibility', 'off', ...
                              'FaceAlpha', 1);
patch(X_EQ(1, :), X_EQ(2, :), 'g', ...
                              'HandleVisibility', 'off', ...
                              'FaceAlpha', .5);

patch(outer_REACH(1, :), outer_REACH(2, :), 'w', ...
                                            'HandleVisibility', 'off', ...
                                            'FaceAlpha', 1);
patch(outer_REACH(1, :), outer_REACH(2, :), 'b',  ...
                                            'HandleVisibility', 'off', ...
                                            'FaceAlpha', .5); % outer ring

patch(inner_REACH(1, :), inner_REACH(2, :), 'w', ...
                                            'HandleVisibility', 'off', ...
                                            'FaceAlpha', 1);
patch(inner_REACH(1, :), inner_REACH(2, :), 'g', ...
                                            'HandleVisibility', 'off', ...
                                            'FaceAlpha', .5); % inner ring


patch(X_EQB(1, :), X_EQB(2, :), 'w', ...
                                'HandleVisibility', 'off', ...
                                'FaceAlpha', 1);
patch(X_EQB(1, :), X_EQB(2, :), 'r', ...
                                'HandleVisibility', 'off', ...
                                'FaceAlpha', .5);

scatter(yeq(1, :), yeq(2, :), 'k', 'filled', ...
                              'HandleVisibility', 'off');
scatter(yeq_back(1, :), yeq_back(2, :), 'k', 'filled', ...
                                        'HandleVisibility', 'off');

xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-1.75, 1.75, -1.75, 1.75]);
xticks([-1 0 1])
yticks([-2 -1 0 1])


Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = F(x, w)
    out = [-x(2) + x(1)*(4 - 4*x(1)^2 - x(2)^2) + w(1); ...
            x(1) + x(2)*(4 - x(1)^2 - 4*x(2)^2) + w(2)];
end 

function out = d(X, Xhat, w)
    x = X(1:2);
    xh = Xhat(1:2);

    out = [-xh(2) + x(1)*(4 - 4*x(1)^2) + w(1); ...
            x(1)  + x(2)*(4 - 4*x(2)^2) + w(2)] + ...
          [- x(1)* l(x(1), x(2), xh(2)); ...
           - x(2)* l(x(2), x(1), xh(1))];
end




function out = l(a, b, c)
    if (a >= 0  &&  0 >= b && b <= -c) || ...
       (a < 0 &&  0 <= b && b >= -c)
   
        out = b^2;
           
    elseif (a >= 0  && c >= 0 && b > -c) || ...
           (a < 0 && c <= 0 && b < -c)  
       
        out = c^2;
        
    elseif (a >= 0  && b > 0 && c < 0)  || ...
           (a < 0 && b < 0 && c > 0)
       
        out = b*c;
    end
end




function out = E(xm)
    global W
    x = xm(1:2); xhat = xm(3:4);
    out = [d(x, xhat, W(:, 1)); ...
           d(xhat, x, W(:, 2))];
end

function out = dG(X, Xhat, wh)
    x = X(1:2);
    xh = Xhat(1:2);
    
    out = [x(2)  - x(1)*(4 - 4*x(1)^2 - l(-x(1), x(2), xh(2))) - wh(1); ...
          -xh(1) - x(2)*(4 - 4*x(2)^2 - l(-x(2), x(1), xh(1))) - wh(2)];   
end

function out = EG(xm)
    global W
    x = xm(1:2); xhat = xm(3:4);
    out = [dG(x, xhat, W(:, 2)); ...
           dG(xhat, x, W(:, 1))];
end

function [inner, outer] = boundary2(set)
    angles = atan2(set(2, :), set(1, :)); % angles of all points in 'set'
    holder_in = [];
    holder_out = [];
    for theta = (-pi: .002: pi)
        % get points in holder with close angle to theta
        tt = min(round(abs(theta - angles), 3)); %approximate closest angle between holder and theta
        [i, j, val] = find(abs(theta - angles) <= tt + .01); % 

        if tt <=.1 && ~isequal(size(i, 2), 0)
            thing = set(:, j);
            [M, I1] = min(vecnorm(thing));
            [M, I2] = max(vecnorm(thing));
            holder_in = [holder_in, thing(:, I1)];
            holder_out = [holder_out, thing(:, I2)];
        end
    end
    inner = unique(holder_in', 'rows')';
    outer = unique(holder_out', 'rows')';
end







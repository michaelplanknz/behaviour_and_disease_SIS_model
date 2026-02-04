function traj = solveODE(ODE_handle, tSpan, IC)

% Solve ODE
[traj.t, Y] = ode45(ODE_handle, tSpan, IC);

% Extract epi variables from ODE solution
traj.SN = Y(:, 1);
traj.SB = Y(:, 2);
traj.B = Y(:, 3);
traj.S = Y(:, 1)+Y(:, 2);



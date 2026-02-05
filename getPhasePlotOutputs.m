function results = getPhasePlotOutputs(par, dx, dxquiv, dyquiv, pQuiv, tManifold, pert, opts)

% FUnction to calculate the quantities needed for the phase plots
% (vector field, nullclines, equilibria, stable/unstable manifolds)

R0 = par.Beta/par.Gamma;

% Evaluate vector field over a grid
sx = 0:dxquiv:1;
by = 0:dyquiv:1;
[SX, BY] = meshgrid(sx, by);
[U, V] = myVectorField(SX, BY, par);

% Scale arrows with a power p of their initial length so long
% vectors don't dominate too much
Z = sqrt(U.^2+V.^2);
U = U.*Z.^(pQuiv-1);
V = V.*Z.^(pQuiv-1);

% Get (B,S) coords for the S nullcline
Snull1 = 0:dx:1;
Bnull1 = 1/par.q * (1 - 1./(R0*Snull1));

% Get (B,S) coords for the B nullcline
Bnull2 = 0:dx:1;
Snull2 = 1 + (par.Tau*Bnull2.^2 - par.Alpha*Bnull2./(1-Bnull2))/par.Chi;

% Get coords for BDFE
Bstar = 0.5 * (1 + [1, -1]*sqrt(1-4*par.Alpha/par.Tau));

% Find endemic equilibria
% First get good initial conditions for root-solving by moving along
% the S nullcline and seeing where it crosses the B hullcline
nx = length(Snull1);
xr0 = [];
leftFlag = 1;
for ix = 1:nx
    % Only consider values in [0, 1]
    if Bnull1(ix) >= 0 & Bnull1(ix) <= 1
        leftFlagPrev = leftFlag;
        % Compute flag variable indicating whether the current
        % point on the S nullclnie is left of the B nullcline
        leftFlag = Snull1(ix) < interp1(Bnull2, Snull2, Bnull1(ix));
        if leftFlag ~= leftFlagPrev
            % Crossing detected - append initial condition to
            % xr0
            xr0 = [xr0, [0.5*(Snull1(ix-1)+Snull1(ix)); 0.5*(Bnull1(ix-1)+Bnull1(ix))] ];
        end
    end
end

% Number of EEs to look for (number of crossing pts up to max of 3):
nEqs = min(3, size(xr0, 2));
EE = nan(2, 3);
Ystable1 = [];
Ystable2 = [];
Yunstable1 = [];
Yunstable2 = [];
for iEq = 1:nEqs
    EE(:, iEq) = fsolve(@(x)my2Dsystem(0, x, par), xr0(:, iEq), opts );
end

% Find stable and unstable manifolds for the saddle equilibrium (if one exists) 
if nEqs > 1
    % Usually this is the 2nd of 3 EEs...
    saddleEq = EE(:, 2);
elseif nEqs == 1 & 1/par.q*(1-1/R0) < Bstar(2)
    % ...but there is also a case where BDFE- is a saddle (regime 6)
    saddleEq = [1; Bstar(2)];
else
    saddleEq = [];
end

if ~isempty(saddleEq)
    % Get Jacobian and calculate its eigenvalues and
    % eigenvectors
    J = getJacobian(saddleEq, par);
    [eV, eD] = eig(J);
    eD = diag(eD);
    iStable = find(real(eD) < 0);
    iUnstable = find(real(eD) > 0);
    eigStable = eV(:, iStable);
    eigUnstable = eV(:, iUnstable);

    % Solve for stable manifrold (run backwards in time)
    ICstable1 = saddleEq + pert*eigStable;
    ICstable2 = saddleEq - pert*eigStable;
    if ICstable1(1) <= 1
        [~, Ystable1] = ode45(@(t, y)(-my2Dsystem(t, y, par)), tManifold, ICstable1);
    end
    if ICstable1(2) <= 1
        [~, Ystable2] = ode45(@(t, y)(-my2Dsystem(t, y, par)), tManifold, ICstable2);                    
    end
    % Solve for unstable manifrold (run forwards in time)
    ICunstable1 = saddleEq + pert*eigUnstable;
    ICunstable2 = saddleEq - pert*eigUnstable;
    [~, Yunstable1] = ode45(@(t, y)(my2Dsystem(t, y, par)), tManifold, ICunstable1);
    [~, Yunstable2] = ode45(@(t, y)(my2Dsystem(t, y, par)), tManifold, ICunstable2);
end

% Truncate manifolds to stay inside [0,1]x[0,1]
Ystable1 = trunc01(Ystable1);
Ystable2 = trunc01(Ystable2);
Yunstable1 = trunc01(Yunstable1);
Yunstable2 = trunc01(Yunstable2);

% Store outputs in results structure
results.sx = sx;
results.by = by;
results.U = U;
results.V = V;
results.Snull1 = Snull1;
results.Bnull1 = Bnull1;
results.Snull2 = Snull2;
results.Bnull2 = Bnull2;
results.Bstar = Bstar;
results.EE = EE; 
results.Ystable1 = Ystable1;
results.Ystable2 = Ystable2;
results.Yunstable1 = Yunstable1;
results.Yunstable2 = Yunstable2;


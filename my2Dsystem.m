function f = my2Dsystem(t, y, par)

% Evaluate the right-hand sides of the 2D system ODEs, by evaluating
% myVectorField at the specified input and returning the ouput as a column
% vector

[dSdt, dBdt] = myVectorField(y(1, :), y(2, :), par);

f = [dSdt; dBdt];


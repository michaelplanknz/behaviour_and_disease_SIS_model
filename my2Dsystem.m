function f = my2Dsystem(t, y, par)

[dSdt, dBdt] = myVectorField(y(1, :), y(2, :), par);

f = [dSdt; dBdt];


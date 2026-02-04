function f = mySteadyState(x, par)

[dSdt, dBdt] = myVectorField(x(1), x(2), par);

f = [dSdt; dBdt];


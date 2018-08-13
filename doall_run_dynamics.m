lengths = [0.1, 0.3];
m = 0.1;
g = 9.81;
I = m * (lengths(1)^2 + lengths(2)^2)/12.0;
stiction_tolerance = 1.0e-3;
relative_tolerance = 1e-3;
mu = 0.25;
y0 = 1.0;
vx0 = -1.0;
vy0 = 0.0;
w0 = 0;
sim_time = 2.0;
nsteps = 200;

penetration_allowance = 1.0e-6;

% Estimate contact stiffness/damping
damping_ratio = 1.0;
k = m*g/penetration_allowance;
omega = sqrt(k/m);
time_scale = 1.0/omega;
d = damping_ratio * time_scale / penetration_allowance;

% Save parameters into a struct
params.lengths = lengths;
params.m = m;
params.I = I;
params.g = g;
params.stiction_tolerance = stiction_tolerance;
params.relative_tolerance = relative_tolerance;
params.k = k;
params.d = d;
params.mu = mu;

box_dynamics = @(t, x) box_xdot(x, params);

x0 = [0; y0; 0; 
      vx0; vy0; w0];

options = odeset('RelTol',1e-6,'Stats','on');
tspan = linspace(0, sim_time, nsteps);
[tt,xx] = ode23(box_dynamics,tspan,x0,options);

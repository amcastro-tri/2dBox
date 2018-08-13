function state_xdot = box_xdot(state_x, params)
% Params
lengths = params.lengths;
m = params.m;
I = params.I;
g = params.g;
k = params.k;
d = params.d;
mu = params.mu;
stiction_tolerance = params.stiction_tolerance;
relative_tolerance = params.relative_tolerance;

ev = relative_tolerance * stiction_tolerance;
ev2 = ev*ev;

% State
q = state_x(1:3);
v = state_x(4:6);
p_WBo = q(1:2);
%theta = q(3);
%v_WBo = v(1:2);
%w = v(3);

p_BoC_W = calc_contact_points(q, lengths);

[Jn, Jt] = calc_jacobians(p_BoC_W);

% Calc signed distance and distance rate.
x = zeros(4,1);
for ic = 1:4
    p_WC = p_WBo + p_BoC_W(:, ic);
    x(ic) = -p_WC(2);   
end
xdot = -Jn * v;

% Normal force
fn = calc_normal_force(x, xdot, k, d);

% Tangent velocities.
vt = Jt * v;

% Tangent (friction) forces.
ft = zeros(4, 1);
for ic=1:4
    vt_ic = vt(ic);
    slip = abs(vt_ic);
    mu_ic = stribeck_friction2(slip, mu, stiction_tolerance);
    
    sign = vt_ic / sqrt(vt_ic^2 + ev2);
    
    ft(ic) = -mu_ic * fn(ic) * sign;    
end

% Genralized forces
tau = [0; -m*g; 0] + Jn'*fn + Jt'*ft;

% Generalized accelerations
vdot = zeros(3,1);

vdot(1:2) = tau(1:2) / m;
vdot(3) = tau(3) / I;

% Output
state_xdot = [v; vdot];

function [fn, dfdx, dfdxdot] = calc_normal_force(x, xdot, k, d)

% Value of the normal force.
kv = k * (1+d*xdot);
kv_p = max(0, kv);
x_p = max(0, x);
fn = kv_p .* x_p;

% Gradient with respect to x, xdot.
Hx = (x >= 0);
Hkv = (kv_p >=0);


dfdx = kv_p.*Hx;
dfdxdot = k*d*Hkv.*x_p;

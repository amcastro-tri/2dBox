function fn = calc_normal_force(x, xdot, k, d)

kv = k * (1+d*xdot);
kv_p = max(0, kv);
x_p = max(0, x);
fn = kv_p .* x_p;

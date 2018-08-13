function [Jn, Jt] = calc_jacobians(p_BoC_W)

% Normal Jacobian.
Jn = zeros(4, 3);
for ic = 1:4
    Jn(ic, :) = [0.0, 1.0, p_BoC_W(1, ic)];
end

% Tangential Jacobian.
Jt = zeros(4, 3);
for ic = 1:4
    Jt(ic, :) = [1.0, 0.0, -p_BoC_W(2, ic)];
end



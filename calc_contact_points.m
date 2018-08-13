function p_BoC_W = calc_contact_points(q, lengths)
%p_WBcm = q(1:2);
theta = q(3);

% The four corners, in the body frame.
p_BoC = zeros(2, 4);
p_BoC(:, 1) = [-lengths(1); -lengths(2)] / 2;
p_BoC(:, 2) = [ lengths(1); -lengths(2)] / 2;
p_BoC(:, 3) = [ lengths(1);  lengths(2)] / 2;
p_BoC(:, 4) = [-lengths(1);  lengths(2)] / 2;

% Rotation matrix
c = cos(theta);
s = sin(theta);
xhat = [c; s];
yhat = [-s; c];
R_WB = [xhat, yhat];

% The four corners, in the world frame.
p_BoC_W = zeros(2, 4);
for ic=1:4
    p_BoC_W(:, ic) = R_WB * p_BoC(:, ic);
end

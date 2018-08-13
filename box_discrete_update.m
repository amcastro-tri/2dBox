function state_x = box_discrete_update(itime, state_x0, params)

% Params
lengths = params.lengths;
m = params.m;
I = params.I;
g = params.g;
k = params.k;
d = params.d;
mu = params.mu;
h = params.h;
stiction_tolerance = params.stiction_tolerance;
relative_tolerance = params.relative_tolerance;

ev = relative_tolerance * stiction_tolerance;
ev2 = ev*ev;

% State
q0 = state_x0(1:3);
v0 = state_x0(4:6);
p_WBo = q0(1:2);
%theta = q(3);
%v_WBo = v(1:2);
%w = v(3);

% Sizes
nc = 4;
nv = 3;

p_BoC_W = calc_contact_points(q0, lengths);

[Jn, Jt] = calc_jacobians(p_BoC_W);


% Calc signed distance at t = t0.
x0 = zeros(4,1);
for ic = 1:4
    p_WC = p_WBo + p_BoC_W(:, ic);
    x0(ic) = -p_WC(2);   
end

% Tangent velocities.
%vt = Jt * v;

% Tangent (friction) forces.
% ft = zeros(4, 1);
% for ic=1:4
%     vt_ic = vt(ic);
%     slip = abs(vt_ic);
%     mu_ic = stribeck_friction2(slip, mu, stiction_tolerance);
%     
%     sign = vt_ic / sqrt(vt_ic^2 + ev2);
%     
%     ft(ic) = -mu_ic * fn(ic) * sign;    
% end

% Genralized forces
tau = [0; -m*g; 0]; % + Jn'*fn + Jt'*ft;


% Newton-Rapshon loop
max_iters = 100;
v = v0;
Gn = zeros(nc, nv);
M = [m, 0, 0;
     0, m, 0;
     0, 0, I];
pstar = M*v0 + h*tau;
 
vn_err = 2*ev;
for it=1:max_iters
    
    % Normal velocities
    vn = Jn*v;
    
    % Penetration distance and rate of change.    
    xdot = -vn;
    x = x0 + h * xdot;
       
    % Normal force and gradients.
    [fn, dfdx, dfdxdot] = calc_normal_force(x, xdot, k, d);
    
    % Check for norm of residual here so we get to compute the forces
    % with the latest velocity update.
    if (vn_err < ev) 
        break;
    end

    % Residual
    R = M*v - pstar - h*Jn'*fn; 

    % Normal forces Jacobian. Gn = dfn/dv.
    for ic=1:nc
        Gn(ic, :) = -(h*dfdx(ic)+dfdxdot(ic))*Jn(ic, :);
    end
    
    % System Jacobian. J = dRdv.
    J = M - h*Jn'*Gn;
    
    if (norm(J-J') > 1.0e-16*norm(J))
        msg = ['Matrix J is not symmetric. J = ' sprintf('\n') sprintf('%f %f %f\n',J)];
        error(msg);
    end
    
    % Velocity update
    dv = -J\R;
    
    % Update velocities.
    v = v + dv;
    
    dvn = Jn * dv;
    vn_err = norm(dvn);
end

if (vn_err > ev)
    % If we are here is because the NR iteration faild. Abort.
    msg = sprintf('NR iteration did not converge.\n It: %d.\n vn_err: %g.\n Time step: %d\n', it, vn_err, itime);
    error(msg);
end

% Update the state
q = q0 + h*v;
state_x = [q; v];



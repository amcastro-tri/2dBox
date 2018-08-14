function [state_x, fn, ft, vn, vt, vn_err, vt_err] = box_discrete_update(itime, state_x0, params)

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
Gt = zeros(nc, nv);
M = [m, 0, 0;
     0, m, 0;
     0, 0, I];
pstar = M*v0 + h*tau;

% Compute rigid Delassus.
Wnn = Jn * (M \ Jn');
mtilde = pinv(Wnn);
 
vn_err = 2*ev;
for it=1:max_iters
    
    % Normal/tangential velocities
    vn = Jn*v;
    vt = Jt*v;
    
    % Penetration distance and rate of change.    
    xdot = -vn;
    x = x0 + h * xdot;       
       
    % Normal force and gradients.
    [fn, dfdx, dfdxdot] = calc_normal_force(x, xdot, k, d);
    
    % Friction forces and gradients.
    [ft, dft_dvt, dft_dfn] = calc_friction_force(vt, fn, params);
    
    % Check for norm of residual here so we get to compute the forces
    % with the latest velocity update.
    if (vn_err < ev && vt_err < ev) 
        break;
    end

    % Residual
    R = M*v - pstar - h*Jn'*fn - h*Jt'*ft;

    % Normal forces Jacobian. Gn = dfn/dv.
    for ic=1:nc
        Gn(ic, :) = -(h*dfdx(ic)+dfdxdot(ic))*Jn(ic, :);
    end
    
    % Friction forces Jacobian. Gt = dft/dv
    for ic=1:nc
        Gt(ic, :) = dft_dvt(ic)*Jt(ic, :) + dft_dfn(ic) * Gn(ic, :);  
    end
    
    % System Jacobian. J = dRdv.
    J = M - h*Jn'*Gn - h*Jt'*Gt;
    
    %if (norm(J-J') > 1.0e-16*norm(J))
    %    msg = ['Matrix J is not symmetric. J = ' sprintf('\n') sprintf('%f %f %f\n',J)];
    %    error(msg);
    %end
    
    % Velocity update
    dv = -J\R;
    
    dvn = Jn * dv;
    dvt = Jt * dv;
    
    % Limit velocity update.
    alpha = 1.0;
    for ic=1:nc
        vt0 = vt(ic);
        
        % We started in the "strong" gradients region. Safe.
        if (abs(vt0) < stiction_tolerance)
            continue
        end
        
        vt1 = vt(ic) + dvt(ic);
        
        % If velocity changes direction.
        if (vt0 * vt1 < 0)
            % Clamp to within the stiction region, keeping the sign.
            vt_alpha = sign(vt0) * stiction_tolerance/2;
            
            alpha_ic = (vt_alpha-vt0)/dvt(ic);
            if (alpha_ic < 0 || alpha_ic > 1)
                error('Alpha is outsude [0, 1]');
            end
            
            alpha = min(alpha, alpha_ic);
        end        
    end
    
    % Update velocities.
    v = v + alpha*dv;
        
    vn_err = norm(dvn);
    vt_err = norm(dvt);
end

if (vn_err > ev || vt_err > ev)
    % If we are here is because the NR iteration failed. Abort.
    msg = sprintf('NR iteration did not converge.\n It: %d.\n vn_err: %g.\n vt_err: %g. \n Time step: %d\n', it, vn_err, vt_error, itime);
    error(msg);
end

% Update the state
q = q0 + h*v;
state_x = [q; v];



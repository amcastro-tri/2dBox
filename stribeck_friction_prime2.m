function mu = stribeck_friction_prime2(v, mus, vs)
xv = v / vs;

mu = zeros(size(xv));
for i=1:length(xv)
    x = xv(i);

    if (x >= 1)
        mu(i) = 0;
    else
        mu(i) = mus * (2*(1-x));
    end    
end

mu = mu / vs;

end
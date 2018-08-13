function mu = stribeck_friction2(v, mus, vs)
xv = v / vs;

mu = zeros(size(xv));
for i=1:length(xv)
    x = xv(i);

    if (x >= 1)
        mu(i) = mus;
    else
        mu(i) = mus * (1-(x-1)^2);
    end    
end

end
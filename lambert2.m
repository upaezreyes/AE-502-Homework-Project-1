function z = lambert2(mu, dt, z0, A, r1, r2); 

error = 1e-8; %tolerance error
n_max = 50000; %number of iterations 

D = 1; %ratio
n = 0; %initial number of iterations
z = z0; 

while abs(D) > error & n <= n_max
    n = n + 1; %adding number of iterations
    C = stumpC(z); 
    S = stumpS(z); 
    y  = r1 + r2 + A.*((z.*S - 1)./(sqrt(C)));
    F = ((y./C).^(3./2)).*S + A.*sqrt(y) - sqrt(mu).*dt;
    if z == 0 
        y0 = r1 + r2 + A.*((z.*(1./6) - 1)./(sqrt(0.5)));
        dFdz = (sqrt(2)./40).*y0.^(3./2) + (A./8).*(sqrt(y0) ...
                + A.*sqrt(1./(2.*y0))); 
    else
        dFdz = ((y./C).^(3./2)).*((1./(2.*z)).*(C ...
               -(3./2).*(S./C)) + (3./4).*((S.^2)./C)) ...
               + (A./8).*(3.*(S./C).*sqrt(y)...
               + A.*sqrt(C./y)); 
    end
    D = F/dFdz; 
    z = z - D; 

end

end


function X = Kepler_Universal(mu, dt, r0, Vr0, X0, d)

n0 = 0; %initial iteration 
n_max = 1000; % max number of iterations 
error = 1e-8;  % error tolerance 
X = X0; %let first iteration of X be X0
D = 1; 

while abs(D) > error & n0 <= n_max
    z = d.*X.^2; %z value
    n0 = n0 + 1; %adding iterations
    C = stumpC(z); %call stumpC function
    S = stumpS(z); %call stumpS function
    
    A = ((r0.*Vr0)./sqrt(mu)); 
    B = (1 - d.*r0); 
    
    F = A.*(X.^2).*C + B.*(X.^3).*S + r0.*X - sqrt(mu).*dt; 

    dFdx = A.*X.*(1 - d.*(X.^2).*S) + (1 - d.*r0).*(X.^2).*C + r0; 
    
    D = F/dFdx; 
    X = X - D; 
end

end


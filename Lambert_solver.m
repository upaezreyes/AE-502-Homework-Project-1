clc 
clear variables 

% 3600 sec for Curtis' book
% 1200 sec for ex given by Prof. Siegfried
dt = 3600; % time
%dt = linspace(0,3600,100); 
mu = 398600; % gravity constant
z0 = 10; %initial z value 

% positions vectors from Curtis's book: 
r1_v = [5000 10000 2100]; % initial position
r2_v = [-14600 2500 7000]; % final position 

% position vectors given by Professor Siegfried
% r1_v = [5644 2830 4170]; % initial position
% r2_v = [-2240 7320 4980]; % final position 

r1 = sqrt(sum(r1_v.^2)); % magnitude of initial position
r2 = sqrt(sum(r2_v.^2)); % magnitude of final position 

r1xr2 = cross(r1_v, r2_v); %cross product of initial & final positions
dot_r1r2 = dot(r1_v, r2_v); %dot product of the initial & final positions
 

%for prograde trajectory: -------------------------------------------------
if r1xr2(3) >= 0
    delta_theta = acosd(dot_r1r2./(r1.*r2));  %change of true anomaly angle
else
    delta_theta = 360 - acosd(dot_r1r2./(r1.*r2)); 
end 

% for retrograde trajectory:
% -----------------------------------------------------------------------
% if r1xr2(3) >= 0
%     delta_theta = 360 - acosd(dot_r1r2./(r1.*r2));  %change of true anomaly angle
% else
%     delta_theta =  acosd(dot_r1r2./(r1.*r2));
% end 



%constant A: 
A = sind(delta_theta).*sqrt((r1.*r2)./(1-cosd(delta_theta))); 

z = lambert2(mu, dt, z0, A, r1, r2); 

% if z < 0 
%     fprintf('The Orbit is a Hyperbola because z is:');
%     disp(z);
% elseif z > 0 
%     fprintf('The Orbit is an Ellipse because z is:');
%     disp(z);
% else 
%     fprintf('The Orbit is an Parabola because z is:');
%     disp(z);
% end

C = stumpC(z); 
S = stumpS(z); 

%calcualte y(y); 
y  = r1 + r2 + A.*((z.*S - 1)./(sqrt(C)));

%calculate the F & G functions 
f = 1 - (y./r1); 
g = A.*sqrt(y./mu); 
gdot = 1 - (y./r2); 

%=============================================================
%calculate the initial and final velocities 
V1 = (1./g).*(r2_v - f.*r1_v); %inital velocity vector
V2 = (1./g).*(gdot.*r2_v - r1_v); %final velocity vector
%==========================================================

fprintf('The initial velocity vector [i j k](km) is:')
disp(V1)
fprintf('The final velocity vector [i j k](km) is:')
disp(V2)

% deltaV = V2 - V1; 
% dV = sqrt(sum(deltaV.^2));
% fprintf('The deltaV is:')
% disp(dV)

%=========================================================================
%Orbital Elements from the initial state position & velocity vectors
%========================================================================
fprintf('Orbital Elments from r0 & v0:')
%1.) distance: was already calculated 
fprintf('\ndistance (km), r = ')
disp(r1)

%2.) speed:
V1_mag = sqrt(sum(V1.^2));
fprintf('Speed (km/s), v =')
disp(V1_mag); 

%3.) radial velocty:  
dot_r1v1 = dot(r1_v, V1); 
Vr = dot_r1v1./r1; 
fprintf('Radial Velocity (km/s), Vr =')
disp(Vr)

%4.) & 5.) specific angular momentum: 
h_v = cross(r1_v, V1); 
h = sqrt(sum(h_v.^2)); 
fprintf('Specific Angular Momentum (km^2/s), h =')
disp(h)

%6.) inclination: 
i = acosd(h_v(3)./h); 
fprintf('Inclination (degrees), i =')
disp(i)

%7.), 8.) & 9.) ascension of the ascending node:
K_v = [0 0 1]; %k unit vector
N_v = cross(K_v, h_v); % N vector
N = sqrt(sum(N_v.^2)); % magnitude of N vector 
if N_v(2) > 0 
    B = acosd(N_v(1)./N); 
else 
    B = 360 - acosd(N_v(1)./N); 
end
fprintf('Ascension of the Ascending Node (degrees) =')
disp(B)

%10.) & 11.) eccentricity: 
e_v = (1./mu).*((V1_mag.^2 - (mu./r1)).*r1_v - r1*Vr.*V1); 
e = sqrt(sum(e_v.^2)); 
fprintf('Eccentricity, e =')
disp(e)

%12.) argument of perigee: 
if e_v(3) >= 0
    w = acosd(dot(N_v./N, e_v./e)); 
else 
    w = 360 - acosd(dot(N_v./N, e_v./e))
end
fprintf('Argument of Perigee (degrees), w =')
disp(w)

%13.) true anomaly:
if Vr >= 0
    theta = acosd(dot(e_v./e, r1_v./r1)); 
else 
    theta = 360 - acosd(dot(e_v./e, r1_v./r1));
end 
fprintf('True Anomaly (degrees), theta =')
disp(theta)

%14.) semi-major axis
a = (h.^2)./(mu.*(1 - e.^2)); 
fprintf('Semi-Major Axis (km), a =')
disp(a)
%hw/proj 1
%question 6: 

clc
clear variables 

%=============================================================
% Oumuamua: 
%=============================================================
r1 = [3.515868886595499e-2,
      -3.162046390773074,
      4.493983111703389]*1.496e8; % km, velocity vector
V1 = [-2.317577766980901e-3,
      9.843360903693031e-3,
      -1.541856855538041e-2]*1731; % km/s, position vector

%=============================================================
% Bosivo:
% =================================================================
% r1 = [7.249472033259724,
%       14.61063037906177,
%       14.24274452216359]*1.496e8; % km, velocity vector
% V1 = [-8.241709369476881e-3,
%       -1.156219024581502e-2,
%       -1.317135977481448e-2]*1731; % km/s, position vector

%=========================================================================
%Orbital Elements from the initial state position & velocity vectors
%========================================================================
mu = 1.327e11; 
fprintf('Orbital Elments from initial state vectors:')

%1.) distance: was already calculated 
r1_mag = sqrt(sum(r1.^2)); % magnitude of position vector
fprintf('\nDistance (km), r = ')
disp(r1_mag)

%2.) speed:
V1_mag = sqrt(sum(V1.^2)); % magnitude of velocity vector 
fprintf('Speed (km/s), v =')
disp(V1_mag); 

%3.) radial velocty:  
dot_r1v1 = dot(r1, V1); 
Vr = dot_r1v1./r1_mag; 
fprintf('Radial Velocity (km/s), Vr =')
disp(Vr)

%4.) & 5.) specific angular momentum: 
h_v = cross(r1, V1); 
h = sqrt(sum(h_v.^2)); 
fprintf('Specific Angular Momentum (km^2/s), h =')
disp(h)

%6.) inclination: 
i = acosd(h_v(3)./h); 
fprintf('Inclination (degrees), i =')
disp(i)

%7.), 8.) & 9.) ascension of the ascending node:
K_v = [0,
       0,
       1]; %k unit vector
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
e_v = (1./mu).*((V1_mag.^2 - (mu./r1_mag)).*r1 - r1_mag*Vr.*V1); 
e = sqrt(sum(e_v.^2)); 
fprintf('Eccentricity, e =')
disp(e)


%12.) argument of perigee: 
if e_v(3) >= 0
    w = acosd(dot(N_v./N, e_v./e)); 
else 
    w = 360 - acosd(dot(N_v./N, e_v./e));
end
fprintf('Argument of Perigee (degrees), w =')
disp(w)

%13.) true anomaly:
if Vr >= 0
    theta = acosd(dot(e_v./e, r1./r1_mag)); 
else 
    theta = 360 - acosd(dot(e_v./e, r1./r1_mag));
end 
fprintf('True Anomaly (degrees), theta =')
disp(theta)

%14.) semi-major axis
a = (h.^2)./(mu.*(1 - e.^2)); 
%a = r1_mag/(2 - ((r1_mag.*V1_mag.^2)/mu)); 
fprintf('Semi-Major Axis, a =')
disp(a)
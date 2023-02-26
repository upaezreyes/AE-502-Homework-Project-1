clc 
clear variables 

mu = 398600; %(au^3/day^2), sun's gravitational constant 
dt = 3600; %(days), time 

%Initial state position vector at epoch (2027-Jan-01
r1Ix = 7000; %(au), rx vector component
r1Iy = -12124; %(au), ry vector component
r1Iz = 0; %(au), rz vector component 
r0_v = [r1Ix r1Iy r1Iz]; %au, inital state position vector

%Initial state velocity vector
v1Ix = 2.6679; %(au/day), vx vector component
v1Iy = 4.6210; %(au/day), vy vector component
v1Iz = 0; %(au/day), vz vector component 
v0_v = [v1Ix v1Iy v1Iz]; %au/day, inital state velocity vector

%calculate magnitude of the position & velocity vectors 
r0 = sqrt((r1Ix.^2) + (r1Iy.^2) + (r1Iz.^2)); %(au), initial position magnit.
v0 = sqrt(v1Ix.^2 + v1Iy.^2 + v1Iz.^2); %(au/day), initial velocity magnit.

%calculate the radial component velocity
Vr0 = (r1Ix.*v1Ix + r1Iy.*v1Iy + r1Iz.*v1Iz)./r0; %(au/day)

%calculate the reciprocal d of the simimajor axis, a
d = (2./r0) - ((v0.^2)./mu); %(1/au), reciprocal of the a

if d < 0 
    fprintf('The Trajectory is a Hyperbola becasue α is:')
    disp(dt)
elseif d > 0 
    fprintf('The Trajectory is an Ellipse becasue α is:')
    disp(d)
else 
    fprintf('The Trajectory is a Parabola because α is:')
    disp(d)
end

X0 = sqrt(mu).*abs(d).*dt; %starting value of X
% for hw/proj 1: 
X = Kepler_Universal(mu, dt, r0, Vr0, X0, d);
fprintf('The Universal Anomaly is: ')
disp(X)

%======================================================================
% for curtis example (example_3_06) 
X2 = Kepler_Universal(398600,3600,10000,3.0752,115.637, -1/19655); 
%disp(X2)
%=======================================================================

%back to X:
z = d.*X.^2; % z value 
C = stumpC(z); % call stumpC function to get C value
S = stumpS(z); % call stumpS function to get S value

%calculate the position vector (rf) using the f and g funtions
f = 1 - ((X.^2)./r0).*C; % f function
g = dt - (1./sqrt(mu)).*(X.^3).*S; % g function


rf1 = f.*r0_v; %(au), 
rf2 = g.*v0_v; %(au), 
rf_v = rf1 + rf2; % au, final state position vector
rf_v_mag = sqrt(sum(rf_v.^2)); % au, magnitude of rf_v
fprintf('The position vector, r(t), [i j](km) is: ') 
disp(rf_v)

%calculate the velocity vector (vf) using the fdot & gdot functions
fdot = (sqrt(mu)./(rf_v_mag.*r0)).*(d.*(X.^3).*S - X); 
gdot = 1 - ((X.^2)./rf_v_mag).*C; 

vf1 = fdot.*r0_v; 
vf2 = gdot.*v0_v; 
vf_v = vf1 + vf2; % au/day, final state velocity vector
fprintf('The velocity vector, v(t), [i j](km/s) is: ')
disp(vf_v)


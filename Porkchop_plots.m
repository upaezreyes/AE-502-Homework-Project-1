%===================================
clc 
clear variables 

%Calling position & velocity vectors from Earth: 
[~,departure] = LoadData('Earth_Borisov.txt'); 
%tdep = Data1(:,:)
depart_t = departure.calDate(:,1); % time in columnw
depart_v = departure.v(:,1:3); % velocity in column
depart_r = departure.r(:,1:3); % position in column

% v_depart = sqrt(sum(departure_v.^2)); %magnitude of depart velocity
% r_depart = sqrt(sum(departure_r.^2)); %magnitude of depart position

%Calling position & velocity vectors from Oumuamau: 
[err,arrival] = LoadData('Borisov.txt'); 
arriv_t = arrival.calDate(:,1); % time in column
arriv_v = arrival.v(:,1:3); % velocity in column
arriv_r = arrival.r(:,1:3); % position in column

% v_arrival = sqrt(sum(arrival_v.^2)); 
% r_arrival = sqrt(sum(arrival_r.^2));
 
%fileID = fopen('results3.dat','w');
depart_length = length(depart_t); 
arrival_length = length(arriv_t); 
total = depart_length*arrival_length; 

x = depart_length; 
y = arrival_length; 

%====================================================
%Oumuamua: 
x_axis = linspace(0,335,335); 
y_axis = linspace(0,519,519); 
%====================================================

%====================================================
%Borisov: 
% x_axis = linspace(0,1278,1278); 
% y_axis = linspace(0,946,946); 
%=================================================

for j = 1:y   %arrival_t   %213:1:730    %arrival
    for i = 1:x   %departure_t  %1:1:212   %departure
        
        %time of fligth: 
        dt = ((arriv_t(j,1) - depart_t(i,1)))*86400 ; % sec, 

        mu = 1.327e11; % km^3/s^2, gravity constant
        z0 = 20; %initial z value
        %take the cross product btw the deparute & arrival position vectors: 
        r1xr2 = cross(depart_r(i,:),arriv_r(j,:));
        %take the dot product bwt the deparutre & arrival position vectors:
        dot_r1r2 = dot(depart_r(i,:),arriv_r(j,:)); 

        r1 = sqrt(sum(depart_r(i,:).^2)); % magnitude of depart. position vectors
        r2 = sqrt(sum(arriv_r(j,:).^2)); % magnitude of arrival position vectors
        
        %=============================================================
        %for prograde trajectory:
        %=============================================================
        if r1xr2(3) >= 0
            %change of true anomaly angle
            delta_theta = real(acosd(dot_r1r2./(r1.*r2))); 
         else
            delta_theta = 360 - real(acosd(dot_r1r2./(r1.*r2))); 
        end
        %===========================================================

        %====================================================
        % Retrograte Trajectory
        %====================================================
%         if r1xr2(3) >= 0
%             delta_theta = 360 - acosd(dot_r1r2./(r1.*r2));  %change of true anomaly angle
%         else
%             delta_theta =  acosd(dot_r1r2./(r1.*r2));
%         end
        %========================================================

        %constant A:
        A = sind(delta_theta).*sqrt((r1.*r2)./(1-cosd(delta_theta)));
        z = lambert2(mu, dt, z0, A, r1, r2);
        C = stumpC(z);
        S = stumpS(z);

        %calcualte y(y);
        y = r1 + r2 + A.*((z.*S - 1)./(sqrt(C)));

        %calculate the F & G functions
        f = 1 - (y./r1);
        g = A.*sqrt(y./mu);
        gdot = 1 - (y./r2);

        % s/c departure velocity: 
        try
         V1 = real((1./g).*(arriv_r(j,:) - f.*depart_r(i,:))); % V1
        catch
         V1  = [30,
                30,
                30];
        end

        C3a = real(V1 - depart_v(i,:)); % deltaV departure vector 
        C3(i,j) = real(sqrt(sum(C3a.^2))); % magnitude of deltaV departure
        C3b = C3.'; % magnitude of deltaV departure

        % s/c arrival velocity:
        try
            V2 = (1./g).*(gdot.*arriv_r(j,:) - depart_r(i,:)); % V2
        catch
             V2 = [30,
                   30,
                   30]; 
        end

        Vinfa = real(V2 - arriv_v(j,:)); % deltaV arrival vector
        Vinf(i,j) = real(sqrt(sum(Vinfa.^2))); % mag. of deltaV arrival
        Vinfb = Vinf.'; % mag. of deltaV arrival
      
        % Total Delta V: 
        DeltaV(i,j) = C3(i,j) + Vinf(i,j);
        DV = DeltaV.'; 

    end
    
end

% contour line plots: ==============================================
% figure
% contour(x_axis, y_axis, DV,[20:2:50],"ShowText","on")
% xlabel('Departure Day')
% ylabel("arrival Day")
% title('Delta V')
% 
% 
% figure
% contour(x_axis,y_axis,C3b, [0:5:20],"ShowText","on")
% xlabel("Departure Day")
% ylabel("Arrival Day")
% title('V1')

% contour color code plots: =========================================
% total delta V: 
figure
contourf(x_axis, y_axis,DV,[20:2:60])
xlabel('Departure Day')
ylabel('arrival Day')
title('Rendez-Vous, ∆V')
h = colorbar; 
ylabel(h,'deltaV')

% departure delta V: 
figure
contourf(x_axis,y_axis,C3b,[0:2:20])
xlabel("Departure Day")
ylabel("Arrival Day")
title('Fly-by, ∆V')
h = colorbar; 
ylabel(h,'Vinf')

% arrival delta V: 
% figure
% contourf(x_axis,y_axis,Vinfb,[0:2:20])
% xlabel("Departure Day")
% ylabel("Arrival Day")
% title('Fly-by, ∆V_arrival')
% h = colorbar; 
% ylabel(h,'Vinf')
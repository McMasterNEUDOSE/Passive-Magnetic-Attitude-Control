% Current best model
% 1-2-3 Euler angles
tic;
clear all;
close all;
% Importing files
m_bar=0.5; % magnetic moment of bar magnet
vol=1e-7;%volume of hysteresis rods per axis
days=1;
tstep=1; % Time step of propagation data

[ggbody_name,attecfname,B_bodyname,velbody_name,tname]=fileselec(days); %selects file names for simulation
ggbody=cell2mat(struct2cell(load(ggbody_name))); % Import gravity gradient torque in body inertial frame
attecfSTK=cell2mat(struct2cell(load(attecfname))); % Import attitude of body inertial frame to ecf
attecf=[attecfSTK(:,4) attecfSTK(:,1) attecfSTK(:,2) attecfSTK(:,3)]; % Q4 is scalar in STK, but Q1 is scalar part in Matlab so it needs to be converted
B_body_in=(1e-9)*cell2mat(struct2cell(load(B_bodyname))); % Import magnetic field in body inertial frame (convert from nT to T)
velbody=(1e3)*cell2mat(struct2cell(load(velbody_name))); % Import velocity in the body's inertial frame, convert from km/s to m/s
time=cell2mat(struct2cell(load(tname)))'; % Import time of simulation in epoch seconds

orbitperiod=3600+32*60+38.54;%from : https://keisan.casio.com/exec/system/1224665242
orbits=time/orbitperiod;


endind=round(1*length(time))-1; % determines how long to loop through final loop for calculating torques
endtime=endind*tstep; %final time for ODE 45 to compute to

% ----------------------------------------------------------------
% Initial Conditions/ initializing

w0=[5 ;5 ;5]*pi/180;
att0=[13.9;104.1;-71.6]*pi/180;
%att0=[25.9;132.1;-71.6]*pi/180; % initial beta=178.6 degrees
att0=angle2quat(att0(1),att0(2),att0(3),'XYZ');
init=[w0; att0'];
J=[6433703.78 0 0; 0 6484881.04 0;0 0 2901645.34]*(1e-9);
% J=[6433703.78 1176.84 -285846.68; 1176.84 6484881.04 18874.53; -259656.12 18874.53 2901645.34]*(1e-9); % Moment of inertia Matrix

%---------------------------------------------------------------------
% Propagation inital
tspan1=[0 tstep];
options = odeset('AbsTol', 1e-3); % Set absolute and relative error tolerances for ode45
[ti, Xi] = ode23t(@(ti,Xi) eqset(Xi,ggbody(1,:),B_body_in(2,:),B_body_in(1,:),velbody(1,:),m_bar,vol), tspan1, init,options); % Solve ODE
wi=Xi(:,1:3); %this is w in the body frame
Qi=Xi(:,4:7);
Qi=quatnormalize(Qi);
[xi,yi,zi]=quat2angle(Qi,'YXZ');

%---------------------------------------------------------------------
%  Second Initial Conditions/ initializing

w0=wi(end,:)';
att0=Qi(end,:)';
init=[w0; att0];

%---------------------------------------------------------------------
% Propagation main

tspan=[tstep endtime];
% options = odeset('RelTol', 1e-7, 'AbsTol', 1e-7); % Set absolute and relative error tolerances for ode45
% at every timstep ODE 45 Generates, which is DIFFERENT from the timestep
% of the data collection from STK as Matlab determines this based on
% convergence, round the time, and use the find functiont to see the time
% from STK it is  closest to. This returns an indec, and that index is
% where we take the disturbance torque values from
[t2, X2] = ode23t(@(t2,X2) eqset(X2,ggbody(find(time==round(t2)),:),B_body_in(find(time==round(t2)),:),B_body_in(find(time==round(t2))-1,:),velbody(find(time==round(t2)),:),m_bar,vol), tspan, init,options); % Solve ODE

w2=X2(:,1:3); %this is w in the body frame
Q2=X2(:,4:7);
Q2=quatnormalize(Q2); 
w=[wi;w2];
Q=[Qi;Q2];
t=[ti;t2];
[x,y,z]=quat2angle(Q,'XYZ');
att=[x y z];
att=att*180/pi;

%-------------------------------------------------------------------
% Calculating pointing error of satellite
% Dot product of unit vector of magenetic field with unit vector of y-axis
% of satellite in satellite fixed body frame, then solve for cosine of the
% pointing error
B_norm=B_body_in./vecnorm(B_body_in,2,2); % normalize the magentic field vector in the body inertial frame
% Q needs to be reduced, becuase there are more timesteps than would've
% been taken from the STK data. So it'll be taken at the timesteps closest
% to the real timesteps from STK, which is what QWsel does
[Q_samp,W_samp]=QWsel(time(1:(endind+1)),t,Q,tstep,w);
B_body_rot=quatrotate(Q_samp,B_body_in(1:(endind+1),:)); % rotate the magnetic field vector from body inertial to body fixed
B_body=quatrotate(Q_samp,B_norm(1:(endind+1),:)); % rotate the magnetic field unit vector from body inertial to body fixed

% this dot product with the y-axis of the satellite, will just be the y
% column of the magnetic field, becuase in the body fixed frame, the y-axis
% is just (0,1,0), so then we directly just do an arc cosine
B_body_unit=B_body./vecnorm(B_body,2,2);  %confirm it is normalized
% What is occuring above is the magnetic field vector is rotated into the
% body fixed axes frame, and in that frame, the y axes which we care about
% is (0,1,0), and I also normalized the magnetic field above, so a
% dotproduct between the two (both of magnitude 1), would just be the y
% component of the magnetic field in the body fixed frame, then take arc
% cosine to find the angle between teh vectors
beta=acos(B_body_unit(:,2)); % pointing error in radians
beta_deg=beta*180/pi; %convert to degrees
pointacc=10*ones(length(time),1);
req=60*ones(length(time),1);


% --------------------------------------------------------------------
% Torque Calculations
% These just make the arrays to store the torque data
Qprop=Q_samp;
Wprop=W_samp;
L=length(Qprop);
ggrot=quatrotate(Qprop,ggbody);
ggrotmag=vecnorm(ggrot,2,2);

hyst_torque=zeros(L,3);
hyst_mag=zeros(L,1);
barmag_torque=zeros(L,3);
barmag_mag=zeros(L,1);
mag_total=zeros(L,3);
magtotal_norm=zeros(L,1);
gyro_torque=zeros(L,3);
gyro_mag=zeros(L,1);
aero_torque=zeros(L,3);
aero_mag=zeros(L,1);


for i=1:L
    %aerodynamic torque
    aero_torque(i,:)=aero(Qprop(i,:),velbody(i,:));
    aero_mag(i)=norm(aero_torque(i,:));

    % magnetic torques
    if i==1
        bi=quatrotate(Qprop(2,:),B_body_in(2,:));
        b0=quatrotate(Qprop(1,:),B_body_in(1,:));
        barmag_torque(1,:)=barmag_torque_calc(b0,m_bar);
        barmag_mag(1)=norm(barmag_torque(1,:));
        hyst_torque(1,:)=hyst_torque_calc(bi,b0,vol);
        hyst_mag(1)=norm(hyst_torque(1,:));
        mag_total(1,:)=barmag_torque(1,:)+hyst_torque(1,:);
        magtotal_norm(1)=norm(mag_total(1,:));
    else
        bi=quatrotate(Qprop(i,:),B_body_in(i,:));
        b0=quatrotate(Qprop(i-1,:),B_body_in(i-1,:));
        barmag_torque(i,:)=barmag_torque_calc(b0,m_bar);
        barmag_mag(i)=norm(barmag_torque(i,:));
        hyst_torque(i,:)=hyst_torque_calc(bi,b0,vol);
        hyst_mag(i)=norm(hyst_torque(i,:));
        mag_total(i,:)=barmag_torque(i,:)+hyst_torque(i,:);
        magtotal_norm(i)=norm(mag_total(i,:));
    end
    gyro_torque(i,:)=gyro_torque_calc(Wprop(i,:));
    gyro_mag(i)=norm(gyro_torque(i,:));
    
   
end

% total disturbance torques
dist=aero_torque+ggrot;
mag_dist=vecnorm(dist,2,2);
max_dist=max(mag_dist);

Bnorm=vecnorm(B_body_rot,2,2);
Bmin=min(Bnorm);

m_min=15*max_dist/(Bmin*sind(10)); % the minimum bar magnet moment needed for stability
point_error=mean(beta_deg(end-86400:end)) % mean over the last day of orbit


% --------------------------------------------------------------------


% Plotting


figure()
plot(orbits(1:(endind+1)),beta_deg);
hold on
plot(orbits,pointacc);
plot(orbits,req,'r');
title('Pointing error of Bar Magnet axis');
ylabel('Beta (degrees)');
xlabel('Orbits');
movegui('northeast');


figure()
plot(t/orbitperiod,w(:,1)*180/pi);
title('Euler Angle Rates vs Time');
ylabel('Angular Velocity (deg/s)');
xlabel('Orbits');
hold on
plot(t/orbitperiod,w(:,3)*180/pi);
plot(t/orbitperiod,w(:,2)*180/pi);
legend('X','Z','Y');
movegui('southeast');
hold off

figure()
semilogy(orbits,aero_mag,'r');
hold on
semilogy(orbits,ggrotmag,'g');
semilogy(orbits,barmag_mag,'c');
semilogy(orbits,hyst_mag,'y');
semilogy(orbits,gyro_mag,'k');
title('Magnitude of Torques acting on the Satellite');
xlabel('Orbits');
ylabel('Norm of Torque (N m)');
legend({'Aerodynamic Torque','Gravity Gradient Torque','Bar Magnet Torque','Hysteresis Rods Torque','Gyroscopic Torque'},'Location','southwest');
movegui('southwest');



%------------------------------------------------------------------------
% saving data for Liam
% Qecef=quatmultiply(attecf,Q_samp);
% Qsave=[time Qecef(:,2) Qecef(:,3) Qecef(:,4) Qecef(:,1)];
% Qsavename='Q_7day.txt';
% fileID = fopen(Qsavename,'w');
% fprintf(fileID,'%e %e %e %e %e\n',Qsave');
% fclose(fileID);
% todeg=180/pi;
% Qecef=quatmultiply(attecf,Q_samp);
% Qsave=[time Qecef(:,2) Qecef(:,3) Qecef(:,4) Qecef(:,1) todeg*W_samp(:,1) todeg*W_samp(:,2) todeg*W_samp(:,3)];
% QWsavename='QW_7day.txt';
% fileID = fopen(QWsavename,'w');
% fprintf(fileID,'%e %e %e %e %e\n',Qsave');
% fclose(fileID);

toc;

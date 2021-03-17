close all;
clear all;
tic;

% Simulation parameters
days=14;
tstep=20;
rho_year = 2022; % Year for air density (2022,2023,2020,2014)
w0=[2;2;2]*pi/180; % nanoracks worst case scenario

% Physical Design Parameters
m = 0.32; % Possible from combination of magnets from KJ
hyst_l = 0.05; % Orignal Length of hysteresis rods in m
hyst_d = 0.001; % Diameter of hysteresis rods in m
nrods = 2; % rods per axis
vol=nrods*0.25*pi*hyst_l*hyst_d^2; % Volume of hysteresis rods

mres=[-0.0046;-0.0046;0.000575]; %residual magnetic moment in A*m^2

% --------------------------------------------------------

% Loading in Data
[data_fname,rho_fname] = file_selec(days,tstep,rho_year);
data=cell2mat(struct2cell(load(data_fname))); 
time=data(:,1); pos=data(:,2:4); vel=data(:,5:7); B_in=data(:,8:10); Bdot=data(:,11:13); 
sun=data(:,14:16).*data(:,17); 
dens=cell2mat(struct2cell(load(rho_fname)));
% --------------------------------------------------------

     
% Define constants
period=3600+32*60+39; % orbital period from https://keisan.casio.com/exec/system/1224665242
mu0=4*pi*1e-7;
Bs=0.3;
Br=6.0618e-4; % from gerhardt dissertation
Hc=0.3381; %from gerhardt dissertation
BsFac=0.5*pi/Bs;
k=(1/Hc)*tan(0.5*pi*Br/Bs);
rad2deg=180/pi;
tspan=[0 days*86400];
% --------------------------------------------------------

% Initial values for ODE solver
% att_init=[15.5;-1.2;-30]*pi/180; %initial beta  0.0555 deg
% att_init=[45;35;-65]*pi/180; %initial beta 45.5022 deg
% att_init=[30;88.8;-104]*pi/180; %initial beta  90 deg
% att_init=[-150;58.8;-102.6]*pi/180; %initial beta 121.3192 deg

% att_init=[15;180;-150]*pi/180; % NON PHYSICAL EULER ANGLES initial beta 178.6576 deg
att_init=[-165.7;1.3;35]*pi/180; %initial beta   178.7483 deg
% att_init=[15;177;-150]*pi/180; %  initial beta 178.1416 deg
% att_init=[15;179.5;-150]*pi/180; %  initial beta   179.0905 deg


att0=angle2quat(att_init(1),att_init(2),att_init(3),'XYZ');

B0=quatrotate(att0,B_in(1,:));
Bdot0 = cross(-w0, B0) + quatrotate(att0, Bdot(1,:));
Bh0 = (2*Bs/pi)*atan(k*((B0/mu0)-sign(Bdot0/mu0)*Hc))';
% Bh0 =[0;0;0];


init=[w0; att0';Bh0];

% ODE solver options
% options=odeset('RelTol',1e-7,'AbsTol',1e-10,'Stats','on');
% options=odeset('RelTol',1e-10,'AbsTol',1e-10,'Stats','on');
% options=odeset('RelTol',1e-11,'AbsTol',1e-11,'Stats','on');
options=odeset('RelTol',1e-12,'AbsTol',1e-12,'Stats','on');
% options=odeset('RelTol',1e-13,'AbsTol',1e-13,'Stats','on');
% options=odeset('RelTol',1e-7,'AbsTol',1e-7,'Stats','on');
% --------------------------------------------------------

% Run simulations
[T,X]=ode113quat(@(t,x) eqset(x,t,time,B_in,Bdot,pos,sun,vel,dens,m,hyst_l,hyst_d,nrods,mres),tspan,init,options);
% [T,X]=ode420(@(t,x) eqset(x,t,time,B_in,Bdot,pos,sun,vel,dens,m,hyst_l,hyst_d,nrods,mres),tspan,init,options);
% [T,X]=odeRK4(@(t,x) eqset(x,t,time,B_in,Bdot,pos,sun,vel,dens,m,hyst_l,hyst_d,nrods,mres),tspan,init,0.1);


% Output data from simulation
w=X(:,1:3); % Angular velocities
q=X(:,4:7); % Attitude quaternions
BH=X(:,8:10); % Hysteresis magnetic field
t_days=T/86400; % Time in days
wdeg=(180/pi)*w; % Convert angular rates to degrees/second

% Saving data
if days<20
    X_save = downsample(X,20);
    X_save_name = strcat('X_save_',num2str(days),'days.mat');
    q_save = downsample(q,20);
    q_save_name = strcat('q_save_',num2str(days),'days.mat');
    w_save = downsample(wdeg,20);
    w_save_name = strcat('w_save_',num2str(days),'days.mat');
    t_save = downsample(T,20);
    t_save_name = strcat('t_save_',num2str(days),'days.mat');
    save(X_save_name,'X_save');
    save(q_save_name,'q_save');
    save(t_save_name,'t_save');
    save(w_save_name,'w_save');
else
    X_save_name = strcat('X_save_',num2str(days),'days.mat');
    X_save = downsample(X,100);
    q_save = downsample(q,100);
    q_save_name = strcat('q_save_',num2str(days),'days.mat');
    w_save = downsample(wdeg,100);
    w_save_name = strcat('w_save_',num2str(days),'days.mat');
    t_save = downsample(T,100);
    t_save_name = strcat('t_save_',num2str(days),'days.mat');
    save(X_save_name,'X_save');
    save(q_save_name,'q_save');
    save(t_save_name,'t_save');
    save(w_save_name,'w_save');
end
% --------------------------------------------------------

% PLOTTING SECTION
% Plot angular rates
w_figname = strcat('w_',num2str(days),'day.png');
w_fig = figure();
plot(t_days,wdeg(:,1),'g','LineWidth',1);
hold on
plot(t_days,wdeg(:,2),'r','LineWidth',1);
plot(t_days,wdeg(:,3),'b','LineWidth',1);
legend('wx','wy','wz','LineWidth',1);
title('Satellite angular rates');
xlabel('Time (days)');
ylabel('Angular Velocity (deg/s)');
set(gca,'fontname','arial')  % Set it to arial
set(gca,'fontsize',13)  % Set it to arial
% ylim([-10 10]);
ax = gca;
ax.LineWidth = 1;
ay = gca;
ay.LineWidth = 1;
fig = gcf;
fig.Color = [1 1 1];
box off
grid 
movegui('northwest');
saveas(w_fig,w_figname);


goal=10*ones(length(T),1); % Goal pointing error
Beta    = calcBeta(T, quatnormalize(q), interp1(time,B_in,T)); % Calculate pointing error angle
% Calculate settling time and steady state error
if days>6
    wind=length(T)/10;
    runav=movmean(Beta,wind);
    under=runav<10;
    undervals=find(under==1);
    if length(undervals)>1
        set_ind=undervals(1);
        set_time=T(set_ind)/86400 % settling time in days
        ss_error=mean(Beta(set_ind:end)) % average steady state error
    end
end

beta_figname = strcat('beta_',num2str(days),'day.png');
beta_fig = figure();
plot(t_days,Beta,'r','LineWidth',1);
hold on
plot(t_days,goal,'b','LineWidth',1);
title('Pointing error with respect to geomagnetic field lines');
xlabel('Time (days)');
ylabel('Pointing Error (deg)');
set(gca,'fontname','arial')  % Set it to arial
set(gca,'fontsize',13)  % Set it to arial
ax = gca;
ax.LineWidth = 1;
ay = gca;
ay.LineWidth = 1;
fig = gcf;
fig.Color = [1 1 1];
box off
grid 
movegui('northeast');
saveas(beta_fig,beta_figname);


Nad    = calcNad(T, quatnormalize(q), interp1(time,pos,T));

figure()
plot(t_days,Nad,'b','LineWidth',1);
title('Long axis angle to nadir');
xlabel('Time (days)');
ylabel('Nadir angle (deg)');
set(gca,'fontname','arial')  % Set it to arial
set(gca,'fontsize',13)  % Set it to arial
ax = gca;
ax.LineWidth = 1;
ay = gca;
ay.LineWidth = 1;
fig = gcf;
fig.Color = [1 1 1];
box off
grid on
movegui('north');

%Energy analyis
% E_k=dot((w.^2),([0.0065 0.0065 0.0029].*ones(length(T),3)),2);
% E_p=-dot(([0 m 0].*ones(length(T),3)),quatrotate(q,interp1(time,B_in,T)),2);
% E_T=E_k+E_p;
% figure()
% plot(t_days,E_k,t_days,E_p,t_days,E_T);
% title('Energy of Satellite');
% xlabel('Time (days)');
% ylabel('Energy (J)');
% legend('Kinetic','Potential','Total');

if days>10
    BH_downs = downsample(BH,10000);
    t_downs = downsample(T,10000);
    w_downs = downsample(w,10000);
    q_downs = downsample(q,10000);
else
    BH_downs = downsample(BH,10);
    t_downs = downsample(T,10);
    w_downs = downsample(w,10);
    q_downs = downsample(q,10);
end



[bar_torque,hyst_torque,gyro_torque,gg_torque,aero_torque,sun_torque,eddy_torque,res_torque]...
    =plotT(t_downs,time,q_downs,w_downs,interp1(time,B_in,t_downs),BH_downs,interp1(time,pos,t_downs),interp1(time,vel,t_downs),interp1(dens(:,1),dens(:,2),t_downs),interp1(time,sun,t_downs),m,hyst_l,hyst_d,nrods);

t_days_downs = (1/86400)*t_downs;

torque_mat = [t_days_downs, vecnorm(bar_torque,2,2),vecnorm(hyst_torque,2,2),vecnorm(gyro_torque,2,2),vecnorm(gg_torque,2,2),vecnorm(aero_torque,2,2)...
   ,vecnorm(sun_torque,2,2),vecnorm(eddy_torque,2,2),vecnorm(res_torque,2,2)];

torque_save_name = strcat('torque_save_',num2str(days),'days.mat');
save(torque_save_name,'torque_mat');

z_torque_mat = [t_days_downs,abs(hyst_torque(:,3)),abs(gg_torque(:,3)),abs(aero_torque(:,3))...
   ,abs(sun_torque(:,3)),abs(eddy_torque(:,3)),abs(res_torque(:,3))];

z_torque_save_name = strcat('z_torque_save_',num2str(days),'days.mat');
save(z_torque_save_name,'z_torque_mat');


ztorque_figname = strcat('ztorque_',num2str(days),'day.png');
ztorque_fig = figure();
semilogy(t_days_downs,abs(hyst_torque(:,3)),'r','LineWidth',1);
hold on
semilogy(t_days_downs,abs(res_torque(:,3)),'g','LineWidth',1);
semilogy(t_days_downs,abs(gg_torque(:,3)),'y','LineWidth',1);
semilogy(t_days_downs,abs(aero_torque(:,3)),'b','LineWidth',1);
semilogy(t_days_downs,abs(sun_torque(:,3)),'c','LineWidth',1);
semilogy(t_days_downs,abs(eddy_torque(:,3)),'k','LineWidth',1);
title('Torques acting about Satellite Z-axis');
xlabel('Time (days)');
ylabel('Torque (N m)');
legend({'||T_H||','||T_R||','||T_{GG}||','||T_D||','||T_{SP}||','||T_{EC}||'},'Location','eastoutside');
set(gca,'fontname','arial')  % Set it to arial
set(gca,'fontsize',13)  % Set it to arial
ax = gca;
ax.LineWidth = 1;
ay = gca;
ay.LineWidth = 1;
fig = gcf;
fig.Color = [1 1 1];
box off
grid on
movegui('southwest');
% saveas(ztorque_fig,ztorque_figname);

torque_figname = strcat('torque_',num2str(days),'day.png');
torque_fig = figure();
semilogy(t_days_downs,vecnorm(bar_torque,2,2),'r','LineWidth',1);
hold on
semilogy(t_days_downs,vecnorm(hyst_torque,2,2),'g','LineWidth',1);
semilogy(t_days_downs,vecnorm(res_torque,2,2),'b','LineWidth',1);
semilogy(t_days_downs,vecnorm(gg_torque,2,2),'m','LineWidth',1);
semilogy(t_days_downs,vecnorm(aero_torque,2,2),'c','LineWidth',1);
semilogy(t_days_downs,vecnorm(sun_torque,2,2),'k','LineWidth',1);
semilogy(t_days_downs,vecnorm(eddy_torque,2,2),'Color',[0.4660 0.6740 0.1880],'LineWidth',1);
title('Norm of External Torques');
xlabel('Time (days)');
ylabel('Torque (N m)');
legend({'||T_B||','||T_H||','||T_R||','||T_{GG}||','||T_D||','||T_{SP}||','||T_{EC}||'},'Location','eastoutside');
set(gca,'fontname','arial')  % Set it to arial
set(gca,'fontsize',13)  % Set it to arial
ax = gca;
ax.LineWidth = 1;
ay = gca;
ay.LineWidth = 1;
fig = gcf;
fig.Color = [1 1 1];
box off
grid on
movegui('southwest');
saveas(torque_fig,torque_figname);


toc;
% % 



% L=length(T);
% H = zeros(L,3);
% Bdot_rot=zeros(L,3);
% for i=1:L
%     Bnow  = quatrotate(q(i,:),lininterp1(time,B_in(:,1:3),T(i)));
%     Bdot_body=quatrotate(q(i,:),interp1(time,Bdot(:,1:3),T(i)));
%     Bdot_rot(i,:)= cross(-w(i,:), Bnow) + Bdot_body;
%     H(i,:)  = Bnow/mu0;
% end
% BinvTan=atan(k*(H-sign(Bdot_rot)*Hc))/BsFac;
% 
% figure()
% plot(t_days,BinvTan(:,1),t_days,BH(:,1))
% legend('invtan','flatley');
% hyst_titlex=strcat('Hysteresis Loop in X');
% hyst_titlez=strcat('Hysteresis Loop in Y');
% 
% figure()
% plot(H(:,1),BH(:,1),'r');
% hold on
% plot(H(:,1),BinvTan(:,1),'b');
% title(hyst_titlex);
% xlabel('Magnetizing Field H (A/m)');
% ylabel('Magnetic Flux Density (T)');
% legend('Flatley','Invtan');
% movegui('southeast');
% figure()
% plot(H(:,2),BH(:,2),'r');
% hold on
% plot(H(:,2),BinvTan(:,2),'b');
% title(hyst_titlez);
% xlabel('Magnetizing Field H (A/m)');
% ylabel('Magnetic Flux Density (T)');
% legend('Flatley','Invtan');
% movegui('south');

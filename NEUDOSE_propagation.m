close all;
tic;

% Simulation parameters
days=365;
tstep=20;
w0=[5;5;5]*pi/180; % nanoracks worst case scenario
% w0=[0.17 ;-0.97 ;2.93]*pi/180; %CSSWE initial velocity
% w0=[0;0;0]; % best case scenario

% Physical Design Parameters
m=0.32; % Magnetic moment of bar magnet in A*m^2
hyst_l = 0.05; % Length of hysteresis rods in m
hyst_d = 0.001; % Diameter of hysteresis rods in m
nrods = 2; % rods per axis
vol=nrods*0.25*pi*hyst_l*hyst_d^2; % Volume of hysteresis rods
mres=[-0.0039;-0.0055;0.0004]; %residual magnetic moment in A*m^2
% --------------------------------------------------------

% Loading in Data
str_days=num2str(days); str_tstep=num2str(tstep); strm=num2str(m); strvol=num2str(vol); 
path='data_';
fname=strcat(path,str_days,'days_',str_tstep,'s.mat');
data=cell2mat(struct2cell(load(fname))); 
time=data(:,1); pos=data(:,2:4); vel=data(:,5:7); B_in=data(:,8:10); Bdot=data(:,11:13); 
sun=data(:,14:16); 
if days<= 28
    dens=cell2mat(struct2cell(load('neudoserho60s_28day.mat')));

elseif days <=60
     dens=cell2mat(struct2cell(load('neudoserho60s_60day.mat')));
elseif days <=120
     dens=cell2mat(struct2cell(load('neudoserho120s_120day.mat')));
elseif days <=180
     dens=cell2mat(struct2cell(load('neudoserho180s_180day.mat')));
elseif days<= 365
    dens=cell2mat(struct2cell(load('neudoserho120s_365day.mat')));
end
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
att_init=[30;90;-104]*pi/180; %initial beta 178.6952 deg
att0=angle2quat(att_init(1),att_init(2),att_init(3),'XYZ');
B0=quatrotate(att0,B_in(1,:));
Bdot0 = cross(-w0, B0) + quatrotate(att0, Bdot(1,:));
Bh0 = (2*Bs/pi)*atan(k*((B0/mu0)-sign(Bdot0/mu0)*Hc))';
init=[w0; att0';Bh0];

% ODE solver options
% options=odeset('RelTol',1e-7,'AbsTol',1e-10,'Stats','on');
% options=odeset('RelTol',1e-7,'AbsTol',1e-9,'Stats','on');
options=odeset('RelTol',1e-7,'AbsTol',1e-7,'Stats','on');
% --------------------------------------------------------

% Run simulations
% [T,X]=ode420(@(t,x) eqset(x,t,time,B_in,Bdot,pos,sun,vel,dens,m,hyst_l,hyst_d,nrods,mres),tspan,init,options);
[T,X]=ode113quat(@(t,x) eqset(x,t,time,B_in,Bdot,pos,sun,vel,dens,m,hyst_l,hyst_d,nrods,mres),tspan,init,options);

% Output data from simulation
w=X(:,1:3); % Angular velocities
q=X(:,4:7); % Attitude quaternions
BH=X(:,8:10); % Hysteresis magnetic field
% Downsample data
X = downsample(X,200);
w=X(:,1:3);
q=X(:,4:7);
BH=X(:,8:10);
T = downsample(T,200);
t_days=T/86400; % Time in days
wdeg=(180/pi)*w; % Convert angular rates to degrees/second

% Saving data
X_save_name = strcat('X_save_',num2str(days),'days.mat');
q_save = downsample(q,100);
q_save_name = strcat('q_save_',num2str(days),'days.mat');
w_save = downsample(wdeg,100);
w_save_name = strcat('w_save_',num2str(days),'days.mat');
t_save = downsample(X,100);
t_save_name = strcat('t_save_',num2str(days),'days.mat');
save(X_save_name,'X');
save(q_save_name,'q_save');
save(t_save_name,'t_save');
save(w_save_name,'w_save');

% --------------------------------------------------------

% PLOTTING SECTION
% Plot angular rates
figure()
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
ylim([-10 10]);
ax = gca;
ax.LineWidth = 1;
ay = gca;
ay.LineWidth = 1;
fig = gcf;
fig.Color = [1 1 1];
box off
grid on
movegui('northwest');


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

figure()
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
fig = gcf;
fig.Color = [1 1 1];
box off
grid on
movegui('northeast');

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

toc;
% % 
[bar_torque,hyst_torque,gyro_torque,gg_torque,aero_torque,sun_torque,eddy_torque,res_torque]...
    =plotT(T,time,q,w,interp1(time,B_in,T),BH,interp1(time,pos,T),interp1(time,vel,T),interp1(dens(:,1),dens(:,2),T),interp1(time,sun,T),m,hyst_l,hyst_d,nrods);

% [bar_torque,hyst_torque,gyro_torque,gg_torque,aero_torque,sun_torque,eddy_torque,res_torque]...
%     =plotT(T,time,q,w,interp1(time,B_in,T),BH,interp1(time,pos,T),interp1(time,vel,T),interp1(dens(:,1),dens(:,2),T),interp1(time,sun,T),m,volx,volz);

figure()
semilogy(t_days,bar_torque,'.','MarkerFaceColor','b');
hold on
semilogy(t_days,hyst_torque,'.','MarkerEdgeColor','r','MarkerFaceColor','r');
semilogy(t_days,gyro_torque,'.','MarkerEdgeColor','m','MarkerFaceColor','m');
semilogy(t_days,res_torque,'.','MarkerEdgeColor','g','MarkerFaceColor','g');
semilogy(t_days,gg_torque,'.','MarkerEdgeColor','c','MarkerFaceColor','c');
semilogy(t_days,aero_torque,'.','MarkerEdgeColor',[0.6350 0.0780 0.1840],'MarkerFaceColor',[0.6350 0.0780 0.1840]);
semilogy(t_days,sun_torque,'.','MarkerEdgeColor','k','MarkerFaceColor','k');
semilogy(t_days,eddy_torque,'.','MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerFaceColor',[0.4940 0.1840 0.5560]);
title('Magnitude of Torques acting on the Satellite');
xlabel('Time (days)');
ylabel('Norm of Torque (N m)');
legend({'||L_B||','||L_H||','||L_{GY}||','||L_R||','||L_{GG}||','||L_D||','||L_{SP}||','||L_{EC}||'},'Location','eastoutside');
set(gca,'fontname','arial')  % Set it to arial
set(gca,'fontsize',13)  % Set it to arial
ax = gca;
ax.LineWidth = 1;
ay = gca;
ay.LineWidth = 1;
fig = gcf;
fig = gcf;
fig.Color = [1 1 1];
box off
grid on
movegui('southwest');



if T(end)>80000
%orbit 15 plotting
    tstart=14/16;
    tend=15/16;
    figure()
    semilogy(t_days,bar_torque,'.','MarkerFaceColor','b');
    hold on
    semilogy(t_days,hyst_torque,'.','MarkerEdgeColor','r','MarkerFaceColor','r');
    semilogy(t_days,gyro_torque,'.','MarkerEdgeColor','m','MarkerFaceColor','m');
    semilogy(t_days,res_torque,'.','MarkerEdgeColor','g','MarkerFaceColor','g');
    semilogy(t_days,gg_torque,'.','MarkerEdgeColor','c','MarkerFaceColor','c');
    semilogy(t_days,aero_torque,'.','MarkerEdgeColor',[0.6350 0.0780 0.1840],'MarkerFaceColor',[0.6350 0.0780 0.1840]);
    semilogy(t_days,sun_torque,'.','MarkerEdgeColor','k','MarkerFaceColor','k');
    semilogy(t_days,eddy_torque,'.','MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerFaceColor',[0.4940 0.1840 0.5560]);
%     title('Magnitude of Torques acting on the Satellite in orbit 15');
    xlim([tstart tend]);
    xlabel('Time (days)');
    ylabel('Norm of Torque (N m)');
    legend({'||L_B||','||L_H||','||L_{GY}||','||L_R||','||L_{GG}||','||L_D||','||L_{SP}||','||L_{EC}||'},'Location','eastoutside');
    set(gca,'fontname','arial')  % Set it to arial
    set(gca,'fontsize',13)  % Set it to arial
    ax = gca;
    ax.LineWidth = 1;
    ay = gca;
    ay.LineWidth = 1;
    fig = gcf;
    fig = gcf;
    fig.Color = [1 1 1];
    box off
    grid on
    movegui('south');

   

end


if days>6
    % Orbit 105
    endt=(105*period)/86400;
    tstart=endt-1/16;
    orbitnumber=round(T(end)/period);
    str_orbit=num2str(orbitnumber);
    figure()
    semilogy(t_days,bar_torque,'.','MarkerFaceColor','b');
    hold on
    semilogy(t_days,hyst_torque,'.','MarkerEdgeColor','r','MarkerFaceColor','r');
    semilogy(t_days,gyro_torque,'.','MarkerEdgeColor','m','MarkerFaceColor','m');
    semilogy(t_days,res_torque,'.','MarkerEdgeColor','g','MarkerFaceColor','g');
    semilogy(t_days,gg_torque,'.','MarkerEdgeColor','c','MarkerFaceColor','c');
    semilogy(t_days,aero_torque,'.','MarkerEdgeColor',[0.6350 0.0780 0.1840],'MarkerFaceColor',[0.6350 0.0780 0.1840]);
    semilogy(t_days,sun_torque,'.','MarkerEdgeColor','k','MarkerFaceColor','k');
    semilogy(t_days,eddy_torque,'.','MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerFaceColor',[0.4940 0.1840 0.5560]);
%     title('Magnitude of Torques during orbit 105');
    xlim([tstart endt]);
    xlabel('Time (days)');
    ylabel('Norm of Torque (N m)');
    legend({'||L_B||','||L_H||','||L_{GY}||','||L_R||','||L_{GG}||','||L_D||','||L_{SP}||','||L_{EC}||'},'Location','eastoutside');
    set(gca,'fontname','arial')  % Set it to arial
    set(gca,'fontsize',13)  % Set it to arial
    ax = gca;
    ax.LineWidth = 1;
    ay = gca;
    ay.LineWidth = 1;
    fig = gcf;
    fig = gcf;
    fig.Color = [1 1 1];
    box off
    grid on
    movegui('southeast');

   
end
% 


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

close all;
tic;

days=10;
tstep=20;
m=0.45;
vol=5e-7;
str_days=num2str(days); str_tstep=num2str(tstep); strm=num2str(m); strvol=num2str(vol);
path='data_';
fname=strcat(path,str_days,'days_',str_tstep,'s.mat');
data=cell2mat(struct2cell(load(fname))); 
time=data(:,1); pos=data(:,2:4); vel=data(:,5:7); B_in=data(:,8:10); Bdot=data(:,11:13); 
sun=data(:,14:16); 
dens=cell2mat(struct2cell(load('neudoserho60s_28day.mat')));

period=3600+32*60+39; % orbital period from https://keisan.casio.com/exec/system/1224665242
mu0=4*pi*1e-7;
Bs=0.3;
Br=6.0618e-4; % from gerhardt dissertation
Hc=0.3381; %from gerhardt dissertation
BsFac=0.5*pi/Bs;
k=(1/Hc)*tan(0.5*pi*Br/Bs);


rad2deg=180/pi;
tspan=[0 days*86400];
w0=[5;5;5]*pi/180; 
att_init=[13.9;104.1;-71.6]*pi/180; 
att0=angle2quat(att_init(1),att_init(2),att_init(3),'XYZ');
B0=quatrotate(att0,B_in(1,:));
Bdot0 = cross(-w0, B0) + quatrotate(att0, Bdot(1,:));
Bh0 = (2*Bs/pi)*atan(k*((B0/mu0)-sign(Bdot0/mu0)*Hc))';
init=[w0; att0';Bh0];


options=odeset('RelTol',1e-7,'AbsTol',1e-10,'Stats','on');
% options=odeset('RelTol',1e-7,'AbsTol',1e-9,'Stats','on');

% [T,X]=ode420(@(t,x) eqset(x,t,time,B_in,Bdot,pos,sun,vel,dens),tspan,init,options);
[T,X]=ode113quat(@(t,x) eqset(x,t,time,B_in,Bdot,pos,sun,vel,dens,m,vol),tspan,init,options);


w=X(:,1:3);
q=X(:,4:7);
BH=X(:,8:10);
t_days=T/86400;
% 
wdeg=(180/pi)*w;

figure()
plot(t_days,wdeg(:,1),t_days,wdeg(:,3),t_days,wdeg(:,2));
legend('wx','wz','wy');
wtitle=strcat('Angular Velocities with m=',strm,' A*m^2, vol=',strvol,' m^3');
title(wtitle);
xlabel('Time (days)');
ylabel('Angular Velocity (deg/s)');
grid on
movegui('northwest');

goal=10*ones(length(T),1);
Btitle=strcat('Pointing error with m=',strm,' A*m^2, vol=',strvol,' m^3');
Beta    = calcBeta(T, quatnormalize(q), interp1(time,B_in,T));
figure()
plot(t_days,Beta);
hold on
plot(t_days,goal,'r');
title(Btitle);
xlabel('Time (days)');
ylabel('Pointing Error (deg)');
grid on
movegui('northeast');


toc;
% % 
[bar_torque,hyst_torque,gyro_torque,gg_torque,aero_torque,sun_torque,eddy_torque,res_torque]...
    =plotT(T,time,q,w,interp1(time,B_in,T),BH,interp1(time,pos,T),interp1(time,vel,T),interp1(dens(:,1),dens(:,2),T),interp1(time,sun,T),m,vol);
% 

figure()
semilogy(t_days,bar_torque,'.');
hold on
semilogy(t_days,hyst_torque,'.');
semilogy(t_days,gyro_torque,'.');
semilogy(t_days,res_torque,'.');
semilogy(t_days,gg_torque,'.');
semilogy(t_days,aero_torque,'.');
semilogy(t_days,sun_torque,'.');
semilogy(t_days,eddy_torque,'.');
title('Magnitude of Torques acting on the Satellite');
xlabel('Time (days)');
ylabel('Norm of Torque (N m)');
legend({'||L_B||','||L_H||','||L_{GY}||','||L_R||','||L_{GG}||','||L_D||','||L_{SP}||','||L_{EC}||'},'Location','eastoutside');
movegui('southwest');
fig = gcf;
fig.Color = [1 1 1];
grid on



if T(end)>80000
%orbit 15 plotting
    tstart=1*15/16;
    tend=1;
    figure()
    semilogy(t_days,bar_torque,'.');
    hold on
    semilogy(t_days,hyst_torque,'.');
    semilogy(t_days,gyro_torque,'.');
    semilogy(t_days,res_torque,'.');
    semilogy(t_days,gg_torque,'.');
    semilogy(t_days,aero_torque,'.');
    semilogy(t_days,sun_torque,'.');
    semilogy(t_days,eddy_torque,'.');
    title('Magnitude of Torques acting on the Satellite in orbit 15');
    xlim([tstart tend]);
    xlabel('Time (days)');
    ylabel('Norm of Torque (N m)');
    legend({'||L_B||','||L_H||','||L_{GY}||','||L_R||','||L_{GG}||','||L_D||','||L_{SP}||','||L_{EC}||'},'Location','eastoutside');
    movegui('south');
    fig = gcf;
    fig.Color = [1 1 1];
    grid on
   

end


if days>6
    % Orbit 105
    endt=(105*period)/86400;
    tstart=endt-1/16;
    orbitnumber=round(T(end)/period);
    str_orbit=num2str(orbitnumber);
    figure()
    semilogy(t_days,bar_torque,'.');
    hold on
    semilogy(t_days,hyst_torque,'.');
    semilogy(t_days,gyro_torque,'.');
    semilogy(t_days,res_torque,'.');
    semilogy(t_days,gg_torque,'.');
    semilogy(t_days,aero_torque,'.');
    semilogy(t_days,sun_torque,'.');
    semilogy(t_days,eddy_torque,'.');
    title('Magnitude of Torques during orbit 105');
    xlim([tstart endt]);
    xlabel('Time (days)');
    ylabel('Norm of Torque (N m)');
    legend({'||L_B||','||L_H||','||L_{GY}||','||L_R||','||L_{GG}||','||L_D||','||L_{SP}||','||L_{EC}||'},'Location','eastoutside');
    movegui('southeast');
    fig = gcf;
    fig.Color = [1 1 1];
    grid on
   
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

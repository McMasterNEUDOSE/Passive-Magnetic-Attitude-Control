%clear all;
close all;
tic;
m=0.45;
vol=1e-7;
volx=2e-7;
volz=1e-7;
days=14;
tstep=20;
str_days=num2str(days); str_tstep=num2str(tstep); str_m=num2str(m); str_vol=num2str(vol);
path='C:\Users\NEUDOSE-2\Documents\Zan Simulations\ECI Propagation\ECI Data\data_';
fname=strcat(path,str_days,'days_',str_tstep,'s.mat');
data=cell2mat(struct2cell(load(fname))); 
time=data(:,1); pos=data(:,2:4); vel=data(:,5:7); B_in=data(:,8:10); Bdot=data(:,11:13); 
sun=data(:,14:16); int=data(:,17); % solar intensity (is there sunlight or not)
sunlight=sun.*int; %the unit distance vector from earth to sun, when sun is shining on the satellite

mu0=4*pi*1e-7;
Bs=0.3;
Br=6.0618e-4; % from gerhardt dissertation
Hc=0.3381; %from gerhardt dissertation
BsFac=0.5*pi/Bs;
k=(1/Hc)*tan(0.5*pi*Br/Bs);


rad2deg=180/pi;
tspan=[0 days*86400];
w0=[5 ;5 ;5]*pi/180; 
att_init=[25.9;132.1;-71.6]*pi/180; % initial beta=178.6 degrees
att0=angle2quat(att_init(1),att_init(2),att_init(3),'XYZ');
B0=quatrotate(att0,B_in(1,:));
Bdot0 = cross(-w0, B0) + quatrotate(att0, Bdot(1,:));
Bh0 = (2*Bs/pi)*atan(k*((B0/mu0)-sign(Bdot0/mu0)*Hc))';
options=odeset('MaxStep',0.01);
init=[w0; att0';Bh0];

[T,X]=ode420(@(t,x) eqset(x,t,time,B_in,m,Bdot,pos,sunlight,vel,volx,volz),tspan,init);%,options);
w=X(:,1:3);
wmag=vecnorm(w,2,2);
q=X(:,4:7);
%qnorm=vecnorm(q,2,2);
BH=X(:,8:10);

% plotting
wtitle=strcat('Angular Velocity with m=',str_m,' A m^2',' hyst vol=',str_vol,' m^3');
t_days=T/86400;
figure()
plot(t_days,rad2deg*w(:,1));
hold on
plot(t_days,rad2deg*w(:,2));
plot(t_days,rad2deg*w(:,3));
plot(t_days,rad2deg*wmag);
title(wtitle)
xlabel('Time (days)');
ylabel('Angular Velocity (deg/s)');
legend('x','y','z','|w|');
movegui('northwest');

Btitle=strcat('Pointing error with m=',str_m,' A m^2',', hyst vol=',str_vol,' m^3');
Beta    = calcBeta(T, quatnormalize(q), interp1(time,B_in,T));
figure()
plot(t_days,Beta);
title(Btitle);
xlabel('Time (days)');
ylabel('Pointing Error (deg)');
movegui('northeast');

[bar_torque,hyst_torque,gyro_torque,gg_torque,res_torque,sun_torque,eddy_torque,aero_torque]...
    =plotT(T,time,q,w,B_in,BH,pos,sunlight,vel,m,volx,volz);
toc;
msize=3;
period=5600;
figure()
semilogy(t_days,bar_torque,'linestyle','none','marker','.','MarkerSize',msize);
hold on
semilogy(t_days,hyst_torque,'linestyle','none','marker','.','MarkerSize',msize);
semilogy(t_days,gyro_torque,'linestyle','none','marker','.','MarkerSize',msize);
semilogy(t_days,gg_torque,'linestyle','none','marker','.','MarkerSize',msize);
semilogy(t_days,res_torque,'linestyle','none','marker','.','MarkerSize',msize);
semilogy(t_days,sun_torque,'linestyle','none','marker','.','MarkerSize',msize);
semilogy(t_days,eddy_torque,'linestyle','none','marker','.','MarkerSize',msize);
semilogy(t_days,aero_torque,'linestyle','none','marker','.','MarkerSize',msize);
title('Magnitude of Torques acting on the Satellite');
xlabel('Time (days)');
ylabel('Norm of Torque (N m)');
legend({'||L_B||','||L_H||','||L_{GY}||','||L_{GG}||','||L_{R}||','||L_{S}||','||L_{EC}||','||L_D||'},'Location','eastoutside');
movegui('southwest');
fig = gcf;
fig.Color = [1 1 1];
grid on


if T(end)>80000
%orbit 16 plotting
    tstart=1*15/16;
    tend=1;
    figure()
    semilogy(t_days,bar_torque,'linestyle','none','marker','.','MarkerSize',msize);
    hold on
    semilogy(t_days,hyst_torque,'linestyle','none','marker','.','MarkerSize',msize);
    semilogy(t_days,gyro_torque,'linestyle','none','marker','.','MarkerSize',msize);
    semilogy(t_days,gg_torque,'linestyle','none','marker','.','MarkerSize',msize);
    semilogy(t_days,res_torque,'linestyle','none','marker','.','MarkerSize',msize);
    semilogy(t_days,sun_torque,'linestyle','none','marker','.','MarkerSize',msize);
    semilogy(t_days,eddy_torque,'linestyle','none','marker','.','MarkerSize',msize);
    semilogy(t_days,aero_torque,'linestyle','none','marker','.','MarkerSize',msize);
    title('Magnitude of Torques acting on the Satellite in orbit 16');
    xlabel('Time (days)');
    ylabel('Norm of Torque (N m)');
    legend({'||L_B||','||L_H||','||L_{GY}||','||L_{GG}||','||L_{R}||','||L_{S}||','||L_{EC}||','||L_D||'},'Location','eastoutside');
    movegui('south');
    fig = gcf;
    fig.Color = [1 1 1];
    grid on
    xlim([tstart tend]);

end

if days>1
    % final orbit plotting
    endday=t_days(end);
    tstart=endday-15/16;
    tend=endday;
    orbnitnumber=round(T(end)/period);
    str_orbit=num2str(orbitnumber);
    fin_orb_title=strcat('Magnitude of Torques during orbit','    ',str_orbit);
    figure()
    semilogy(t_days,bar_torque,'linestyle','none','marker','.','MarkerSize',msize);
    hold on
    semilogy(t_days,hyst_torque,'linestyle','none','marker','.','MarkerSize',msize);
    semilogy(t_days,gyro_torque,'linestyle','none','marker','.','MarkerSize',msize);
    semilogy(t_days,gg_torque,'linestyle','none','marker','.','MarkerSize',msize);
    semilogy(t_days,res_torque,'linestyle','none','marker','.','MarkerSize',msize);
    semilogy(t_days,sun_torque,'linestyle','none','marker','.','MarkerSize',msize);
    semilogy(t_days,eddy_torque,'linestyle','none','marker','.','MarkerSize',msize);
    semilogy(t_days,aero_torque,'linestyle','none','marker','.','MarkerSize',msize);
    title(fin_orb_title);
    xlabel('Time (days)');
    ylabel('Norm of Torque (N m)');
    legend({'||L_B||','||L_H||','||L_{GY}||','||L_{GG}||','||L_{R}||','||L_{S}||','||L_{EC}||','||L_D||'},'Location','eastoutside');
    movegui('southeast');
    fig = gcf;
    fig.Color = [1 1 1];
    grid on
    xlim([tstart tend]);
end
% if T(end)>80000
% %orbit 15 plotting
%     endtime=(86400-period):-20:(84000-period);
%     indmat=find(round(T)==endtime); 
%     ind=indmat(1);
%     indmat2=find(round(T)==86400); 
%     ind2=indmat2(1);
%     orbitnumber=round(endtime/period);
%     str_orbit=num2str(orbitnumber);
%     fin_orb_title=strcat('Magnitude of Torques during orbit','    ',str_orbit);
%     fin_orbtime=(T(ind:ind2)-endtime)/60;
%     figure()
%     semilogy(fin_orbtime,bar_torque(ind:ind2),'linestyle','none','marker','.','MarkerSize',msize);
%     hold on
%     semilogy(fin_orbtime,hyst_torque(ind:ind2),'linestyle','none','marker','.','MarkerSize',msize);
%     semilogy(fin_orbtime,gyro_torque(ind:ind2),'linestyle','none','marker','.','MarkerSize',msize);
%     semilogy(fin_orbtime,gg_torque(ind:ind2),'linestyle','none','marker','.','MarkerSize',msize);
%     semilogy(fin_orbtime,res_torque(ind:ind2),'linestyle','none','marker','.','MarkerSize',msize);
%     semilogy(fin_orbtime,eddy_torque(ind:ind2),'linestyle','none','marker','.','MarkerSize',msize);
%     title(fin_orb_title);
%     xlabel('Time (minutes)');
%     ylabel('Norm of Torque (N m)');
%     legend({'||L_B||','||L_H||','||L_{GY}||','||L_{GG}||','||L_{R}||','||L_{EC}||'},'Location','eastoutside');
%     movegui('south');
%     fig = gcf;
%     fig.Color = [1 1 1];
%     grid on
% end
% 
% if days>1
%     orbitsT=round(T(end)/period);
%     L=length(T);
%     ind=L-round(0.5*L/orbitsT);
%     orbitnumber=orbitsT-1;
%     str_orbit=num2str(orbitnumber);
%     fin_orb_title=strcat('Magnitude of Torques during orbit',' ',str_orbit);
%     fin_orbtime=(T(ind:end)-T(ind))/60;
%     figure()
%     semilogy(fin_orbtime,bar_torque(ind:end),'linestyle','none','marker','.','MarkerSize',msize);
%     hold on
%     semilogy(fin_orbtime,hyst_torque(ind:end),'linestyle','none','marker','.','MarkerSize',msize);
%     semilogy(fin_orbtime,gyro_torque(ind:end),'linestyle','none','marker','.','MarkerSize',msize);
%     semilogy(fin_orbtime,gg_torque(ind:end),'linestyle','none','marker','.','MarkerSize',msize);
%     semilogy(fin_orbtime,res_torque(ind:end),'linestyle','none','marker','.','MarkerSize',msize);
%     semilogy(fin_orbtime,eddy_torque(ind:end),'linestyle','none','marker','.','MarkerSize',msize);
%     title(fin_orb_title);
%     xlabel('Time (minutes)');
%     ylabel('Norm of Torque (N m)');
%     legend({'||L_B||','||L_H||','||L_{GY}||','||L_{GG}||','||L_{R}||','||L_{EC}||'},'Location','eastoutside');
%     movegui('southeast');
%     fig = gcf;
%     fig.Color = [1 1 1];
%     grid on
% end
% mu0=4*pi*1e-7;
% H = zeros(length(T),3);
% for i=1:length(T)
%     Bnow    = quatrotate(q(i,:),interp1(time,B_in(:,1:3),T(i)));
%     H(i,:)  = Bnow/mu0;
% end
% hyst_titlex=strcat('Hysteresis Loop in X hyst vol=',str_vol,' m^3');
% hyst_titlez=strcat('Hysteresis Loop in Z hyst vol=',str_vol,' m^3');
% 
% figure()
% plot(H(:,1),BH(:,1),'r');
% title(hyst_titlex);
% xlabel('Magnetizing Field H (A/m)');
% ylabel('Magnetic Flux Density (T)');
% movegui('southeast');
% figure()
% plot(H(:,3),BH(:,3),'r');
% title(hyst_titlez);
% xlabel('Magnetizing Field H (A/m)');
% ylabel('Magnetic Flux Density (T)');
% movegui('southwest');
% figure()
% plot(H(:,1),BH(:,1),'r');
% title(hyst_titlex);
% xlabel('Magnetizing Field H (A/m)');
% ylabel('Magnetic Flux Density (T)');
% movegui('south');
% xlim([30 33]);
% figure()
% plot(H(:,1),BH(:,1),'r');
% title(hyst_titlex);
% xlabel('Magnetizing Field H (A/m)');
% ylabel('Magnetic Flux Density (T)');
% movegui('north');
% xlim([-10 -5]);




function Beta = calcBeta(t, q, B)
y=[0, 1, 0];                    % CHANGE THIS to reflect alignment axis 
Beta = zeros(1, length(t));

for  i = 1:length(t)
    mag         = quatrotate(q(i,:),B(i,:)./norm(B(i,:)));
    Beta(:,i)   = atan2d(norm(cross(mag,y)),dot(mag,y));
end
end


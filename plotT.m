function [bar_torque,hyst_torque,gyro_torque,gg_torque,res_torque,sun_torque,eddy_torque,aero_torque]...
    =plotT(T,time,q,w,B_in,Bhyst,eci,sun,vel,m,volx,volz)
% This function is used to calculate torques for plotting
% Constants
L=length(T);
J=[6433703.78 0 0; 0 6484881.04 0;0 0 2901645.34]*(1e-9);
mu0=4*pi*1e-7;
mue= (6.67408e-11)*(5.972e24); %earth gravitational parameter
rho= 3.89e-12 ;% from http://www.braeunig.us/space/atmos.htm assuming mean solar activity
cd=2.4; % drag coeffecient from paper
rd=(1e-2)*([5;5;11.35]-[4.99;5.44;11.28]);% distance from centre of mass to geometric centre 
S=[227;227;100]*1e-4;% area vector on each face, turned from cm^2 to m^2
diagS=[227 0 0; 0 227 0; 0 0 100]*1e-4;
cr=0.8; %reflectivity of satellite (from gerhardt paper)
Ps=4.5e-6; % solar radiation pressure at earth (from gerhardt paper)
ki = diag([147.3*(2/3),147.3*(2/3),49.3]); %for eddy currents from gerhardt paper, function of A, so 2U is 2/3 3U A
% Initializing arrays for storing torques
B_samp=zeros(L,3);
bar_torque=zeros(L,1); %magnitude of torque
hyst_torque=zeros(L,1);
gyro_torque=zeros(L,1);
res_torque=zeros(L,1);
sun_torque=zeros(L,1);
eddy_torque=zeros(L,1);
aero_torque=zeros(L,1);
bar_mom=[0;m;0];
mres=[-0.0059;-0.0083;0.0004]; %residual magnetic moment from gerhardt paper needs to measured for neudose



for i=1:L
    R= quat2rotm(q(i,:));
    B_samp(i,:)=interp1(time,B_in,T(i))*R;
    bar_torque(i)=norm(cross(bar_mom,B_samp(i,:)));
    hyst_mom=[volx*Bhyst(i,1)/mu0;0;volz*Bhyst(i,3)/mu0];
    hyst_torque(i)=norm(cross(hyst_mom,B_samp(i,:)));
    % gyroscopic stiffness
    wi=w(i,:);
    Wskew=[0 -wi(3) wi(2); wi(3) 0 -wi(1); -wi(2) wi(1) 0];
    gyro_torque(i)=norm(Wskew*J*wi');
    % gravity gradient shite
    pos=interp1(time,eci,T(i))*R+rd';
    coeff=3*mue/((norm(pos))^5);
    gg_torque(i)=norm(coeff*cross(pos',J*pos'));
    % aerodynamic torque
    v=interp1(time,vel,T(i))*R;
    aero_torque(i)=norm(0.5*rho*cd*dot(S,v)*cross(v,rd));
    % residual
    res_torque(i)=norm(cross(mres,B_samp(i,:)));
    % Eddy torque
    B_body=B_samp(i,:);
    eddy_torque(i)= norm(cross(sum(dot(ki,diag(B_body/norm(B_body))))*(cross(wi,B_body)), B_body));
    % Solar torque
      suni=interp1(time,sun,T(i));
      sun_torque(i)=norm(cross(rd,(-Ps*cr*diagS*suni')));
    
end
end
    
    
    

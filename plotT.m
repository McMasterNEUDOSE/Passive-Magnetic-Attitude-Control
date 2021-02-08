function [bar_torque,hyst_torque,gyro_torque,gg_torque,aero_torque,sun_torque,eddy_torque,res_torque]...
    =plotT(T,time,q,w,B_in,Bhyst,eci,vel,dens,sun,m,hyst_l,hyst_d,nrods)
% function [bar_torque,hyst_torque,gyro_torque,gg_torque,aero_torque,sun_torque,eddy_torque,res_torque]...
%     =plotT(T,time,q,w,B_in,Bhyst,eci,vel,dens,sun,m,volx,volz)
% This function is used to calculate torques for plotting

vol = nrods*0.25*pi*hyst_l*hyst_d^2;
HystVol=[vol vol 0];
Nd  = 1/((4*hyst_l)/(sqrt(pi)*hyst_d) + 2);% demag factor

% Constants
L=length(T);
% magnetic parameters
mu0=4*pi*1e-7;
mhyst = (HystVol/mu0).*(Bhyst-B_in)./(1-Nd);
% mhyst = (HystVol/mu0).*Bhyst;
% mres=[-0.0039 -0.0055 0.0004].*ones(L,3); %residual magnetic moment from gerhardt paper scaled by 2/3
mres=[-0.0046 -0.0046 0.000575].*ones(L,3); %residual magnetic moment in A*m^2

bar_mom=[0 0 m].*ones(L,3);
% Gravity Gradient Parameters
mue= (6.67408e-11)*(5.972e24); %earth gravitational parameter
% J= diag([0.0064, 0.0065, 0.0029]);
J= diag([0.0102 , 0.0104, 0.0046]); 
% Aerodynamic Parameters
cd=2.4; % drag coeffecient from paper
rd=[-0.00244 0.00025 -0.00432].*ones(L,3);% distance from centre of mass to geometric centre 
S=[0.02 0.02 0.02].*ones(L,3);% area vector on each face, turned from cm^2 to m^2
% Solar torque parameters
diagS=diag([0.02,0.02,0.01]);
cr=0.8; %reflectivity of satellite (from gerhardt paper)
Ps=4.5e-6; % solar radiation pressure at earth (from gerhardt paper)
% eddy current shite
ki = diag([98.2,98.2,49.3]); %eddy currents from gerhardt paper scaled by 2/3 to go from 3U to 2U
% Initializing arrays for storing torques
B_samp=quatrotate(q,B_in);
v=quatrotate(q,vel);
sunrot=quatrotate(q,sun);
posrot=quatrotate(q,eci);

bar_torque=cross(bar_mom,B_samp);
hyst_torque=cross(mhyst,B_samp);
res_torque=cross(mres,B_samp);
sun_torque=cross(rd,(-Ps*cr*[0.03 0.03 0.01].*sunrot));
aero_torque=0.5*cd*dens.*dot(S,v,2).*cross(v,rd);

gg_torque=3*mue*cross(posrot,[0.0222 0.0218 0.005].*posrot)./vecnorm(posrot,2,2).^5;

eddy_torque=zeros(L,3);
gyro_torque=zeros(L,3);

for i=1:L
    wi=w(i,:);
    Wskew=[0 -wi(3) wi(2); wi(3) 0 -wi(1); -wi(2) wi(1) 0];
    gyro_torque(i,:)=Wskew*J*wi';
    B_body=B_samp(i,:);
    eddy_torque(i,:)= cross(sum(dot(ki,diag(B_body/norm(B_body))))*(cross(wi,B_body)), B_body);  
end
end
    
    
    

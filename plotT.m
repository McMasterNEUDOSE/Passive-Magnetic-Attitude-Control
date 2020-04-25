function [bar_torque,hyst_torque,gyro_torque,gg_torque,aero_torque,sun_torque,eddy_torque,res_torque]...
    =plotT(T,time,q,w,B_in,Bhyst,eci,vel,dens,sun,m,vol)
% This function is used to calculate torques for plotting

% hystH   = 95 * 1e-3;    % dimensions, m
% hystDiam= 1 * 1e-3;
% Nd  = 1/((4*hystH)/(sqrt(pi)*hystDiam) + 2);% demag factor
HystVol=[vol 0 vol];
% Constants
L=length(T);
% magnetic parameters
mu0=4*pi*1e-7;
% mhyst = (HystVol/mu0).*Bhyst./(1-Nd);
mhyst = (HystVol/mu0).*Bhyst;
mres=[-0.0039 -0.0055 0.0004].*ones(L,3); %residual magnetic moment from gerhardt paper scaled by 2/3
bar_mom=[0 m 0].*ones(L,3);
% Gravity Gradient Parameters
mue= (6.67408e-11)*(5.972e24); %earth gravitational parameter
J= diag([0.0064, 0.0065, 0.0029]);
% Aerodynamic Parameters
cd=2.4; % drag coeffecient from paper
rd=[0.0001 -0.0044 0.0007].*ones(L,3);% distance from centre of mass to geometric centre 
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

bar_torque=vecnorm(cross(bar_mom,B_samp),2,2);
hyst_torque=vecnorm(cross(mhyst,B_samp),2,2);
res_torque=vecnorm(cross(mres,B_samp),2,2);
sun_torque=vecnorm(cross(rd,(-Ps*cr*[0.03 0.03 0.01].*sunrot)),2,2);
aero_torque=vecnorm(dens.*dot(S,v,2).*cross(v,rd),2,2);

gg_torque=vecnorm(3*mue*cross(posrot,[0.0222 0.0218 0.005].*posrot)./vecnorm(posrot,2,2).^5,2,2);

eddy_torque=zeros(L,1);
gyro_torque=zeros(L,1);

for i=1:L
    wi=w(i,:);
    Wskew=[0 -wi(3) wi(2); wi(3) 0 -wi(1); -wi(2) wi(1) 0];
    gyro_torque(i)=norm(Wskew*J*wi');
    B_body=B_samp(i,:);
    eddy_torque(i)= norm(cross(sum(dot(ki,diag(B_body/norm(B_body))))*(cross(wi,B_body)), B_body));  
end
end
    
    
    

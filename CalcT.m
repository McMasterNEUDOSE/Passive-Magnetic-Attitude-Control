function T=CalcT(B,BH,pos,sun,w,vel,rho,m,vol)

% magnetic shit
mbar=[0;m;0];
mres=[-0.0039;-0.0055;0.0004]; %residual magnetic moment from gerhardt paper scaled by 2/3
mu0=4*pi*1e-7;
HystVol=[vol 0 vol];

% hystH   = 95 * 1e-3;    % dimensions, m
% hystDiam= 1 * 1e-3;
% Nd  = 1/((4*hystH)/(sqrt(pi)*hystDiam) + 2);% demag factor

% mhyst = (HystVol/mu0).*(BH'-B)./(1-Nd);
mhyst = (HystVol/mu0).*(BH'-B);


Mtot=mbar+mhyst'-mres;
% Gravity Gradient parameters
J= diag([0.0064, 0.0064, 0.0029]);
mue= (6.67408e-11)*(5.972e24); %earth gravitational parameter
r_mag=norm(pos); coeff=3*mue/(r_mag^5); %constants for gravity gradient
% Aerodynamic parameters
cd=2.4; % drag coeffecient from paper
rd=[0.0001;-0.0044;0.0007];% distance from centre of mass to geometric centre 
rdmag=norm(rd);
rdunit=rd/rdmag;
vmag=norm(vel);
vunit=vel/vmag;
pos=pos+rd';
S=[0.02;0.02;0.01];% area vector on each face
% Solar torque parameters
diagS=diag([0.02,0.02,0.01]);
cr=0.8; %reflectivity of satellite (from gerhardt paper)
Ps=4.5e-6; % solar radiation pressure at earth (from gerhardt paper)
% Eddy current parameters
ki = diag([98.2,98.2,49.3]); %eddy currents from gerhardt paper scaled by 2/3 to go from 3U to 2U
% torques
mag_torque=cross(Mtot,B');
gg_torque=coeff*cross(pos',J*pos');
aero_torque=rho*(vmag^2)*rdmag*abs(dot(S,vunit))*cross(vunit,rdunit)';
sun_torque=cross(rd,(-Ps*cr*diagS*sun'));
eddy_torque= cross(sum(dot(ki,diag(B/norm(B))))*(cross(w,B)), B);
T=mag_torque+gg_torque+sun_torque+eddy_torque'+aero_torque;
end
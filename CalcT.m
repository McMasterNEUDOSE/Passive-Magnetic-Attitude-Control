function T=CalcT(B,BH,w,pos,sun,vel,rho,m,hyst_l,hyst_d,nrods,mres)
% function T=CalcT(B,BH,w,pos,sun,vel,m,hyst_l,hyst_d,nrods,mres)
% function T=CalcT(B,BH,w,pos,sun,m,hyst_l,hyst_d,nrods,mres)

% General Constants
S=(1e-4)*[200;200;100];% area vector on each face
% S=[0.02;0.02;0.01];% area vector on each face

rd=[-0.00244;0.00025;-0.00432];% distance from centre of mass to geometric centre 

% Bar Magnet 
mbar=[0;0;m];

% Hysteresis rods
mu0=4*pi*1e-7;
vol = nrods*0.25*pi*hyst_l*hyst_d^2;
HystVol=[vol vol 0];
Nd  = 1/((4*hyst_l)/(sqrt(pi)*hyst_d) + 2);% demag factor
mhyst = (HystVol/mu0).*(BH'-B)./(1-Nd);

% Gravity Gradient torque parameters
J= diag([0.0102 , 0.0104, 0.0046]); 
mue= (6.67408e-11)*(5.972e24); %earth gravitational parameter
r_mag=norm(pos); coeff=3*mue/(r_mag^5); %constants for gravity gradient

% Solar torque parameters
cr=0.8; %reflectivity of satellite (from gerhardt paper)
Ps=4.5e-6; % solar radiation pressure at earth (from gerhardt paper)
diagS=diag([0.02,0.02,0.01]);

%Aerodynamic torque parameters
cd=2.4; % drag coeffecient from paper
% cd=2.2; % drag coeffecient 

rdmag=norm(rd);
rdunit=rd/rdmag;
vmag=norm(vel);
vunit=vel/vmag;


% Eddy torque stuff
ki = [98.2;98.2;49.3]; %eddy currents from gerhardt paper scaled by 2/3 to go from 3U to 2U

% ----------------------------------------
% Torques
bar_torque = cross(mbar,B');
hyst_torque = cross(mhyst',B');
res_torque = cross(mres,B');
sun_torque=cross(rd,(-Ps*cr*diagS*sun'));
gg_torque=coeff*cross(pos',J*pos');
eddy_torque= dot(B/norm(B),ki)*cross(cross(w,B),B)';
aero_torque=0.5*cd*rho*(vmag^2)*rdmag*abs(dot(S,vunit))*cross(vunit,rdunit')';
% aero_torque=0.5*rho*(vmag^2)*cd*rdmag*abs(dot(S,vunit))*cross(vunit,rdunit)';
% aero_torque=0.5*rho*cd*abs(dot(vel,S))*cross(vel,rd)';
% aero_torque=cross(0.5*rho*cd*abs(dot(vel,S))*vel,rd)';
% aero_torque=cross(rd,(aero_coeff*vunit'));
% aero_torque=cross(rd,(aero_coeff*vunit'));
% aero_torque = aero_coeff*aero_T_vector;




T=bar_torque+hyst_torque+res_torque+sun_torque+gg_torque+eddy_torque+aero_torque;
% T=bar_torque+hyst_torque+res_torque+sun_torque+gg_torque+eddy_torque;
end

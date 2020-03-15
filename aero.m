function torque=aero(Q,V)
% Q is the current attitude, V is the current velocity
% rotate V into the body fixed from body inertial frame
rho= 3.89e-12 ;% from http://www.braeunig.us/space/atmos.htm assuming mean solar activity
cd=2.4; % drag coeffecient from paper, needs to be updated
% centre of mass: (4.99,5.44,11.28) cm, while geometric centre is:
% (5,5,11.35)
rd=(1e-2)*([5;5;11.35]-[4.99;5.44;11.28]);% distance from centre of mass to geometric centre 
S=[227;227;100]*1e-4;% area vector on each face, turned from cm^2 to m^2
V=quatrotate(Q,V); % rotate velocity vector into body frame
k=0.5*rho*cd*dot(S,V);
torque=k*cross(V,rd);
end
function torque=calcT(B_body,BH,pos,sun,w,vel,m,volx,volz)
mu0=4*pi*1e-7;
mbar=[0;m;0];
J=[6433703.78 0 0; 0 6484881.04 0;0 0 2901645.34]*(1e-9);
mue= (6.67408e-11)*(5.972e24); %earth gravitational parameter
r_mag=norm(pos); coeff=3*mue/(r_mag^5); %constants for gravity gradient
rho= 3.89e-12 ;% from http://www.braeunig.us/space/atmos.htm assuming mean solar activity
cd=2.4; % drag coeffecient from paper
rd=(1e-2)*([5;5;11.35]-[4.99;5.44;11.28]);% distance from centre of mass to geometric centre 
pos=pos+rd';
S=[227;227;100]*1e-4;% area vector on each face, turned from cm^2 to m^2
diagS=[227 0 0; 0 227 0; 0 0 100]*1e-4;
cr=0.8; %reflectivity of satellite (from gerhardt paper)
Ps=4.5e-6; % solar radiation pressure at earth (from gerhardt paper)
ki = diag([147.3*(2/3),147.3*(2/3),49.3]); %for eddy currents from gerhardt paper, function of A, so 2U is 2/3 3U A
mhyst=[volx*BH(1)/mu0;0;volz*BH(3)/mu0];
mres=[-0.0059;-0.0083;0.0004]; %residual magnetic moment from gerhardt paper needs to measured for neudose
bar_torque=cross(mbar,B_body)';
hyst_torque=cross(mhyst,B_body)';
gg_torque=coeff*cross(pos',J*pos');
res_torque=cross(mres,B_body)';
sun_torque=cross(rd,(-Ps*cr*diagS*sun'));
eddy_torque= cross(sum(dot(ki,diag(B_body/norm(B_body))))*(cross(w,B_body)), B_body);
aero_torque=0.5*rho*cd*dot(S,vel)*cross(vel,rd)';
torque=bar_torque+hyst_torque+gg_torque+res_torque+sun_torque+eddy_torque'+aero_torque;
end
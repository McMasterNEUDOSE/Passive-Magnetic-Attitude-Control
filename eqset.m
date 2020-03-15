function dyneq=eqset(x,gg,B_inertiali,B_inertial0,velocity,m_bar,vol)
%q is current attitude, w is angular velocity, T is external torque
% gg is external torque in body inertial frame
% velocity is orbital velocity of satellite
w=[x(1);x(2);x(3)];
q=[x(4);x(5);x(6);x(7)];
gg=quatrotate(q',gg);
W=[0 w'];
B_bodyi=quatrotate(q',B_inertiali);
B_body0=quatrotate(q',B_inertial0);

mag_torque=magtorque(B_bodyi',B_body0',m_bar,vol);
aero_torque=aero(q',velocity);

torque=gg'+mag_torque'+aero_torque';

% old moment of inertia
% Moment of inertia matrix from Juan
%J=[6433703.78 1176.84 -285846.68; 1176.84 6484881.04 18874.53; -259656.12 18874.53 2901645.34]*(1e-9);
J=[6433703.78 0 0; 0 6484881.04 0;0 0 2901645.34]*(1e-9);


Wskew=[0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
Wdot = J\(torque-Wskew*J*w); %Equation 2.1 from gerhardt dissertation

q0=q(1);
q1=q(2);
q2=q(3);
q3=q(4);
Qdot=[q0 -q1 -q2 -q3; q1 q0 -q3 q2; q2 q3 q0 -q1;q3 -q2 q1 q0]*W'; %Equation 8.4 from gerhardt dissertation
dyneq=[Wdot; Qdot];
end

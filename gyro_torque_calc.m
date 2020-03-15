function torque=gyro_torque_calc(w)
% Moment of inertia matrix from Juan
%J=[6433703.78 1176.84 -285846.68; 1176.84 6484881.04 18874.53; -259656.12 18874.53 2901645.34]*(1e-9);
J=[6433703.78 0 0; 0 6484881.04 0;0 0 2901645.34]*(1e-9);
Wskew=[0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
torque=-Wskew*J*w';
end
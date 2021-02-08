clear all;
% This script calculates thet initial error angle between the
% bar magnet and the magnetic field
B0 = [-7.22723018278003e-07,-9.06296946062727e-06,3.25701853359517e-05];
B0_norm = B0/norm(B0);
B0_ang = 90*B0_norm;

% Below are some initial attitudes 
% att_init=[30;88.8;-104]*pi/180; %initial beta  90 deg
% att_init=[-170.7;1.3;35]*pi/180; %initial beta     173.7513 deg
% att_init=[-166.4;1.3;35]*pi/180; %initial beta   178.0493 deg
% att_init=[15;177;-150]*pi/180; %  initial beta 178.1416 deg
% att_init=[15;179.5;-150]*pi/180; %  initial beta   179.0905 deg
% att_init=[-150;58.8;-102.6]*pi/180; %initial beta 121.3192 deg
% att_init=[15.5;-1.2;-30]*pi/180; %initial beta  0.0555 deg
att_init=[-165.7;1.3;35]*pi/180; %initial beta   178.7483 deg

att0=angle2quat(att_init(1),att_init(2),att_init(3),'XYZ');
Beta_init = Beta_check(att0,B0)

function Beta = Beta_check(q, B)
z=[0, 0, 1];                    % CHANGE THIS to reflect alignment axis 

mag    = quatrotate(q,B/norm(B));
Beta   = atan2d(norm(cross(mag,z)),dot(mag,z));
end

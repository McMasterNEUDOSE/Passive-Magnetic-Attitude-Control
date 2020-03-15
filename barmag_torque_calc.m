function torque=barmag_torque_calc(B_body,mbar)
m=[0;mbar;0];
torque=cross(m,B_body')';
end
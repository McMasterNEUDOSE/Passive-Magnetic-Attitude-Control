function torque=magtorque(B_bodyi,B_body0,m_bar,vol)
% HYSTERESIS STUFF GOES HERE
% qconj rotates the bar magnet axes from the body fixed frame to the body
% inertial frame
mu0=4*pi*1e-7;
H=B_bodyi/mu0;
% Hysteresis Rods:
%From paper, L/D ~ 100-300
Bx=B_bodyi(1);
By=B_bodyi(2);
Bz=B_bodyi(3);

Bs=0.3;
Br=6.0618e-4; % from gerhardt dissertation
Hc=0.3381; %from gerhardt dissertation
mu_hyst=1.5e4;

k=(1/Hc)*tan(pi*Br/(2*Bs));
signx=(B_body0(1)-B_bodyi(1))/abs((B_body0(1)-B_bodyi(1)));
signz=(B_body0(3)-B_bodyi(3))/abs((B_body0(3)-B_bodyi(3)));
Bx=(2*Bs/pi)*atan(k*(H(1)+(signx*Hc))); 
Bz=(2*Bs/pi)*atan(k*(H(3)+(signz*Hc)));    
% INVESTIGATE EQUATION 2.7 FROM GERHARDT DISSERTATION
mhystx=vol*Bx/mu0;
mhystz=vol*Bz/mu0;
m=[mhystx ;m_bar; mhystz]; % Magnetic moment of bar magnet only
% m=quatrotate(qconj,m');
torque=cross(m,B_bodyi');
end

% function torque=magtorque(B_bodyi,B_body0,m_bar,vol)
% % HYSTERESIS STUFF GOES HERE
% % qconj rotates the bar magnet axes from the body fixed frame to the body
% % inertial frame
% mu0=4*pi*1e-7;
% H=B_bodyi/mu0;
% % Hysteresis Rods:
% %From paper, L/D ~ 100-300
% Bx=B_bodyi(1);
% By=B_bodyi(2);
% Bz=B_bodyi(3);
% 
% Bs=0.3;
% Br=6.0618e-4; % from gerhardt dissertation
% Hc=0.3381; %from gerhardt dissertation
% mu_hyst=1.5e4;
% 
% k=(1/Hc)*tan(pi*Br/(2*Bs)); % equation 8.16 from gerhardt dissertation
% if B_bodyi(1)>B_body0(1)
%     Bx=(2*Bs/pi)*atan(k*(H(1)-Hc));   
% elseif B_bodyi(1)<B_body0(1)
%     Bx=(2*Bs/pi)*atan(k*(H(1)+Hc));
% end    
% if B_bodyi(3)>B_body0(3)
%     Bz=(2*Bs/pi)*atan(k*(H(3)-Hc));    
% elseif B_bodyi(3)<B_body0(3)
%     Bz=(2*Bs/pi)*atan(k*(H(3)+Hc));
% end
% % INVESTIGATE EQUATION 2.7 FROM GERHARDT DISSERTATION
% mhystx=vol*Bx/mu0;
% mhystz=vol*Bz/mu0;
% m=[mhystx ;m_bar; mhystz]; % Magnetic moment of bar magnet only
% % m=quatrotate(qconj,m');
% torque=cross(m,B_bodyi');
% end
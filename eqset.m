function dyneq=eqset(x,t,time,B,Bdot,eci,sunlight,v,dens,m,hyst_l,hyst_d,nrods,mres)
w=x(1:3);
q=x(4:7)';
BH=x(8:10);

% Interpolations/rotations
rho=lininterp1(dens(:,1),dens(:,2),t);
B_body=quatrotate(q,lininterp1(time,B,t));
Bdot=quatrotate(q,lininterp1(time,Bdot,t));
sun=quatrotate(q,lininterp1(time,sunlight,t));
pos=quatrotate(q,lininterp1(time,eci,t));
vel=quatrotate(q,lininterp1(time,v,t));
Bdot= (cross(-w, B_body) + Bdot)'; % Using the transport theorem

% constants
mu0=4*pi*1e-7;
Bs=0.3;
Br=6.0618e-4; % from gerhardt dissertation
Hc=0.3381; %from gerhardt dissertation
BsFac=0.5*pi/Bs;
k=(1/Hc)*tan(0.5*pi*Br/Bs);


Jvec=[0.0102 ; 0.0104 ; 0.0046]; % Moment of inertia
w_vec=[(Jvec(3)-Jvec(2))*w(2)*w(3);(Jvec(1)-Jvec(3))*w(1)*w(3);(Jvec(2)-Jvec(1))*w(2)*w(1)];

% Setting differential equations
T=CalcT(B_body,BH,w,pos,sun,vel,rho,m,hyst_l,hyst_d,nrods,mres); % External torques
Wdot=(T-w_vec)./Jvec;
W=[0 w'];

Qdot=0.5*quatmultiply(q,W);
dBH = (k/BsFac)*(cos(BsFac*BH).^2).*(Bdot/mu0).*...
    (((2*Hc)^-1*((B_body'/mu0)-(k^-1)*tan(BsFac*BH)+sign((Bdot/mu0))*Hc)).^2);
dyneq = [Wdot;Qdot';dBH];
end

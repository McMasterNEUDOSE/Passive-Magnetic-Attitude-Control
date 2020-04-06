function dBh=flatley(Bdot,BH,B_body)
mu0=4*pi*1e-7;
Bs=0.3;
Br=6.0618e-4; % from gerhardt dissertation
Hc=0.3381; %from gerhardt dissertation
BsFac=0.5*pi/Bs;
k=(1/Hc)*tan(0.5*pi*Br/Bs);
dBh = (k/BsFac)*(cos(BsFac*BH).^2).*(Bdot/mu0).*...
    (((2*Hc)^-1*((B_body/mu0)-(k^-1)*tan(BsFac*BH)+sign((Bdot/mu0))*Hc)).^2);
end
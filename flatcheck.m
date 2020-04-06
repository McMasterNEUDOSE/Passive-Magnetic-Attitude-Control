function BH=flatcheck(B_body,Bdot,B_flat)
mu0=4*pi*1e-7;
Bs=0.3;
Br=6.0618e-4; % from gerhardt dissertation
Hc=0.3381; %from gerhardt dissertation
BsFac=0.5*pi/Bs;
k=(1/Hc)*tan(0.5*pi*Br/Bs);
H=B_body/mu0;
BHHi= (2*Bs/pi)*atan(k*(H+Hc));
BHLo= (2*Bs/pi)*atan(k*(H-Hc));
BH=min(max(B_flat,BHLo),BHHi);
% BHx=min(max(B_flat(1)),BHLo(1)),BHHi(1));
% BHy=min(max(B_flat(2),BHLo(2)),BHHi(2));
% BHz=min(max(B_flat(3),BHLo(3)),BHHi(3));
% BH=[BHx BHy BHz];
end
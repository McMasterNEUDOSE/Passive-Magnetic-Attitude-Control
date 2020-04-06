function dyneq=eqset(x,t,time,B,m,Bdot,eci,sunlight,vel,volx,volz)
w=[x(1);x(2);x(3)];
q=[x(4) x(5) x(6) x(7)];
% BH_flat=[x(8) x(9) x(10)]; % ADD IN BOUNDING FUNCTION WITH INVTAN
BH=[x(8) x(9) x(10)]; % ADD IN BOUNDING FUNCTION WITH INVTAN
BH_flat=BH;
R= quat2rotm(q); 
B_body=interp1(time,B,t)*R;
Bdot=interp1(time,Bdot,t)*R;
pos=interp1(time,eci,t)*R;
sun=interp1(time,sunlight,t)*R;
vel=interp1(time,vel,t)*R;
Bdot= cross(-w, B_body) + Bdot;
BH=flatcheck(B_body,Bdot,BH_flat);

J=[6433703.78 0 0; 0 6484881.04 0;0 0 2901645.34]*(1e-9);
Wskew=[0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
Wdot = J\(calcT(B_body,BH,pos,sun,w',vel,m,volx,volz)-Wskew*J*w); %Equation 2.1 from gerhardt dissertation
W=[0 w'];
Qdot=0.5*quatmultiply(q,W); %Equation 8.4 from gerhardt dissertation
dBH=flatley(Bdot,BH,B_body);
dyneq=[Wdot; Qdot';dBH'];
end
function [Q_new,W_new]=QWsel(t_sim,t_ode45,Q,tstep,W)
% this function generates an attitude quaternion matrix, that will have
% equal amounts of values as there are timesteps in the STK simulation,
% so each magnetic field vector can be rotated into the body fixed frame
% ---------------------------------------------------------------------
% t_sim is the time vector imported which gives the time values for the STK
% simulation
% t_ode45 is the tiem vector that is calculated wehn using ODE45
% Q is the attitude quaternion vector solved for using ODE45
% tstep is the timestep of the simulation
% ---------------------------------------------------------------------
% these if statements determine how many digits to round the time vector to
% depending on the time step
if tstep==0.1
    dig=3;
elseif tstep==1
    dig=2;
end

t_round=round(t_ode45,dig);
L=length(t_sim);
indices=zeros(L,1);
Q_new=zeros(L,4);
W_new=zeros(L,3);
for i=1:L
    diff=t_round-t_sim(i);
    absdiff=abs(diff);
    val=find(absdiff==min(absdiff));
    if length(val>1)
        val=val(1);
    end
    indices(i,1)=val;
    Q_new(i,:)=Q(indices(i),:);
    W_new(i,:)=W(indices(i),:);
end

end

    
    

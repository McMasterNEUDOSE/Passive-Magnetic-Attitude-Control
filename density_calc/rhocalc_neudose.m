% dat=table2array(load('NeudoseLLa.mat'));
close all;
days=180;
timestep=180;
% lat=dat(:,1);
% lon=dat(:,2);
% alt=dat(:,3)*1e3; %convert to m
% alt = alt*1e3;
L=length(alt);
f107Average=(168.5)*ones(L,1);
f107Daily=(128.7)*ones(L,1);
magneticIndex=48*ones(L,1);
year=2022*ones(L,1);
dayofyear=zeros(L,1);
UTsec=zeros(L,1);
totday=172;
totsec=16*3600;
for i=1:L
    if totsec<86400
        totsec=totsec+timestep;
    else 
        totsec=0;
        totday=totday+1;
        
    end
if totday > 365
    totday = totday - 365;
end
dayofyear(i)=totday;
UTsec(i)=totsec;
end
            
[T rho] = atmosnrlmsise00(alt,lat,lon, year, dayofyear, UTsec, f107Average, f107Daily, magneticIndex);
t=0:timestep:(days*86400);
rho_air=[t' rho(3:end,6)];
rho_air(end,2)=rho_air(end-2,2);
rho_air(end-1,2)=rho_air(end-2,2);

save_name = strcat('neudoserho',num2str(timestep),'s_',num2str(days),'day.mat');
save(save_name,'rho_air');
figure()
plot(rho_air(:,2));
figure()
semilogy(rho_air(:,2));

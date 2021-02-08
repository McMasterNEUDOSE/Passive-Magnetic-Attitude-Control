close all;
clear all;
days=28;

if days==28
    LLA_data = readmatrix('Neudose_LLA_Position_28day.csv');
    timestep=60;
elseif days ==60
    LLA_data = readmatrix('Neudose_LLA_Position_60day.csv');
    timestep=60;
elseif days ==120
    LLA_data = readmatrix('Neudose_LLA_Position_120day.csv');
    timestep=120;
elseif days ==365
    LLA_data = readmatrix('Neudose_LLA_Position_365day.csv');
    timestep=120;
end
    
lat = LLA_data(:,2);
lon = LLA_data(:,3);
alt = (1e3)*LLA_data(:,4); % Altitude is in km convert to m

    
L=length(alt);
% f107Average=(168.5)*ones(L,1); % CSSWE Value
% f107Average=(146.5)*ones(L,1); % 2014 Average Value (worst recent year)
% f107Average=(69.7)*ones(L,1); % June 2020 Value : https://www.spaceweather.gc.ca/solarflux/sx-5-mavg-en.php
f107Average=(87.73)*ones(L,1); % June 2022 predicted Value : https://www.swpc.noaa.gov/products/solar-cycle-progression
% f107Average=(110.4)*ones(L,1); % June 2023 predicted Value : https://www.swpc.noaa.gov/products/solar-cycle-progression



% f107Daily=(128.7)*ones(L,1); % CSSWE Value
% f107Daily=(146.4)*ones(L,1); % 2014 Average Value (worst recent year)
% f107Daily=(67.6)*ones(L,1); %June 21 2020 Value : https://spaceweather.gc.ca/solarflux/sx-5-flux-en.php
f107Daily=(87.8)*ones(L,1); % June 2022 predicted Value : https://www.swpc.noaa.gov/products/solar-cycle-progression
% f107Daily=(109.9)*ones(L,1); % June 2023 predicted Value : https://www.swpc.noaa.gov/products/solar-cycle-progression




magneticIndex=48*ones(L,1);
% magneticIndex=15*ones(L,1);

year=2020*ones(L,1);
dayofyear=zeros(L,1);
UTsec=zeros(L,1);
totday=173; %June 21 2020 day of year
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
            
[T ,rho] = atmosnrlmsise00(alt,lat,lon, year, dayofyear, UTsec, f107Average, f107Daily, magneticIndex);
t=0:timestep:(days*86400);
rho_air=[t' rho(:,6)];

save_name = strcat('neudoserho',num2str(timestep),'s_',num2str(days),'day.mat');
save(save_name,'rho_air');
figure()
plot(rho_air(:,2));
figure()
semilogy(rho_air(:,2));

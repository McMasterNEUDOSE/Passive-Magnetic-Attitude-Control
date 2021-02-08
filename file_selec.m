function [data_fname,rho_fname] = file_selec(days,tstep,rho_year)
str_days=num2str(days); str_tstep=num2str(tstep);  
path='data_';
data_fname=strcat(path,str_days,'days_',str_tstep,'s.mat');

if rho_year == 2022
    rho_folder = '2022 Solar Activity Density/';
elseif rho_year == 2020
    rho_folder = '2020 Solar Activity Density/';
elseif rho_year == 2023
    rho_folder = '2023 Solar Activity Density/';
elseif rho_year == 2014
    rho_folder = '2014 Solar Activity Density/';
end
if days<= 28
    dens='neudoserho60s_28day.mat';
elseif days <=60
     dens='neudoserho60s_60day.mat';
elseif days <=120
     dens='neudoserho120s_120day.mat';
elseif days <=180
     dens='neudoserho180s_180day.mat';
elseif days<= 365
    dens='neudoserho120s_365day.mat';
end

rho_fname = strcat(rho_folder,dens);

end
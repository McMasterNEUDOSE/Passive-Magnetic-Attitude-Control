clear all; close all; 
tic;
days=365; % Number of days
tstep = 60; % desired timestep for attitude file


q_fname=strcat('q_save_',num2str(days),'days.mat'); % File name of quaternions save 
w_fname=strcat('w_save_',num2str(days),'days.mat'); % File name of angular velocity save
t_fname=strcat('t_save_',num2str(days),'days.mat'); % Filename of time save


q_downsamp=cell2mat(struct2cell(load(q_fname))); % Import quaternions
w_downsamp=(180/pi)*cell2mat(struct2cell(load(w_fname))); % Import angular velocities, comvert to deg/s
t_downsamp=cell2mat(struct2cell(load(t_fname))); % Import time

t_end = t_downsamp(end); %Ending time, final value in downsampled time array
t_STK = 0:tstep:t_end; % Time vector for STK
L = round(t_end/tstep); % number of points

q_interp = zeros(L,4); % Initialize interpolated quaternions
w_interp = zeros(L,3); % Initialize interpolated angular velocites
for i=1:L
    t_curr = t_STK(i);
    q_interp(i,:)= lininterp1(t_downsamp,q_downsamp,t_curr); % Interpolate
    w_interp(i,:)= lininterp1(t_downsamp,w_downsamp,t_curr);
end

% -------------------------------------
% TESTING STUFF
%  t_STK = 0:tstep:1000;
%  q_interp = zeros(1000,4);
%  w_interp = zeros(1000,3);
% 
% for i=1:1001
%     t_curr = t_STK(i);
%     q_interp(i,:)= lininterp1(t_downsamp,q_downsamp,t_curr);
%     w_interp(i,:)= lininterp1(t_downsamp,w_downsamp,t_curr);
% end
% -------------------------------------

q_STK = [q_interp(:,2:4) q_interp(:,1)]; % Convert from Scalar firt to scalar last
w_STK = w_interp;

% Saving attitude file
a_fname = strcat('attitude_',num2str(days),'days_',num2str(tstep),'s.txt'); % Set file name of attitude save file

a_save_STK = [t_STK' q_STK w_STK]; % Maxe combined matrix of time, quaternions, angular velocities
writematrix(a_save_STK,a_fname,'delimiter','space'); % Save matric to space delimted text file


toc;



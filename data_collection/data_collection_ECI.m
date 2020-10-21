clear all;
close all;
tic
app = actxserver('STK11.application');
app.visible = 2;
root = app.Personality2; 

days=180;
timestep=60; %Timestep of data acquisition in seconds
start_time='21 Jun 2022 16:00:00.000';
stop_time='18 Dec 2022 16:00:00.000';
scenario = root.Children.New('eScenario','Data_Collection');
scenario.SetTimePeriod(start_time,stop_time);
scenario.StartTime = start_time;
scenario.StopTime = stop_time;
root.ExecuteCommand('Animate * Reset');
str_days=strcat(num2str(days),'days');
str_tstep=strcat(num2str(timestep),'s');
savename=strcat('data_',str_days,'_',str_tstep,'.mat');

% Adding the satellite to the scenario
satellite = scenario.Children.New('eSatellite','Neudose'); 


% IAgSatellite satellite: Satellite object
model = satellite.VO.Model;
model.ModelData.Filename = 'STKData\VO\Models\Space\cubesat_2u.dae';
orbitmarker = model.OrbitMarker;
orbitmarker.SetMarkerImageFile('C:\Program Files\AGI\STK 11\STKData\VO\Markers\Satellite.ppm');
orbitmarker.MarkerData.IsTransparent = true;
orbitmarker.PixelSize = 18;
orbitmarker.OrientationMode = 'eVOMarkerOrientationFollowDirection';

%Keplerian orbit defined (try to find just circular next time)
keplerian = satellite.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical');
keplerian.SizeShapeType = 'eSizeShapeAltitude';
keplerian.LocationType = 'eLocationTrueAnomaly';
keplerian.Orientation.AscNodeType = 'eAscNodeLAN';
keplerian.SizeShape.PerigeeAltitude = 401.4;
keplerian.SizeShape.ApogeeAltitude = 408;
keplerian.Orientation.Inclination = 51.64;
%Same parameters as the iss https://heavens-above.com/orbit.aspx?satid=25544
keplerian.Orientation.ArgOfPerigee = 223.6999;
keplerian.Orientation.AscNode.Value = 115.2267;
keplerian.Location.Value = 136.3668;
satellite.Propagator.InitialState.Representation.Assign(keplerian);
satellite.Propagator.Propagate;
% IAgSatellite satellite: Satellite object
% satellite.Attitude.('eprofileECIVelNadir');
% basic = satellite.Attitude.Basic;
root.ExecuteCommand('SetAttitude */Satellite/Neudose Profile InertFix Euler 0 0 0 123');

% %below determines the relative rates at which it is spinning per axis
% basic.Profile.Body.AssignXYZ(1,1,1)
% % Below determines which axes it is spinning about
% basic.Profile.Inertial.AssignXYZ(1,1,1);
% basic.Profile.Rate = 0;  % rev/se

root.UnitPreferences.Item('DateFormat').SetCurrentUnit('EpSec');

%--------------------------------------------------------------------------------

% % ECI Magnetic Field in nT
mag_elems_eci = {'B Field ECI-X';'B Field ECI-Y'; 'B Field ECI-Z'}; 
satDP_ecI=satellite.DataProviders.GetDataPrvTimeVarFromPath('SEET Magnetic Field').Exec(scenario.StartTime,scenario.StopTime,timestep); 
bx_eci = cell2mat(satDP_ecI.DataSets.Item(cast(16,'int32')).GetValues); 
by_eci = cell2mat(satDP_ecI.DataSets.Item(cast(17,'int32')).GetValues); 
bz_eci = cell2mat(satDP_ecI.DataSets.Item(cast(18,'int32')).GetValues);
magdata_eci = [bx_eci by_eci bz_eci]*1e-9; %convert to nT
bdot_eci=diff(magdata_eci)/timestep;
bdot_eci=[bdot_eci; bdot_eci(end,:)/timestep];


root.UnitPreferences.Item('DateFormat').SetCurrentUnit('EpSec'); 
pos = satellite.DataProviders.Item('Vectors(Body)').Group.Item('Position').Exec(scenario.StartTime,scenario.StopTime,timestep);
vel = satellite.DataProviders.Item('Vectors(Body)').Group.Item('Velocity').Exec(scenario.StartTime,scenario.StopTime,timestep);
time = cell2mat(pos.DataSets.GetDataSetByName('Time').GetValues); 
x = cell2mat(pos.DataSets.GetDataSetByName('x').GetValues); 
y = cell2mat(pos.DataSets.GetDataSetByName('y').GetValues); 
z = cell2mat(pos.DataSets.GetDataSetByName('z').GetValues); 
vx = cell2mat(vel.DataSets.GetDataSetByName('x').GetValues); 
vy = cell2mat(vel.DataSets.GetDataSetByName('y').GetValues); 
vz = cell2mat(vel.DataSets.GetDataSetByName('z').GetValues);

pos_eci=[x y z]*1e3; %convert to m from km
vel_eci=[vx vy vz]*1e3; %covnert to m/s from km/s



% sun stuff
root.UnitPreferences.Item('DateFormat').SetCurrentUnit('EpSec'); 
inertial = satellite.DataProviders.Item('Vectors(Inertial)').Group.Item('Sun').Exec(scenario.StartTime,scenario.StopTime,timestep);
s_x = cell2mat(inertial.DataSets.GetDataSetByName('x').GetValues); 
s_y = cell2mat(inertial.DataSets.GetDataSetByName('y').GetValues); 
s_z = cell2mat(inertial.DataSets.GetDataSetByName('z').GetValues); 
sun_pos=[s_x s_y s_z]*1e3; %convert to m from km
sun_unit=sun_pos./vecnorm(sun_pos,2,2);

% Solar intensity
root.UnitPreferences.Item('DateFormat').SetCurrentUnit('EpSec'); 
intensity = satellite.DataProviders.Item('Solar Intensity').Exec(scenario.StartTime,scenario.StopTime,timestep);
int=cell2mat(intensity.DataSets.GetDataSetByName('Intensity').GetValues);
int=int/100;

data=[time pos_eci vel_eci magdata_eci bdot_eci sun_unit int];
save(savename,'data');
toc;
% 

# Passive-Magnetic-Attitude-Control

#ECEF_basic_propagation.m is the main code that propagates the attitude, it also plots the pointing error, angular velocities and all modelled torques on the satellite
#eqset.m Sets up the 10 coupled ODES needed to propagate the attitude, angular velocity and hysteresis model (flatley). This utilizes many of the other functions
#ode113quat and ode420 are modified versions of ode113  and ode45 which include a normalization step for the quaternions
#CalcT caluclates all torques on the satellite
#calcBeta caluclates the pointing error of the satellite
#PlotT calculates all the torques for plotting

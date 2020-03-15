# Passive-Magnetic-Attitude-Control

#ECEF_basic_propagation.m is the main code that propagates the attitude, it also plots the pointing error, angular velocities and all modelled torques on the satellite
#aero.m caluclates the aerodynamic torque
#barmag_torque_calc.m calculates the permanent magnet torque for plotting
#equset.m Sets up the 7 coupled ODES needed to propagate the attitude and the angular velocity. This utilizes many of the other functions
#fileselec.m selects the data files to be used depending on how long the simulation will porpagate for. Currently there is only data for 1,2,4,5,7,10,14 days of propagation. NOTE these diorectories are for my computer, and would need to be changed to run on another computer
#gyro_torque_calc.m calculated the gyroscopic stifness torque for plotting
#hyst_torque_calc.m calculates the hysteresis torque for plotting
#magtorque.m is used to calculaute the total magnetic torque for attitude propagation
#QWsel is used in order to deal with the difference in time steps used byt the ODE solver and the data from STK

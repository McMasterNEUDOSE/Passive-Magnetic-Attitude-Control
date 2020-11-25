import time

"""
reference : https://help.agi.com/stk/11.0.1/Content/stk/importfiles-01.htm#quatAngVels
and https://help.agi.com/stk/11.0.1/Content/stk/example_AttitudeFile4.htm
"""

days = 365 # Number of days in simulation
tstep = 60 # Timestep in minutes
n_points = days*86400*days +1 # Total number of points (starts at 0 so add 1)

a_fname ='attitude_'+str(days)+'days_'+str(tstep)+'s.txt' # Naming convention for the txt files

# Header below
header="""stk.v.11.7.0
BEGIN Attitude
NumberOfAttitudePoints %s
ScenarioEpoch          21 Jun 2022 16:00:00.00000
InterpolationOrder     1
CentralBody            Earth
CoordinateAxes         Inertial
AttitudeTimeQuatAngVels
%s
END Attitude"""


time1=time.time()

def add_header(header,a_fname,n_points):
    STK_fname="STK_"+a_fname.strip("txt")+"a" #create attitude file name
    with open(a_fname,'r') as f:
        data=f.read()
    with open(STK_fname,'w') as f:
        print(header%(n_points,data),file=f) #Print data and n points to new file
        f.close()

def num_points(filename):
    with open(filename,'r') as f:
        data=f.readlines()
        length=len(data)
        return length

def main(a_fname,n_points,header):
    add_header(header,a_fname,n_points)
    time2=time.time()
    runtime=abs(time2-time1)
    print("Runtime: "+str(runtime))

if __name__ == "__main__":
    main(a_fname,n_points,header)




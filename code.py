import matplotlib.pyplot as plt
import numpy as np


T_air = 25
density = 7800 #kg/m^3
c_p = 460 #J/kg*k spec heat of steel
k = 45 #W/mK, thermal conductivity of wall
h = 150 #W/(m^2 K), The convective heat transfer of air/steel surface
thickness = 0.050 #m
dt = 1 #seconds




def node_temps(nodes, time_period):
    '''Point 1 - plot temp of 10 nodal points on a graph, and print max temp of
    inner and outer nodes'''
    temps_next = np.zeros((time_period+1,nodes))
    temps_next[0,:] = 25
    time = np.arange(0, time_period+1, 1)
    T_fur = 25 +475*(1 - np.e**(-time/180))
    
    dx = thickness/nodes
    A = k/dx #alpha in notes
    B = (density*c_p*dx)/(2*dt) #beta in notes    
    
    constants = mat_constants(nodes, A, B)
    
    for j in range(0, time_period, 1): # j is the current time step
        temps_curr = current_temps(nodes, j, temps_next, T_fur, B)
        temps_next[j+1] = np.linalg.solve(constants, temps_curr)
   
    return temps_next, time




def mat_constants(n, A, B):
    '''initialises the constants matrix at size n'''
    
    x = np.zeros((n,n))
    
    x[0,0] = -A-B-h
    x[0,1] = A
    
    for i in range(1, n-1):
        x[i, i-1] = A
        x[i, i] = -2*(A+B)
        x[i, i+1] = A
        
    x[n-1, n-2] = A
    x[n-1, n-1] = -A-B-h
    return x





def current_temps(n, j, temps_next, T_fur, B):
    '''Gets the temperatures at the current time steps'''
    
    temps_curr = np.zeros(n) 
    
    temps_curr[0] = -B*temps_next[j,0] - (h*T_fur[j+1])
    
    for k in range(1, n-1): #k is the node number, j is the time step
        temps_curr[k] = -2*B*temps_next[j,k]
        
    temps_curr[n-1] = -B*temps_next[j,n-1] - (h*T_air)    
      
    return temps_curr




def plot_graphs(data, nodes, time, time_period):
    '''plotting graphs'''
        
    for column in range(0, nodes):
        T_data = data[:,column]
        plt.plot(time, T_data)
    
    plt.title(f'Temps over t={time_period}')
    plt.xlabel('Time (s)')
    plt.ylabel('Temp (degrees C)')
    plt.show()   
    
    return None




#point 1
'''
data, time = node_temps(10, 2000)
print(f'Max inner temp = {data[-1,0]}')
print(f'Max outer temp = {data[-1,9]}')
plot_graphs(data, 10, time, 2000)  #getting data from node temps function then handing it to plot function
data1, time = node_temps(10, 5000) 
plot_graphs(data1, 10, time, 5000)
'''


# point 2, supposed to overlay them but then can't see the difference so could 
# just do them 1 at a time and arrange them side by side on the report maybe?
'''
data1, time = node_temps(10, 5000)
#data2, time = node_temps(50, 5000)  
#data3, time = node_temps(100, 5000) 
plot_graphs(data1, 10, time, 5000)
#plot_graphs(data2, 50, time, 5000)
#plot_graphs(data3, 100, time, 5000)
'''


#point 3
'''
data, time = node_temps(50, 5000)
data = data[::100,:]
time = time[::100]
print(f'Max inner temp = {data[-1,0]}')
print(f'Max outer temp = {data[-1,49]}')
plot_graphs(data, 50, time, 5000)
'''


#point 4
'''
temp_grad = np.zeros((5001, 48))
data, time = node_temps(50, 5000)
for j in range(0, 5000):
    for node in range( 1, 49): 
        temp_grad[j, node-1] = (-data[j+1, node-1] + data[j+1, node+1])/(2*0.005) 

temp_grad = temp_grad[0:-2] #This line and next as gradient goes to inifinity elsewise
time = time[0:-2]  # try commenting out to see what i mean
for column in range(0, 48):
    T_data = temp_grad[:,column]
    plt.plot(time, T_data)

plt.title(f'Thermal gradient over t=5000')
plt.xlabel('Time (s)')
plt.ylabel('Thermal gradient')
plt.show()

print(f'max value is {np.max(temp_grad)}')
'''


# point 5
# change the value of h to 150,000
'''
data, time = node_temps(50, 5000)
data = data[::100,:]
time = time[::100]
print(f'Max inner temp = {data[-1,0]}')
print(f'Max outer temp = {data[-1,49]}')
plot_graphs(data, 50, time, 5000)
'''

def node_temps6(nodes, time_period, h_in, h_out):
    '''Point 1 - plot temp of 10 nodal points on a graph, and print max temp of
    inner and outer nodes'''
    temps_next = np.zeros((time_period+1,nodes))
    temps_next[0,:] = 25
    time = np.arange(0, time_period+1, 1)
    T_fur = 25 +475*(1 - np.e**(-time/180))
    
    dx = thickness/nodes
    A = k/dx #alpha in notes
    B = (density*c_p*dx)/(2*dt) #beta in notes    
    
    constants = mat_constants6(nodes, A, B, h_in, h_out)
    
    for j in range(0, time_period, 1): # j is the current time step
        temps_curr = current_temps6(nodes, j, temps_next, T_fur, B, h_in, h_out)
        temps_next[j+1] = np.linalg.solve(constants, temps_curr)
   
    return temps_next, time


def current_temps6(n, j, temps_next, T_fur, B, h_in, h_out):
    '''Gets the temperatures at the current time steps'''
    
    temps_curr = np.zeros(n) 
    
    temps_curr[0] = -B*temps_next[j,0] - (h_in*T_fur[j+1])
    
    for k in range(1, n-1): #k is the node number, j is the time step
        temps_curr[k] = -2*B*temps_next[j,k]
        
    temps_curr[n-1] = -B*temps_next[j,n-1] - (h_out*T_air)    
      
    return temps_curr



def mat_constants6(n, A, B, h_in, h_out):
    '''initialises the constants matrix at size n'''
    
    x = np.zeros((n,n))
    
    x[0,0] = -A-B-h_in
    x[0,1] = A
    
    for i in range(1, n-1):
        x[i, i-1] = A
        x[i, i] = -2*(A+B)
        x[i, i+1] = A
        
    x[n-1, n-2] = A
    x[n-1, n-1] = -A-B-h_out
    return x


# Part 6
'''
h_in =150
h_out = 1621   #some h value I determined by guess and check
data, time = node_temps6(50, 5000, h_in, h_out)
#plot_graphs(data, 50, time, 5000)


temp_grad = np.zeros((5001, 48))
for j in range(0, 5000):
    for node in range( 1, 49): 
        temp_grad[j, node-1] = (-data[j+1, node-1] + data[j+1, node+1])/(2*0.005) 

temp_grad = temp_grad[:-2, :] # same as checkpoint four, removing infinite part of grad
time = time[:-2]
        
plt.plot(time, temp_grad)        
plt.title(f'Thermal gradient over t=5000 with fan')
plt.xlabel('Time (s)')
plt.ylabel('Thermal gradient')
plt.show()

print(f'max temp inner node is {np.max(data)}')
print(f'max gradient is {np.max(temp_grad)}')
'''
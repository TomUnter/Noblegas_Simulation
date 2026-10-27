#Velocity Verlet + Temperature regulator + Cutoff List


import numpy as np
from scipy import constants
import time
#array with normaldistributed velocities for every atom
def get_velocities(mu, sigma, N):
#mu = 0
#sigma = 1
    v = np.zeros((N, 3))
    N = np.arange(0, N)
    for i in N:
        v[i] = np.random.normal(mu, sigma, 3)
   
    return v

#NxN array with distances between every atom
#Index as list to indicate atoms in a cutoff range
def distance(x, N):
    dist = np.zeros((N, N))
    Index = []
    for i in range(len(x)):
        Index.append([])
        for j in range(len(x)):
            dist[i,j] = np.sqrt((x[i,0]-x[j,0])**2+(x[i,1]-x[j,1])**2+(x[i,2]-x[j,2])**2)
            if dist[i, j] < 1 and dist[i, j] > 0:
                Index[i].append(j)
    return dist, Index

#LJ Force Vector is a Nx3 Vector with Lennard Jones Forces for every atom in x,y,z direction
#distance Matrix and its Index list are computed every x iterations since atoms do not move too significantly for small time intervals 
def LJ(dist, x, N, Index):
    eps = 1
    sig = 0.34
    A = 4*eps*sig**12
    B = 4*eps*sig**6
    Force_Vector_new = np.zeros((N, 3))
    for i in range(0, len(Index)):
        for k in Index[i]:
            Force_Vector_new[i] +=  (12*eps*(sig**12)/(dist[i,k]**13)-6*eps*(sig**6)/(dist[i,k]**7))*((x[i]-x[k])/np.sqrt(np.sum((x[i]-x[k])**2)))
            
    return Force_Vector_new

#Pressure calculation in respect to the equation in the lecture

def Pressure(dist, N, T):
    eps = 1
    sig = 0.34
    dV = np.zeros(())
    for i in range(len(dist)):
        for j in range(i+1, len(dist)):
            dV = np.sum((12*eps*(sig**12)/(dist[i,j]**13)-6*eps*(sig**6)/(dist[i,j]**7)))
    Pressure = (1/1000)*(N*constants.R*10**(-3)*T-1/3*dV)
    return Pressure

#Velocity Verlet differential equation using the lennard jones force vector.
#

def Verlet(x_new, x, Force_Vector, Force_Vector_new, deltaT, v_new):
    m = 39.94
    
    x_new2 = x_new + v_new*deltaT + (Force_Vector/2/m)*deltaT**2
    v_new = v_new + ((Force_Vector_new + Force_Vector)/2/m)*deltaT
    
    return x_new2, x_new, v_new

#Simple and not well implemented boundary conditions.
# The 10**12 is a fail safe if the speed of the particles exceeds for some reason the lennard jones part for repulsion
#But the lennard jones repulsion term is enough

def LJ_boundary(x_new, N):
    eps = 1
    sig = 0.34
    x_dist_low = np.zeros((N, 1))
    x_dist_up = np.zeros((N, 1))
    y_dist_low = np.zeros((N, 1))
    y_dist_up = np.zeros((N, 1))
    z_dist_low = np.zeros((N, 1))
    z_dist_up = np.zeros((N, 1))
    #boundary_dist = np.zeros((10, 3))
    Force_Vector_boundary = np.zeros((N, 3))
    for i in range(0, N):
        #boundary_dist[i] = 
        if (x_new[i, 0] > 10):
            Force_Vector_boundary[i, 0] = -10**12
        elif x_new[i, 0] < 0:
            Force_Vector_boundary[i, 0] = 10**12
        else:
            x_dist_low[i] = x_new[i, 0]-0
            x_dist_up[i] = np.absolute(10 - x_dist_low[i])
            Force_Vector_boundary[i, 0] = 12*eps*(sig**12/x_dist_low[i]**13)-12*eps*(sig**12/x_dist_up[i]**13)
            
        if (x_new[i, 1] > 10):
            Force_Vector_boundary[i, 1] = -10**12
        elif x_new[i, 1] < 0:
            Force_Vector_boundary[i, 1] = 10**12
        else:
            y_dist_low[i] = x_new[i, 1]-0
            y_dist_up[i] = np.absolute(10 - y_dist_low[i])
            Force_Vector_boundary[i, 1] = 12*eps*(sig**12/y_dist_low[i]**13)-12*eps*(sig**12/y_dist_up[i]**13)
            
        if (x_new[i, 2] > 10):
            Force_Vector_boundary[i, 2] = -10**12
        elif x_new[i, 2] < 0:
            Force_Vector_boundary[i, 2] = 10**12
        else:
            z_dist_low[i] = x_new[i, 2]-0
            z_dist_up[i] = np.absolute(10 - z_dist_low[i])
            Force_Vector_boundary[i, 2] = 12*eps*(sig**12/z_dist_low[i]**13)-12*eps*(sig**12/z_dist_up[i]**13)
        
        
        
    return Force_Vector_boundary 

#write_file writes the velocities in a txt- and the coordinates in a pdb.file
def write_file(FILE1, FILE2, coords, vels, n_iter):

    '''

    FILE1 = Coordinaten File (Endung .pdb)
    FILE2 = Velocity File (Endung .txt, o.ä.)
    coords = Nx3 Frame der akuellen Coordinaten
    vels = Nx3 Frame der aktuellen Velocities
    n_iter = Nummer der aktuellen Iteration
    '''

    round_ = 2

    FILE1.write('Model {}\n'.format(n_iter)) #Beginn eines Frames
    
    for i in range(len(coords)):
        #ATOM  #ID   #NAME     #RESN      #RESID          #X    #Y     #Z
        coord_str = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format('ATOM', int(i),
                                                                                                                        'PAR', ' ', 'ARG',
                                                                                                                        ' ', int(i), ' ',
                                                                                                                        np.round(coords[i][0], round_)*10,
                                                                                                                        np.round(coords[i][1], round_)*10,
                                                                                                                        np.round(coords[i][2], round_)*10,
                                                                                                                         1.0, 1.0, ' ', ' ')


        FILE1.write(coord_str) #Schreibe Koordinaten
        FILE2.write('{}\t{}\t{}\n'.format(vels[i][0], vels[i][1], vels[i][2])) #Schreibe Geschwindigkeiten

    FILE1.write('TER\nENDMDL\n') #Ende eines Files


#AtomSim is the Simulation

def AtomSim(N, delL, deltaT):
    import numpy as np
    FILE1 = open('coords.pdb', 'w')
    FILE2 = open('vels.txt', 'w')
     
##################################################################################################
#Initialisierung
##################################################################################################
    start = time.time()  
    #from random import sample
    x = np.zeros((N+20, 3)) #empty coordinate array
    Anzahl = np.arange(0, N+20) #range for for loop
    Gitter = np.arange(1, delL, 0.5).tolist() #Gitterpunkte auf denen Atome liegen dürfen start:0, ende: delL aus Angabe , 0.5=5A Abstand
    n = 3 #for x,y,z coordinates
    for i in Anzahl:  
        x[i,:] = np.random.choice(Gitter, n, replace=False) #Random koordinaten Nx3 for every direction and atom
    x_new = np.unique(x, axis=0) #erkenne duplicates
    random_indices = np.random.choice(len(x_new), size=N, replace=False) #Indizes ohne replacement um unique rows zu erhalten
    x = x_new[random_indices, :] # erstellung uniquer coordinates
    v_new = get_velocities(0,0.15, N) #erstellung der start velocities standard deviation is calculated beforehand for given T0
    write_file(FILE1, FILE2, x, v_new, 1) #schreibe startkoordinaten in file
##################################################################################################
#1. Run
##################################################################################################   
    m = 39.94
    #Force_Vector_boundary = np.zeros((10, 3))
    #Force_Vector = np.zeros((10, 3))
    dist, Index = distance(x, N) #calculate distance between all atoms
    #Force_Vector_boundary, x_new = LJ_boundary(x)
    Force_Vector = LJ(dist, x, N, Index)#force between all atoms
    Force_Vector_boundary = LJ_boundary(x, N) #boundary force to keep atoms in pseudo box
    Force_Vector = Force_Vector + Force_Vector_boundary #addition of both forces
    x_new1 = x + v_new*deltaT+(Force_Vector/2*m)*deltaT**2 #first calculation step of coordinates according to Verlet
    #T0 = np.sum(m*v_new**2)/3/N/(constants.R*10**(-3))
    T0 = 100 #Starttemperature in correlation to start velocity
    write_file(FILE1, FILE2, x_new1, v_new, 2) 
    count = 1
##################################################################################################
# Iteration
##################################################################################################
    for i in range(0,10000000): 
        if i % 20 == 0:
            dist, Index = distance(x_new1, N) #calculate distance between all atoms 
        Force_Vector_boundary = LJ_boundary(x_new1, N) 
        #print(Force_Vector_boundary)
        Force_Vector_new = LJ(dist, x_new1, N, Index)  #force between all atoms
        #print(Force_Vector)
        Force_Vector_new = Force_Vector_boundary + Force_Vector_new 
        x_new1, x, v_new = Verlet(x_new1, x, Force_Vector, Force_Vector_new, deltaT, v_new) #Velocity Verlet 
        Force_Vector = Force_Vector_new #Set current Force Vector to one step behind (necessary for Velocity Verlet)
        #T = np.sum(m*v_new**2)/(3*N*constants.R)
        T = np.sum(m*v_new**2)*(10**3)/(3*N*(constants.R)) #Berendson Temperature regulation
        #print(T, T0)
        v_new = np.sqrt(1+(deltaT/1)*((T0/T)-1))*v_new #new speed with Temperature regulation
        #T1 = np.sum(m*v_new**2)*(10**3)/(3*N*(constants.R))
        #print(T, T1, T0)
        #print(v_new)
        #P = Pressure(dist, N, T)
        count += 1
        
      
        if count == 1500: #write only every 1500 entries and adjust temperature for downregulation 
            write_file(FILE1, FILE2, x_new1, v_new, i+3)
            if T0 <= 300 and T0 > 30:
                T0 = T0 - 0.5
            print('Iteration: ', i)
            print(T,T0)
            count = 0
    end = time.time() 
    print('elapsed time: ', end - start)
    FILE1.close()
    FILE2.close()
AtomSim(150, 9, 10**(-3))       

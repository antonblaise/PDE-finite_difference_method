# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 03:08:08 2021

@author: OMEN
"""

import numpy as np
import pandas as pd


#%% Functions

def meanError(matrix1, matrix2):
    if np.shape(matrix1) == np.shape(matrix2):  # if both matrices have the same size
        error = 0
        count = 0
        for i in range(np.shape(matrix1)[0]):
            for j in range(np.shape(matrix1)[1]):
                if matrix1[i,j] != 0 and matrix2[i,j] != 0:  # Only when both values are NOT zero
                    error += abs(matrix1[i,j] - matrix2[i,j])
                    count += 1
        avgError = error/count
        return avgError
    else:
        print("\nBoth matrices must be equal in size.\n")
    
def analytical(x,t):
    uxt = np.exp( -t ) * np.cos(pi*( x-0.5 ))
    return uxt

def crankNicolson(u):
    
    for j in range(np.shape(u)[0]-1):
        
        A = np.zeros((np.shape(u)[1]-1, np.shape(u)[1]-1))
        
        b = np.zeros(np.shape(A)[1])
        
        coeff = [-r, 2+(2*r), -r]
        
        # Build matrix A
        for m in range(np.shape(A)[0]):
            n = m-1
            if m == 0:
                A[m,0] = coeff[1]
                A[m,1] = coeff[2]
            elif m == np.shape(A)[0]-1:
                A[m,m-1] = coeff[0]
                A[m,m] = coeff[1]
            else:
                for member in coeff:
                    A[m,n] = member
                    n += 1
    
        # Build matrix b
        
        for i in range(np.shape(u)[1]-2):
            if i != 0:
                b[i-1] = r*u[j,i-1] + (2-(2*r))*u[j,i] + r*u[j,i+1]
        
        # Ax = b
        y = np.linalg.solve(A, b)
        
        # register values of y into the consequent row of u
        
        
        for i in range(np.shape(u)[1]-1):
            if i != 0:
                u[j+1,i] = '{0:.5f}'.format(y[i-1])
    
    return u

#%% Display title

print("SSCE2393-36 NUMERICAL METHODS\n\nGroup Assignment\n")
print("\nPartial Differential Equations: Finite Difference Method\n")


#%% Define values

pi = np.pi

h = [0.2, 0.5, 0.1]
k = [0.25, 0.2, 0.1]

part = ['a', 'b', 'c']
    
max_x = 2
max_t = 2 
    


#%% Calculation

for count in range(3):
    
    #define these values according to current values of h and k
    r = k[count]/((pi*h[count])**2)
    
    x = np.linspace(0, 2, num=int(max_x/h[count])+1)
    t = np.linspace(0, 2, num=int(max_t/k[count])+1)
    
    print("\n("+str(part[count])+") h = "+str(h[count])+"; k = "+str(k[count])+"\n")
    
    #Initiate the BIG matrix u(i,j)
    u = np.zeros((int(max_t/k[count])+1, int(max_x/h[count])+1), dtype=float)

    # Insert values for when t = 0, make all values into 5 d.p.
    u[0,:] = [ '{0:.5f}'.format(np.cos( pi * (x[i]-0.5) )) for i in range(len(x))]
    
    init_u = u
    
    u = crankNicolson(u)
    
    real_u = np.zeros(np.shape(u))
    
    for j in range(np.size(real_u, axis=0)):
        for i in range(np.size(real_u, axis=1)-1):
            real_u[j,i] = '{0:.5f}'.format(analytical(x[i], t[j]))


#%% Error calculation

    error = meanError(u, real_u)
        
    print("\nError = "+str(error))


#%% Export calculated result
    
    df = pd.DataFrame(u)
    real_df = pd.DataFrame(real_u)
    
    # result_Q1a_calculated_Explicit.csv
    
    FILENAME = "result_Q1" + str(part[count])
    C_FILENAME = FILENAME + "_calculated_CrankNic.csv"
    A_FILENAME = FILENAME + "_analytical.csv"
    
    fileLabel = "Q1 "
    C_fileLabel = fileLabel + "("+str(part[count])+")"
    A_fileLabel = "Analytical " + C_fileLabel
    
    df.to_csv(C_FILENAME, index_label=(C_fileLabel))
    
    print("\nCalculation complete!\n\nCalculation result has been exported as "+C_FILENAME+" into the current directory.\n")
        
    real_df.to_csv(A_FILENAME, index_label=(A_fileLabel))
    
    print("\nAnalytical solution has been exported as "+A_FILENAME+" into the current directory.\n")

#%%      
# =============================================================================
#      PLOTTING
# =============================================================================
    
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    X, T = np.meshgrid(x, t)
    error = '{0:.5f}'.format(error) # Change to 5 d.p.
    
    # Colours
    myColor = ['blue', 'green', 'purple']
    
    fileName = "Q1"+str(part[count])
    title = "("+str(part[count])+") h = "+str(h[count])+"; k = "+str(k[count])

    print("\nVisualising data. Please wait ... \n")

    # Calculated
    plot1 = plt.figure(1)
    
    ax = plot1.add_subplot(111, projection="3d")
    
    # axis labels
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_zlabel("u(x,t)")
    
    C_fileName = fileName + "_calculated_CrankNic.png"
    C_title = title + " - Calculated (Crank Nicolson)"
    
    ax.text2D(0.15, 0.95, C_title+"\n     Error = "+str(error), transform=ax.transAxes)
    
    ax.plot_wireframe(X, T, np.array(u), color=myColor[count])
    
    plt.savefig(C_fileName)
    
    plt.show()
    plt.close()
    
    #%% Grid plot
    
    df = pd.DataFrame(u)

    fig, ax = plt.subplots(figsize=(20,15))
    
    gridTitle = "\n("+str(part[count])+") h = "+str(h[count])+"; k = "+str(k[count])+" - Calculated\n\n Crank Nicolson method\n"
    
    sns.heatmap(u, annot=u).set_title(gridTitle+" (Error = "+str(error)+")\n", fontsize=20)
    
    plt.savefig("GRID_"+C_fileName)
    
    plt.show()
    plt.close()
   
    #%%
    
    # Analytical
    plot2 = plt.figure(2)
    
    ax = plot2.add_subplot(111, projection="3d")
    
    # axis labels
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_zlabel("u(x,t)")
    
    A_fileName = fileName + "_analytical.png"
    A_title = title + " - Analytical"
    
    ax.text2D(0.15, 0.95, A_title+"\n ", transform=ax.transAxes)
    
    ax.plot_wireframe(X, T, np.array(real_u), color=myColor[count])
    
    plt.savefig(A_fileName)
    
    plt.show()
    plt.close()

    #%% Grid plot
    
    real_df = pd.DataFrame(real_u)

    fig, ax = plt.subplots(figsize=(20,15))
    
    gridTitle = "\n("+str(part[count])+") h = "+str(h[count])+"; k = "+str(k[count])+" - analytical\n"
    
    sns.heatmap(real_u, annot=real_u).set_title(gridTitle+"\n", fontsize=20)
    
    plt.savefig("GRID_"+A_fileName)
    
    plt.show()
    plt.close()

#%% Divider
    
    print("\nDone!\n")
    
    print('\n______________________________________________________________________________\n')


#%% Done plotting

print("\nAll plots are saved as PNG to the current directory.\n")

print("\n ~ ~ ~ T H E   E N D ~ ~ ~ \n")

print("\nTHANK YOU AND ALL THE BEST ! ! !\n")


while True:
    pass




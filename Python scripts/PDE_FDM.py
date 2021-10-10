# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 00:42:58 2021

@author: Antonius
"""

import numpy as np
import pandas as pd

#%% Preset

repeat = input() #take user input

pi = np.pi

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

while repeat.lower() == 'y':
    print("SSCE2393-36 NUMERICAL METHODS\n\nGroup Assignment\n")
    print("\nPartial Differential Equations: Finite Difference Method\n")
    
    
    
    
    #%% Choose a, b or c
    
    print("\nPlease choose an option below:\n")
    print("(a) h = 0.2, k = 0.25\n(b) h = 0.5, k = 0.2\n(c) h = 0.1, k = 0.1")
    option = str(input("> "))
    
    correctChoices = ['a','b','c']
    
    # fail safe
    while option.lower() not in correctChoices:
        option = str(input("> "))
        
    if option.lower() == 'a':
        h = 0.2
        k = 0.25
    elif option.lower() == 'b':
        h = 0.5
        k = 0.2
    elif option.lower() == 'c':
        h = 0.1
        k = 0.1
        
    
    
    #%% Define values
    
    r = k/((pi*h)**2)
    
    max_x = 2
    max_t = 2
    
    x = np.linspace(0, 2, num=int(max_x/h)+1)
    t = np.linspace(0, 2, num=int(max_t/k)+1)
    
    
    #%% Create and configure the BIG matrix u(i,j)
    
    u = np.zeros((int(max_t/k)+1, int(max_x/h)+1), dtype=float)
    
    # Insert values for when t = 0, make all values into 5 d.p.
    u[0,:] = [ '{0:.5f}'.format(np.cos( pi * (x[i]-0.5) )) for i in range(len(x))]
    
    
    #%% Calculations and interation (choose a method)
    
    print("\nPlease choose a method for calculation")
    print("\n1. Explicit method\n2. Crank Nicolson method\n")
    mtd = str(input("> "))
    
    # fail safe
    while mtd != '1' and mtd != '2':
        mtd = str(input("> "))
    
    if mtd == '1':
        for j in range(np.size(u, axis=0)-1):
            for i in range(np.size(u, axis=1)-1):
                if i != 0:
                    u[j+1,i] = '{0:.5f}'.format( r*u[j,i-1] + (1-2*r)*u[j,i] + r*u[j,i+1] )
    elif mtd == '2':
        u = crankNicolson(u)


    #%% Calculate the analytical solution and find the man error
    
    real_u = np.zeros(np.shape(u))
    
    for j in range(np.size(real_u, axis=0)):
        for i in range(np.size(real_u, axis=1)-1):
            real_u[j,i] = '{0:.5f}'.format(analytical(x[i], t[j]))

    
    #%% Get error
    
    error = meanError(u, real_u)
    
    print("\nError = "+str(error))
    
    #%% Export calculated result
    
    df = pd.DataFrame(u)
    real_df = pd.DataFrame(real_u)
    
    # result_Q1a_calculated_Explicit.csv
    
    FILENAME = "result_Q1" + option.lower()
    C_FILENAME = FILENAME + "_calculated_"
    A_FILENAME = FILENAME + "_analytical.csv"
    
    fileLabel = "Q1 "
    C_fileLabel = fileLabel + "("+option.lower()+")"
    A_fileLabel = "Analytical " + C_fileLabel
    
    if mtd == '1':
        C_FILENAME += "Explicit.csv"
    elif mtd == '2':
        C_FILENAME += "CrankNic.csv"
        
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
    
    # Colour
    if option.lower() == 'a':
        myColor = "blue"
    elif option.lower() == 'b':
        myColor = "green"
    elif option.lower() == 'c':
        myColor = "purple"
    
    # Initiate file name and title
    fileName = "Q1"+option.lower()
    title = "("+option.lower()+") h = "+str(h)+"; k = "+str(k)
    
    print("\nWhich one do you want to visualise?\n")
    print("1. Calculated result\n2. Analytical result\n3. Both\n4. None")
    
    choices = ['1', '2', '3', '4']
    toPlot = str(input("> "))
    
    # fail safe
    while toPlot not in choices:
        toPlot = str(input("> "))
    
    if toPlot != '4':
        print("\nVisualising data. Please wait ... \n")
    
    #%% ONE
    
    # title
    if toPlot == '1':  # Calculated
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        
        # axis labels
        ax.set_xlabel("x")
        ax.set_ylabel("t")
        ax.set_zlabel("u(x,t)")
        
        fileName += "_calculated"
        title += " - Calculated"
        
        if mtd == '1':
            title += " (Explicit)"
            fileName += "_Explicit.png"
        elif mtd == '2':
            title += " (Crank Nicolson)"
            fileName += "_CrankNic.png"
        
        ax.text2D(0.15, 0.95, title+"\n     Error = "+str(error), transform=ax.transAxes)
        
        ax.plot_wireframe(X, T, np.array(u), color=myColor)
        
        plt.savefig(fileName)
        
        plt.show()
        
        #%% Grid plot
        
        df = pd.DataFrame(u)
    
        fig, ax = plt.subplots(figsize=(20,15))
        
        gridTitle = "\n("+option.lower()+") h = "+str(h)+"; k = "+str(k)+" - Calculated\n\n"
        
        if mtd == '1':
            gridTitle += " Explicit method\n"
        elif mtd == '2':
            gridTitle += " Crank Nicolson method\n"
        
        sns.heatmap(u, annot=u).set_title(gridTitle+" (Error = "+str(error)+")\n", fontsize=20)
        
        plt.savefig("GRID_"+fileName)
        
        plt.show()
        
        #%% TWO
        
    elif toPlot == '2':  # Analytical
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        
        # axis labels
        ax.set_xlabel("x")
        ax.set_ylabel("t")
        ax.set_zlabel("u(x,t)")
        
        fileName += "_analytical.png"
        title += " - Analytical"
        
        ax.text2D(0.15, 0.95, title+"\n", transform=ax.transAxes)
        
        ax.plot_wireframe(X, T, np.array(real_u), color=myColor)
        
        plt.savefig(fileName)
        
        plt.show()
    
    #%% Grid plot
        
        real_df = pd.DataFrame(real_u)
    
        fig, ax = plt.subplots(figsize=(20,15))
        
        gridTitle = "\n("+option.lower()+") h = "+str(h)+"; k = "+str(k)+" - Analytical\n\n"
        
        if mtd == '1':
            gridTitle += " Explicit method\n"
        elif mtd == '2':
            gridTitle += " Crank Nicolson method\n"
        
        sns.heatmap(real_u, annot=real_u).set_title(gridTitle+"\n", fontsize=20)
        
        plt.savefig("GRID_"+fileName)
        
        plt.show()
        
    #%% THREE
    
    elif toPlot == '3':  # Both
        
        # Calculated
        plot1 = plt.figure(1)
        
        ax = plot1.add_subplot(111, projection="3d")
        
        # axis labels
        ax.set_xlabel("x")
        ax.set_ylabel("t")
        ax.set_zlabel("u(x,t)")
        
        C_fileName = fileName + "_calculated"
        C_title = title + " - Calculated"
        
        if mtd == '1':
            C_title += " (Explicit)"
            C_fileName += "_Explicit.png"
        elif mtd == '2':
            C_title += " (Crank Nicolson)"
            C_fileName += "_CrankNic.png"
        
        ax.text2D(0.15, 0.95, C_title+"\n     Error = "+str(error), transform=ax.transAxes)
        
        ax.plot_wireframe(X, T, np.array(u), color=myColor)
        
        plt.savefig(C_fileName)
        
        plt.show()
        
        #%% Grid plot
        
        df = pd.DataFrame(u)
    
        fig, ax = plt.subplots(figsize=(20,15))
        
        gridTitle = "\n("+option.lower()+") h = "+str(h)+"; k = "+str(k)+" - Calculated\n\n"
        
        if mtd == '1':
            gridTitle += " Explicit method\n"
        elif mtd == '2':
            gridTitle += " Crank Nicolson method\n"
        
        sns.heatmap(u, annot=u).set_title(gridTitle+" (Error = "+str(error)+")\n", fontsize=20)
        
        plt.savefig("GRID_"+C_fileName)
        
        plt.show()
       
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
        
        
        ax.text2D(0.15, 0.95, A_title+"\n", transform=ax.transAxes)
        
        ax.plot_wireframe(X, T, np.array(real_u), color=myColor)
        
        plt.savefig(A_fileName)
        
        plt.show()
    
        #%% Grid plot
        
        real_df = pd.DataFrame(real_u)
    
        fig, ax = plt.subplots(figsize=(20,15))
        
        gridTitle = "\n("+option.lower()+") h = "+str(h)+"; k = "+str(k)+" - analytical\n\n"
        
        sns.heatmap(real_u, annot=real_u).set_title(gridTitle+"\n", fontsize=20)
        
        plt.savefig("GRID_"+A_fileName)
        
        plt.show()
    
    #%% FOUR
    
    elif toPlot == '4':  # None
        pass
    
    if toPlot != '4':
        print("\nAll plots are saved as PNG to the current directory.\n")
        print("\n* * * NOTE * * * ")
        print("\nPlease be informed that the grid plot shows limited number of decimal places due to space constrains.\n")
        print("For full review of data in 5 decimal places, please use the CSV files generated.\nThank you for understanding!\n")
    
    #%% Ending
    
    print("\n ~ ~ ~ T H E   E N D ~ ~ ~ \n")
    print("\nWould you like to repeat? [y/n]\n")
    repeat = str(input("> "))
    
    # fail safe
    while repeat.lower() != 'y' and repeat.lower() != 'n':
        repeat = str(input("[y/n]> "))

    if repeat.lower() == 'y':
        print('\n______________________________________________________________________________\n\n\n')



#%% If user wants to stop and leave

print("\nTHANK YOU AND ALL THE BEST ! ! !\n")

while True:
    pass









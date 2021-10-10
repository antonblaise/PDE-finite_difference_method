
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

pi = np.pi

def analytical(x,t):
    uxt = np.exp( -t ) * np.cos(pi*( x-0.5 ))
    return uxt

resolution = 1000

x = np.linspace(0, 2, resolution)
t = np.linspace(0, 2, resolution)

full_u = np.zeros([np.size(x), np.size(t)])

X, T = np.meshgrid(x, t)
    
for j in range(np.size(full_u, axis=0)):
    for i in range(np.size(full_u, axis=1)-1):
        full_u[j,i] = analytical(x[i], t[j])
    
#%% Surface plot
        
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

# axis labels
ax.set_xlabel("x")
ax.set_ylabel("t")
ax.set_zlabel("u(x,t)")


ax.text2D(0.15, 0.95, "Actual surface\n", transform=ax.transAxes)

ax.plot_wireframe(X, T, np.array(full_u), color='red')

plt.show()

#%%

full_df = pd.DataFrame(full_u)
    
fig, ax = plt.subplots(figsize=(20,15))

sns.heatmap(full_u).set_title("Actual grid\n", fontsize=20)

plt.show()






















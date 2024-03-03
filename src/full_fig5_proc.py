
#PATH = "~/Git/A-Statistical-Approach-to-Quantum-Mechanics/GPU/"
#FILE_PATH = PATH + "data/data_fig9"
totalConfigurations = 10000
configurations_to_use = 20 
fig, ax = plt.subplots()
FILE_PATH = "data/data_fig5.csv"
ax.set_title(f"fig5: Probability Density")
Monte_Carlo_Iterations = 1360
radius = 4 
resolution = radius*0.05
tot_counts = np.zeros(np.arange(2*radius/resolution+1).shape[0]-1)
MCPerHistogram = 10
for i in range(int(Monte_Carlo_Iterations/MCPerHistogram)):
    data = pd.read_csv(FILE_PATH, sep=";", skiprows=np.delete(np.arange(Monte_Carlo_Iterations), np.arange(i*MCPerHistogram, (i+1)*MCPerHistogram, 1))).to_numpy()
    counts, bins = np.histogram(data.flatten(), np.arange(2*radius/resolution+1)*resolution-radius)
    tot_counts = tot_counts + counts
psi = tot_counts/tot_counts.sum()/(bins[1]-bins[0])
ax.scatter(bins[:-1]+resolution*0.5, psi, linewidths=1, marker="x")
x = np.linspace(-radius, radius, 100)
ax.plot(x, 3.14**(-1/2)*np.exp(-x**2), linestyle="--", color="r", linewidth=0.5)
ax.plot(x, 0.59*np.exp(-1.1*x**2), linestyle="--", color="b", linewidth=0.5)
fig.set_size_inches(11.69, 8.27)
fig.savefig("fig5.pdf", dpi=500)
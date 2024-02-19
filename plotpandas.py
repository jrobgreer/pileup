import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
#data = pd.read_csv("pileupcorrection50perc.csv")
data = pd.read_csv("calib140224.csv")
print(data)

data = data[data['Area']>0]
#data = data[data['Area']<1e6]
#data = data[data['Timestamp']<5000000]

# plt.scatter(data['Area'], data['Timestamp'], s=10, c=data['Chi2']/data['NDF'])

# plt.colorbar()
# plt.xlabel("Area")
# plt.ylabel("Timestamp")
# plt.grid()
# plt.minorticks_on()
# plt.show()

# plt.hist2d(data['Area'], data['Timestamp'], bins=(500,5000), norm=mpl.colors.LogNorm())
# plt.show()
# upper_range = 3e6


plt.close()

# us 
print(data['Area'])
plt.hist(data['Area'], bins=1000, label='us', histtype='step')

# for i in [1, 10, 100, 150, 200, 500, 5000]:
#     plt.hist(data['Area'][data['Timestamp']<i*1000], bins=np.linspace(0,upper_range, 100), label=f'<{i}us', histtype='step')

#plt.yscale('log')
#plt.hist(data['Area'][data['Timestamp']>500], bins=np.linspace(0,upper_range, 500), label='>500', histtype='step')
plt.xlabel("Pulse Area")
plt.ylabel("Counts")
plt.legend()

plt.show()
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
# data = pd.read_csv("pileupcorrection50perc.csv")
data = pd.read_csv("soilfit_EXTRADATA.csv")
print(data['EventID'].unique())

# data = data[data['Area'] > 0]


def guo_fit(x, par):
    '''par[0] - A
       par[1] - t0
       par[2] - theta1
       par[3] - theta2'''
    return par[0]*(np.exp(-(x-par[1])/par[2]) - np.exp(-(x-par[1])/par[3]))


pulse = np.zeros(1031)
for i in data['EventID'].unique():
    waveform_record = data[data['EventID'] == i]
    pulse = np.zeros(1031)
    for idx, row in waveform_record.iterrows():
        # print(idx)
        # print(row['Par0'])
        # print(pulse)

        pulse_to_add = np.array(guo_fit(np.linspace(1, 1030, 1031), [
                                row['Par0'], row['Par1'], row['Par2'], row['Par3']]))
        pulse_to_add[:int(row['Par1'])] = 0

        pulse = pulse+pulse_to_add

        print("Rem Pulse area: ", row['Rem Pulse Area'])
        print("Chi2: ", row['Chi2'])
        print("NDF: ", row['NDF'])
        print("Chi2/NDF: ", row['Chi2']/row['NDF'])
        print("-------------------------------")

        plt.plot(np.linspace(1, 1030, 1031), pulse)
        plt.ylim(0, 15000)
        plt.title("Waveform : {}".format(i))
        plt.show()
        plt.close()
    #          i for i in [waveform_record['Par0']]]))
    # plt.show()
    # plt.close()
# data = data[data['Area']<1e6]
# data = data[data['Timestamp']<5000000]

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
# print(data['Area'])
# plt.hist(data['Area'], bins=1000, label='us', histtype='step')

# for i in [1, 10, 100, 150, 200, 500, 5000]:
#     plt.hist(data['Area'][data['Timestamp']<i*1000], bins=np.linspace(0,upper_range, 100), label=f'<{i}us', histtype='step')

# plt.yscale('log')
# plt.hist(data['Area'][data['Timestamp']>500], bins=np.linspace(0,upper_range, 500), label='>500', histtype='step')
# plt.xlabel("Pulse Area")
# plt.ylabel("Counts")
# plt.legend()

# plt.show()

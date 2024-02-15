import spectranalysis as sp
import getwave as gw
from itertools import chain
import pandas as pd

pulse_collection, pulse_timestamps = gw.get_pulse_collection('soilfromoutside150224_15mins_GENON.dat', baseline=0.1)
# plt.close()
# plt.hist(pulse_timestamps)
# plt.show()
areas = []
timestamps = []
chi2 = []
ndf = []
print("EXPECT ", int(len(pulse_collection)), " PULSES")
input("Press to continue")

for pulse_idx in range(int(len(pulse_collection))):

    print("_______________________________________________________________________________")
    print(" PULSE INDEX ", pulse_idx)
    print("_______________________________________________________________________________")
    pulse = sp.Pulse(pulse_collection, pulse_timestamps, pulse_idx)

    #pulse.butter_lowpass_filtfilt(cutoff=15e6, fs=250e6, plotting=False) #25e6 was
    pulse.get_peaks2(min_dist_between_peaks=20, gradient_threshold=10)
    pulse.fit2()

    areas.append(pulse.areas)
    timestamps.append(pulse.true_timestamps)
    chi2.append(pulse.chi2)
    ndf.append(pulse.ndf)
    
    # Get the first load of peaks
    # try:
    #     pulse.get_peaks(plotting=False, threshold=500)
    # except:
    #     pulse.peak_heights = []
    #     pulse.peak_times = []
    #     print("No peaks found")
        
    
    # # As long as there is a peak to be fit to, fit 
    # while len(pulse.peak_heights) >=1:
    #     print("Peaks avaiable: ", len(pulse.peak_heights))
    #     pulse.fit(prefit_window=10, postfit_window=100)
    #     areas.append(pulse.area)
    #     timestamps.append(pulse.timestamp)

    
    #end = time.time()
    #elapsed = end - st

    #time_rem = (int(len(pulse_collection)) - pulse_idx) * elapsed / 60

    # print("Estimated time remaining : ", time_rem, " mins")
  

print(timestamps)
print(areas)
flat_timestamps = list(chain.from_iterable(timestamps))
flat_areas = list(chain.from_iterable(areas))
flat_chi2 = list(chain.from_iterable(chi2))
flat_ndf = list(chain.from_iterable(ndf))
print(flat_timestamps)
print(flat_areas)

# plt.close()
# plt.hist(np.array(flat_timestamps), bins=1000)
# plt.xlabel("Timestamps")
# plt.show()

# plt.close()
# plt.hist(np.array(flat_areas), bins=1000)
# plt.xlabel("Areas")
# plt.show()


# plt.close()
# plt.scatter(np.array(flat_areas), np.array(flat_timestamps))
# plt.xlabel("Areas")
# plt.ylabel("Timestamps")
# plt.show()



df = pd.DataFrame({'Area': np.array(flat_areas), 'Timestamp': np.array(flat_timestamps), 'Chi2': np.array(flat_chi2), 'NDF': np.array(flat_ndf)})
# Save to CSV file
csv_filename = 'soil_GeneratorON_15mins.csv'
df.to_csv(csv_filename, index=False)

# plt.close()
# plt.hist(areas, bins=1000)
# plt.show()
    

    # for idx, ph in enumerate(pulse.peak_heights):
    #     pulse.fit(prefit_window=10, postfit_window=100)
    #     if len(pulse.peak_heights)==0:
    #         break
    #     ##pulse.butter_lowpass_filtfilt(cutoff=1e6, fs=250e6, plotting=True)
    #     else:
    #         pulse.get_peaks(plotting=True, threshold=200)
    
    # plt.grid()
    # plt.minorticks_on()
    # plt.legend()
    # plt.show()



# Code needs to fit to tails of preamp pulses and deconvolve pileup 
# Should be able to timestamp events within a record


import numpy as np 
import ROOT
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter
from scipy.signal import butter, filtfilt
from scipy.fft import fft, fftfreq
from scipy.signal import argrelextrema
import getwave as gw
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import ROOT   
import array
import gc
import pandas as pd
import time
# COMPASS cannot timestamp events
# Changed setup to manually timestamp

# For neutron generator setup, we have an external trigger based upon when the neutron pulse occurs
# The TTL pulse from the generator goes into a gate generator
# The gate produced is sent into TRG-IN on DT5730
# Then set up 

class Pulse: 
    
    def __init__(self, pulse_collection, pulse_timestamps, index):
        self.record = pulse_collection[index]
        # Currently a crude timestamp based on the timestamp of the first trigger - need to add timestamping for events within a record - TTT + the offset of the extra pulse
        self.timestamp = pulse_timestamps[index]
        self.time = np.linspace(0,1029,1030)
        self.record_length = len(pulse_collection[index])


    def fft(self, plotting=False):
        '''Produce FFT spectrum - optional plotting'''
        self.yf = fft(pulse)
        self.xf = fftfreq(self.record_length, 4e-9) # 4ns between samples, should probably make this more general 
        
        if plotting == True:
            plt.plot(self.xf, np.abs(self.yf))
            plt.show()
            plt.close()

    def butter_lowpass_filtfilt(self, cutoff, fs, order=5, plotting=False):
        '''Low pass filter'''
        
        def butter_lowpass(cutoff, fs, order=5):
            nyq = 0.5 * fs
            normal_cutoff = cutoff / nyq
            b, a = butter(order, normal_cutoff, btype='low', analog=False)
            return b, a
            
        b, a = butter_lowpass(cutoff, fs, order=order)
        y = filtfilt(b, a, self.record)

        if plotting == True:
            # Plot pre and post filter waveforms
            #plt.plot(self.record, label='Original')
            self.record = y
            plt.plot(self.record, label=f'Filtered @ {cutoff} Hz')
            #plt.legend()
            # plt.show()
            # plt.close() 

        self.record = y
    
    def plot(self):
        plt.plot(self.record, label='Original')
        #plt.show()

    
    def get_peaks(self, plotting=False, interpolation_length=10000, minimum_distance=10, threshold=500, diff_scaling=10, minimum_time=0):
        '''Differentiate and locate peaks based on gradient, second differential, and minimum distance apart
           Disgusting code, needs tidying, but it works most of the time'''

        # Calculate gradient
        grad = []

        for i in range(len(self.record)-1):

            # Factor of diff_scaling used for visibility in plot
            # 8 ns per sample
            grad.append(diff_scaling*(self.record[i+1]-self.record[i])/8)

        grad = np.array(grad)

        # Calculate gradient of gradient
        gradgrad = []
        
        for i in range(len(grad)-1):
            gradgrad.append(diff_scaling*(grad[i+1]-grad[i])/8)

        gradgrad = np.array(gradgrad)

        # Interpolate the gradient, to get values closer to zero for steep changes - show maxima
        fine_grad = interp1d(np.linspace(0,1029,1029), grad)
        fine_gradgrad = interp1d(np.linspace(0,1029,1028), gradgrad)
        pulse_interp = interp1d(np.linspace(0,1030, 1030), self.record)
        time_interp = interp1d(np.linspace(0,1030,1030), self.time)

        # Create arrays
        fine_grad = fine_grad(np.linspace(0,1029, interpolation_length))
        fine_gradgrad = fine_gradgrad(np.linspace(0,1029, interpolation_length))
        pulse_interp = pulse_interp(np.linspace(0,1029, interpolation_length))
        time_interp = time_interp(np.linspace(0,1029, interpolation_length))
        #print(time_interp)

        print("Making event mask")

        # Check the gradient is negative and small - ie we have reached a peak, and heading back down after
        # Also check if the pulse value is above threshold
        # Check second diff is negative, we are expecting gradient to get more negative
        #print(time_interp)
        #print(minimum_time)
        #print(time_interp[time_interp>minimum_time])
        
        # plt.close()
        # plt.axvline(x=minimum_time, label='Cutoff Time', linestyle='-')
        # plt.plot(self.record, label='Original')
        # plt.plot(grad, label='Gradient')
        # plt.plot(gradgrad, label='Double Gradient')
        # ##plt.scatter(peak_times, peak_heights, color='black', marker='*', s=100, label='Peaks')
        # plt.legend()
        # plt.show()
        # plt.close()

        event_mask = (fine_grad<0) & (fine_grad>-500) & (pulse_interp>threshold) & (fine_gradgrad > -100) & (time_interp > minimum_time)
        
        #plt.scatter(np.linspace(0,1029, interpolation_length)[event_mask], pulse_interp[event_mask], label='FULLMASK', marker='+', color='r', s=100)

        # Take first instance in found peaks to be THE peak, can then be subtracted from rest of them
        peak_times = []
        peak_heights = []

        # Take first peak as a peak
        #print(np.linspace(0,1029, interpolation_length)[event_mask])
        
        peak_t = np.linspace(0,1029, interpolation_length)[event_mask][0]
        peak_v = pulse_interp[event_mask][0]

        peak_times.append(peak_t)
        peak_heights.append(peak_v)
      
        for peak_idx in range(len(np.linspace(0,1029, interpolation_length)[event_mask])):
            
            if np.linspace(0,1029, interpolation_length)[event_mask][peak_idx]>peak_t+minimum_distance:
                peak_t = np.linspace(0,1029, interpolation_length)[event_mask][peak_idx]
                peak_times.append(peak_t)

                peak_v = pulse_interp[event_mask][peak_idx]
                peak_heights.append(peak_v)

            else:
                continue

        if plotting==True:
            plt.axvline(x=minimum_time, label='Cutoff Time', linestyle='-')
            plt.plot(self.record, label='Original')
            plt.plot(grad, label='Gradient')
            plt.plot(gradgrad, label='Double Gradient')
            plt.scatter(peak_times, peak_heights, color='black', marker='*', s=100, label='Peaks')
            plt.legend()
            plt.show()
            plt.close()

        self.peak_times = peak_times
        self.peak_heights = peak_heights

    
    def fit(self, prefit_window=200, postfit_window=500):

        if len(self.peak_times) >= 1:
      

            # # Convert numpy arrays to PyROOT arrays
            x_array = self.time[self.time>self.peak_times[0]-prefit_window]
            y_array = self.record[self.time>self.peak_times[0]-prefit_window]

            # plt.close()
            # plt.plot(x_array, y_array)
            # plt.show()

            # Create a TGraph object
            graph = ROOT.TGraph(len(self.time[self.time>self.peak_times[0]-prefit_window]), x_array, y_array)

            def guo_fit(x,par):
                '''par[0] - A
                   par[1] - t0
                   par[2] - theta1
                   par[3] - theta2'''
                return par[0]*(np.exp(-(x[0]-par[1])/par[2]) - np.exp(-(x[0]-par[1])/par[3]))
            


            # Create a TF1 object for the fit function
            #fit_function = ROOT.TF1("landau", 'landau', self.peak_times[0]-prefit_window, self.peak_times[0]+postfit_window)
            fit_function = ROOT.TF1('guo_fit', guo_fit, self.peak_times[0]-prefit_window, self.peak_times[0]+postfit_window, 4)
            fit_function.SetParameters(2*self.peak_heights[0], self.peak_times[0], 50, 20)

            #Perform the fit
            fit_result = graph.Fit(fit_function, "RS")
            print(fit_result.IsValid())
            if fit_result.IsValid() ==False:
                print("FIT FAILURE")

                canvas = ROOT.TCanvas("canvas", "Guo Fit", 800, 600)
                graph.Draw("AP")
                
                fit_function.Draw("same")

                canvas.Update()
                canvas.Draw()
                
                # check if at end of pulse and fitting not possible in some way, if this is the case,
                # probably not necessary to throw away whole pulse

                input("Press Enter to close the canvas...")
            else:
                print("FIT SUCCESSFUL")

            #fitted_pulse = np.array([fit_function.GetParameter(0) * ROOT.TMath.Landau(t, fit_function.GetParameter(1), fit_function.GetParameter(2)) for t in self.time[self.time>self.peak_times[0]-prefit_window]])
            fitted_pulse = np.array([guo_fit([t], [fit_function.GetParameter(i) for i in range(fit_function.GetNpar())]) for t in self.time[self.time>self.peak_times[0]-prefit_window]])

            area = np.sum(fitted_pulse)
            self.area= area

            # EDIT THIS TO CORRECT THE TIMESTAMP
            self.timestamp = self.timestamp
            print("AREA: ", area)

            # plt.close()
            # plt.plot(self.time, self.record, label='Original full pulse')
            # self.record[self.time>self.peak_times[0]-prefit_window] = self.record[self.time>self.peak_times[0]-prefit_window] - fitted_pulse
            # plt.plot(self.time[self.time>self.peak_times[0]], self.record[self.time>self.peak_times[0]], label='Original pulse - fit region', linestyle='-')
            # plt.plot(self.time[self.time>self.peak_times[0]-prefit_window], fitted_pulse, label='Fit pulse')
            # plt.plot(self.time, self.record, label='Subtracted fit pulse')

            # plt.legend()
            # plt.show()

            # plt.close()

            # Arbitrary offset of 50 now before another peak is allowed
            try:
                self.get_peaks(minimum_time=self.peak_times[0]+50, threshold=500)
            except:
                self.peak_heights = []
                self.peak_times = []
                print("No peaks found")


            # Now find area of this fitted exponential from the maximum onwards 

            # Then look at doing 2 pulse situations where we need to start
            # iterative subtraction

            # Fit to first, find area, subtract from full waveform
            # Fit to next, find area, subtract from full waveform
            # Fit to next ... etc.
            # In the end the 'waveform' should be flat and we have extracted
            # N areas where N is number of peaks found
        else:
            print("Insufficient peaks")
            






        
    


    



pulse_collection, pulse_timestamps = gw.get_pulse_collection('testdataGENERATOR.dat')


areas = []
timestamps = []
print("EXPECT ", int(len(pulse_collection)), " PULSES")
input("Press to ocntinue")

for pulse_idx in range(int(0.03*len(pulse_collection))):

    
    st = time.time()

    print("_______________________________________________________________________________")
    print(" PULSE INDEX ", pulse_idx)
    print("_______________________________________________________________________________")
    pulse = Pulse(pulse_collection, pulse_timestamps, pulse_idx)
    #pulse.plot()
    
    pulse.butter_lowpass_filtfilt(cutoff=22e6, fs=250e6, plotting=False) #25e6 was
    
    # Get the first load of peaks
    try:
        pulse.get_peaks(plotting=False, threshold=500)
    except:
        pulse.peak_heights = []
        pulse.peak_times = []
        print("No peaks found")
        
    
    # As long as there is a peak to be fit to, fit 
    while len(pulse.peak_heights) >=1:
        print("Peaks avaiable: ", len(pulse.peak_heights))
        pulse.fit(prefit_window=10, postfit_window=100)
        areas.append(pulse.area)
        timestamps.append(pulse.timestamp)

    
    end = time.time()
    elapsed = end - st

    time_rem = (int(0.02*len(pulse_collection)) - pulse_idx) * elapsed / 60

    print("Estimated time remaining : ", time_rem, " mins")
    #input("Look look look")



df = pd.DataFrame({'Area': areas, 'Timestamp': timestamps})
# Save to CSV file
csv_filename = 'pileupcorrection.csv'
df.to_csv(csv_filename, index=False)

plt.close()
plt.hist(areas, bins=1000)
plt.show()
    

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


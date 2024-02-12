from gimmedatwave.gimmedatwave import gimmedatwave as gdw
import matplotlib.pyplot as plt
import numpy as np

#file = 'testdataGENERATOR.dat'



def get_pulse_collection(file, digitizer_family = gdw.DigitizerFamily.X725):
    parser = gdw.Parser(file, digitizer_family=digitizer_family)
    print(parser.n_entries)

    pulse_collection = np.empty((parser.n_entries, 1030))
    pulse_timestamps = np.empty(parser.n_entries)
    print(pulse_collection)
    for i in range(parser.n_entries):
        event = parser.get_event(i)
        pulse_collection[i] = -1*(0.1*16384 + event.record - 16384)
        pulse_timestamps[i] = event.header.trigger_time_tag 
        #plt.plot(8e-9*(np.linspace(0, 1030*4, 1030)+event.header.trigger_time_tag)/1e-6, -1*(-16384+(event.record)+0.1*16384))
        #times.append(8e-9*event.header.trigger_time_tag)

    #plt.show()  
    return pulse_collection, pulse_timestamps



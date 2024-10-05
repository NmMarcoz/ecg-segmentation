##PORTE DO CÓDIGO DE SEGMENTAÇÃO ECG
##AUTOR DO ALGORTIMO: JONATHAN ARAUJO QUEIROZ
#AUTOR DO PORTE: MARCOS RAFAEL NOGUEIRA MOREIRA

import pandas as pd
import numpy as np
import detectionPeaks as dt
import wfdb as wf
#import ecg_plot as ecg
import matplotlib.pyplot as plt

#record = wf.rdsamp('16265')
#df =pd.DataFrame(record[0], columns=record[1]['sig_name'])
#df.to_csv('16265.csv')
df_csv = pd.read_csv('./data/16265.csv')

signal = df_csv._get_label_or_level_values('ECG1').tolist()
signalY = df_csv._get_label_or_level_values('ECG2').tolist()

signalXY = [signal, signalY]
#signal = signal[5000:100000]
fs = round(81*0.7,2)
#fs = round(20*0.2, 2)
print(f"fs: {fs}")
theta = 0.4
lbda = 0.6
lbdaP = 0.15
thetaQRS = 0.15
lbdaQRS = 0.15
#---------------------pré processamento-----------------------#
def min_max_scaling(signal):
    """Normalize a signal using min-max scaling."""
    return (signal - signal.min()) / (signal.max() - signal.min())
def signal_extract(signal, window_size_percentage):
    window_size  = int(len(signal)*window_size_percentage)
    start_index = int(len(signal)/2) - window_size//2
    end_index = start_index + window_size
    return signal[start_index:end_index]

#RETIRAR A MÉDIA
def signal_mean(signal):
   signal_mean = pd.array(signal)
   signal_mean = signal_mean - signal_mean.mean()
   return signal_mean


def signal_std(signal):
    signalPd = pd.Series(signal)
    signal_std = signalPd.std()
    normalized_signal = signalPd / signal_std
    return normalized_signal

windowedSignal = signal_extract(signal, 0.01)
signal_mean = signal_mean(windowedSignal)
signal_std = signal_std(signal_mean)
final_signal = signal_std
print(final_signal)

peaks = dt.dtPeaks(final_signal, [0,60], fs, 0 )

qrs_amplitude = peaks[1]
qrs_index = peaks[2]
delay = peaks[3]
peaks_array = pd.array(peaks[0], int)

print(f"qrsIndex: {qrs_index}")
print(f"qrsamplitude: {qrs_amplitude}")
#print(qrs_amplitude)

#-------------------------CICLO-----------------------------------#
# if(qrs_index[1] < theta*fs):
#     qrs_index[1] = []
#     qrs_amplitude[1] = []


final_signal = final_signal[0:20]
# signal = signal[0:175]
# peaks[0] = peaks[0][0:175]

timeAxis = np.arange(len(peaks[0])) / fs
timeAxisNormalSignal = np.arange(len(signal))/fs

# plt.figure(figsize=(10, 6))
# plt.plot(timeAxis, peaks[0])
#     #plt.plot(timeAxisNormalSignal, signal)
# plt.xlabel('Tempo (s)')
# plt.ylabel('Batimentos')
# plt.title('ECG Signal')
# plt.grid(True)  # Add a grid for better readability
# plt.show()

def plotSignal(signal, title):
    signal = signal[0:500]
    timeAxis = np.arange(len(signal)) / fs
    plt.figure(figsize=(10, 6))
    plt.plot(timeAxis, signal)
        # plt.plot(timeAxisNormalSignal, signal)
    plt.xlabel('Tempo (s)')
    plt.ylabel('Batimentos')
    plt.title(title)
    plt.grid(True)  # Add a grid for better readability
    plt.show()

plotSignal(peaks[0], "com detecção de pico")
plotSignal(signal, "sem detecção de pico")
plotSignal(qrs_amplitude, "qrs amplitude")



##PORTE DO CÓDIGO DE SEGMENTAÇÃO ECG
##AUTOR DO ALGORTIMO: JONATHAN ARAUJO QUEIROZ
#AUTOR DO PORTE: MARCOS RAFAEL NOGUEIRA MOREIRA

import pandas as pd
import wfdb as wf
import pygwalker as pyg


# record = wf.rdsamp('16265')
# df =pd.DataFrame(record[0], columns=record[1]['sig_name'])
# df.to_csv('16265.csv')
df_csv = pd.read_csv('16265.csv')

#print(df_csv._get_label_or_level_values('ECG1'))
#-----------------------PARAMETROS------------------------#
df_array = df_csv._get_label_or_level_values('ECG1').tolist()
print(type(df_array))
df_array = df_array[5000:100000]
fs = round(81*0.7,2)
theta = 0.4
lbda = 0.6
lbdaP = 0.15
thetaQRS = 0.15
lbdaQRS = 0.15







# def signal_extract_pd(signal, window_size_percentage):
#     signal = pd.read_csv('16265.csv')
#     window_size = int(signal.)
#

#---------------------pré processamento-----------------------#
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

def detect_peaks(signal, threshold):
    signalPd = pd.Series(signal)
    rolling_max = signalPd.rolling(window = 2).max()
    peaks = (signalPd >  rolling_max) & (signalPd > threshold)
    return peaks

windowedSignal = signal_extract(df_array, 0.1)
signal_mean = signal_mean(windowedSignal)
signal_std = signal_std(signal_mean)
final_signal = signal_std
print(final_signal)



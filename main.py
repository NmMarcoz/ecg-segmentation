##PORTE DO CÓDIGO DE SEGMENTAÇÃO ECG
##AUTOR DO ALGORTIMO: JONATHAN ARAUJO QUEIROZ
#AUTOR DO PORTE: MARCOS RAFAEL NOGUEIRA MOREIRA
import sys
import os
import pandas as pd
import numpy as np
from matplotlib.backends.backend_pgf import PdfPages
from datetime import datetime
from scipy import signal as sg
import time  # Importa o módulo time

import detectionPeaks as dt
#import wfdb as wf
#import ecg_plot as ecg
import matplotlib.pyplot as plt
import matplotlib as mt

#plt.rcParams['pgf.texsystem'] = 'xelatex'

#os.environ["PATH"] += os.pathsep + '/usr/local/texlive/2018/bin/x86_64-darwin'
print("recursion limite", sys.getrecursionlimit())

# record = wf.rdsamp('16265')
# df =pd.DataFrame(record[0], columns=record[1]['sig_name'])
# df.to_csv('16265.csv')RuntimeError: 'xelatex' not found; install it or change rcParams['pgf.texsystem'] to an available TeX implementation

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

start_index = 0
end_index = 0

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



def bandpass_filter(signal, fs, lowcut=0.5, highcut=50.0, order=2):
    """
    Aplica um filtro passa-banda Butterworth ao sinal
    
    Parameters:
    signal: array-like - Sinal de entrada
    fs: float - Frequência de amostragem
    lowcut: float - Frequência de corte inferior (Hz)
    highcut: float - Frequência de corte superior (Hz)
    order: int - Ordem do filtro
    """
    nyquist = 0.5 * fs
    # Garante que as frequências normalizadas estejam entre 0 e 1
    low = min(lowcut / nyquist, 0.99)
    high = min(highcut / nyquist, 0.99)
    
    # Verifica se as frequências são válidas
    if low >= high:
        raise ValueError("lowcut deve ser menor que highcut")
    if low <= 0:
        low = 0.01
        
    b, a = sg.butter(order, [low, high], btype='band')
    return sg.filtfilt(b, a, signal)

windowedSignal = signal_extract(signal, 0.05)
filtered_signal = bandpass_filter(windowedSignal, fs, lowcut=0.5, highcut=40.0, order=4)
signal_mean = signal_mean(filtered_signal)
signal_std = signal_std(signal_mean)
final_signal = signal_std
#print(final_signal)

# Normalize both signals before SNR calculation
originalSignal = min_max_scaling(np.array(signal[0:len(final_signal)]))
filtredSignal = min_max_scaling(np.array(final_signal))
noise = originalSignal - filtredSignal
snr = 10 * np.log10(np.sum(filtredSignal**2) / np.sum(noise**2))

signal_power = np.sum(filtredSignal**2)
noise_power = np.sum(noise**2)
print(f"Signal power: {signal_power}")
print(f"Noise power: {noise_power}")
print(f"Power ratio: {signal_power/noise_power}")
print(f"SNR (dB): {snr}")

print("processando...")
start_time = time.time()  # Marca o início do tempo
peaks = dt.dtPeaks(final_signal, [0,60], fs, 0)
end_time = time.time()  # Marca o fim do tempo
print("processado!")

# Calcula o tempo de execução
execution_time = end_time - start_time
print(f"Tempo de execução da função dtPeaks: {execution_time:.4f} segundos")

ecgH = peaks[0]
qrs_amplitude = peaks[1]
qrs_index = peaks[2]
delay = peaks[3]
peaks_array = pd.array(peaks[0], int)

#print(f"qrsIndex: {qrs_index}")
#print(f"qrsamplitude: {qrs_amplitude}")
#print(qrs_amplitude)

#TODO -> tenho que fazer essa parte do ciclo, e da segmentação.
#-------------------------CICLO-----------------------------------#
if qrs_index[1] < theta  *fs:
    qrs_index[1] = []
    qrs_amplitude[1] = []

if qrs_index[-1] > len(signal) - (0.6*fs):
    qrs_index[-1] = []
    qrs_amplitude[-1] = []
# Batimentos
B = np.zeros((len(qrs_amplitude), int(fs + 1)))

P = np.zeros((len(qrs_amplitude), round(lbdaP * fs) + 1))

QRS = np.zeros((len(qrs_amplitude), round(thetaQRS * fs) + round(lbdaQRS * fs) + 1))

T = np.zeros((len(qrs_amplitude), round(0.3 * fs) + 1))

#print("qrs index", qrs_index[0])

#segmentos
for i in range(len(qrs_amplitude)):
    #batimentos
    #print("qrs index", qrs_index[i])
    B[i] = signal[qrs_index[i] - round(theta*fs):qrs_index[i] + round(lbda*fs)]

    start_P = qrs_index[i] - round(theta * fs)
    end_P = start_P + round(lbdaP * fs)
    P[i] = signal[start_P:end_P + 1]

    start_QRS = qrs_index[i] - round(thetaQRS * fs)
    end_QRS = qrs_index[i] + round(lbdaQRS * fs)
    QRS[i] = signal[start_QRS:end_QRS + 1]

    # Onda T
    start_T = (qrs_index[i] + round(lbda * fs)) - round(0.3 * fs)
    end_T = qrs_index[i] + round(lbda * fs)
    T[i] = signal[start_T:end_T + 1]
#---------------------------------------------------------------#
B = B.T
P = P.T
QRS = QRS.T
T = T.T

# signal = signal[0:175]
# peaks[0] = peaks[0][0:175]

timeAxis = np.arange(len(peaks[0])) / fs
timeAxisNormalSignal = np.arange(len(signal))/fs


def saveToPng(figs):
    date = datetime.today().isoformat(timespec='seconds')

    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    axs = axs.ravel()  
    
    # Plot each figure in its own subplot
    for i, fig_item in enumerate(figs):
        if isinstance(fig_item, plt.Figure):
            data = fig_item.axes[0].lines[0].get_data()
            axs[i].plot(data[0], data[1])
            axs[i].set_title(fig_item.axes[0].get_title())
            axs[i].set_xlabel(fig_item.axes[0].get_xlabel())
            axs[i].set_ylabel(fig_item.axes[0].get_ylabel())
            axs[i].grid(True)
            fig_item.savefig(f'./graficos/individual_fig_{i}.png', dpi = 300, bbox_inches='tight')
            plt.close(fig_item)  
    
    plt.tight_layout()
    

    plt.savefig(f'./graficos/ecg_analysis_{date}.png', dpi=300, bbox_inches='tight')
    plt.close()


def plotSignal(signal, title, reduce = 0, invert = False):

    timeAxis = []
    if(reduce != 0):
        signal = signal[0:reduce]
        timeAxis = np.arange(reduce) / fs
    else:
        timeAxis = np.arange(len(signal)) / fs
    plot = plt.figure(figsize=(10, 6))
    plt.subplot(2,1,1)
    # plt.plot(timeAxis, signal)
    if(invert):
        plt.plot(signal, timeAxis)
    else:
        plt.plot(timeAxis,signal)
    # plt.plot(timeAxisNormalSignal, signal)
    plt.xlabel('Tempo (s)')
    plt.ylabel('Sinal')
    plt.title(title)
    plt.grid(True)  # Add a grid for better readability
    #plt.show()
    return plot

#
#
signal = plotSignal(signal, "SINAL (ORIGINAL)", 500)
final_signal = plotSignal(final_signal, "SINAL (FILTRADO)", 500)
ecgProcessado = plotSignal(ecgH, "SINAL (FILTRADO E PROCESSADO)", 500)
#batimentos = plotSignal(B[0], "BATIMENTOS", 300)
complexoQrs = plotSignal(QRS, "COMPLEXO QRS")
ondaT = plotSignal(T, "ONDA T")
ondaP = plotSignal(P, "ONDA P")


saveToPng([signal, final_signal, ecgProcessado, complexoQrs, ondaT, ondaP])
#plot_segments(B, P, QRS, T, fs)

def plot_ecg_segments(B, P, QRS, T, fs):
    nf = 16
    fig = plt.figure(figsize=(15, 10))
    
    # Plot dos batimentos
    ax1 = plt.subplot(2,1,1)
    for i in range(B.shape[1]):  # Plot cada batimento
        ax1.plot(np.arange(len(B[:,i]))/fs, B[:,i], 'b-', alpha=0.1)
    ax1.plot(np.arange(len(B[:,0]))/fs, np.mean(B, axis=1), 'r-', linewidth=2)  # Média
    ax1.set_xlabel('Time (s)', fontsize=nf)
    ax1.set_ylabel('Amplitude (mV)', fontsize=nf)
    ax1.set_title('Heartbeats', fontsize=nf)
    ax1.grid(True)
    
    # Plot da onda P
    ax2 = plt.subplot(2,3,4)
    for i in range(P.shape[1]):
        ax2.plot(np.arange(len(P[:,i]))/fs, P[:,i], 'b-', alpha=0.1)
    ax2.plot(np.arange(len(P[:,0]))/fs, np.mean(P, axis=1), 'r-', linewidth=2)
    ax2.set_xlabel('Time (s)', fontsize=nf)
    ax2.set_ylabel('Amplitude (mV)', fontsize=nf)
    ax2.set_title('P Wave', fontsize=nf)
    ax2.grid(True)
    
    # Plot do complexo QRS
    ax3 = plt.subplot(2,3,5)
    for i in range(QRS.shape[1]):
        ax3.plot(np.arange(len(QRS[:,i]))/fs, QRS[:,i], 'b-', alpha=0.1)
    ax3.plot(np.arange(len(QRS[:,0]))/fs, np.mean(QRS, axis=1), 'r-', linewidth=2)
    ax3.set_xlabel('Time (s)', fontsize=nf)
    ax3.set_ylabel('Amplitude (mV)', fontsize=nf)
    ax3.set_title('QRS Complex', fontsize=nf)
    ax3.grid(True)
    
    # Plot da onda T
    ax4 = plt.subplot(2,3,6)
    for i in range(T.shape[1]):
        ax4.plot(np.arange(len(T[:,i]))/fs, T[:,i], 'b-', alpha=0.1)
    ax4.plot(np.arange(len(T[:,0]))/fs, np.mean(T, axis=1), 'r-', linewidth=2)
    ax4.set_xlabel('Time (s)', fontsize=nf)
    ax4.set_ylabel('Amplitude (mV)', fontsize=nf)
    ax4.set_title('T Wave', fontsize=nf)
    ax4.grid(True)
    
    plt.tight_layout()
    return fig

# Uso:
fig = plot_ecg_segments(B, P, QRS, T, fs)
plt.savefig('./graficos/ecg_segments.png', dpi=300, bbox_inches='tight')
plt.close()


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
from ecgdetectors import Detectors
import detectionPeaks as dt
#import wfdb as wf
#import ecg_plot as ecg
import matplotlib.pyplot as plt
import matplotlib as mt

sampleMap = {
    "FA": 256,
    "SINUS": 256
}

WINDOW_SIZE = 12000
CURRENT_SAMPLE = "FA"

#plt.rcParams['pgf.texsystem'] = 'xelatex'

#os.environ["PATH"] += os.pathsep + '/usr/local/texlive/2018/bin/x86_64-darwin'
print("recursion limite", sys.getrecursionlimit())

#record = wf.rdsamp('16265')
# record = wf.rdsamp("04015")
# df =pd.DataFrame(record[0], columns=record[1]['sig_name'])
# #df.to_csv('16265.csv')
# df.to_csv("04015.csv")
df_csv = pd.read_csv('./data/04015.csv')

signal = df_csv._get_label_or_level_values('ECG1').tolist()
signalY = df_csv._get_label_or_level_values('ECG2').tolist()

signalXY = [signal, signalY]
#signal = signal[5000:100000]
#fs = round(81*0.7,2)
fs = sampleMap.get(CURRENT_SAMPLE)
detectors = Detectors(fs)
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
    """Normaliza o sinal dividindo pelo desvio padrão."""
    return np.array(signal) / np.std(signal)



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

def getSnr (noise, signal, title):
    originalNoise = min_max_scaling(np.array(noise[0:len(signal)]))
    filtredSignal = min_max_scaling(np.array(signal))
    noise = originalNoise - filtredSignal
    snr = 10 * np.log10(np.sum(filtredSignal**2) / np.sum(noise**2))
    signal_power = np.sum(filtredSignal**2)
    noise_power = np.sum(noise**2)
    print(f"Signal power: {signal_power}")
    print(f"Noise power: {noise_power}")
    print(f"Power ratio: {signal_power/noise_power}")
    print(f"SNR (dB) {title}: {snr}")


start = int(round(len(signal) * 0.01))
end = len(signal) - int(round(len(signal) * 0.01))
windowedSignal = signal[start:end] #tira 1% do começo e do fim
#filtered_signal = bandpass_filter(windowedSignal, fs, lowcut=0.5, highcut=40.0, order=4)
signal_mean = signal_mean(windowedSignal)
#signal_std = signal_std(signal_mean)
final_signal =  ((windowedSignal - np.mean(windowedSignal)) / np.std(windowedSignal))
#final_signal = signal_std
#print(final_signal)

# Normalize both signals before SNR calculation
getSnr(signal, final_signal, "pré processamento")

print("processando...")
start_time = time.time()  # Marca o início do tempo
peaks = dt.dtPeaks(final_signal[0:WINDOW_SIZE], [0,60], fs, 0)
end_time = time.time()  # Marca o fim do tempo
print("processado!")

# Calcula o tempo de execução
execution_time = end_time - start_time
print(f"Tempo de execução da função dtPeaks: {execution_time:.4f} segundos")
#rr_peaks = detectors.pan_tompkins_detector(final_signal)
ecgH = peaks[0]
qrs_amplitude = peaks[1]
qrs_index = peaks[2]
delay = peaks[3]
peaks_array = pd.array(peaks[0], int)

getSnr(final_signal, ecgH, "processado")

#print(f"qrsIndex: {qrs_index}")
#print(f"qrsamplitude: {qrs_amplitude}")
#print(qrs_amplitude)

#-------------------------CICLO-----------------------------------#
if qrs_index[0] < round(theta * fs):
    qrs_index = np.delete(qrs_index, 0)
    qrs_amplitude = np.delete(qrs_amplitude, 0)

if qrs_index[-1] > (len(signal) - round(lbda * fs)):
    qrs_index = np.delete(qrs_index, -1)
    qrs_amplitude = np.delete(qrs_amplitude, -1)


# Batimentos
B = np.zeros((len(qrs_amplitude), round(theta * fs) + round(lbda * fs) + 1))
P = np.zeros((len(qrs_amplitude), round(lbdaP * fs) + 1))
QRS = np.zeros((len(qrs_amplitude), round(thetaQRS * fs) + round(lbdaQRS * fs) + 1))
T = np.zeros((len(qrs_amplitude), round(0.3 * fs) + 1))

for i in range(len(qrs_amplitude)):
    # Batimentos
    start_B = qrs_index[i] - round(theta * fs)
    end_B = qrs_index[i] + round(lbda * fs)
    B[i] = signal[int(start_B) : int(end_B) + 1]  # +1 para incluir o end_B

    # Onda P
    start_P = qrs_index[i] - round(theta * fs)
    end_P = start_P + round(lbdaP * fs)
    P[i] = signal[int(start_P) : int(end_P) + 1]

    # Complexo QRS
    start_QRS = qrs_index[i] - round(thetaQRS * fs)
    end_QRS = qrs_index[i] + round(lbdaQRS * fs)
    QRS[i] = signal[int(start_QRS) : int(end_QRS) + 1]

    # Onda T
    start_T = (qrs_index[i] + round(lbda * fs)) - round(0.3 * fs)
    end_T = qrs_index[i] + round(lbda * fs)
    T[i] = signal[int(start_T) : int(end_T) + 1]

B = B.T
P = P.T
QRS = QRS.T
T = T.T

# signal = signal[0:175]
# peaks[0] = peaks[0][0:175]

timeAxis = np.arange(len(peaks[0])) / fs
timeAxisNormalSignal = np.arange(len(signal))/fs

def saveIndivualsToPng(figs):
    date = datetime.today().isoformat(timespec='seconds')
    for i, fig_item in enumerate(figs):
        if isinstance(fig_item, plt.Figure):
            data = fig_item.axes[0].lines[0].get_data()
            fig_item.savefig(f'./graficos/{CURRENT_SAMPLE}/individual_fig_{i}.png', dpi = 300, bbox_inches = 'tight')
        plt.tight_layout()
        plt.close()

def saveToPng(figs):
    date = datetime.today().isoformat(timespec='seconds')

    fig, axs = plt.subplots(2, 2, figsize=(15, 5))
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
            #fig_item.savefig(f'./graficos/{CURRENT_SAMPLE}/individual_fig_{i}.png', dpi = 300, bbox_inches='tight')
            plt.close(fig_item)  
    
    plt.tight_layout()
    
    plt.savefig(f'./graficos/{CURRENT_SAMPLE}/ecg_analysis_{date}.png', dpi=300, bbox_inches='tight')
    plt.close()


def plotSignal(signal, title, reduce = 0, invert = False):


    timeAxis = []
    if(reduce != 0):
        signal = signal[0:reduce]
        timeAxis = np.arange(len(signal)) / fs
    else:
        timeAxis = np.arange(len(signal)) / fs
    plot = plt.figure(figsize=(10, 6))
    
    # plt.plot(timeAxis, signal)
    if(invert):
        plt.plot(signal, timeAxis)
    else:
        plt.plot(timeAxis,signal)
    # plt.plot(timeAxisNormalSignal, signal)
    plt.xlabel('Amostra')
    plt.ylabel('Sinal (mv)')
    plt.title(title)
    plt.grid(True)  # Add a grid for better readability
    #plt.show()
    return plot

def plot_ecg_segments(B, P, QRS, T, fs):
    nf = 16
    fig = plt.figure(figsize=(15, 10))
    
    # Configurações gerais
    plt.rcParams.update({'font.size': nf})
    
    # Full beats plot (top)
    ax1 = plt.subplot(2,1,1)
    for beat in B.T:  # Plotar cada batimento separadamente
        ax1.plot(beat, linewidth=1)
    ax1.set_xlabel('Samples', fontsize=nf)
    ax1.set_ylabel('mV', fontsize=nf)
    ax1.set_title('Healthy beats (50000)', fontsize=nf)
    ax1.tick_params(labelsize=nf)
    
    # P wave plot
    ax2 = plt.subplot(2,3,4)
    for p_wave in P.T:
        ax2.plot(p_wave, linewidth=1)
    ax2.set_xlabel('Samples', fontsize=nf)
    ax2.set_ylabel('mV', fontsize=nf)
    ax2.set_title('P waves (50000)', fontsize=nf)
    ax2.tick_params(labelsize=nf)
    
    # QRS complex plot
    ax3 = plt.subplot(2,3,5)
    for qrs in QRS.T:
        ax3.plot(qrs, linewidth=1)
    ax3.set_xlabel('Samples', fontsize=nf)
    ax3.set_ylabel('mV', fontsize=nf)
    ax3.set_title('QRS complexes (50000)', fontsize=nf)
    ax3.tick_params(labelsize=nf)
    
    # T wave plot
    ax4 = plt.subplot(2,3,6)
    for t_wave in T.T:
        ax4.plot(t_wave, linewidth=1)
    ax4.set_xlabel('Samples', fontsize=nf)
    ax4.set_ylabel('mV', fontsize=nf)
    ax4.set_title('T waves (50000)', fontsize=nf)
    ax4.tick_params(labelsize=nf)
    
    plt.tight_layout()
    return fig

#
#
#rr_peaks_figure = plotSignal(signal[0:3000], "PICOS RR")
signal_figure = plotSignal(signal, "SINAL (ORIGINAL)", 500)
final_signal_figure = plotSignal(final_signal, "SINAL (PRE PROCESSADO)", 500)
ecgProcessado = plotSignal(ecgH, "SINAL (FILTRADO E PROCESSADO)", 500)
#batimentos = plotSignal(B[0], "BATIMENTOS", 300)
complexoQrs = plotSignal(QRS, "COMPLEXO QRS", WINDOW_SIZE)
ondaT = plotSignal(T, "ONDA T", WINDOW_SIZE)
ondaP = plotSignal(P, "ONDA P", WINDOW_SIZE)


saveIndivualsToPng([signal_figure, final_signal_figure, ecgProcessado, complexoQrs, ondaT, ondaP])
saveToPng([signal_figure, final_signal_figure, ecgProcessado])
#plot_segments(B, P, QRS, T, fs)


# Uso:
fig = plot_ecg_segments(B[0:WINDOW_SIZE], P[0:WINDOW_SIZE], QRS[0:WINDOW_SIZE], T[0:WINDOW_SIZE], fs)
plt.savefig(f'./graficos/{CURRENT_SAMPLE}/ecg_segments{datetime.now()}.png', dpi=300, bbox_inches='tight')
plt.close()


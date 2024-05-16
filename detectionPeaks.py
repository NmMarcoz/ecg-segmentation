#initialize
import pandas as pd
import numpy as np
import scipy as sc




def dtPeaks(ecg,fs,gr,qrsAmpRaw,qrsIRaw,delay,method):
    qrsC = []
    qrsI = []
    SIG_LEV = 0
    noisC = []
    noisI = []
    delay = 0
    skip = 0
    notNois = 0
    selectedRR = []
    mSelectedRR = []
    meanRR = 0
    qrsIRaw = []
    qrsAmpRaw = []
    serBack = 0
    testM = 0
    siglBuf = []
    noislBuf = []
    thrsBuf = []
    siglBuf1 = []
    noislBuf1 = []
    # t = <<-- aqui tenho que botar um intervalo de 5 minutos. Ainda não tenho muita ideia de como funciona.

    if (fs == 200):
        #PRIMEIRO PASSO
        # Primeiro filtro: Passa-baixa
        b = [1, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 1]
        a = [1,-2,1]
        hL = sc.lfiter(b,a,[1, np.zeros((1,12))])
        ecgL = np.convolve(ecg, hL)
        ecgL = ecgL/np.max(ecgL)

        #SEGUNDO PASSO
        # Segundo filtro: Passa-alta

        b = [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32,-32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        a = [1,-1]

        hH = sc.lfilter(b,a, [1, np.zeros((1,32))])
        ecgH = np.convolve(ecgL, hH)
        ecgH = ecgH/ np.max((( np.abs(ecgH))))

    else:
        # Filtro Passa-banda
        # Filtro passa-banda para o cancelamento de artefatos de baixa e alta frequência

        f1 = 5
        f2 = 15
        n = 3
        rp = 1
        rs = 20
        ## Filtro elípitico de 2nd ordem de tempo contínuo

        [b,a] = sc.signal.butter(2,[f1,f2]/(fs/2))
        [b,a] = sc.signal.ellip(3,1,20,[f1,f2]/(fs/2))

        ecgH = sc.lfilter(b,a, ecg)
        ecgH = ecgH/ np.max(np.abs(ecgH))




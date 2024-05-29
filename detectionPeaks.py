#initialize
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
        f1Normal = f1 / (fs/2)
        f2Normal = f2 / (fs/2)
        b,a = sc.signal.butter(2,[f1Normal,f2Normal], "bandpass")
       # b,a = sc.signal.ellip(3,1,20,[f1,f2]/(fs/2))

        ecgH = sc.lfilter(b,a, ecg)
        ecgH = ecgH/ np.max(np.abs(ecgH))

    hD = [-1,-2,0,2,1]*(1/8)
    ecgD = np.convolve(ecgH, hD)
    ecgD = ecgD/np.max(ecgD)

    ecgS = ecgD**2

    ecgM = np.convolve(ecgS, np.ones(1, round(0.150*fs))/round(0.150*fs))

    [pks, locs] = sc.signal.find_peaks(ecgM, round(0.2*fs))
    #parei aqui

    thrSig = (np.max(ecgM[:2*fs]))/3 #25% do sinal
    thrNoise = (np.mean(ecgM[:2*fs]))/2
    sigLevel = thrSig
    noiseLevel = thrNoise

    #--------------------------------------------

    thrSig1 = (np.max(ecgH[:2*fs]))/3
    thrNoise1 = (np.max(ecgH[:2*fs]))/2
    sigLevel1 = thrSig1
    noiseLevel1 = thrNoise1

    #--------------------------------------------
    y_i = np.zeros_like(locs, dtype=float)
    x_i = np.zeros_like(locs, dtype=int)
    for i in pks:
        if(locs(i)-round(0.150*fs)>=1) and locs(i)<=len(ecgH):
            y_i[i], x_i[i] = np.max(ecgH[locs(i)-round(0.150*fs):locs(i)])
        else:
            if i ==1:
                y_i[i], x_i[i] = np.max(ecgH[1:locs(i)])
                serBack = 1
            elif locs(i)>= len(ecgH):
                y_i[i], x_i[i] = np.max(ecgH[locs(i)-round(0.150*fs):locs(i)])

    #### update the heart_rate (two heart rate mans one of the most recent and the other selected
        if(len(qrsC))>=9:
            diffRR = np.diff(qrsI[-8:-1]) #calculando o intervalo RR
            meanRR = np.mean(diffRR) ## média dos 8 números anteriores
            comp = qrsI[-1]-qrsI(-2) # Últimos RR
            if(comp <= 0.92*meanRR or comp >= 1.16*meanRR):
                thrSig = 0.5*(thrSig)
                thrSig1 = 0.5*(thrSig1)
            else:
                mSelectedRR = meanRR ## média das últimas batidas regulares


        #calculo da média das últimas 8 ondas R para ter certeza que o QRS
        # não está perdido

        if (testM):
            if(locs[i] - qrsI[-1] >= round(1.66*testM)):
                pksTemp, locsTemp = np.max(ecgM(qrsI[-1+round*0.200*fs:locs[i]-round(0.200*fs)]))
                locsTemp = qrsI[-1] + round(0.200*fs) + locsTemp -1

                if pksTemp > thrNoise:
                    qrsC = [qrsC, pksTemp]
                    qrsI = [qrsI, locsTemp]

                    if locsTemp <= len(ecgH):
                        y_i_t, x_i_t = np.max(ecgH[locsTemp - round(0.150*fs):locsTemp])
                    else:
                        y_i_t, x_i_t = np.max(ecgH[locsTemp - round(0.150*fs): -1])
                    if y_i_t > thrNoise:
                        qrsIRaw = [qrsIRaw, locsTemp-round(0.150*fs) + (x_i_t -1)]
                        qrsAmpRaw = [qrsAmpRaw, y_i_t]
                        sigLevel1 = 0.25*y_i_t + 0.75*sigLevel1
                    notNois = 1
                    sigLevel = 0.25*pksTemp + 0.75*sigLevel
                else:
                    notNois = 0


        ## Detectando ruído nos picos QRS
        if(pks[i] >= thrSig):
            if(len(qrsC)>=3):
                if(locs[i] - qrsI[-1]<= round(0.3600*fs)):
                    slope1 = np.mean(np.diff(ecgM[locs[i]-round(0.075*fs):locs[i]]))
                    slope2 = np.mean(np.diff(ecgM[qrsI[-1]-round(0.075*fs):qrsI[-1]]))

                    if np.abs(slope1)<= np.abs(0.5*slope2):
                        noisC= [noisC, pks[i]]
                        noisI = [noisI, locs[i]]
                        skip =1 # indentificação da onda T

                        noiseLevel1 = 0.125*y_i + 0.875*noiseLevel1
                        noiseLevel = 0.125*pks[i] + 0.875*noiseLevel
                    else:
                        skip = 0
            if skip == 0:
                qrsC = [qrsC, pks[i]]
                qrsI = [qrsI, locs[i]]

            if y_i >= thrSig:
                if serBack:
                    qrsIRaw = [qrsIRaw, x_i] #salvando o index do bandpass
                else:
                    qrsIRaw = [qrsIRaw, locs[i]-round(0.150*fs)+(x_i-1)]

                qrsAmpRaw = [qrsAmpRaw, y_i]
                sigLevel1 = 0.125*y_i + 0.875*sigLevel1

            sigLevel = 0.125*pks[i] + 0.875*sigLevel

        elif (thrNoise <= pks[i] and pks[i] < thrSig):
            noiseLevel1 = 0.125*y_i + 0.875*noiseLevel1
            noiseLevel = 0.125*pks[i] + 0.875*noiseLevel
        elif (pks[i]< thrNoise):
            noisC = [noisC, pks[i]]
            noisI = [noisI, locs[i]]

            noiseLevel1 = 0.125*y_i + 0.875*noiseLevel1

            noiseLevel = 0.125*pks[i] + 0.875*noiseLevel

        if (not(noiseLevel == 0) or not(sigLevel == 0)):
            thrSig = noiseLevel + 0.25*(np.abs(sigLevel - noiseLevel))
            thrNoise = 0.5*(thrSig)

        if(not(noiseLevel1 == 0) or not(sigLevel1 == 0)):
            thrSig1 = noiseLevel1 + 0.25*(np.abs(sigLevel - noiseLevel1))
            thrNoise1 = 0.5*(thrSig1)


        #-----------------------------
        # take a track of thresholds of smoothed signal
        siglBuf = [siglBuf, sigLevel]
        noislBuf = [noislBuf, noiseLevel]
        thrsBuf = [thrsBuf, thrSig]

        #take a track of thresholds of filtered signal
        siglBuf1 = [siglBuf1, sigLevel1]
        noislBuf1 = [noislBuf1, noiseLevel1]
        thrsBuf1 = [thrsBuf1, thrSig1]

        skip = 0
        notNois = 0
        serBack = 0
























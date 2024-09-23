#initialize
import numpy as np
import scipy as sc




#def dtPeaks(ecg,fs,gr,qrsAmpRaw,qrsIRaw,delay,method):
def dtPeaks(ecg, min, fs, flagFigura):
    qrsC = []
    qrsI = []
    SIG_LEV = 0
    noisC = []
    noisI = []
    delay = 0
    skip = 0
    notNois = 0
    selectedRR = []
    mSelectedRR = 0
    meanRR = 0
    qrsIRaw = []
    qrsAmpRaw = []
    serBack = 0
    testM = 0
    siglBuf = []
    noislBuf = []
    thrsBuf = []
    thrsBuf1 = []
    siglBuf1 = []
    noislBuf1 = []
    # t = <<-- aqui tenho que botar um intervalo de 5 minutos. Ainda não tenho muita ideia de como funciona.
    ts = 1/fs
    fc = 15
    fc1 = 5
    wn = fc/(fs/2)
    wn1 = fc1/(fs/2)
    ordF = 60

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

        hH = sc.signal.lfilter(b,a, [1, np.zeros((1,32))])
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
        #b = sc.signal.butter(2,[f1Normal,f2Normal], "bandpass")
       # b,a = sc.signal.ellip(3,1,20,[f1,f2]/(fs/2))
        ecgH = sc.signal.lfilter(b,a, ecg)
        ecgH = ecgH/ np.max(np.abs(ecgH))

    hD = np.array([-1,-2,0,2,1])
    hD = hD/8
    ecgD = np.convolve(ecgH, hD)
    ecgD = ecgD/np.max(ecgD)

    ecgS = ecgD**2
    window_size = round(0.150*fs)
    ecgM = np.convolve(ecgS, np.ones(window_size)/round(0.150*fs))
    #print(ecgM)
    print(len(ecgM))
    print(len(ecgH))

    data = sc.signal.find_peaks(ecgM, distance = round(0.2*fs))
    #data = sc.signal.find_peaks(ecgM, height=2)
    #pks = PICOS
    #pks = sc.signal.find_peaks(ecgM, distance =round(0.2*fs))[0]
    pks = sc.signal.find_peaks(ecgM, distance = round(0.2*fs))[0]
    #locs = MOMENTOS DE PICO

    locs = data[0]

    #pks = sc.signal.find_peaks(ecgM, round(0.2*fs))
    print(f"pks: {pks}")
    print(f"locs: {locs}" )
    print(f"data: {data[0]}")
    print(f"data: {data[1]}")
    sliceCut = int(2*fs)
    #thrSig = (np.max(ecgM[:sliceCut]))/3 #25% do sinal
    thrSig = np.max(ecgM[1:int(2*fs)]*1/3)
    thrNoise = np.mean(ecgM[:int(2*fs)]*0.5)
    #thrNoise = (np.mean(ecgM[:sliceCut]))/2
    sigLevel = thrSig
    noiseLevel = thrNoise

    #--------------------------------------------

    thrSig1 = (np.max(ecgH[:sliceCut]))/3
    thrNoise1 = (np.mean(ecgH[:sliceCut]))/2
    sigLevel1 = thrSig1
    noiseLevel1 = thrNoise1


    #--------------------------------------------
    #y_i = np.zeros_like(locs, dtype=float)
    #x_i = np.zeros_like(locs, dtype=int)
    y_i = 0
    x_i = 0
    y_i_t = 0
    x_i_t = 0
    pksTemp = 0
    print(f"noiseLevel1: {noiseLevel1}")
    print(f"thsig1: {thrSig1}")

    for i in range(len(locs)):
        #locate the corresponding peak in the filtered signal
        print(f"thrnoise = {thrNoise}")
        print(f"thsig1: {thrSig1}")
        print(f"y_i: {y_i}")
        print(f"locs[i]: {locs[i]}")
        print(f"tamanho do ecgH: {len(ecgH)}")
        if(locs[i]-round(0.150*fs) >=1) and locs[i]<=len(ecgH):
            print("entrou na condição do locs i< ecgH")
            y_i = np.max(ecgH[locs[i]-round(0.150*fs):locs[i]])
            x_i = np.argmax(ecgH[locs[i] - round(0.150*fs):locs[i]])
        else:
            if i ==1:
                y_i = np.max(ecgH[1:locs[i]])
                x_i = np.argmax(ecgH[1:locs[i]])
                serBack = 1
            elif locs[i]>= len(ecgH):
                y_i  = np.max(ecgH[locs[i]-round(0.150*fs):])
                x_i=  np.argmax(ecgH[locs[i]-round(0.150*fs):])

    #update the heart_rate (two heart rate mans one of the most recent and the other selected
        if(len(qrsC)>=9):
            print("entrou no tamanho do qrsC >= 9")
            diffRR = np.diff(qrsI[-8:]) #calculando o intervalo RR
            meanRR = np.mean(diffRR) ## média dos 8 números anteriores
            comp = qrsI[-1]-qrsI[-2] # Últimos RR
            # print(f"comp: {comp}")
            # print(f"meanRR: {meanRR}")
            # print(f"0.92 * meanRR: {meanRR*0.92}")
            # print(f"1.16 * meanRR: {meanRR*1.16}")
            if(comp <= 0.92*meanRR or comp >= 1.16*meanRR):
                print("e deu update no thrSig * 0.5")
                thrSig *= 0.5
                #thrNoise = 0.5*(thrSig) # nao sei se isso aqui <- é pra ficar
                thrSig1 *= 0.5
            else:
                mSelectedRR = meanRR ## média das últimas batidas regulares


        #calculo da média das últimas 8 ondas R para ter certeza que o QRS
        # não está perdido
        if (mSelectedRR):
            print("mselectedRR")
            testM = mSelectedRR
        elif (meanRR and mSelectedRR == 0):
            testM = meanRR
        else:
            testM = 0

        if (testM):
            # print(f"teste M: {testM}")
            # print(f"qrsI : {qrsI}")
            # print(f"pksTemp : {pksTemp}")
            # print(f"valor de thrNoise : {thrNoise}")
            # print(f"valor de thrSig: {thrSig}")
            # print(f"valor do y_i_t: {y_i_t}")
            # print(f"valor de pks[i]: {pks[i]}")
            if(locs[i] - qrsI[-1] >= round(1.66*testM)):
                print("TO MALUCO")
                pksTemp= np.max(ecgM[qrsI[-1] + round(0.200*fs): locs[i] - round(0.200*fs)])
                locsTemp = np.argmax(ecgM[qrsI[-1] + round(0.200*fs): locs[i] - round(0.200*fs)])
                locsTemp = qrsI[-1] + round(0.200*fs) + locsTemp -1
                print(f"y_i_t = {y_i_t}")
                print(f"pksTemp: {pksTemp}")
                print(f"thrNoise: {thrNoise}")
                if pksTemp > thrNoise:
                    print("pksTemp > thrNoise")
                    #qrsC = [qrsC, pksTemp]
                    #qrsI = [qrsI, locsTemp]
                    qrsC = np.append(qrsC, pksTemp)
                    qrsI = np.append(qrsI, locsTemp)

                    if locsTemp <= len(ecgH):
                        print("locs temp menor que ech")
                        y_i_t= np.max(ecgH[int(locsTemp - round(0.150*fs)):int(locsTemp)])
                        x_i_t = np.argmax(ecgH[int(locsTemp - round(0.150*fs)):int(locsTemp)])
                    else:
                        print("locs maior que ech")
                        y_i_t = np.max(ecgH[locsTemp - round(0.150*fs): -1])
                        x_i_t = np.argmax(ecgH[locsTemp - round(0.150*fs): -1])
                if y_i_t > thrNoise1:
                    print("y_i_t > thrNoise")
                    qrsIRaw = [qrsIRaw, locsTemp-round(0.150*fs) + (x_i_t -1)]
                    qrsAmpRaw = [qrsAmpRaw, y_i_t]
                    sigLevel1 = 0.25*y_i_t + 0.75*sigLevel1
                notNois = 1
                sigLevel = 0.25*pksTemp + 0.75*sigLevel
            else:
                notNois = 0


        ## Detectando ruído nos picos QRS
        print(f"PKS I: {pks[i]}")
        if(pks[i] >= thrSig):
            print("pks i maior que thrSig")
            #print(f"valor de y_i: {y_i}")
            if(len(qrsC)>=3):
                print("SESSAO DO QRS C > 3")
                print("len de qrsC maior que 3")
                if(locs[i] - qrsI[-1]<= round(0.3600*fs)):
                    print("locs i < qrsI CHECK")
                    slope1 = np.mean(np.diff(ecgM[locs[i]-round(0.075*fs):locs[i]]))
                    slope2 = np.mean(np.diff(ecgM[qrsI[-1]-round(0.075*fs):qrsI[-1]]))

                    if np.abs(slope1)<= np.abs(0.5*slope2):
                        print("ABS <= SLOPE2 CHECK")
                        noisC= [noisC, pks[i]]
                        noisI = [noisI, locs[i]]
                        skip =1 # indentificação da onda T

                        noiseLevel1 *= 0.125*y_i + 0.875
                        print(f"noiseLevel1 2nd vez: {noiseLevel1}")
                        noiseLevel = 0.125*pks[i] + 0.875*noiseLevel
                        print(f"pks I: {pks[i]}")
                    else:
                        skip = 0
            print("SESSAO TERMINADA")
            print("----------------")
            if skip == 0:
                #qrsC = np.concatenate((qrsC, pks))
                print("skipou")
                qrsC.append(pks[i])
                qrsI.append(locs[i])
                #qrsI = np.concatenate((qrsI, locs))

            #if np.all(y_i > thrSig):

            if y_i >= thrSig1:
                print("Y_I MAIOR QUE TRHSIG")
                if serBack:
                    qrsIRaw = [qrsIRaw, x_i] #salvando o index do bandpass
                else:
                    qrsIRaw = [qrsIRaw, locs[i]-round(0.150*fs)+(x_i-1)]

                #qrsAmpRaw = [qrsAmpRaw, y_i]
                qrsAmpRaw.append(y_i)
                print(f"valor do qrsAmpRaw {qrsAmpRaw}")
                sigLevel1 = 0.125*y_i + 0.875*sigLevel1

            sigLevel = 0.125*pks[i] + 0.875*sigLevel

        elif (thrNoise <= pks[i] and pks[i] < thrSig):
            print("ZUBUMAFU")
            noiseLevel1 = 0.125*y_i + 0.875*noiseLevel1
            #print(f"noiseLevel1 terceira vez: {noiseLevel1}")
            noiseLevel = 0.125*pks[i] + 0.875*noiseLevel
        elif (pks[i]< thrNoise):
            noisC = [noisC, pks[i]]
            noisI = [noisI, locs[i]]

            noiseLevel1 = 0.125*y_i + 0.875*noiseLevel1
            #print(f"noiseLevel1 quarta vez: {noiseLevel1}")

            noiseLevel = 0.125*pks[i] + 0.875*noiseLevel

        if (noiseLevel != 0) or (sigLevel != 0):
            print("ajustou thrnoise aqui")
            print(f"thrnoise antes do ajuste {thrNoise}")

            thrSig = noiseLevel + 0.25*(np.abs(sigLevel - noiseLevel))
            thrNoise = 0.5*(thrSig)
            print(f"thrnoise depois do ajuste {thrNoise}")
            #print(f"tipo do noiselevel 1: {noiseLevel1}")
            #print(f"tipo do noiselevel : {noiseLevel}")

        if(noiseLevel1 != 0) or (sigLevel1 != 0):
            print(".............")
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

    return [ecgH,  qrsAmpRaw, qrsIRaw, delay]






















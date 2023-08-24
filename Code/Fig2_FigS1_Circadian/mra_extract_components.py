from pylab import*
import pandas as pd
import pyboat
import os
import pyboat.plotting as pl
import pywt
import numpy as np
from modwt import*
import xlsxwriter

# CHANGE LINES 130-132 AND 151 !

def plotDWT(t, signal, NumLevels, c_wtmra, dt):
    f, axarr = plt.subplots(NumLevels+3, sharex=True, figsize=(6/1.5, 12/1.5))
    #print( signal )
    c_max = ceil(max(signal.flatten()))
    for i, j in enumerate(c_wtmra):     # 7 elements 6
        axarr[i].plot(t, j)
        #axarr[i].set_ylim(-c_max, c_max)
        axarr[i].set_xticks(arange(0, t[-1]+24., 24.))
        #axarr[i].set_yticks([-c_max, 0, c_max])
        axarr[i].set_xlim(0, t[-1] - t[-1] % 24 + 24)
        if i == NumLevels:
            axarr[i].set_title("Wavelet Smooth", fontsize=10)     # Last One is the smootheness coefficient of the last detail which corresponds to the trend of the signal
        else:
            axarr[i].set_title("D"+str(i+1) + r"$\rightarrow$ ["+str(round((2.**(i+1)*dt), 2))+"h, "+ str(round((2.**(i+2)*dt), 2)) +"h[",  fontsize=10)
    axarr[i+1].plot(t, signal, "k-", linewidth=2, label="Signal")
    axarr[i+1].plot(t, sum(c_wtmra, axis=0), linestyle="--", color="red", label="$ D_{\sum} + S$")
    axarr[i+1].set_ylim(-c_max, c_max)
    axarr[i+1].set_yticks([-c_max, 0, c_max])
    axarr[i+1].plot(t, j, color="blue", label=("D"+str(i+1)))
    leg1 = axarr[i+1].legend(loc="upper right", prop={'size':10}, fancybox=True)
    leg1.get_frame().set_alpha(0.5)


    axarr[i+2].plot(t, signal, color="k", linestyle="-", label="Signal")
    axarr[i+2].plot(t, signal - j, color="gray", linestyle="-", label="Signal - S")
    axarr[i+2].plot(t, c_wtmra[i-3], color="blue", label = "D"+str(NumLevels-2), linestyle="--")
    leg2 = axarr[i+2].legend(loc="upper right", prop={'size':10}, fancybox=True)
    leg2.get_frame().set_alpha(0.5)
    axarr[i+2].set_ylim(-c_max, c_max)
    axarr[i+2].set_yticks([-c_max, 0, c_max])
    axarr[i+2].set_xlabel("Time in h")
    f.tight_layout()

def PlotMRA(t, signal, MRA, dt=30./60, NumLevels=7, t_edge=1):
    f, axarr = plt.subplots(NumLevels+2, sharex=True, figsize=(6/1.7, 12/1.7))
    c_max = ceil(max(signal.flatten()))
    
    i_edge = round(t_edge/dt)
    ii_edge = round(t_edge/dt) #defined by CE, original version ii_edge = i_edge
    
    OverallEnergy = sum([sum(n[i_edge:-ii_edge]**2) for n in MRA[:-1]])/len(t[i_edge:-ii_edge])
    
    for i, j in enumerate(MRA):     # 7 elements 6
        axarr[i].plot(t, j, color="gray", linestyle="--")
        
        # calculate energy in detail j
        DetailEnergy = sum(j[i_edge:-ii_edge]**2) / len(t[i_edge:-ii_edge]) / OverallEnergy
        
        axarr[i].plot(t[i_edge:-ii_edge], j[i_edge:-ii_edge], label = str(round(DetailEnergy*100,2)) + r"% in " + r"["+str(round((2.**(i+1)*dt)))+"h, "+ str(round((2.**(i+2)*dt))) +"h[")
        axarr[i].set_xticks(arange(0, t[-1]+24., 24.))
        axarr[i].set_xlim(0, t[-1])
        axarr[i].set_ylabel("D$_{"+str(i+1) +"}$")
        leg2 = axarr[i].legend(loc="upper right", prop={'size':10}, fancybox=True, handlelength=0, handletextpad=0)
        leg2.get_frame().set_alpha(0.5)
    axarr[i+1].plot(t, signal, "k-", linewidth=2, label="Signal")
    axarr[i+1].plot(t, sum(MRA, axis=0), linestyle="--", label="$ D_{\sum} + S$")
    axarr[i+1].set_ylim(-c_max, c_max)
    axarr[i+1].set_yticks([-c_max, 0, c_max])
    leg1 = axarr[i+1].legend(loc="upper right", prop={'size':8}, fancybox=True, ncol=2)
    leg1.get_frame().set_alpha(0.5)

    axarr[i+1].set_xlabel("Time in h")
    f.tight_layout()

def PlotMRA_CoarseGrained(t, signal, MRA, dt=30./60, NumLevels=7, t_edge=1): #before t.edge = 12.
    f, axarr = plt.subplots(5, sharex=True, figsize=(6/1.7, 12/2.5))
    c_max = ceil(max(signal.flatten()))
    
    i_edge = round(t_edge/dt) #original version has one t_edge/i_edge and applies it to both sides
    ii_edge = round(t_edge/dt) #defined by CE, original version ii_edge = i_edge
    
    
    OverallEnergy = sum([sum(n[i_edge:-ii_edge]**2) for n in MRA[:-1]])/len(t[i_edge:-ii_edge])
    D1Energy = sum(MRA[0][i_edge:-ii_edge]**2) / len(t[i_edge:-ii_edge]) / OverallEnergy
    D2Energy = sum(MRA[1][i_edge:-ii_edge]**2) / len(t[i_edge:-ii_edge]) / OverallEnergy
    D3Energy = sum(MRA[2][i_edge:-ii_edge]**2) / len(t[i_edge:-ii_edge]) / OverallEnergy
    D4Energy = sum(MRA[3][i_edge:-ii_edge]**2) / len(t[i_edge:-ii_edge]) / OverallEnergy
    D5Energy = sum(MRA[4][i_edge:-ii_edge]**2) / len(t[i_edge:-ii_edge]) / OverallEnergy
    D6Energy = sum(MRA[5][i_edge:-ii_edge]**2) / len(t[i_edge:-ii_edge]) / OverallEnergy
    D7Energy = sum(MRA[6][i_edge:-ii_edge]**2) / len(t[i_edge:-ii_edge]) / OverallEnergy
    D8Energy = sum(MRA[7][i_edge:-ii_edge]**2) / len(t[i_edge:-ii_edge]) / OverallEnergy

     
    D12  = MRA[0] + MRA[1] #Noise
    D34  = MRA[2] + MRA[3] #Ultradian
    D5   = MRA[4] #Circadian
    D678 = MRA[5] + MRA[6] + MRA[7] #Infradian
    
    c_label  = ["Noise", "Ultradian", "Circadian", "Infradian"]
    test  = ["Noise", "Ultradian", "Circadian", "Infradian"]
    c_energy = [D1Energy+D2Energy, D3Energy+D4Energy, D5Energy, D6Energy+D7Energy+D8Energy]
    
    for i, j in enumerate([D12, D34, D5, D678]):
        axarr[i].plot(t, j, color="gray", linestyle="--")
        
        c_maxtmp = max(abs(j[i_edge:-ii_edge])) 
        axarr[i].plot(t[i_edge:-ii_edge], j[i_edge:-ii_edge], label = str(round(c_energy[i]*100,2)))
        axarr[i].set_xticks(arange(0, t[-1]+24., 24.))
        axarr[i].set_xlim(0, t[-1])
        axarr[i].set_ylim(-c_maxtmp, c_maxtmp) 
        axarr[i].set_ylabel(c_label[i])
        leg2 = axarr[i].legend(loc="upper right", prop={'size':10}, fancybox=True, handlelength=0, handletextpad=0)
        leg2.get_frame().set_alpha(0.5)
        
        finalvalues = str(round(c_energy[i],4))
        print(col,finalvalues)
        worksheet.write(col+1,i+1,finalvalues) 

        
    axarr[4].plot(t, signal, "k-", linewidth=2, label="Signal")
    axarr[4].plot(t, sum(MRA, axis=0), linestyle="--", label="$ D_{\sum} + S$")
    axarr[4].set_ylim(-c_max, c_max)
    axarr[4].set_yticks([-c_max, 0, c_max])
    leg1 = axarr[4].legend(loc="upper right", prop={'size':8}, fancybox=True, ncol=2)
    leg1.get_frame().set_alpha(0.5)

    axarr[4].set_xlabel("Time in h")
    f.tight_layout()
                      
# Load data
path='/Users/carolinector/Nextcloud/Manuscripts/CircadianClockTNBC/' # !! specify path for where excel sheet with MRA values is supposed to be stored !!
sheet = "MRA_Values_All_Celllines_CutOff_48h" # !! change to name of the excel sheet for the MRA values !!
Data = pd.read_excel("./Manuscript_RawData_LumiCycle_all_celllines.xlsx")#,sep=",", skipinitialspace=True) #!! put this excel file where the code is !!
#print (Data.columns)

c_label  = ["Noise", "Ultradian", "Circadian", "Infradian"]

workbook=xlsxwriter.Workbook(path+sheet+'.xlsx')
worksheet=workbook.add_worksheet()

for j in range(len(c_label)):
    worksheet.write(0,j+1,c_label[j])

for y in range(len(Data.columns)):
    worksheet.write(y,0,Data.columns[y])
worksheet.write(0,0,'Id') 
    
#@NICA: here I create a second workbook for the detrended data (calculated in line 180)
sheet2 = "Detrended_LumiCycle_All_Celllines_CutOff_48h_All" #excelfile for detrended data

workbook2=xlsxwriter.Workbook(sheet2+'.xlsx')
worksheet2=workbook2.add_worksheet()

# Set pyboat parameters
dt = 10./60         # 10 min in hours
T_cut_off = 48      # cut off period

# ((len(Data.columns))-1): 
    
for ii in range((len(Data.columns))-1): 

    # Filter Data based on time
    Data = Data.loc[Data["Time"] < 121]
    
    t           = Data["Time"].values

    BiolData    = Data.drop("Time", axis=1)

    signal      = BiolData.iloc[:,ii].values
    
    columns = BiolData.shape[1]

    col=ii
    
    #addition of column names (=sample names) to second workbook
    worksheet2.write(0,ii,BiolData.columns[ii])
    
    # CWT parameters
    tested_periods = 200         # increase this to make the period grid 
    periods = np.linspace(2*dt,T_cut_off, tested_periods)
    time_unit = 'h'
    
    #initialize wavelet analyzer #update 15.12.2022
    wAn = pyboat.WAnalyzer(periods, dt, time_unit_label=time_unit) #update 15.12.2022
    
    # # sinc filter detrending
    trend  = wAn.sinc_smooth(signal, T_c=T_cut_off) #update 15.12.2022
    #trend  = pyboat.sinc_smooth(signal, T_cut_off, dt) #old
    detr_signal = signal - trend
    #print(detr_signal)

    for jj in range(len(t)):
        worksheet2.write(jj+1, col, detr_signal[jj])
    
    # plot trend and detrended signal
    plot(t, signal, color="k", linestyle = "-")
    plot(t, detr_signal, color="gray", linestyle = "-")
    plot(t, trend, color="gray", linestyle = "-")
    
    # CWT parameters
    #tested_periods = 200         # increase this to make the period grid 
    #periods = np.linspace(2*dt,T_cut_off, tested_periods)
    #time_unit = 'h'
    
    # wavelet analysis
    #modulus, wlet = pyboat.compute_spectrum(detr_signal, dt, periods) #old
    modulus, wlet = wAn.compute_spectrum(detr_signal) #update 15.12.2022
    
    # get ridge
    #ridge = pyboat.get_maxRidge_ys(modulus) #old
    #ridge_results = pyboat.eval_ridge(ridge, wlet, detr_signal, periods, t) #old
    
    ridge_results = wAn.get_maxRidge() #update 15.12.2022, #uncomment this one when no ridge threshold or smoothing should be used
    #ridge_results = wAn.get_maxRidge(power_thresh=0) #smoothing_wsize=20)  #uncomment this one if you want to set a power threshold filter and/or you want to smooth the ridge
    wAn.draw_Ridge()                    
    
    # Plot ridge
    #print(ridge_results)
    #ridge_results.to_csv("./Ridge_"+str(BiolData.columns[ii])+".csv")        # saves ridge results as csv
    #ridge_results.to_excel("./Ridge_"+str(BiolData.columns[ii])+".xlsx")     # # saves ridge results as xlsx/xls
    
    # plot results
    ax_sig, ax_spec = pl.mk_signal_modulus_ax(time_unit)
    pl.plot_signal_modulus((ax_sig, ax_spec), t, detr_signal, modulus, periods)
    pl.draw_Wavelet_ridge(ax_spec, ridge_results)
     
    #
    # Start DWT
    #
    
    # downsample
    t      = t[::3]
    signal = detr_signal[::3]
    
    # DWT parameters
    c_wavelet = "db20"
    NumLevels = 7
    
    # downsample t
    #t      = t[::3]
    
    #signal = detr_signal
    MaxLength = len(array(t))
    
    LastOne = len(t)
    t = array(t).copy()[:LastOne]
    level = pywt.swt_max_level(LastOne)
    c_modwt = modwt(signal, c_wavelet, NumLevels)
    c_wtmra = modwtmra(c_modwt, c_wavelet)
    
    plotDWT(t, signal, NumLevels, c_wtmra, dt=30./60)
    
    
    fig=plt.figure()
    PlotMRA_CoarseGrained(t, signal, c_wtmra, dt=30./60, NumLevels=7, t_edge=1)
    plt.savefig("MRA_"+str(BiolData.columns[ii])+".svg")
    #plt.savefig("MRA_"+str(BiolData.columns[ii])+".png")
    #plt.savefig("Last_figure_"+name+".svg")
    print(BiolData.columns[ii])
    

#workbook.close()
workbook2.close()

    
    

# TOMAS DAQ data analysis 
# An AA inc. product
import os
import re
import sys
import numpy as np              
import tkfilebrowser
from pathlib import Path
from scipy import spatial
from scipy import optimize
from lvm_read import read
import matplotlib.pyplot as plt  
import matplotlib
import matplotlib.tri as tri
from mpl_smithchart import SmithAxes
from tkinter import filedialog as fd
from tkinter import *
from tkinter import StringVar, OptionMenu
from tkinter.messagebox import showinfo

#matplotlib.use('module://matplotlib-backend-kitty')

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def PlotFile(selection):
    ToRead = ch0filepath
    if not ToRead:
        print("Error, no file selected")
        sys.exit(0)
    if selection == 1:
        #header
        data = read(ToRead)
        keys = data[0]['Channel names']
        for i,key in enumerate(keys):
            x = data[0]['data'][:,0]
            y = data[0]['data'][:,i]
            plt.plot(x,y,label=key)
            plt.xlabel("time (s)")
            plt.ylabel("Voltage")
            plt.legend()
        plt.show()
    
def SelectSignals(ToRead,convert,GasType):
    win = tk.Toplevel()
    win.wm_title("Select Signals")
    # Create a listbox
    listbox = Listbox(win, width=40, height=10, selectmode=MULTIPLE)
     
    # Inserting the listbox items
    if not ToRead:
        print("Error, no file selected")
        sys.exit(0)
    data = read(ToRead)
    keys = data[0]['Channel names']
    for i,key in enumerate(keys):
        if key=="X_Value" or key=="Comment":
                continue
        listbox.insert(i, key)
     
    # Function for printing the
    # selected listbox value(s)
    def selected_item():
        # Traverse the tuple returned by
        # curselection method and print
        # corresponding value(s) in the listbox
        ToPlot = []
        unit = "Volts"
        for i in listbox.curselection():
            ToPlot.append(listbox.get(i))
        lookuptable = np.array(data[0]['Channel names'][:-1])
        plt.clf()
        for key in ToPlot:
            i = np.where(key == lookuptable)[0]
            x = data[0]['data'][:,0]
            y = data[0]['data'][:,i]
            if convert:
                unit,y = Convert(key,y,GasType)
            plt.plot(x,y,label=key)
            plt.xlabel("time (s)")
            plt.ylabel(unit)
            plt.legend()
        plt.show()

     
    # Create a button widget and
    # map the command parameter to
    # selected_item function
    btn = Button(win, text='Plot Selected', command=selected_item)
     
    # Placing the button and listbox
    btn.pack(side='bottom')
    listbox.pack()

def Convert(key,y,GasType):
    unit = "Not known"
    if key == "ECF" or key == "ECR":
        y *= 600
        unit = "Watts"
    if key == "ICF":
        y = 10**(((y - 2.196569)/0.0257915 + 70)/10 - 3) #Using 35MHz calibration
        unit = "Watts"
    if key == "ICR":
        y = 10**(((y-2.188172)/0.02554517 + 70)/10 - 3)
        unit = "Watts"
    if key == "Penning":
        y = 10**((y*1.25) - 9.75)/1000
        if GasType == "H":
            y *= 2.4
        elif GasType == "He":
            y *= 5.9
        elif GasType == "Ar":
            y *= 0.8
        unit = "mbar"
    if key == "Baratron":
        y = (y-1)*0.125*0.1
        unit = "mbar"
    if key == "ICPhase":
        y = 190.31-95.57214*y
        unit = "degrees (°)"

    return unit,y 

def IfCa(FolderLocation):
    DirCas = []
    ICps = []
    ICss = []
    ICas = []
    EICps = []
    EICss = []
    EICas = []

    dirList = os.listdir(FolderLocation)
    # Get Folders location and their names to infer Ca
    for DirName in dirList:
        CaString = ''
        for letter in DirName:
            if letter == 'A':
                continue
            if letter == 'P':
                break
            CaString += letter
        Ca = int(CaString)
        DirCas.append(Ca)
        # Get file inside folder
        FolderPath = FolderLocation + "/" + DirName
        FileNames = os.listdir(FolderPath) #there are pkl and lvm files
        for file in FileNames:
            fileExt = file.split('.')[-1]
            if fileExt == 'LVM':
                FileName = file
        FilePath = Path(FolderPath + "/" + FileName)

        data = read(FilePath)
        i_cp = int(np.where('ICp 1' == np.array(data[0]['Channel names']))[0][0])
        i_ca = int(np.where('ICa 1' == np.array(data[0]['Channel names']))[0][0])
        i_cs = int(np.where('ICs 1' == np.array(data[0]['Channel names']))[0][0])

        ICp = data[0]['data'][:,i_cp]
        ICp = ICp[np.where(ICp>-30)]
        EICp = np.std(ICp[10:-10])
        ICp = ICp.mean()
        ICa = data[0]['data'][:,i_ca]
        ICa = ICa[np.where(ICa>-30)]
        EICa = np.std(ICa[10:-10])
        ICa = ICa.mean()
        ICs = data[0]['data'][:,i_cs]
        ICs = ICs[np.where(ICs>-30)]
        EICs = np.std(ICs[10:-10])
        ICs = ICs.mean()

        x = data[0]['data'][:,0] #Ik denk tijd

        ICps.append(ICp)
        EICps.append(EICp)
        ICss.append(ICs)
        EICss.append(EICs)
        ICas.append(ICa)
        EICas.append(EICa)

    plt.errorbar(DirCas,ICss,EICss,fmt='o',linewidth=1,capsize=6,label='ICs')
    plt.errorbar(DirCas,ICps,EICps,fmt='o',linewidth=1,capsize=6,label='ICp')
    plt.errorbar(DirCas,ICas,EICas,fmt='o',linewidth=1,capsize=6,label='ICa')
    plt.legend()
    plt.title('Measured currents i.f.o prematching Capacitor step value')
    plt.xlabel('Ca (step)')
    plt.ylabel('I (dBm)')
    plt.savefig('figures/IfCa_25MHz.pdf')
    plt.show()

###########################
# Gui for selecting data: #
###########################
import tkinter as tk
from tkinter import filedialog

def open_Ch0_file_dialog():
    global ch0filepath
    ch0filepath = filedialog.askopenfilename(title="Select one DAQ file",initialdir=".." ,filetypes=[("text files", "*.txt"), ("All files", "*.*")])
    if ch0filepath:
        selected_file_ch0_label.config(text=f"Selected File: {ch0filepath}")

def open_Ch1_file_dialog():
    global ch1filepath
    ch1filepath = filedialog.askdirectory(title="select a folder containing multiple DAQ files")
    if ch1filepath:
        selected_file_ch1_label.config(text=f"Selected folder: {ch1filepath}")

global dirs
dirs = []
def get_directories():
    dirs.append(tkfilebrowser.askopendirnames())
    return dirs

def CloseDialog():
    root.destroy()

root = tk.Tk()
root.title("I vs Ca")

top = Frame(root)
bottom = Frame(root)

open_button = tk.Button(root, text="Frequency folder", command=open_Ch1_file_dialog)
open_button.pack(pady=20)

selected_file_ch1_label = tk.Label(root, text="Selected Folder:")
selected_file_ch1_label.pack()

CaTraj = tk.Button(root, text="Plot I vs Ca", command= lambda: IfCa(ch1filepath))

CaTraj.pack(pady=10)

done_button = tk.Button(root, text="Done", command=CloseDialog)
done_button.pack(pady=10)

root.mainloop()

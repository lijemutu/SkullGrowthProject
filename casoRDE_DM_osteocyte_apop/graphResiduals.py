# %%
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import numpy as np

def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(data, key=alphanum_key)

def fieldParser(fileLog,logFolder):
    """
    Given a file log extract the last residual per time step and returns a numpy array
        fileLog must be a string like "DUxFinalRes_","aFinalRes_"
        logFolder must be a string "logs/" or "."
    """
    for logFileS in sorted_alphanumeric(os.listdir(logFolder)):
        if logFileS.startswith(fileLog):
            if logFileS == fileLog+"0":
                logFile = pd.read_csv("{folder}{file}" .format(folder=logFolder,file=logFileS) ,sep='\t',header=None,names=['TimeStep','Residual'])
            #print(logFileS)
            logFileTemporal = pd.read_csv("{folder}{file}" .format(folder=logFolder,file=logFileS) ,sep='\t',header=None,names=['TimeStep','{}' .format(logFileS)])
            logFile = pd.concat([logFile,logFileTemporal],axis=1)
    field = []
    for index,row in logFile.iterrows():
        field.append(row.dropna()[-1])

    return np.array(field)
    
# %%
dux = fieldParser("DUxFinalRes_","logs/")
duy = fieldParser("DUyFinalRes_","logs/")
duz = fieldParser("DUzFinalRes_","logs/")

du = dux +duy +duz /3

a = fieldParser("aFinalRes_","logs/")
h = fieldParser("hFinalRes_","logs/")
o = fieldParser("oFinalRes_","logs/")
b = fieldParser("bFinalRes_","logs/")
ocy = fieldParser("ocyFinalRes_","logs/")

# %%
plt.plot(du,label="DU")
plt.plot(a,label="a")
plt.plot(h,label="h")
plt.plot(o,label="o")
plt.plot(b,label="b")
plt.plot(ocy,label="ocy")
plt.legend()

#plt.ylim(0, 0.01)
plt.show()


# %%

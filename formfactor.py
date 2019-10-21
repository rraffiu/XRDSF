import pandas as pd
from numpy import *
import sys

url ="https://raw.github.com/rraffiu/XRDSF/master/expdata.csv"

df = pd.read_csv(url,index_col=0)
df.index=df.index.str.strip()

def eaff(elem,k):
    try:
        el=df.loc[elem]
    except:
        sys.exit("ERROR: Incorrect symbol or no data available for this element.")
    el = el.str.strip()
    a1=float(el.loc['a1'])
    a2=float(el.loc['a2'])
    a3=float(el.loc['a3'])
    a4=float(el.loc['a4'])
    b1=float(el.loc['b1'])
    b2=float(el.loc['b2'])
    b3=float(el.loc['b3'])
    b4=float(el.loc['b4'])
    c=float(el.loc['c'])
    return a1*exp(-b1*(k/(4*pi))**2)+a2*exp(-b2*(k/(4*pi))**2) \
           +a3*exp(-b3*(k/(4*pi))**2)+a4*exp(-b4*(k/(4*pi))**2)+c

def pdf(k,k0=0.0,sig=0.05):
    return exp(-power((k-k0)/sig,2)/2)

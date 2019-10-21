import pandas as pd
from numpy import *
import sys

urlpc = "https://raw.github.com/rraffiu/XRDSF/master/point.csv"
df = pd.read_csv(urlpc,index_col=2)

urlea ="https://raw.github.com/rraffiu/XRDSF/master/expdata.csv"
df = pd.read_csv(urlea,index_col=0)
df.index=df.index.str.strip()

def eaff(sym,k):
    try:
        el=df.loc[sym]
    except:
        sys.exit("ERROR: Incorrect symbol ("+sym+") or no data available for this element.")
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

def pcff(sym):
    try:
        el=df.loc[sym]
    except:
        sys.exit("ERROR: Incorrect symbol ("+sym+") or no data available for this element.")
    return float(df.loc[sym,'AtomicNumber'])


def gpdf(k,k0=0.0,sig=0.05):
    return exp(-power((k-k0)/sig,2)/2)

def lpdf(k,k0=0.0,sig=0.05):
    return  sig/(2*pi)*(1/((k-k0)**2+(0.5*sig)**2))

# -*- coding: utf-8 -*-
# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
from dash.dependencies import Input, Output
from dash import dash_table
from dash import dcc
from dash import html
import dash_daq as daq
from flask import request
import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import lxml.etree as etree
from xml.dom import minidom
from pathlib import Path
import webbrowser
from tkinter import filedialog
from tkinter import *

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
colors = {'background': '#FFFFFF','text': '#3ead84'}

# File path to run folder
def loadRun():
    root = Tk()
    root.withdraw()
    path = filedialog.askdirectory()
    #path = '/Users/davindersandhu/Desktop/MiSeqDemo'
    #path = '/Users/davindersandhu/Desktop/rmgtest/4/20200827_FS10000333_11_BPA73105-2411'                 #iSeq
    #ath = '/Users/davindersandhu/Desktop/rmgtest/3/200907_NS500523_0453_AHK2TFBGXG' #NextSeq 500
    #path = '/Users/davindersandhu/Desktop/rmgtest/NovaSeq_S4_Illumina_DNA_PCR-free_v1.5_chemistry_-195476289'   #NovaSeq
    # Execution command on path
    cmd = '/Users/davindersandhu/PycharmProjects/sav4mac/interop-1.3.2-Darwin-AppleClang/bin/dumptext '+path
    so = os.popen(cmd).read().strip().split('\n')
    cmd1 = '/Users/davindersandhu/PycharmProjects/sav4mac/interop-1.3.2-Darwin-AppleClang/bin/summary '+path
    #print(path, cmd)
    so1 = os.popen(cmd1).read().strip().split('\n')
    cmd2 = '/Users/davindersandhu/PycharmProjects/sav4mac/interop-1.3.2-Darwin-AppleClang/bin/imaging_table '+path
    so2 = os.popen(cmd2).read().strip().split('\n')
    return (so,so1,so2,path)
(so,so1,so2,path)=loadRun()

def rp(path):
    path = path+'/RunParameters.xml'
    myfile = Path(path)
    if myfile.is_file():
        tree = etree.parse(path)
        #root = tree.getroot()
        root = tree.findall('*')
        #print(root)
        output=[]
        for child in root:
            output.append(str(child.tag)+': '+str(child.text))
        return(output)
    else:
        return([])

def readlengths(path):
    xmldoc = minidom.parse(path + '/RunInfo.xml')
    itemlist = xmldoc.getElementsByTagName('Read')
    itemlist.pop()
    returnlist = []
    for s in itemlist:
        returnlist.append(s.attributes['NumCycles'].value)
    x = []
    for i in range(0, len(returnlist)):
        if i == 0:
            x.append(int(returnlist[i]))
        else:
            x.append(int(returnlist[i]) + int(x[i - 1]))
    return(x)

def q30bycycle(so):

    #Q30 by cycle plot generation
    start = 0
    end = 0
    for i in range(0,len(so)):
        if so[i] == '# Q2030,1':
            start=i
            break
    for i in range(start+5,len(so)):
        if so[i][0]=='#':
            end = i
            break


    for i in range(start,end):
        if so[i][0]=="L":
            start=i
            break

    y=pd.DataFrame(so[start+1:end])
    t = (y[0].str.split(",",expand=True))
    t.columns=so[start].split(",")
    t=t.astype('int64')
    t['Q30Percent'] = 100* t['Q30']/t['Total']
    t['Surface'] =((t.Tile <2000)).map({True:'Top',False:'Bottom'})
    med = t.groupby(['Cycle'])[['Q30Percent']].mean()
    med.name='Median'
    t = t.merge(med,on='Cycle')
    t = t.sort_values('Cycle')
    #print(t)
    fig = go.Figure()
    hlines = readlengths(path)
    for i in hlines:
        fig.add_vline(x=i, line_width=1, line_dash="dash", line_color="black")
    fig.add_trace(go.Box(x=t.Cycle, y=t.Q30Percent_x,line=dict(color='blue')))
    fig.add_trace(go.Scatter(x=t.Cycle,y=t.Q30Percent_y,name='Mean Q30', mode='lines'))
    fig.update_layout(showlegend=False, title='Q30 by Cycle', xaxis={'title':'Cycle'}, yaxis={'title':'%>= Q30'})
    fig.update_xaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_yaxes(rangemode="tozero")
    return fig

def percbasebycycle(so):
    #cmd = '/Users/davindersandhu/Desktop/InterOp-1.1.15-Darwin-AppleClang/bin/dumptext '+ path
    #so = os.popen(cmd).read().strip().split('\n')

    start = 0
    end = 0

    for i in range(0,len(so)):
        if so[i] == '# CorrectedInt,1':
            start=i
            break
    for i in range(start+5,len(so)):
        if so[i][0]=='#':
            end = i
            break


    for i in range(start,end):
        if so[i][0]=="L":
            start=i
            break

    y=pd.DataFrame(so[start+1:end])
    t = (y[0].str.split(",",expand=True))
    t.columns=so[start].split(",")
    t=t.astype('float64')
    t = t.sort_values(['Cycle', 'Tile'])
    #t.to_csv('test0.csv', sep='\t')
    #print(t.columns)
    t['CountTotal']=t['CalledCount_A']+t['CalledCount_C']+t['CalledCount_G']+t['CalledCount_T']
    t['A_perc']=100*t['CalledCount_A']/t['CountTotal']
    t['C_perc']=100*t['CalledCount_C']/t['CountTotal']
    t['G_perc']=100*t['CalledCount_G']/t['CountTotal']
    t['T_perc']=100*t['CalledCount_T']/t['CountTotal']


    #print(t.columns)
    t['Surface'] =((t.Tile <2000)).map({True:'Top',False:'Bottom'})

    medPA = t.groupby(['Cycle'])[['A_perc']].mean()
    medPA.name='MedianPA'
    #print(medPA)
    t= t.merge(medPA,on='Cycle')

    medPC = t.groupby(['Cycle'])[['C_perc']].mean()
    medPC.name='MedianPC'
    t= t.merge(medPC,on='Cycle')

    #print(t.groupby(['Cycle'])[['G_perc']].mean())
    #print(t.Cycle, type(t.Cycle))
    #print(t.G_perc, type(t.G_perc))
    medPG = t.groupby(['Cycle'])[['G_perc']].mean()
    medPG.name='MedianPG'
    #print(medPG)
    t= t.merge(medPG,on='Cycle')

    medPT = t.groupby(['Cycle'])[['T_perc']].mean()
    medPT.name='MedianPT'
    t = t.merge(medPT,on='Cycle')
    #print(t)
    fig= go.Figure()
    t = t.sort_values(['Cycle','Tile'])
    #print(t)
    #print(t.columns)
    #t.to_csv('test1.csv',sep='\t')
    hlines = readlengths(path)
    for i in hlines:
        fig.add_vline(x=i, line_width=1, line_dash="dash", line_color="black")
    fig.add_trace(go.Scatter(x=t.Cycle,y=t.A_perc_y, name='A', mode='lines',line=dict(color='firebrick')))
    fig.add_trace(go.Scatter(x=t.Cycle,y=t.C_perc_y, name='C', mode='lines',line=dict(color='forestgreen')))
    fig.add_trace(go.Scatter(x=t.Cycle,y=t.G_perc_y, name='G', mode='lines',line=dict(color='royalblue')))
    fig.add_trace(go.Scatter(x=t.Cycle,y=t.T_perc_y, name='T', mode='lines',line=dict(color='black')))
    fig.update_layout(showlegend=True, title='%Base by Cycle', xaxis={'title':'Cycle'}, yaxis={'title':'% Base'})
    fig.update_xaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    return fig

def fwhmbycycle(so):
    start = 0
    end = 0
    for i in range(0,len(so)):
        if so[i] == '# Extraction,1':
            start=i
            break
    for i in range(start+5,len(so)):
        if so[i][0]=='#':
            end = i
            break

    for i in range(start,end):
        if so[i][0]=="L":
            start=i
            break

    y=pd.DataFrame(so[start+1:end])
    t = (y[0].str.split(",",expand=True))

    t.columns=so[start].split(",")
    t=t.astype('float64')
    t = t.sort_values("Cycle")
    fig=go.Figure()
    if len(t.columns)==12:
        t.Focus_A=pd.to_numeric(t.Focus_A)
        t.Focus_C=pd.to_numeric(t.Focus_C)
        t.Focus_G=pd.to_numeric(t.Focus_G)
        t.Focus_T=pd.to_numeric(t.Focus_T)

        medFA = t.groupby(['Cycle'])[['Focus_A']].apply(np.average)
        medFA.name='MedianFA'
        t= t.join(medFA,on=['Cycle'])

        medFC = t.groupby(['Cycle'])[['Focus_C']].apply(np.average)
        medFC.name='MedianFC'
        t= t.join(medFC,on=['Cycle'])

        medFG = t.groupby(['Cycle'])[['Focus_G']].apply(np.average)
        medFG.name='MedianFG'
        t= t.join(medFG,on=['Cycle'])

        medFT = t.groupby(['Cycle'])[['Focus_T']].apply(np.average)
        medFT.name='MedianFT'
        t= t.join(medFT,on=['Cycle'])

        fig.add_trace(go.Scatter(x=t.Cycle,y=t.MedianFA,name='A',mode='lines',line=dict(color='firebrick')))
        fig.add_trace(go.Scatter(x=t.Cycle,y=t.MedianFC,name='C',mode='lines',line=dict(color='forestgreen')))
        fig.add_trace(go.Scatter(x=t.Cycle,y=t.MedianFG,name='G',mode='lines',line=dict(color='royalblue')))
        fig.add_trace(go.Scatter(x=t.Cycle,y=t.MedianFT,name='T',mode='lines',line=dict(color='black')))

    elif len(t.columns)==8 and t.columns[6]=='Focus_RED':
        t.Focus_RED=pd.to_numeric(t.Focus_RED)
        t.Focus_GREEN=pd.to_numeric(t.Focus_GREEN)

        medFR = t.groupby(['Cycle'])[['Focus_RED']].apply(np.average)
        medFR.name='MedianFR'
        t= t.join(medFR,on=['Cycle'])

        medFG = t.groupby(['Cycle'])[['Focus_GREEN']].apply(np.average)
        medFG.name='MedianFG'
        t= t.join(medFG,on=['Cycle'])
        fig.add_trace(go.Scatter(x=t.Cycle,y=t.MedianFR,name='Red',mode='lines',line=dict(color='red')))
        fig.add_trace(go.Scatter(x=t.Cycle,y=t.MedianFG,name='Green',mode='lines',line=dict(color='green')))

    elif len(t.columns)==8 and t.columns[6]=='Focus_Red':
        t.Focus_Red=pd.to_numeric(t.Focus_Red)
        t.Focus_Green=pd.to_numeric(t.Focus_Green)

        medFR = t.groupby(['Cycle'])[['Focus_Red']].apply(np.average)
        medFR.name='MedianFR'
        t= t.join(medFR,on=['Cycle'])

        medFG = t.groupby(['Cycle'])[['Focus_Green']].apply(np.average)
        medFG.name='MedianFG'
        t= t.join(medFG,on=['Cycle'])
        fig.add_trace(go.Scatter(x=t.Cycle,y=t.MedianFR,name='Red',mode='lines',line=dict(color='red')))
        fig.add_trace(go.Scatter(x=t.Cycle,y=t.MedianFG,name='Green',mode='lines',line=dict(color='green')))


    if len(t.columns)==8 and t.columns[6]=='Focus_1':
        t.Focus_1=pd.to_numeric(t.Focus_1)
        t.Focus_2=pd.to_numeric(t.Focus_2)

        medF1 = t.groupby(['Cycle'])[['Focus_1']].apply(np.average)
        medF1.name='MedianF1'
        t= t.join(medF1,on=['Cycle'])

        medF2 = t.groupby(['Cycle'])[['Focus_2']].apply(np.average)
        medF2.name='MedianF2'
        t= t.join(medF2,on=['Cycle'])


        fig.add_trace(go.Scatter(x=t.Cycle,y=t.MedianF1,name='Red',mode='lines',line=dict(color='red')))
        fig.add_trace(go.Scatter(x=t.Cycle,y=t.MedianF2,name='Green',mode='lines',line=dict(color='green')))
    hlines = readlengths(path)
    for i in hlines:
        fig.add_vline(x=i, line_width=1, line_dash="dash", line_color="black")
    fig.update_layout(showlegend=True, title='FWHM by Cycle', xaxis={'title':'Cycle'}, yaxis={'title':'FWHM'})
    fig.update_layout(yaxis=dict(range=[0,6]))
    fig.update_xaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    return fig

# so = os.popen(cmd).read().strip().split('\n')
#
#
def intensitycycle(so):
    # Intensity by cycle
    start3 = 0
    end3 = 0
    for i in range(0,len(so)):
        if so[i] == '# Extraction,1':
            start3=i
    for i in range(start3+5, len(so)):
        if so[i][0]=='#':
            #print(so[i], so[i][0], i)
            end3 = i
            break
    #
    x3 = so[start3+4:end3]
    #x3.to_csv('x3.csv',sep='\t')
    y3=pd.DataFrame(x3)

    t3 = (y3[0].str.split(",",expand=True))

    #t3.to_csv('t4.csv', sep='\t')
    #print(so[start3+3].split(","))
    colList=(so[start3 + 3].split(","))
    t3.columns=colList

    t3=t3.astype('float64')
    if len(t3.columns)==8:
        t3.columns = map(str.lower, t3.columns)
        if 'maxintensity_red' in t3.columns:
            meanRed = t3.groupby(['cycle'])[['maxintensity_red']].mean()
            t3 = t3.merge(meanRed,on='cycle')

            meanGreen = t3.groupby(['cycle'])[['maxintensity_green']].mean()
            t3 = t3.merge(meanGreen,on='cycle')

            fig3 = go.Figure()
            fig3.add_trace(
                go.Scatter(x=t3.cycle, y=t3.maxintensity_red_y, mode='lines', name='Red', line=dict(color='firebrick')))
            fig3.add_trace(
                go.Scatter(x=t3.cycle, y=t3.maxintensity_green_y, mode='lines', name='Green', line=dict(color='forestgreen')))
            hlines = readlengths(path)
            for i in hlines:
                fig3.add_vline(x=i, line_width=1, line_dash="dash", line_color="black")
            fig3.update_layout(showlegend=True, title='Intensity by Cycle', xaxis={'title': 'Cycle'},
                               yaxis={'title': 'Intensity'})
            fig3.update_xaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
            fig3.update_yaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
            fig3.update_yaxes(rangemode="tozero")
            return (fig3)
        elif 'maxintensity_1' in t3.columns:
            meanRed = t3.groupby(['cycle'])[['maxintensity_1']].mean()
            t3 = t3.merge(meanRed, on='cycle')

            meanGreen = t3.groupby(['cycle'])[['maxintensity_2']].mean()
            t3 = t3.merge(meanGreen, on='cycle')

            fig3 = go.Figure()
            fig3.add_trace(
                go.Scatter(x=t3.cycle, y=t3.maxintensity_1_y, mode='lines', name='Red', line=dict(color='firebrick')))
            fig3.add_trace(
                go.Scatter(x=t3.cycle, y=t3.maxintensity_2_y, mode='lines', name='Green',
                           line=dict(color='forestgreen')))
            hlines = readlengths(path)
            for i in hlines:
                fig3.add_vline(x=i, line_width=1, line_dash="dash", line_color="black")
            fig3.update_layout(showlegend=True, title='Intensity by Cycle', xaxis={'title': 'Cycle'},
                               yaxis={'title': 'Intensity'})
            fig3.update_xaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
            fig3.update_yaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
            fig3.update_yaxes(rangemode="tozero")
            return (fig3)

    else:
        medA = t3.groupby(['Cycle'])[['MaxIntensity_A']].mean()
        medA.name='MedianA'
        t3= t3.merge(medA,on='Cycle')
        #
        medC = t3.groupby(['Cycle'])[['MaxIntensity_C']].mean()
        medC.name='MedianC'
        t3= t3.merge(medC,on='Cycle')
        #
        medG = t3.groupby(['Cycle'])[['MaxIntensity_G']].mean()
        medG.name='MedianG'
        t3= t3.merge(medG,on='Cycle')
        #
        medT = t3.groupby(['Cycle'])[['MaxIntensity_T']].mean()
        medT.name='MedianT'
        t3= t3.merge(medT,on='Cycle')
        #
        #print(t3.columns)
        fig3 = go.Figure()
        fig3.add_trace(go.Scatter(x=t3.Cycle,y=t3.MaxIntensity_A_y, mode='lines',name='A',line=dict(color='firebrick')))
        fig3.add_trace(go.Scatter(x=t3.Cycle,y=t3.MaxIntensity_C_y, mode='lines',name='C',line=dict(color='forestgreen')))
        fig3.add_trace(go.Scatter(x=t3.Cycle,y=t3.MaxIntensity_G_y, mode='lines',name='G',line=dict(color='royalblue')))
        fig3.add_trace(go.Scatter(x=t3.Cycle,y=t3.MaxIntensity_T_y, mode='lines',name='T',line=dict(color='black')))
        hlines = readlengths(path)
        for i in hlines:
            fig3.add_vline(x=i, line_width=1, line_dash="dash", line_color="black")
        fig3.update_layout(showlegend=True, title='Intensity by Cycle', xaxis={'title':'Cycle'}, yaxis={'title':'Intensity'})
        fig3.update_xaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
        fig3.update_yaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
        fig3.update_yaxes(rangemode="tozero")
        return (fig3)

def errbycycle(so):

    start = 0
    error = 0
    end = 0
    for i in range(0,len(so)):
        if so[i].find('# Error') !=-1:
            start=i
            error = 1
            break
    for i in range(start+5,len(so)):
        if so[i][0]=='#':
            end = i
            break


    for i in range(start,end):
        if so[i][0]=="L":
            start=i
            break
    y=pd.DataFrame(so[start+1:end])
    t = (y[0].str.split(",",expand=True))
    t.columns=so[start].split(",")
    t=t.astype('float64')
    #print(start, t.columns)
    fig = go.Figure()
    hlines = readlengths(path)
    for i in hlines:
        fig.add_vline(x=i, line_width=1, line_dash="dash", line_color="black")
    if error == 1:
        fig.add_trace(go.Box(x=t.Cycle, y=t.ErrorRate,line=dict(color='blue')))
    else:
        fig.add_trace(go.Box(x=t.Cycle, y=[0],line=dict(color='blue')))
    fig.update_layout(showlegend=False, title='Error Rate by Cycle', xaxis={'title':'Cycle'}, yaxis={'title':'% Error Rate'})
    fig.update_xaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_yaxes(rangemode="tozero")
    return fig

def snrbycycle(so):

    start = 0
    end = 0
    for i in range(0,len(so)):
        if so[i].find('# CorrectedInt') !=-1:
            start=i
            break
    for i in range(start+5,len(so)):
        if so[i][0]=='#':
            end = i
            break


    for i in range(start,end):
        if so[i][0]=="L":
            start=i
            break

    y=pd.DataFrame(so[start+1:end])
    t = (y[0].str.split(",",expand=True))
    t.columns=so[start].split(",")
    t=t.astype('float64')
    #print(start, t.columns)
    fig = go.Figure()
    hlines = readlengths(path)
    for i in hlines:
        fig.add_vline(x=i, line_width=1, line_dash="dash", line_color="black")
    fig.add_trace(go.Box(x=t.Cycle, y=t.SignalToNoise,line=dict(color='blue')))
    fig.update_layout(showlegend=False, title='SNR by Cycle', xaxis={'title':'Cycle'}, yaxis={'title':'Signal to Noise Ratio'})
    fig.update_xaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_yaxes(rangemode="tozero")
    return fig

#####################################
######### Data by lane ##############
#####################################

def clusterlane(so2):
    start = 2
    end = len(so2)
    datacd = so2[3:end]
    y = pd.DataFrame(datacd)
    t = (y[0].str.split(",", expand=True))
    t1 = t.iloc[:, :10]
    t1.columns = so2[start].split(",")[:10]
    t1 = t1.astype('float64')
    rawcd = t1.groupby(['Lane'])[['Density(k/mm2)']].mean()
    t1 = t1.merge(rawcd, on='Lane')
    pfcd = t1.groupby(['Lane'])[['Density Pf(k/mm2)']].mean()
    t1 = t1.merge(pfcd, on='Lane')
    fig = go.Figure()
    fig.add_trace(go.Box(x=t1.Lane, y=t1['Density(k/mm2)_x'], line=dict(color='blue'),name='Raw Cluster Density'))
    fig.add_trace(go.Box(x=t1.Lane, y=t1['Density Pf(k/mm2)_x'], line=dict(color='green'), name = 'Cluster Density PF'))
    fig.update_layout(showlegend=False, title='Cluster Density by Lane', xaxis={'title':'Lane'}, yaxis={'title':'Cluster Density (k/mm^2)'})
    fig.update_xaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_yaxes(rangemode="tozero")
    return(fig)

def percentlane(so2):
    start = 2
    end = len(so2)
    datacd = so2[3:end]
    y = pd.DataFrame(datacd)
    t = (y[0].str.split(",", expand=True))
    t1 = t.iloc[:, :10]
    t1.columns = so2[start].split(",")[:10]
    t1 = t1.astype('float64')
    percpf = t1.groupby(['Lane'])[['% Pass Filter']].mean()
    t1 = t1.merge(percpf, on='Lane')
    fig = go.Figure()
    fig.add_trace(go.Box(x=t1.Lane, y=t1['% Pass Filter_x'], line=dict(color='green'), name = '% Pass Filter'))
    fig.update_layout(showlegend=False, title='% Pass Filter by Lane', xaxis={'title':'Lane'}, yaxis={'title':'% Pass Filter'})
    fig.update_xaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_layout(yaxis_range=[0,100])
    return(fig)

def alignedlane(so2):
    start = 2
    end = len(so2)
    datacd = so2[3:end]
    y = pd.DataFrame(datacd)
    t = (y[0].str.split(",", expand=True))
    t1 = t.iloc[:, :11]
    t1.columns = so2[start].split(",")[:11]
    print(t1.columns)
    print(t1['% Aligned'][0])
    for i in t1['% Aligned']:
        #print(t1['% Aligned'][i])
        if [i]=='nan':
            t1['% Aligned'][i]=0

    t1 = t1.astype('float64')
    percpf = t1.groupby(['Lane'])[['% Aligned']].mean()
    t1 = t1.merge(percpf, on='Lane')
    fig = go.Figure()
    fig.add_trace(go.Box(x=t1.Lane, y=t1['% Aligned_x'], line=dict(color='green'), name = '% Pass Filter'))
    fig.update_layout(showlegend=False, title='% Aligned by Lane', xaxis={'title':'Lane'}, yaxis={'title':'% Aligned to PhiX'})
    fig.update_xaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_layout(yaxis_range=[0,100])
    return(fig)

def prephasingLane(so):
    start = 0
    end = 0
    for i in range(0, len(so)):
        if so[i].find('# Tile') != -1:
            start = i
            break
    for i in range(start + 5, len(so)):
        if so[i][0] == '#':
            end = i
            break

    for i in range(start, end):
        if so[i][0] == "L":
            start = i
            break

    if end ==0:
        end = len(so)
    print(so[start])
    y = pd.DataFrame(so[start + 3:end])
    print(y)
    t = (y[0].str.split(",", expand=True))
    print(type(t))
    t.names = so[start+2].split()
    print(t)
    t = t.astype('float64')
    # print(start, t.columns)
    fig = go.Figure()
    fig.add_trace(go.Box(x=t[0], y=t[8], line=dict(color='green')))
    fig.update_layout(showlegend=False, title='Prephasing by Lane', xaxis={'title': 'Lane'},
                      yaxis={'title': 'Legacy Prephasing Rate'})
    fig.update_xaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_layout(yaxis_range=[0,1])
    return fig

def phasingLane(so):
    start = 0
    end = 0
    for i in range(0, len(so)):
        if so[i].find('# Tile') != -1:
            start = i
            break
    for i in range(start + 5, len(so)):
        if so[i][0] == '#':
            end = i
            break

    for i in range(start, end):
        if so[i][0] == "L":
            start = i
            break

    if end ==0:
        end = len(so)
    print(so[start])
    y = pd.DataFrame(so[start + 3:end])
    print(y)
    t = (y[0].str.split(",", expand=True))
    print(type(t))
    t.names = so[start+2].split()
    print(t)
    t = t.astype('float64')
    # print(start, t.columns)
    fig = go.Figure()
    fig.add_trace(go.Box(x=t[0], y=t[9], line=dict(color='green')))
    fig.update_layout(showlegend=False, title='Phasing by Lane', xaxis={'title': 'Lane'},
                      yaxis={'title': 'Legacy Phasing Rate'})
    fig.update_xaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
    fig.update_layout(yaxis_range=[0,1])
    return fig

def occupiedlane(so2):
    start = 2
    end = len(so2)
    datacd = so2[3:end]
    y = pd.DataFrame(datacd)
    t = (y[0].str.split(",", expand=True))
    t1 = t.iloc[:, :11]
    listcol = so2[start].split(",")[:11]
    print(t1.columns)
    if "% Occupied" in listcol:
        #print(t1['% Aligned'][0])
        t1 = t1.astype('float64')
        fig = go.Figure()
        fig.add_trace(go.Box(x=t1.iloc[:,0], y=t1.iloc[: , -1], line=dict(color='green'), name='% Pass Filter'))
        fig.update_layout(showlegend=False, title='% Occupied by Lane', xaxis={'title': 'Lane'},
                          yaxis={'title': '% Occupied'})
        fig.update_xaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
        fig.update_yaxes(showline=True, linewidth=1, linecolor='grey', mirror=True)
        fig.update_yaxes(rangemode="tozero")
        return (fig)

def summary(so1):

    #cmd = '/Users/davindersandhu/Desktop/InterOp-1.1.15-Darwin-AppleClang/bin/summary '+path
    #print(path, cmd)
    #so = os.popen(cmd).read().strip().split('\n')
    t1s=0
    t1e=0
    #print(len(so))
    for i in range(0,len(so1)):
        #print(so[i])
        if len(so1[i])>0 and so1[i][0]=='L':
            t1s=i
            break
    for i in range(t1s,len(so1)):
        if len(so1[i]) and so1[i][0]=='T':
            t1e=i
            break
    #print(t1s,t1e)
    y1 = pd.DataFrame(so1[t1s+1:t1e])
    #print(y1)
    t1=y1[0].str.split(",",expand=True)
    t1.columns=so1[t1s].split(",")

    t2s=t1e
    t2e=t1e
    for i in range(t1e,len(so1)):
        if len(so1[i])>0 and so1[i][0]=='R':
            t2s=i
            break
    for i in range(t2s,len(so1)):
        if len(so1[i])>0 and so1[i][0]=='E':
            t2e=i
            break
    y2 = pd.DataFrame(so1[t2s:t2e])
    t2=y2[0].str.split(",",expand=True)
    t2.columns=so1[t2s+1].split(",")
    #print(t2)

    return(t1,t2)

@app.callback(
    Output('stop-button-output-1', 'children'),
    Input('my-stop-button-1', 'n_clicks')
)
def update_output(n_clicks):
    if n_clicks >0:
        func = request.environ.get('werkzeug.server.shutdown')
        if func is None:
            raise RuntimeError('Not running with the Werkzeug Server')
        func()



# Dash
app.title = 'webSAVvy'
app.layout = html.Div(children=[
    daq.StopButton(
        id='my-stop-button-1',
        n_clicks=0
    ),
    html.Div(id='stop-button-output-1'),
    html.Div([
        html.H1(
            children='webSAV',
            style={
                'textAlign': 'center',
                'color': colors['text']
            }
        ),
        html.Div(
            children="a browser-based approach to Sequencing Analysis Viewer" ,
            style={
                'textAlign': 'center',
                'color': 'black'
            }
        ),
    html.Div(id='output-data-upload'),
        html.Hr(
            style={
                'height':'2px'
            }
        ),
        html.Div(
            children="Current folder: "+path,
            style={
                'textAlign': 'left',
                'color': 'black'
            }
        ),
        html.A(html.Button('Load run folder'),href='/'),
        #html.Button('Load run folder', id='load-run'),
        html.Hr(
            style={
                'height':'2px'
            }
        )
    ]),
    html.Div([
        dcc.Tabs([
            dcc.Tab(label='Data by Cycle', children=[
                html.Div([
                    html.Div([
                        dcc.Graph(
                            id='intensitybycycle',
                            figure=intensitycycle(so),
                            config={'displaylogo': False}
                        )
                    ], className='six columns'),
                    html.Div([
                        dcc.Graph(
                            id='%Q30bycycle',
                            figure=q30bycycle(so),
                            config={'displaylogo':False}
                        )
                    ],className='six columns')]),
                html.Div([
                    html.Div([
                        dcc.Graph(
                            id='%Basebycycle',
                            figure=percbasebycycle(so),
                            config={'displaylogo':False}
                        )
                    ],className='six columns'),
                    html.Div([
                        dcc.Graph(
                            id='FWHMbycycle',
                            figure=fwhmbycycle(so),
                            config={'displaylogo':False}
                        )
                    ],className='six columns')
                ]),
                html.Div([
                    html.Div([
                        dcc.Graph(
                            id='Errbycycle',
                            figure=errbycycle(so),
                            config={'displaylogo':False}
                        )
                    ],className='six columns'),
                    html.Div([
                            dcc.Graph(
                            id='SNRbycycle',
                            figure=snrbycycle(so),
                            config={'displaylogo':False}
                        )
                    ],className='six columns')
                ])
            ]),
            #Data by lane tab
            dcc.Tab(label='Data by Lane', children=[
                html.Div([
                    html.Div([
                        dcc.Graph(
                            id='ClusterLane',
                            figure=clusterlane(so2),
                            config={'displaylogo':False}
                        )
                    ], className='six columns'),
                    html.Div([
                        dcc.Graph(
                            id='pfLane',
                            figure=percentlane(so2),
                            config={'displaylogo':False}
                        )
                    ],className='six columns')]),
                html.Div([
                    html.Div([
                        dcc.Graph(
                            id='alignedLane',
                            figure=alignedlane(so2),
                            config={'displaylogo':False}
                        )
                    ],className='six columns'),
                    html.Div([
                        dcc.Graph(
                            id='phasingLane',
                            figure=phasingLane(so),
                            config={'displaylogo':False}
                        )
                    ],className='six columns')
                ]),
                html.Div([
                    html.Div([
                        dcc.Graph(
                            id='prephasingLane',
                            figure=prephasingLane(so),
                            config={'displaylogo':False}
                        )
                    ],className='six columns'),
                    html.Div([
                        dcc.Graph(
                            id='occLane',
                            figure=occupiedlane(so2),
                            config={'displaylogo':False}
                        )

                    ],className='six columns')
                ])
            ]),
            #dcc.Tab(label='Flow Cell Chart', children=[
            #    html.Div([
            #        html.Div([

            #        ], className='six columns'),
            #        html.Div([

            #       ],className='six columns')]),
            #   html.Div([
            #        html.Div([

            #        ],className='six columns'),
            #        html.Div([

            #        ],className='six columns')
            #    ])
            #]),
            dcc.Tab(label='Summary', children=[
                html.Hr(
                    style={
                      'height':'2px'
                    }
                ),
                html.Div([
                    dash_table.DataTable(
                        id='table',
                        columns=[{"name": i, "id": i} for i in summary(so1)[0].columns],
                        data=summary(so1)[0].to_dict('records'),
                    ),

                ]),
                html.Hr(
                    style={
                      'height':'2px'
                    }
                ),
                html.Div([
                    dash_table.DataTable(
                        id='table2',
                        columns=[{"name": i, "id": i} for i in summary(so1)[1].columns],
                        data=summary(so1)[1].to_dict('records'),
                        style_header = {'display': 'none'}
                    )
                ])
            ])
        ])
    ])
])
webbrowser.open('http://127.0.0.1:8050/', new=2)

if __name__ == '__main__':
    app.run_server(debug=False)


# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 13:21:22 2020

@author: mu_hp
"""
import numpy as np
import csv
import pandas as pd
from sympy import *

df = pd.read_csv (r'Dados.csv')
print (df)

precipt = df["PrecipitacaoTotal"].tolist()
mes = df['Mês'].tolist()

print(precipt)
print('\n\n\n')
print(mes)

# Para prever a função a um valor específico, comentar a proxima linha e descomentar a debaixo
#x = symbols("x")
x = 71

coeficientes = []
pontos = 60

for indice in range(pontos):
    L  = 1
    L1 = 1
    L2 = 1
    for j in range(pontos):
        if indice != j:
            L1 *= (x - mes[j])
            L2 *= (mes[indice]-mes[j])
            L = L1/L2
    coeficientes.append(L)

print ('\n\n\n')
#print ('L59 =', coeficientes[58])

pn = 0

for i in range(len(coeficientes)):
    pn += precipt[i]*coeficientes[i]

print ('P60(x) = ', pn)
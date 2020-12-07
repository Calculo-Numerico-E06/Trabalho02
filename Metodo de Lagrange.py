# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 13:21:22 2020

Autores: Felipe Murakami Caldas Tourinho; Gabriel Biscardi; Murilo Henrique Pasini Trevisan
"""
import numpy as np
import csv
import pandas as pd
from sympy import *
import sympy


# Leitura dos valores do arquivo csv e montagem da tabela
df = pd.read_csv (r'Dados.csv')
print (df)

# Conversão das colunas da tabela em listas indexadas
precipt = df["PrecipitacaoTotal"].tolist()
mes = df['Mês'].tolist()

# Visualização dos resultados das listas
#print(precipt)
#print('\n\n\n')
#print(mes)

"""

 Metodo de Lagrange

"""

# Para prever a função a um valor específico, comentar a proxima linha e descomentar a debaixo
x = symbols("x")
#x = 7

# Coeficientes Lk do metodo de lagrange
coeficientes = []

# Numero de pontos a serem utilizados da tabela
pontos = 60

# Laço de leitura dos valores da coluna dos meses, e formação dos coeficientes Lk
for indice in range(pontos):
    L  = 1
    L1 = 1
    L2 = 1
    for j in range(pontos):
        if indice != j:
            L1 *= (x - mes[j])
            L2 *= (mes[indice]-mes[j])
            Lk = L1/L2
    coeficientes.append(Lk)

# Impressão dos primeiros coeficientes para verificação
#print ('\n\n\n')
#print ('L0 =', coeficientes[0])
#print ('L1 =', coeficientes[1])
#print ('L2 =', coeficientes[2])

# Laço de repetição para realizar o somatório dos Lk com seus respectivos Yk
pnl = 0
for k in range(len(coeficientes)):
    pnl += precipt[k]*coeficientes[k]

# Impressão do polinomio de interpolacao de Lagrange, utilizando a biblioteca sympy para simplificar os termos em x
print ('\n\n\nPolinomio de interpolacao de Lagrange')
print ('\nP(',str(pontos-1),')=(x) = ', sympy.expand(pnl))

"""

Metodo dos minimos quadrados

"""
"""
F(x) = C0 + C1*x + C2*x^2 + ...
"""
# Montagem do vetor de resultados ( usando para poder usar o .dot da matriz sem precisar do precipt todo)
Y = np.zeros(pontos)
for i in range(pontos):
    Y[i] = precipt[i]

# Definicao do grau do polinomio de minimos quadrados
grau = 30

# Criando os vetores Hk da matriz
Hk = np.zeros((grau+1, pontos))
for i in range(grau+1):
    for j in range(pontos):
        Hk[i][j] = pow(mes[j], i)

#print ('\n\n\nMatriz H',Hk)

# Montando o sistema Ax = b para solução dos coeficientes Ck
A = np.zeros((grau+1, grau+1))
b = np.zeros(grau+1)
for i in range(grau+1):
    for j in range(grau+1):
        A[i][j] = Hk[i].dot(Hk[j])
    b[i] = Hk[i].dot(Y)
    
#print('\n\n\nMatriz A',A)
#print('\n\n\nMatriz b',b)

# Resolvendo o sistema Ax = b:
Pnmmq = np.linalg.solve(A,b)

#print('\n\n\n',Pnmmq)

#for i in range(len(Pnmmq)):
#    print("\nC" + str(i) + "=", Pnmmq[i])

#Montagem do polinomio a partir dos coeficientes encontrados
F = 0
for i in range(len(Pnmmq)):
    F += Pnmmq[i]*pow(x, i)
    
# Impressao do polinomio no terminal
print ('\n\n\nPolinomio do minimos quadrados de grau =', grau)
print ('\nP'+ str(grau) + '(x)=', F)
    

"""

Metodo dos minimos quadrados trigonometrico

"""
"""
F(x) = a0*cos(0x) + b0*sen(0x) + a1*cos(1x) + b1*sen(1x) + ... + ak*cos(kx) + bk*sen(kx) 

"""

# Os valores Yk se mantém os mesmos do anterior assim como o grau do polinomio e o numero de pontos

ordem = 60

# Uma matriz diagonal que terá ordem*2 + 1 linhas e colunas, a diagonal é preenchida com N para o primeiro e N/2 no restante
D = np.zeros((2*ordem + 1, 2* ordem + 1))
for i in range(2*ordem+1):
    for j in range(2*ordem+1):
        if i == 0 and j == 0:
            D[i][j] = pontos
        elif i == j:
            D[i][j] = pontos/2
        else:
            D[i][j] = 0

#print ('\n\n\nD =',D)

# Divisao da amostra em um intervalo de 0 a 2pi
Xt = np.zeros(pontos)
for i in range(pontos):
    Xt[i] = i*2*(np.pi)/pontos

# Vetor de senos e cossenos de x multiplicado pelo indice
V = np.zeros((pontos, 2*ordem+1))
for j in range(2*ordem+1):
    for i in range(pontos):
        if j == 0:
            V[i][j] = 1
            
        elif j%2 == 0:
            V[i][j] = np.cos(j*Xt[i])
            
        else:
            V[i][j] = np.sin(j*Xt[i])
           
#print('   1        cos1x           sen1x         cos2x         sen2x'   )
#print('V =',V)
Vt = V.T
#print(Vt)

# Vetor de resultados b
R = np.zeros(2*ordem+1)
for i in range(2*ordem+1):
    R[i] = Vt[i].dot(Y)

# Solução do sistema linear
#print('R =', R)
Ptr = np.linalg.solve(D,R)


# Soma dos termos da função
T = 0
#print('Ptr =',Ptr)
for i in range(len(Ptr)):
        if i == 0:
            T += Ptr[i] 
        elif i%2 == 0:
            T += Ptr[i]*cos(i*x)
        else:
            T += Ptr[i]*sin(i*x)


# Impressão da função no terminal
print('\n\n\nPolinomio trigonometrico de ordem '+str(ordem)+' com '+str(pontos)+' pontos\n' )
print(T)      
            


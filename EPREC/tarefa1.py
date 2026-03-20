import numpy as np
from matplotlib import pyplot as plt
import math

# Funções da tarefa a:
# Integração de Crank

def CriarMatriz(N,M): # cria vetores a serem utilizados
    a = np.zeros(N-1) # Vetor Inicial (sem fronteiras)
    q = np.zeros(N-1) # Vetor b que recebe a com complementos
    F1 = np.zeros(M+1) # Fronteira 0
    F2 = np.zeros(M+1) # Fronteira 1
    return a, q, F1, F2


def LDLt_Decomp(type,lbd,N): # Cria matriz decomposta ou dada uma matriz, ela a decompoe.
    if type == "Tridiagonal": # Tridiagonal Simetrica
        a = 2*(1+lbd)  # N-1 Diagonal Maior
        b = -lbd  # N-2 Diagonal Menor
        L_lista = np.zeros(N-2) # l1 l2 l3 l4...
        D_lista = np.zeros(N-1) # d1 d2 d3 d4 d5...
        D_lista[0] = a
        for i in range(N-2):
            L_lista[i] = b/D_lista[i]
            D_lista[i+1] = a - b*L_lista[i]
        return D_lista, L_lista



def Find_f_1_and_2(p,deltaX,t,t2,x): # As forcantes para todos os testes
    h = deltaX
    rt = (10*(1+math.cos(5*t)))
    rt2 = (10*(1+math.cos(5*t2)))
    gh = (1/h)
    if p+h/2 >= x >= p-h/2:
        f = rt*gh
        f2 = rt2*gh
    else:
        f = 0
        f2 = 0

    return f, f2

def CrankMethod(lbd,a,f,f2,i,deltaT,N): # Crank classico, usa matriz "a"(sem aditivos) e coloca aditivos
    if i == 0:
        s = 2*(1-lbd)*a[i] + lbd*(0+a[i+1]) + (deltaT)*(f + f2) + (lbd)*0
    elif i == N-2:
        s = 2*(1-lbd)*a[i] + lbd*(a[i-1]+0) + (deltaT)*(f + f2) + (lbd)*0
    elif i not in (0, N-2):
        s = 2*(1-lbd)*a[i] + lbd*(a[i-1]+a[i+1]) + (deltaT)*(f + f2)
    return s

def solveLDLt(type,q,N,L_lista,D_lista,lbd):  # Resolve LDLt para...
    if type == "Triadiagonal": #... a tridiagonal simetrica
        Y = np.zeros(N-1)
        Y[0] = q[0]
        for i in range(1, N-1):
            Y[i] = q[i] - L_lista[i-1]*Y[i-1]
    
            X = np.zeros(N-1)
            X[-1] = Y[-1]/D_lista[-1]
        for i in range(2, N):
            X[-i] = (Y[-i] - (-lbd)*X[-(i-1)])/D_lista[-i]
        return X

    
def plotting(lista,comment,teste,N,episolon,tipo): # Plota a distribuicao de calor de alguma lista dada
    y=lista
    x=np.linspace(0,1,len(lista))
    plt.figure(figsize=(7,7))
    plt.plot(x,y)
    plt.suptitle(comment)
    if teste == "f":
        plt.title('U(T,x) do Teste {0}, com N={1}, com  ε={2}, GS {3}'.format(teste, N,episolon ,tipo))
    elif teste == "a":
        plt.title('U(T,x) do Teste {0}, com 4 fontes, {1}'.format(teste,tipo))
    else:
        plt.title('U(T,x) do Teste {0}, com N={1}, GS {2}'.format(teste, N, tipo))
    plt.xticks(np.arange(0, 1, 0.1))
    plt.grid()
    if teste in ("b","c","d"):
        plt.ylim(0,25)
    elif teste in ("a"):
        plt.ylim(-2,2)
    elif teste in ("e","f"):
        plt.ylim(0,140)
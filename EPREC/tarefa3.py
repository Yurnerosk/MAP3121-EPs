import numpy as np
import math


def Produto(A,B): # Produto interno entre A e B
    bar=0
    indice=len(A)
    for i in range(indice):
        bar=bar + A[i]*B[i]
    return bar

def Schmidt(store,tipo,qual,given): #Gram Schmidt

    K=len(store) # NUMERO DE FONTES
    P=len(store[0]) # TAMANHO DA BARRA

    u = np.empty([K,P]) #vetores ortogonalizados(q nao unitarios)
    
    e = np.empty([K,P]) #e = q unitario

    print("\n--- Gran Schmidt " + tipo +" em processo... --- \n")
        
    if qual == "c":  # Schmidt Classico
        i=0
        while i<K:

            if i==0:
                u[i] = store[i]

                g=0
                while g<P:
                    raiz = math.sqrt(Produto(u[i],u[i]))
                    e[i,g] = u[i,g] / raiz
                    g=g+1

            else:
                segura=0
                
                g=0
                while g<i:
                    segura = segura + Proj(u[g],store[i])   #Proj_u1(v1)
                    g=g+1
                    
                u[i] = store[i] - segura
                
                g=0
                while g<P:
                    raiz = math.sqrt(Produto(u[i],u[i]))
                    e[i,g] = u[i,g] / raiz
                    g=g+1
            i = i+1
    
    
    elif qual == "m": # Schmidt Modificado
        i=0
        while i<K:
    
            if i==0:
                u[i] = store[i]
                
                g=0
                while g<P:
                    raiz = math.sqrt(Produto(u[i],u[i]))
                    e[i,g] = u[i,g] / raiz
                    g=g+1
            
            else:
                g=0
                while g<i:

                    if g==0:
                        u[i] = store[i] - Proj(u[0],store[i])
                    else:
                        
                        d=0
                        while d<i:  
                            u[i] = u[i] -  Proj(u[d],u[i])
                            d=d+1
                    g=g+1
                
                g=0
                while g<P:
                    raiz = math.sqrt(Produto(u[i],u[i]))
                    e[i,g] = u[i,g] / raiz
                    g=g+1
            i = i+1
        
    print("\n--- Gram Schmidt Finalizado! --- \n")
    
    R = np.empty([K,K]) # R
        
    for i in range(K):
        for j in range(K):
            if j>i:
                break
            else:
                R[i][j] = Produto(e[j],store[i])

    print("QR fatorado com sucesso! ")
    return R, e

def Proj(u,a): #Proj a em u
    up = Produto(u,a)
    down = Produto(u,u)
    proj = (up/down)*np.array(u)
    
    return proj

def Resolve_QR(Q,R,given,store): #Pega Q, R, e retorna coeficientes prontos
    S=len(Q)
    e1 = np.empty(S) #e1 = Qt*b
    
    K=len(store) # Fontes
    x = np.empty(K) # Solucao
    
    for i in range(S):
        e1[i]=Produto(Q[i],given)
    
    i=len(R)-1
    while i>-1:

        if i == len(R)-1:
            x[i] = e1[i]/R[i,i]
        else:
            add = 0
            
            j=K-1
            while j>i:
                add = add + x[j]*R[j,i]
                j=j-1
            x[i] = (e1[i] - add)/R[i,i]
        i=i-1
    print("QR resolvido com sucesso!")
    return x # c o e f i c i e n t e s    p r o n t o s

def Veritas(Q): #Quanto maior, menor a ortogonalidade
    S=len(Q)
    highest = 0
    
    for i in range(S):
        for j in range(S):
            if i!=j:
                result = abs(Produto(Q[i], Q[j]))
                if result > highest:
                    highest = result
                    
    return highest

def Erro(store, nf, x, given, N, deltaX, teste, tipo, episolon):
    alpha = store[:]
    
    i=0
    while i<nf:
        alpha[i] = x[i]*store[i]
        i=i+1
    beta=np.zeros(N-1)
    
    i=0
    while i<nf:
        beta = beta + alpha[i]
        i=i+1
        
    comment = "Solução obtida com os coeficientes"
    
    import tarefa1 as T1
    
    T1.plotting(beta,comment,teste,N,episolon,tipo)
    gamma =  given - beta
    
    delta = 0
    
    i=0
    while i<(N-1):
        delta = delta + gamma[i]**2
        i=i+1
        
    E = math.sqrt(deltaX * delta)
    return E
import numpy as np
import math

teste = str(input("Digite teste (a/b/c/d/e/f) ... : ")) # Escolher teste

print("Teste " + teste + " será realizado ...")


qual = str(input("GRAM SCHMIDT Clássico (c) ou Modificado (m) :")) # Escolhe qual estilo de Gram Schmidt
if qual == "c":
    tipo = "clássico"
elif qual == "m":
    tipo = "modificado"
    
print("Gram Schmidt " + tipo + " será realizado ...")

if teste == "a":
    switch = "on" # manual choice
    if switch == "on":
        A = [[1,2,-2,2,0,1], # esta transverso aqui, mas usei do jeito certo no programa
             [2,3,0,11,2,-3],
             [-1,9,4,5,-1,1],
             [3,16,-1,3,4,7]
             ]
    elif switch == "off":
        A = [[1,0,1,0,-1,1], # opcional, para testar
             [1,1,1,1,1,1],
             [1,1,1,1,1,2],
             [1,1,1,1,1,3],
             ]
    b = np.array([1,0,1,0,-1,1])
    store = np.array(A)
    given = b
    
    import tarefa3 as T3
    
    R, Q = T3.Schmidt(store,tipo,qual,given) # Decompoe store e given em R e Q
    
    x = T3.Resolve_QR(Q,R,given,store) # Usa Q, R para achar x

    ax = np.zeros(len(store[0]))
    for i in range(len(store[0])): # coluna
        for j in range(len(store)):
            ax[i]= A[j][i] * x[j] + ax[i]
    
    DIST = given - ax
    dist = math.sqrt(T3.Produto(DIST,DIST))
    print("Distancia entre b e ax :")
    print(dist)
    
    print(" ")
    print("Verificando Ortogonalidade :")
    print("Valor Máximo = ",T3.Veritas(Q))    
    
    import tarefa1 as T1
    
    episolon = 0
    N=0
    comment = "Solução dada pelo enunciado"
    T1.plotting(given,comment,teste,N,episolon,tipo)
    
    comment = "Solução obtida pelos coeficientes"
    T1.plotting(ax,comment,teste,N,episolon,tipo)

if teste == "b": # predefinido, sem ruido
    N=128
    nf = 4
    position = [0.15, 0.3, 0.7, 0.8]
    
if teste in ("c","d","e","f"): # de arquivo TXT, sem ruido (c) ou com ruido (d)
    N=int(input("Digite N: (128/ 256/ 512/ 1024/ 2048) ... : "))

    
    if teste in ("c","d"):
        file1 = open("teste1.txt","r") 
        episolon = 0.01 # para o ruido
        nf = 10
    elif teste in ("e","f"):
        file1 = open("teste2.txt","r")
        episolon = 0
        if teste == "f":
            ask = int(input("Digite 1 para episolon = 0.01  ou Digite 2 para episolon = 0.00001  : "))
            if ask ==1:
                episolon = 0.01
            elif ask ==2:
                episolon = 0.00001
        print("Escolhido episolon = ", episolon)
        nf = 25
        
    f=file1.readlines()
    
    import tarefa2 as T2
    
    position ,given = T2.separar(N,f,teste,episolon,tipo)
    
if teste != "a":                                                                        ### integração
    (lbd, deltaX, deltaT, M) = (N, 1/N, 1/N, N)
    
    print("--- Teste {0}: N={1} com {2} fonte(s), nas posicoes {3}. --- \n".format(teste, N, nf, position))
    print("Calculado U(T,x)... 0%")
    for posicao in range(nf):
        p = position[posicao] # posicao da fonte
        
        import tarefa1 as T1
        
        D_lista, L_lista = T1.LDLt_Decomp("Tridiagonal",lbd,N) # matrizes ja foram decompostas
        a, q, F1, F2= T1.CriarMatriz(N,M) # Matrizes vazias para serem usadas
    
        for k in range(0, M): # O primeiro resolve para T = 0
            t = (k)*deltaT #k
            t2 = (k+1)*deltaT #k+1
            for i in range(0, N-1):  # 0 para N-2
                x = (i+1)*deltaX #i
                    
                f, f2 = T1.Find_f_1_and_2(p,deltaX,t,t2,x) #

                q[i] = T1.CrankMethod(lbd,a,f,f2,i,deltaT,N)
            a = T1.solveLDLt("Triadiagonal",q,N,L_lista,D_lista,lbd)
    
        if posicao==0: # Comeco a guardar os valores de Uk(T,x) na matriz "store"
            store=[a]
        else:
            store=np.append(store,[a],axis=0)
            
        print("Calculado U{1}(T,x)... {0}%".format(float(100*((posicao+1)/nf)),posicao)) ### integração
    print("\n--- Fontes Concluídas --- \n")
    if teste == "b":
        given=2.3*np.array(store[0])+3.7*np.array(store[1])+0.3*np.array(store[2])+4.2*np.array(store[3])
        comment = "Solução dada pelo enunciado"
        episolon = 0
        N=0
        T1.plotting(given,comment,teste,N,episolon,tipo)
    
    print("\n--- Iniciando Fatoração QR --- \n")
    
    import tarefa3 as T3
    
    R, Q= T3.Schmidt(store,tipo,qual,given)
    
    x = T3.Resolve_QR(Q,R,given,store) # Montado Rx=Qt*b
    
    print(" ")
    print("Seus coeficientes são:")
    for i in range(nf):
        print("a{0} = {1}".format(i,x[i]))
    
    
    if teste in ('c','d','e','f'):
        E=T3.Erro(store, nf, x, given, N, deltaX, teste, tipo, episolon)
        print("\nSeu Erro Quadrático:")
        print("E = {0}".format(E))
    
    print(" ")
    print("Verificando Ortogonalidade :")
    print("Valor Máximo = ",T3.Veritas(Q))    
    
    
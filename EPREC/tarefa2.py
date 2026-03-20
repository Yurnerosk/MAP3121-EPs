from random import random


def reader(lista): # Abre o arquivo TXT. Usado pelo separar()
    holder=""
    filtrado=[]
    retirar = " \n\t" # nem sei pra q serve esse t?
    for linha in lista:
        for letra in linha:
            if letra not in retirar:
                
                holder = holder + letra
                
            else:
                if holder != "":
                    filtrado.extend([holder])
                    holder = ""
                else:
                    pass
    if holder != "":
        filtrado.extend([holder])
    else:
        pass
    return filtrado

def separar(N,f,teste,episolon,tipo): # Escolhe qual parte do TXT deve pegar, dado N .
    step = int(2048/N)
    count = 0
    distr=[] # Parte N-1 da distribuicao
    pos=[] # Parte das posicoes
    while count <= 2048:
        if count == 0 or count == 2048:
            pass
        else:
            distr.extend([reader(f[count + 1])])
        count = count + step
    count=0
    for item in distr: # converter para float e talvez colocar um ruido (d)
        for conteudo in item:
            if teste in ("c","e"): # sem ruido
                distr[count]=float(conteudo)
            if teste in ("d","f"): # com ruido
                r=(random()-0.5)*2
                distr[count]= ( 1 + r*episolon ) * float(conteudo)
        count=count+1
    pos.extend(reader(f[0]))
    i=0
    for item in pos:
        pos[i]=float(item)
        i=i+1
    comment = "Solução do arquivo"
    
    import tarefa1 as T1
    
    T1.plotting(distr,comment,teste,N,episolon,tipo)
    return pos, distr

from Bio import SeqIO

def refinarsequencia(tamanhoMin,limiar, registros):
    #matriz em que as sequencias refinadas serao salvas
    tudo=[]

    for registro in registros:

       #Os nucliotideos que  estiverem acimar do limiar serao salvos aqui
        parte=[]
        media=0
        i=0
       #verifica da nucleotidio
        for nucleotideo in registro.letter_annotations['phred_quality']:
           # print("padrao")
           # print(registro[i])

            if(nucleotideo>=limiar):
                #print("refinada")
                #print(registro[i])
                #print(i)
                #localizacao.append(i)
                media= media+nucleotideo
                parte.append(registro[i])

            i = i + 1
       #  print("registro %i" % count)
       #  print(parte)
       # print(localizacao)
       #Adciona a sequencia a matriz caso a media da qualidade esteja acima do limiar
       # e a sequencia seja grande o suficiente
        if(parte.len()>tamanhoMin and  media/i>=limiar):
            tudo.append(parte)
        #print("refinou o %i registro" % count)


    return tudo




# transforma a lista em sequencia para melhor manipulacao
registros = list (SeqIO.parse("C:\\Users\\ricar\Downloads\\1-200613_S7_L001_R1_001.fastq\\1-200613_S7_L001_R1_001.fastq","fastq"))
count=0
#print(len(registros))
#for registro in registros:
#   print("regristro numero %i \n" %count)
#    print("%s \n %s \n" %(registro.id, registro.letter_annotations['phred_quality']))
#    count+=1
#    if count== 10:
#        break
refinados = refinarsequencia(60,20,registros)


for registro in registros:

   print("\n regristro numero %i \n" %count)
   print("%s \n  %s \n " %(registro.id, registro.seq))
   print("qualidade do registro \n")
   print(registro.letter_annotations['phred_quality'])
   print("%s registro refinado" )
   print(refinados[count])

   count+=1
   if count== 2:
        break

#print("numero total é %s" % count)

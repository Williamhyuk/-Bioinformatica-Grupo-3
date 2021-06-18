from Bio import SeqIO
import numpy as np
from itertools import product
import os, glob

folder_path = 'bd/'
fasta_paths = glob.glob(os.path.join(folder_path, '*.fasta'))
basededatos=[]
for fasta_path in fasta_paths:
    #print(fasta_path)
    for seq_record in SeqIO.parse(fasta_path, "fasta"):
        #print(seq_record.id)
        basededatos.append(str(seq_record.seq.upper()))
        #print(str(seq_record.seq.upper()))
        #print()

#print(basededatos)
sequences = SeqIO.parse("pruebas/Q4G0N8.fasta", "fasta")
for record in sequences:
   query = str(record.seq.upper())

print(query)

AMINOACID_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
BLOSUM62 = {
    ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
    ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
    ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
    ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
    ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
    ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
    ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
    ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
    ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
    ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
    ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
    ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
    ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
    ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
    ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
    ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
    ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
    ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
    ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
    ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
    ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
    ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
    ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
    ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
    ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
    ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
    ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
    ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
    ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
    ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
    ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
    ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
    ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
    ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
    ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
    ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
    ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
    ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
    ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
    ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
    ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
    ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
    ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
    ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
    ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
    ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
    ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
    ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
    ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
    ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
    ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
    ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
    ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
    ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
    ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
    ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
    ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
    ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
    ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
    ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
    ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
    ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
    ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
    ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
    ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
    ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
    ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
    ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
    ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4
}

def comparar(a,b):
   if (a, b) in BLOSUM62.keys():
        return BLOSUM62[(a, b)]
   else:
        return BLOSUM62[(b, a)]

def score(a,b):
  return comparar(a[0],b[0])+comparar(a[1],b[1])+comparar(a[2],b[2])

def BLAST(cad1,DB):
  #cad1="EYYPMVPGTYIVTITWGGQNIGRSPFEVKVGTECGNQKVRAWGPGLEGGVVGKSADFVVE"
  #DB=["PQLPITNFSRDWQSGRALGALVDSCAEYYPMVPDSWDASKPVTNAREAMQQADDWLGIPQ","VITPEEIVDPNVDEHSVMTYLSQFPKAKLKPGAPLRPKLNPKKARAYGPGIEPTGNMVKK","RAEFTVETRSAGQGEVLVYVEDPAGHQEEAKVTANNDKNRTFSVWYVPEVTGTHKVTVLF","AGQHIAKSPFEVYVDKSQGDASKVTAQGPGLEPSGNIANKTTYFEIFTAGAGTGEVEVVI","QDPMGQKGTVEPQLEARGDSTYRCSYQPTMEGVHTVHVTFAGVPIPRSPYTVTVGQACNP","SACRAVGRGLQPKGVRVKETADFKVYTKGAGSGELKVTVKGPKGEERVKQKDLGDGVYGF"]
  lencad1=len(cad1)
  partes=[[cad1[i:i+3]] for i in range(0,lencad1-2)]

  # vECINOS

  posibilidades=[''.join(j) for j in  product('CSTPAGNDEQHRKMILVFYW', repeat=3)]
  print("PARTES",partes)
  for i in partes:
    for j in posibilidades:
      if ((comparar(i[0][0], j[0])+comparar(i[0][1], j[1])+comparar(i[0][2], j[2]))>12):
        scor=comparar(i[0][0], j[0])+comparar(i[0][1], j[1])+comparar(i[0][2], j[2])
        #print('score=',scor)
        i.append(j)
    i.remove(i[0])
  print("PARTES",partes)

  #BLAST

  for k in range(len(partes)-1): # Recorre el vector con los vecinos [[abc,aeb,acc..],[],[]]
    for kk in partes[k]: # Recorre el subvector [abc,aeb,acc..]
      for i in DB: # Recorre ["","",""]
        for j in range(len(i)-2): #Recorre "ABCDEFGHI" DESDE "A" HASTA "G" 
          if kk[0]==i[j]:
            if kk[1]==i[j+1]:
              if kk[2]==i[j+2]:
                print("----------------")
                print("KK",kk)
                print("CAD1",cad1[k:k+3])
                print("i",i[j:j+3])
                ttt=score(cad1[k:k+3],i[j:j+3])
                print("score",ttt)
                izqQuery=k
                derQuery=k+2
                izqDB=j
                derDB=j+2
                print("iq",izqQuery)
                print("dq",derQuery)
                print("id",izqDB)
                print("dd",derDB)
                while (ttt>15):#22 PARA PROTINAS Y 20 PARA DNA
                  beneficioDer=0
                  beneficioIzq=0

                  if izqQuery>0 and izqDB>0:
                    beneficioIzq=comparar(cad1[izqQuery-1],i[izqDB-1])
                  else:
                    beneficioIzq=-1000
                  if derQuery<len(cad1)-1 and derDB<len(i)-1:
                    #print("d",derecha)
                    #print("kk",cad1)
                    #print("i",i)
                    beneficioDer=comparar(cad1[derQuery+1],i[derDB+1])
                  else:
                    beneficioDer=-1000
                  
                  if (beneficioDer>beneficioIzq):
                    derQuery+=1
                    derDB+=1
                    ttt+=beneficioDer
                  elif (beneficioDer<beneficioIzq):
                    izqQuery-=1
                    izqDB-=1
                    ttt+=beneficioIzq
                  else:
                    break
                    
                  #print("d1",derecha1)
                  #print("i1",izquierda1)
                  #print("i2",izquierda2)
                  #print("d2",derecha2)
                  

                #print("i1",izquierda1)
                #print("d1",derecha1)
                
                #print("i2",izquierda2)
                #print("d2",derecha2)

                print("cad1",cad1[izqQuery:derQuery+1])
                print("cadi",i[izqDB:derDB+1])
                print("scoreFinal",ttt)


BLAST(query,basededatos)

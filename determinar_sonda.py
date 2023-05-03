from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO

# Lê as sequências a partir de um arquivo fasta
records = SeqIO.parse("/home/jpedro/ProjetosBioinfo/alvos_para_qPCR_multiplex/arquivo_fasta.fasta", "fasta")

# Salva as sequências em um arquivo temporário
SeqIO.write(records, "temp.fasta", "fasta")

# Define a linha de comando do ClustalW
clustalw_cline = ClustalwCommandline("clustalw2", infile="temp.fasta")

# Executa o clustalw2
stdout, stderr = clustalw_cline()

# Lê o resultado do alinhamento
alinhamento = AlignIO.read("temp.aln", "clustal")

from fuzzywuzzy import fuzz
#import argparse
#import PySimpleGui

seq1 = alinhamento[0].seq
seq2 = alinhamento[1].seq
regioes_divergentes = []
similaridade_sonda = float(input('similaridade da sonda: '))

for i in range(0, len(seq1), 30):
    sub_lista1 = "".join(seq1[i:i+30])
    sub_lista2 = "".join(seq2[i:i+30])
    similaridade_seq = fuzz.ratio(sub_lista1, sub_lista2) #retorna a similaridade entre as duas listas pulando de 2 em 2
    if similaridade_seq < similaridade_sonda:
        regioes_divergentes.append([sub_lista2, similaridade_seq])

nucleotideos_conservados = []

for var1, var2 in zip(seq1, seq2):

    posicao_var1 = seq1.index(var1)
    posicao_var2 = seq2.index(var2)

    if var1 == var2 and ((seq1[posicao_var1:posicao_var1+20]) == (seq2[posicao_var2:posicao_var2+20])):
        nucleotideos_conservados.append(var1)

regioes_conservadas = []

for i in range(0, len(seq2), 20):
    if str(seq2[i:i+20]) in regioes_divergentes:
        continue    
    regioes_conservadas.append(str(seq2[i:i+20])) # Depois tentar i:i+100, por exemplo, para desenhar os primers...
    # Ou continuar encontrando as regiões de 20 em 20 para depois expandir elas.

# print(regioes_conservadas) 

print('----------------------------------------\n----------------------------------------\n----------------------------------------\n----------------------------------------\n')

print('Sondas: ')
for j in regioes_divergentes:
    print(f' |A sequência é [{j[0]}]|\n |Sua similaridade é {j[1]}|\n--------------------')

# Adicionar as coordenadas junto ao retorno da sequência e da similaridade (posição do nucleotídeo inicial:posição do nucleotídeo final)
# Criar região Conservada-Divergente-Conservada 
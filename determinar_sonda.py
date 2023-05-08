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
#import argparse || estudar esta biblioteca
#import PySimpleGui || estudar esta biblioteca

seq1 = alinhamento[0].seq
seq2 = alinhamento[1].seq
regioes_divergentes = []
similaridade_sonda = float(input('similaridade da sonda: '))

regioes_conservadas = []
for i in range(0, len(seq1), 30):
    sub_lista1 = "".join(seq1[i:i+30])
    sub_lista2 = "".join(seq2[i:i+30])
    similaridade_seq = fuzz.ratio(sub_lista1, sub_lista2) #retorna a similaridade entre as duas listas pulando de 2 em 2
    if similaridade_seq == 100:
        regioes_conservadas.append(sub_lista2)
    elif similaridade_seq <= similaridade_sonda:
        regioes_divergentes.append([sub_lista2, similaridade_seq])

# Depois tentar i:i+100, por exemplo, para desenhar os primers...
# Ou continuar encontrando as regiões de 20 em 20 para depois expandir elas.

print(regioes_conservadas) 

print('\n----------------------------------------\n----------------------------------------\n----------------------------------------\n----------------------------------------\n')

print('Sondas: ')
for j in regioes_divergentes:
    print(f' |A sequência é [{j[0]}]|\n |Sua similaridade é {j[1]}|\n--------------------')

# Adicionar as coordenadas junto ao retorno da sequência e da similaridade (posição do nucleotídeo inicial:posição do nucleotídeo final)
# Criar região Conservada-Divergente-Conservada 
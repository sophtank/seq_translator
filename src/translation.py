#!/bin/usr/env python3

import argparse, csv, os
import pandas as pd
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
required = parser.add_argument_group("Required Arguments")
required.add_argument('-i','--input', type=str, help = 'path to intended input file')
required.add_argument('-o', '--output', type=str, help='path to intended output TEXT file')
parser.add_argument('-p', '--pandas', type=str, help = "path to intended output GRAPHIC file")
args = parser.parse_args()
file_in = open(args.input)
lines = file_in.readlines()

header = []
dna = []
rna = []

#creating two separate lists - a list of headers and a list of sequences because input is .fasta and not string literal
for line in lines:
    line = line.strip()

    if line[0] == '>':
        header.append(line)
    if line[0] != '>':
        dna.append(line)

#going protein translation
rna = [r.replace('T', 'U') for r in dna]
print(rna)

#dictionary of all amino acid codons
codons = {
        "UUU" : "F", "UUC" : "F", "UUA" : "L", "UUG" : "L",
        "CUU" : "L", "CUC" : "L", "CUA" : "L", "CUG" : "L",
        "AUU" : "I", "AUC" : "I", "AUA" : "I", "AUG" : "M",
        "GUU" : "V", "GUC" : "V", "GUA" : "V", "GUG" : "V",
        "UCU" : "S", "UCC" : "S", "UCA" : "S", "UCG" : "S",
        "CCU" : "P", "CCC" : "P", "CCA" : "P", "CCG" : "P",
        "ACU" : "T", "ACC" : "T", "ACA" : "T", "ACG" : "T",
        "GCU" : "A", "GCC" : "A", "GCA" : "A", "GCG" : "A",
        "UAU" : "Y", "UAC" : "Y", "UAA" : "STOP", "UAG" : "STOP",
        "CAU" : "H", "CAC" : "H", "CAA" : "Q", "CAG" : "Q",
        "AAU" : "N", "AAC" : "N", "AAA" : "K", "AAG" : "K",
        "GAU" : "D", "GAC" : "D", "GAA" : "E", "GAG" : "E",
        "UGU" : "C", "UGC" : "C", "UGA" : "STOP", "UGG" : "W",
        "CGU" : "R", "CGC" : "R", "CGA" : "R", "CGG" : "R",
        "AGU" : "S", "AGC" : "S", "AGA" : "R", "AGG" : "R",
        "GGU" : "G", "GGC" : "G", "GGA" : "G", "GGG" : "G",
    }


protein = []

def diction(rna):
    #finding the longest unit in the list, finding the length of that to dictate range in the second for loop
    longest = int(max(enumerate(rna), key=lambda x: len(x[1]))[0])
    x = len(rna[longest])
    for i in range(0, len(rna)):
        current_string = rna[i]
        for i in range (0,x,3):
            temp = current_string[i:i+3]
            new = codons[temp]
            if new != "STOP":
                protein.append(new)
            if new == "STOP":
                protein.append(new)
                break
    totalprotein = ''.join(protein)
    global listform 
    listform = totalprotein.split("STOP")

diction(rna)


with open(args.output, 'w') as reader:
        for i in range(0, len(header)):
            reader.write(header[i] + "\n")
            reader.write(listform[i] + "\n")

for i in range(0,len(dna)):
    #AT and GC counts
    AT = dna[i].count("A") + dna[i].count("T")
    GC = dna[i].count("G") + dna[i].count("C")
    total = AT + GC
    ATC = (AT/total) * 100
    GTC = (GC/total) * 100
    GTC = round(GTC, 2)
    ATC = round(ATC, 2)

print(ATC)
print(GTC)
row = ["Content" , "Count"]
row1 = ["AT Content" , ATC]
row2 = ["GC Content", GTC]

os.remove('/home/stank/seq_translator/ATGC.csv')
with open('ATGC.csv', 'a') as csv_file:
   writer = csv.writer(csv_file)
   writer.writerow(row)
   writer.writerow(row1)
   writer.writerow(row2)


df = pd.read_csv('/home/stank/seq_translator/ATGC.csv')

content = df["Content"]
amount = df["Count"]
colors = ["#D88C9A" , "#F2D0A9"]
explode = (0.1 , 0)
plt.pie(amount , labels=content , explode=explode, colors=colors,autopct='%1.1f%%',shadow=True, startangle=140)
plt.title("AT and GC Content by Percentage")
plt.savefig(args.pandas)
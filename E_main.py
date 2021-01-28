# -*- coding: utf-8 -*-

import E_bruteforce as brute
import E_greedy as greed
import E_randomized as rand

'''DNA = ['tagtggtcttttgagtgtagatctgaagggaaagtatttccaccagttcggggtcacccagcagggcagggtgacttaat',
'cgcgactcggcgctcacagttatcgcacgtttagaccaaaacggagttggatccgaaactggagtttaatcggagtcctt',
'gttacttgtgagcctggttagacccgaaatataattgttggctgcatagcggagctgacatacgagtaggggaaatgcgt',
'aacatcaggctttgattaaacaatttaagcacgtaaatccgaattgacctgatgacaatacggaacatgccggctccggg',
'accaccggataggctgcttattaggtccaaaaggtagtatcgtaataatggctcagccatgtcaatgtgcggcattccac',
'tagattcgaatcgatcgtgtttctccctctgtgggttaacgaggggtccgaccttgctcgcatgtgccgaacttgtaccc',
'gaaatggttcggtgcgatatcaggccgttctcttaacttggcggtgcagatccgaacgtctctggaggggtcgtgcgcta',
'atgtatactagacattctaacgctcgcttattggcggagaccatttgctccactacaagaggctactgtgtagatccgta',
'ttcttacacccttctttagatccaaacctgttggcgccatcttcttttcgagtccttgtacctccatttgctctgatgac',
'ctacctatgtaaaacaacatctactaacgtagtccggtctttcctgatctgccctaacctacaggtcgatccgaaattcg']'''

DNA = ['tagcatcgc',
    'cctagaacg',
    'ggtgcatag']

k=3

print('\nBruteForceMotifSearch: ')
test_motif = brute.BruteForceMotifSearch(DNA, k)
print('>tablica pozycji: ', test_motif)
print('>motywy w kolejnych sekwencjach: ')
for i, index in enumerate(test_motif):
    print(DNA[i][index : index+3])


print('\nGreedyMotifSearch: ')
test_motif, test_score = greed.GreedyMotifSearch(DNA, k)   
print('>tablica pozycji: ', test_motif)
print('>score: ', test_score)
print('>motywy w kolejnych sekwencjach: ')
for i, index in enumerate(test_motif):
    print(DNA[i][index : index+k])

print('\nRandomizedMotifSearch: ')       
score, test_motif = rand.RandomizedMotifSearch(DNA, k)
print('>motif: ', test_motif, '\n>score: ', score)
        


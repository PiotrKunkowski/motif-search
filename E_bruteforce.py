# -*- coding: utf-8 -*-

import E_score as sc
import itertools as it

def BruteForceMotifSearch(DNA, k):
    '''przeszukiwanie metodą brutalnej siły. iteruje po wszystkich możliwosciach.
    zwraca macierz pozycji maksymalizujących score'''
    assert k <= len(DNA[0]), "k-mer nie może być dłuższy niż łańcuch, w którym go szukamy"
    positions = [0 for x in range(len(DNA))]
    best_score = 0
    comb_list=[x for x in range(len(DNA[0])-k+1)]
    combined = list(it.product(comb_list, repeat = len(DNA))) 
    #funkcja product generuje wszystkie możliwe kombinacje podanych list
    #repeat mówi ile razy mamy powtórzyć operację
    for a in combined:
        if (sc.score(a, DNA, k) > best_score):
            best_score = sc.score(a, DNA, k)
            positions = a
    assert(positions != []), "macierz pozycji nie może być pusta"
    return positions

if __name__ == '__main__':
    print("____________TEST WYDAJNOŚCI FUNKCJI BRUTE FORCE____________")
    test_DNA = ['tagccatgc',
    'cctagaacg',
    'ggtgcatag']
    
    test_motif = BruteForceMotifSearch(test_DNA, 3)
    print('>tablica pozycji: ', test_motif)
    print('>motywy w kolejnych sekwencjach: ')
    for i, index in enumerate(test_motif):
        print(test_DNA[i][index : index+3])
    
    assert (test_motif == (0,2,6)), "funkcja zle wyszukuje motywy"
    assert (len(test_motif) == 3), "motyw musi mieć okreslona dlugosc"
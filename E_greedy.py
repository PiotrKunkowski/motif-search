# -*- coding: utf-8 -*-

import E_score as sc

def GreedyMotifSearch(DNA, k):
    motifs = [0 for x in range(len(DNA))]
    best_motifs = [0 for x in range(len(DNA))]
    best_score = sc.score(best_motifs, DNA, k)
    for j in range(0, len(DNA[0]) - k):
        #iteruje po k-merach w pierwszym łańcuchu
        motifs[0] = j 
        for i in range(1, len(DNA)):
            #szuka najbardziej prawdopodobnego k-meru w każdej sekwencji po kolei (zaczynając od drugiej)
            my_profile = sc.profile(motifs, DNA, k)
            mot, motifs[i] = sc.fit_motif_to_profile(DNA[i], my_profile, k)
            
            if sc.score(motifs, DNA, k) > best_score:
                best_motifs = []
                best_score = sc.score(motifs, DNA, k)
                for motif in motifs:
                    #best motifs nadpisywał się
                    best_motifs.append(motif)
    return best_motifs, best_score


if __name__ == '__main__':
    print("____________TEST WYDAJNOŚCI FUNKCJI GREEDY____________")
    test_DNA = ['tagcatcgc',
    'cctagaacg',
    'ggtgcatag']
    
    test_motif, test_score = GreedyMotifSearch(test_DNA, 3)
    
    print('>tablica pozycji: ', test_motif)
    print('score: ', test_score)
    print('>motywy w kolejnych sekwencjach: ')
    for i, index in enumerate(test_motif):
        print(test_DNA[i][index : index+3])
    
    assert (test_motif == [3,0,4]), "funkcja zle wyszukuje motywy"
    assert (len(test_motif) == 3), "motyw musi mieć okreslona dlugosc"
    
    
    
    
    
    
    

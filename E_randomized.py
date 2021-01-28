# -*- coding: utf-8 -*-

import random
import E_score as sc

def random_starting_positions(DNA, k):
    '''zwraca losowe pozycje startowe i k-mery na tych miejscach'''
    pos_tab = []
    for seq in DNA:
        start_point = random.randrange(len(seq)-k)
        pos_tab.append(start_point)
    print('losowe pozycje poczatkowe: ', pos_tab)
    assert(len(pos_tab) > 0), "tablica pozycji początkowych nie może być pusta"
    return pos_tab


def RandomizedMotifSearch(DNA, k):
    '''Randomized Motif Search opiera się na losowaniu macierzy pozycji początkowych i budowaniu 
    na jej podstawie macierzy profilu. Znajduje motyw najlepiej pasujący do danej macierzy profilu
    i porównuje score wylosowanego motywu z tym zbudowanym na podstawie profilu.
    '''
    my_motif = [] #motyw znaleziony na podstawie macierzy profilu, uzupełniany niżej
    positions = random_starting_positions(DNA, k)
    max_score = sc.score(positions, DNA, k)
    best_motif = positions
    index_list = []
    assert(positions != []), "tablica pozycji początkowych nie może być pusta"
    assert(best_motif == positions), "wylosowany motyw jest póki co najlepszym motywem"
    while True:
        my_profile = sc.profile(positions, DNA, k) #buduje profil
        for sekwencja in DNA:
            #iteruje po łańcuchach DNA i dodaje do listy indeksy
            fit_motif, index = sc.fit_motif_to_profile(sekwencja, my_profile,k)
            my_motif.append(fit_motif)
            index_list.append(index)
        my_score = sc.score(index_list, DNA, k)
        if my_score > max_score:
            return my_score, my_motif
        else:
            best_motif=[]
            for i, index in enumerate(positions):
                best_motif.append(DNA[i][index : index+k])
            return max_score, best_motif


if __name__ == '__main__':
    print("____________TEST WYDAJNOŚCI FUNKCJI RANDOMIZED____________")
    test_DNA = ['tagccatgc',
    'cctagaacg',
    'ggtgcatag']
    
    for i in range(0, 50):
        test_score, test_motif = RandomizedMotifSearch(test_DNA, 3)
        print(test_score, test_motif)

    assert (len(test_motif) == 3), "motyw musi mieć okreslona dlugosc"

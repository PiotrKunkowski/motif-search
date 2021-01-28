# -*- coding: utf-8 -*-  

def profile(motif, DNA, k):
    '''wraca profil prawdopodobienstwa
    motif - motyw o dlugosci k, DNA - macierz sekwencji nukleotydowych, k - dlugosc meru'''
    profile_matrix = []
    for x in range(0, k):  
        #iteruje po długoci meru
        nuc = [0, 0, 0, 0]     #lista częstoci występowania poszczególnych nukleotydów
        for y in range(0, len(motif)):    
            #układa profil iterując po łańcuchach
            if DNA[y][x + motif[y]] == 'a':   #jeśli w sekwencji nr y na miejscu x + indeks nukleotydem jest adenina
                nuc[0] += 1                   #wartosc A w liscie zwieksza sie o 1
            if DNA[y][x + motif[y]] == 'c':
                nuc[1] += 1
            if DNA[y][x + motif[y]] == 'g':
                nuc[2] += 1
            if DNA[y][x + motif[y]] == 't':
                nuc[3] += 1
        for i in range(0, len(nuc)):            
            #iteruje po długości listy
            nuc[i] = (nuc[i] + 1) / len(motif)  
        profile_matrix.append(nuc)
    return profile_matrix

def fit_motif_to_profile(seq, profile, k):
    '''na podstawie profilu znajduje mer, który najbardziej pasuje do tego profilu.
    macierz profilu opiera się na funkcji prawdopodobieństwa'''
    max_probab = 0  
    for i in range(0, len(seq) - k):  
        # iteruje po sekwencjach
        prob = 1  #wartosc prawdopodbienstwa, póki co niska
        motif = seq[i:i + k]  # Przypisanie aktualnie sprawdzanego meru
        for j, nuc in enumerate(motif):  
            # iteruje po merach
            #funkcja enumerate służy nadaniu indeksów
            if nuc == 'a':
                prob = prob * profile[j][0]
            elif nuc == 'c':
                prob = prob * profile[j][1]
            elif nuc == 'g':
                prob = prob * profile[j][2]
            else:  
                #sprawdzanie znaku oraz zwiekszanie prawdopodobienstwa o odpowiednia wartosc
                prob = prob * profile[j][3]
        if prob > max_probab:  
            max_probab = prob
            best_motif = motif  
            index = i
    return best_motif, index  # Zwracany jest indeks o najwyzszym prawdopodobienstwie i motyw


def score(s, DNA, k):
    """ 
    s - tablica pozycji poczatkowych, DNA - tablica stringów (sekwencji), k - dlugosc meru
    funkcja oblicza score dla k-meru 
    """
    my_score = 0
    for i in range(k):
        # iteruje po pozycjach w sekwencji
        nuc_dictionary =	{   #słownik nukleotydów, tu będzie przechowywana macierz profilu
                "a": 0,
                "c": 0,
                "t": 0,
                "g": 0
                }
        for j, sval in enumerate(s):
            # iteruje po łańcuchach DNA
            base = DNA[j][sval+i] 
            nuc_dictionary[base] += 1
        my_score += max(nuc_dictionary.values())
    return (my_score)

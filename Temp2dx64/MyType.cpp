#ifndef MY_TYPE
#define MY_TYPE

#include <vector>

// Struktura reprezentująca węzeł
struct Node {
    double x, y, t, CR; // współrzędne i inne parametry węzła
    int status;         // status węzła
};

// Struktura reprezentująca element skończony
struct Element {
    int nop[4];        // Indeksy węzłów elementu
    int Npov;          // Liczba kontaktowych powierzchni elementu (0-1-2)
    int aPov[2];       // Lokalne numery kontaktowych powierzchni elementu
};

// Struktura reprezentująca siatkę MES w 2D
struct Gr2d {
    int nh;            // Liczba węzłów
    int ne;            // Liczba elementów
    int nbn;           // Liczba węzłów na brzegu
    int ncn;           // Liczba węzłów kontaktowych
    int nhPov;         // Liczba kontaktowych powierzchni siatki

    std::vector<Element> EL; // Wskaźnik na elementy siatki
    std::vector<Node> ND;    // Wskaźnik na węzły siatki
};

// Struktura dla współrzędnych lokalnych
struct CorL {
    double N, E; // Współrzędne lokalne
};

// Struktura dla granic siatki MES
struct Gran {
    int N_p;                     // Liczba punktów całkowania po powierzchni granicy
    std::vector<CorL> P;         // Punkty całkowania (współrzędne lokalne)
    std::vector<double> W;       // Wagi całkowania
    std::vector<std::vector<double>> N1, N2, Nf; // Funkcje kształtu węzłów
    std::vector<char> UZEL;      // Lokalne numery węzłów granicy
};

// Struktura reprezentująca element (ELEM)
struct Elem {
    int nbn;                     // Liczba węzłów na brzegu
    int nbnp;                    // Liczba punktów brzegowych
    int N_p;                     // Liczba punktów całkowania
    std::vector<std::vector<double>> N1, N2, Nf; // Funkcje kształtu
    std::vector<CorL> P;         // Punkty całkowania
    std::vector<CorL> L;         // Współrzędne węzłów elementu
    std::vector<double> W;       // Wagi punktów całkowania
    Gran Sf[4];                  // Granice elementu (max 4)
};

#endif
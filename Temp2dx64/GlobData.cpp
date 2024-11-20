#ifndef  GLOB_DATA
#define  GLOB_DATA

#include <vector>
#include "MyType.cpp"

class GlobData {
public:
    // Globalne zmienne
    double mTbegin;           // Początkowa temperatura
    double mTime;             // Czas procesu
    double mdTime;            // Początkowa wartość przyrostu czasu
    double mTau;              // Bieżący czas
    double mT_otoczenia;      // Temperatura otoczenia
    double mAlfa;             // Współczynnik wymiany ciepła
    double mH0;               // Wysokość przekroju
    double mB0;               // Szerokość przekroju
    int mNhH;                 // Liczba węzłów na wysokości
    int mNhB;                 // Liczba węzłów na szerokości
    int mLDA;                 // Szerokość pasma macierzy MES
    double mC;                // Pojemność cieplna
    double mK;                // Współczynnik przewodzenia ciepła
    double mR;                // Gęstość

    Gr2d mGr;                 // Siatka MES
    Elem mEL4;             // Dane użytego elementu skończonego

    // Kontrolne punkty
    int mContrPoints[9];
    double mcpX[9], mcpY[9];

    // Macierze i wektory
    double est[4][4];         // Macierz bieżącego elementu
    double r[4];              // Wektor obciążeń bieżącego elementu
    std::vector<double> mB;   // Globalny wektor obciążeń
    std::vector<std::vector<double>> mA; // Globalna macierz
    std::vector<double> mX;   // Globalny wektor niewiadomych
};

#endif
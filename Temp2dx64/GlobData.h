#ifndef GLOB_DATA_H
#define GLOB_DATA_H

#include <vector>
#include <array>
#include <Eigen/Sparse>
#include "MyType.h"

struct GlobData {

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
    std::array<double, 9> mcpX; // Współrzędne X punktów kontrolnych
    std::array<double, 9> mcpY; // Współrzędne Y punktów kontrolnych
    std::array<int, 9> mContrPoints; // Indeksy punktów kontrolnych w węzłach

    // Macierze i wektory
    double est[4][4];         // Macierz bieżącego elementu
    double r[4];              // Wektor obciążeń bieżącego elementu
    // Macierze i wektory
    Eigen::SparseMatrix<double> mA; // Globalna macierz równań (rzadka)
    Eigen::VectorXd mB;            // Globalny wektor obciążeń
    Eigen::VectorXd mX;            // Globalny wektor niewiadomych
};

#endif

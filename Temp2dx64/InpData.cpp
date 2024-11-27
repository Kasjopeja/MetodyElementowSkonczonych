#include <iostream>
#include <fstream>
#include <vector>
#include "InpData.h"

extern GlobData data;

void InpData() {
    std::ifstream inFile("C:\\Users\\Lenovo\\CLionProjects\\Metody Elementow Skonczonych\\Temp2dx64\\dane.txt");

    if (!inFile.is_open()) {
        std::cerr << "Blad: Nie udalo sie otworzyc pliku 'dane.txt'" << std::endl;
        return;
    }

    std::string temp; // Zmienna tymczasowa do odczytu opisów

    // Odczyt danych wraz z pominięciem opisów
    inFile >> data.mTbegin; getline(inFile, temp); // Odczyt wartości i pominięcie opisu
    inFile >> data.mTime; getline(inFile, temp);
    inFile >> data.mdTime; getline(inFile, temp);
    inFile >> data.mT_otoczenia; getline(inFile, temp);
    inFile >> data.mAlfa; getline(inFile, temp);
    inFile >> data.mH0; getline(inFile, temp);
    inFile >> data.mB0; getline(inFile, temp);
    inFile >> data.mNhH; getline(inFile, temp);
    inFile >> data.mNhB; getline(inFile, temp);
    inFile >> data.mC; getline(inFile, temp);
    inFile >> data.mK; getline(inFile, temp);
    inFile >> data.mR; getline(inFile, temp);

    inFile.close();
    std::cout << "Dane wczytane poprawnie z pliku 'indata.t2d'" << std::endl;
}

void ALLOCATE_Matrix() {
    int mLDA = 0;
    int NeMaxB = 0;

    // Iterowanie przez elementy w siatce
    for (size_t NEL = 0; NEL < data.mGr.EL.size(); ++NEL) {
        const auto &element = data.mGr.EL[NEL];
        std::vector<int> nk(data.mEL4.nbn);

        // Przypisywanie numerów węzłów
        for (int i = 0; i < data.mEL4.nbn; ++i) {
            nk[i] = element.nop[i];
        }

        // Obliczanie szerokości pasma macierzy
        for (int i = 0; i < data.mEL4.nbn; ++i) {
            int ii = nk[i];
            for (int j = 0; j < data.mEL4.nbn; ++j) {
                int jj = nk[j];
                int jB = jj - ii + 1;
                if (jB >= mLDA) {
                    mLDA = jB;
                    NeMaxB = NEL;
                }
            }
        }
    }

    data.mLDA = mLDA;

    try {
        // Alokacja macierzy i wektorów
        data.mA.resize(data.mLDA, std::vector<double>(data.mGr.nh, 0.0));
        data.mB.resize(data.mGr.nh, 0.0);
        data.mX.resize(data.mGr.nh, 0.0);
    } catch (const std::bad_alloc &e) {
        std::cerr << "Błąd alokacji pamięci: " << e.what() << std::endl;
    }
}
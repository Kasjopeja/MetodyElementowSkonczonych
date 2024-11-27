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

    // Zainicjalizowanie macierzy i wektorów
    data.mLDA = 0; // Początkowo zerowa, ustalana później
    data.mB.resize(data.mNhH * data.mNhB, 0.0); // Wektor globalny
    data.mA.resize(data.mNhH * data.mNhB, std::vector<double>(data.mNhH * data.mNhB, 0.0)); // Macierz globalna
    data.mX.resize(data.mNhH * data.mNhB, 0.0); // Wektor niewiadomych

    inFile.close();
    std::cout << "Dane wczytane poprawnie z pliku 'indata.t2d'" << std::endl;
}

#ifndef AUXILIAR_HPP
#define AUXILIAR_HPP

// Bibliotecas padrão C++
#include <iostream>				// criação de Input/output files
#include <fstream>				// criação de classes de Input/output files
#include <vector>				// cria e manipula vetores
#include <cmath>				// antigo <math.h> - libera funções matemáticas, ex. pow
#include <string>				// manipulação de strings
#include <sstream>
#include <chrono>				// manipulação do tempo



std::vector < double > LeComp(std::string);
std::vector < std::vector < double > > LeMatriz(std::string);
void EscreveMatriz(std::vector < std::vector < double > >, std::string, std::string);
void MostraMatrz(std::vector < std::vector < double > >);
void MostraVetor(std::vector < double >);
void MostraVetorInt(std::vector < int >);
void MostraTempo(double);
void EscreveLog(std::string, std::string, double, double, double, double, double, unsigned int = 0);

#endif // AUXILIAR_HPP
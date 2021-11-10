#ifndef MATRIZ_HPP
#define MATRIZ_HPP

// Bibliotecas padrão C++
#include <iostream>				// criação de Input/output files
#include <fstream>				// criação de classes de Input/output files
#include <vector>				// cria e manipula vetores
#include <cstdlib>				// gera numeros randomicos, habilita o system("PAUSE")
#include <cmath>				// antigo <math.h> - libera funções matemáticas, ex. pow
#include <functional>			// transfere funções entre objetos
#include <limits>				// Limites das variáveis, usado em Random
#include <algorithm>			// função find em EnxameParticulas.cpp
#include <ctime>				// antigo <time.h> - Biblioteca de manipulações do tempo
#include <iomanip>				// manipulação de Input/output files
#include <locale>				// utilitários de localização
#include <string>				// manipulação de strings
#include <omp.h>		//*\\	// parelelismo, usado em Enxame de Partícula
#include <random>				// manipulação de números randômicos
#include <chrono>				// manipulação do tempo


template <class T>
class Matriz
{
    public:
        Matriz(int , int);
        std::vector < std::vector < T > > M;

};

class MatrizOperacoes
{
    public:
        std::vector < std::vector <double> > R;
        std::vector < std::vector <double> > T;
        void MultMatriz(std::vector < std::vector <double> >&  , std::vector < std::vector <double> >&  );
        void MatrizTransposta(std::vector < std::vector <double> >& );
};



#endif // MATRIZ_HPP

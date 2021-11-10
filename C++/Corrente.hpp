#ifndef CORRENTE_HPP
#define CORRENTE_HPP

// Bibliotecas internas
#include "PropsMX.hpp"


class Corrente : public PropsMX
{
    public:
        
        // Variáveis
        std::string nome;
        std::vector < bool > descartadas;
            
        // funções
        Corrente(std::string, double, double, double, bool = 1, std::vector < double > = {1});
        void RemoveNulas(std::vector < double >);
        void BalancoFases(bool = 1);
        void MarcaComoSolvente();
        
};
#endif // CORRENTE_HPP

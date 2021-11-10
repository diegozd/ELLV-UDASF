#ifndef MISTURADOR_HPP
#define MISTURADOR_HPP

#include "Corrente.hpp"

class Misturador 
{
    public:

        // funções
        Misturador();
        Corrente mistura(Corrente&, Corrente&, std::string, double = 0, double = 0);
        std::vector <int> BuscaRepetidas(Corrente&, Corrente&);
        Corrente BalancoMassa(Corrente&, Corrente&);


};

#endif // MISTURADOR_HPP
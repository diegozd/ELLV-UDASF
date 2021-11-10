#ifndef SPLITFASE_HPP
#define SPLITFASE_HPP

#include "Corrente.hpp"

class SplitFase 
{
    public:

        // funções
        SplitFase();
        std::vector < Corrente > SepFases(Corrente&);
        Corrente SepSolvente(Corrente&, std::string);
};

#endif // SPLITFASE_HPP
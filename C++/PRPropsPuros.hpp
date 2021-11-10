#ifndef PRPROPSPUROS_HPP
#define PRPROPSPUROS_HPP

#include <math.h> 

// Bibliotecas internas
// #include "Matriz.hpp"
#include "MGProps.hpp"

class PRPropsPuros : public MGProps
{
    public:
        // Estruturas
        typedef struct 
        {
            double A;
            double B;
            double Tr;
            double a;
            double a_alpha;
            double b;
            double alpha;
            double kappa;
            double Tc;
            double Pc;
            double w;
            double q;
        } CtsPRStr;

        // Variáveis
        std::vector < double > TbPR;
        std::vector < double > TbPR_viaG;
        unsigned int maxiter = 100;
        
        // Funções
        PRPropsPuros(std::string);
        void CalcTbPR(double);
        void CalcTbPR_viaG(double);
        double FZPR(double, double, double);
        double FTPR(double, double, unsigned int);
        double FTPR_viaG(double, double, unsigned int);
        double HRPR(double, double, unsigned int, bool);
        double SRPR(double, double, unsigned int, bool);
        double GRPR(double, double, unsigned int, bool);
        CtsPRStr calcCtesPR(double, double,  unsigned int);
        double CalcRho(double, double, unsigned int);
        double FZLPR(double, double, double);
        double FZVPR(double, double, double);
};

#endif // PRPROPSPUROS_HPP

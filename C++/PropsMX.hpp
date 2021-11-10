#ifndef PROPSMX_HPP
#define PROPSMX_HPP

// Bibliotecas internas
// #include "Matriz.hpp"
#include "PRPropsPuros.hpp"
#include "auxiliar.hpp"
#include <math.h> 
#include <Eigen/Dense>
// #include <stdlib.h>     /* abs */

class PropsMX : public PRPropsPuros
{
    public:

        // estruturas
        typedef struct 
        {
            double tpd;                     //tpd da fase teste
            bool EstFis;                    //estado fisico da fase
            double MMw;                     //Massa molecular média da fase
            double roh;                     //Massa especifica da fase por PR @ T e P da corrente [kg/m³]
            double Z;                       //Fator de compressibildiade da fase por PR @ T e P da corrente
            double vazM;                    //Vazão mássica da fase [kg/s]   
            double vazN;                    //Vazão molar da fase [kmol/s]
            double beta;                    //fração molar da fase na corrente    
            std::vector < double > y;       //composição da fase 
            std::vector < double > phi;     //coeficiente de fugacidade de cada componente na fase     
            std::vector < double > vazNi;   //vazão molar de cada componente na fase 
            std::vector < double > vazMi;   //vazão massica de cada componente na fase  
        } FaseStr;
        
        typedef struct
        {
            double T;                                   //Temperatura da mistura [K]
            double P;                                   //Presão da mistura [bar]
            double vazM;                                //Vazão mássica da corente [kg/s] 
            double vazN;                                //Vazão molar da corente [kmol/s]
            double MMw;                                 //Massa molar média da mistura kg/kmol
            double roh;                                 //Massa especifica média das fases da corrente calculada por PR @ T e P da corrente [kg/m³]
            std::vector < double > z;                   //Composição global da mistura
            std::vector < double > vazNi;               //vazão molar de cada componente global da mistura [kg/s]
            std::vector < double > vazMi;               //vazão massica de cada componente  global da mistura [kmol/s]
            std::vector < FaseStr > F;                  //vetor de fases do flash
            std::vector < std::vector < double > > K;   // matriz de K=yi/ykeyF onde tem os componentes nas linhas e as fases nas colunas
        } FlshStr;

        typedef struct 
        {
            std::vector < double > E;
            double Q;
        } Qstr;

        typedef struct
        {
            std::vector < std::vector < double > > aij;
            std::vector < CtsPRStr > Cprs;      
        } PurosCtsStr;

        typedef struct
        {
            PurosCtsStr PCts;
            double A;
            double B;
            double a;
            double b; 
            double q;   
        } MxCtsStr;


        // Variáveis
        FlshStr FLSH;                //objeto que contem o resultado do equilíbrio de fases da corrente 
        double Taij;
        unsigned int warning = 0;
        bool Jamostrei = 1;

        // Funções
        PropsMX(std::string);
        PurosCtsStr CalcCtesPuros(double, double, bool = 0);
        MxCtsStr CalcCtesMx(double, double, std::vector < double >, bool = 0);
        double ZPRTP(double, double, bool, std::vector < double >);
        std::vector < double > phiMxi(double, double, bool, std::vector < double >, bool = 0);
        FaseStr tpd(bool, std::vector < double >, bool, FlshStr);
        std::vector < FaseStr > AnaliseEstabilidade(FlshStr, bool = 0); 
        Qstr CalcQ(FlshStr);
        void FlashELLV(double = 0, double = 0, std::vector < double > = {1});

        //Funções teste não utilizadas
        double ResRachfordRice(double, std::vector < double >, std::vector < double >);
        FlshStr OrvalhoTELV(double , std::vector < double > );
        FlshStr BolhaTELV(double , std::vector < double > );
        FlshStr FlashLV(double , double , std::vector < double > );
        void Envelope(double);

        
        
        
        
};

#endif // PROPSMX_HPP
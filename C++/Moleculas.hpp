#ifndef MOLECULAS_HPP
#define MOLECULAS_HPP

// Bibliotecas padrão C++
#include <iostream>				// criação de Input/output files
#include <fstream>				// criação de classes de Input/output files
#include <vector>				// cria e manipula vetores
#include <cmath>				// antigo <math.h> - libera funções matemáticas, ex. pow
#include <string>				// manipulação de strings
#include <sstream>
#include <stdlib.h>     /* abs */

// Bibliotecas internas
#include "Matriz.hpp"

class Moleculas
{

	private:
		unsigned int tamanhoLigas;
		unsigned int tamanhoAtribMSOL;

		std::vector < std::vector < std::vector <long long int> > > MSOL;
		std::vector < std::vector < std::vector <int> > > Ligas;
		std::vector < std::vector <double> > MSOL_expandida;

		MatrizOperacoes T_Mjoback_SOL;
		MatrizOperacoes Mjoback_moleculas;
		MatrizOperacoes T_Mest;
		MatrizOperacoes Matomos_moleculas;
		MatrizOperacoes T_Matomos_moleculas;
		MatrizOperacoes Mmw;

    public:
		
		// estruturas
        typedef struct 
        {
			bool SOLV;			//flag que indica se esse componente é solvente
            double MMw;         //Massa molecular de uma molecula kg/kmol
            double Tm;          //Normal melting point (Tm [K])
            double Tb;          //Normal boiling point (Tb [K])
            double Tc;          //Critical temperature (Tc [K])
            double Pc;          //Critical pressure (Pc [bar])
            double Vc;          //Critical Volum (Vc [cm³/mol])
            double Gf;          //Standard Gibbs energy at 298 K (Gf [kJ/mol])
            double Hf;          //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
            double Hv;          //Standard enthalpy of vaporization at 298 K (Hv [kJ/mol])
            double Hfus;        //Standard enthalpy of fusion (Hfus [kJ/mol])
            double VmolRckt;    //Vol molar (kmol/m3) na Tref=298K volume molar rackett, pag 407 perry.
            double w;           //Fator acênctrico
            double rho;         //Massa especifica da molécula calculada por PR @ T e P da corrente [kg/m³]
            // double vazM;        //Vazão mássica da molécula global na corente [kg/s] 
            // double vazN;        //Vazão molar da molécula global na corente [kmol/s]
            std::vector < std::vector <int> > Ligas;
            // std::vector < std::vector <long long int> > MSOL;
            std::vector <double> MSOL_expandida;
			std::vector <double> MGgroupsPO;
			std::vector <double> Matomos_moleculas;
        } MoleStr;

		//variaveis
		unsigned int nMol;
		unsigned int nGroupsPO;
		std::vector <MoleStr> mol;

		//funcoes
		Moleculas(std::string);
		void LeLigas(std::string);
		void LeMSOL_expandida(std::string);
		void LeMSOL(std::string);
};

#endif // MOLECULAS_HPP

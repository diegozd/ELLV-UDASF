// Bibliotecas padrão C++
#include <iostream>				// criação de Input/output files
#include <fstream>				// criação de classes de Input/output files
#include <vector>				// cria e manipula vetores
#include <cmath>				// antigo <math.h> - libera funções matemáticas, ex. pow
#include <string>				// manipulação de strings
#include <chrono>				// manipulação do tempo
#include <sstream>
#include <filesystem>
#include <unistd.h>
#include <stdlib.h>     /* abs */
#include <filesystem>
#include <cstring>


// Bibliotecas internas
#include "Matriz.hpp"			// classe de Matrizes
#include "Moleculas.hpp"		// classe da geração molecular
#include "MGProps.hpp"
#include "PRPropsPuros.hpp"
#include "PropsMX.hpp"
#include "Corrente.hpp"
#include "auxiliar.hpp"
#include "Misturador.hpp"
#include "SplitFase.hpp"

// Definindo variáveis globais
#define _USE_MATH_DEFINES


int main(int argc, char** argv)
{
	//inicia contagem do tempo
	auto start = std::chrono::high_resolution_clock::now();
	
	auto tmp = std::filesystem::current_path();
	std::string ProjPath = tmp;
	ProjPath += "/run/";

	// mensagens de boas vindas
	std::cout << std::endl;
	std::cout << "ELLV para moleculas de Reconstrucao molecular" << std::endl;
	std::cout << "Diego Telles Fernandes" << std::endl;
	std::cout << "PEQ - COPPE - UFRJ" << std::endl;
	std::cout << "Trabalho da disciplina Termodinamica de Solucoes (COQ712) " << std::endl;
	std::cout << "Compilado em " << __DATE__ << " as " << __TIME__ << std::endl << std::endl;

	// Dia da execução do algoritimo
	std::chrono::system_clock::time_point today = std::chrono::system_clock::now();
	std::time_t tt;
	tt = std::chrono::system_clock::to_time_t ( today );
	std::cout << ctime(&tt);
	std::cout << std::endl;

	// Definindo valores padrão para variáveis
	// double Tatm = 25 +273.15;	// K
	double TC = 60;				// C
	double T = TC + 273.15;		// K
	double Patm = 1.01325;		// bar
	double Pbarkg = 1.01972;	// (kgf/cm2g)/bar
	double P = Patm;			// bar
	double Pkg = 0 ;			// kgf/cm2g
	double RSOM = 4;			// RSO massica desejada
	std::string nomeRV = "RV";
	std::string nomeSOLV = "C3";

	// Lendo flags de entrada
	for (int i = 1 ; i < argc ; i++)
	{
		// leitura da Temperatura (T) se fornecida pelo usuário
		if(strcmp(argv[i], "-T") == 0) T = std::stod(argv[i+1]);

		// leitura da Temperatura (T) se fornecida pelo usuário
		if(strcmp(argv[i], "-TC") == 0) T = std::stod(argv[i+1]) + 273.15;

		// leitura da Pressão (P) se fornecida pelo usuário
		if(strcmp(argv[i], "-P") == 0) Pkg = std::stod(argv[i+1]);

		// leitura da RSO se fornecida pelo usuário
		if(strcmp(argv[i], "-RSO") == 0) RSOM = std::stod(argv[i+1]);

		// leitura da Corrente de RV se fornecida pelo usuário
		if(strcmp(argv[i], "-RV") == 0) nomeRV = argv[i+1];

		// leitura da Corrente de RV se fornecida pelo usuário
		if(strcmp(argv[i], "-SOLV") == 0) nomeSOLV = argv[i+1];
	}

	//ajustando pressão
	P = (Pkg+1.03323)/Pbarkg;

	// Definindo Razão Solvente Óleo molar
	double fSOLV = 15;				// fator de conversão RSO massica em molar para Solvente Generico
	// double fC7 = 7.5;				// fator de conversão RSO massica em molar para Solvente C7
	double RSOSOLV = fSOLV*RSOM;	// RSO molar Solvente Generico
	// double RSOC7M = 1;				// RSO massica desejada para C7
	// double RSOC7 = fC7*RSOC7M;		// RSO molar C3

	//Criando corrente RV
	std::string pathRV = ProjPath + nomeRV + "/";
	double vazNRV = 1;		//kmol/s
	Corrente RV(pathRV, T, P, vazNRV);
	std::cout << " Objeto " << RV.nome << " construido com sucesso" << std::endl << std::endl;


	// comentário geral inicia aqui
	/**/

	// Criando corrente Solvente Generico
	std::string pathSOLV = ProjPath + nomeSOLV + "/";
	double vazNSOLV = vazNRV*RSOSOLV;	//kmol/s
	Corrente SOLV(pathSOLV, T, P, vazNSOLV, 0);
	fSOLV = RV.FLSH.MMw/SOLV.FLSH.MMw;
	RSOSOLV = fSOLV*RSOM;
	SOLV.FLSH.vazN = vazNRV*RSOSOLV;
	SOLV.BalancoFases();
	SOLV.MarcaComoSolvente();
	std::cout << " Objeto " << SOLV.nome << " construido com sucesso" << std::endl;
	std::cout << " Real RSO SOLV:" << SOLV.FLSH.vazM/RV.FLSH.vazM  << std::endl;
	std::cout << std::endl;

	// //Criando corrente C7
	// std::string nomeC7 = "C7";
	// std::string pathC7 = ProjPath + nomeC7 + "/";
	// double vazNC7 = vazNRV*RSOC7;	//kmol/s
	// Corrente C7(pathC7, Tatm, Patm, vazNC7, 0);
	// fC7 = RV.FLSH.MMw/C7.FLSH.MMw;
	// RSOC7 = fC7*RSOC7M;
	// C7.FLSH.vazN = vazNRV*RSOC7;
	// C7.BalancoFases();
	// C7.MarcaComoSolvente();
	// std::cout << " Objeto " << C7.nome << " construido com sucesso" << std::endl;
	// std::cout << " Real RSO C7:" << C7.FLSH.vazM/RV.FLSH.vazM  << std::endl;
	// std::cout << std::endl;

	//Criando o objeto misturador
	Misturador MX;

	// //Criando corrente de mistura C7+RV
	// Corrente RVC7_MX = MX.mistura(C7,RV,"RVC7_MX",Tatm);
	// std::cout << " Objeto " << RVC7_MX.nome << " construido com sucesso" << std::endl;
	// std::cout << std::endl;

	//Criando corrente de mistura C3+RV
	Corrente RVSOLV_MX = MX.mistura(SOLV,RV,"RVSOLV_MX");
	std::cout << " Objeto " << RVSOLV_MX.nome << " construido com sucesso" << std::endl;
	std::cout << std::endl;
	// if (warning > 0) msg.append(RVSOLV_MX.nome);

	//Cria objeto separador de fases
	SplitFase SEP;

	// criando vetor de correntes para mistura RVSolvente
	std::vector < Corrente > VCRVSOLV = SEP.SepFases(RVSOLV_MX);
	std::cout << " Objeto com vetor de correntes a partir de " << RVSOLV_MX.nome << " foi construido com sucesso" << std::endl;
	
	
	//Avaliando rendimento UDASF
	unsigned int iOdes = 0;
	double minMMw = 1e10;
	for (unsigned int j = 0 ; j < VCRVSOLV.size() ; j++)
	{
		if (VCRVSOLV[j].FLSH.MMw < minMMw) 
		{
			iOdes = j;
			minMMw = VCRVSOLV[j].FLSH.MMw;
		}
	}
	
	Corrente ODES = SEP.SepSolvente(VCRVSOLV[iOdes], "ODES");
	Corrente RASF = ODES;
	unsigned int iRASF1 = 0, iRASF2 = 0, iRASF3 = 0;
	std::cout << "VCRVSOLV.size():" << VCRVSOLV.size() << std::endl;
	if (VCRVSOLV.size() == 1) std::cout << "NAO FORMOU DUAS FASES COM SOLV" << std::endl;
	else if (VCRVSOLV.size() == 2) RASF = SEP.SepSolvente(VCRVSOLV[1-iOdes], "RASF");
	else if (VCRVSOLV.size() == 3)
	{
		if (iOdes == 0) { iRASF1 = 1; iRASF2 = 2;}
		else if (iOdes == 1) { iRASF1 = 0; iRASF2 = 2;}
		else if (iOdes == 2) { iRASF1 = 0; iRASF2 = 1;}
		Corrente RASF_1 = SEP.SepSolvente(VCRVSOLV[iRASF1], "RASF_1");
		Corrente RASF_2 = SEP.SepSolvente(VCRVSOLV[iRASF2], "RASF_2");
		RASF = MX.mistura(RASF_1, RASF_2, "RASF");
	}
	else if (VCRVSOLV.size() == 4)
	{
		if (iOdes == 0) { iRASF1 = 1; iRASF2 = 2; iRASF3 = 3;}
		else if (iOdes == 1) { iRASF1 = 0; iRASF2 = 2; iRASF3 = 3;}
		else if (iOdes == 2) { iRASF1 = 0; iRASF2 = 1; iRASF3 = 3;}
		else { iRASF1 = 0; iRASF2 = 1; iRASF3 = 2;}
		Corrente RASF_1 = SEP.SepSolvente(VCRVSOLV[iRASF1], "RASF_1");
		Corrente RASF_2 = SEP.SepSolvente(VCRVSOLV[iRASF2], "RASF_2");
		Corrente RASF_3 = SEP.SepSolvente(VCRVSOLV[iRASF3], "RASF_3");
		Corrente RASF_m = MX.mistura(RASF_1, RASF_2, "RASF_m");
		RASF = MX.mistura(RASF_m, RASF_3, "RASF");
	}
	else std::cout << "RV+SOLVENTE FORMOU MAIS DE 4 FASES" << std::endl;

	// //avaliando teor de asfaltenos
	// unsigned int iSolvC7 = 0;
	// minMMw = 1e10;
	// for (unsigned int j = 0 ; j < RVC7_MX.FLSH.F.size() ; j++)
	// {
	// 	if (RVC7_MX.FLSH.F[j].MMw < minMMw) 
	// 	{
	// 		iSolvC7 = j;
	// 		minMMw = RVC7_MX.FLSH.F[j].MMw;
	// 	}
	// }
	// std::cout << std::endl;
	// std::cout << std::endl;
	// std::cout << "####    TEOR DE ASFALTENOS    ####" << std::endl;
	// if (RVC7_MX.FLSH.F.size() == 1) std::cout << "NAO FORMOU DUAS FASES COM C7" << std::endl;
	// else if (RVC7_MX.FLSH.F.size() == 2)
	// {
	// 	std::cout << "SOLUVEIS EM C7:   " << 100*(RVC7_MX.FLSH.F[iSolvC7].vazM - RVC7_MX.FLSH.F[iSolvC7].vazMi[0])/RV.FLSH.vazM << " %" << std::endl;
	// 	std::cout << "INSOLUVEIS EM C7: " << 100*(RVC7_MX.FLSH.F[1-iSolvC7].vazM - RVC7_MX.FLSH.F[1-iSolvC7].vazMi[0])/RV.FLSH.vazM << " %" << std::endl;
	// }
	// else if (RVC7_MX.FLSH.F.size() == 3)
	// {
	// 	unsigned int Insol1 = 0, Insol2 = 0;
	// 	if (iSolvC7 == 0) { Insol1 = 1; Insol2 = 2;}
	// 	else if (iSolvC7 == 1) { Insol1 = 0; Insol2 = 2;}
	// 	else if (iSolvC7 == 2) { Insol1 = 0; Insol2 = 1;}
	// 	std::cout << "SOLUVEIS EM C7:   " << 100*(RVC7_MX.FLSH.F[iSolvC7].vazM - RVC7_MX.FLSH.F[iSolvC7].vazMi[0])/RV.FLSH.vazM << " %" << std::endl;
	// 	std::cout << "INSOLUVEIS EM C7: " << 100*((RVC7_MX.FLSH.F[Insol1].vazM - RVC7_MX.FLSH.F[Insol1].vazMi[0])+(RVC7_MX.FLSH.F[Insol2].vazM - RVC7_MX.FLSH.F[Insol2].vazMi[0]))/RV.FLSH.vazM << " %" << std::endl;
	// }
	// else std::cout << "RV+C7 FORMOU MAIS DE 3 FASES" << std::endl;
	// std::cout << std::endl;

	std::cout << std::endl;
	std::cout << std::setprecision(2) << std::fixed;
	std::cout << "###############################################" << std::endl;
	std::cout << "####       RENDIMENTOS DESASFALTACAO       ####" << std::endl;
	std::cout << "###############################################" << std::endl;
	std::cout << std::endl;

	if (RVSOLV_MX.warning > 0)
	{
		std::cout << "####    ERRO:    ####" << std::endl;
		std::cout << "Problema na formação de duas fases." << std::endl  << std::endl;
	}

	std::cout << "T: " << T << " K  =>  " << T-273.15 << " C" << std::endl;
    std::cout << "P: " << P << " bar =>  " << P * Pbarkg - 1.03323 << " kgf/cm2g" << std::endl;
	std::cout << std::endl;

	std::cout << "OLEO:      " << RV.nome << std::endl;
	std::cout << "SOLVENTE:  " << SOLV.nome << std::endl;
	std::cout << "RSO:       " << SOLV.FLSH.vazM/RV.FLSH.vazM  << std::endl;
	std::cout << "rend ODES: " << 100*(ODES.FLSH.vazM)/RV.FLSH.vazM << " %" << std::endl;
	std::cout << "rend RASF: " << 100*(RASF.FLSH.vazM)/RV.FLSH.vazM << " %" << std::endl;
	std::cout << std::endl;

	EscreveLog(RV.nome, SOLV.nome, T, P , SOLV.FLSH.vazM/RV.FLSH.vazM, 100*(ODES.FLSH.vazM)/RV.FLSH.vazM, 100*(RASF.FLSH.vazM)/RV.FLSH.vazM, RVSOLV_MX.warning);
	
	//Final do comentário geral
	/**/

	// Escrevendo propriedades da mistura
	// RV.CalcTbPR(Patm);
	// std::cout << "mol" << "\t\t" << "MMW" << "\t\t" << "TbPR[K]" << "\t\t" << "TbMG[K]" << "\t\t" << "Tc[K]" << "\t\t" << "Pc[bar]" << "\t\t" << "w" << std::endl;
	// for (unsigned int i = 0 ; i < RV.mol.size() ; i++)
	// {
	// 	std::cout << std::setprecision(1) << std::fixed;
	// 	std::cout << i  << "\t\t" << RV.mol[i].MMw << "\t\t" << RV.TbPR[i] << "\t\t" << RV.mol[i].Tb << "\t\t" << RV.mol[i].Tc << "\t\t" ;
	// 	std::cout << std::setprecision(2) << std::fixed;
	// 	std::cout <<  RV.mol[i].Pc << "\t\t" ;
	// 	std::cout << std::setprecision(4) << std::fixed;
	// 	std::cout << RV.mol[i].w  << std::endl;
	// }

	// RV.Envelope(Patm);
	



	// // salva matriz kij
	// RVSOLV_MX.GroupKij(T);
	// EscreveMatriz(RVSOLV_MX.kij, pathRV, "kij_SOLV.csv");
	// std::cout << std::setprecision(4) << std::fixed;
	// for (unsigned int i = 0 ; i < RV.FLSH.z.size() ; i++) std::cout << RV.mol[i].MMw << "\t" << RV.FLSH.z[i] << std::endl;




	//finaliza contagem do tempo e mostra
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = duration_cast <std::chrono::microseconds> (stop - start);
	double TempoDouble = (double) duration.count() ;
	MostraTempo(TempoDouble);

	std::cout << std::endl << "FIM!" << std::endl << std::endl;
    return 0;
}
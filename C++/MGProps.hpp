#ifndef MGPROPS_HPP
#define MGPROPS_HPP

#include <math.h> 

// Bibliotecas internas
// #include "Matriz.hpp"
#include "Moleculas.hpp"

class MGProps : public Moleculas
{
    private:
        // Variáveis
        int nProps;
        int MaxElementoGeral;
        unsigned int maxLigas_e;
        bool aviso;
        
        std::vector < int > MaxElemento;
        std::vector < double > ParamAjustMG;
        std::vector < std::vector < double > > _MGgroupsPO;
        std::vector < std::vector < int > > nElementos;
        std::vector < std::vector < double > > ParMGPO;
        std::vector < std::vector < double > > PropsMGPuros;
        std::vector < int > SUB_UNIFAC;
        std::vector < double > FA;
	    std::vector < std::vector < double > > sub_FA;
        std::vector < std::vector < double > > MG_UNIFAC;
        std::vector < std::vector < double > > MG_FA;
        std::vector < std::vector < double > > parametros_SW;
        std::vector < std::vector < double > > Grupos_alfa_Mol;
        std::vector < std::vector < double > > Akl;
        std::vector < std::vector < double > > Bkl;
        std::vector < std::vector < std::vector < double > > > MolecPETROSIM;

        MatrizOperacoes somaContriMG;
        MatrizOperacoes Molec_UNIFAC;
	    MatrizOperacoes Molec_FA;
        MatrizOperacoes SumNiwi;
    
    public:
        //variáveis
        double R;
        double RPa;
        double Tref;
        std::vector < std::vector < double > > kij;

        // Funções
        MGProps(std::string, bool = 1);
        bool CH3(int, int);
        bool CH2(int, int);
        bool CH(int, int);
        bool aCH(int, int);
        bool aCfaC(int, int);
        bool aCfnC(int, int);
        bool nCH2(int, int);
        bool nCH(int, int);
        bool aC_CH3(int, int);
        bool aC_CH2(int, int);
        bool aC_CH(int, int);
        bool aC_C(int, int);
        bool aN(int, int);
        bool CH2SH(int, int);
        bool CH2NH2(int, int);
        bool aCqq (int, int);
        bool aCS (int, int);
        bool aCN (int, int);
        bool aCO (int, int);
        bool SH (int, int);
        bool NH2 (int, int);
        bool OH (int, int);
        bool CH2CO (int, int);
        bool CH3CO (int, int);
        bool CHCO (int, int);  
        bool aC_CO (int, int);    
        bool aldeido_CO (int, int);
        bool aldeido_aC_CO (int, int);
    

        bool pmul_e(int, int, int);
        bool tal_e(int, int, int);
        bool toel_e(int, int);
        bool CH2pertenceKO(int, int);
        bool CHpertenceKO(int, int);
        void LendoTipos();
        void BuscaUltimoC();
        int ContaLigaAnel(int, int);
        int AnelPred(int, int);
        bool aCpertence_aCCH2aC(int, int);
        std::vector< std::vector < int > > GCH2COCH2(int);
        std::vector< std::vector < int > > GCH3CO(int);
        std::vector< std::vector < int > > GCHCO(int);
        std::vector< std::vector < int > > GaC_CO(int);
        std::vector< std::vector < int > > GaCCH2aC(int);
        std::vector< std::vector < int > > Galdeido_aC_CO(int);
        std::vector<int> BuscaLigas_e(int, int);
        void ContaMG();
        void CalcPropsPuros();
        void GroupUNIFAC();
	    void GroupKij(double);
	    void Fator_Ac();
        void IniciaGrupos_alfa_Mol(unsigned int);
};



#endif // MGPROPS_HPP

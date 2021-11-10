#include "SplitFase.hpp"

SplitFase::SplitFase() { }


std::vector < Corrente > SplitFase::SepFases(Corrente &IN)
{
    std::vector < Corrente > VCs;
    Corrente AuxC = IN; 
    AuxC.nome = "";
    AuxC.FLSH.vazN = 0;
    AuxC.FLSH.vazM = 0;
    AuxC.FLSH.z.clear();
    AuxC.FLSH.vazNi.clear();
    AuxC.FLSH.vazMi.clear();
    AuxC.FLSH.F.clear();
    AuxC.FLSH.K.clear();

    for (unsigned int j = 0 ; j < IN.FLSH.F.size() ; j++)
    {   
        AuxC.nome = IN.nome + "_F" + std::to_string(j);
        AuxC.FLSH.vazN = IN.FLSH.F[j].vazN;
        AuxC.FLSH.vazM = IN.FLSH.F[j].vazM;
        AuxC.FLSH.MMw = IN.FLSH.F[j].MMw;
        AuxC.FLSH.roh = IN.FLSH.F[j].roh;
        AuxC.FLSH.z = IN.FLSH.F[j].y;
        AuxC.FLSH.vazNi = IN.FLSH.F[j].vazNi;
        AuxC.FLSH.vazMi = IN.FLSH.F[j].vazMi;
        AuxC.FLSH.F.push_back(IN.FLSH.F[j]);

        AuxC.FLSH.K.resize(AuxC.FLSH.z.size());
        for (unsigned int i = 0 ; i < AuxC.FLSH.z.size() ; i++) AuxC.FLSH.K[i].resize(1,1);

        VCs.push_back(AuxC);
    }
    return VCs;
}

Corrente SplitFase::SepSolvente(Corrente &IN, std::string nome)
{
    Corrente AuxC = IN;
    std::vector < double > SumY(IN.FLSH.F.size(), 0);
    AuxC.nome = nome;
    AuxC.FLSH.vazN = 0;
    AuxC.FLSH.vazM = 0;
    AuxC.nMol = 0;
    AuxC.mol.clear();
    AuxC.FLSH.z.clear();
    AuxC.FLSH.vazNi.clear();
    AuxC.FLSH.vazMi.clear();
    AuxC.FLSH.F.clear();
    AuxC.FLSH.K.clear();

    std::cout << std::endl << std::endl << "####    CORRENTE " << AuxC.nome << "    ####" << std::endl;

    for (unsigned int i = 0 ; i < IN.nMol ; i++)
    {
        if (IN.mol[i].SOLV == 0)
        {//molecula não é solvente
            AuxC.mol.push_back(IN.mol[i]);
            AuxC.nMol++;

            AuxC.FLSH.vazNi.push_back(IN.FLSH.vazNi[i]);
            AuxC.FLSH.vazMi.push_back(IN.FLSH.vazMi[i]);
            AuxC.FLSH.vazN += IN.FLSH.vazNi[i];
            AuxC.FLSH.vazM += IN.FLSH.vazMi[i];
        }

    }
    AuxC.FLSH.MMw = AuxC.FLSH.vazM/AuxC.FLSH.vazN;

    for (unsigned int i = 0 ; i < AuxC.nMol ; i++) AuxC.FLSH.z.push_back(AuxC.FLSH.vazNi[i]/AuxC.FLSH.vazN);

    AuxC.FlashELLV();
    AuxC.BalancoFases();

    return AuxC;
}
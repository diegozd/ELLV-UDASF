#include "Misturador.hpp"

Misturador::Misturador() { }

Corrente Misturador::mistura(Corrente &C1, Corrente &C2, std::string nome, double T, double P)
{   
    Corrente OUT = BalancoMassa(C1, C2);
    OUT.nome = nome;
    std::cout << std::endl << std::endl << "####    CORRENTE " << OUT.nome << "    ####" << std::endl;
    if (T == 0) T = (C1.FLSH.T * C1.FLSH.vazM + C2.FLSH.T * C2.FLSH.vazM)/OUT.FLSH.vazM;
    if (P == 0)
    {
        if (C1.FLSH.P < C2.FLSH.P) P = C1.FLSH.P;
        else P = C2.FLSH.P;
    }
    // OUT.Jamostrei = 0;   //mostrar kij
    OUT.FlashELLV(T,P);
    OUT.BalancoFases();

    return OUT;
}

std::vector <int> Misturador::BuscaRepetidas(Corrente &C1, Corrente &C2)
{
    
    unsigned int marcador = 0;
    unsigned int nMolC2 = C2.mol.size();

    unsigned int nMolDescartadas = 0; //numero de moleculas repetidas
    std::vector <int>  descartadas(nMolC2, -1);

    for (unsigned int i = 0 ; i < C1.mol.size() ; i++)
    {//varrendo moléculas de C1
        
        for (unsigned int j = 0 ; j < C2.mol.size() ; j++)
        {//varendo moléculas de C2
            
            if (descartadas[j] < 0)
            {//molecula ainda não descartada
                for (unsigned int k = 0 ; k < C1.mol[i].MSOL_expandida.size() ; k++)
                {// varrendo colunas de MSOL_expandida
                    
                    if(C1.mol[i].MSOL_expandida[k] == C2.mol[j].MSOL_expandida[k])
                    {//se aquele grupo em C1 foi igual ao mesmo mgrupo em C2 aumenta o contador
                        marcador++;
                    }
                    
                }
                if ((marcador == C1.mol[i].MSOL_expandida.size()) && (C1.mol[i].Tc == C2.mol[i].Tc))
                {//significa que a molécula i é toda igual a n
                    
                    descartadas[j] = i;
                    nMolDescartadas++;
                }
                marcador = 0;
            }
        }
    }
    return descartadas;
}

Corrente Misturador::BalancoMassa(Corrente &C1, Corrente &C2)
{
    Corrente OUT = C1;
    OUT.FLSH.vazM = C1.FLSH.vazM + C2.FLSH.vazM;
    OUT.FLSH.vazN = C1.FLSH.vazN + C2.FLSH.vazN;

    std::vector <int>  descartadas = BuscaRepetidas(C1, C2);

    for (unsigned int i = 0 ; i < descartadas.size() ; i++)
    {

        if (descartadas[i] < 0)
        {//molecula não será descartada
            OUT.mol.push_back(C2.mol[i]);
            OUT.FLSH.vazNi.push_back(C2.FLSH.vazNi[i]);
            OUT.FLSH.vazMi.push_back(C2.FLSH.vazMi[i]);
            OUT.nMol++;
        }
        else
        {// molecula será descartada
            OUT.FLSH.vazNi[descartadas[i]] += C2.FLSH.vazNi[i];
            OUT.FLSH.vazMi[descartadas[i]] += C2.FLSH.vazMi[i];
            
            // OUT.mol[descartadas[i]].vazM += C2.mol[i].vazM;
            // OUT.mol[descartadas[i]].vazN += C2.mol[i].vazN;
        }
    }
    OUT.FLSH.z.resize(OUT.mol.size());
    for (unsigned int i = 0 ; i < OUT.FLSH.z.size() ; i++) OUT.FLSH.z[i] = OUT.FLSH.vazNi[i]/OUT.FLSH.vazN;

    return OUT;
}
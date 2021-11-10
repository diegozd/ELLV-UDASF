#include "Corrente.hpp"

Corrente::Corrente(std::string path, double T, double P, double _vazN, bool aviso, std::vector < double > z) : PropsMX(path)
{
    std::vector < double > _z;
    unsigned int inicio = path.rfind("/", path.size()-2) + 1;
    unsigned int tamanho = path.size()-inicio-1;
    nome = path.substr(inicio,tamanho);
    
    std::cout << std::endl << std::endl << "####    CORRENTE " << nome << "    ####" << std::endl;
    if(z.size() != nMol) _z = LeComp(path);
    else _z = z;
    FLSH.T = T;
    FLSH.P = P;
    FLSH.vazN = _vazN;
    RemoveNulas(_z);    //remove moleculas de composição significativa e monta a composição da fase FLSH.z
    FlashELLV();
    BalancoFases(aviso);

    std::string path_descartadas = path + "descartadas.csv";
	std::ofstream fout1(path_descartadas);
	for (unsigned int i = 0 ; i < descartadas.size() ; i++)
	{
		fout1 << descartadas[i];
		fout1 << std::endl;
	}
    fout1.close();
}

void Corrente::RemoveNulas(std::vector < double > zin)
{
    std::vector < double > zout;
    descartadas.resize(zin.size(),0);
    std::vector < MoleStr > molOut;
    double tol = 1e-8;
    double SumZ = 0;

    // colocando em zout somente moléculas que possuem composição maior que a tolerancia
    for (unsigned int i = 0 ; i < zin.size() ; i++)
    {
        if (zin[i] >= tol)
        {
            zout.push_back(zin[i]);
            molOut.push_back(mol[i]);
            SumZ += zin[i];
        }
        else
        {
            descartadas[i] = 1;
            nMol--;
        }
    }

    // normalizando
    for (unsigned int i = 0 ; i < zout.size() ; i++) zout[i] = zout[i]/SumZ;

    FLSH.z = zout;
    mol = molOut;
}

void Corrente::BalancoFases(bool aviso)
{
    FLSH.vazNi.resize(nMol);
    FLSH.vazMi.resize(nMol);
    double SumV = 0;
    if (aviso == 1) std::cout << std::setprecision(2) << std::fixed;
    if (aviso == 1) std::cout << "    PROPRIEDADES DAS FASES (j) E DA MISTURA (mix):" << std::endl;
    FLSH.vazM = 0;

    for (unsigned int i = 0 ; i < nMol ; i++)
    {    
        FLSH.vazNi[i] = FLSH.z[i]*FLSH.vazN;
        FLSH.vazMi[i] = FLSH.vazNi[i]*mol[i].MMw;
        FLSH.vazM += FLSH.vazMi[i];
    }

    for (unsigned int j = 0 ; j < FLSH.F.size() ; j++)
    {
        FLSH.F[j].Z = ZPRTP(FLSH.T, FLSH.P, FLSH.F[j].EstFis, FLSH.F[j].y);
        FLSH.F[j].vazN = FLSH.F[j].beta*FLSH.vazN;
        FLSH.F[j].vazM = 0;
        FLSH.F[j].vazNi.resize(nMol);
        FLSH.F[j].vazMi.resize(nMol);
        for (unsigned int i = 0 ; i < nMol ; i++)
        {
            FLSH.F[j].vazNi[i] = FLSH.F[j].y[i]*FLSH.F[j].vazN;
            FLSH.F[j].vazMi[i] = FLSH.F[j].vazNi[i]*mol[i].MMw;
            FLSH.F[j].vazM += FLSH.F[j].vazMi[i];
        }
        FLSH.F[j].MMw = FLSH.F[j].vazM/FLSH.F[j].vazN;
        FLSH.F[j].roh = FLSH.P*FLSH.F[j].MMw/(FLSH.F[j].Z*R*FLSH.T)/1000;

        SumV += FLSH.F[j].vazM/FLSH.F[j].roh;

        if (aviso == 1) std::cout << "j: " << j << " MMw:" << FLSH.F[j].MMw << "\t roh:" << FLSH.F[j].roh << "\t VazN:" << FLSH.F[j].vazN << "\t VazM:" << FLSH.F[j].vazM << std::endl;
    }
    FLSH.MMw = FLSH.vazM/FLSH.vazN;
    FLSH.roh = FLSH.vazM/SumV;
    if (aviso == 1) std::cout << "mix: MMw:" << FLSH.MMw << "\t roh:" << FLSH.roh << "\t VazN:" << FLSH.vazN << "\t VazM:" << FLSH.vazM << std::endl;
    if (aviso == 1) std::cout << std::defaultfloat;
}

void Corrente::MarcaComoSolvente()
{
    for (unsigned int i = 0 ; i < nMol ; i++)
    {
        mol[i].SOLV = 1;
    }
}
#include "PropsMX.hpp"

// Constructor da classe
PropsMX::PropsMX(std::string path) : PRPropsPuros(path)
{
    Taij = 0;
}

PropsMX::PurosCtsStr PropsMX::CalcCtesPuros(double T, double P, bool aviso)
{
    PurosCtsStr rsult;
    CtsPRStr Ci;
    CtsPRStr Cj;

    if ((T/Taij > 1.1) || (T/Taij < 0.9))
    {
        // std::cout << "Precisei recalcular aij" << std::endl;
        GroupKij(T);
        Taij = T;
    }

    // //  bz + mCyC6 + água
    // kij[0][1] = 1.6742e-03;
    // kij[0][2] = 0.5;
    // kij[1][2] = 0.5;
    // kij[1][0] = kij[0][1];
    // kij[2][0] = kij[0][2];
    // kij[2][1] = kij[1][2];

    //  bz + nC6 + água
    // kij[0][1] = 7e-03;
    // kij[0][2] = 0.5;
    // kij[1][2] = 0.5;
    // kij[1][0] = kij[0][1];
    // kij[2][0] = kij[0][2];
    // kij[2][1] = kij[1][2];

    // bz + nC6
    // kij[0][1] = 7e-03;
    // kij[1][0] = kij[0][1];


    // kij GLP
    // kij.clear();
    // kij = LeMatriz("/home/diegozd/Documentos/PEQ/COQ712-TermoSol/ELLV/run/GLP/kij_petrosim.csv");
    
    // mostra kij
    std::cout << std::setprecision(4) << std::fixed;
    if (Jamostrei == 0)
    {
        std::cout << "kij:" <<  std::endl;
        for(unsigned int i = 0 ; i < kij.size() ; i++)
        {
            for(unsigned int j = 0 ; j < kij[i].size() ; j++)
            {
                std::cout << kij[i][j] << "\t";
            }
            std::cout << std::endl;
        } 
        Jamostrei = 1;
    }
    

    rsult.aij.resize(nMol);
    for (unsigned int i = 0 ; i < rsult.aij.size() ; i++) 
    {
        rsult.aij[i].resize(nMol);
        Ci = calcCtesPR(T,P,i);
        rsult.Cprs.push_back(Ci);
        for (unsigned int j = 0 ; j < rsult.aij[i].size() ; j++) 
        {
            Cj = calcCtesPR(T,P,j);
            rsult.aij[i][j] = pow((Ci.a_alpha*Cj.a_alpha),0.5)*(1 - kij[i][j]);
            // if ((rsult.aij[i][j] < 0) || (isnan(rsult.aij[i][j])) || (isinf(rsult.aij[i][j])))  std::cout << " Ci.a_alpha:" << Ci.a_alpha << " Cj.a_alpha:" << Cj.a_alpha << " kij[" << i << "][" << j << "]:" << kij[i][j] << std::endl;
        }
    }

    return rsult;
}

PropsMX::MxCtsStr PropsMX::CalcCtesMx(double T, double P, std::vector < double > z, bool aviso)
{
    MxCtsStr result;
    result.PCts = CalcCtesPuros(T, P, aviso);

    result.a = 0;
    result.b = 0;
    for (unsigned int i = 0 ; i < result.PCts.aij.size() ; i++) 
    {
        for (unsigned int j = 0 ; j < result.PCts.aij[i].size() ; j++) 
        {
            // if ((z[i] >= 0) && (z[j] >= 0)) 
            result.a += z[i]*z[j]*result.PCts.aij[i][j];
        }
        result.b += z[i]*result.PCts.Cprs[i].b;
    }
    result.A = result.a*P/(pow((R*T),2));
    result.B = result.b*P/(R*T);
    if (result.b != 0) result.q = result.a/(result.b*R*T);
    else result.q = 0;

    // if ((isnan(result.A)) || isnan(result.B) || isnan(result.q) || result.A <=0 ) 
    // std::cout << " A:" << result.A << " B:" << result.B << " q:" << result.q << " a:" << result.a << " b:" << result.b << " T:" << T << " P:" << P << std::endl;

    return result;
}

double PropsMX::ZPRTP(double T, double P, bool EstFis_z, std::vector < double > z)
{
    MxCtsStr cts = CalcCtesMx(T, P, z);

    double Z = 0;
    
    // if (EstFis_z == 1) Z = FZPR(1, cts.A, cts.B);    //vapor
    // else Z = FZPR(cts.B, cts.A, cts.B);    //liquido
    if (EstFis_z == 1) Z = FZVPR(cts.A, cts.B, cts.q);    //vapor
    else Z = FZLPR(cts.A, cts.B, cts.q);    //liquido

    return Z;
}

std::vector < double > PropsMX::phiMxi(double T, double P, bool EstFis_z, std::vector < double > z, bool aviso)
{
    MxCtsStr cts = CalcCtesMx(T, P, z, aviso);  
    std::vector < double > phimx;
    double Zl = 0;
    double Zv = 0;
    double Z = 0;
    double I = 0;
    double qi_barra = 0;
    double ai_barra = 0;
    double bi = 0;
    double Szaij = 0;
    double sig = 1+pow(2,0.5);
    double eps = 1-pow(2,0.5);

    Zv = FZVPR(cts.A, cts.B, cts.q);    //vapor
    Zl = FZLPR(cts.A, cts.B, cts.q);    //liquido
    
    // if ((isnan(Zv)) || isnan(Zl) || isinf(Zv)|| isinf(Zl)) std::cout << "Zv:" << Zv << " Zl:" << Zl << " A:" << cts.A << " B:" << cts.B << " q:" << cts.q << std::endl;

    if (Zv == Zl) 
    {
        if (Zv > 0.5) 
        {
            if (EstFis_z == 0)
            {
                std::cout << "sem raiz de líquido Z:" << Zv << " A:" << cts.A << " B:" << cts.B << " a:" << cts.a << " b:" << cts.b << std::endl;
                Zl = cts.B; 
            } 

        }
        else 
        {
            if (EstFis_z == 1) 
            {
                std::cout << "sem raiz de vapor Z = " << Zl << std::endl;
                Zv = 1;
            }
        }
    }
    

    if (EstFis_z == 1) Z = Zv;      //vapor
    else Z = Zl;                    //liquido

    phimx.resize(nMol);

    for (unsigned int i = 0 ; i < phimx.size() ; i++)
    {    
        // if((cts.b <= 0) || (cts.a == 0) || (fail == 1))
        // if (fail == 1)
        // {
        //     phimx[i] = failPhiLiq;
        // }
        // else if (fail == 2)
        // {
        //     phimx[i] = 1;
        // }
        // else 
        if((cts.b <= 0) || (cts.a == 0)  )
        {
            // std::cout << "  i:" << i << " bi:" <<  bi << " b:" << cts.b <<  " B:" << cts.B << " Szaij:" << Szaij << " a:" << cts.a << " A:" << cts.A <<" q:" << cts.q << " Z:" << Z <<  " ai_barra:" << ai_barra <<  " qi_barra:" << qi_barra << " I:" << I << " phimx[" << i << "]:" <<  phimx[i] << " log(phimx):" << log(phimx[i]) << std::endl;
            if (EstFis_z == 1) phimx[i] = 1;
            else phimx[i] = 0+1E-10;
        } 
        else
        {
            if ((Z-cts.B) <= 0) Z = cts.B+1E-10;

            bi = cts.PCts.Cprs[i].b;
            Szaij = 0;
            for (unsigned int j = 0 ; j < phimx.size() ; j++) Szaij += z[j]*cts.PCts.aij[i][j];
            ai_barra = 2*Szaij-cts.a;
            qi_barra = cts.q*(1 + ai_barra/cts.a - bi/cts.b);               //no livro do van ness é qi_barra = q*(1 + ai_barra/a + bi/b); na apostila da materia é qi_barra = q*(1 + ai_barra/a - bi/b);
            I = (1/(sig-eps))*log((Z+sig*cts.B)/(Z+eps*cts.B));
            phimx[i] = exp((bi/cts.b)*(Z-1) - log(Z-cts.B) - qi_barra*I);
        }
        // std::cout << "  i:" << i << " bi:" <<  bi << " b:" << cts.b <<  " B:" << cts.B << " Szaij:" << Szaij << " a:" << cts.a << " A:" << cts.A <<" q:" << cts.q << " Z:" << Z <<  " ai_barra:" << ai_barra <<  " qi_barra:" << qi_barra << " I:" << I << " phimx[" << i << "]:" <<  phimx[i] << " log(phimx):" << log(phimx[i]) << std::endl;
    }

    return phimx;
}

PropsMX::FaseStr PropsMX::tpd(bool EstFis_y, std::vector < double > y0, bool EstFis_z, FlshStr bulk)
{

    FaseStr FTest;
    FTest.y.resize(y0.size());
    std::vector < double > phi_z = phiMxi(bulk.T,bulk.P,EstFis_z,bulk.z);
    std::vector < double > Y(y0.size());
    std::vector < double > h(y0.size());
    
    double tol = 1e-8;
    double SerroQ = 0;
    double SerroQ_old = 1e10;
    bool PrimeiraVez = 0;
    bool suave = 0;
    double SumY = 0;
    double K = 0;
    unsigned int iter = 0;


    for (unsigned int i = 0 ; i < y0.size() ; i++) h[i] = log(bulk.z[i]) + log(phi_z[i]);

    FTest.y = y0;
    do
    {
        FTest.phi = phiMxi(bulk.T,bulk.P,EstFis_y,FTest.y);
        SumY = 0;
        SerroQ = 0;

        for (unsigned int i = 0 ; i < y0.size() ; i++) 
        {
            Y[i] = exp(h[i] - log(FTest.phi[i]));
            SumY += Y[i];
        }
        for (unsigned int i = 0 ; i < y0.size() ; i++)
        {
            y0[i] = Y[i]/SumY;
            SerroQ += pow((FTest.y[i]-y0[i]),2);
            if (suave == 0) 
            {
                FTest.y[i] = y0[i];
                if (iter >= 50) suave = 1;
            }
        }
        if (suave == 1)
        {

            if ((SerroQ_old < SerroQ) && (PrimeiraVez == 1)) SerroQ = 0;
            else 
            {
                for (unsigned int i = 0 ; i < y0.size() ; i++)
                {
                    // std::cout << "FTest.y[" << i << "]:" << y0[i] << std::endl;
                    FTest.y[i] = (y0[i] + FTest.y[i])/2;
                }
                // std::cout << "iter:"<< iter << " SerroQ:" << SerroQ << " SerroQ_old:" << SerroQ_old << " SumY:" << SumY << " tpd:" << -log(SumY)*R*bulk.T << std::endl;
                PrimeiraVez = 1;
            }
            
        }

        if ((isnan(SumY)) || (isinf(SumY)))
        {
            std::cout << "deu nan em tpd" << std::endl;
            SerroQ = tol;
            SumY = 1;
        }
        if (SerroQ_old == SerroQ) suave = 1;
        // if (iter >= 100) std::cout << "SerroQ:" << SerroQ << std::endl; 
        
        SerroQ_old = SerroQ;
    
        iter++;
        if (iter >= maxiter*2) std::cout << "tpd atingiu maxiter. Sauve:" << suave << std::endl;
    } while ((SerroQ > tol) && (iter < maxiter*2));

    K = -log(SumY);
    FTest.tpd = K*R*bulk.T;
    FTest.EstFis = EstFis_y;
    FTest.beta = 0;

    // std::cout << "tpd:" << FTest.tpd <<  std::endl; 

    return FTest;
}

std::vector < PropsMX::FaseStr > PropsMX::AnaliseEstabilidade(FlshStr bulk, bool EstFis_z)
{
    FaseStr FTest;
    std::vector < FaseStr > PosFasEst;
    unsigned int n = 0;
    unsigned int purosTPDnega = 0;
    bool EstFis_y = EstFis_z;
    std::vector < double > XL(bulk.z.size());
    std::vector < double > y(bulk.z.size(),0);
    std::vector < double > yXL(bulk.z.size(),0);
    double tol = 1e-10;
    double tolx = 1e-08;
    int origem1My = 0;
    std::vector < double > yQueDeuOrigemMy;
    std::vector < double > AddMy(bulk.z.size()+1);
    std::vector < std::vector < double > > My;     //vetor que concentra todas as possiveis fases líquidas

    //montando a fase referencia
    FTest.y = bulk.z;
    FTest = tpd(EstFis_y, bulk.z, EstFis_z, bulk);
    // std::cout << "    AnaliseEstabilidade-referencia tpd:" << FTest.tpd << " EstFis:" << FTest.EstFis << " y[1]:" << FTest.y[1] << std::endl;
    // FTest.phi = phiMxi(bulk.T,bulk.P,FTest.EstFis,FTest.y);
    // FTest.tpd = 0;
    if (EstFis_z == 0) 
    {
        for (unsigned int i = 0 ; i < FTest.y.size() ; i++) AddMy[i] = FTest.y[i];
        AddMy[bulk.z.size()] = FTest.tpd;
        n++;
        origem1My++;
        My.push_back(AddMy);
        yQueDeuOrigemMy.push_back(0);
    }
    PosFasEst.push_back(FTest);

    // testando estabilidade da mistura na fase oposta
    EstFis_y = 1 - EstFis_z;
    FTest = tpd(EstFis_y, bulk.z, EstFis_z, bulk);
    // std::cout << "    AnaliseEstabilidade-opsta tpd:" << FTest.tpd << " EstFis:" << FTest.EstFis << " y[1]:" << FTest.y[1] << std::endl;
    if (FTest.tpd < -tol) PosFasEst.push_back(FTest);

    // // testando estabilidade da mistura na fase vapor
    // EstFis_y = 1;
    // FTest = tpd(EstFis_y, bulk.z, EstFis_z, bulk);
    // std::cout << "    AnaliseEstabilidade-0 FTest.tpd vapor = " << FTest.tpd << std::endl;
    // if (FTest.tpd < -tol) PosFasEst.push_back(FTest);

    // testando estabilidade da mistura na fase líquida
    // EstFis_y = 0;
    // FTest = tpd(EstFis_y, bulk.z, EstFis_z, bulk);
    // std::cout << "    AnaliseEstabilidade-1 FTest.tpd líquido = " << FTest.tpd << std::endl;
    // if (FTest.tpd < -tol) 
    // {
    //     for (unsigned int i = 0 ; i < FTest.y.size() ; i++) AddMy[i] = FTest.y[i];
    //     AddMy[bulk.z.size()] = FTest.tpd;
    //     n++;
    //     origem1My++;
    //     My.push_back(AddMy);
    //     yQueDeuOrigemMy.push_back(0);
    // }

    EstFis_y = 0;
    // testando se cada um dos componentes puros pode fazer uma fase líquida
	for (unsigned int i = 0 ; i < bulk.z.size() ; i++) 
	{
		if (i == 0) y[i] = 1;
		else 
		{
			y[i-1] = 0;
			y[i] = 1;
		}
        FTest = tpd(EstFis_y, y, EstFis_z, bulk);
        // std::cout << "    AnaliseEstabilidade-líquido puro "<< i <<" tpd:" << FTest.tpd << " EstFis:" << FTest.EstFis << " y[1]:" << FTest.y[1] <<std::endl;
		if (FTest.tpd < -tol) 
        {
            XL[i] = 1;
            for (unsigned int j = 0 ; j < bulk.z.size() ; j++) AddMy[j] = FTest.y[j];
            AddMy[bulk.z.size()] = FTest.tpd;
            purosTPDnega++;
            My.push_back(AddMy);
            yQueDeuOrigemMy.push_back(i);
        }
	}

    // juntando todas os puros que quiseram fazer outra fase numa mistura só
    for (unsigned int i = 0 ; i < XL.size() ; i++)
    {
        if(XL[i] == 1) yXL[i] = (double) 1/purosTPDnega;
        else yXL[i] = 0;
    }

    FTest = tpd(EstFis_y, yXL, EstFis_z, bulk);
    if (FTest.tpd  < -tol) 
    {
        for (unsigned int j = 0 ; j < bulk.z.size() ; j++) AddMy[j] = FTest.y[j];
        AddMy[bulk.z.size()] = FTest.tpd;
        n++;
        My.push_back(AddMy);
        yQueDeuOrigemMy.push_back(-1);
    }

    n += purosTPDnega;
    
    // vendo se ha diferença entre as composições geradas por cada fase valida
    bool Npodei = 0;                // variavel que indica se a mistura j pode ser representada pela i
    unsigned int ncabecas = 0;      //variavel que indica o numero total de misturas cabeça
    std::vector < int > validas(n,-1);
    for (unsigned int i = 0 ; i < n ; i++)
    {// varrendo correntes validas
        // std::cout << "My[" << i << "][tpd]:" << My[i][bulk.z.size()] << std::endl;
        if ((ncabecas == 0) || (validas[i] < 0))
        {// ainda não há cabeças OU a mistura não pode ser representada por uma cabeça
            validas[i] = i;
            ncabecas++;

            for (unsigned int j = i+1 ; j < n ; j++) 
            {// varrendo as proximas misturas 
                Npodei = 0;
                for (unsigned int k = 0 ; k < My[j].size() ; k++)
                {
                    // std::cout << "Erro y:" << My[i][k] << " " << My[j][k] << " " << pow((My[i][k] - My[j][k]),2) << std::endl;
                    // if (abs(My[i][k] - My[j][k]) > tolx) 
                    if (pow((My[i][k] - My[j][k]),2) > tolx) 
                    {
                        Npodei = 1;
                        k = My[j].size();
                    }
                }
                if (Npodei == 0) validas[j] = i; //mistura j pode ser representada pela i
            }
        }
    }

    unsigned int checkcabecas = 0;
    int compnowarning = 0;
    for (unsigned int i = 0 ; i < validas.size() ; i++) 
    {
        // std::cout << "validas[" << i << "]:" << validas[i] << std::endl;
        compnowarning = i;
        if (validas[i] == compnowarning)
        {
            checkcabecas++;
            if (yQueDeuOrigemMy[i] == -1) y = yXL;              //quem deu origem esse TPD foi a mistura de todos os puros que queria fazer outra fase
            else if ((i == 0) && (origem1My == 1)) y = bulk.z;  //quem deu origem ao primeiro my foi composição bulk liquido
            else                                                //quem deu origem a esse TPD foi algum puro
            {                                                   
                for (unsigned int j = 0 ; j < bulk.z.size() ; j++) 
                {
                    if (j == yQueDeuOrigemMy[i]) y[j] = 1;
                    else y[j] = 0;
                }
            }
            FTest = tpd(EstFis_y, y, EstFis_z, bulk);
		    if (FTest.tpd < -tol) PosFasEst.push_back(FTest);
        }
    }

    // MostraVetorInt(validas);
    // std::string path = "/home/diegozd/Documentos/PEQ/COQ712-TermoSol/ELLV/run/RV/";
    // std::string filename = "My.csv";
    // EscreveMatriz(My,path,filename);

    return PosFasEst;
}

void PropsMX::FlashELLV(double T, double P, std::vector < double > z)
{

    // std::cout << std::endl << "      Inicio FlashELLV: " << std::endl;

    FlshStr _FLSH;
    FlshStr BFLS;
    Qstr Qnew;
    Qstr BestQ;
    int keyF = 0;
    unsigned int iter = 0;
    unsigned int itery;
    unsigned int iterAlpha;
    unsigned int Tentativas = 0;
    int inativeF = -1;
    int inativeFold = 0;
    bool convergiu = 1;
    bool convergiuQ = 1;
    bool NaoPode = 0;
    double MinTPD = 10000;
    double MaxBeta;
    double SumBeta;
    double alpha = 1;
    double SumZ;
    double SumYquad;
    double minBetaNeg;
    double lastBQ = 0;
    std::vector < double > SumZiEiphiIJ;
    std::vector < double > SumY;
    std::vector < std::vector < double > > y_old;
    std::vector < double > SumBetaK;
    std::vector < FaseStr > VFina;                        // vetor de FaseStrs inativas
    std::vector < int > repetidas;
    
    if (T == 0) 
    {
        T = FLSH.T;
        _FLSH.T = FLSH.T;
    }
    else _FLSH.T = T;
    if (P == 0) 
    {
        P = FLSH.P;
        _FLSH.P = FLSH.P;
    }
    else _FLSH.P = P;
    if(z.size() != nMol) 
    {
        z = FLSH.z;
        _FLSH.z = FLSH.z;
    }
    else _FLSH.z = z;
    
    Taij = 0;
    // analise de estabilidade propondo fase líquida
    _FLSH.F = AnaliseEstabilidade(_FLSH);

    // exibindo informações sobre AnaliseEstabilidade
    // std::cout << "Nfases = " << _FLSH.F.size() << std::endl;
    // for (unsigned int i = 0 ; i < _FLSH.F.size() ; i++) 
    // {
    //     std::cout << "    tpd = " << _FLSH.F[i].tpd << " _FLSH.F[" << i << "].EstFisFT_y = " << _FLSH.F[i].EstFis << " y = " ;
    //     SumZ = 0;
    //     for(unsigned int j = 0 ; j < _FLSH.F[i].y.size() ; j++) 
    //     {
    //         if (j < 5 ) std::cout << _FLSH.F[i].y[j] << " ";
    //         else if (j == 5) std::cout << "... ";
    //         SumZ += _FLSH.F[i].y[j];
    //     }
    //     std::cout  << " Sum = " << SumZ << std::endl;
    // }
    

    SumBetaK.resize(z.size());
    _FLSH.K.resize(z.size());
    SumZiEiphiIJ.resize(_FLSH.F.size());
    Eigen::MatrixXd g(_FLSH.F.size(),1);
    Eigen::MatrixXd dBeta(_FLSH.F.size(),1);
    Eigen::MatrixXd EiBeta(_FLSH.F.size(),1);
    Eigen::MatrixXd EiBetaNew(_FLSH.F.size(),1);
    Eigen::MatrixXd H(_FLSH.F.size(),_FLSH.F.size());
    
    //estimativa inicial de beta (homogênea)
    // std::cout << std::endl << " estimativa inicial de beta:" ;
    for (unsigned int j = 0 ; j < _FLSH.F.size() ; j++) 
    {//varrendo as fases
        _FLSH.F[j].beta = (double) 1/_FLSH.F.size();
        if (_FLSH.F[j].tpd < MinTPD) 
        {
            MinTPD = _FLSH.F[j].tpd;
            // keyF = j;
        }
        // std::cout << _FLSH.F[j].beta << " ";
    }
    // std::cout << std::endl;

    // i) avaliar Q
    BFLS = _FLSH;
    BestQ = CalcQ(_FLSH);

    do
    {
        g.resize(_FLSH.F.size(),1);
        dBeta.resize(_FLSH.F.size(),1);
        EiBeta.resize(_FLSH.F.size(),1);
        EiBetaNew.resize(_FLSH.F.size(),1);
        H.resize(_FLSH.F.size(),_FLSH.F.size());
        SumY.resize(_FLSH.F.size());
        y_old.resize(_FLSH.F.size());

        // std::cout << std::endl << std::endl << "            ITER:" << iter << std::endl; 

        // ii) calcular grandiente e Hessiana e resolver o problema de newton
        for (unsigned int j = 0 ; j < _FLSH.F.size() ; j++)
        {
            EiBeta(j,0) = _FLSH.F[j].beta;
            SumZiEiphiIJ[j] = 0;
            for (unsigned int i = 0 ; i < z.size() ; i++) SumZiEiphiIJ[j] += z[i]/(BestQ.E[i]*_FLSH.F[j].phi[i]);
            g(j,0) = 1 - SumZiEiphiIJ[j];
        }
        // std::cout << "g:\n" << g << std::endl << std::endl;

        for (unsigned int j = 0 ; j < _FLSH.F.size() ; j++)
        {
            for (unsigned int k = 0 ; k < _FLSH.F.size() ; k++)
            {
                H(j,k) = 0;
                for (unsigned int i = 0 ; i < z.size() ; i++) 
                {
                    H(j,k) += z[i]/(pow(BestQ.E[i],2)*_FLSH.F[j].phi[i]*_FLSH.F[k].phi[i]);
                }
                if (k == j) H(j,k) += 1e-10;
            }
        }
        // std::cout << "H:\n" << H << std::endl<< std::endl;

        // iii) calcular beta novo e ajustar alpha se algum beta < 0
        alpha = 1;
        dBeta = -H.inverse()*g;
        EiBetaNew = EiBeta + alpha*dBeta;
        // std::cout << "H.inverse:\n" << H.inverse() << std::endl;
        // std::cout << std::endl << "EiBeta:\n" << EiBeta << std::endl;
        // std::cout << std::endl << "dBeta:\n" << dBeta << std::endl;
        // std::cout << std::endl << "EiBetaNew:\n" << EiBetaNew << std::endl;
        // std::cout << std::endl << "check:\n" << H*dBeta + g << std::endl;
        // std::cout << std::endl;
        
        minBetaNeg = 100;
        //ajustando alpha caso algum beta seja negativo
        for (unsigned int k = 0 ; k < _FLSH.F.size() ; k++) 
        {
            if ((EiBetaNew(k,0) < 0) && (EiBetaNew(k,0) < minBetaNeg))
            {
                alpha = -0.999*EiBeta(k,0)/dBeta(k,0);
                minBetaNeg = EiBetaNew(k,0);
            } 

        }
        // std::cout << "alpha0:" << alpha << std::endl;
        
        if (alpha != 1) alpha = 0.999*alpha;
        // double alpha0 = alpha;
        // double dalpha, Ql;
        double Qvelho;

        iterAlpha = 0;
        convergiuQ = 1;
        do
        {
            
            EiBetaNew = EiBeta + alpha*dBeta;
            // std::cout << iterAlpha << " EiBetaNew:\n" << EiBetaNew << std::endl;
            
            MaxBeta = 0;
            inativeF = -1;
            SumBeta = EiBetaNew.sum();
            for (unsigned int j = 0 ; j < _FLSH.F.size() ; j++) 
            {
                // normalizando beta
                // _FLSH.F[j].beta = EiBetaNew(j,0);
                _FLSH.F[j].beta = EiBetaNew(j,0)/SumBeta;
                // std::cout << _FLSH.F[j].beta << " ";

                if (_FLSH.F[j].beta > MaxBeta) 
                {
                    MaxBeta = _FLSH.F[j].beta;
                    keyF = j;
                }
                if (_FLSH.F[j].beta <= 1e-8) 
                {
                    _FLSH.F[j].beta = 0;
                    inativeF = j;
                }
            }

            // iv) recalcula Q
            Qvelho = Qnew.Q;
            Qnew = CalcQ(_FLSH);

            if ((Qnew.Q <= BestQ.Q) || (Qnew.Q == Qvelho))
            {
                convergiuQ = 0;
                BestQ = Qnew;
                BFLS = _FLSH;
            }
            else 
            {
            //     if (iterAlpha < 30) 
            //     {
                    // alpha0 = alpha;       // alpha que gereou Qnew
                    alpha = alpha/2;      // proximo alpha
                // }
                // else
                // {
                //     //primeira vez que entra aqui, alpha0 tem valor que gerou Qdalpha
                    
                //     dalpha = (alpha-alpha0);
                //     Ql = (Qnew.Q - Qdalpha)/dalpha;
                //     alpha0 = alpha; 
                //     alpha = alpha - Ql/Qdalpha;
                //     std::cout << "alpha:" << alpha << " alpha0:" << alpha0 << " dalpha:" << dalpha << " Ql:" << Ql << " Qnew.Q:" << Qnew.Q << std::endl;
                // }
            }
            
            iterAlpha++;
            if (iterAlpha >= 100) 
            {
                std::cout << "Aingiu limite iteracoes em Q na função ELLV" << std::endl;
                convergiuQ = 0;
            }
            if (inativeF >= 0) convergiuQ = 0;
            if (_FLSH.F.size() == 1) convergiuQ = 0;
            

        } while (convergiuQ == 1);

        // mostra infos na tel
        // std::cout << "iterAlpha:" << iterAlpha << " alpha:" << alpha << std::endl;
        // std::cout << "Qnew:" << Qnew.Q << " BestQ:" << BestQ.Q << std::endl;
        // std::cout << "Beta:";
        // for (unsigned int j = 0 ; j < _FLSH.F.size() ; j++) std::cout << _FLSH.F[j].beta << " ";
        // std::cout << std::endl;
        
        if (inativeF >= 0) 
        {// pegou o cara inativo
            // std::cout << "inativeF:" << inativeF << " g:(inativeF,0):" << g(inativeF,0) << std::endl;
            if (g(inativeF,0) < 0)
            {
                VFina.push_back(_FLSH.F[inativeF]);
                // std::cout << "g < 0 => VFina.size():" << VFina.size() << std::endl;
            }
            _FLSH.F.erase(_FLSH.F.begin()+inativeF);
            if(keyF > inativeF) keyF--;
            // std::cout << "FASE " << inativeF << " RETIRADA  =>  Nfases:" << _FLSH.F.size() << " keyF:" << keyF <<  std::endl;
        }

        // convergindo composição das fases 
        itery = 0;
        do
        {
            // std::cout << "      itery:" << itery << " keyF:" << keyF << std::endl;
            for (unsigned int i = 0 ; i < z.size() ; i++) 
            {
                SumBetaK[i] = 0;
                _FLSH.K[i].resize(_FLSH.F.size());
                for (unsigned int j = 0 ; j < _FLSH.F.size() ; j++) 
                {
                    _FLSH.K[i][j] = _FLSH.F[keyF].phi[i]/_FLSH.F[j].phi[i];
                    if (isnan(_FLSH.K[i][j])) _FLSH.K[i][j] = 1;
                    else SumBetaK[i] += _FLSH.F[j].beta*(_FLSH.K[i][j] - 1);
                }
            }
            
            SumYquad = 0;
            for (unsigned int j = 0 ; j < _FLSH.F.size() ; j++) 
            {
                y_old[j] = _FLSH.F[j].y;
                SumY[j] = 0;
                for (unsigned int i = 0 ; i < z.size() ; i++) 
                {
                    _FLSH.F[j].y[i] = z[i]*_FLSH.K[i][j]/(1 + SumBetaK[i]);
                    if (isnan(_FLSH.F[j].y[i])) _FLSH.F[j].y[i] = 0;
                    SumY[j] += _FLSH.F[j].y[i];
                }
                for (unsigned int i = 0 ; i < z.size() ; i++) 
                {
                    _FLSH.F[j].y[i] = _FLSH.F[j].y[i]/SumY[j];
                    SumYquad += pow((_FLSH.F[j].y[i] - y_old[j][i]),2);
                }
                _FLSH.F[j].phi = phiMxi(T,P,_FLSH.F[j].EstFis,_FLSH.F[j].y);
            }
            itery++;
        } while ((itery < 20) && (SumYquad > 1e-4));
        // std::cout << "itery:" << itery << " SumYquad:" << SumYquad << std::endl;

        Qnew = CalcQ(_FLSH);
        // std::cout << "iter:" << iter << " Qnew:" << Qnew.Q << " BestQ:" << BestQ.Q << std::endl;

        if (Qnew.Q <= BestQ.Q)
        {
            BestQ = Qnew;
            BFLS = _FLSH;
        }

        // iv) após a convergência checa se alguma fse desativada tem gradiente negativo, se tiver reativa a fase
        convergiu = ((abs(lastBQ - BestQ.Q) > 1e-10) && (_FLSH.F.size() > 1));
        if ((inativeF >= 0) && (convergiu == 0))
        {
            BestQ = Qnew;
            BFLS = _FLSH;
            convergiu = 1;
        }
        if ((VFina.size() > 0) && (inativeFold == 0)) inativeFold = 1;
        else if (inativeFold == 1)
        {
            inativeFold = 0;
            // std::cout << "RETOMANDO A FASE COM g < 0" << std::endl;
            _FLSH.F.push_back(VFina[0]);
            VFina.erase(VFina.begin()+0);
            if(keyF > inativeF) keyF++;
            convergiu = 1;
        }
        // std::cout << "convergiu:" << convergiu << " Qdiff:" << abs(lastBQ - BestQ.Q) << std::endl;
        
        lastBQ = BestQ.Q;
        iter++;

        if (convergiu == 0)
        {
            //ajustando composição e beta para uma unica fase
            if (BFLS.F.size() == 1)
            {
                BFLS.F[0].y = z;
                BFLS.F[0].phi = phiMxi(T,P,BFLS.F[0].EstFis,BFLS.F[0].y);
                BFLS.F[0].beta = 1;
                BFLS.F[0].Z = ZPRTP(T,P,BFLS.F[0].EstFis,BFLS.F[0].y);
            }

            //removendo nulas e normalizando
            SumY.resize(BFLS.F.size());
            for (unsigned int j = 0 ; j < BFLS.F.size() ; j++) 
            {
                SumY[j] = 0;
                for (unsigned int i = 0 ; i < z.size() ; i++) 
                {   
                    // removendo nulas
                    if (BFLS.F[j].y[i] > 1e-10)
                    {
                        SumY[j] += BFLS.F[j].y[i];
                    }
                    else BFLS.F[j].y[i] = 0;
                }  
                SumZ = 0;
                for (unsigned int i = 0 ; i < z.size() ; i++) 
                {
                    //normalizando
                    BFLS.F[j].y[i] = BFLS.F[j].y[i]/SumY[j];
                    SumZ += BFLS.F[j].y[i];
                }
            }

            
            // verificando se as fases se repetem
            repetidas.resize(BFLS.F.size(),-1);
            for (unsigned int j = 0 ; j < BFLS.F.size() ; j++) 
            {
                SumY[j] = 0;
                if(repetidas[j] < 0)
                {
                    for (unsigned int k = j+1 ; k < BFLS.F.size() ; k++)
                    {
                        if (BFLS.F[j].EstFis == BFLS.F[k].EstFis)
                        {
                            NaoPode = 0;
                            for (unsigned int i = 0 ; i < z.size() ; i++) 
                            {
                                // std::cout << j << " " << k << " " << BFLS.F[j].y[i] << " " << BFLS.F[k].y[i] <<  " " <<  abs(BFLS.F[j].y[i] - BFLS.F[k].y[i]) << std::endl;   
                                if (pow((BFLS.F[j].y[i] - BFLS.F[k].y[i]),2) > 1e-03)
                                {
                                    NaoPode = 1;
                                    i = z.size();
                                }
                            }
                            if (NaoPode == 0) repetidas[k] = j;
                        }
                    } 
                }
            }

            //removendo fases repetidas
            for (unsigned int j = BFLS.F.size()-1 ; j > 0 ; j--)
            {
                if (repetidas[j] >= 0) 
                {
                    std::cout << " ##########  Fase " << j << " eh igual a fase " << repetidas[j] <<  " ########" <<  std::endl;
                    warning++;
                    BFLS.F[repetidas[j]].beta += BFLS.F[j].beta;
                    BFLS.F.erase(BFLS.F.begin()+j);
                }
            } 

            //removendo fases nulas
            for (unsigned int j = BFLS.F.size()-1 ; j > 0 ; j--)
            {
                if (BFLS.F[j].beta == 0) 
                {
                    BFLS.F.erase(BFLS.F.begin()+j);
                }
            } 

            //Verificando estabilidade
            _FLSH.F.clear();
            _FLSH.F.resize(BFLS.F.size());
            _FLSH.T = BFLS.T;
            _FLSH.P = BFLS.P;
            _FLSH.z = BFLS.F[0].y;
            for (unsigned int j = 0 ; j < BFLS.F.size() ; j++) 
            {
                _FLSH.F[j] = tpd(BFLS.F[j].EstFis, BFLS.F[j].y, BFLS.F[0].EstFis, _FLSH);
                _FLSH.F[j].beta = BFLS.F[j].beta;
                BFLS.F[j].tpd = _FLSH.F[j].tpd;
            }
            // _FLSH.z = BFLS.z;
            _FLSH.F = AnaliseEstabilidade(_FLSH, BFLS.F[0].EstFis); 
            

            // _FLSH.F.erase(_FLSH.F.begin()+0);
            // if (_FLSH.F.size() > 1) _FLSH.F.erase(_FLSH.F.begin()+0);
            for (unsigned int j = 0 ; j < _FLSH.F.size() ; j++) _FLSH.F[j].beta = (double) 1.0/_FLSH.F.size();
            if (BFLS.F.size() < _FLSH.F.size())
            // if (_FLSH.F.size() > 0)
            {
                // exibindo informações sobre AnaliseEstabilidade
                std::cout << "Nfases = " << _FLSH.F.size() << std::endl;
                for (unsigned int i = 0 ; i < _FLSH.F.size() ; i++) 
                {
                    std::cout << "    tpd = " << _FLSH.F[i].tpd << " _FLSH.F[" << i << "].EstFisFT_y = " << _FLSH.F[i].EstFis << " y = " ;
                    SumZ = 0;
                    for(unsigned int j = 0 ; j < _FLSH.F[i].y.size() ; j++) 
                    {
                        if (j < 5 ) std::cout << _FLSH.F[i].y[j] << " ";
                        else if (j == 5) std::cout << "... ";
                        SumZ += _FLSH.F[i].y[j];
                    }
                    std::cout  << " Sum = " << SumZ << std::endl;
                }

                if (Tentativas < 0) 
                {
                    std::cout << "Tentativa:" << Tentativas << std::endl;
                    for (unsigned int j = 0 ; j < BFLS.F.size() ; j++) std::cout << "BFLS.F[" << j << "].beta:" << BFLS.F[j].beta << std::endl;
                    convergiu = 1;      //força NÃO convergência
                    Tentativas++;
                }
                else
                {
                    std::cout << "Aviso de instabilidade. Verificar tpds:" << std::endl;
                    for (unsigned int j = 0 ; j < BFLS.F.size() ; j++) std::cout << "BFLS.F[" << j << "].tpd = " << BFLS.F[j].tpd << "\t EstFis:" << BFLS.F[j].EstFis << "\t beta:" << BFLS.F[j].beta << "\t y[0]:" << BFLS.F[j].y[0] <<  std::endl;
                    for (unsigned int j = 0 ; j < _FLSH.F.size() ; j++) std::cout << "_FLSH.F[" << j << "].tpd = " << _FLSH.F[j].tpd << "\t EstFis:" << _FLSH.F[j].EstFis << "\t beta:" << _FLSH.F[j].beta << "\t y[0]:" << _FLSH.F[j].y[0] << std::endl;
                    std::cout << std::endl;
                    
                }
            }
        }

        if(iter >= 1000) std::cout << "Atingiu maximo de iteracoes em ELLV" << std::endl;
        
    } while ((convergiu == 1) && (iter < 1000));

    
    std::cout << std::endl << "    RESULTADOS ELLV:    " << std::endl;
    std::cout << std::setprecision(1) << std::fixed;
    std::cout << "T:      " << BFLS.T << " K" << std::endl;
    std::cout << "P:      " << BFLS.P << " bar" << std::endl;
    std::cout << "Nmol:   " << z.size() << std::endl;
    std::cout << "Nfases: " << BFLS.F.size() << std::endl;

    std::cout << "EstFis (0->Liquido ; 1->Vapor):" << std::endl << "    ";
    for (unsigned int j = 0 ; j < BFLS.F.size() ; j++) std::cout << BFLS.F[j].EstFis << " \t\t";
    std::cout << std::endl;
    
    std::cout << std::setprecision(6);
    std::cout << "Z :" << std::endl << "    ";
    for (unsigned int j = 0 ; j < BFLS.F.size() ; j++) 
    {
        BFLS.F[j].Z = ZPRTP(BFLS.T, BFLS.P, BFLS.F[j].EstFis, BFLS.F[j].y);
        std::cout << BFLS.F[j].Z << " \t";
    }
    std::cout << std::endl;

    std::cout << "beta:" << std::endl << "    ";
    for (unsigned int j = 0 ; j < BFLS.F.size() ; j++) std::cout << BFLS.F[j].beta << " \t";
    std::cout << std::endl;

    std::cout << "y:" << std::endl;
    for (unsigned int i = 0 ; i < z.size() ; i++) 
    {
        if (i < 10) std::cout << "0" << i << ": ";
        else std::cout << i << ": ";
        for (unsigned int j = 0 ; j < BFLS.F.size() ; j++) 
        {   
            if(i == 0) SumY[j] = 0;
            std::cout << BFLS.F[j].y[i] << " \t";
            SumY[j] += BFLS.F[j].y[i];
        }
        std::cout << std::endl; 
    }
    std::cout << "Sum y:"<< std::endl; 
    for (unsigned int j = 0 ; j < BFLS.F.size() ; j++) std::cout << SumY[j] << " \t";

    std::cout << std::endl << std::endl;
    std::cout << std::defaultfloat;

    Taij = 0;
    BFLS.vazN = FLSH.vazN;
    FLSH = BFLS;
}

PropsMX::Qstr PropsMX::CalcQ(FlshStr _FLSH)
{
    Qstr out;
    out.E.resize(_FLSH.z.size());
    double SumZiLnE = 0;
    double MaxE = 0;
    double SumBeta = 0;
    std::vector < bool > deuinf(_FLSH.z.size(),0);

    for (unsigned int j = 0 ; j < _FLSH.F.size() ; j++) SumBeta += _FLSH.F[j].beta;

    for (unsigned int i = 0 ; i < _FLSH.z.size() ; i++)
    {//varrendo moleculas
        out.E[i] = 0;
        for (unsigned int k = 0 ; k < _FLSH.F.size() ; k++) 
        {//varrendo fases
            out.E[i] += _FLSH.F[k].beta/_FLSH.F[k].phi[i];
            if (isinf(out.E[i])) 
            {
                deuinf[i] = 1;
                k = _FLSH.F.size();
            }
            else if (out.E[i] > MaxE) MaxE = out.E[i];
        }
        SumZiLnE += _FLSH.z[i]*log(out.E[i]);
    }

    for (unsigned int i = 0 ; i < _FLSH.z.size() ; i++) 
    {//varrendo moleculas
        if (deuinf[i] == 1) out.E[i] = 100*MaxE;
    }
        
    out.Q = SumBeta - SumZiLnE;

    if (isnan(out.Q))
    {
        std::cout << "Deu nan em Q. SumBeta:" << SumBeta << " SumZiLnE:" << SumZiLnE << std::endl;
        // for ()
        // {

        // }
    }

    return out;
}


void PropsMX::Envelope(double P)
{
    FlshStr bolha;
    FlshStr orvalho;
    std::vector <double> z(2);
    
    std::cout << std::endl;
    std::cout << "pt" << "\t\t" << "z" << "\t\t" << "Tbolha[K]" << "\t\t" << "Torbvalho[K]"<< "\t\t" << "x1"  << "\t\t" << "y1" << std::endl;
    unsigned int npts = 21;
    for (unsigned int i = 0 ; i < npts ; i++)
    {
        z[0] = (double)i*1/(npts-1);
        z[1] = 1-z[0];

        bolha = BolhaTELV(P, z);
        orvalho = OrvalhoTELV(P, z);

        std::cout << i << "\t\t" << z[0] << "\t\t" << bolha.T << "\t\t" << orvalho.T << "\t\t" << bolha.F[0].y[0] << "\t\t" << orvalho.F[1].y[0] << std::endl;
    }
}

double PropsMX::ResRachfordRice(double beta, std::vector < double > z, std::vector < double > K)
{
    double residuo = 0;
    for (unsigned int i = 0 ; i <  z.size() ; i++) 
    {
        residuo += z[i]*(K[i] - 1)/(1 + beta*(K[i]-1));
    }
    return residuo;
}

PropsMX::FlshStr PropsMX::FlashLV(double T, double P, std::vector < double > z)
{
    FlshStr result;
    double FB = 0;
    double FBdB = 0;
    double FlB = 0;
    double dB = 0.00001;
    double e = 0.000001;
    double SumX, SumY;
    unsigned int iter = 0;
    std::vector < double > K(z.size());
    
    result.K.resize(z.size());
    result.F.resize(2);         // 2 fases 0-V e 1-L
    result.F[0].EstFis = 1;     // fase 0 no estado vapor
    result.F[1].EstFis = 0;     // fase 1 no estado liquido
    result.z = z;               // composição global da mistura
    result.T = T;
    result.P = P;

    // estimativas iniciais
    result.F[0].beta = 0.5;
    result.F[1].beta = 1 - result.F[0].beta;

    result.F[0].y = z;          // composição da fase vapor igual a dacarga como estimativa inicial
    result.F[1].y = z;          // composição da fase liquida igual a dacarga como estimativa inicial

    result.F[0].phi.resize(z.size(),1);
    result.F[1].phi = phiMxi(T,P,0,result.F[1].y);
    
    do
    {
        for (unsigned int i = 0 ; i < z.size() ; i++) 
        {
            // if (result.F[0].phi[i] == 0) K[i] = 100000;
            // else 
            K[i] = result.F[1].phi[i]/result.F[0].phi[i];
        }

        FB = ResRachfordRice(result.F[0].beta, z, K);
        FBdB = ResRachfordRice(result.F[0].beta+dB, z, K);
        FlB = (FBdB - FB)/dB;
        result.F[0].beta = result.F[0].beta - FB/FlB;

        SumX = 0;
        SumY = 0;
        for (unsigned int i = 0 ; i < z.size() ; i++)
        {
            result.F[1].y[i] = z[i]/(1 + result.F[0].beta*(K[i]-1));
            result.F[0].y[i] = K[i]*result.F[1].y[i];

            SumX += result.F[1].y[i];
            SumY += result.F[0].y[i];
        }

        result.F[0].phi = phiMxi(T,P,1,result.F[0].y);
        result.F[1].phi = phiMxi(T,P,0,result.F[1].y);
        
        iter++;
        if (iter >= maxiter) std::cout << "FlashLV atingiu maxiter " << std::endl;
        // std::cout << std::endl;
        // std::cout << "iter = " << iter <<  std::endl;
        // std::cout << "FB = " << FB << " " << abs(FB) << " " << FlB << " " << (abs(FB) > e) << " " << (iter < 50) << " " << (abs(FB) > e && (iter < 50))<< std::endl;
        // std::cout << "beta = " << result.F[0].beta << std::endl;
        // std::cout << "y[0] = " << result.F[0].y[0] << " " << result.F[1].y[0] << std::endl;
        // std::cout << "y[1] = " << result.F[0].y[1] << " " << result.F[1].y[1] << std::endl;
        // std::cout << "Sum = " << SumY << " " << SumX << std::endl;
        // std::cout << "phi[0] = " << result.F[0].phi[0] << " " << result.F[1].phi[0] << std::endl;
        // std::cout << "K[i] = " << K[0] << " " << K[1] << std::endl;
        // std::cout << std::endl;

    } while (abs(FB) > e && (iter < maxiter));

    if (result.F[0].beta > 1)
    {
        result.F[0].beta = 1;
        result.F[0].y = z;          // composição da fase vapor igual a dacarga como estimativa inicial
        for (unsigned int i = 0 ; i < result.F[1].y.size() ; i++) result.F[1].y[i] = 0;
    }
    else if (result.F[0].beta < 0)
    {
        result.F[0].beta = 0;
        for (unsigned int i = 0 ; i < result.F[0].y.size() ; i++) result.F[0].y[i] = 0;
        result.F[1].y = z;          // composição da fase líquida igual a dacarga como estimativa inicial
    }

    result.F[0].Z = ZPRTP(T,P,1,result.F[0].y);
    result.F[1].Z = ZPRTP(T,P,0,result.F[1].y);

    result.F[1].beta = 1 - result.F[0].beta;

    for (unsigned int i = 0 ; i < z.size() ; i++) 
    {
        result.K[i].resize(1);
        result.K[i][0] = K[i];
    }
    
    return result;
}

PropsMX::FlshStr PropsMX::OrvalhoTELV(double P, std::vector < double > z)
{
    FlshStr result;
    double T0 = 0;
    double FT = 0;
    double FTdT = 0;
    double FlT = 0;
    double dT = 0.01;
    double e = 0.00001;
    double Sumx = 0;
    unsigned int iter = 0;
    std::vector < double > K(z.size());
    std::vector < double > KdT(z.size());
    std::vector < double > phiVdT;
    std::vector < double > phiLdT;

    CalcTbPR(P);
    unsigned int n = z.size();
    for (unsigned int i = 0; i < z.size() ; i++)
    {
        if (isnan(TbPR[i])) n--;
        else T0 +=  TbPR[i];
        // std::cout << "TbPR[" << i << "] = " << TbPR[i] << std::endl;
    }
    T0 = T0/n;
    
    result.K.resize(z.size());
    
    result.F.resize(2);         // 2 fases 0-V e 1-L
    result.F[0].EstFis = 1;     // fase 0 no estado vapor
    result.F[1].EstFis = 0;     // fase 1 no estado liquido
    result.z = z;               // composição global da mistura
    result.P = P;
    result.F[0].beta = 1;
    result.F[1].beta = 1 - result.F[0].beta;
    result.F[0].y = z;          // composição da fase vapor igual a dacarga como estimativa inicial

    // estimativas iniciais
    result.T = T0;
    result.F[1].y = z;          // composição da fase liquida igual a dacarga como estimativa inicial
  
    result.F[0].phi = phiMxi(result.T,P,1,result.F[0].y);
    result.F[1].phi = phiMxi(result.T,P,0,result.F[1].y);

    // std::cout << std::endl << "vai entrar no DO T = " << T0 << std::endl;
   
    do
    {
        for (unsigned int i = 0 ; i < z.size() ; i++) K[i] = result.F[1].phi[i]/result.F[0].phi[i];
        FT = ResRachfordRice(result.F[0].beta, z, K);
        phiVdT = phiMxi(result.T+dT,P,1,result.F[0].y);
        phiLdT = phiMxi(result.T+dT,P,0,result.F[1].y);
        for (unsigned int i = 0 ; i < z.size() ; i++) KdT[i] = phiLdT[i]/phiVdT[i];
        FTdT = ResRachfordRice(result.F[0].beta, z, KdT);
        FlT = (FTdT - FT)/dT;
        result.T = result.T - FT/FlT;

        Sumx = 0;
        for (unsigned int i = 0 ; i < z.size() ; i++) 
        {
            result.F[1].y[i] = z[i]/(1 + result.F[0].beta*(K[i]-1));
            Sumx += result.F[1].y[i];
        }


        // std::cout << "iter = " << iter << std::endl;
        // std::cout << "FT = " << FT << " " << abs(FT) << " " << FlT <<  " " << (abs(FT) > e) <<  std::endl;
        // std::cout << "T = " << result.T << std::endl;
        // std::cout << "y[0] = " << result.F[0].y[0] << " " << result.F[1].y[0] << std::endl;
        // std::cout << "y[1] = " << result.F[0].y[1] << " " << result.F[1].y[1] << std::endl;
        // std::cout << "Sum_y = " << result.F[0].y[0]+result.F[0].y[1] << " " << Sumx << std::endl;
        // std::cout << "Zintesuss = " << ZPRTP(result.T,P,1,result.F[0].y) << " " << ZPRTP(result.T,P,0,result.F[1].y) << std::endl;
        // // MxCtsStr ctsv = CalcCtesMx(result.T, P, result.F[0].y);
        // // MxCtsStr ctsl = CalcCtesMx(result.T, P, result.F[1].y);
        // // std::cout << "ZNR = " << FZPR(1,ctsv.A,ctsv.B) << " " << FZPR(0,ctsl.A,ctsl.B) << std::endl;
        // std::cout << "phi[0] = " << result.F[0].phi[0] << " " << result.F[1].phi[0] << std::endl;
        // std::cout << "phi[1] = " << result.F[0].phi[1] << " " << result.F[1].phi[1] << std::endl;
        // std::cout << "K[i] = " << K[0] << " " << K[1] << std::endl;
        // std::cout << std::endl;

        if (isnan(result.T)) iter = maxiter;
        else
        {
            result.F[0].phi = phiMxi(result.T,P,1,result.F[0].y);
            result.F[1].phi = phiMxi(result.T,P,0,result.F[1].y);

            iter++;
        }
        

    } while ((abs(FT) > e) && (iter < maxiter));

    result.F[0].Z = ZPRTP(result.T,P,1,result.F[0].y);
    result.F[1].Z = ZPRTP(result.T,P,0,result.F[1].y);

    for (unsigned int i = 0 ; i < z.size() ; i++) 
    {
        result.K[i].resize(1);
        result.K[i][0] = K[i];
    }
    
    return result;
}

PropsMX::FlshStr PropsMX::BolhaTELV(double P, std::vector < double > z)
{
    FlshStr result;
    double T0 = 0;
    double FT = 0;
    double FTdT = 0;
    double FlT = 0;
    double dT = 0.01;
    double e = 0.00001;
    double Sumy = 0;
    unsigned int iter = 0;
    std::vector < double > K(z.size());
    std::vector < double > KdT(z.size());
    std::vector < double > phiVdT;
    std::vector < double > phiLdT;

    // CalcTbPR(P);
    for (unsigned int i = 0; i < z.size() ; i++)
    {
        T0 +=  mol[i].Tb;
        // T0 +=  TbPR[i];
    }
    T0 = 0.9*T0/z.size();

    result.K.resize(z.size());
    
    result.F.resize(2);         // 2 fases 0-V e 1-L
    result.F[0].EstFis = 1;     // fase 0 no estado vapor
    result.F[1].EstFis = 0;     // fase 1 no estado liquido
    result.z = z;               // composição global da mistura
    result.P = P;
    result.F[0].beta = 0;
    result.F[1].beta = 1 - result.F[0].beta;
    result.F[1].y = z;          // composição da fase liquida igual a dacarga como estimativa inicial
    
    // estimativas iniciais
    result.T = T0;
    result.F[0].y = z;          // composição da fase vapor igual a dacarga como estimativa inicial
    
    result.F[0].phi = phiMxi(result.T,P,1,result.F[0].y);
    result.F[1].phi = phiMxi(result.T,P,0,result.F[1].y);
    
    // std::cout << std::endl << "vai entrar no DO T = " << T0 << std::endl;
    do
    {
        for (unsigned int i = 0 ; i < z.size() ; i++) K[i] = result.F[1].phi[i]/result.F[0].phi[i];
        
        // std::cout << "K[i] = " << K[0] << " " << K[1] << std::endl;
        FT = ResRachfordRice(result.F[0].beta, z, K);
        phiVdT = phiMxi(result.T+dT,P,1,result.F[0].y);
        phiLdT = phiMxi(result.T+dT,P,0,result.F[1].y);
        for (unsigned int i = 0 ; i < z.size() ; i++) KdT[i] = phiLdT[i]/phiVdT[i];

        // std::cout << "KdT[i] = " << KdT[0] << " " << KdT[1] << std::endl;

        FTdT = ResRachfordRice(result.F[0].beta, z, KdT);
        FlT = (FTdT - FT)/dT;
        result.T = result.T - FT/FlT;

        Sumy = 0;
        for (unsigned int i = 0 ; i < z.size() ; i++) 
        {
            result.F[0].y[i] = K[i]*z[i]/(1 + result.F[0].beta*(K[i]-1));
            Sumy += result.F[0].y[i];
        }

        // std::cout << "iter = " << iter << std::endl;
        // std::cout << "FT = " << FT << " " << abs(FT) << " " << FlT <<  " " << (abs(FT) > e) <<  std::endl;
        // std::cout << "T = " << result.T << std::endl;
        // std::cout << "y[0] = " << result.F[0].y[0] << " " << result.F[1].y[0] << std::endl;
        // std::cout << "y[1] = " << result.F[0].y[1] << " " << result.F[1].y[1] << std::endl;
        // std::cout << "Sum_y = " << Sumy <<  " " << result.F[1].y[0] + result.F[1].y[1] << std::endl;
        // std::cout << "phi[0] = " << result.F[0].phi[0] << " " << result.F[1].phi[0] << std::endl;
        // std::cout << "phi[1] = " << result.F[0].phi[1] << " " << result.F[1].phi[1] << std::endl;
        // std::cout << "K[i] = " << K[0] << " " << K[1] << std::endl;
        // std::cout << std::endl;

        if (isnan(result.T)) iter = maxiter;
        else
        {
            result.F[0].phi = phiMxi(result.T,P,1,result.F[0].y);
            result.F[1].phi = phiMxi(result.T,P,0,result.F[1].y);

            iter++;
        }

    } while ((abs(FT) > e) && (iter < maxiter));

    // std::cout << "Bolha DONE" << std::endl;

    result.F[0].Z = ZPRTP(result.T,P,1,result.F[0].y);
    result.F[1].Z = ZPRTP(result.T,P,0,result.F[1].y);

    for (unsigned int i = 0 ; i < z.size() ; i++) 
    {
        result.K[i].resize(1);
        result.K[i][0] = K[i];
    }
    
    return result;
}
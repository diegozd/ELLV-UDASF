#include "auxiliar.hpp"

std::vector < double > LeComp(std::string path)
{
	std::vector < double > z; 
    std::string filename = "Composicao.csv";
	path += filename;
	std::vector < std::vector < double > > M = LeMatriz(path); 
	double sumz = 0;
	double tol =  1e-6;

	z.resize(M.size());
	for (unsigned int i = 0; i < M.size(); i++)
	{
		
		z[i] = M[i][1];
		sumz += z[i];
	}

	if (abs(1 - sumz) > tol)
	{
		for (unsigned int i = 0; i < M.size(); i++) z[i] = z[i]/sumz;
	}

	return z;
}


std::vector < std::vector < double > > LeMatriz(std::string path)
{
	std::vector < std::vector < double > > M; 
	std::ifstream file(path);

	//Msg de Erro de abertura
	//if (!file.is_open()) throw std::runtime_error("Erro: Nao foi possivel abrir o arquivo de parametros P1.csv");
	if (!file.is_open())
	{
		std::cout << "Erro: nao foi possivel abrir " << path << std::endl;
		exit(EXIT_FAILURE);
	}

	//variaveis auxiliares
	std::string line, colum;
	unsigned int c = 0, l = 0;

	// encontrando a dimens�o da matriz de par�metros
	if (file.good())
	{
		//lendo as linhas
		while (std::getline(file, line))
		{
			// se estiver na primeira lilnha
			if (l == 0)
			{
				// Create a stringstream from line
				std::stringstream ss(line);

				//conta as colunas
				while (std::getline(ss, colum, ',')) c++;
			}

			//conta as linhas 
			l++;
		}
	}

	//voltando ao come�o do arrquivo
	file.clear();
	file.seekg(0);


	M.resize(l);
	for (unsigned int i = 0; i < l; i++)
	{
		M[i].resize(c);
		//lendo a linha
		std::getline(file, line, '\n');
		std::stringstream  iss(line);

		for (unsigned int j = 0; j < c; j++)
		{
			//lendo a coluna
			std::getline(iss, colum, ',');
			M[i][j] = std::stod(colum);
		}
	}
	// Fecha arquivo
	file.close();
	return M;
}

void EscreveMatriz(std::vector < std::vector < double > > M, std::string path, std::string filename)
{
	path += filename;
	std::ofstream fout1(path);
	fout1.precision(4);
	for (unsigned int i = 0 ; i < M.size() ; i++)
	{
		for (unsigned int j = 0 ; j < M[i].size() ; j++)
		{//varre nucleo
			
			fout1 << M[i][j];
			if (j < (M[i].size() - 1)) fout1 << ",";
			// std::cout << M[i][j] << ",";
		}
		fout1 << std::endl;
		// std::cout << std::endl;
	}
}

void MostraMatrz(std::vector < std::vector < double > > M)
{
	for (unsigned int i = 0 ; i < M.size() ; i++)
	{
		for (unsigned int j = 0 ; j < M[i].size() ; j++) std::cout << M[i][j] << " ";
		std::cout << std::endl;
	}
}

void MostraVetorInt(std::vector < int >  V)
{
	for (unsigned int i = 0 ; i < V.size() ; i++) std::cout << V[i] << " ";
	std::cout << std::endl;
}

void MostraVetor(std::vector < double >  V)
{
	for (unsigned int i = 0 ; i < V.size() ; i++) std::cout << V[i] << " ";
	std::cout << std::endl;
}

void MostraTempo(double TempoDouble)
{
	double auxtempo = 0;
	double auxtempo2 = 0;
	double auxtempo3 = 0;

	TempoDouble = TempoDouble/(1e+06);

	std::cout << std::endl << "Tempo de execucao: ";
	if (TempoDouble > 3600)
	{
		auxtempo = TempoDouble/3600;
		auxtempo3 = TempoDouble - auxtempo*3600;
		auxtempo2 = auxtempo3/60;
		auxtempo3 = TempoDouble - auxtempo*3600 - auxtempo2*60;
		std::cout << auxtempo << "h " << auxtempo2 << "min " << auxtempo3 << "s" << std::endl;
	}
	else if (TempoDouble > 60)
	{
		auxtempo = TempoDouble/60;
		auxtempo2 = TempoDouble - auxtempo*60;
		std::cout << auxtempo << "min " << auxtempo2 << "s" << std::endl;
	}
	else
	{
		std::cout << TempoDouble << "s" << std::endl;
	}
}

void EscreveLog(std::string RV, std::string SOLV, double T, double P , double RSO, double ODES, double RASF, unsigned int warning)
{
	double TC = T-273.15;
	double Pkg = P * 1.01972 -1.03323;
	int Tin = TC;
	int Pin = P;
	int RSOin = RSO;

	std::string path = "./run/resultados/" + RV + "_-_" + SOLV + "-R"+ std::to_string(RSOin) + "-T" + std::to_string(Tin) + "-P" + std::to_string(Pin) + ".txt";
	std::ofstream log(path, std::ios_base::app | std::ios_base::out);

	// Dia da execução do algoritimo
	std::chrono::system_clock::time_point today = std::chrono::system_clock::now();
	std::time_t tt;
	tt = std::chrono::system_clock::to_time_t ( today );
	
	log << "\n\n";
	log << "############################################### \n";
	log << "####       RENDIMENTOS DESASFALTACAO       #### \n";
	log << "############################################### \n";
	log << "\n";
	log << ctime(&tt);
	log.precision(4);
	log << "\n";

	if (warning > 0)
	{
		log << "####    ERRO:    ####\n";
		log << "Problema na formação de duas fases.\n \n";
	}

	log << "T: " << T << " K  =>  " << TC << " C \n";
	log << "P: " << P << " bar => " << Pkg <<" kgf/cm2g\n";
	log << "\n";
	log << "OLEO:      " << RV << "\n";
	log << "SOLVENTE:  " << SOLV << "\n";
	log << "RSO:       " << RSO  << "\n";
	log << "rend ODES: " << ODES << " % \n";
	log << "rend RASF: " << RASF << " % \n";

}
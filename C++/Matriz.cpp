#include "Matriz.hpp"


template <class T> Matriz<T>::Matriz(int l,int c)
{
    M.resize(l);
    for ( int i = 0 ; i < l ; i++)
    {
        M[i].resize(c);
    }
}

template class Matriz<int>;
template class Matriz<double>;

void MatrizOperacoes::MultMatriz(std::vector < std::vector <double> > &M1 , std::vector < std::vector <double> > &M2)
{
    if (M1[0].size() != M2.size())
    {
        std::cout << "Matrizes nao multiplicaveis" << std::endl;
    }
    else
    {
		R.resize(M1.size());
        for (unsigned int i = 0 ; i < M1.size() ; i++)
        {
            R[i].resize(M2[0].size());
        }
    }

    for (unsigned int i = 0 ; i < R.size() ; i++)
    {
        for(unsigned int j = 0 ; j < R[0].size() ; j++)
        {
            R[i][j] = 0;
            for(unsigned int k = 0 ; k < M2.size();k++)
            {
                R[i][j] += M1[i][k]*M2[k][j];
            }
        }
    }
}

void MatrizOperacoes::MatrizTransposta(std::vector < std::vector <double> > &M1)
{
	T.resize(M1[0].size());
	for (unsigned int i = 0 ; i < M1[0].size() ; i++)
	{
		T[i].resize(M1.size());
	}

	for(unsigned int i =0 ; i < M1.size() ; i++)
	{
		for(unsigned int j = 0; j < M1[0].size() ; j++)
		{
			T[j][i] = M1[i][j];
		}
	}
}


#include <iostream>
#include <cstdlib> 
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <ctime>
#include <cstdio>

#define HSIZE 5
#define VSIZE 5

#define RANK 1
/*
 *	Aleksander Balicki - 220989;
 * 	Dominika Rogozińska - 221094
 *	Pracownia 2.20
 *	Porównanie eliminacji Gaussa: bez wyboru, z wyborem wierszy, kolumn, pełnym
 *	kompilujemy:
 *		g++ program.cpp
 */

using namespace std;
typedef vector<double> row;
typedef vector< vector<double> > matrix;

matrix gen_random(int cols, int rows, int rank)
{
	matrix m(rows);
	if(rank > min(cols, rows)|| rank == 0)
		fprintf(stderr, "Rank wrong\n");

	int lin_dep = min(cols, rows)- rank;

	for(int rz = 0; rz < rows; ++rz)
	{
		row r0(cols);
		m[rz] = r0;
	}

	for(int cz = 0; cz < cols; ++cz)
		m[0][cz] = (double(rand())/ RAND_MAX);

	for(int r = 1; r < rows; r++)
	{
		int scalar;
		do
			scalar = 10 * double(rand())/ RAND_MAX;
		while(scalar == 0);

		for(int c = 0; c < cols; ++c)
		{
			if(r == c && lin_dep <= 0)
			{
				double old;
				do
				{
					old = m[r][c];
					m[r][c]+= double(rand())/ RAND_MAX; 
				}while(m[r][c] == old);
			}
			else if(r == c)
			{
				m[r][c] = m[0][c];
				lin_dep--;
			}
			else
				m[r][c] = m[0][c];
			m[r][c] *= scalar;
		}
	}
	return m;	
}



void print(matrix m)
{	
	for(int i = 0; i < m[0].size(); ++i)
	{
		for(int j = 0; j < m.size(); ++j)
			cout << " " << m[i][j];
		cout << endl;
	}
}

void random_1_1000(matrix &m)
{
	for(int i = 0; i < m[0].size(); ++i)
		for(int j = 0; j < m.size(); ++j)
			m[j][i] = -1. + ((double) rand() / RAND_MAX)*2;

	int h = ((double) rand()/RAND_MAX)*HSIZE;
	int v = ((double) rand()/RAND_MAX)*VSIZE;
	m[h][v] += -999 + ((double) rand() / RAND_MAX)*1998
	h = ((double) rand()/RAND_MAX)*HSIZE;
	v = ((double) rand()/RAND_MAX)*VSIZE;
	m[h][v] += -999 + ((double) rand() / RAND_MAX)*1998
	h = ((double) rand()/RAND_MAX)*HSIZE;
	v = ((double) rand()/RAND_MAX)*VSIZE;
	m[h][v] += -999 + ((double) rand() / RAND_MAX)*1998
}

void hilbert(matrix &m)
{
	for(int i = 0; i < m[0].size(); ++i)
		for(int j = 0; j < m.size(); ++j)
			m[j][i] = 1./(i + j + 1);
}

void pei(matrix &m, double d)
{
	for(int i = 0; i < m[0].size(); ++i)
		for(int j = 0; j < m.size(); ++j)
			if(i == j)
				m[j][i] = d;
			else
				m[j][i] = 1;
}

int compute_nonempty(matrix m)
{
	int nonempty = 0;

	for(int i = 0; i < m[0].size(); ++i)
	{
		int current = 0;
		for(int j=0; j < m.size(); ++j)
			if(m[i][j]!= 0)
				current = 1;

		nonempty += current;
	}
	return nonempty;
}

int find_max_in_col(matrix m, int k)
{
	double max = m[k][k];
	int max_row = k;
	for(int i = k+1; i < m.size(); ++i)
	{
		if(fabs(m[i][k])> fabs(max))
		{
			max = m[i][k];
			max_row = i;
		}
	}
	return max_row;
}

vector<int> find_max_in_subm(matrix m, int k)
{
	double max = m[k][k];
	vector<int> row_col(2);
	row_col[0] = k;
	row_col[1] = k;
	for(int i = k+1; i < m.size(); ++i)
	{
		for(int j = k+1; j < m.size(); ++j)
		{
			if(fabs(m[i][j])> fabs(max))
			{
				max = m[i][j];
				row_col[0] = i;
				row_col[1] = j;
			}
		}
	}
	return row_col;
}

void switch_rows(matrix &m , int r1, int r2)
{
	double temp;
	for(int i = 0; i < m[r1].size(); i++)
	{
		temp = m[r1][i];
		m[r1][i] = m[r2][i];
		m[r2][i] = temp;
	}
}

void switch_cols(matrix &m , int c1, int c2)
{
	double temp;
	for(int i = 0; i < m.size(); i++)
	{
		temp = m[i][c1];
		m[i][c1] = m[i][c2];
		m[i][c2] = temp;
	}
}

void subtract_row(matrix &m, int from, int what, double times)
{
	for(int i = 0; i < m[0].size(); i++)
		m[from][i] -= times * m[what][i];
}


int gauss_rank_col(matrix m)
{
	cout << "Gaussian with choice from column" << endl;
	cout << "Before:" << endl;
	print(m);
	cout << endl;
	for(int k = 0; k < m[0].size()- 1; ++k)
	{
		int max_row = find_max_in_col(m, k);
		switch_rows(m, k, max_row);
		if(m[k][k] != 0)
			for(int i = k + 1; i < m.size(); i++)
			{
				subtract_row(m, i, k, m[i][k]/m[k][k]);
				m[i][k] = 0.0;
			}
	}
	cout << "After:" << endl;
	print(m);
	cout << endl;
	return compute_nonempty(m);
}

int gauss(matrix m)
{
	cout << "Gaussian without choices" << endl;
	cout << "Before:" << endl;
	print(m);
	cout << endl;
	for(int k = 0; k < m.size()- 1; ++k)
	{
		if(m[k][k] != 0)
			for(int i = k + 1; i < m.size(); i++)
			{
				subtract_row(m, i, k, m[i][k]/m[k][k]);
				m[i][k] = 0.0;
			}
	}
	cout << "After:" << endl;
	print(m);
	cout << endl;
	return compute_nonempty(m);
}

int gauss_rank_full(matrix m)
{
	cout << "Gaussian with full choice" << endl;
	cout << "Before:" << endl;
	print(m);
	for(int k = 0; k < m[0].size()- 1; ++k)
	{
		vector<int> row_col = find_max_in_subm(m, k);
		switch_rows(m, row_col.front(), k);
		switch_cols(m, row_col.back(), k);
		if(m[k][k]!= 0)
			for(int i = k + 1; i < m.size(); i++)
			{
				subtract_row(m, i, k, m[i][k]/ m[k][k]);
				m[i][k] = 0;
			}
	}
	cout << "After:" << endl;
	print(m);
	return compute_nonempty(m);
}


void printrow(row r)
{
	for(int i = 0; i < r.size(); ++ i)
		cout << " " << r[i];
}


int main(int argc,char **argv)
{    
	int rank;
	srand(time(NULL));
	
	matrix m = gen_random(HSIZE, VSIZE, RANK);
	
	cout << endl;
	rank = gauss_rank_col(m);
	cout << endl;

	cout << endl;
	rank = gauss_rank_full(m);
	cout << endl;
	
	cout << endl;
	rank = gauss(m);
	cout << endl;
}

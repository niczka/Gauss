#include <iostream>
#include <cstdlib> 
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <ctime>
#include <cstdio>

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
typedef vector<row> matrix;

void print(const row &r)
{
	for(int i = 0; i < r.size(); ++ i)
		cout << " " << r[i];
	cout << endl;
}

void print(const matrix &m)
{	
	for(int i = 0; i < m[0].size(); ++i)
		print(m[i]);
}

matrix new_matrix(int n)
{
	matrix m(n);
	for (int i = 0; i < n; i++)
		m[ i ] = row(n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			m[ i ][ j ] = 0.0;
	
	return m;
}

void random_1_1000(matrix &m)
{
	for(int i = 0; i < m[0].size(); ++i)
		for(int j = 0; j < m.size(); ++j)
			m[j][i] = -1. + ((double) rand() / RAND_MAX)*2;

	int h = ((double) rand()/RAND_MAX)*m[0].size();
	int v = ((double) rand()/RAND_MAX)*m.size();
	m[h][v] += -999 + ((double) rand() / RAND_MAX)*1998;
	h = ((double) rand()/RAND_MAX)*m[0].size();
	v = ((double) rand()/RAND_MAX)*m.size();
	m[h][v] += -999 + ((double) rand() / RAND_MAX)*1998;
	h = ((double) rand()/RAND_MAX)*m[0].size();
	v = ((double) rand()/RAND_MAX)*m.size();
	m[h][v] += -999 + ((double) rand() / RAND_MAX)*1998;
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
			m[j][i] = (i == j) ? d : 1;
}

void dominating(matrix &m)
{
	for(int i = 0; i < m[0].size(); ++i)
		for(int j = 0; j < m.size(); ++j)
			m[j][i] = (i == j) ? 0 : -50 +  ((double)rand() / RAND_MAX) * 100;

	for(int i = 0; i < m[0].size(); ++i)
	{
		int sum = 0;
		for(int j = 0; j < m.size(); ++j)
			sum += fabs(m[i][j]);
		m[i][i] = sum + 1;
	}
}

int find_max_in_col(matrix &m, int k)
{
	double max = m[k][k];
	int max_row = k;
	for(int i = k+1; i < m.size(); ++i)
	{
		if(fabs(m[i][k]) > fabs(max))
		{
			max = m[i][k];
			max_row = i;
		}
	}
	return max_row;
}

int find_max_in_row(matrix &m, int k)
{
	double max = m[k][k];
	int max_col = k;
	for(int i = k+1; i < m.size(); ++i)
	{
		if(fabs(m[k][i]) > fabs(max))
		{
			max = m[k][i];
			max_col = i;
		}
	}
	return max_col;
}

vector<int> find_max_in_subm(matrix &m, int k)
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

void switch_rows(matrix &m, int r1, int r2)
{
	double temp;
	for(int i = 0; i < m[r1].size(); i++)
	{
		temp = m[r1][i];
		m[r1][i] = m[r2][i];
		m[r2][i] = temp;
	}
}

void switch_cols(matrix &m, int c1, int c2)
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


matrix gauss_rank_col(matrix m)
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
	return m;
}

matrix gauss_rank_row(matrix m)
{
	cout << "Gaussian with choice from row" << endl;
	cout << "Before:" << endl;
	print(m);
	cout << endl;
	for(int k = 0; k < m[0].size()- 1; ++k)
	{
		int max_col = find_max_in_row(m, k);
		switch_cols(m, k, max_col);
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
	return m;
}

matrix gauss(matrix m)
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
	return m;
}

matrix gauss_rank_full(matrix m)
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
				m[i][k] = 0.0;
			}
	}
	cout << "After:" << endl;
	print(m);
	return m;
}

int main(int argc,char **argv)
{    
	int rank;
	srand(time(NULL));
	
	matrix m, mg, mgf, mgc;
	m = new_matrix(5);
	hilbert(m);
	
	cout << endl;
	mgc = gauss_rank_col(m);
	cout << endl;

	cout << endl;
	mgf = gauss_rank_full(m);
	cout << endl;
	
	cout << endl;
	mg = gauss(m);
	cout << endl;
	
	cout << endl;
	mg = gauss_rank_row(m);
	cout << endl;
}

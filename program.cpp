#include <iostream>
#include <cstdlib> 
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <ctime>
#include <cstdio>
#include <iomanip>
/*
 *	Aleksander Balicki - 220989;
 * 	Dominika Rogozińska - 221094
 *	Pracownia 2.20
 *	Porównanie eliminacji Gaussa: bez wyboru, z wyborem wierszy, kolumn, pełnym
 *	kompilujemy:
 *		g++ program.cpp >> out.csv
 */

using namespace std;
typedef vector<double> row;
typedef vector<row> matrix;
typedef vector<int> vars;
typedef vector<double> result;

#define PRECISION 20

typedef struct
{
	vars v;
	result r;
} solution;

typedef struct
{
	vars v;
	matrix m;
	row b;
} ematrix;

void print(const vars &v)
{
	for(int i = 0; i < v.size(); ++ i)
		cout << " x" << v[i];
	cout << endl;
}

void print(const row &r)
{
	for(int i = 0; i < r.size(); ++ i)
		cout << "\t" << fixed << setprecision (9) << r[i];
}

void print(const matrix &m, const row &b)
{	
	for(int i = 0; i < m[0].size(); ++i)
		cout << "\t" << "x" << i;
	
	cout << endl;

	for(int i = 0; i < m[0].size(); ++i)
	{
		print(m[i]);
		cout << "\t=\t" << fixed << setprecision (9) << b[i] << endl;
	}
}

void print(const ematrix &e)
{
	for(int i = 0; i < e.v.size(); ++i)
		cout << "\t" << "x" << e.v[i];
	cout << endl;

	for(int i = 0; i < e.m[0].size(); ++i)
	{
		print(e.m[i]); 
		cout << "\t=\t" << fixed << setprecision (9) << e.b[i] << endl;
	}
}

void print(const solution &s)
{
	for(int i = 0; i < s.v.size(); i++)
		cout << "\tx" << s.v[i] << " = " << fixed << setprecision (30) << s.r[i] << endl;
}

ematrix new_ematrix(matrix m, vars v, row b)
{
	ematrix e;
	e.m = m;
	e.v = v;
	e.b = b;
	return e;
}

solution new_solution(vars v, result r)
{
	solution s;
	s.v = v;
	s.r = r;
	return s;
}

vars new_vars(int n)
{
	vars v(n);
	for(int i = 0; i < n; i++)
		v[i]=i;
	return v;
}


void random_row(row &r)
{
	for (int i = 0; i < r.size(); i++)
		r[i] = -500 +  ((double)rand() / RAND_MAX) * 1000;
}

row new_random_row(int n)
{
	row r(n);
	random_row(r);
	return r;
}

matrix new_matrix(int n)
{
	matrix m(n);
	for (int i = 0; i < n; i++)
	{
		m[ i ] = row(n);
		for (int j = 0; j < n; j++)
			m[ i ][ j ] = 0.0;
	}
	
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

void switch_rows(matrix &m, int r1, int r2, row &b)
{
	double temp;
	for(int i = 0; i < m[r1].size(); i++)
	{
		temp = m[r1][i];
		m[r1][i] = m[r2][i];
		m[r2][i] = temp;
	}
	temp = b[r1];
	b[r1] = b[r2];
	b[r2] = temp;
}

void switch_cols(matrix &m, int c1, int c2, vars &v)
{
	double temp;
	for(int i = 0; i < m.size(); i++)
	{
		temp = m[i][c1];
		m[i][c1] = m[i][c2];
		m[i][c2] = temp;
	}

	int itemp = v[c1];
	v[c1] = v[c2];
	v[c2] = itemp;
}

void subtract_row(matrix &m, int from, int what, double times, row &b)
{
	for(int i = 0; i < m[0].size(); i++)
		m[from][i] -= times * m[what][i];

	b[from] -= times * b[what];
}


ematrix gauss_rank_col(matrix m, row b)
{
	vars v = new_vars(m.size());
	for(int k = 0; k < m[0].size()- 1; ++k)
	{
		int max_row = find_max_in_col(m, k);
		switch_rows(m, k, max_row, b);
		
		if(m[k][k] != 0)
			for(int i = k + 1; i < m.size(); i++)
			{
				subtract_row(m, i, k, m[i][k]/m[k][k], b);
				m[i][k] = 0.0;
			}
	}
	return new_ematrix(m, v, b);
}

ematrix gauss_rank_row(matrix m, row b)
{
	vars v = new_vars(m.size());
	for(int k = 0; k < m[0].size()- 1; ++k)
	{
		int max_col = find_max_in_row(m, k);
		switch_cols(m, k, max_col, v);

		if(m[k][k] != 0)
			for(int i = k + 1; i < m.size(); i++)
			{
				subtract_row(m, i, k, m[i][k]/m[k][k], b);
				m[i][k] = 0.0;
			}
	}
	return new_ematrix(m, v, b);
}

ematrix gauss(matrix m, row b)
{
	vars v = new_vars(m.size());
	for(int k = 0; k < m.size()- 1; ++k)
	{
		if(m[k][k] != 0)
			for(int i = k + 1; i < m.size(); i++)
			{
				subtract_row(m, i, k, m[i][k]/m[k][k], b);
				m[i][k] = 0.0;
			}
	}
	return new_ematrix(m, v, b);
}

ematrix gauss_rank_full(matrix m, row b)
{
	vars v = new_vars(m.size());
	for(int k = 0; k < m[0].size()- 1; ++k)
	{
		vector<int> row_col = find_max_in_subm(m, k);
		switch_rows(m, row_col.front(), k, b);
		switch_cols(m, row_col.back(), k, v);


		if(m[k][k]!= 0)
			for(int i = k + 1; i < m.size(); i++)
			{
				subtract_row(m, i, k, m[i][k]/ m[k][k], b);
				m[i][k] = 0.0;
			}
	}
	return new_ematrix(m, v, b);
}

solution back_substitution(const ematrix &e)
{
	matrix m = e.m;
	result r(m.size());
	for(int i = (m.size()) - 1; i >= 0; i--)
	{
    		r[i] = e.b[i];
    		for (int j = i + 1; j < m.size(); j++)
      			r[i] -= m[i][j] * r[j];
    		r[i] /= m[i][i];
  	}
	return new_solution(e.v, r);
}

solution sort(const solution &s)
{
	solution s1;
	s1.v = new_vars(s.v.size());
	row r(s.v.size());
	s1.r = r;

	for(int i = 0; i < s.v.size(); i++)
		s1.r[s.v[i]] = s.r[i];
	return s1;
}

result mult(matrix m, row b)
{
	result r(b.size());
	for (int i = 0; i < r.size(); i++)
	{
		r[i] = 0;
		for (int j = 0; j < b.size(); j++)
			r[i] += m[i][j] * b[j];
	}
	return r;
}

row subtract(result a, result b)
{
	row r(a.size());
	for(int i = 0; i < r.size(); i++)
		r[i] = a[i] - b[i];
	return r;
}

double norm_row_inf(row x)
{
	double max = fabs(x[0]);
	for(int i = 1; i < x.size(); i++)
		if (fabs(x[i]) > max)
			max = fabs(x[i]);
	return max;
}

double eq(matrix m, row b, solution x)
{
	return norm_row_inf(subtract(mult(m, x.r), b));
}



int main(int argc,char **argv)
{    
	int rank;
	srand(time(NULL));

	solution s;
	matrix m;
	ematrix mg, mgf, mgc, mgr;
	row b;
	
	cout << "Size | Gaussian | Gaussian with choice from column | Gaussian with full choice | Gaussian with choice from row" << endl;
	
	cout << "Hilbert matrix" << endl;
	for (int i = 3; i <= 100; i++)
	{
		m = new_matrix(i);
		hilbert(m);
		
		b = new_random_row(i);
		
		mg = gauss(m, b);
		s = sort(back_substitution(mg));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << " ;";
		
		mgc = gauss_rank_col(m, b);
		s = sort(back_substitution(mgc));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << " ;";

		mgf = gauss_rank_full(m, b);
		s = sort(back_substitution(mgf));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << " ;";
		
		mgr = gauss_rank_row(m, b);
		s = sort(back_substitution(mgr));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << endl;
	}



	cout << "PEI matrix" << endl;
	for (int i = 3; i <= 100; i++)
	{
		m = new_matrix(i);
		pei(m, 5);
		
		b = new_random_row(i);		
		
		mg = gauss(m, b);
		s = sort(back_substitution(mg));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << " ;";
		
		mgc = gauss_rank_col(m, b);
		s = sort(back_substitution(mgc));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << " ;";

		mgf = gauss_rank_full(m, b);
		s = sort(back_substitution(mgf));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << " ;";
		
		mgr = gauss_rank_row(m, b);
		s = sort(back_substitution(mgr));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << endl;
	}	
	
	cout << "Random matrix" << endl;
	for (int i = 3; i <= 100; i++)
	{
		m = new_matrix(i);
		random_1_1000(m);
		
		b = new_random_row(i);
		
		mg = gauss(m, b);
		s = sort(back_substitution(mg));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << " ;";
		
		mgc = gauss_rank_col(m, b);
		s = sort(back_substitution(mgc));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << " ;";

		mgf = gauss_rank_full(m, b);
		s = sort(back_substitution(mgf));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << " ;";
		
		mgr = gauss_rank_row(m, b);
		s = sort(back_substitution(mgr));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << endl;
	}	
	
	
	cout << "Dominating diagonal matrix" << endl;
	for (int i = 3; i <= 100; i++)
	{
		m = new_matrix(i);
		dominating(m);
		
		b = new_random_row(i);
		
		mg = gauss(m, b);
		s = sort(back_substitution(mg));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << " ;";
		
		mgc = gauss_rank_col(m, b);
		s = sort(back_substitution(mgc));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << " ;";

		mgf = gauss_rank_full(m, b);
		s = sort(back_substitution(mgf));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << " ;";
		
		mgr = gauss_rank_row(m, b);
		s = sort(back_substitution(mgr));
		cout << i << " " << eq(m, b, s) << fixed << setprecision (PRECISION) << endl;
	}	
	/*
	solution s1;
	matrix m1;
	ematrix mg1, mgf1, mgc1, mgr1;
	m1 = new_matrix(5);
	m1[0][0] = 2;
	m1[0][1] = 1;
	m1[0][2] = 1;
	m1[0][3] = 2;
	m1[0][4] = 0;
	m1[1][0] = 1;
	m1[1][1] = 1;
	m1[1][2] = 1;
	m1[1][3] = 1;
	m1[1][4] = 1;
	m1[2][0] = 2;
	m1[2][1] = 0;
	m1[2][2] = 3;
	m1[2][3] = 1;
	m1[2][4] = 0;
	m1[3][0] = 0;
	m1[3][1] = 1;
	m1[3][2] = 0;
	m1[3][3] = 5;
	m1[3][4] = 0;
	m1[4][0] = 1;
	m1[4][1] = 0;
	m1[4][2] = 3;
	m1[4][3] = 1;
	m1[4][4] = 4;
	
	row b1(5);
	b1[0] = 21;
	b1[1] = 15;
	b1[2] = 21;
	b1[3] = 14;
	b1[4] = 20;

	hilbert(m1);

	cout << endl;
	cout << "Gaussian with choice from column" << endl;
	cout << "Before:" << endl;
	print(m1,b1);
	cout << endl;
	cout << "After:" << endl;
	mgc1 = gauss_rank_col(m1, b1);
	print(mgc1);
	cout << "Variable order:" << endl;
	print(mgc1.v);
	cout << "Solution:" << endl;
	s1 = back_substitution(mgc1);
	print(s1);
	cout << "Sorted solution" << endl;
	s1 = sort(back_substitution(mgc1));
	print(s1);
	cout << eq(m1, b1, s1);

	cout << endl;
	cout << "Gaussian with full choice" << endl;
	cout << "Before:" << endl;
	print(m1,b1);
	cout << endl;
	mgf1 = gauss_rank_full(m1, b1);
	cout << "After:" << endl;
	print(mgf1);
	cout << "Variable order:" << endl;
	print(mgf1.v);
	cout << "Solution:" << endl;
	s1 = back_substitution(mgf1);
	print(s1);
	cout << "Sorted solution" << endl;
	s1 = sort(back_substitution(mgf1));
	print(s1);
	cout << eq(m1, b1, s1);
	
	cout << endl;
	cout << "Gaussian without choices" << endl;
	cout << "Before:" << endl;
	print(m1,b1);
	cout << endl;
	mg1 = gauss(m1, b1);
	cout << "After:" << endl;
	print(mg1);
	cout << "Variable order:" << endl;
	print(mg1.v);
	cout << "Solution:" << endl;
	s1 = back_substitution(mg1);
	print(s1);
	cout << "Sorted solution" << endl;
	s1 = sort(back_substitution(mg1));
	print(s1);
	cout << eq(m1, b1, s1);

	
	cout << endl;
	cout << "Gaussian with choice from row" << endl;
	cout << "Before:" << endl;
	print(m1,b1);
	cout << endl;
	mgr1 = gauss_rank_row(m1, b1);
	cout << "After:" << endl;
	print(mgr1);
	cout << "Variable order:" << endl;
	print(mgr1.v);
	cout << "Solution:" << endl;
	s1 = back_substitution(mgr1);
	print(s1);
	cout << "Sorted solution" << endl;
	s1 = sort(back_substitution(mgr1));
	print(s1);
	cout << eq(m1, b1, s1);*/

}

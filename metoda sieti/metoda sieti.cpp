// metoda sieti.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <chrono>

using namespace std;

//napisane vo Visual Studio 2015 - (kompilátor: Visual C++ 14.0)

vector<double> TDMA(const int N, vector<vector<double>> matrix, vector<double> vec) //pre tridiagonalne matice (Thomas algorithm)
{
	int i;
	vector<double> x;
	double help;
	x.resize(N);

	for (i = 1; i < N; i++)
	{
		help = matrix[i][i - 1] / matrix[i - 1][i - 1];
		matrix[i][i - 1] -= matrix[i - 1][i - 1] * help;
		matrix[i][i] -= matrix[i - 1][i] * help;;
		vec[i] -= vec[i - 1] * help;
	}

	x[N - 1] = vec[N - 1] / matrix[N - 1][N - 1];

	for (i = N - 2; i > -1; i--)
		x[i] = (vec[i] - x[i + 1] * matrix[i][i + 1]) / matrix[i][i];

	return x;
}

int main()
{
	vector<vector<double>> M_matrix, riesenie;
	vector<double> u_old, u_new, b_vector, temp;
	int i, j, k, m;
	double h, mi;

	//konštanty
	const double pi = 3.14159265358979323;
	const int n = 100;
	const double c = 4.99*(10E-6);
	double tau;
	int time;

	//nastavenie rozmeru poli
	riesenie.resize(3);
	u_new.resize(n);
	u_old.resize(n);
	M_matrix.resize(n);
	temp.resize(n);
	b_vector.resize(n);
	for (i = 0; i < 3; i++) riesenie[i].resize(n + 2);
	for (i = 0; i < n; i++) M_matrix[i].resize(n);

	//nastavenie presnosti vypisu
	cout.precision(10);
	cout.setf(std::ios::fixed, std::ios::floatfield);

	//vypocet h
	h = 1.0 / double(n + 1);

	//zapis deliacich bodov do pola rieseni
	for (j = 0; j < n + 2; j++)
		riesenie[0][j] = j*h;

#pragma region tau=0.1 explicitná schéma

	//inicializacia merania casu, vyuziva chrono kniznicu
	auto zaciatok_explicit = chrono::high_resolution_clock::now();

	tau = 0.1;
	time = 10000;

	//vypocet mi
	mi = (c*tau) / pow(h, 2);

	if (tau >(pow(h, 2) / (2 * c))) cout << "nie je splnena nutna podmienka pre stabilitu riesenia (tau <= h*h/2*c)" << endl;

	//naplnenie matice M
	for (i = 0; i < n; i++)
		M_matrix[i][i] = 1 - mi * 2;

	for (i = 0; i < n - 1; i++)
		M_matrix[i][i + 1] = mi * 1;

	for (i = 1; i < n; i++)
		M_matrix[i][i - 1] = mi* 1;

	//okrajove podmienky
	b_vector[0] = 0 * mi;
	b_vector[n - 1] = 1 * mi;

	//pociatocne podmienky - inicializacia u_old
	
	for (i = 0; i < n; i++)
		u_old[i] = sin(2.0 * pi*riesenie[0][i + 1]) + riesenie[0][i + 1];

	//kontrolne vypisy

	cout << "tau \t" << tau << endl;
	cout << "c \t" << c << endl;
	cout << "h \t" << h << endl;
	cout << "mi \t" << mi << endl;
	cout << "h*h/2*c \t" << pow(h, 2) / (2 * c) << endl;

	//hlavny cyklus
	for (i = 0; i < time; i++)
	{
		for (k = 0; k < n; k++)
		{
			u_new[k] = 0;
			if (k == 0)
				for (m = k; m < k + 2; m++)
					u_new[k] += M_matrix[k][m] * u_old[m];
			else if (k == n - 1)
				for (m = k - 1; m < k + 1; m++)
					u_new[k] += M_matrix[k][m] * u_old[m];
			else
				for (m = k - 1; m < k + 2; m++)
					u_new[k] += M_matrix[k][m] * u_old[m];
		}

		u_new[0] += b_vector[0];
		u_new[n - 1] += b_vector[n - 1];

		for (j = 0; j < n; j++) u_old[j] = u_new[j];
	}

	for (j = 0; j < n; j++)
		riesenie[2][j + 1] = u_new[j];

	auto koniec_explicit = chrono::high_resolution_clock::now();
	cout << "explicitna schema (tau = 0.1) - vypoctovy cas: " << chrono::duration_cast<chrono::nanoseconds>(koniec_explicit - zaciatok_explicit).count()*10E-10 << "s" << std::endl;

	cout << endl;

#pragma endregion

#pragma region tau=1 explicitná schéma - pre tau = 1 je nestabilná

	//inicializacia merania casu, vyuziva chrono kniznicu
	auto zaciatok_implicit = chrono::high_resolution_clock::now();

	tau = 1;
	time = 1000;

	//vypocet mi
	mi = (c*tau) / pow(h, 2);

	if (tau >(pow(h, 2) / (2 * c))) cout << "nie je splnena nutna podmienka pre stabilitu riesenia (tau <= h*h/2*c)" << endl;

	//naplnenie matice M
	for (i = 0; i < n; i++)
		M_matrix[i][i] = 1 - mi * 2;

	for (i = 0; i < n - 1; i++)
		M_matrix[i][i + 1] = mi * 1;

	for (i = 1; i < n; i++)
		M_matrix[i][i - 1] = mi * 1;

	//okrajove podmienky
	b_vector[0] = 0 * mi;
	b_vector[n - 1] = 1 * mi;

	//pociatocne podmienky - inicializacia u_old

	for (i = 0; i < n; i++)
		u_old[i] = sin(2.0 * pi*riesenie[0][i + 1]) + riesenie[0][i + 1];

	//kontrolne vypisy

	cout << "tau \t" << tau << endl;
	cout << "c \t" << c << endl;
	cout << "h \t" << h << endl;
	cout << "mi \t" << mi << endl;
	cout << "h*h/2*c \t" << pow(h, 2) / (2 * c) << endl;

	//hlavny cyklus
	for (i = 0; i < time; i++)
	{
		for (k = 0; k < n; k++)
		{
			u_new[k] = 0;
			if (k == 0)
				for (m = k; m < k + 2; m++)
					u_new[k] += M_matrix[k][m] * u_old[m];
			else if (k == n - 1)
				for (m = k - 1; m < k + 1; m++)
					u_new[k] += M_matrix[k][m] * u_old[m];
			else
				for (m = k - 1; m < k + 2; m++)
					u_new[k] += M_matrix[k][m] * u_old[m];
		}

		u_new[0] += b_vector[0];
		u_new[n - 1] += b_vector[n - 1];

		for (j = 0; j < n; j++) u_old[j] = u_new[j];
	}

	for (j = 0; j < n; j++)
		riesenie[1][j + 1] = u_new[j];

	auto koniec_implicit = chrono::high_resolution_clock::now();
	cout << "implicitna schema (tau = 1) - vypoctovy cas: " << chrono::duration_cast<chrono::nanoseconds>(koniec_implicit - zaciatok_implicit).count()*10E-10 << "s" << std::endl;

#pragma endregion

#pragma region tau=1 implicitná schéma - stabilná vždy!!!!!

	////inicializacia merania casu, vyuziva chrono kniznicu
	//auto zaciatok_implicit = chrono::high_resolution_clock::now();

	//tau = 1;
	//time = 1000;

	////vypocet mi
	//mi = (c*tau) / pow(h, 2);

	////naplnenie matice M
	//for (i = 0; i < n; i++)
	//	M_matrix[i][i] = 1 + mi * 2;

	//for (i = 0; i < n - 1; i++)
	//	M_matrix[i][i + 1] = -mi * 1;

	//for (i = 1; i < n; i++)
	//	M_matrix[i][i - 1] = -mi * 1;

	////okrajove podmienky
	//b_vector[0] = 0 * mi;
	//b_vector[n - 1] = 1 * mi;

	////pociatocne podmienky - inicializacia u_old

	//for (i = 0; i < n; i++)
	//	u_old[i] = sin(2.0 * pi*riesenie[0][i + 1]) + riesenie[0][i + 1];

	////kontrolne vypisy

	//cout << "tau \t" << tau << endl;
	//cout << "c \t" << c << endl;
	//cout << "h \t" << h << endl;
	//cout << "mi \t" << mi << endl;
	//cout << "h*h/2*c \t" << pow(h, 2) / (2 * c) << endl;

	////hlavny cyklus
	//for (i = 0; i < time; i++)
	//{
	//	u_old[0] += b_vector[0];
	//	u_old[n - 1] += b_vector[n - 1];
	//	u_old.swap(TDMA(n, M_matrix, u_old));
	//}

	//for (j = 0; j < n; j++)
	//	riesenie[1][j + 1] = u_old[j];

	//auto koniec_implicit = chrono::high_resolution_clock::now();
	//cout << "implicitna schema (tau = 1) - vypoctovy cas: " << chrono::duration_cast<chrono::nanoseconds>(koniec_implicit - zaciatok_implicit).count()*10E-10 << "s" << std::endl;

#pragma endregion

	riesenie[1][0] = b_vector[0];
	riesenie[1][n + 1] = 1;

	riesenie[2][0] = b_vector[0];
	riesenie[2][n + 1] = 1;

	ofstream outfile1; //zápis do suboru csv na odovzdanie ulohy
	outfile1.open("priloha2.csv");

	//nastavenie presnosti vypisu
	outfile1.precision(12);
	outfile1.setf(std::ios::fixed, std::ios::floatfield);

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < n + 2; j++)
		{
			if (j < n + 1) outfile1 << riesenie[i][j] << ",";
			else outfile1 << riesenie[i][j];
		}
		outfile1 << "\n";
	}
	outfile1.close();

	system("PAUSE");

    return 0;
}


#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <climits>
#include <iostream>
#include <vector>
#include <unordered_set>
#include "myrand.h"

#define N 10000 // population size
#define nWeek 4 // number of week
#define nTeam 6 // number of team
#define nField 2 // number of field
#define nTimeSlot nWeek*10 // number of timeslot
#define l nTimeSlot*nField

using namespace std;
void getIndexFromMatchId(int*, int*, int);
void exitWithError(int);
int evaluate(int*);
int eval(int*);

MyRand myrand;

// each team has a match with every other team
int nMatch = nTeam * (nTeam - 1) / 2;
vector<int> matches[nTeam];

// date preference [Team][Date]
int preferenceArray[10][10] =
{
	{ -1, -2, 1, -1, 1, 0, 0, 0, 2, 0, },
	{ -2, -1, 2, 0, -1, 0, 0, 1, 1, 0, },
	{ 0, 1, -2, -1, -1, 0, 0, 1, 2, 0, },
	{ 2, 1, 0, 0, -1, 1, 0, 0, -2, -1, },
	{ -1, 2, 0, 0, -2, 1, -1, 0, 0, 1, },
	{ 0, -2, 0, -1, 0, 2, 0, 1, -1, 1, },
	{ 0, 1, 2, 0, -2, 0, 1, -1, 0, -1, },
	{ 0, 0, -2, 0, 1, 1, -1, 0, -1, 2, },
	{ 2, -1, 0, 0, 0, 1, -1, -2, 0, 1, },
	{ -2, 0, 1, 0, 0, 2, -1, 1, 0, -1, },
};
int main()
{

	int *chromosomes[N];
	int *chromos[N];

	int mid = 0;
	for (int i = 0; i < nTeam; ++i)
	{
		for (int j = i + 1; j < nTeam; ++j)
		{
			matches[i].push_back(mid);
			matches[j].push_back(mid);
			mid++;
		}
	}

	for (int i = 0; i < N; i++)
	{
		chromosomes[i] = new int[l];
		chromos[i] = new int[nMatch];
		myrand.uniformArray(chromosomes[i], l, 0, l-1);
		for (int j = 0; j < l; ++j)
		{
			if (chromosomes[i][j] < nMatch)
			{
				chromos[i][chromosomes[i][j]] = j;
			}
		}
	}

	int temp[20] = {
		0, 45, 45, 45,
		1, 45, 45, 45,
		2, 45, 45, 45,
		45, 45, 45, 45,
		45, 45, 45, 45,
	};

	int oldGood = 0, good = 0;
	for (int i = 0; i < N; ++i)
	{
		if (evaluate(chromosomes[i]) != INT_MIN)
		{
			oldGood++;
		}
		if (eval(chromos[i]) > -900) good++;
	}

	cout << "oldGood: " << oldGood << " " << "good: " << good << endl;;


	for (int i = 0; i < N; i++)
	{
		delete chromosomes[i];
		delete chromos[i];
	}

	system("pause");
	return 0;
}

void getIndexFromMatchId(int* a, int* b, int matchId)
{
	*a = 0;
	*b = 0;
	int i = nTeam - 1;
	int j = 1;
	while (matchId >= 0)
	{
		if (matchId < i)
		{
			*b = matchId + j;
			return;
		}

		j++;
		*a = *a + 1;
		matchId -= i;
		i--;
	}

	cerr << "This should never happen" << endl;
	exitWithError(-1);
}

void exitWithError(int ret) 
{
	system("pause");
	exit(ret);
}

int evaluate(int* chromosome)
{
	int fitness = 0;

	for (int w = 0; w < nWeek; ++w)
	{
		unordered_set<int> b2bSet, b2bSetLastDay;
		for (int d = 0; d < 5; ++d) // 5 days per week
		{
			int base = (w * 20 + d * 4);
			int t1, t2, t3, t4;
			for (int i = 0; i < 4; ++i)
			{
				if (chromosome[base + i] >= nMatch) continue;
				getIndexFromMatchId(&t1, &t2, chromosome[base + i]);
				fitness += preferenceArray[t1][d * 2 + i / 2];
				fitness += preferenceArray[t2][d * 2 + i / 2];

				for (int j = i + 1; j < 4; ++j) 
				{
					if (chromosome[base + j] >= nMatch) continue;
					getIndexFromMatchId(&t3, &t4, chromosome[base + j]);
					// do not allow game from the same day
					if (t1 == t3 || t1 == t4 || t2 == t3 || t2 == t4) return INT_MIN;
				}

				for (int dComp = d + 1; dComp < d + 3 && dComp < 5; ++dComp)
				{
					for (int k = 0; k < 4; ++k)
					{
						if (chromosome[w * 20 + dComp * 4 + k] >= nMatch) continue;
						getIndexFromMatchId(&t3, &t4, chromosome[w * 20 + dComp * 4 + k]);
						if (t1 == t3 || t1 == t4 || t2 == t3 || t2 == t4)
						{
							if (dComp - d == 1)
							{
								fitness -= 3;
								if (t1 == t3 || t1 == t4) 
								{
									if (b2bSetLastDay.find(t1) != b2bSetLastDay.end()) return INT_MIN;
									else b2bSet.insert(t1);
								}
								if (t2 == t3 || t2 == t4)
								{
									if (b2bSetLastDay.find(t2) != b2bSetLastDay.end()) return INT_MIN;
									else b2bSet.insert(t2);
								}
							}
							else if (dComp - d == 2) fitness -= 1;
						}
					}
				}
			}
			b2bSetLastDay = b2bSet;
			b2bSet.clear();
		}
	}

	return fitness;
}

int eval(int* chromos)
{
	int days[nWeek * 5] = {0};
	int fitness[nTeam] = { 0 };

	for (int i = 0; i < nTeam; i++)
	{
		for (int j = 0; j < nTeam - 1; j++)
		{
			fitness[i] += preferenceArray[i][(chromos[matches[i][j]] % 20) / nField];
			if (days[chromos[matches[i][j]] / (2 * nField)]) {
				fitness[i] -= 1000;
				//cout << "same day collision" << endl;
			}
			else days[chromos[matches[i][j]] / 2 / nField] = 1;
		}
		
		for (int j = 0; j < nWeek; ++j)
		{
			if (days[j * 5] && days[j * 5 + 1]) fitness[i] -= 3;
			for (int k = 2; k < 5; ++k)
			{
				if (days[j * 5 + k - 2] && days[j * 5 + k]) fitness[i] -= 1;
				if (days[j * 5 + k - 1] && days[j * 5 + k]) fitness[i] -= 3;
				if (days[j * 5 + k - 2] && days[j * 5 + k - 1] && days[j * 5 + k]) fitness[i] -= 1000;
				//cout << "b2b2b" << endl;
			}
		}
		memset(days, 0, sizeof(int) * nWeek * 5);
	}

	int f = 0;
	for (int i = 0; i < nTeam; ++i)
	{
		f += fitness[i];
	}

	return f;

}
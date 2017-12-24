#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <climits>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
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
void start(int**, int**);
void orderXO(int* a, int* b, int* c, int* d);
void partiallyMappedXO(int *a, int *b, int *c, int *d);

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
	int *chromos[N];
	int *chromosBuffer[N];

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
		chromos[i] = new int[nMatch];
		chromosBuffer[i] = new int[nMatch];
		myrand.uniformArray(chromos[i], nMatch, 0, l-1);
	}

	int good = 0;
	for (int i = 0; i < N; ++i)
	{
		if (eval(chromos[i]) > -900) good++;
	}

	cout << "good: " << good << endl;

	start(chromos, chromosPool);

	for (int i = 0; i < N; i++)
	{
		delete chromos[i];
		delete chromosPool[i];
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

struct sort_pred {
	bool operator()(const pair<int*, int> &left, const pair<int*, int> &right) {
		return left.second < right.second;
	}
};

void start(int** chromos, int** chromosBuffer)
{
	int order[N];
	vector<pair<int*, int>> vec;

	while (1)
	{
		myrand.uniformArray(order, N, 0, N);
		for (int i = 0; i < N; i += 2)
		{
			partiallyMappedXO(chromos[order[i]], chromos[order[i+1]], chromosBuffer[i], chromosBuffer[i+1]);
		}

		for (int i = 0; i < N; ++i)
		{
			vec.push_back(pair<int*, int>(chromos[i], eval(chromos[i])));
			vec.push_back(pair<int*, int>(chromosBuffer[i], eval(chromosBuffer[i])));
		}

		sort(vec.begin(), vec.end(), sort_pred());

	}
}

void orderXO(int* a, int* b, int* c, int* d)
{

}

void partiallyMappedXO(int *a, int *b, int *c, int *d)
{
	unordered_map<int, int> ah, bh;
	for (int i = 0; i < nMatch; ++i)
	{
		ah[a[i]] = i;
		bh[b[i]] = i;
	}

	memset(c, -1, sizeof(int) * nMatch);
	memset(d, -1, sizeof(int) * nMatch);

	// generate the random segment
	int head, tail;
	head = myrand.uniformInt(0, nMatch);
	tail = myrand.uniformInt(0, nMatch);

	// make sure that head > tail
	if (head > tail)
	{
		head = head ^ tail;
		tail = head ^ tail;
		head = head ^ tail;
	}

	for (int i = head; i <= tail; ++i)
	{
		c[i] = a[i];
		if (bh[a[i]] < head || bh[a[i]] > tail)
		{
			int s = b[i];
			while (ah[s] >= head && ah[s] <= tail)
			{
				s = b[ah[s]];
			}
			c[bh[a[i]]] = s;
		}

		d[i] = b[i];
		if (ah[b[i]] < head || ah[b[i]] > tail)
		{
			int s = a[i];
			while (bh[s] >= head && bh[s] <= tail)
			{
				s = a[bh[s]];
			}
			d[ah[b[i]]] = s;
		}
	}

	for (int i = 0; i < nMatch; ++i)
	{
		if (c[i] == -1)
		{
			c[i] = b[i];
		}

		if (d[i] == -1)
		{
			d[i] = a[i];
		}
	}

}
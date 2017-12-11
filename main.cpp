#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <climits>
#include "myrand.h"

#define N 100 // population size 

using namespace std;
void getIndexFromMatchId(int*, int*, int, int);
void exitWithError(int);

MyRand myrand;

// 10 time slots per week
int nWeek = 8;
int nTeam = 10;
int nField = 2;
int nTimeSlot = nWeek * 10 * 2;

int l = nTimeSlot * nField;

// each team has a match with every other team
int nMatch = nTeam * (nTeam - 1) / 2;

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

	for (int i = 0; i < N; i++)
	{
		chromosomes[i] = new int[l];
		myrand.uniformArray(chromosomes[i], 80, 0, 79);
	}




	for (int i = 0; i < N; i++)
	{
		delete chromosomes[i];
	}

	system("pause");
	return 0;
}

void getIndexFromMatchId(int* a, int* b, int matchId, int nTeam)
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
		for (int t = 0; t < 10; ++t)
		{
			
		}
	}

}
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <climits>
#include <iostream>
#include <vector>
#include <cstring>
#include <unordered_set>
#include <unordered_map>
#include "myrand.h"

#define N 1000 // population size
#define nWeek 4 // number of week
#define nTeam 6 // number of team
#define nField 2 // number of field
#define nTimeSlot nWeek*10 // number of timeslot
#define l nTimeSlot*nField
#define unbalancedPenalty nTeam
#define gpc 0.7

using namespace std;
void getIndexFromMatchId(int*, int*, int);
void exitWithError(int);
int evaluate(int*);
int eval(int*, int* contribution = 0, bool log = false);
int start(int**, int**);
void orderXO(int* a, int* b, int* c, int* d);
void partiallyMappedXO(int *a, int *b, int *c, int *d);
void partiallyMappedXOT(int *a, int *b, int *c, int *d);
void cycleXO(int *a, int *b, int *c, int *d);
void baselineRandom(int generation);
void baselineGreedy(int generation);
void baselineGreedyNFE(int nfe);
void greedyInitialization(int ** chromos);

MyRand myrand;

// each team has a match with every other team
int nMatch = nTeam * (nTeam - 1) / 2;
vector<int> matches[nTeam];

// date preference [Team][Date]
int preferenceArray[17][10] =
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
	{ 1, 0, -1, 0, -1, 2, -2, 0, 0, 1, },
	{ -2, -1, 0, 1, 0, 0, 1, 0, -1, 2, },
	{ 1, 0, -1, -1, 0, 0, 1, 2, -2, 0, },
	{ 0, -1, 1, 0, 1, 0, 0, 2, -1, -2, },
	{ 0, 0, 0, -1, 2, 1, -1, 0, -2, 1, },
	{ 1, -1, -2, 0, 0, 0, 2, 1, 0, -1, },
	{ 1, -1, -2, 0, 0, 0, 2, 1, 0, -1, },
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

	//greedyInitialization(chromos);

	int good = 0;
	for (int i = 0; i < N; ++i)
	{
		if (eval(chromos[i]) > -900) good++;
	}

	int generation = start(chromos, chromosBuffer);

	//baselineGreedy(100);
	baselineGreedyNFE(generation * N * 2);

	for (int i = 0; i < N; i++)
	{
		delete chromos[i];
		delete chromosBuffer[i];
	}

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
	cout << "Press Enter to exit";
	getc(stdin);
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

int eval(int* chromos, int* contribution, bool log)
{
	int days[nWeek * 5] = {0};
	int fitness[nTeam] = { 0 };

	if (contribution != 0)
	{
		memset(contribution, 0, sizeof(int) * nMatch);
	}

	for (int i = 0; i < nTeam; i++)
	{
		for (int j = 0; j < nTeam - 1; j++)
		{
			fitness[i] += preferenceArray[i][(chromos[matches[i][j]] % 20) / nField];
			if (contribution != 0)
			{
				contribution[matches[i][j]] += preferenceArray[i][(chromos[matches[i][j]] % 20) / nField];
			}
			if (days[chromos[matches[i][j]] / (2 * nField)]) {
				fitness[i] -= 1000;
				if (contribution != 0)
				{
					contribution[matches[i][j]] -= 1000;
				}
			}
			else days[chromos[matches[i][j]] / 2 / nField] = 1;
		}

		if (contribution != 0)
		{
			for (int j = 0; j < nTeam - 1; j++)
			{
				int day = chromos[matches[i][j]] / 2 / nField;
				if (((day - 2) / 5) == (day / 5) && days[day] && days[day - 2])
				{
					contribution[matches[i][j]] -= 2;
				}
				if (((day - 1) / 5) == (day / 5) && days[day] && days[day - 1])
				{
					contribution[matches[i][j]] -= 1;
				}
				if (((day + 1) / 5) == (day / 5) && days[day] && days[day + 1])
				{
					contribution[matches[i][j]] -= 1;
				}
				if (((day + 2) / 5) == (day / 5) && days[day] && days[day + 2])
				{
					contribution[matches[i][j]] -= 2;
				}
			}
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
			if (!days[j * 5] && !days[j * 5 + 1] && !days[j * 5 + 2] && !days[j * 5 + 3] && !days[j * 5 + 4])
			{
				fitness[i] -= 3;
			}
		}
		memset(days, 0, sizeof(int) * nWeek * 5);
	}

	int f = 0;
	int min = INT_MAX;
	for (int i = 0; i < nTeam; ++i)
	{
		if(log) cout << fitness[i] << " ";
		f += fitness[i];
		if (fitness[i] < min) min = fitness[i];
	}
	if (log) cout << endl;

	return f + unbalancedPenalty * min;

}

struct sort_pred {
	bool operator()(const pair<int*, int> &left, const pair<int*, int> &right) {
		return left.second > right.second;
	}
};

int start(int** chromos, int** chromosBuffer)
{
	int order[N];
	vector<pair<int*, int>> vec;

	int generation = 1;

	while (1)
	{
		myrand.uniformArray(order, N, 0, N-1);
		for (int i = 0; i < N; i += 2)
		{
			// if ((i % (N/50)) == 0) cout << "*";
			if (myrand.flip(gpc)) partiallyMappedXOT(chromos[order[i]], chromos[order[i+1]], chromosBuffer[i], chromosBuffer[i+1]);
			else
			{
				chromosBuffer[i][0] = -1;
				chromosBuffer[i+1][0] = -1;
			}
		}

		for (int i = 0; i < N; ++i)
		{
			vec.push_back(pair<int*, int>(chromos[i], eval(chromos[i])));
			if (chromosBuffer[i][0] != -1) vec.push_back(pair<int*, int>(chromosBuffer[i], eval(chromosBuffer[i])));
			else
			{
				vec.push_back(pair<int*, int>(chromosBuffer[i], INT_MIN));
			}
		}

		sort(vec.begin(), vec.end(), sort_pred());

		int sum = 0;
		int count = 0;

		for (int i = 0; i < N; ++i)
		{
			if (vec[i].second > -500)
			{
				count++;
			}
			sum += vec[i].second;
			chromos[i] = vec[i].first;
			chromosBuffer[i] = vec[i + N].first;
		}

		//cout << "Generation: " << generation << endl;
		//cout << "the average fitness of this generation is: " << float(float(sum) / float(N)) << endl;
		//cout << "valid chromosome in this generation is: " << count << endl;
		//cout << "the best fitness value of this generation is: " << vec[0].second << endl;

		cout << generation << ", ";
		cout << float(float(sum) / float(N)) << ", ";
		cout << count << ", ";
		cout << vec[0].second << endl;

		generation++;

		if ((float(vec[0].second) - float(float(sum) / float(N))) < 0.0001)
		{
			cout << "Converged at generation:" << generation << endl;
			eval(vec[0].first, 0, true);
			break;
		}

		vec.clear();
	}
	return generation;
}

void orderXO(int* a, int* b, int* c, int* d)
{
	unordered_map<int, int> ah, bh, ch, dh;
	for (int i = 0; i < nMatch; ++i)
	{
		ah[a[i]] = i;
		bh[b[i]] = i;
	}

	memcpy(c, a, sizeof(int)* nMatch);
	memcpy(d, b, sizeof(int)* nMatch);

	// generate the random segment
	int head, tail;
	head = myrand.uniformInt(0, nMatch - 1);
	tail = myrand.uniformInt(0, nMatch - 1);

	// make sure that head > tail
	if (head > tail)
	{
		head = head ^ tail;
		tail = head ^ tail;
		head = head ^ tail;
	}

	int idx;
	int idxa = (tail + 1) % nMatch;
	int	idxb = idxa;

	for (int i = 0; i < nMatch - (tail - head + 1); ++i)
	{
		idx = (tail + i + 1) % nMatch;
		while (ah.find(b[idxb]) != ah.end() && ah[b[idxb]] >= head && ah[b[idxb]] <= tail)
		{
			idxb++;
			idxb = idxb % nMatch;
		}
		c[idx] = b[idxb];

		while (bh.find(a[idxa]) != bh.end() && bh[a[idxa]] >= head && bh[a[idxa]] <= tail)
		{
			idxa++;
			idxa = idxa % nMatch;
		}
		d[idx] = a[idxa];

		idxa++;
		idxa = idxa % nMatch;
		idxb++;
		idxb = idxb % nMatch;

	}
}

void partiallyMappedXO(int *a, int *b, int *c, int *d)
{
	if (a == b || a == c || a == d || b == c || b == d || c == d)
	{
		cout << "duplicated chromosome pointer" << endl;
		getc(stdin);
	}

	unordered_map<int, int> ah, bh;
	for (int i = 0; i < nMatch; ++i)
	{
		ah[a[i]] = i;
		bh[b[i]] = i;
	}

	memset(c, -1, sizeof(int)* nMatch);
	memset(d, -1, sizeof(int)* nMatch);

	// generate the random segment
	int head, tail;
	head = myrand.uniformInt(0, nMatch-1);
	tail = myrand.uniformInt(0, nMatch-1);

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
		if (bh.find(a[i]) != bh.end() && (bh[a[i]] < head || bh[a[i]] > tail))
		{
			int s = b[i];
			while (ah.find(s) != ah.end() && ah[s] >= head && ah[s] <= tail)
			{
				s = b[ah[s]];
			}
			c[bh[a[i]]] = s;
		}

		d[i] = b[i];
		if (ah.find(b[i]) != ah.end() && (ah[b[i]] < head || ah[b[i]] > tail))
		{
			int s = a[i];
			while (bh.find(s) != bh.end() && bh[s] >= head && bh[s] <= tail)
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

	ah.clear();
	bh.clear();
	for (int i = 0; i < nMatch; ++i)
	{
		ah[c[i]] = i;
		bh[d[i]] = i;
	}
}

void partiallyMappedXOT(int *a, int *b, int *c, int *d)
{
	if (a == b || a == c || a == d || b == c || b == d || c == d)
	{
		cout << "duplicated chromosome pointer" << endl;
		getc(stdin);
	}

	unordered_set<int> chromosToChange;
	unordered_map<int, int> ah, bh;
	int *contributionA = new int[nMatch];
	int *contributionB = new int[nMatch];
	float *pc = new float[nMatch];
	memset(pc, 0, sizeof(float) * nMatch);

	eval(a, contributionA);
	eval(b, contributionB);

	for (int i = 0; i < nMatch; ++i)
	{
		ah[a[i]] = i;
		bh[b[i]] = i;
		if (contributionA[i] < -500)
		{
			pc[i] += 0.25;
		}
		else if (contributionA[i] < 0)
		{
			pc[i] += 0.15;
		}
		else pc[i] += 0.1;

		if (contributionB[i] < -500)
		{
			pc[i] += 0.25;
		}
		else if (contributionB[i] < 0)
		{
			pc[i] += 0.15;
		}
		else pc[i] += 0.1;
	}

	memset(c, -1, sizeof(int)* nMatch);
	memset(d, -1, sizeof(int)* nMatch);

	for (int i = 0; i < nMatch; ++i)
	{
		if (myrand.flip(pc[i]))
		{
			chromosToChange.insert(i);
		}
	}

	for (int i : chromosToChange)
	{
		c[i] = a[i];
		if (bh.find(a[i]) != bh.end() && chromosToChange.find(bh[a[i]]) == chromosToChange.end())
		{
			int s = b[i];
			while (ah.find(s) != ah.end() && chromosToChange.find(ah[s]) != chromosToChange.end())
			{
				s = b[ah[s]];
			}
			c[bh[a[i]]] = s;
		}

		d[i] = b[i];
		if (ah.find(b[i]) != ah.end() && chromosToChange.find(ah[b[i]]) == chromosToChange.end())
		{
			int s = a[i];
			while (bh.find(s) != bh.end() && chromosToChange.find(bh[s]) != chromosToChange.end())
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

	ah.clear();
	bh.clear();
	for (int i = 0; i < nMatch; ++i)
	{
		ah[c[i]] = i;
		bh[d[i]] = i;
	}

	delete[] contributionA;
	delete[] contributionB;
	delete[] pc;
}

void cycleXO(int *a, int *b, int *c, int *d)
{
	if (a == b || a == c || a == d || b == c || b == d || c == d)
	{
		cout << "duplicated chromosome pointer" << endl;
		getc(stdin);
	}

	unordered_map<int, int> ah, bh;
	for (int i = 0; i < nMatch; ++i)
	{
		ah[a[i]] = i;
		bh[b[i]] = i;
	}

	memcpy(c, a, sizeof(int) * nMatch);
	memcpy(d, b, sizeof(int) * nMatch);

	bool backward = true;
	int pivot = myrand.uniformInt(0, nMatch - 1);
	int idx = pivot;

	while (1)
	{
		c[idx] = b[idx];
		d[idx] = a[idx];

		if (ah.find(b[idx]) == ah.end() || ah[b[idx]] == pivot)
		{
			if (ah.find(b[idx]) != ah.end())
			{
				backward = false;
			}
			break;
		}

		else idx = ah[b[idx]];
	}

	idx = pivot;
	while (backward)
	{
		if (bh.find(a[idx]) == bh.end())
		{
			break;
		}

		idx = bh[a[idx]];

		c[idx] = b[idx];
		d[idx] = a[idx];

	}

}

void baselineRandom(int generation)
{
	int *chromos, *best;
	chromos = new int[nMatch];
	best = new int[nMatch];
	int bestFitness = INT_MIN, f;
	for (int i = 0; i < generation * N; ++i)
	{
		myrand.uniformArray(chromos, nMatch, 0, l-1);
		f = eval(chromos);
		if (f > bestFitness)
		{
			memcpy(best, chromos, sizeof(int)* nMatch);
			bestFitness = f;
		}
	}

	cout << "Best fitness using pure random baseline: " << bestFitness << endl;

	delete [] chromos;
	delete [] best;
}

void baselineGreedy(int generation)
{
	int * chromos = new int[nMatch];
	int * buffer = new int[nMatch];
	int best = INT_MIN;

	for (int i = 0; i < generation; ++i)
	{
		if ((i % (generation / 50)) == 0) cout << "*";

		myrand.uniformArray(chromos, nMatch, 0, l-1);

		unordered_map<int, int> h;

		for (int i = 0; i < nMatch; ++i)
		{
			h[chromos[i]] = i;
		}

		while (1) 
		{
			int jIdx, kIdx, temp;
			int currentBestFitness = eval(chromos);
			int f = 0;
			for (int j = 0; j < nMatch; ++j)
			{
				for (int k = 0; k < l; ++k)
				{
					memcpy(buffer, chromos, sizeof(int)* nMatch);

					if (buffer[j] == k) continue;
					temp = buffer[j];
					buffer[j] = k;
					if (h.find(k) != h.end())
					{
						buffer[h[k]] = temp;
					}

					f = eval(buffer);
					if (f > currentBestFitness)
					{
						jIdx = j;
						kIdx = k;
						currentBestFitness = f;
					}
				}
			}

			if (eval(chromos) == currentBestFitness)
			{
				if (currentBestFitness > best) 	best = currentBestFitness;
				break;
			}
			else
			{
				temp = chromos[jIdx];
				chromos[jIdx] = kIdx;
				h.erase(temp);
				if (h.find(kIdx) != h.end())
				{
					chromos[h[kIdx]] = temp;
					h[temp] = h[kIdx];
				}
				h[kIdx] = jIdx;
			}
		}
	}

	cout << endl << "Best fitness using greedy baseline: " << best << endl;

	delete [] chromos;
	delete [] buffer;
}

void baselineGreedyNFE(int nfe)
{
	int * chromos = new int[nMatch];
	int * buffer = new int[nMatch];
	int best = INT_MIN;
	int count = nfe;
	int loops = 0;

	while (1)
	{
		loops++;
		myrand.uniformArray(chromos, nMatch, 0, l - 1);

		unordered_map<int, int> h;

		for (int i = 0; i < nMatch; ++i)
		{
			h[chromos[i]] = i;
		}

		while (1)
		{
			int jIdx, kIdx, temp;
			int currentBestFitness = eval(chromos);
			int f = 0;
			for (int j = 0; j < nMatch; ++j)
			{
				for (int k = 0; k < l; ++k)
				{
					memcpy(buffer, chromos, sizeof(int)* nMatch);

					if (buffer[j] == k) continue;
					temp = buffer[j];
					buffer[j] = k;
					if (h.find(k) != h.end())
					{
						buffer[h[k]] = temp;
					}

					f = eval(buffer);
					count--;
					if ((count % (nfe / 50)) == 0) cerr << "*";
					if (f > currentBestFitness)
					{
						jIdx = j;
						kIdx = k;
						currentBestFitness = f;
					}
				}
			}

			if (eval(chromos) == currentBestFitness)
			{
				if (currentBestFitness > best) 	best = currentBestFitness;
				break;
			}
			else
			{
				temp = chromos[jIdx];
				chromos[jIdx] = kIdx;
				h.erase(temp);
				if (h.find(kIdx) != h.end())
				{
					chromos[h[kIdx]] = temp;
					h[temp] = h[kIdx];
				}
				h[kIdx] = jIdx;
			}
		}
		if (count < 0) break;
	}

	cout << endl << "Best fitness using greedy NFE baseline: " << best << endl;
	cout << "Greedy NFE tested " << loops << " chromosomes" << endl;

	delete[] chromos;
	delete[] buffer;
}

void greedyInitialization(int ** chromos)
{
	int *buffer = new int[nMatch];

	for (int i = 0; i < N; ++i)
	{
		if ((i % (N / 50)) == 0) cout << "*";
		myrand.uniformArray(chromos[i], nMatch, 0, l-1);

		unordered_map<int, int> h;

		for (int j = 0; j < nMatch; ++j)
		{
			h[chromos[i][j]] = j;
		}

		int currentBestFitness = eval(chromos[i]);
		while (1)
		{
			int jIdx, kIdx, temp;
			int f = 0;
			for (int j = 0; j < nMatch; ++j)
			{
				for (int k = j+1; k < nMatch; ++k)
				{
					memcpy(buffer, chromos[i], sizeof(int)* nMatch);

					buffer[j] ^= buffer[k];
					buffer[k] ^= buffer[j];
					buffer[j] ^= buffer[k];

					f = eval(buffer);
					if (f > currentBestFitness)
					{
						jIdx = j;
						kIdx = k;
						currentBestFitness = f;
					}
				}
			}

			if (eval(chromos[i]) == currentBestFitness)
			{
				break;
			}
			else
			{
				chromos[i][jIdx] ^= chromos[i][kIdx];
				chromos[i][kIdx] ^= chromos[i][jIdx];
				chromos[i][jIdx] ^= chromos[i][kIdx];
			}
		}
	}

	cout << endl;

	delete[] buffer;
}

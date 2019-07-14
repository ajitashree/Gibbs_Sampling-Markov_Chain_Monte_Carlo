#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <ctime>
#include <cassert>
#include <cstdlib>
#include <map>
#include <ctime>
#include "Factor.cpp"

using namespace std;
double time1 = 0.0;
double LikeLihood = 0;
map <char, int> char2int;
map <int, char> int2char;
vector <string> correctWords;
long totalChars = 0, charsSame = 0; 
long totalWords = 0, wordsSame = 0;
double ocr[1000][10], transition[10][10];
char charaters[10] = {'e', 't', 'a', 'o', 'i', 'n', 's', 'h', 'r', 'd'};

void computeTable(long N, Factor *f, vector <long> list, long task, vector <long> vec, long imgID)
{
	if(N == 0)
	{
		double p = 0;
		if (task == 1) p = ocr[imgID][list[0]];
		else if (task == 2) p = transition[list[0]][list[1]];
		else p = (list[0] == list[1]) ? log(5) : 0; // DANGER

		f->table.push_back(make_pair(list, p));
		return;
	}
	else
	{
		for (long i = 0; i < 10; i++)
		{
			list.push_back(i);
			computeTable(N - 1, f, list, task, vec, imgID);
			list.pop_back();
		}
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void updateMatrix(long** matrix, vector <long> image, long size, long begin, 
	vector <Factor*> &initial_facts)
{
	vector <long> vec, tmp;
	for(long i = begin; i < size; i++)
	{
		for (long j = i + 1; j < size; j++)
		{
			if(image[i] == image[j])
			{
				vec.push_back(i);
				vec.push_back(j);

				matrix[i][j] = 1;
				matrix[j][i] = 1;

				Factor* f = createFactor(vec, 10);
				computeTable(vec.size(), f, tmp, 3, vec, 0);
				initial_facts.push_back(f);

				vec.pop_back();
				vec.pop_back();
			}
		}
	}
}
void updatePairs(long** matrix, vector <long> image, long size1, 
	vector <Factor*> &initial_facts)
{
	vector <long> vec, tmp;
	for(long i = 0; i < size1; i++)
	{
		for (long j = size1; j < image.size(); j++)
		{
			if(image[i] == image[j])
			{
				vec.push_back(i);
				vec.push_back(j);

				matrix[i][j] = 1;
				matrix[j][i] = 1;

				Factor* f = createFactor(vec, 10);
				computeTable(vec.size(), f, tmp, 4, vec, 0);
				initial_facts.push_back(f);

				vec.pop_back();
				vec.pop_back();
			}
		}
	}
}
void updateTans(long** matrix, long size, long begin, vector <Factor*> &initial_facts)
{
	vector <long> vec, tmp;
	for (long i = begin; i < size - 1; i++)
	{
		vec.push_back(i);
		vec.push_back(i + 1);

		matrix[i][i + 1] = 1;
		matrix[i + 1][i] = 1;

		Factor* f = createFactor(vec, 10);
		computeTable(vec.size(), f, tmp, 2, vec, 0);
		initial_facts.push_back(f);

		vec.pop_back();
		vec.pop_back();
	}
}

void initialize_matrix(long** matrix, vector <long> image, long size1, long size2, 
	long category, vector <Factor*> &initial_facts)
{
	if (category >= 1)
	{
		vector <long> vec, tmp;
		for (long i = 0; i < image.size(); i++)
		{
			vec.push_back(i);
			Factor *f = createFactor(vec, 10);
			computeTable(vec.size(), f, tmp, 1, vec, image[i]);
			initial_facts.push_back(f);
			vec.pop_back();
		}
	}

	if (category >= 2)
	{
		updateTans(matrix, size1, 0, initial_facts);
		updateTans(matrix, image.size(), size1, initial_facts);
	}

	if (category >= 3)
	{
		updateMatrix(matrix, image, size1, 0, initial_facts);
		updateMatrix(matrix, image, image.size(), size1, initial_facts);
	}

	if (category >= 4)
		updatePairs(matrix, image, size1, initial_facts);

}
long getMaxval(long imageID_index, vector <vector <long> > samples_gen, long trueVal)
{
	long assign;
	long maxVal = -1;
	map <long, long> val;

	for (auto sample : samples_gen)
	{
		if (val.find(sample[imageID_index]) == val.end()) 
			val[sample[imageID_index]] = 0;

		val[sample[imageID_index]]++;
	}
	for (auto v : val)
	{
		if (v.second > maxVal)
		{
			maxVal = v.second;
			assign = v.first;
		}
	}
	LikeLihood += log(val[trueVal]/(double)samples_gen.size()); // Likelihood Computation....
	return assign;
}


string prediction(vector <vector <long> > samples_gen, vector <long> image, string finalWord)
{
	string str = "";
	vector <long> vec_str;
	long start_s = clock();
	for (long i = 0; i < image.size(); i++)
	{
		vec_str.push_back(getMaxval(i, samples_gen, char2int[finalWord[i]]));
	}
	for (auto i : vec_str)
		str += int2char[i];
	long stop_s = clock();
	time1 += 1.0 * (stop_s - start_s) / CLOCKS_PER_SEC;
	return str;
}



bool issubset(long id, vector <long> scp)
{
	for (long i = 0; i < scp.size(); i++)
	{
		if (id == scp[i]) return true;
	}
	return false;
}

vector <long> getSample(vector <long> sample, vector <Factor*> &initial_facts, long turn)
{
	double Z = 0;
	long assign;
	vector <Factor*> fact;

	for (auto f : initial_facts)
	{
		if (issubset(turn, f->scope)) fact.push_back(f);
	}

	for (long f = 0; f < fact.size(); f++)
	{
		for (long s = 0; s < sample.size(); s++)
		{
			if (s != turn && issubset(s, fact[f]->scope))
			{
				fact[f] = reduction(fact[f], s, sample[s], 10);
			}
		}
	}

	Factor *Fact_mul = compute_psi(fact, 10);

	
	Fact_mul = normalize(Fact_mul);
	
	double r = ((double) rand() / (RAND_MAX));
	double CF = 0;

	for (long i = 0; i < Fact_mul->tableSize; i++)
	{
		CF += exp(Fact_mul->table[i].second);
		assign = Fact_mul->table[i].first[0];

		if (r <= CF) break;			
	}
	sample[turn] = assign;
	return sample;
}



// 1: OCR , 2 : OCR + Trans, 3 : OCR + Trans + Skip, 4 : OCR + Trans + Skip + Parskip
void Modelling(vector <long> image, long size1, long size2, string word1, string word2, long category,
	long burn_insamples, long gen_samples)
{
	vector <vector <long> > samples_gen;
	vector <Factor*> initial_facts;

	long size = image.size();
	long** matrix = new long*[size];

	for (long i  = 0; i < size; i++)
	{
		matrix[i] = new long[size];
		for (long j = 0; j < size; j++)
			matrix[i][j] = 0;
	}

	vector <long> samples;
	for (long i = 0; i < image.size(); i++)
		samples.push_back(0);

	initialize_matrix(matrix, image, size1, size2, category, initial_facts);

	
	long turn = 0;
	while(burn_insamples--)
	{
		samples = getSample(samples, initial_facts, turn);
		turn = (turn + 1) % image.size();
	}
	
	while(gen_samples--)
	{
		samples = getSample(samples, initial_facts, turn);
		samples_gen.push_back(samples);
		turn = (turn + 1) % image.size();
	}
	
	string str = prediction(samples_gen, image, word1 + word2);
	string str1 = str.substr(0, size1);
	string str2 = str.substr(size1, size2);

	if (str1 == word1) wordsSame++;
	if (str2 == word2) wordsSame++;

	for (long i = 0; i < str1.size(); i++)
	{
		if (str1[i] == word1[i]) charsSame++;
	}
	for (long j = 0; j < str2.size(); j++)
	{
		if (str2[j] == word2[j]) charsSame++;
	}
}

void printAccuracy()
{
	double charAcc, wordAcc, LL;
	long total = totalWords;
	charAcc = (double)charsSame/totalChars * 100;
	wordAcc = (double)wordsSame/totalWords * 100;
	LL = (double)LikeLihood/total;

	/*
	cout << endl;
	cout << "Char Wise Accuracy : ";
	cout << (double)charsSame/totalChars * 100 << endl;

	cout << "Word Wise Accurcay : ";
	cout << (double)wordsSame/totalWords * 100 << endl;

	cout << "Average Log LikeLihood : ";
	cout << (double)LikeLihood/total<< endl;
	cout << "Time Taken : ";
	cout << time1 << endl;
	cout << endl;*/
	printf("%lf %lf %lf %lfs\n\n", charAcc, wordAcc, LL, time1);
}

//RUN :  ./a.out category fileindex
int main(int argc, char const *argv[])
{

	char ch, ch1;
	long n, index = 0;
	double prob;
	long burn_insamples = 1000, gen_samples = 10000;
	string str, imgID1, imgID2, space;

	long category = atoi(argv[1]);
	long fileindex = atoi(argv[2]);

	ifstream ocrDat, transDat, dataDat, wordsDat;
	ocrDat.open("./OCRdataset-2/potentials/ocr.dat");
	transDat.open("./OCRdataset-2/potentials/trans.dat");

	if (fileindex == 1)
	{
		dataDat.open("./OCRdataset-2/data/data-tree.dat");
		wordsDat.open("./OCRdataset-2/data/truth-tree.dat");
	}
	else if (fileindex == 2)
	{
		dataDat.open("./OCRdataset-2/data/data-loops.dat");
		wordsDat.open("./OCRdataset-2/data/truth-loops.dat");
	}
	else if(fileindex == 3)
	{
		dataDat.open("./OCRdataset-2/data/data-treeWS.dat");
		wordsDat.open("./OCRdataset-2/data/truth-treeWS.dat");
	}
	else if (fileindex == 4)
	{
		dataDat.open("./OCRdataset-2/data/data-loopsWS.dat");
		wordsDat.open("./OCRdataset-2/data/truth-loopsWS.dat");
	}	
	
	for (int i = 0; i < 10; i++)
	{
		int2char[i] = charaters[i];
		char2int[charaters[i]] = i;
	}
	for (long i = 0; i < 10000; i++)
	{
		ocrDat >> n >> ch >> prob;
		ocr[n][char2int[ch]] = log(prob);
	}
	for (long i = 0; i < 100; i++)
	{
		transDat >> ch >> ch1 >> prob;
		transition[char2int[ch]][char2int[ch1]] = log(prob);
	}
	while(getline(wordsDat, str))
	{
		if (str.size() != 0)
		{
			correctWords.push_back(str);
			totalChars += str.size();
		}
	}

	totalWords = correctWords.size();

	for (int i = 0; i < 10; i++)
	{
		int2char[i] = charaters[i];
		char2int[charaters[i]] = i;
	}

	while(getline(dataDat, imgID1))
	{
		vector <long> image1, image2, image;
		stringstream img1(imgID1);

		while(img1 >> n)
		{
			image1.push_back(n);
			image.push_back(n);
		}

		getline(dataDat, imgID1);
		stringstream img2(imgID1);

		while(img2 >> n)
		{
			image2.push_back(n);
			image.push_back(n);
		}

		Modelling(image, image1.size(), image2.size(), correctWords[index], correctWords[index + 1], 
			category, burn_insamples, gen_samples);
		getline(dataDat, space);
		index += 2;
	}
	printAccuracy();
}
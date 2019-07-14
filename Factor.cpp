#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <cassert>
#include <algorithm>
#include <limits.h>
#include <map>
#include <set>
using namespace std;

struct Factor
{
	vector <long> scope; // index of 3-digit image ID's
	map <long, long> scpMap;
	long tableSize;
	vector <pair< vector <long>,double> > table;

};

void initialize_table(long N, Factor *f, vector <long> list, long dist_values, long val)
{
	if(N <= 0)
	{
		f->table.push_back(make_pair(list, val));
		return;
	}
	else
	{
		for (long i = 0; i < dist_values; i++)
		{
			list.push_back(i);
			initialize_table(N - 1, f, list, dist_values, val);
			list.pop_back();
		}
	}

}

Factor* createFactor(vector <long> _scope, long dist_values)
{
	Factor *newFactor = new Factor;
	newFactor->scope = _scope;
	newFactor->tableSize = pow(dist_values, _scope.size());
	return newFactor;
}

Factor *getFactorCPD(map <vector <long>, double> CPD, vector <long> _scope, map <long, long> _scpMap)
{
	Factor *newFactor = new Factor;
	newFactor->scope = _scope;
	newFactor->tableSize = CPD.size();
	newFactor->scpMap = _scpMap;

	for (auto i : CPD)
	{
		newFactor->table.push_back(make_pair(i.first, i.second));
	}
	return newFactor;
}

Factor *normalize(Factor *f)
{
	double Z = 0;
	for (long i = 0; i < f->tableSize; i++)
		Z += exp(f->table[i].second);

	if (Z != 0)
	{
		for(long i = 0; i < f->tableSize; i++)
			f->table[i].second = log(exp(f->table[i].second)/Z);
	}

	return f;
}

void makeZero(Factor *f)
{
	for (long i = 0; i < f->table.size(); i++)
		f->table[i].second = 0.0;
}

Factor *copyFactor(Factor *fact)
{
	Factor *newFact = createFactor(fact->scope, 10);
	for (long f = 0; f < fact->table.size(); f++)
	{
		newFact->table.push_back(fact->table[f]);
	}
	newFact->tableSize = newFact->table.size();
	return newFact;
}

vector<pair<long, long> > intersect_scope(vector <long> s1, vector <long> s2, vector <long > &union_Scope)
{
	vector <pair<long, long> > intersect;

	union_Scope = s1;
	for (long i = 0; i < s1.size(); i++)
	{
		for (long j = 0; j < s2.size(); j++)
		{
			if (s1[i] == s2[j]) 
			{
				intersect.push_back(make_pair(i, j));
				break;
			}
		}

	}
	bool flag;
	for (auto e : s2)
	{
		flag = false;
		for (auto e1 : s1)
		{
			if (e1 == e)
			{
				flag = true;
				break;
			}
		}
		if (!flag) 
			union_Scope.push_back(e);
	}
	return intersect;
}

bool matched(vector <long> assign1, vector <long> assign2, vector <long> s1, vector <long> s2)
{
	for (long i = 0; i < s1.size(); i++)
	{
		for (long j = 0; j < s2.size(); j++)
		{
			if (s1[i] == s2[j])
			{
				if (assign1[i] != assign2[j]) return false;
			}
		}
	}
	return true;
}

bool isSub(long id, vector <long> scp)
{
	for (long i = 0; i < scp.size(); i++)
	{
		if (id == scp[i]) return true;
	}
	return false;
}

void initialize_tablenew(long N, vector <vector <long> > &assignment, vector <long> list, long origN, 
		Factor *f)
{
	if(N <= 0)
	{
		f->table.push_back(make_pair(list, 0.0));
		return;
	}
	else
	{
		for (long i = 0; i < assignment[origN - N].size(); i++)
		{
			list.push_back(assignment[origN - N][i]);
			initialize_tablenew(N - 1, assignment, list, origN, f);
			list.pop_back();
		}
	}

}

Factor *mutiplynew(Factor* f1, Factor* f2, long dist_values)
{
	long N = 0;

	vector <long> vec, union_Scope, attrs;
	vector <vector <long> > assignment;
	vector <pair<long, long> > intersect;
	map <long, long> sMap;

	intersect = intersect_scope(f1->scope, f2->scope, union_Scope); 

	assert((intersect.size() != 0));

	N = union_Scope.size();

	for (long u = 0; u < union_Scope.size(); u++)
	{
		if (isSub(union_Scope[u], f1->scope))
		{
			long dist = f1->scpMap[union_Scope[u]];
			sMap[union_Scope[u]] = dist;
			for (long d = 0; d < dist; d++)
				attrs.push_back(d);
		}
		else
		{
			long dist = f2->scpMap[union_Scope[u]];
			sMap[union_Scope[u]] = dist;
			for (long d = 0; d < dist; d++)
				attrs.push_back(d);
		}
		assignment.push_back(attrs);
		attrs.clear();
	}
	Factor *tmp = createFactor(union_Scope, dist_values);
	tmp->scpMap = sMap;

	initialize_tablenew(N, assignment, vec, N, tmp);

	//cout << union_Scope.size() << endl;
	for (long u = 0; u < tmp->table.size(); u++)
	{
		for (long r = 0; r < f1->table.size(); r++)
		{
			if (matched(tmp->table[u].first, f1->table[r].first, union_Scope, f1->scope))
			{
				for (long n = 0; n < f2->table.size(); n++)
				{
					if (matched(tmp->table[u].first, f2->table[n].first, union_Scope, f2->scope))
					{
						//tmp->table[u].second = f1->table[r].second * f2->table[n].second; 
						tmp->table[u].second = f1->table[r].second + f2->table[n].second; 
						break;
					}
				}
				break;
			}
		}
	}
	tmp->tableSize = tmp->table.size();
	return tmp;
}

Factor *compute_psinew(vector <Factor*> list, long dist_values)
{
	Factor *psi = list[0];
	for(long i = 1; i < list.size(); i++)
	{
		psi = mutiplynew(psi, list[i], dist_values);
	}
	return psi;
}

// create a map and then multiply...............
Factor *mutiply(Factor* f1, Factor* f2, long dist_values)
{
	vector <long> vec, union_Scope;
	vector <pair<long, long> > intersect;
	intersect = intersect_scope(f1->scope, f2->scope, union_Scope); 

	assert((intersect.size() != 0));

	Factor *tmp = createFactor(union_Scope, dist_values);
	initialize_table(union_Scope.size(), tmp, vec, dist_values, 0);

	//cout << union_Scope.size() << endl;
	for (long u = 0; u < tmp->tableSize; u++)
	{
		for (long r = 0; r < f1->tableSize; r++)
		{
			if (matched(tmp->table[u].first, f1->table[r].first, union_Scope, f1->scope))
			{
				for (long n = 0; n < f2->tableSize; n++)
				{
					if (matched(tmp->table[u].first, f2->table[n].first, union_Scope, f2->scope))
					{
						//tmp->table[u].second = f1->table[r].second * f2->table[n].second; 
						tmp->table[u].second = f1->table[r].second + f2->table[n].second; 
						break;
					}
				}
				break;
			}
		}
	}

	return tmp;
}



bool is_subset(vector <long> scp1, vector <long> scp2)
{
	bool flag;
	for (long i = 0; i < scp1.size(); i++)
	{
		flag = false;
		for (long j = 0; j < scp2.size(); j++)
		{
			if (scp1[i] == scp2[j])
			{
				flag = true;
				break;
			}
		}
		if (!flag) return false;
	}
	return true;
}

double div(double a, double b)
{
	//if (a == 0 && b == 0) return 0;
	//return a/b;
	return a - b;
}

Factor *Divide(Factor* f1, Factor* f2, long dist_values)
{
	vector <long> vec;
	assert(is_subset(f2->scope, f1->scope));

	Factor *tmp = createFactor(f1->scope, dist_values);
	initialize_table(f1->scope.size(), tmp, vec, dist_values, 0);

	double prob = 0;
	for (long i = 0; i < tmp->tableSize; i++)
	{
		for (long j = 0; j < f2->tableSize; j++)
		{
			if (matched(tmp->table[i].first, f2->table[j].first, f1->scope, f2->scope))
			{	
				tmp->table[i].second = div(f1->table[i].second,f2->table[j].second);
				break;
			}
		}
	}
	return tmp;
}

Factor *reduction(Factor* f, long scopeID, long scopeVal, long dist_values)
{
	vector <long> vec;
	long index = -1, tableSize = 0;
	map <long, long> scpmappr = f->scpMap;
	vector <long> scope = f->scope;

	for (long i = 0; i < scope.size(); i++)
	{
		if (scope[i] == scopeID)
		{
			index = i;
			break;
		}
	}
	scope.erase(scope.begin() + index);
	scpmappr.erase(scopeID);
	Factor *tmp = createFactor(scope, dist_values);

	for (long r = 0; r < f->table.size(); r++)
	{
		if (f->table[r].first[index] == scopeVal)
		{
			vector <long> vec = f->table[r].first;
			vec.erase(vec.begin() + index);
			tmp->table.push_back(make_pair(vec, f->table[r].second));
		}
	}
	tmp->tableSize = tmp->table.size();
	tmp->scpMap = scpmappr;
	return tmp;
}


Factor *marginalize(Factor* f, long scopeID, long dist_values)
{
	vector <long> v;
	long index = -1;
	vector <long> scope = f->scope;

	for (long i = 0; i < scope.size(); i++)
	{
		if (scope[i] == scopeID)
		{
			index = i;
			break;
		}
	}
	scope.erase(scope.begin() + index);
	Factor *tmp = createFactor(scope, dist_values);

	initialize_table(scope.size(), tmp, v, dist_values, 0);

	for (long r = 0; r < f->tableSize; r++)
	{
		vector <long> vec = f->table[r].first;
		vec.erase(vec.begin() + index);

		for (long n = 0; n < tmp->table.size(); n++)
		{
			if (tmp->table[n].first == vec)
			{
				//tmp->table[n].second += (f->table[r].second);
				tmp->table[n].second += exp(f->table[r].second);
				break;
			}
		}

	}
	// new....
	for (long n = 0; n < tmp->table.size(); n++)
	{
		tmp->table[n].second = log(tmp->table[n].second);
	}
	tmp->tableSize = tmp->table.size();
	return tmp;

}

double maxi(double a, double b)
{
	if (a >= b) return a;
	return b;
}

Factor *max_marginalize(Factor* f, long scopeID, long dist_values)
{
	vector <long> v;
	long index = -1;
	vector <long> scope = f->scope;

	for (long i = 0; i < scope.size(); i++)
	{
		if (scope[i] == scopeID)
		{
			index = i;
			break;
		}
	}
	scope.erase(scope.begin() + index);
	Factor *tmp = createFactor(scope, dist_values);

	initialize_table(scope.size(), tmp, v, dist_values, LONG_MIN);

	vector <double> val;

	for (long r = 0; r < f->tableSize; r++)
	{
		vector <long> vec = f->table[r].first;
		vec.erase(vec.begin() + index);

		for (long n = 0; n < tmp->table.size(); n++)
		{

			if (tmp->table[n].first == vec)
			{
				double prob = tmp->table[n].second;
				tmp->table[n].second = maxi(prob, f->table[r].second);
				break;
			}
		}

	}

	tmp->tableSize = tmp->table.size();
	return tmp;

}

vector <long> removeVars(vector <long> psi_scope, vector <long> tao_scope)
{
	vector <long> diff;

	bool flag;
	for (auto t : psi_scope)
	{
		flag = false;
		for (long p = 0; p < tao_scope.size(); p++)
		{
			if (tao_scope[p] == t)
			{
				flag = true;
				break;
			}
		}
		if (!flag) diff.push_back(t);

	}
	return diff;
}

Factor* sumOutVariables(Factor* f, vector <long> psi_scope, vector <long> tao_scope, long dist_values)
{
	Factor *tmp = f;
	vector <long> varToRemove = removeVars(psi_scope, tao_scope);

	if (varToRemove.size() == 0)
		return f;

	for (auto v : varToRemove)
		tmp = marginalize(tmp, v, dist_values);
	return tmp;
}

Factor* max_sumOutVariables(Factor* f, vector <long> psi_scope, vector <long> tao_scope, long dist_values)
{
	Factor *tmp = f;
	vector <long> varToRemove = removeVars(psi_scope, tao_scope);

	if (varToRemove.size() == 0)
		return f;

	for (auto v : varToRemove)
		tmp = max_marginalize(tmp, v, dist_values);
	return tmp;
}

Factor *compute_psi(vector <Factor*> list, long dist_values)
{
	Factor *psi = list[0];
	for(long i = 1; i < list.size(); i++)
	{
		psi = mutiply(psi, list[i], dist_values);
	}
	return psi;
}

Factor *compute_tao(Factor* f, vector <long> scopeID, long dist_values)
{
	Factor *tao = f;
	for (auto s : scopeID)
	{
		tao = marginalize(tao, s, dist_values);
	}
	return tao;
}

void printFactor(Factor *f)
{
	cout << endl;
	for (auto e : f->scope)
		cout << e << " ";
	cout << "prob";
	cout << endl;
	for (long r = 0; r < f->table.size(); r++)
	{
		for (auto e : f->table[r].first)
			cout << e << " ";
		cout << f->table[r].second << endl;
	}
	cout << endl;

}

/*
int main(int argc, char const *argv[])
{
	Factor *fact1 = NULL, *fact2 = NULL, *fact3 = NULL, *fact4 = NULL;
	
	vector <long> scope;
	map <long, long> mapVal;
	scope.push_back(1);
	scope.push_back(2);
	scope.push_back(3);

	mapVal[1] = 2;
	mapVal[2] = 2;
	mapVal[3] = 3;

	fact1 = createFactor(scope, 2);
	fact2 = createFactor(scope, 2);

	fact1->scpMap = mapVal;
	
	vector <long> vec;
	vec.push_back(0);
	vec.push_back(0);
	vec.push_back(0);
	fact1->table.push_back(make_pair(vec, 1));
	vec.clear();

	vec.push_back(0);
	vec.push_back(0);
	vec.push_back(1);
	fact1->table.push_back(make_pair(vec, 2));
	vec.clear();

	vec.push_back(0);
	vec.push_back(0);
	vec.push_back(2);
	fact1->table.push_back(make_pair(vec, 2));
	vec.clear();

	vec.push_back(0);
	vec.push_back(1);
	vec.push_back(0);
	fact1->table.push_back(make_pair(vec, 3));
	vec.clear();
	vec.push_back(0);
	vec.push_back(1);
	vec.push_back(1);
	fact1->table.push_back(make_pair(vec, 2));
	vec.clear();

	vec.push_back(0);
	vec.push_back(1);
	vec.push_back(2);
	fact1->table.push_back(make_pair(vec, 4));
	vec.clear();

	vec.push_back(1);
	vec.push_back(0);
	vec.push_back(0);
	fact1->table.push_back(make_pair(vec, 2));
	vec.clear();

	vec.push_back(1);
	vec.push_back(0);
	vec.push_back(1);
	fact1->table.push_back(make_pair(vec, 5));
	vec.clear();


	vec.push_back(1);
	vec.push_back(0);
	vec.push_back(2);
	fact1->table.push_back(make_pair(vec, 2));
	vec.clear();


	vec.push_back(1);
	vec.push_back(1);
	vec.push_back(0);
	fact1->table.push_back(make_pair(vec, 2));
	vec.clear();


	vec.push_back(1);
	vec.push_back(1);
	vec.push_back(1);
	fact1->table.push_back(make_pair(vec, 7));
	vec.clear();


	vec.push_back(1);
	vec.push_back(1);
	vec.push_back(2);
	fact1->table.push_back(make_pair(vec, 2));
	vec.clear();


	printFactor(fact1);

	fact2 = reduction(fact1, 2, 1, 2);
	printFactor(fact2);
	//fact3 = marginalize(fact1, 2, 2);
	fact3 = mutiplynew(fact1, fact2, 2);
	fact4 = mutiply(fact1, fact2, 2);
	printFactor(fact3);
	printFactor(fact4);

	
	vector <Factor*> fact4;
	fact4.push_back(fact1);
	fact4.push_back(fact2);
	fact4.push_back(fact3);
	printFactor(fact1);
	printFactor(fact2);
	printFactor(fact3);
	Factor *fact5 = reduction(fact3, 1, 0, 2);
	printFactor(fact5);
	
	Factor *fact5 = max_marginalize(fact3, 3, 2);
	printFactor(fact5);
	
	Factor* Fpsi = compute_psi(fact4, 2);
	//printFactor(Fpsi);

	vector <long> s;
	s.push_back(1);
	s.push_back(2);
	Factor* nn = compute_tao(Fpsi, s, 2);
	//printFactor(nn);

	cout << "fact1" << endl;
	printFactor(fact1);
	cout << "fact2" << endl;
	printFactor(fact2);
	Factor *fact5 = Divide(fact1, fact2, 2);
	printFactor(fact5);

	vector <long> v, b;
	v.push_back(1);
	v.push_back(2);
	v.push_back(3);
	b.push_back(2);
	b.push_back(1);
	Factor* dd = sumOutVariables(fact5, v, b, 2);
	//Factor* dd = reduction(fact1, 2, 1, 2);
	printFactor(dd);
	return 0;
}

*/
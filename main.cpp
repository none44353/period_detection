#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <time.h>

using namespace std;


#include "ssummary.h"
#include "BF.h"

string datapath[5] = {"formatted00.dat",
                "formatted01.dat",
                "formatted02.dat",
                "formatted03.dat",
                "formatted04.dat"};
ifstream fin;
/*FILE *fin[5] = {fopen("formatted00.dat", "rb"), 
				fopen("formatted01.dat", "rb"), 
				fopen("formatted02.dat", "rb"), 
				fopen("formatted03.dat", "rb"), 
				fopen("formatted04.dat", "rb")};*/

const int M = 1e8;

pair <uint64_t, double> Read()
{   
    static bool isfirstRead = true;
    static int curFinID = 0;
    static double offset = 0;
    static double lastT = 0;

	double t; uint64_t s;
    if (isfirstRead) {
        isfirstRead = false;
        fin.open(datapath[curFinID], std :: ios :: binary);
    }

	if (fin.eof()) {
        fin.close();
        fin.open(datapath[++curFinID], std :: ios :: binary);
    }
    fin.read((char*)&t, sizeof(double));
    fin.read((char*)&s, sizeof(uint64_t));

    
    t += offset;
    if(t < lastT) {
        offset += lastT - t;
        t += (lastT - t);
    }
    lastT = t;

	return make_pair(s, t * 1e308); //假设时间戳的单位是s，现在的效果是每0.05s有1W个item流过
}

int main() {
    for (int i = 0; i < M; ++i) {
        auto e = Read();
       // printf("%llu\n", e.first);
    }
	return 0;
}
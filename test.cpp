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
#include "detector.h"

string datapath[5] = {"../../usr/share/dataset/CAIDA2018/dataset/130000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130400.dat"};

/*                
string datapath[5] = {"../formatted00.dat",
                "../formatted01.dat",
                "../formatted02.dat",
                "../formatted03.dat",
                "../formatted04.dat"};*/
ifstream fin;

const int M = 1e8;
const double unitTime = 0.01;

pair <uint64_t, double> Read()
{   
    static bool isfirstRead = true;
    static int curFinID = 0;
    static double offset = 0;
    static double lastT = 0;

	double t; uint64_t s, _s;
    if (isfirstRead) {
        isfirstRead = false;
        fin.open(datapath[curFinID], std :: ios :: binary);
    }

	if (fin.eof()) {
        cerr << "?????????" << endl;
        fin.close();
        fin.open(datapath[++curFinID], std :: ios :: binary);
    }
    fin.read((char*)&s, sizeof(uint64_t)); //srcip(4)+dstip(4) 
    fin.read((char*)&_s, 5);//srcport destport protcol

    fin.read((char*)&t, sizeof(double));

    
    t += offset;
    if(t < lastT) {
        offset += lastT - t;
        t += (lastT - t);
    }
    lastT = t;

  /*  t *= 1e100;
    t *= 1e100;
    t *= 1e100;
    t *= 1e8;*/
	return make_pair(s, t); //假设时间戳的单位是s，现在的效果是每0.05s有1W个item流过
}

map <uint64_t, double> startTime_predict;
map <uint64_t, double> endTime_predict;
int mem;

int main() {//包含元素个数
    Alg* durationTimeDetector = new Alg(12800, 10, 0.6); //L * 24 / 1024 = 300KB
    Alg* intervalDetector = new Alg(12800, 10, 0.6);  //L * 24 / 1024 = 300KB
    bloomfliter* timeStamp = new bloomfliter(2, 64000); //t * M * 8 / 1024 = 100KB
    bloomfliter* startTime = new bloomfliter(2, 64000); //t * M * 8 / 1024 = 100KB

    for (int i = 0; i < M; ++i) {
        auto e = Read(); 
        
        uint64_t id = e.first;
        double curTime = e.second;
        
        if (i && i % 500000 == 0) {
            cerr << i << ' ' << setiosflags(ios::fixed) << curTime << ' '<< id << endl;
        }
    }
    
    fin.close();
	return 0;
}
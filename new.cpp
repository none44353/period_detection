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
#include <set>
#include <string>
#include <time.h>

using namespace std;

#define Range 500
#define delta 0.5

#include "ssummary.h"
#include "BF.h"
#include "detector.h"
#include "detectorHG.h"

string datapath[60] = {"./130000.dat"};

/*
string datapath[60] = {"../../usr/share/dataset/CAIDA2018/dataset/130000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130600.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130900.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131600.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131900.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132600.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132900.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133600.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133900.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134600.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134900.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135600.dat", 
                "../../usr/share/dataset/CAIDA2018/dataset/135700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135900.dat"};*/

ifstream fin;

const int M = 2e5;
const double P = 0.8;
const double unitTime = 0.001;

pair <uint64_t, double> Read()//新CAIDA
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
        fin.close();
        fin.open(datapath[++curFinID], std :: ios :: binary);
        if (curFinID > 60) {
            fin.close();
            exit(0);
        }
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

	return make_pair(s, t * 2); 
}

pair <uint64_t, double> input[M + 7];

bool check_near(double c, double x) {
    if (fabs(x - c) > Range * 2) return false;
    return fabs(x - c) <= delta * max(fabs(c), 0.001);
}

struct getGT{
    map <uint64_t, int> count;
    map <uint64_t, set <double> > table;
    
    void init() {
        count.clear();
        table.clear();
    }

    void insert(uint64_t id, double key) {
        if (count.find(id) == count.end()) 
            count[id] = 0, table[id].clear();
        count[id]++, table[id].insert(key);
    }

    pair <bool, double> query(const uint64_t& id) {
        if (count.find(id) == count.end()) return make_pair(false, -1);
        int n = count[id]; 
        if (n < 10) return make_pair(false, -1);
        
        int x = n * P;
        set <double> :: iterator itl = table[id].begin();
        set <double> :: iterator itr = table[id].begin();

        for (int i = 1; i < x; ++i) ++itr;
        for (; itr != table[id].end(); ++itl, ++itr) {
            double xl = *itl;
            double xr = *itr;
            double c = xl / (1 - delta);
          //  printf("check [%.6lf %.6lf] %.6lf\n", xl, xr, c);
            if (c > 0.001) {
                if (c * (1 + delta) >= xr) return make_pair(true, c);
            }
            else {
                c = xl + 0.001 * delta;
                if (c + 0.001 * delta >= xr) return make_pair(true, c);
            }
        }

        return make_pair(false, -1);
    }

    double queryPercentage(const uint64_t& id, double c) {
        set <double> :: iterator it = table[id].begin();
        int cnt = 0;
        for (; it != table[id].end(); ++it)
            if (check_near(c, *it)) ++cnt;
        return (double)cnt / count[id];     
    }

}intervalGT;


map <uint64_t, double> timeStamp;
map <uint64_t, double> startTime;
map <uint64_t, int> counter;
map <uint64_t, pair <bool, double> > intervalAnswer;
map <uint64_t, bool> check;

void GroundTruth() {
    for (int i = 0; i < M; ++i) {
        auto e = input[i]; 

        uint64_t id = e.first;
       // if (id != 17213472216328176400ULL) continue;
        double curTime = e.second;
        double lastTime = 0;
        
        if (timeStamp.find(id) != timeStamp.end()) lastTime = timeStamp[id];
        if (curTime - lastTime <= unitTime); //连续
        else {//新的预测
            double st = 0;
            if (startTime.find(id) != startTime.end()) st = startTime[id]; 

            if (counter[id] >= 5) intervalGT.insert(id, curTime - st);
            
            counter[id] = 0;
            startTime[id] = curTime;
        }
        
        if (counter.find(id) == counter.end()) counter[id] = 0;
        ++counter[id];
        timeStamp[id] = curTime;
    }

    intervalAnswer.clear();
    for (int i = 0; i < M; ++i){
        auto e = input[i]; 
        uint64_t id = e.first;
        
       // if (id != 17019272071912997107ULL) continue;
      //  if (id != 17213472216328176400ULL) continue;
        if (intervalAnswer.find(id) != intervalAnswer.end()) continue;
        pair <bool, double> result = intervalGT.query(id);
        intervalAnswer[id] = result;
      //  printf("c %.6lf per %.6lf\n", result.second, intervalGT.queryPercentage(id, result.second));
    }
}

long double precision[2], recall[2], centre_precision[2];
int cnt;


void OURS() {
    bloomfliter* ntimeStamp = new bloomfliter(4, 32768); //t * M * 8 / 1024 = 1024KB
    bloomfliter* nstartTime = new bloomfliter(4, 32768); //t * M * 8 / 1024 = 1024KB
    Alg* intervalDetector = new Alg(6400, 10, P);  //L * 24 / 1024 = 150KB
    AlgHG* HGDetector = new AlgHG(4, 200000, 10, 0.4);  //L * 24 / 1024 = 150KB

    for (int i = 0; i < M; ++i) {
        auto e = input[i]; 

        uint64_t id = e.first;
        //if (id != 2051345630196722313ULL) continue;
       // if (id != 17213472216328176400ULL) continue;
      //  if (id != 17019272071912997107ULL) continue;
        double curTime = e.second;
        double lastTime = ntimeStamp -> query(id);

        if (curTime - lastTime <= unitTime); //连续
        else {//新的预测
            double st = nstartTime -> query(id);

            if (counter[id] >= 5) intervalDetector -> insert(id, curTime - st);
           // if (counter[id] >= 5) printf("#%.6lf\n", curTime - st);
            if (counter[id] >= 5) HGDetector -> insert(id, curTime - st);
            counter[id] = 0;
            nstartTime -> insert(id, curTime);
        }
        
        if (counter.find(id) == counter.end()) counter[id] = 0;
        ++counter[id];
        ntimeStamp -> insert(id, curTime);
    }

/*
17019272071912997107
14973527361033788018
14553577667936389904*/
    check.clear();
    for (int i = 0; i < M; ++i) {
        auto e = input[i]; 
        uint64_t id = e.first;
       // if (id != 17019272071912997107ULL) continue;
        
        if (check.find(id) != check.end()) continue;
        check[id] = true; ++cnt;
        pair <bool, double> guess = HGDetector -> query(id);
        pair <bool, double> result = intervalAnswer[id];
        if (guess.first) {
            precision[1]++;
          //  printf("%llu\n", id);
            if (result.first) precision[0]++;           
        }

        if (result.first) {
            recall[1]++;
            if (guess.first) recall[0]++;
           // else printf("%llu\n", id);
        }

        if (guess.first && result.first) {
            double p = intervalGT.queryPercentage(id, guess.second);
            centre_precision[1]++;
            centre_precision[0] += (p >= P ? 1 : 0);
        }
    }
}

//17213472216328176400 
//1006749987763913488 
int main() {
    for (int i = 0; i < M + 1; ++i) input[i] = Read();
    
    GroundTruth(); puts("calc GT");
    OURS(); puts("calc OURS");
    printf("%d\n", cnt);
    printf("precision %.lf %.6lf\n", (double)precision[1], (double)(precision[0] / precision[1]));
    printf("recall %.lf %.6lf\n", (double)recall[1], (double)(recall[0] / recall[1]));
    printf("centre_precision %.lf %.6lf\n", (double)centre_precision[1], (double)(centre_precision[0] / centre_precision[1]));

	return 0;
}
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

//string datapath[60] = {"./130000.dat"};


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
                "../../usr/share/dataset/CAIDA2018/dataset/135900.dat"};

ifstream fin;

const int M = 2e7;

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
    return fabs(x - c) < delta * max(fabs(c), 0.001);
}

const double threshold = 1;
const int P = 10;

struct getGT{
    map <uint64_t, double> sum;
    map <uint64_t, double> sum2;
    map <uint64_t, int> count;
    map <uint64_t, set <double> > table;
    void init() {
        sum.clear();
        sum2.clear();
        count.clear();
        table.clear();
    }
    void insert(uint64_t id, double key) {
        if (sum.find(id) == sum.end()) 
            sum[id] = sum2[id] = count[id] = 0, table[id].clear();
        sum[id] += key, sum2[id] += key * key, count[id]++, table[id].insert(key);
    }

    pair <bool, double> query(const uint64_t& id) {
        if (sum.find(id) == sum.end()) return make_pair(false, -1);
        int n = count[id]; 
        if (n < 10) return make_pair(false, -1);
        int x = n / P; //去除最大最小P%元素 
      //  x = 0;
        double s = sum[id], s2 = sum2[id];
  
        set <double> :: iterator it = table[id].begin();
        for (int i = 0; i < x; ++i) 
            --n, s -= *it, s2 -= (*it) * (*it), ++it;
        set <double> :: reverse_iterator rit = table[id].rbegin();
        for (int i = 0; i < x; ++i)
            --n, s -= *rit,s2 -= (*it) * (*it), ++rit;

        double mn = *it, mx = *rit;
        
        double mean = s / n;
        /*
        double std2 = (double)s2 / n - mean * mean;
        double std = sqrt(std2);
        double res = std / mean;*/

       // printf("query %llu %d %.5lf %.5lf\n", id, n, mean, res);

        //if (res <= threshold)
          //  return make_pair(true, mean);

        if (check_near(mean, mn) && check_near(mean, mx))
            return make_pair(true, mean); 
        return make_pair(false, mean);
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
        double curTime = e.second;
        double lastTime = 0;
        
        if (timeStamp.find(id) != timeStamp.end()) lastTime = timeStamp[id];
        if (curTime - lastTime <= unitTime); //连续
        else {//新的预测
            double st = 0;
            if (startTime.find(id) != startTime.end()) st = startTime[id]; 

            if (counter[id] >= 5) {
                intervalGT.insert(id, curTime - st);
               // printf("insert %llu %.4lf\n", id, curTime - st);
            }
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
        
        if (intervalAnswer.find(id) != intervalAnswer.end()) continue;
        pair <bool, double> result = intervalGT.query(id);
        intervalAnswer[id] = result;
    }
}

long double precision[2], recall[2], are[2];
int cnt;


void OURS() {
    bloomfliter* ntimeStamp = new bloomfliter(4, 32768); //t * M * 8 / 1024 = 1024KB
    bloomfliter* nstartTime = new bloomfliter(4, 32768); //t * M * 8 / 1024 = 1024KB
    Alg* intervalDetector = new Alg(6400, 10, 0.2);  //L * 24 / 1024 = 150KB
    AlgHG* HGDetector = new AlgHG(4, 1600, 10, 0.2);  //L * 24 / 1024 = 150KB

    for (int i = 0; i < M; ++i) {
        auto e = input[i]; 

        uint64_t id = e.first;
        //if (id != 17213472216328176400ULL) continue;
        double curTime = e.second;
        double lastTime = ntimeStamp -> query(id);

        if (curTime - lastTime <= unitTime); //连续
        else {//新的预测
            double st = nstartTime -> query(id);

            if (counter[id] >= 5) intervalDetector -> insert(id, curTime - st);
            if (counter[id] >= 5) HGDetector -> insert(id, curTime - st);
            counter[id] = 0;
            nstartTime -> insert(id, curTime);
        }
        
        if (counter.find(id) == counter.end()) counter[id] = 0;
        ++counter[id];
        ntimeStamp -> insert(id, curTime);
    }

    check.clear();
    int Mp = 0, Lp = 0;
    int pre = 0, are_n = 0;
    double are_sum = 0;
    for (int i = 0; i < M; ++i){
        auto e = input[i]; 
        uint64_t id = e.first;
        
        if (check.find(id) != check.end()) continue;
        check[id] = true; ++cnt;
        pair <bool, double> guess = HGDetector -> query(id);
        pair <bool, double> result = intervalAnswer[id];
        if (guess.first) {
            precision[1]++;
            if (result.first) precision[0]++;           
        }

        if (result.first) {
            recall[1]++;
            if (guess.first) recall[0]++;
        }

        if (guess.first && result.first) {
            are[1]++;
            are[0] += fabs(guess.second - result.second) / result.second;
        }

        if (!guess.first) continue;
        double p_ours = intervalGT.queryPercentage(id, guess.second);
        double p_GT = intervalGT.queryPercentage(id, result.second);
        if (p_ours >= p_GT) ++Mp;
        else ++Lp;

        if (p_ours >= 0.4) ++pre;
        else ++are_n, are_sum += (0.4 - p_ours) / 0.4;

    }
    printf("ours vs mean: %d %d %.6lf\n", Mp, Lp, (double)Mp / (Mp + Lp));
    printf("ours: precision=%.6lf are=%.6lf %d\n", (double)pre / (Mp + Lp), (double)are_sum / are_n, are_n);
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
    printf("are %.lf %.6lf\n", (double)are[1], (double)(are[0] / are[1]));

	return 0;
}
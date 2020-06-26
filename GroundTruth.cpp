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

#include "ssummary.h"
#include "BF.h"
#include "detector.h"

/*
string datapath[50] = {"../DataSet/CAIDA/formatted00.dat",
                "../DataSet/CAIDA/formatted01.dat",
                "../DataSet/CAIDA/formatted02.dat",
                "../DataSet/CAIDA/formatted03.dat",
                "../DataSet/CAIDA/formatted04.dat"}; //旧CAID*/

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

//const double unitTime = 0.01;
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


/*pair <uint64_t, double> Read() //旧CAIDA
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

    t *= 1e308;
	return make_pair(s, t); //假设时间戳的单位是s，现在的效果是每1s有50W个item流过
}*/
#define delta 0.01
bool check_near(double c, double x) {
    return fabs(x - c) < delta;
}

struct getGT{
    map <uint64_t, double> sum;
    map <uint64_t, int> count;
    map <uint64_t, set <double> > table;
    void init() {
        sum.clear();
        count.clear();
        table.clear();
    }
    void insert(uint64_t id, double key) {
       // printf("insert %llu %.5lf\n", id, key);
        if (sum.find(id) == sum.end()) 
            sum[id] = count[id] = 0, table[id].clear();
        sum[id] += key, count[id]++, table[id].insert(key);
    }

    pair <bool, double> query(const uint64_t& id) {
        if (sum.find(id) == sum.end()) return make_pair(false, -1);
        int n = count[id]; double s = sum[id];
        if (n < 10) return make_pair(false, -1);
        int x = n / 5; //去除最大最小20%元素
        double mn, mx;
        set <double> :: iterator it = table[id].begin();
        for (int i = 0; i < x; ++i) 
            --n, s -= *it, ++it;
        mn = *it;
        set <double> :: reverse_iterator rit = table[id].rbegin();
        for (int i = 0; i < x; ++i)
            --n, s -= *rit, ++rit;
        mx = *rit;
        double mean = s / n;
        if (check_near(mean, mn) && check_near(mean, mx))
            return make_pair(true, mean);
        return make_pair(false, mean);
    }
}duration, interval;

map <uint64_t, double> timeStamp;
map <uint64_t, double> startTime;
map <uint64_t, double> startTime_predict;
map <uint64_t, double> endTime_predict;
map <uint64_t, int> counter;
map <uint64_t, int> counter_sum;
map <uint64_t, int> counter_num;
map <uint64_t, double> counter_sum2;
int n[2], precision[2];
double aae[2];

void check(double p, double x, int type) {
    if (fabs(p - x) < delta) {
        ++precision[type];
        aae[type] += fabs(p - x);
    }
}

int main() {
    for (int i = 0; i < M + 1; ++i) {
        auto e = Read(); 

        uint64_t id = e.first;
        double curTime = e.second;
        //if (id != 5941329809716782921UL) continue;
        //if (id != 10780215907711968843UL) continue;
        //printf("%d %llu %.6lf\n", i, id, curTime);
        double lastTime = 0;
        if (timeStamp.find(id) != timeStamp.end()) lastTime = timeStamp[id];
        if (curTime - lastTime <= unitTime); //连续
        else {//验证预测 、新的预测
            double st = 0;
            if (startTime.find(id) != startTime.end()) st = startTime[id]; 

       //     printf("%.6lf %d %.6lf\n", lastTime, counter[id], curTime);

            double startp = -1; //下一个周期开始时间的预测
            if (startTime_predict.find(id) != startTime_predict.end()) startp = startTime_predict[id];
            if (startp > 0) check(startp, curTime, 0);

            if (counter[id] >= 10) {
                duration.insert(id, lastTime - st); 
                interval.insert(id, curTime - st);
            }
            startTime[id] = curTime;

           // if (i >= 50000) {
                //新的预测
                //关于结束时间的预测
                if (counter[id] >= 10) {
                    //关于开始时间的预测
                    auto it = interval.query(id);
                    if (it.first) {
                        ++n[0];
                        startp = curTime + it.second; 
                        startTime_predict[id] = startp;
                    } 
                    else 
                        if (startTime_predict.find(id) != startTime_predict.end()) 
                        startTime_predict.erase(startTime_predict.find(id));
                }
                else 
                    if (startTime_predict.find(id) != startTime_predict.end()) 
                    startTime_predict.erase(startTime_predict.find(id));
            //}

            if (counter[id] != 0) {
                if (counter_num.find(id) == counter_num.end()) counter_num[id] = counter_sum[id] = counter_sum2[id] = 0;

                counter_num[id]++, counter_sum[id] += counter[id], counter_sum2[id] += (double)counter[id] * counter[id];
            }

            counter[id] = 0;
        }
        
        if (counter.find(id) == counter.end()) counter[id] = 0;
        ++counter[id];
        timeStamp[id] = curTime;

        if (i && i % 50000 == 0) {
          //  cerr << i << ' ' << setiosflags(ios::fixed) << curTime << ' '<< id << endl;
            printf("M = %d\n", i);
            cerr << i << endl;
            for (int k = 0; k < 1; ++k) {
                printf("n=%d precision=%.6lf aae=%.6lf\n", n[k], (double)precision[k] / n[k], aae[k] / precision[k]);
                cerr << "n=" << n[k] << setprecision(6) << setiosflags(ios::fixed)
                     << " precision=" << (double)precision[k] / n[k] 
                     << " aae=" << aae[k] / precision[k] << endl;
            }
        }
    }

    freopen("statistic.csv", "w", stdout);
    auto its2 = counter_sum2.begin();
    for (auto it = counter_num.begin(), its = counter_sum.begin(); it != counter_num.end(); ++it, ++its, ++its2) {
        double mean = (double)its -> second / it -> second;
        int n = it -> second;
        double std2 = (double)its2 -> second / n - mean * mean;
        double std = sqrt(std2);

        printf("%llu, %.2lf, %.2lf, %.2lf, %d\n", it -> first, mean, std, std / mean, n);
    }

    fin.close();
	return 0;
}
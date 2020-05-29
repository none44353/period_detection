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
                
string datapath[60] = {"/usr/share/dataset/CAIDA2018/dataset/130000.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130100.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130200.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130300.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130400.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130500.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130600.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130700.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130800.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130900.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131000.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131100.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131200.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131300.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131400.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131500.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131600.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131700.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131800.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131900.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132000.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132100.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132200.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132300.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132400.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132500.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132600.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132700.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132800.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132900.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133000.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133100.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133200.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133300.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133400.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133500.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133600.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133700.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133800.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133900.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134000.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134100.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134200.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134300.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134400.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134500.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134600.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134700.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134800.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134900.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135000.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135100.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135200.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135300.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135400.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135500.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135600.dat", 
                "/usr/share/dataset/CAIDA2018/dataset/135700.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135800.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135900.dat"};

ifstream fin;

const int M = 2e9;
const double unitTime = 0.01;

/*
string datapath[50] = {"../formatted00.dat",
                "../formatted01.dat",
                "../formatted02.dat",
                "../formatted03.dat",
                "../formatted04.dat"};*/

/*pair <uint64_t, double> Read()
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

map <uint64_t, double> startTime_predict;
map <uint64_t, double> endTime_predict;
map <uint64_t, int> counter;

int n[2], precision[2];
double aae[2];

void check(double p, double x, int type) {
    if (fabs(p - x) < delta) {
        ++precision[type];
        aae[type] += fabs(p - x);
    }
}

int mem;

int main() {//包含元素个数
    freopen("result.txt", "w", stdout);
    cerr << "150KB" << endl;
    Alg* durationTimeDetector = new Alg(6400, 10, 0.6); //L * 24 / 1024 = 150KB
    Alg* intervalDetector = new Alg(6400, 10, 0.6);  //L * 24 / 1024 = 150KB
    bloomfliter* timeStamp = new bloomfliter(2, 64000); //t * M * 8 / 1024 = 100KB
    bloomfliter* startTime = new bloomfliter(2, 64000); //t * M * 8 / 1024 = 100KB
        
    for (int i = 0; i < M + 1; ++i) {
        auto e = Read(); 

        uint64_t id = e.first;
        double curTime = e.second;
        double lastTime = timeStamp -> query(id);
        if (curTime - lastTime <= unitTime); //连续
        else {//验证预测 、新的预测
            double st = startTime -> query(id); 

            //验证预测
            double endp = -1; //上一个周期结束时间的预测
            if (endTime_predict.find(id) != endTime_predict.end()) endp = endTime_predict[id];
            if (endp > 0) check(endp, lastTime, 1);

            double startp = -1; //下一个周期开始时间的预测
            if (startTime_predict.find(id) != startTime_predict.end()) startp = startTime_predict[id];
            if (startp > 0) check(startp, curTime, 0);

            if (counter[id] >= 20) {
                durationTimeDetector -> insert(id, lastTime - st); 
                intervalDetector -> insert(id, curTime - st);
            }
            startTime -> insert(id, curTime);

            if (i >= 500000000) {
                //新的预测
                //关于结束时间的预测
                if (counter[id] >= 20) {
                    auto it = durationTimeDetector -> query(id);
                    if (it.first) {
                        ++n[1];
                        endp = curTime + it.second; 
                        endTime_predict[id] = endp;
                    }
                    else 
                        if (endTime_predict.find(id) != endTime_predict.end())
                            endTime_predict.erase(endTime_predict.find(id));
                    
                    //关于开始时间的预测
                    it = intervalDetector -> query(id);
                    if (it.first) {
                        ++n[0];
                        startp = curTime + it.second; 
                        startTime_predict[id] = startp;
                    } 
                    else 
                        if (startTime_predict.find(id) != startTime_predict.end()) 
                        startTime_predict.erase(startTime_predict.find(id));
                }
            }
                
            counter[id] = 0;
        }
        
        ++counter[id];
        timeStamp -> insert(id, curTime);

        if (i && i % 50000000 == 0) {
            cerr << '?' << endl;
          //  cerr << i << ' ' << setiosflags(ios::fixed) << curTime << ' '<< id << endl;
            printf("M = %d\n", i);
            cerr << i << endl;
            for (int k = 0; k < 2; ++k) {
                printf("n=%d precision=%.6lf aae=%.6lf\n", n[k], (double)precision[k] / n[k], aae[k] / precision[k]);
                cerr << "n=" << n[k] << setprecision(6) << setiosflags(ios::fixed)
                     << " precision=" << (double)precision[k] / n[k] 
                     << " aae=" << aae[k] / precision[k] << endl;
            }
        }
    }
    
    fin.close();
	return 0;
}

//
//  main.cpp
//  JSP-PSO
//
//  Created by 林理露 on 16/4/25.
//  Copyright © 2016年 林理露. All rights reserved.
//
//  使用中级元启发式算法解决JSP问题
//  使用PSO粒子群优化算法，多次迭代使其收敛于最优解
//

#include <iostream>
#include <cmath>
#include <climits>
#include <queue>
#include <ctime>
#include <cassert>
#ifndef DEBAG
#define DEBAG 0
#endif
using namespace std;
//参数
#define Vmax 1      //最高速度
#define Xmax 5      //最大范围
#define Wmax 2.2    //最大惯性
#define Wmin 0.4    //最小惯性
#define Itermax 600 //最高迭代次数
#define Pmax 50     //种群规模
#define K 5         //最高不变解跳出次数

class step{
public:
    //工序
    int j_NO;//对应工件序号
    int m_NO;//对应机器序号
    int duration;//完成工序所需时间
    step():j_NO(0),m_NO(0),duration(0){}
    step(int j_NO,int m_NO,int t):j_NO(j_NO),m_NO(m_NO),duration(t){}
};

class pairs{
public:
    //用于把连续粒子离散化
    double part_x;//粒子位置
    int idx;//对应下标
    pairs():part_x(-1),idx(-1){}
    pairs(double part_x,int idx):part_x(part_x),idx(idx){}
    
    friend bool operator<(const pairs &a,const pairs &b){
        return (a.part_x > b.part_x);
    }
};

class particle{
public:
    //粒子
    int n,p,size;//粒子内部信息容量
    int value;//粒子的价值(对应make_span的长短)
    double* x;//位置
    double* v;//速度
    double* Pbest;//个人所到达过的佳点
    int Pb_value;//个人最佳点价值
    int* operations;//用于存储解码出的规划
    
    particle():n(0),p(0),size(0),x(nullptr),v(nullptr),\
    Pbest(nullptr),value(INT_MAX),Pb_value(INT_MAX){}
    
    particle(int n,int p):n(n),p(p),size(n*p),value(INT_MAX),Pb_value(INT_MAX){
        x = new double[size];
        v = new double[size];
        Pbest = new double[size];
        operations = new int[size];
    }
    
    particle(const particle& r){
        n = r.n, p = r.p, size = r.size;
        value = r.value;
        Pb_value = r.Pb_value;
        x = new double[size];
        v = new double[size];
        Pbest = new double[size];
        operations = new int[size];
        for (int i = 0; i < size; ++i) {
            x[i] = r.x[i];
            v[i] = r.v[i];
            Pbest[i] = r.Pbest[i];
            operations[i] = r.operations[i];
        }
    }
    
    void init(){
        //粒子初始化
        for (int i = 0; i < size; ++i) {//随机生成
            x[i] = 2*Xmax*rand()/double(RAND_MAX)-Xmax;
            v[i] = 2*Vmax*rand()/double(RAND_MAX)-Vmax;
        }
    }
    
    particle& operator=(const particle& r){
        n = r.n,p = r.p,size = r.size;
        value = r.value;
        Pb_value = r.Pb_value;
        x = new double[size];
        v = new double[size];
        Pbest = new double[size];
        operations = new int[size];
        for (int i = 0; i < size; ++i) {
            x[i] = r.x[i];
            v[i] = r.v[i];
            Pbest[i] = r.Pbest[i];
            operations[i] = r.operations[i];
        }
        return *this;
    }
    
    bool operator<(const particle& r){
        return (value < r.value);
    }
    
    void ROV() {
        //ROV规则 离散化(解码)
        priority_queue<pairs> q;
        for (int i = 0; i < p; ++i) {
            for (int j = 0; j < n; ++j) {
                q.push(pairs(x[i*n+j],j));
            }
            for (int j = 0; j < n; ++j) {
                operations[i*n+q.top().idx] = j;
                q.pop();
            }
        }
    }
    
    ~particle(){
//        if (x != nullptr)
            delete[] x;
//        if (v != nullptr)
            delete[] v;
//        if (Pbest != nullptr)
            delete[] Pbest;
//        if (operations != nullptr)
            delete[] operations;
    }
        
};


int main(int argc, const char * argv[]) {
    //主函数
//    freopen("/users/LLL/Desktop/la12.in", "r", stdin);
    //基础输入
    clock_t t1 = clock();
    int n = 0,m = 0,p = 0,size = 0,size2 = 0;//n个工件,m台机器,p道工序,矩阵规模size
    char tp;
    cin>>n>>tp>>m>>tp>>p>>tp;
    size = n*p,size2 = n*m;
    srand (int(time(NULL)));
    //内存分配
    step** steps = new step*[n+1];
    particle *parts = new particle[Pmax];
    int* Tm = new int[size2];//代表每台机器上的累计加工时间
    int* Tp = new int[size];//代表每个工件上的累计加工时间
    int* idx = new int[n];//代表每个工件当前进行到的工序数
    
    //工序、时间 信息输入
    for (int i = 0; i< n; ++i) {
        steps[i] = new step[p];//顺带分配内存
        for (int j = 0 ; j < p; ++j) {
            cin>>steps[i][j].duration;
        }
    }
    
    if (DEBAG) {
        cout<<"duration"<<endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < p; ++j) {
                cout<<steps[i][j].duration<<" ";
            }
            cout<<endl;
        }
    }
    
    //工序、机器 信息输入
    for (int i = 0; i< n; ++i) {
        for (int j = 0 ; j < p; ++j) {
            cin>>steps[i][j].m_NO;
            steps[i][j].j_NO = i;//顺带初始化其它
        }
    }
    
    if (DEBAG) {
        cout<<"machines:"<<endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < p; ++j) {
                cout<<steps[i][j].m_NO<<" ";
            }
            cout<<endl;
        }
    }
    
    //参数初始化
    double w = Wmin,c1 = 2,c2 = 2;
    particle Gbest(n,p);
    double make_span_old = 0,make_span_new = 0;
    
    //粒子群初始化
    for (int i = 0; i < Pmax; ++i) {
        parts[i] = particle(n,p);
        parts[i].init();
    }
    
//    if (DEBAG) {
//        cout<<"after init:"<<endl;
//        for (int i = 0; i < Pmax; ++i) {
//            cout<<"n: "<<parts[i].n<<" p: "<<parts[i].p<<" size: "<<parts[i].size<<endl;
//            cout<<"x:"<<endl;
//            for (int j = 0; j < size; ++j) {
//                cout<<parts[i].x[j]<<" ";
//            }
//            cout<<endl<<"v:"<<endl;
//            for (int j = 0; j < size; ++j) {
//                cout<<parts[i].v[j]<<" ";
//            }
//            cout<<endl;
//        }
//    }
    
    int iter_times = 0,ctr = 0;
    bool updated = false;
    for (iter_times = 0; iter_times < Itermax; ++iter_times) {
        updated = false;//用于判断解是否更新过
        particle* min = &parts[0];
        for (int i = 0; i < Pmax; ++i) {
            //解码
            parts[i].ROV();//离散化
            if (DEBAG) {
                cout<<"After ROV"<<endl;
                for (int j = 0; j < size; ++j) {
                    cout<<parts[i].operations[j]<<" ";
                }
                cout<<endl;
            }
            
            //初始化
            for (int j = 0; j < size2; ++j) {Tm[j] = 0;}
            for (int j = 0 ; j < size; ++j) {Tp[j] = 0;}
            for (int j = 0 ; j < n; ++j) {idx[j] = 0;}
            int max = INT_MIN;
            //计算目标函数值
            for (int j = 0; j < size; ++j) {
                int NO = parts[i].operations[j];
                int T = steps[NO][idx[NO]].duration;
                int R = (Tp[NO] > Tm[steps[NO][idx[NO]].m_NO-1] ? \
                Tp[NO] : Tm[steps[NO][idx[NO]].m_NO-1]);
                Tp[NO] = R+T;
                Tm[steps[NO][idx[NO]].m_NO-1] = R+T;
                ++idx[NO];
            }
            for (int j = 0 ; j < size2; ++j) {
                if (Tm[j] > max) {
                    max = Tm[j];
                }
            }
            if (DEBAG) {
                cout<<"max: "<<max<<endl;
            }
            parts[i].value = max;
            assert(max > 0);
            //搜索局部最优
            if (parts[i].value < parts[i].Pb_value) {
                parts[i].Pb_value = parts[i].value;
                for (int j = 0; j < size; ++j) {
                    parts[i].Pbest[j] = parts[i].x[j];
                }
            }
            //搜索全局最优
            if (parts[i].value < (*min).value) {
                min = &parts[i];
            }
        }
        
        //更新全局最优
       
        if (min->value < Gbest.value) {
            if (DEBAG) {
                cout<<"updated:"<<Gbest.value<<"->"<<min->value<<endl;
            }
            updated = true;
            Gbest = *min;
        }else{
            if (DEBAG) {
                cout<<"abandoned:"<<"x "<<min->value<<endl;
            }
        }
        
        if (DEBAG) {
            cout<<"Values:"<<endl;
            for (int i = 0; i < Pmax; ++i) {
                assert(parts[i].value > 0);
                cout<<parts[i].value<<" ";
            }
            cout<<endl;
        }
        if (DEBAG) {
            cout<<"Gbest:"<<endl;
            for (int i = 0; i < size; ++i) {
                cout<<Gbest.x[i]<<" ";
            }
            cout<<"Gbest value: "<<Gbest.value<<endl;
            cout<<endl;
        }
        
        //粒子迭代
        for (int k = 0; k < Pmax; ++k) {
            for (int i = 0; i < size; ++i) {
                parts[k].v[i] = w*parts[k].v[i] + \
                c1*rand()/double(RAND_MAX)*(parts[k].Pbest[i]-parts[k].x[i])+\
                c2*rand()/double(RAND_MAX)*(Gbest.x[i]-parts[k].x[i]);
                parts[k].x[i] += parts[k].v[i];
                //限制在范围内 引入随机化因素
                if (fabs(parts[k].v[i]) > Vmax) {
                    parts[k].v[i] = (1-0.1*rand()/double(RAND_MAX))* \
                    (parts[k].v[i]>0 ? Vmax : -Vmax);
                }
                if (fabs(parts[k].x[i]) > Xmax) {
                    parts[k].x[i] = (1-0.1*rand()/double(RAND_MAX))* \
                    (parts[k].x[i]>0 ? Xmax : -Xmax);
                }
            }
        }
        
        
        //动态参数更新
        w = Wmin + (Wmax - Wmin)/double(n) * iter_times;
        
        make_span_new = Gbest.value;
        
        //判定是否继续迭代(连续K次解都不变则跳出)
        if (updated && fabs(make_span_new - make_span_old) < 0.000000001 && ++ctr > K) {
            break;
        }
        make_span_old = make_span_new;
    }
    
    //输出
//    if (DEBAG) {
        cout<<"iterate times: "<<iter_times<<endl;
//    }
    clock_t t2 = clock();
    float diff = ((float)(t2 - t1) / 1000000.0F );
    printf("%f\n",diff);
    if (DEBAG) {
        cout<<"new schdule: "<<endl;
        for (int j = 0; j < size; ++j) {
            cout<<Gbest.operations[j]<<" ";
        }
        cout<<endl;
    }
    cout<<make_span_new<<endl;
    
    if (DEBAG) {
        cout<<"Pbests"<<endl;
        for (int i = 0; i < Pmax ; ++i) {
            cout<<parts[i].Pb_value<<" ";
        }
        cout<<endl;
    }
    
    //尾处理
    delete[] parts;
    for (int i=0; i<n; ++i) {
        delete[] steps[i];
    }
    return 0;
}

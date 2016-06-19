//
//  main.cpp
//  JSP-HHA
//
//  Created by 林理露 on 16/5/28.
//  Copyright © 2016年 林理露. All rights reserved.
//
//  使用高级超启发式算法解决JSP问题
//  用高层GA框架选择低阶启发式规则，从而寻找最优解
//

#include <iostream>
#include <climits>
#include <cmath>
#include <ctime>
#include <vector>
#include <list>
#include <cassert>
#ifndef DEBAG
#define DEBAG 0
#endif
using namespace std;
//参数
#define Pc 0.65     //基因重组概率
#define Pm 0.005    //基因突变概率
#define Gmax 600    //最多代数
#define K 400       //最多解不变跳出次数
#define Cmax 50     //染色体数量
#define Hctr 6      //启发式规则数量
#define EV 15       //环境变量 越高环境越好 种群生存的概率越大
class step;

//启发式规则列表
//用于判断 工序s1 是否优先于 工序s2
bool FA(step* s1,step* s2);           //先到先得
bool EFT(step* s1,step* s2);          //最早完工时间
bool LU(step* s1,step* s2);           //最低使用率优先
bool MA(step* s1,step* s2);           //最高空闲率
bool SPT(step* s1,step* s2);          //最短加工时间
bool TIS(step* s1,step* s2);          //在机器中所处时间最长

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
    //用于辅助同一机器上工序排序
    step* stp;//存储该工序的地址
    int idx;//存储该工序在对应工件中的位置
    pairs():stp(nullptr),idx(INT_MAX){}
    pairs(step* st,int idx):stp(st),idx(idx){}
    pairs(const pairs& pr){
        stp = pr.stp;
        idx = pr.idx;
    }
    
    pairs& operator=(const pairs& pr){
        stp = pr.stp;
        idx = pr.idx;
        return *this;
    }
    
    friend bool operator==(const pairs& p1,const pairs& p2){
        return (p1.stp == p2.stp);
    }
};

class chromosome{
public:
    //染色体
    int n,p;//染色体容量信息
    double fitness;//适应度(= 1/(obj(g)+1))
    int ms;//工作时长
    int* DNA;//存储对应的基因编码
    //用于存储可能会用到的启发式规则
    bool (*hf[Hctr])(step*,step*)={FA,EFT,LU,MA,SPT,TIS};
    chromosome():n(0),p(0),fitness(0),ms(0),DNA(nullptr){}
    chromosome(int n,int p):n(n),p(p),ms(0),fitness(0){
        DNA = new int[p];
    }
    chromosome(const chromosome& cr){
        assert(cr.n > 0),assert(cr.p > 0);
        n = cr.n, p = cr.p;
        fitness = cr.fitness,ms = cr.ms;
        DNA = new int[p];
        for (int i = 0; i < p; ++i) {
            DNA[i] = cr.DNA[i];
        }
    }
    
    void init(int n,int p){
        //随机分配算法
        assert(n > 0),assert(p > 0);
        this->n = n, this->p = p;
        DNA = new int[p];
        for (int i = 0; i < p; ++i) {
            DNA[i] = rand()%Hctr;
        }
    }
    
    friend void cross_over(chromosome c1,chromosome c2){
        //交叉操作
        if (DEBAG) {
            cout<<"p"<<c1.p<<" "<<c2.p<<endl;
        }
        int a = rand()%(c1.p), b = rand()%(c2.p);
        int begin = a < b ? a : b,end = a > b ? a : b;
        int temp = 0;
        for (int i = begin; i <= end; ++i) {
            temp = c1.DNA[i];
            c1.DNA[i] = c2.DNA[i];
            c2.DNA[i] = temp;
        }
    }
    
    chromosome& operator=(const chromosome& cr){
        assert(cr.n > 0),assert(cr.p > 0);
        n = cr.n, p = cr.p;
        fitness = cr.fitness,ms = cr.ms;
        DNA = new int[p];
        for (int i = 0; i < p; ++i) {
            DNA[i] = cr.DNA[i];
        }
        return *this;
    }
    
    void mutate(double avg){
        //变异操作
        for (int i = 0; i < p; ++i) {
            if (rand()/double(RAND_MAX) < Pm*(1+EV*(fitness-avg))){
                //适应度越差 变异率越高
                DNA[i] = rand()%Hctr;
            }
        }
    }
    
    ~chromosome(){
        delete[] DNA;
    }
};

//多个全局变量 方便启发式规则访问
int *Tm = nullptr,*Tp = nullptr,*idxs = nullptr;
step** steps = nullptr;
list<pairs>* same_machine = nullptr;
chromosome* chromo = nullptr;
int *sum_j_time = nullptr;
int main(int argc, const char * argv[]) {
    //主函数
//    freopen("/users/LLL/Desktop/ft10.in", "r", stdin);
//    freopen("/users/LLL/Desktop/ft06.out", "w", stdout);
    //基础输入
    int n = 0,m = 0,p = 0;
    clock_t t1 = clock();
    char tp;
    cin>>n>>tp>>m>>tp>>p>>tp;
    srand(int(time(NULL)));
    
    //内存分配
    steps = new step*[n];//用于存储工序信息
    same_machine = new list<pairs>[m];//用于存储同一机器上的工序
    chromo = new chromosome[Cmax];//存储染色体
    sum_j_time = new int[n];//存储对应工件所有工序所用总时间
    
    //工序、时间 信息输入
    for (int i = 0; i< n; ++i) {
        int sum_j = 0;
        steps[i] = new step[p];//顺带分配内存
        for (int j = 0 ; j < p; ++j) {
            cin>>steps[i][j].duration;
            sum_j += steps[i][j].duration;
        }
        sum_j_time[i] = sum_j;
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
        for (int j = 0 ; j < m; ++j) {
            cin>>steps[i][j].m_NO;
            steps[i][j].j_NO = i;//顺带初始化其它
            same_machine[steps[i][j].m_NO-1].\
            push_back(pairs(&steps[i][j],j));
        }
    }
    
    if (DEBAG) {
        cout<<"machines:"<<endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                cout<<steps[i][j].m_NO<<" ";
            }
            cout<<endl;
        }
    }
    
    if (DEBAG) {
        cout<<"same machine:"<<endl;
        for (int i = 0; i < p; ++i) {
            cout<<"machine "<<i<<endl;
            for(auto t : same_machine[i]){
                cout<<t.stp->duration<<'\t';
            }
            cout<<endl;
        }
    }

    //初始化染色体
    for (int i = 0; i < Cmax; ++i) {
        chromo[i].init(n,p);
    }
    
    if (DEBAG) {
        cout<<"DNAs"<<endl;
        for (int i = 0; i < Cmax; ++i) {
            for (int j = 0; j < p; ++j) {
                cout<<chromo[i].DNA[j]<<" ";
            }
            cout<<endl;
        }
    }
    //GA开始
    vector<chromosome> match;
    int iter_times = 0,ctr = 0,make_span_old = 0,make_span_new = 0;
    Tm = new int[m],Tp = new int[n];//存储机器、工件加工用时
    idxs = new int[n];//存储当前工件进行到的工序号
    bool *selected = new bool[m];//用于防止重复选择，陷入在一处
    list<pairs>* same_machine_temp = new list<pairs>[m];
    double sum = 0;
    for (iter_times = 0; iter_times < Gmax; ++iter_times) {
        //计算适应值
        for (int k = 0; k < Cmax; ++k) {
            //初始化Tm,Tp,idxs,selected,same_machine_temp
            for (int i = 0; i < m; ++i) {
                same_machine_temp[i] = same_machine[i];
            }
            if (DEBAG) {
                for (int i = 0; i < m; ++i) {
                    cout<<"same_machine_temp:"<<endl;
                    cout<<"machine "<<i<<": ";
                    for(auto s: same_machine_temp[i]){
                        cout<<s.stp->duration<<" ";
                    }
                    cout<<endl;
                }
            }
            
            for (int i = 0; i < n; ++i) {
                Tp[i] = 0;
                idxs[i] = 0;
            }
            
            for (int i = 0; i < m; ++i) {
                Tm[i] = 0;
                selected[i] = false;
            }
                
            //计算适应度
            while (1) {
                int min = INT_MAX,min_idx = -1;
                int max = INT_MIN,min2 = INT_MAX;
                //找到当前 最早、最晚 完工的机器
                for (int i = 0; i < m; ++i) {
                    if (Tm[i] < min && !selected[i]) {
                        min = Tm[i];
                        min_idx = i;
                    }
                    if (Tm[i] > max) {
                        max = Tm[i];
                    }
                }
                for (int i = 0; i < m; ++i) {
                    if (Tm[i] > min && Tm[i] < min2) {
                        min2 = Tm[i];
                    }
                }
                
                bool is_em = true;
                for (int i = 0; i < Cmax; ++i) {
                    if (!same_machine_temp[i].empty()) {
                        is_em = false;
                    }
                }
                
                if (!same_machine_temp[min_idx].empty() && min_idx != -1) {
                    if (DEBAG) {
                        cout<<"current machine: "<<min_idx+1<<endl;
                        cout<<"before:"<<endl;
                        cout<<"Tm:"<<endl;
                        for (int i = 0; i < m; ++i) {
                            cout<<Tm[i]<<" ";
                        }
                        cout<<endl;
                        
                        cout<<"Tp:"<<endl;
                        for (int i = 0; i < n; ++i) {
                            cout<<Tp[i]<<" ";
                        }
                        cout<<endl;
                        
                        cout<<"idxs:"<<endl;
                        for (int i = 0; i < n; ++i) {
                            cout<<idxs[i]<<" ";
                        }
                        cout<<endl;
                    }
                    
                    //该机器还有工序可以做
                    step* todo = nullptr;
                    for(auto &s : same_machine_temp[min_idx]){
                        if (DEBAG) {
                            cout<<"s.idx "<<s.idx<<" s.stp->j_NO "<<s.stp->j_NO<<endl;
                        }
                        if (s.idx == idxs[s.stp->j_NO] ) {
                            //是当前的前沿工序 且 根据其对应启发式规则更佳
                            if (todo == nullptr || \
                                chromo[k].hf[chromo[k].DNA[s.stp->m_NO-1]](s.stp,todo)){
                                if (DEBAG) {
                                    cout<<"selected->"<<s.stp->j_NO<<endl;
                                }
                                todo = s.stp;
                            }
                        }
                    }
                    if (todo == nullptr) {
                        //该机器剩下的工序都不是前沿工序
                        selected[min_idx] = true;
                        if (min2 == INT_MAX) {
                            min2 = min;
                        }
                        assert(min2 != INT_MAX);
                        for (int i = 0; i < m; ++i) {
                            if (i == min_idx) {
                                Tm[i] = min2;
                            }
                        }
                    }else{
                        //执行该工序
                        for (int i = 0; i < m; ++i) {
                            selected[i] = false;
                        }
                        assert(todo != nullptr);
                        int T = todo->duration;
                        int R = (Tp[todo->j_NO] > Tm[todo->m_NO-1] ? \
                                 Tp[todo->j_NO] : Tm[todo->m_NO-1]);
                        Tp[todo->j_NO] = R + T;
                        Tm[todo->m_NO-1] = R + T;
                        ++idxs[todo->j_NO];
                        //从same_machine_temp中删掉已完成的工序
                        same_machine_temp[min_idx].remove(pairs(todo,-1));
                        if (DEBAG) {
                            cout<<"list:"<<min_idx+1<<endl;
                            for(auto s : same_machine_temp[min_idx]){
                                cout<<s.stp->duration<<" ";
                            }
                            cout<<endl;
                            cout<<"list size:"<<same_machine_temp[min_idx].size()<<endl;
                        }
                    }
                    if (DEBAG) {
                        cout<<"after:"<<endl;
                        cout<<"Tm:"<<endl;
                        for (int i = 0; i < m; ++i) {
                            cout<<Tm[i]<<" ";
                        }
                        cout<<endl;
                        
                        cout<<"Tp:"<<endl;
                        for (int i = 0; i < n; ++i) {
                            cout<<Tp[i]<<" ";
                        }
                        cout<<endl;
                        
                        cout<<"idxs:"<<endl;
                        for (int i = 0; i < n; ++i) {
                            cout<<idxs[i]<<" ";
                        }
                        cout<<endl;
                    }
                }else if(!is_em && min_idx != -1){
                    //该机器已无可用工序，需要转移
                    selected[min_idx] = true;
                    if (min2 == INT_MAX) {
                        min2 = min;
                    }
                    for (int i = 0; i < m; ++i) {
                        if (i == min_idx) {
                            Tm[i] = min2;
                        }
                    }
                }else{
                    for (int i = 0; i < m; ++i) {
                        selected[i] = false;
                    }
                    //更新适应度
                    if (DEBAG) {
                        cout<<"max:"<<max<<endl;
                        cout<<"k:"<<k<<endl;
                    }
                    assert(max >= 22);
                    chromo[k].ms = max;
                    chromo[k].fitness = (1/(double(max)+1));
                    break;
                }
            }
        }
        
        //轮盘赌选入匹配集
        if(DEBAG){
            cout<<"last gen"<<endl;
            for (int i = 0; i < Cmax; ++i) {
                cout<<chromo[i].ms<<" ";
            }
            cout<<endl;
        }
        if (DEBAG) {
            cout<<"last fit"<<endl;
            for (int i = 0; i < Cmax; ++i) {
                cout<< chromo[i].fitness<<" ";
            }
            cout<<endl;
        }
        sum = 0;
        for (int i = 0; i < Cmax; ++i) {
            sum += chromo[i].fitness;
        }
        double avg_fit = sum/Cmax;//计算出平均适应度
        if (DEBAG) {
            cout<<"sum:"<<sum<<endl;
        }
        for (int i = 0; i < Cmax; ++i) {
            if (rand()/double(RAND_MAX) < \
                EV*(Gmax-iter_times*0.66)/Gmax*chromo[i].fitness/sum) {
                //随时间推移，环境变差，生存率变低
                match.push_back(chromo[i]);
            }
        }
        
        if (DEBAG) {
            cout<<"matches"<<endl;
            for(auto m : match){
                cout<<m.ms<<" ";
            }
            cout<<endl;
        }
        //变异 交叉
        for(auto &s : match){
            s.mutate(avg_fit);
            if (rand()/double(RAND_MAX) < Pc) {
                cross_over(match[rand()%match.size()],match[rand()%match.size()]);
            }
        }
        
        //死光了则结束
        if (match.size() == 0) {
            break;
        }
        //产生下一代染色体
        int j = 0;
        //老一代全部收获
        for (int j = 0; j < match.size(); ++j) {
            chromo[j] = match[j];
        }
        //空缺处随机选取补全
        if (DEBAG) {
            cout<<"match size:"<<match.size()<<endl;
        }
        for (int i = j; i < Cmax; ++i) {
            chromo[i] = match[rand()%match.size()];
        }
        
        int max = INT_MIN;
        for (int i = 0; i < Cmax; ++i) {
            if (chromo[i].ms > max) {
                max = chromo[i].ms;
            }
        }
        make_span_new = max;
        
        if (make_span_old == make_span_new && ++ctr > K) {
            break;
        }
        make_span_old = make_span_new;
        match.clear();
    }
    
    //输出
//    if (DEBAG) {
        cout<<"iter_times: "<<iter_times<<endl;
//    }
    clock_t t2 = clock();
    cout<<"time:"<<endl;
    float diff = ((float)(t2 - t1) / 1000000.0F );
    printf("%f\n",diff);
    cout<<make_span_new<<endl;
    
    //尾处理
    delete[] Tm;
    delete[] Tp;
    delete[] idxs;
    delete[] chromo;
    delete[] same_machine;
    delete[] same_machine_temp;
    for (int i=0; i<n; ++i) {
        delete[] steps[i];
    }
    return 0;
}

//启发式规则列表
//用于判断 工序s1 是否优先于 工序s2
bool FA(step* s1,step* s2){
    //先到先得
    return (Tp[s1->j_NO] < Tp[s2->j_NO]);
}

bool EFT(step* s1,step* s2){
    //最早完工时间
    return (sum_j_time[s1->j_NO] < sum_j_time[s2->j_NO]);
}

bool LU(step* s1,step* s2){
    //最低使用率优先
    return (Tp[s1->j_NO]/double(sum_j_time[s2->j_NO]) < \
            Tp[s2->j_NO]/double(sum_j_time[s2->j_NO]));
}

bool MA(step* s1,step* s2){
    //最高空闲率
    return ((sum_j_time[s1->j_NO]-s1->duration-Tp[s1->j_NO])/(sum_j_time[s1->j_NO]) > \
     (sum_j_time[s2->j_NO]-s2->duration-Tp[s2->j_NO])/(sum_j_time[s1->j_NO]));
}

bool SPT(step* s1,step* s2){
    //最短加工时间
    return (s1->duration < s2->duration);
}

bool TIS(step* s1,step* s2){
    //在机器中所处时间最长
    return ((sum_j_time[s1->j_NO]-Tp[s1->j_NO]) < \
            (sum_j_time[s2->j_NO]-Tp[s2->j_NO]));
}

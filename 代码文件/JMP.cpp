//
//  main.cpp
//  JSP
//
//  Created by 林理露 on 16/3/7.
//  Copyright © 2016年 林理露. All rights reserved.
//
//  使用初级启发式规则解决JSP问题
//  使用前沿贪心方法，即只考虑前沿工序，辅以繁忙度因素
//  使其在合法的情况下，开工时间最小
//

#include <iostream>
#include <deque>
#include <cstdlib>
#include <climits>
#ifndef DEBAG
#define DEBAG 0
#endif
using namespace std;

class step{
public:
    //工序
    int m_NO;//对应机器序号
    int j_NO;//对应工件序号
    double busy;//繁忙度
    int duration;//完成工序所需时间
    int min_t;//最早开工时间
    step():j_NO(0),m_NO(0),busy(0),duration(0),min_t(INT_MAX){}
    step(int j_NO,int m_NO,int t):j_NO(j_NO),m_NO(m_NO),busy(0),duration(t),min_t(0){}
};

class job{
public:
    //工件
    int NO;//序号
    int time_sum;//总时间
    double job_busy;//工件繁忙度
    deque<step> steps;//工件工序集合
    job():NO(0),time_sum(0),job_busy(0){}
};

class machine{
public:
    //机器
    int NO;//编号
    int time_sum;//总时间
    int to_be_handle;//未处理的工序数
    double machine_busy;//机器繁忙度
    machine():NO(0),time_sum(0),to_be_handle(0),machine_busy(0){}
};

int compare(const void* a,const void*b);
int main(int argc, const char * argv[]) {
    //主函数
//    freopen("/users/LLL/Desktop/test1.in", "r", stdin);
    //数量输入
    clock_t t1 = clock();
    
    int n=0,m=0,p=0;
    char tmp;
//    cout<<"Please enter the essential data:"<<endl;
//    cout<<"Enter the number of the job and machine"<<endl;
    cin>>n>>tmp>>m>>tmp>>p>>tmp;
    
    //初始化
    job* jobs = new job[n+1];
    machine* machines = new machine[m+1];
    int *buffer[n+1];
    for (int i=1; i<=n; ++i) {//输入缓冲区将时间先存储起来
        buffer[i] = new int[p+1];//随后同步使用
    }
    for (int i=1; i<= n; ++i) {
        for (int j = 1 ; j <= p; ++j) {
            cin>>buffer[i][j];
        }
    }
    for (int i=1; i<=n; ++i) {
        jobs[i].NO = i;
    }
    for (int i=1; i<=m; ++i) {
        machines[i].NO = i;
    }
    
    //工序、时间 信息输入
//    cout<<"Now enter steps with time."<<endl;
    for (int i=1; i<=n; ++i) {
        for (int j=1; j<=p; ++j) {
            int NO=0,t=buffer[i][j];
            cin>>NO;
            jobs[i].steps.push_back(step(i,NO,t));
            jobs[i].time_sum += t;
            machines[NO].time_sum += t;
            ++machines[NO].to_be_handle;
        }
    }
    
    if (DEBAG) {
        cout<<"job info"<<endl;
        for (int i=1; i<=n; ++i) {
            cout<<i<<" ";
            for (int j=1; j<=p; ++j) {
                step& current =jobs[i].steps[j-1];
                cout<<current.m_NO<<" "<<current.duration<<" ";
            }
            cout<<jobs[i].time_sum<<" ";
            cout<<endl;
        }
        cout<<"machine info"<<endl;
        for (int i=1; i<=m; ++i) {
            cout<<i<<" ";
            cout<<machines[i].time_sum<<" "<<machines[i].to_be_handle<<" "<<endl;
        }
    }
    
    //启发式调度算法开始
    int make_span=0;//调度总时间
    step *last_j = new step[n+1];//用于存储同工件上一个工序
    step *last_m = new step[m+1];//用于存储同机器上一个工序
    step* front_steps = new step[n+1];//前沿工序声明
    for (int i = 0; i <= n; ++i) {
        last_j[i].j_NO = i;
        last_j[i].min_t = 0;
        last_j[i].duration = 0;
    }
    for (int i = 0; i <= m; ++i) {
        last_m[i].m_NO = i;
        last_m[i].min_t = 0;
        last_m[i].duration = 0;
    }
    
    while (1) {
        //获取每个工件的前沿工序
        int is_empty = 1;
        if (DEBAG) {
            cout<<"min_t"<<endl;
        }
        for (int i=1; i<=n; ++i) {
            if (jobs[i].steps.size() > 0) {//还有剩下的工序则提出来一个
                is_empty = 0;
                int j = last_j[i].min_t+last_j[i].duration;//更新min_t
                int m = last_m[jobs[i].steps.front().m_NO].min_t+\
                        last_m[jobs[i].steps.front().m_NO].duration;
                jobs[i].steps.front().min_t = (j > m ? j : m);
                front_steps[i] = jobs[i].steps.front();
                if (DEBAG) {
                    cout<<i<<" "<<front_steps[i].min_t<<" ";
                }
            }else{
                front_steps[i] = step();//没有则置空
            }
        }
        if (DEBAG) {
            cout<<endl<<"front steps"<<endl;
            for (int i=1; i<=n; ++i) {
                cout<<front_steps[i].m_NO<<" ";
            }
            cout<<endl;
        }
        if (is_empty) {
            break;//跳出条件:若工序全部执行完则跳出
        }
    
        //繁忙度更新
        for (int i = 1; i <= n; ++i) {
            //工件繁忙度计算
            //工件繁忙度 ＝ 属于该工件的未处理的工序数 ＊ 未处理的工序所需时间和
            jobs[i].job_busy = jobs[i].steps.size() * jobs[i].time_sum;
        }
        
        for (int i = 1; i <= m; ++i) {
            //机器繁忙度计算
            //工件繁忙度 ＝ 属于该机器的的未处理的工序数 ＊ 未处理的工序所需时间和
            machines[i].machine_busy = machines[i].to_be_handle * machines[i].time_sum;
        }
        
        for (int i = 1; i <= n; ++i) {
            //前沿工序繁忙度计算
            //前沿工序繁忙度 = 工件繁忙度 + 机器繁忙度
            front_steps[i].busy = jobs[front_steps[i].j_NO].job_busy + \
            machines[front_steps[i].m_NO].machine_busy;
        }
        if (DEBAG) {
            cout<<"busy:"<<endl;
            for (int i=1; i<=n; ++i) {
                cout<<front_steps[i].busy<<" ";
            }
            cout<<endl;
        }
        
        qsort(&front_steps[1], n, sizeof(step), compare);//按给定三元组排序
        step& first = front_steps[1];
        
        if (DEBAG) {
            cout<<"sorted front steps"<<endl;
            for (int i=1; i<=n; ++i) {
                cout<<front_steps[i].m_NO<<" ";
            }
            cout<<endl;
            cout<<"first: "<<first.j_NO<<endl;
        }
        
        if (first.j_NO != 0) {//不是空的情况下
            //选中该工序执行工作
            //存储信息
            last_j[first.j_NO] = jobs[first.j_NO].steps.front();
            last_m[first.m_NO] = jobs[first.j_NO].steps.front();
            //更新时间、事务信息
            jobs[first.j_NO].steps.pop_front();
            jobs[first.j_NO].time_sum -= first.duration;
            --machines[first.m_NO].to_be_handle;
            machines[first.m_NO].time_sum -= first.duration;
        }
    }
    
    for (int i = 1 ; i <= n; ++i) {
        if (last_j[i].min_t+last_j[i].duration > make_span) {
            make_span = last_j[i].min_t+last_j[i].duration;
        }
    }
    
    //输出
    clock_t t2 = clock();
    cout<<"time:"<<endl;
    float diff = ((float)(t2 - t1) / 1000000.0F );
    printf("%f\n",diff);
    cout<<make_span<<endl;
    
    //尾处理
    delete[] jobs;
    delete[] machines;
    delete[] front_steps;
    for (int i=1; i<=n; ++i) {
        delete [] buffer[i];
    }

    return 0;
}

int compare(const void* a,const void*b)
{
    step *s1=(step*)a,*s2=(step*)b;
    //使用如下三元组来确定优先级
    //<min_t,busy/duration,j_NO>
    if (s1->j_NO == 0) {//用于把空组排至最后
        return 1;
    }
    if (s2->j_NO == 0) {
        return -1;
    }
    if (s1->min_t < s2->min_t)
    {
        return -1;
    }
    else if (s1->min_t == s2->min_t)
    {
        if (s1->busy/s1->duration < s2->busy/s2->duration)
        {
            return -1;
        }
        else if (s1->busy/s1->duration == s2->busy/s2->duration)
        {
            if (s1->j_NO < s2->j_NO)
            {//编号不存在相等的可能性
                return -1;
            }
            else
            {
                return 1;
            }
        }
        else
        {
            return 1;
        }
    }
    else
    {
        return 1;
    }
}
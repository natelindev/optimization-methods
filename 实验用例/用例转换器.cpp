//
//  main.cpp
//  TransfromLA
//
//  Created by 林理露 on 16/6/4.
//  Copyright © 2016年 林理露. All rights reserved.
//

#include <iostream>
using namespace std;

int a[100][100],b[100][100];
int main(int argc, const char * argv[]) {
    freopen("/users/LLL/Desktop/la122.in", "r", stdin);
    freopen("/users/LLL/Desktop/la12.in", "w", stdout);
    int n=0,m=0;
    cin>>n>>m;
    cout<<n<<'p'<<m<<'m'<<m<<'j'<<endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <m; ++j) {
            cin>>a[i][j]>>b[i][j];
        }
    }
    cout<<endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <m; ++j) {
            cout<<b[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <m; ++j) {
            cout<<a[i][j]+1<<" ";
        }
        cout<<endl;
    }
    return 0;
}

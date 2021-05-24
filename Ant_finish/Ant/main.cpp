#include <iostream>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <cstdio>
#include <ctime>
#include <cmath>

using namespace std;

int just_check[101][101];//用于检测随机生成城市坐标时有无重复

//蚁群算法参数
int ant_num = 50;                        //蚂蚁数量
                                         //遵循蚁群算法的经验，种群数量取10~50
int alpha = 1;                           //信息素重要程度因子
int beta = 5;                            //启发函数重要程度因子
double rho = 0.1;                           //信息素挥发因子
int Q = 1;                               //常系数
int iter = 0;                            //迭代次数初值
int iter_max = 500;                      //最大迭代次
double** Eta;//启发矩阵，两个城市之间的启发函数值
double** Tau;//信息素矩阵，两个城市之间路径上的信息素的浓度
int** Table;//路径记录表,每一行表示第i只蚂蚁的路径

double** cdistance;
int** Route_best;//各代最佳路径
double* Length_best;//各代最佳路径的长度
double* Length_ave;//各代路径的平均长度
//蚁群算法参数

//临时参数
int* tabu;//对城市的访问情况,值为1表示已经访问过了,值为0表示没有访问
int* start;//每个蚂蚁的起点城市
//临时参数

class Cities
{
public:
    Cities(int n = 0);//构造函数设置城市数量，为坐标的变长数组分配空间
    ~Cities()
    {
        delete[] clocal_x;
        delete[] clocal_y;
    }
    void input_random();//录入城市的位置，坐标随机生成
    void input();//录入城市的位置，手动录入
    void output();

    int cnum;//城市数量
    int* clocal_x;//变长数组，用于存储城市的位置，相当于一张地图中的x,y坐标，x,y的取值为[0,100]
    int* clocal_y;

private:
};
Cities::Cities(int n)
{
    cnum = n;
    clocal_x = new int[n];
    clocal_y = new int[n];
    cdistance = new double*[n];
    for(int i = 0; i < cnum; i++)
    {
        cdistance[i] = new double[n];
    }
}
void Cities::input_random()
{
    if(cnum > 0)
    {
        srand((unsigned)time(0));
        for(int i = 0; i < 101; i++)
        {
            memset(just_check[i], 0, sizeof(int) * 101);
        }

        for(int i = 0; i < cnum; i++)
        {
            int tmp_x, tmp_y;
            tmp_x = rand() % 100;
            tmp_y = rand() % 100;
            if(just_check[tmp_x][tmp_y] == 1)
            {
                i--;
            }
            else
            {
                clocal_x[i] = tmp_x;
                clocal_y[i] = tmp_y;
                just_check[tmp_x][tmp_y] = 1;
            }
        }
    }

}
void Cities::input()
{
    if(cnum > 0)
    {
        for(int i = 0; i < 101; i++)
        {
            memset(just_check[i], 0, sizeof(int) * 101);
        }

        for(int i = 0; i < cnum; i++)
        {
            int tmp_x, tmp_y;
            scanf("%d%d", &tmp_x, &tmp_y);
            if(just_check[tmp_x][tmp_y] == 1)
            {
                i--;
            }
            else
            {
                clocal_x[i] = tmp_x;
                clocal_y[i] = tmp_y;
                just_check[tmp_x][tmp_y] = 1;
            }
        }
    }
}
void Cities::output()
{
    for(int i = 0; i < cnum; i++)
    {
        printf("city %d: (%d,%d)\n", i + 1, (clocal_x[i]), (clocal_y[i]));
    }
}

void init(Cities& c1);
void search_route(int n, int id);//对应每一只蚂蚁的寻路过程
int choose(int n, int current_city);//选择下一步要去的城市，并返回这个城市编号
void cal_results(int n, int iter_w);//整理每一次迭代中的寻路结果，整理出此次迭代的最优解信息并更新信息素浓度
void free_space(Cities& c1);

int main()
{
    //城市信息录入
    int city_num, i_mode;
    printf("请输入旅行商问题中城市的数量:\n");
    scanf("%d", &city_num);

    Cities c1(city_num);
    printf("请选择城市坐标信息的生成方式(城市的x,y坐标取值都在[0,100]):\n");
    printf("1.随机生成\n");
    printf("2.键盘录入(相同的(x,y)坐标值不计入)\n");
    scanf("%d", &i_mode);
    if(i_mode == 1)
    {
        c1.input_random();
    }
    else if(i_mode == 2)
    {
        c1.input();
    }

    printf("城市坐标如下:\n");
    c1.output();

    system("pause");
    //城市信息录入结束
    int n = c1.cnum;//城市数量
    init(c1);//初始化蚁群算法的参数

    ofstream out_local("City_local.txt");
    ofstream out_route("Route_best.txt");
    out_local << city_num << endl;
    for(int i = 0; i < n; i++)
    {
        out_local << c1.clocal_x[i] << " " << c1.clocal_y[i] << endl;
    }
    out_local.close();
    //输出蚁群算法找到的最短路径
    printf("开始执行蚁群算法\n");
    clock_t t1 = clock();
    while(iter < iter_max)
    {
        for(int i = 0; i < ant_num; i++)//逐个蚂蚁
        {
            start[i] = rand() % n;//随机产生各个蚂蚁的起点城市
            Table[i][0] = start[i];//构建解空间
        }

        for(int i = 0; i < ant_num; i++)//逐个蚂蚁路径选择
        {
            memset(tabu, 0, sizeof(int) * n);//初始化禁忌表
            search_route(n, i);
        }
        cal_results(n, iter);

        printf("第%d次迭代:\n", iter + 1);
        printf("此次迭代的路径平均距离为: %f\n此次迭代的最优路线如下\n", Length_ave[iter]);
        for(int i = 0; i < n - 1; i++)
        {
            printf("city %d -> city %d\n", Route_best[iter][i] + 1, Route_best[iter][i + 1] + 1);
        }
        printf("city %d -> city %d\n", Route_best[iter][n - 1] + 1, Route_best[iter][0] + 1);
        for(int i = 0; i < n; i++)
        {
            out_route << Route_best[iter][i] + 1 << endl;
        }
        iter++;
        for(int i = 0; i < ant_num; i++)//更新路径记录表
            memset(Table[i], 0, sizeof(int) * n);
    }

    clock_t t2 = clock();
    printf("蚁群算法执行结束，耗时%d ms\n", (int)(t2 - t1));

    int min_index = 0;
    for(int i = 0; i < iter_max; i++)
    {
        if(Length_best[i] < Length_best[min_index])
            min_index = i;
    }
    double min_len =Length_best[min_index];
    printf("最短路径的长度为:%f\n", min_len);
    printf("找到的最短路径如下:\n");
    printf("Starting city: %d, Ending city: %d\n", Route_best[min_index][0] + 1,
        Route_best[min_index][n - 1] + 1);
    for(int i = 0; i < n - 1; i++)
    {
        printf("city %d -> city %d\n", (Route_best[min_index][i]) + 1, (Route_best[min_index][i + 1]) + 1);
    }
    printf("city %d -> city %d\n", (int)((Route_best[min_index][n - 1]) + 1), (int)((Route_best[min_index][0]) + 1));
    //输出蚁群算法找到的最短路径
    free_space(c1);
    printf("算法到此结束\n");
    return 0;
}
void init(Cities& c1)
{
    int n = c1.cnum;                         //城市数量
    //初始化蚁群算法的参数
    //城市距离
    cdistance = new double*[n];
    Eta = new double*[n];
    Tau = new double*[n];
    for(int i = 0; i < n; i++)
    {
        cdistance[i] = new double[n];
        Eta[i] = new double[n];
        Tau[i] = new double[n];
    }
    //城市距离
    //启发函数(启发矩阵)
    //信息素矩阵
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            cdistance[i][j] = sqrt(pow((c1.clocal_x[i] - c1.clocal_x[j]), 2) +
                             pow((c1.clocal_y[i] - c1.clocal_y[j]), 2));
            if(i != j)
            {
                Eta[i][j] = 1 / cdistance[i][j];//城市点不重复的情况下距离不可能为0
            }
            else
            {
                Eta[i][j] = 1;
            }
            Tau[i][j] = 1;
        }
    }
    //城市距离
    //启发函数(启发矩阵)
    //信息素矩阵

    //路径记录表
    Table = new int*[ant_num];
    for(int i = 0; i < ant_num; i++)
    {
        Table[i] = new int[n];
    }
    for(int i = 0; i < ant_num; i++)
    {
        memset(Table[i], 0, sizeof(int) * n);
    }
    //路径记录表

    //各代最佳路径
    Route_best = new int*[iter_max];
    for(int i = 0; i < iter_max; i++)
    {
        Route_best[i] = new int[n];
    }
    for(int i = 0; i < iter_max; i++)
    {
        memset(Route_best[i], 0, sizeof(int) * n);
    }
    //各代最佳路径

    //各代最佳路径的长度
    Length_best = new double[iter_max];
    memset(Length_best, 0, sizeof(double) * iter_max);
    //各代最佳路径的长度

    //各代路径的平均长度
    Length_ave = new double[iter_max];
    memset(Length_ave, 0, sizeof(double) * iter_max);
    //各代路径的平均长度

    //初始化蚁群算法的参数
    tabu = new int[n];
    start = new int[ant_num];
}
void search_route(int n, int id)
{
    int i;
    int current_city;
    for(i = 0; i < n - 1; i++)//第id只蚂蚁逐个城市路径选择，起点城市已经经过一次了
    {
        tabu[Table[id][i]] = 1;//更新禁忌表
        current_city = Table[id][i];//获取当前所在城市
        //选择下一个访问城市
        int result = choose(n, current_city);
        Table[id][i + 1] = result;
        //选择下一个访问城市
    }
}
int choose(int n, int current_city)
{
    //计算城市间转移概率
    double p_transfer[n];
    memset(p_transfer, 0, sizeof(double) * n);
    double sum = 0;

    for(int i = 0; i < n; i++)
    {
        if(tabu[i] == 0)//当前城市和禁忌表中城市的转移概率是0
        {
            p_transfer[i] = pow(Tau[current_city][i], alpha) * pow(Eta[current_city][i], beta);
            sum += p_transfer[i];
        }
    }
    for(int i = 0; i < n; i++)
    {
        p_transfer[i] = p_transfer[i] / sum;
    }
    //计算城市间转移概率

    ////轮盘赌法选择下一个访问城市
    srand((unsigned)time(0));

    double pc[n];
    int target = 0;
    for(int i = 0; i < n; i++)
    {
        if(tabu[i] == 0)
        {
            target = i;//target的初值是1个没有去过的城市，为了避免轮盘赌找不到
        }
    }
    memset(pc, 0, sizeof(pc));
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < i + 1; j++)
            pc[i] += p_transfer[j];
    }
    double number = rand() % (999 + 1) / (double)(999 + 1);//生成一个0~1之间的小数，
    //并且精确到小数点后三位

    for(int i = 0; i < n; i++)
    {
//        int pc_tmp = pc[i] * 1000;
//        int n_tmp = number * 1000;
//        if((fabs(pc[i]) >= fabs(number)) && tabu[i] == 0);
        if((pc[i] >= number) && tabu[i] == 0)
        {
            target = i;
            break;
        }

    }
    return target;
    //轮盘赌法选择下一个访问城市
}
void cal_results(int n, int iter_w)//整理每一次迭代中的寻路结果，整理出此次迭代的最优解信息
{
    //计算每只蚂蚁的路径距离
    int i, j;
    double ant_len[ant_num];
    memset(ant_len, 0, sizeof(double)*ant_num);
    for(i = 0; i < ant_num; i++)
    {
        for(j = 0; j < n - 1; j++)
        {
            ant_len[i] += cdistance[Table[i][j]][Table[i][j + 1]];
        }
        ant_len[i] += cdistance[Table[i][j]][Table[i][0]];
    }
    //计算每只蚂蚁的路径距离

    //计算各自最短路径距离及平均距离
    int min_index = 0;
    for(i = 1; i < ant_num; i++)
    {
        if(ant_len[i] < ant_len[min_index])
            min_index = i;
    }
    Length_best[iter] = ant_len[min_index];

    double asum = 0;
    for(i = 0; i < ant_num; i++)
    {
        asum += ant_len[i];
    }
    Length_ave[iter_w] = asum / ant_num;
    for(i = 0; i < n; i++)
    {
        Route_best[iter_w][i] = Table[min_index][i];
    }
    //计算各自最短路径距离及平均距离

    //信息素浓度更新
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            Tau[i][j] = (1 - rho) * Tau[i][j];
        }
    }
    for(i = 0; i < ant_num; i++)//逐个蚂蚁计算
    {
        for(j = 0; j < n - 1; j++)//逐个城市计算
        {
            Tau[Table[i][j]][Table[i][j + 1]] += Q / ant_len[i];
        }
        Tau[Table[i][j]][Table[i][0]] += Q / ant_len[i];
    }
    //信息素浓度更新
}

void free_space(Cities& c1)
{
    int n = c1.cnum;
    for(int i = 0; i < n; i++)
    {
        delete[] Eta[i];
        delete[] Tau[i];
        delete[] cdistance[i];
    }
    delete[] Eta;
    delete[] Tau;
    delete[] cdistance;

    for(int i = 0; i < ant_num; i++)
    {
        delete[] Table[i];
    }
    delete[] Table;

    for(int i = 0; i < iter_max; i++)
    {
        delete[] Route_best[i];
    }
    delete[] Route_best;

    delete[] Length_best;
    delete[] Length_ave;

    delete[] tabu;
    delete[] start;
}

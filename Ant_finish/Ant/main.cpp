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

int just_check[101][101];//���ڼ��������ɳ�������ʱ�����ظ�

//��Ⱥ�㷨����
int ant_num = 50;                        //��������
                                         //��ѭ��Ⱥ�㷨�ľ��飬��Ⱥ����ȡ10~50
int alpha = 1;                           //��Ϣ����Ҫ�̶�����
int beta = 5;                            //����������Ҫ�̶�����
double rho = 0.1;                           //��Ϣ�ػӷ�����
int Q = 1;                               //��ϵ��
int iter = 0;                            //����������ֵ
int iter_max = 500;                      //��������
double** Eta;//����������������֮�����������ֵ
double** Tau;//��Ϣ�ؾ�����������֮��·���ϵ���Ϣ�ص�Ũ��
int** Table;//·����¼��,ÿһ�б�ʾ��iֻ���ϵ�·��

double** cdistance;
int** Route_best;//�������·��
double* Length_best;//�������·���ĳ���
double* Length_ave;//����·����ƽ������
//��Ⱥ�㷨����

//��ʱ����
int* tabu;//�Գ��еķ������,ֵΪ1��ʾ�Ѿ����ʹ���,ֵΪ0��ʾû�з���
int* start;//ÿ�����ϵ�������
//��ʱ����

class Cities
{
public:
    Cities(int n = 0);//���캯�����ó���������Ϊ����ı䳤�������ռ�
    ~Cities()
    {
        delete[] clocal_x;
        delete[] clocal_y;
    }
    void input_random();//¼����е�λ�ã������������
    void input();//¼����е�λ�ã��ֶ�¼��
    void output();

    int cnum;//��������
    int* clocal_x;//�䳤���飬���ڴ洢���е�λ�ã��൱��һ�ŵ�ͼ�е�x,y���꣬x,y��ȡֵΪ[0,100]
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
void search_route(int n, int id);//��Ӧÿһֻ���ϵ�Ѱ·����
int choose(int n, int current_city);//ѡ����һ��Ҫȥ�ĳ��У�������������б��
void cal_results(int n, int iter_w);//����ÿһ�ε����е�Ѱ·�����������˴ε��������Ž���Ϣ��������Ϣ��Ũ��
void free_space(Cities& c1);

int main()
{
    //������Ϣ¼��
    int city_num, i_mode;
    printf("�����������������г��е�����:\n");
    scanf("%d", &city_num);

    Cities c1(city_num);
    printf("��ѡ�����������Ϣ�����ɷ�ʽ(���е�x,y����ȡֵ����[0,100]):\n");
    printf("1.�������\n");
    printf("2.����¼��(��ͬ��(x,y)����ֵ������)\n");
    scanf("%d", &i_mode);
    if(i_mode == 1)
    {
        c1.input_random();
    }
    else if(i_mode == 2)
    {
        c1.input();
    }

    printf("������������:\n");
    c1.output();

    system("pause");
    //������Ϣ¼�����
    int n = c1.cnum;//��������
    init(c1);//��ʼ����Ⱥ�㷨�Ĳ���

    ofstream out_local("City_local.txt");
    ofstream out_route("Route_best.txt");
    out_local << city_num << endl;
    for(int i = 0; i < n; i++)
    {
        out_local << c1.clocal_x[i] << " " << c1.clocal_y[i] << endl;
    }
    out_local.close();
    //�����Ⱥ�㷨�ҵ������·��
    printf("��ʼִ����Ⱥ�㷨\n");
    clock_t t1 = clock();
    while(iter < iter_max)
    {
        for(int i = 0; i < ant_num; i++)//�������
        {
            start[i] = rand() % n;//��������������ϵ�������
            Table[i][0] = start[i];//������ռ�
        }

        for(int i = 0; i < ant_num; i++)//�������·��ѡ��
        {
            memset(tabu, 0, sizeof(int) * n);//��ʼ�����ɱ�
            search_route(n, i);
        }
        cal_results(n, iter);

        printf("��%d�ε���:\n", iter + 1);
        printf("�˴ε�����·��ƽ������Ϊ: %f\n�˴ε���������·������\n", Length_ave[iter]);
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
        for(int i = 0; i < ant_num; i++)//����·����¼��
            memset(Table[i], 0, sizeof(int) * n);
    }

    clock_t t2 = clock();
    printf("��Ⱥ�㷨ִ�н�������ʱ%d ms\n", (int)(t2 - t1));

    int min_index = 0;
    for(int i = 0; i < iter_max; i++)
    {
        if(Length_best[i] < Length_best[min_index])
            min_index = i;
    }
    double min_len =Length_best[min_index];
    printf("���·���ĳ���Ϊ:%f\n", min_len);
    printf("�ҵ������·������:\n");
    printf("Starting city: %d, Ending city: %d\n", Route_best[min_index][0] + 1,
        Route_best[min_index][n - 1] + 1);
    for(int i = 0; i < n - 1; i++)
    {
        printf("city %d -> city %d\n", (Route_best[min_index][i]) + 1, (Route_best[min_index][i + 1]) + 1);
    }
    printf("city %d -> city %d\n", (int)((Route_best[min_index][n - 1]) + 1), (int)((Route_best[min_index][0]) + 1));
    //�����Ⱥ�㷨�ҵ������·��
    free_space(c1);
    printf("�㷨���˽���\n");
    return 0;
}
void init(Cities& c1)
{
    int n = c1.cnum;                         //��������
    //��ʼ����Ⱥ�㷨�Ĳ���
    //���о���
    cdistance = new double*[n];
    Eta = new double*[n];
    Tau = new double*[n];
    for(int i = 0; i < n; i++)
    {
        cdistance[i] = new double[n];
        Eta[i] = new double[n];
        Tau[i] = new double[n];
    }
    //���о���
    //��������(��������)
    //��Ϣ�ؾ���
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            cdistance[i][j] = sqrt(pow((c1.clocal_x[i] - c1.clocal_x[j]), 2) +
                             pow((c1.clocal_y[i] - c1.clocal_y[j]), 2));
            if(i != j)
            {
                Eta[i][j] = 1 / cdistance[i][j];//���е㲻�ظ�������¾��벻����Ϊ0
            }
            else
            {
                Eta[i][j] = 1;
            }
            Tau[i][j] = 1;
        }
    }
    //���о���
    //��������(��������)
    //��Ϣ�ؾ���

    //·����¼��
    Table = new int*[ant_num];
    for(int i = 0; i < ant_num; i++)
    {
        Table[i] = new int[n];
    }
    for(int i = 0; i < ant_num; i++)
    {
        memset(Table[i], 0, sizeof(int) * n);
    }
    //·����¼��

    //�������·��
    Route_best = new int*[iter_max];
    for(int i = 0; i < iter_max; i++)
    {
        Route_best[i] = new int[n];
    }
    for(int i = 0; i < iter_max; i++)
    {
        memset(Route_best[i], 0, sizeof(int) * n);
    }
    //�������·��

    //�������·���ĳ���
    Length_best = new double[iter_max];
    memset(Length_best, 0, sizeof(double) * iter_max);
    //�������·���ĳ���

    //����·����ƽ������
    Length_ave = new double[iter_max];
    memset(Length_ave, 0, sizeof(double) * iter_max);
    //����·����ƽ������

    //��ʼ����Ⱥ�㷨�Ĳ���
    tabu = new int[n];
    start = new int[ant_num];
}
void search_route(int n, int id)
{
    int i;
    int current_city;
    for(i = 0; i < n - 1; i++)//��idֻ�����������·��ѡ���������Ѿ�����һ����
    {
        tabu[Table[id][i]] = 1;//���½��ɱ�
        current_city = Table[id][i];//��ȡ��ǰ���ڳ���
        //ѡ����һ�����ʳ���
        int result = choose(n, current_city);
        Table[id][i + 1] = result;
        //ѡ����һ�����ʳ���
    }
}
int choose(int n, int current_city)
{
    //������м�ת�Ƹ���
    double p_transfer[n];
    memset(p_transfer, 0, sizeof(double) * n);
    double sum = 0;

    for(int i = 0; i < n; i++)
    {
        if(tabu[i] == 0)//��ǰ���кͽ��ɱ��г��е�ת�Ƹ�����0
        {
            p_transfer[i] = pow(Tau[current_city][i], alpha) * pow(Eta[current_city][i], beta);
            sum += p_transfer[i];
        }
    }
    for(int i = 0; i < n; i++)
    {
        p_transfer[i] = p_transfer[i] / sum;
    }
    //������м�ת�Ƹ���

    ////���̶ķ�ѡ����һ�����ʳ���
    srand((unsigned)time(0));

    double pc[n];
    int target = 0;
    for(int i = 0; i < n; i++)
    {
        if(tabu[i] == 0)
        {
            target = i;//target�ĳ�ֵ��1��û��ȥ���ĳ��У�Ϊ�˱������̶��Ҳ���
        }
    }
    memset(pc, 0, sizeof(pc));
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < i + 1; j++)
            pc[i] += p_transfer[j];
    }
    double number = rand() % (999 + 1) / (double)(999 + 1);//����һ��0~1֮���С����
    //���Ҿ�ȷ��С�������λ

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
    //���̶ķ�ѡ����һ�����ʳ���
}
void cal_results(int n, int iter_w)//����ÿһ�ε����е�Ѱ·�����������˴ε��������Ž���Ϣ
{
    //����ÿֻ���ϵ�·������
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
    //����ÿֻ���ϵ�·������

    //����������·�����뼰ƽ������
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
    //����������·�����뼰ƽ������

    //��Ϣ��Ũ�ȸ���
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            Tau[i][j] = (1 - rho) * Tau[i][j];
        }
    }
    for(i = 0; i < ant_num; i++)//������ϼ���
    {
        for(j = 0; j < n - 1; j++)//������м���
        {
            Tau[Table[i][j]][Table[i][j + 1]] += Q / ant_len[i];
        }
        Tau[Table[i][j]][Table[i][0]] += Q / ant_len[i];
    }
    //��Ϣ��Ũ�ȸ���
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

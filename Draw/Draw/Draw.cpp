#include <iostream>
#include <cstring>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <fstream>
#include "engine.h"
using namespace std;

#define ITER_MAX 500
//迭代次数固定
int main()
{
	Engine* ep;
	mxArray* city_x = NULL;
	mxArray* city_y = NULL;//用来存储城市坐标
	if ((ep = engOpen("")) == NULL)
	{
		printf("Engine Fail");
	}

	ifstream in_route, in_local;
	in_route.open("Route_best.txt", ios::in);
	in_local.open("City_local.txt", ios::in);

	if (!(in_route.is_open())) {
		cout << "route" << endl;
	}
	if (!(in_local.is_open())) {
		cout << "local" << endl;
	}

	int city_num, i = 0, j;
	char alpha;
	char c_num[99];
	//读入City_local.txt文件中第一行的数字
	while (in_local.get(alpha))
	{
		if (alpha == '\n')
			break;
		else
		{
			c_num[i] = alpha;
			i++;
		}
	}
	c_num[i] = '\0';
	//读入City_local.txt文件中第一行的数字

	//得到城市数目
	city_num = atoi(c_num);
	//得到城市数目

	double* local_x = NULL;//存放所有城市的x坐标
	double* local_y = NULL;//存放所有城市的x坐标
	local_x = new double[city_num];
	local_y = new double[city_num];//为了配合matlab，此处用来double型而不是int
	//读取城市坐标
	j = 0;
	for (i = 0; i < city_num; i++)
	{
		int tmp;
		while (in_local.get(alpha))//循环的一次只读取一行
		{
			if (alpha == '\n')
			{
				c_num[j] = '\0';
				tmp = atoi(c_num);
				local_y[i] = (double)tmp;
				j = 0;
				memset(c_num, 0, sizeof(c_num));
				break;
			}
			else if (alpha == ' ')
			{
				c_num[j] = '\0';
				tmp = atoi(c_num);
				local_x[i] = (double)tmp;
				j = 0;
				memset(c_num, 0, sizeof(c_num));
			}
			else
			{
				c_num[j] = alpha;
				j++;
			}
		}
	}
	in_local.close();
	//读取城市坐标

	//城市坐标的数组->mxArray转储
	city_x = mxCreateDoubleMatrix(1, city_num, mxREAL);
	city_y = mxCreateDoubleMatrix(1, city_num, mxREAL);
	//城市坐标的数组->mxArray转储

	memcpy((void*)mxGetPr(city_x), (void*)local_x, sizeof(double) * city_num);
	memcpy((void*)mxGetPr(city_y), (void*)local_y, sizeof(double) * city_num);
	/*for (int i = 0; i < city_num; i++)
	{
		cout << *(mxGetPr(city_x) + i) << " ";
		cout << *(mxGetPr(city_y) + i) << endl;
	}*/
	engPutVariable(ep, "local_x", city_x);
	engPutVariable(ep, "local_y", city_y);

	//设置坐标轴信息
	engEvalString(ep, "hold on");
	engEvalString(ep, "set(gca, 'Xlim', [0,100])");
	engEvalString(ep, "set(gca, 'XTick', [0:10:100])");
	engEvalString(ep, "set(gca, 'XTickLabel', [0:10:100])");
	engEvalString(ep, "set(gca, 'Ylim', [0,100])");
	engEvalString(ep, "set(gca, 'YTick', [0:10:100])");
	engEvalString(ep, "set(gca, 'YTickLabel', [0:10:100])");
	engEvalString(ep, "xlabel('x')");
	engEvalString(ep, "ylabel('y')");
	engEvalString(ep, "title('城市坐标图')");

	//设置坐标轴信息

	//初步绘制散点图(以城市坐标绘制的散点图)
	engEvalString(ep, "scatter(local_x, local_y, 'markerfacecolor', [ 1, 0, 0 ] )");
	engEvalString(ep, "hold on");
	engEvalString(ep, "pause(8)");
	//初步绘制散点图(以城市坐标绘制的散点图)



	//读取每一次迭代的最优路线,并调用MATLAB引擎绘制函数
	double* route = NULL;//存放路径
	route = new double[city_num];


	mxArray* best_route = NULL;
	for (i = 1; i <= ITER_MAX; i++)
	{
		int my_count = 0, tmp;//为了配合matlab下标从1开始的规则
		j = 0;
		while (in_route.get(alpha))
		{
			if (alpha == '\n')
			{
				c_num[j] = '\0';
				tmp = atoi(c_num);
				route[my_count] = (double)tmp;

				j = 0;
				memset(c_num, 0, sizeof(c_num));
				my_count++;
				if (my_count == city_num)
				{
					my_count = 0;

					//设置坐标轴信息
					engEvalString(ep, "hold on");
					engEvalString(ep, "set(gca, 'Xlim', [0,100])");
					engEvalString(ep, "set(gca, 'XTick', [0:10:100])");
					engEvalString(ep, "set(gca, 'XTickLabel', [0:10:100])");
					engEvalString(ep, "xlabel('x')");
					engEvalString(ep, "ylabel('y')");
					engEvalString(ep, "title('TSP旅行商问题迭代最优解')");
					engEvalString(ep, "set(gca, 'Ylim', [0,100])");
					engEvalString(ep, "set(gca, 'YTick', [0:10:100])");
					engEvalString(ep, "set(gca, 'YTickLabel', [0:10:100])");
					//设置坐标轴信息

					//初步绘制散点图(以城市坐标绘制的散点图)
					engEvalString(ep, "scatter(local_x, local_y, 'markerfacecolor', [ 1, 0, 0 ] )");
					engEvalString(ep, "hold on");
					//初步绘制散点图(以城市坐标绘制的散点图)

					best_route = mxCreateDoubleMatrix(1, city_num, mxREAL);
					memcpy((void*)mxGetPr(best_route), (void*)route, sizeof(double) * city_num);
					engPutVariable(ep, "route", best_route);

					engEvalString(ep, "for k = 1:(length(local_x));plot([local_x(1, route(k)), local_x(1, route(k+1))], [local_y(1, route(k)), local_y(1,route(k + 1))], 'color', 'black', 'Linewidth', 1, 'LineStyle', '-');k = k + 1;hold on;end;");
					engEvalString(ep, "plot([local_x(1, route(length(local_x))), local_x(1, route(1))], [local_y(1, route(length(local_x))), local_y(1, route(1))], 'color', 'black', 'Linewidth', 1, 'LineStyle', '-');");
					engEvalString(ep, "pause(2);");
					if (i != ITER_MAX)
					{
						engEvalString(ep, "clf");
					}
					else
					{
						engEvalString(ep, "pause(10);");
					}

					memset(route, 0, sizeof(double) * city_num);//一次迭代结束，清空路线表
					break;
				}

			}
			else
			{
				c_num[j] = alpha;
				j++;
			}
		}
		//本次迭代结束
	}
	in_route.close();
	//读取每一次迭代的最优路线,并调用MATLAB引擎绘制函数

	delete[] local_x;
	delete[] local_y;
	delete[] route;
	engClose(ep);
	return 0;
}

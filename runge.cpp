//Крамер Герман Владиславович
//численное решение дифференциального уравнения
//с использовнием классического метода Рунге-Кутты 4 порядка

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
using std::vector;

const int x0 = 0; //начало отрезка
const int x_end = 5; //конец отрезка
const double EPS = 0.001; //точность

const vector<double> Yi = { 0.0, 3.0, -9.0, -8.0, 0.0 };    //массив начальный условий y0, y1, y2, y3, y4 соответственно

double fx5(const vector<double>& Yi)    //исходное уравнение, приведенное к виду y5 = f(x,y)
{
    return -15 * Yi[4] - 90 * Yi[3] - 270 * Yi[2] - 405 * Yi[1] - 243 * Yi[0];
} 

vector<double> find_k(double step, double yn); //прототипы функций
vector<double> runge_kutta(double step);


int main()
{
    double step = 0.1; //шаг разбиения отрезка

    vector<double> answerY1;
    vector<double> answerY2; // массивый для хранения двух последовательных вычислений
                             // используются для проверки достижения нужной точности
    double n1;
    double n2;
    do
    {
        step /= 2;
        answerY1 = runge_kutta(step);
        step /= 2;
        answerY2 = runge_kutta(step); //вычисление с уменьшенным вдвое шагом
        n1 = answerY1.back();
        n2 = answerY2.back();

    } while (abs(n2 - n1) > EPS); //n2 и n1 значения вычисленной функции на конце отрезка

    double x = 0;
    
    std::ofstream fout;
    fout.open("./answer.txt");
    std::cout << "Answer is:\n";
    fout << "Answer is:\n";    
    fout.width(10);
    for (auto i : answerY2)
    {   
        std::cout << std::setw(10) << x << "\t" << std::setw(10) << i << '\n';
        fout << std::setw(10) << x << "\t" << std::setw(10) << i << '\n';
        x += step;

    }
    fout.close();
}


vector<double> find_k(double step, double yn) //функция нахождения коэффициентов к для Рунге-Кутты
{
    vector<double> k(4);   
    k[0] = yn;
    k[1] = yn + step * k[0] / 2;
    k[2] = yn + step * k[1] / 2;
    k[3] = yn + step * k[2];
    return k;
}

vector<double> runge_kutta(double step) //функция реализации метода Рунге-Кутты
{
    double y5, y4, y3, y2, y1; //переменные для хранения численного значения конкретной производной
    vector<double> y0; // массив для хранения вычисляемых значений искомой функции
        
    y4 = Yi[4];
    y3 = Yi[3];
    y2 = Yi[2];
    y1 = Yi[1];  
    y0.push_back(Yi[0]); // динамической добавление вычисляемых значений искомой функции

    vector<double> for_fx5 = { y0[0], y1, y2, y3, y4 }; //массив для передачи в исходную функцию y5 = f(x,y)
    y5 = fx5(for_fx5); //вычисление y5 в 0

    vector<double> k = find_k(step, y5);
    int index = 1;

    while (index <= int(abs(x_end - x0) / step))
    {           
        y4 = y4 + step / 6 * (k[0] + k[1] + k[2] + k[3]);   //последовательное вычисление соответственных производных для дальнейшего вычисления самой функции
        k = find_k(step, y4);
        y3 = y3 + step / 6 * (k[0] + k[1] + k[2] + k[3]);
        k = find_k(step, y3);
        y2 = y2 + step / 6 * (k[0] + k[1] + k[2] + k[3]);
        k = find_k(step, y2);
        y1 = y1 + step / 6 * (k[0] + k[1] + k[2] + k[3]);
        k = find_k(step, y1);

        y0.push_back(y0[index - 1] + step / 6 * (k[0] + k[1] + k[2] + k[3])); // динамической добавление вычисляемых значений искомой функции

        for_fx5 = { y0[index], y1, y2, y3, y4 };
        y5 = fx5(for_fx5);  //вычисление y5 в следующей точке
        k = find_k(step, y5);
        index++;
    }

    return y0;
}

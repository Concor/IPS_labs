#include <iostream>
#include <math.h>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <thread>
#include <mutex>

# define M_PI 3.14159265358979323846 /* pi, #include <math.h> */

using namespace std;
using namespace std::chrono;

double a = -1, b = 1, computedInt = 2*M_PI;

// функци€ интеграла
double f(double x) {
	return (8 / (2+2*x*x));
}
// вычисление интеграла без дополнительных инструкций
void computeIntegral(double N) {
	const double width = (b - a) / N;
	double trapezoidal_integral = 0;
	high_resolution_clock::time_point start = high_resolution_clock::now();
	for (int step = 0; step < N; step++) {
		const double x1 = a + step * width;
		const double x2 = a + (step + 1) * width;
		trapezoidal_integral += 0.5 * (x2 - x1) * (f(x1) + f(x2));
	}
	high_resolution_clock::time_point finish = high_resolution_clock::now();
	duration<double> duration = (finish - start);
	cout << "Duration is: " << duration.count() << " seconds" << endl;
	cout << "EPS is: " << setprecision(10) << abs(computedInt - trapezoidal_integral) << endl;
}
// вычисление интеграла с распараллеливанием на 8 потоков с помощью pramga
void computeIntegralParallel(double N) {
	const double width = (b - a) / N;
	double trapezoidal_integral = 0;
	high_resolution_clock::time_point start = high_resolution_clock::now();
	#pragma loop(hint_parallel(4))
	for (int step = 0; step < N; step++) {
		const double x1 = a + step * width;
		const double x2 = a + (step + 1) * width;
		trapezoidal_integral += 0.5 * (x2 - x1) * (f(x1) + f(x2));
	}
	high_resolution_clock::time_point finish = high_resolution_clock::now();
	duration<double> duration = (finish - start);
	cout << "Duration is: " << duration.count() << " seconds" << endl;
	cout << "EPS is: " << setprecision(10) << abs(computedInt - trapezoidal_integral) << endl;
}

// глобальные вспомогательные переменные
mutex MUTEX_VAR;
double sum_all = 0.0;
//вспомогательна€ функци€ дл€ ручного создани€ потоков
void INT(double start, double end, const double width) {
	double trapezoidal_integral = 0;
	for (int step = start; step < end; step++) {
		const double x1 = a + step * width;
		const double x2 = a + (step + 1) * width;
		trapezoidal_integral += 0.5 * (x2 - x1) * (f(x1) + f(x2));
	}
	MUTEX_VAR.lock();
	sum_all += trapezoidal_integral;
	MUTEX_VAR.unlock();
}
// создание тредов вручную
void computeIntegralThreads(double N) {
	const double width = (b - a) / N;

	high_resolution_clock::time_point start = high_resolution_clock::now();
	int d = N / 4;
	thread t1(INT, 0*d, 1*d, width);
	thread t2(INT, 1*d, 2*d, width);
	thread t3(INT, 2*d, 3*d, width);
	thread t4(INT, 3*d, 4*d, width);
	t1.join();
	t2.join();
	t3.join();
	t4.join();

	high_resolution_clock::time_point finish = high_resolution_clock::now();
	duration<double> duration = (finish - start);
	cout << "Duration is: " << duration.count() << " seconds" << endl;
	cout << "EPS is: " << setprecision(10) << abs(computedInt - sum_all) << endl;
	sum_all = 0;
}



int main() {
	double N0 = 100;
	double N1 = 1000;
	double N2 = 10000;
	double N3 = 100000;
	double N4 = 1000000;
	double N5 = 10000000;
	cout << computedInt << endl;
	cout << "-----------\nNo Parallel\n-----------" << endl;
	computeIntegral(N0);
	computeIntegral(N1);
	computeIntegral(N2);
	computeIntegral(N3);
	computeIntegral(N4);
	computeIntegral(N5);
	cout << "-----------\nWith Parallel\n-----------" << endl;
	computeIntegralParallel(N0);
	computeIntegralParallel(N1);
	computeIntegralParallel(N2);
	computeIntegralParallel(N3);
	computeIntegralParallel(N4);
	computeIntegralParallel(N5);
	cout << "-----------\nWith threads\n-----------" << endl;
	computeIntegralThreads(N0);
	computeIntegralThreads(N1);
	computeIntegralThreads(N2);
	computeIntegralThreads(N3);
	computeIntegralThreads(N4);
	computeIntegralThreads(N5);
}
#include <stdio.h>
#include <ctime>
//#include <cilk/cilk.h>
//#include <cilk/reducer_opadd.h>
#include <chrono>

using namespace std;

// количество строк в исходной квадратной матрице
const int MATRIX_SIZE = 1500;

/// Функция InitMatrix() заполняет переданную в качестве 
/// параметра квадратную матрицу случайными значениями
/// matrix - исходная матрица СЛАУ
void InitMatrix(double** matrix)
{
	for (int i = 0; i < MATRIX_SIZE; ++i) {
		matrix[i] = new double[MATRIX_SIZE + 1];
	}

	for (int i = 0; i < MATRIX_SIZE; ++i) {
		for (int j = 0; j <= MATRIX_SIZE; ++j) {
			matrix[i][j] = rand() % 2500 + 1;
		}
	}
}

void DeInitMatrix(double** matrix) {
	for (int i = 0; i < MATRIX_SIZE; ++i) {
		delete[]matrix[i];
	}
}

/// Функция SerialGaussMethod() решает СЛАУ методом Гаусса 
/// matrix - исходная матрица коэффиициентов уравнений, входящих в СЛАУ,
/// последний столбей матрицы - значения правых частей уравнений
/// rows - количество строк в исходной матрице
/// result - массив ответов СЛАУ
chrono::duration<double> SerialGaussMethod(double** matrix, const int rows, double* result)
{
	int k;
	double koef;
	
	// прямой ход метода Гаусса
	chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
	for (k = 0; k < rows; ++k) {
		for (int i = k + 1; i < rows; ++i) {
			koef = -matrix[i][k] / matrix[k][k];
			for (int j = k; j <= rows; ++j) {
				matrix[i][j] += koef * matrix[k][j];
			}
		}
	}
	chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();

	// обратный ход метода Гаусса
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];
	for (k = rows - 2; k >= 0; --k) {
		result[k] = matrix[k][rows];
		for (int j = k + 1; j < rows; ++j) {
			result[k] -= matrix[k][j] * result[j];
		}
		result[k] /= matrix[k][k];
	}

	return chrono::duration<double>(t2 - t1);
}

chrono::duration<double> ParallelGaussMethod(double** matrix, const int rows, double* result)
{
	int k;
	double koef;

	// прямой ход метода Гаусса
	chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
	for (k = 0; k < rows; ++k) {
		for (int i = k + 1; i < rows; ++i) {
			koef = -matrix[i][k] / matrix[k][k];
#pragma omp parallel for
			for (int j = k; j <= rows; ++j) {
				matrix[i][j] += koef * matrix[k][j];
			}
		}
	}
	chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();

	// обратный ход метода Гаусса
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];
	for (k = rows - 2; k >= 0; --k) {
		result[k] = matrix[k][rows];
#pragma omp parallel for
		for (int j = k + 1; j < rows; ++j) {
			result[k] -= matrix[k][j] * result[j];
		}
		result[k] /= matrix[k][k];
	}

	return chrono::duration<double>(t2 - t1);
}

#define MATRIX_IN_USE 2

int main()
{
	srand((unsigned)time(0));

	int i;
#if MATRIX_IN_USE == 1
	const int test_matrix_lines = 4;
	double** test_matrix = new double* [test_matrix_lines];

	for (i = 0; i < test_matrix_lines; ++i) {
		test_matrix[i] = new double[test_matrix_lines + 1];
	}
	
	double* result = new double[test_matrix_lines];
	// инициализация тестовой матрицы
	test_matrix[0][0] = 2; test_matrix[0][1] = 5;  test_matrix[0][2] = 4;  test_matrix[0][3] = 1;  test_matrix[0][4] = 20;
	test_matrix[1][0] = 1; test_matrix[1][1] = 3;  test_matrix[1][2] = 2;  test_matrix[1][3] = 1;  test_matrix[1][4] = 11;
	test_matrix[2][0] = 2; test_matrix[2][1] = 10; test_matrix[2][2] = 9;  test_matrix[2][3] = 7;  test_matrix[2][4] = 40;
	test_matrix[3][0] = 3; test_matrix[3][1] = 8;  test_matrix[3][2] = 9;  test_matrix[3][3] = 2;  test_matrix[3][4] = 37;

	chrono::duration<double> elapsed_time = SerialGaussMethod(test_matrix, test_matrix_lines, result);

	for (i = 0; i < test_matrix_lines; ++i) {
		delete[]test_matrix[i];
	}

	printf("Solution (elapsed time = %.8f seconds):\n", elapsed_time.count());
	for (i = 0; i < test_matrix_lines; ++i) {
		printf("x(%d) = %lf\n", i, result[i]);
	}
	delete[] result;
#elif MATRIX_IN_USE == 2
	double** matrix = new double* [MATRIX_SIZE];
	InitMatrix(matrix);
	double* result = new double[MATRIX_SIZE];

	chrono::duration<double> elapsed_time = ParallelGaussMethod(matrix, MATRIX_SIZE, result);

	printf("Solution (elapsed time = %.8f seconds):\n", elapsed_time.count());
	for (i = 0; i < MATRIX_SIZE; ++i) {
		printf("x(%d) = %lf\n", i, result[i]);
	}

	delete[] result;
	DeInitMatrix(matrix);
	delete[] matrix;
#endif



	return 0;
}
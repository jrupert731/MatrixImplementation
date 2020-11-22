#include <iostream>
#include <fstream>

using namespace std;

class matrix {
private:
	int n;
	double *data;


public:
	//constructor
	matrix() {
		int i = 0;
		ifstream file;
		file.open("data.dat");
		file >> n;

		double a = 0;
		data = new double[n*n];
		while (file >> a) {
			data[i] = a;
			i++;
		}
		file.close();
	}

	//constructor for 0 matrix
	matrix(int n) {
		this->n = n;
		data = new double[n*n];
		for (int i = 0; i < n*n; i++) {
			data[i] = 0;
		}
	}
	
	//destructor
	~matrix() {
		delete[] data;
	}

	//copy constructor
	matrix(matrix &orig) {
		n = orig.n;
		data = new double[n];
		for (int i = 0; i < n*n; i++) {
			data[i] = orig.data[i];
		}
	}

	//move constructor
	matrix(matrix &&orig) {
		n = orig.n;
		data = orig.data;
		orig.data = nullptr;
	}

	//operator =
	matrix operator=(matrix &orig) {
		matrix copy = orig;
		swap(n, copy.n);
		swap(data, copy.data);
		return *this;
	}


	matrix sqr() const {
		matrix sqr(n);

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					sqr.data[n * j + i] += data[k * n + i] * data[j * n + k];
					//cout << data[k * n + i] << " x " << data[j * n + k] << " = " << sqr[n * j + i] << endl;
				}
				//cout << endl;
			}
		}
		return sqr;
	}



	matrix cube() const {
		matrix sqr(n);

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					sqr.data[n * j + i] += data[k * n + i] * data[j * n + k];
				}
			}
		}


		matrix cube(n);

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					cube.data[n * j + i] += sqr.data[k * n + i] * data[j * n + k];
					//cout << sqr[k * n + i] << " x " << data[j * n + k] << " = " << cube[n * j + i] << endl;
				}
				//cout << endl;
			}
		}
		return cube;
	}

	friend ostream &operator <<(ostream &os, const matrix &mat) {
		for (int i = 0; i < mat.n; i++){
			for (int j = 0; j < mat.n; j++) {
				os << mat.data[i * mat.n + j] << " ";
			}
			os << endl;
		}
		return os;
	}
};



int main() {
	matrix A;
	matrix B = A.sqr();
	matrix C = A.cube();

	cout << A << endl;
	cout << B << endl;
	cout << C << endl;

	return 0;
}
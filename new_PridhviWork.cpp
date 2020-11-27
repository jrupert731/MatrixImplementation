#define infile "data.dat"
#define debugging 1
#define read_from_cin 0
#define minimum_precision 1e-12
#include <memory.h>
#include <string.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
using std::cout;
using std::endl;
using std::ifstream;

void print_matrix(double* matrix, int n) {
  int counter = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++, counter++) {
      cout << matrix[counter];
      if (j != (n - 1)) {
        cout << "\t";
      }
    }
    if (i != (n - 1)) cout << "\n";
  }
  cout << endl;
}

void print_vector(double* vector, int n) {
  for (int i = 0; i < n; i++) {
    cout << vector[i];
    if (i != n - 1) {
      cout << "\t";
    }
  }
  cout << endl;
}

void print_answer(double* vector, int n) {
  for (int i = 0; i < n; i++) {
    cout << vector[i];
    if (i != n - 1) {
      cout << "\n";
    }
  }
  cout << endl;
}

void print_both(double*& m, double*& v, int& n) {
  cout << "Input Matrix"
       << "\n";
  print_matrix(m, n);
  cout << "Answer Vector"
       << "\n";
  print_vector(v, n);
}
using std::cin;
void get_input(int& n, double*& matrix, double*& answer_vector) {
  if (read_from_cin) {
    // ifstream dat_file(infile);
    // dat_file >> n;
    cin >> n;
    matrix = new double[n * n];
    for (int i = 0; i < (n * n); i++) {
      // dat_file >> matrix[i];
      cin >> matrix[i];
    }
    answer_vector = new double[n];
    for (int i = 0; i < n; i++) {
      // dat_file >> answer_vector[i];
      cin >> answer_vector[i];
    }
  } else {
    ifstream dat_file(infile);
    dat_file >> n;
    // cin >> n;
    matrix = new double[n * n];
    for (int i = 0; i < (n * n); i++) {
      dat_file >> matrix[i];
      // cin >> matrix[i];
    }
    answer_vector = new double[n];
    for (int i = 0; i < n; i++) {
      dat_file >> answer_vector[i];
      // cin >> answer_vector[i];
    }
  }
  if (debugging) {
    print_both(matrix, answer_vector, n);
  }
}

bool run_checks(int& n, double*& a, double*& v) {
  bool answer = (n > 0) && (a != nullptr) && (v != nullptr);
  if (answer && debugging) cout << "Input Basic Checks Passed." << endl;
  return answer;
}

bool row_is_zeros(const double* row, const int& n) {
  bool contains_only_zeros = true;
  for (int j = 0; j < n && contains_only_zeros; j++) {
    contains_only_zeros = contains_only_zeros && (row[j] == 0);
  }
  return contains_only_zeros;
}

int get_index(const int n, const int i, const int j) { return i * n + j; }
void get_nonzero_pivot(int& n, double*& m, double*& v, int& i) {
  int j = i + 1;
  while ((m[get_index(n, j, i)] == 0) && j < n) j++;
  j--;
  double scratch_space[n * n];
  double ans;
  memcpy(scratch_space, m + i * n, sizeof(double) * n);
  memcpy(m + i * n, m + j * n, sizeof(double) * n);
  memcpy(m + j * n, scratch_space, sizeof(double) * n);
  ans = v[i];
  v[i] = v[j];
  v[j] = ans;
  // delete[] scratch_space;
}

void partial_pivot(int& n, double*& m, double*& v, int i) {
  double max_value = 0;
  int max_index = i;
  for (int j = i; j < n; j++) {
    double val = fabs(m[get_index(n, j, i)]);
    if (val > max_value) {
      max_value = val;
      max_index = j;
    }
  }
  if (max_index != i) {
    // cout << "----- Partial Pivot ------" << endl;
    // print_both(m, v, n);
    double* scratch_space = new double[n];
    memcpy(scratch_space, m + i * n, sizeof(double) * n);
    memcpy(m + i * n, m + max_index * n, sizeof(double) * n);
    memcpy(m + max_index * n, scratch_space, sizeof(double) * n);
    double ans = v[i];
    v[i] = v[max_index];
    v[max_index] = ans;
    delete[] scratch_space;
    // print_both(m, v, n);
    // cout << "------ End Partial Pivot ------ " << endl;
  }
}

void convert_row_echelon_form(int& n, double*& m, double*& v) {
  for (int i = 0; i < n; i++) {
    // iterate through the pivots
    partial_pivot(n, m, v, i);
    // get_nonzero_pivot(n, m, v, i);
    int index_of_pivot = get_index(n, i, i);
    double pivot = m[index_of_pivot];
    if (pivot == 0) {
      // pivot = m[index_of_pivot];
      // cout << "unsolvable" << endl;
      // return;
      // pivot = m[index_of_pivot];
    }
    for (int j = i + 1; j < n; j++) {
      // iterate through the rows underneath each pivot
      int index_below_pivot = get_index(n, j, i);
      double& below_pivot = m[index_below_pivot];
      if (below_pivot != 0 && pivot != 0) {
        double reduction_factor = -below_pivot / pivot;
        // cout<<reduction_factor<<endl;
        m[get_index(n, j, i)] = 0;
        for (int k = i + 1; k < n; k++) {
          // iterate through each col in row to rescale it...
          double& rescaled = m[get_index(n, j, k)];
          double& orig = m[get_index(n, i, k)];
          rescaled = rescaled + reduction_factor * orig;
        }
        double& answer = v[j];
        double& orig_answer = v[i];
        answer = answer + reduction_factor * orig_answer;
      }
    }
  }
  // #if 0
  // Move the 0ed out rows to the end.
  int number_displaced = 0;
  int i = 0;
  double* scratch_space = new double[n];
  double ans;
  int last_row_with_zeros = n - 1;
  while (i < n) {
    bool contains_only_zeros = row_is_zeros(m + i * n, n);
    // contains_only_zeros = contains_only_zeros && (v[i] == 0);
    // unnecessary because if left side of eq is 0, then unless bad system
    // supplied, so will right side
    if (contains_only_zeros) {
      number_displaced += 1;

      int curr_row = i + 1;
      bool check_zeroes = true;
      while (check_zeroes == true && curr_row < n) {
        check_zeroes = row_is_zeros(curr_row * n + m, n);
        if (check_zeroes == false) {
          last_row_with_zeros = curr_row;
        } else {
          curr_row++;
        }
      }

      if (!check_zeroes) {
        int& row_being_displaced = curr_row;
        memcpy(scratch_space, m + i * n, sizeof(double) * n);
        memcpy(m + i * n, m + row_being_displaced * n, sizeof(double) * n);
        memcpy(m + row_being_displaced * n, scratch_space, sizeof(double) * n);
        ans = v[i];
        v[i] = v[row_being_displaced];
        v[row_being_displaced] = ans;
      } else {
        // The rest of the rows are 0, so it's okay to exit...
        i = n;
      }
    } else {
      i++;
    }
  }
  delete[] scratch_space;
}
void reduce_row_echelon_form(int& n, double*& m, double*& v) {
  // assume input is in row_echelon_form.;
  for (int i = n - 1; i >= 0; i--) {
    if (m[get_index(n, i, i)] != 0) {
      v[i] /= m[get_index(n, i, i)];
      m[get_index(n, i, i)] = 1;
    }
    if (fabs(v[i]) < minimum_precision) v[i] = 0;
    for (int j = 0; j < i; j++) {
      v[j] -= v[i] * m[get_index(n, j, i)];
      m[get_index(n, j, i)] = 0;
    }
  }
}
int main() {
  int n;
  double* matrix = nullptr;
  double* answer_vector = nullptr;
  get_input(n, matrix, answer_vector);
  if (run_checks(n, matrix, answer_vector)) {
    convert_row_echelon_form(n, matrix, answer_vector);
    if (debugging) print_both(matrix, answer_vector, n);
    reduce_row_echelon_form(n, matrix, answer_vector);
    if (debugging) cout << "----------------" << endl;
    if (debugging) print_both(matrix, answer_vector, n);
    print_answer(answer_vector, n);
    // here lie dragons
  } else {
    cout << "Invald data found, quitting." << endl;
  }
  delete[] matrix;
  delete[] answer_vector;
}
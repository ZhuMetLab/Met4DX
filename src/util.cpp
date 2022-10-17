#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector find_greater_equal_than(NumericVector data, NumericVector values) {
    int value_size = values.length();
    int in_size = data.length();

    int idx = 0;

    NumericVector res(value_size);
    cout << "value_size: " << value_size << endl;
    cout << "in_size: " << in_size << endl;
    for (int i = 0; i < value_size; i++) {
        while (idx < in_size && data[idx] < values[i])
            idx++;
        res[i] = idx + 1;
    }
    return(res);
}

// [[Rcpp::export]]
NumericVector find_dense_min(NumericVector values, int i_max) {
    int i, num_values = values.length();
    NumericVector res(2);

    for (i = i_max; i > 0; i--)
        if (values[i-1] > values[i])
            break;
    res[0] = i + 1;

    for (i = i_max; i < num_values-1; i++)
        if (values[i+1] > values[i])
            break;
    res[1] = i + 1;

    return(res);
}

// [[Rcpp::export]]
NumericVector rect_unique(NumericMatrix data, NumericVector order, NumericVector diff) {
    int i, j, io, jo;
    int num_row = data.nrow();
    int num_diff = diff.length();
    bool is_unique;

    if (2*num_diff != data.ncol()) {
        Rcout << "Error: number of columns in data must be twice the number of differences" << endl;
        throw std::exception();
        return NumericVector();
    }

    NumericVector keep(num_row);

    for (i = 0; i < num_row; i++) {
        io = order[i];
        keep[io] = 1;
        for (j = 0; j < i; j++) {
            jo = order[j];
            is_unique = 0;
            for (int k = 0; k < num_diff; k++) {
                is_unique = is_unique || data(io, 2*k) - data(jo, 2*k+1) > diff[k] || data(jo, 2*k) - data(io, 2*k+1) > diff[k];
            }
            if (keep[jo] && !is_unique) {
                keep[io] = 0;
                break;
            }
        }
    }

    return(keep);
}

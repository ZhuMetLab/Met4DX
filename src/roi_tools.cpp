#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
bool check_continuous_points_above_thr(NumericVector v, int i_start, int num, double thr, int n_skip) {
  int cnt = 0;
  bool res = false;
  int n_skip_cur = 0;
  int n_skip_pre = 0;
  for (int i = i_start; i < v.length(); i++) {
    if (v[i] > thr) {
      cnt++;
      n_skip_cur = 0;
      n_skip_pre = 0;
    } else {
      if (cnt > 0) {
        n_skip_pre = n_skip_cur;
        n_skip_cur++;
      } else {
        n_skip_cur = n_skip + 1;
      }
      if (n_skip_cur > n_skip) {
        cnt = 0;
        n_skip_cur = 0;
        n_skip_pre = 0;
      } else {
        if (n_skip_pre < n_skip) {
          cnt++;
        } else {
          cnt = cnt - n_skip_cur + 1;
          n_skip_pre = 0;
        }
      }
    }
    if (cnt >= num) {
      return(true);
    }
  }

  return(res);
}

// [[Rcpp::export]]
NumericVector get_continuous_points_above_thr_idx(NumericVector v, int i_start, int num, double thr, int n_skip) {
  int cnt = 0;
  int start_idx = 0;
  int end_idx = 0;
  int nv = v.length();
  int n_skip_cur = 0;
  int n_skip_pre = 0;
  NumericVector res(nv);
  for (int i = i_start; i < nv; i++) {
    res[i] = false;
    if (v[i] > thr) {
      cnt++;
      n_skip_cur = 0;
      n_skip_pre = 0;
      if (cnt == 1) {
        start_idx = i;
      } else {
        end_idx = i;
      }
    } else {
      if (cnt > 0) {
        n_skip_pre = n_skip_cur;
        n_skip_cur++;
      } else {
        n_skip_cur = n_skip + 1;
      }
      if (n_skip_cur > n_skip) {
        cnt = 0;
        n_skip_cur = 0;
        n_skip_pre = 0;
      } else {
        if (n_skip_pre < n_skip) {
          cnt++;
        } else {
          cnt = cnt - n_skip_cur + 1;
          n_skip_pre = 0;
        }
      }
    }

    if ((cnt == 0 || i == nv - 1) && (end_idx - start_idx + 1) >= num) {
      for (int j = start_idx; j <= end_idx; j++) {
        res[j] = true;
      }
      start_idx = 0;
      end_idx = 0;
    }
  }
  return(res);
}

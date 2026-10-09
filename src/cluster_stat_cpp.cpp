#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <cmath>
#include <limits>
#include <algorithm>

namespace {

class UnionFind {
public:
    std::vector<int> parent;
    std::vector<int> size;

    explicit UnionFind(int n) : parent(n), size(n, 1) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    int find(int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }

    void unite(int a, int b) {
        a = find(a);
        b = find(b);

        if (a == b) {
            return;
        }

        if (size[a] < size[b]) {
            std::swap(a, b);
        }

        parent[b] = a;
        size[a] += size[b];
    }
};

} // anonymous namespace

// [[Rcpp::export]]
SEXP cluster_stat_cpp(const Rcpp::NumericVector& data,
                      double nnodes,
                      double tcrit) {

    if (!std::isfinite(nnodes) ||
        nnodes < 1.0 ||
        nnodes != std::floor(nnodes) ||
        nnodes > std::numeric_limits<int>::max()) {
        Rcpp::stop("nnodes must be a positive integer.");
    }

    if (!std::isfinite(tcrit) || tcrit < 0.0) {
        Rcpp::stop("tcrit must be finite and nonnegative.");
    }

    const int n = static_cast<int>(nnodes);

    // Divide before multiplying to avoid unnecessary integer overflow.
    R_xlen_t a = static_cast<R_xlen_t>(n);
    R_xlen_t b = static_cast<R_xlen_t>(n - 1);

    if (a % 2 == 0) {
        a /= 2;
    } else {
        b /= 2;
    }

    if (b != 0 &&
        a > std::numeric_limits<R_xlen_t>::max() / b) {
        Rcpp::stop("nnodes is too large.");
    }

    const R_xlen_t expected = a * b;

    if (data.size() != expected) {
        Rcpp::stop(
            "length(data) must equal nnodes * (nnodes - 1) / 2."
        );
    }

    UnionFind uf(n);
    const double* values = data.begin();

    // Pass 1: connect suprathreshold edges.
    //
    // Match R's column-major upper.tri() ordering:
    // (1,2), (1,3), (2,3), (1,4), (2,4), (3,4), ...
    R_xlen_t k = 0;

    for (int col = 1; col < n; ++col) {
        Rcpp::checkUserInterrupt();

        for (int row = 0; row < col; ++row, ++k) {
            const double value = values[k];

            if (!std::isfinite(value)) {
                Rcpp::stop(
                    "data must contain finite t-statistics; "
                    "handle missing or undefined values before calling."
                );
            }

            const double magnitude = std::abs(value);

            // Zero values remain absent even when tcrit == 0.
            if (magnitude > 0.0 && magnitude >= tcrit) {
                uf.unite(row, col);
            }
        }
    }

    // Assign output rows by first encountered node.
    std::vector<int> root(n);
    std::vector<int> output_row(n, -1);
    int ncomponents = 0;

    for (int node = 0; node < n; ++node) {
        root[node] = uf.find(node);

        if (uf.size[root[node]] >= 3 &&
            output_row[root[node]] == -1) {
            output_row[root[node]] = ncomponents++;
        }
    }

    if (ncomponents == 0) {
        return Rcpp::NumericVector::create(0.0, 0.0);
    }

    std::vector<double> edge_count(ncomponents, 0.0);
    std::vector<long double> edge_mass(ncomponents, 0.0L);

    // Pass 2: sum retained edges once, without an adjacency matrix.
    k = 0;

    for (int col = 1; col < n; ++col) {
        Rcpp::checkUserInterrupt();

        for (int row = 0; row < col; ++row, ++k) {
            const double magnitude = std::abs(values[k]);

            if (magnitude > 0.0 && magnitude >= tcrit) {
                const int out = output_row[root[row]];

                if (out >= 0) {
                    edge_count[out] += 1.0;
                    edge_mass[out] +=
                        static_cast<long double>(magnitude - tcrit);
                }
            }
        }
    }

    Rcpp::NumericMatrix result(ncomponents, 2);

    for (int component = 0; component < ncomponents; ++component) {
        result(component, 0) = edge_count[component];
        result(component, 1) =
            static_cast<double>(edge_mass[component]);
    }

    return result;
}
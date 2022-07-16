#include <cstdio>
#include <algorithm>
#include <vector>
using namespace std;

// We assume that vEBTree has size of 2^(2^k) - otherwise implementation is bit more messy.

class vEBTree {
public:
    int u_bits, _min, _max;
    vEBTree *summary;
    vEBTree **clusters;
    int low(int x) {
        return x & ((1 << (u_bits/2)) - 1);
    }
    int high(int x) {
        return (x - low(x))>>(u_bits/2);
    }
    int index(int i, int j) {
        return (i << (u_bits/2)) + j;
    }
    vEBTree(int _u = 1) : u_bits(_u), _min(-1), _max(-1) {
        if (_u == 1) {
            summary = NULL; clusters = NULL;
            return;
        }
        int sub_u = u_bits/2;
        summary = new vEBTree(sub_u);
        clusters = new vEBTree*[1<<sub_u];
    };
    void insert(int x);
    void remove(int x);
    int successor(int x) {
        return 0;
    };
    int predecessor(int x) {
        return 0;
    };
    void print() {
        int x = _min;
        if (x == -1) {
            printf("EMPTY\n"); return;
        }
        else {
            printf("%d",x);
            while (x != -1) {
                x = successor(x);
                printf("-> %d", x);
            }
        }
    }
};

void vEBTree::insert(int x){
    if (u_bits == 1) {
        if (x == 0) {
            _min = 0; _max = (_max==-1 ? 0 : -1);
        }
        else {
            _max = 1; _min = (_min==-1 ? 1 : -1);
        }
    }
    else {
        if (_min == -1) {
            _min = _max = x;
            return;
        }
        if (x < _min) {
            swap(x, _min);
            return insert(x);
        }
        if (x > _max) _max = x;
        int hi = high(x), lo = low(x);
        if (clusters[hi] == NULL) {
            clusters[hi] = new vEBTree(u_bits/2);
            summary -> insert(hi);
        }
        clusters[hi] -> insert(lo);
    }
}

void vEBTree::remove(int x) {
    if (u_bits == 1) {
        if (x == 0) {
            if (_max == 0) _min = _max = -1;
            else if (_max == 1) _min = 1;
        }
        else {
            if (_min == 0) _max = 0;
            else if (_min == 1) _min = _max = -1;
        }
    }
    else {
        if (x == _min) {
            int i = summary -> _min;
            if (i == -1) {
                _min = _max = -1;
                return;
            }
            x = _min = index(i, clusters[i]->_min);
        }
        clusters[high(x)] -> remove(low(x));
        if (clusters[high(x)] -> _min == -1) {
            summary -> remove(high(x));
        }
        if (x == _max) {
            if (summary -> _max == -1) _max = _min;
            else {
                int i = summary -> _max;
                _max = index(i, clusters[i] -> _max);
            }
        }
    }
}

int main() {
    vEBTree *vEB = new vEBTree(16);
    vEB -> insert(2);
    vEB -> insert(3);
    vEB -> insert(4);
    vEB -> insert(5);
    vEB -> insert(7);
    vEB -> insert(14);
    vEB -> insert(15);

    printf("%d\n", vEB -> _min);
    printf("%d\n", vEB -> _max);
    printf("%d\n", (vEB -> predecessor(9)));
    printf("%d\n", (vEB -> successor(9)));

    vEB -> remove(2);
    vEB -> remove(15);
    vEB -> remove(7);

    printf("%d\n", vEB -> _min);
    printf("%d\n", vEB -> _max);

    printf("%d\n", (vEB -> predecessor(9)));
    printf("%d\n", (vEB -> successor(9)));

    return 0;
}
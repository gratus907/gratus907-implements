/*
 * Implementation of Rabin-Karp hashing (Supports multiple-hash)
 * Tested : BOJ 1786
 */
#include <vector>
#include <iostream>
#include <cassert>
#define ll long long
using namespace std;

vector <ll> RabinKarp(string &s, ll window, ll mod, ll p) {
    vector <ll> hash_values;
    ll cur = 0;
    ll hash_offset = 1;
    for (int i = 0; i < window; i++) {
        cur *= p;
        cur += s[i];
        cur %= mod;
        if (i > 0) hash_offset = (hash_offset * p)%mod;
    }
    hash_values.push_back(cur);
    for (int i = window; i < s.length(); i++) {
        cur = (mod + cur - (s[i-window] * hash_offset)%mod)%mod;
        cur = (cur * p + s[i]) % mod;
        hash_values.push_back(cur);
    }
    return hash_values;
}

vector <pair<ll, ll>> configs = {{1000000007, 31}, {998244353, 37}, {1000000009, 59}};
vector <bool> matches;
vector <vector<ll>> hashes;
vector <ll> pattern_hashes;

int main() {
    string T, P;
    getline(cin, T);
    getline(cin, P);
    if (T.length() < P.length()) {
        cout << 0 << '\n';
        return 0;
    }
    matches.resize(T.length(), true);
    for (auto p : configs) {
        hashes.push_back(RabinKarp(T, P.length(), p.first, p.second));
        pattern_hashes.push_back(RabinKarp(P, P.length(), p.first, p.second)[0]);
    }
    for (int j = 0; j < T.length()-P.length()+1; j++) {
        for (ll i = 0; i < configs.size(); i++) {
            matches[j] = matches[j] & (hashes[i][j] == pattern_hashes[i]);
        }
    }
    cout << count(matches.begin(), matches.end(), true) << '\n';
    for (int i = 0; i < T.length()-P.length()+1; i++) {
        if (matches[i]) cout << i+1 << ' ';
    }
    cout << '\n';
}
#include <stdio.h>
#include <string.h>

#define wipe(v, x) (memset((v), (x), sizeof(v)))
#define MAXN 100
#define MAXL 20

int N, A;
char heroi[MAXN][MAXL+1];
char ident[MAXN][MAXL+1];

int st[MAXN];
int inc[MAXN][MAXN];

int adj[MAXN][MAXN], nadj[MAXN];
int map[MAXN];

int marc[MAXN];
char buf[MAXL+1];

int dfs_p(int p) {
	int i, viz;
	if (marc[p]) {
		return 0;
	}
	marc[p] = 1;
	for (i = 0; i < nadj[p]; i++) {
		viz = adj[p][i];
		if (map[viz] == -1 || dfs_p(map[viz])) {
			map[viz] = p;
			return 1;
		}
	}
	return 0;
}

int cycle_p(int p, int dest) {
	int i, no, viz;
	if (marc[p]) {
		return 0;
	}
	marc[p] = 1;
	no = map[p];
	if (no == dest) {
		return 1;
	}
	for (i = 0; i < nadj[no]; i++) {
		viz = adj[no][i];
		if (cycle_p(viz, dest)) {
			return 1;
		}
	}
	return 0;
}

int main() {
	int i, j, k, M, a, b;
	int fail, viz;
	scanf(" %d %d", &N, &A);
	/* read names */

	for (i = 0; i < N; i++) {
		scanf(" %s", heroi[i]);
	}
	for (i = 0; i < N; i++) {
		scanf(" %s", ident[i]);
	}

	/* build graph */
	wipe(inc, 0);
	for (k = 0; k < A; k++) {
		scanf(" %d", &M);
		for (j = 0; j < M; j++) {
			scanf(" %s", buf);
			for (i = 0; i < N; i++) {
				if (strcmp(heroi[i], buf) == 0) {
					break;
				}
			}
			if (i == N) {
				for (i = 0; i < N; i++) {
					if (strcmp(ident[i], buf) == 0) {
						break;
					}
				}
				i = ~i;
			}
			st[j] = i;
		}
		for (i = 0; i < M; i++) {
			for (j = 0; j < M; j++) {
				a = st[i];
				b = st[j];
				if ((a ^ b) < 0) {
					if (a < 0) {
						inc[b][~a] = 1;
					} else {
						inc[a][~b] = 1;
					}
				}
			}
		}
	}

	wipe(nadj, 0);
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (!inc[i][j]) {
				adj[i][nadj[i]++] = j;
			}
		}
	}

	/* matching */
	fail = 0;
	wipe(map, -1);
	for (i = 0; i < N; i++) {
		wipe(marc, 0);
		if (!dfs_p(i)) {
			fail = 1;
			break;
		}
	}
	if (fail) {
		printf("IMPOSSIVEL\n");
		return 0;
	}

	/* possible matchings */
	for (i = 0; i < N; i++) {
		printf("%s:", heroi[i]);
		for (j = 0; j < nadj[i]; j++) {
			viz = adj[i][j];
			wipe(marc, 0);
			if (cycle_p(viz, i)) {
				printf(" %s", ident[viz]);
			}
		}
		printf("\n");
	}
	
	return 0;
}

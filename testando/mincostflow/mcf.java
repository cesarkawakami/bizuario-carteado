class Edge
{
    private int pv, pw, pcap, pflow;
    Edge(int v, int w, int cap)
    { pv = v; pw = w; pcap = cap; pflow = 0; }
    int v() { return pv; }
    int w() { return pw; }
    int cap() { return pcap; }
    int flow() { return pflow; }
    boolean from (int v)
    { return pv == v; }
    int other(int v)
    { return from(v) ? pw : pv; }
    int capRto(int v)
    { return from(v) ? pflow : pcap - pflow; }
    void addflowRto(int v, int d)
    { pflow += from(v) ? -d : d; }
}

class NetworkMinCost
{ 
    private Network G; private int s, t, valid;
    private Edge[] st;
    private int[] mark, phi;
    private int ST(int v) { return st[v].other(v); }
    private void dfsR(Edge x, int v);
    // See Exercise 22.117
    private int phiR(int v)
    {
	if (mark[v] == valid) return phi[v];
	phi[v] = phiR(ST(v)) - st[v].costRto(v);
	mark[v] = valid;
	return phi[v];
    }  
    private int lca(int v, int w)
    { 
	mark[v] = ++valid; mark[w] = valid;
	while (v != w)
	    {
		if (v != t) v = ST(v);
		if (v != t && mark[v] == valid) return v;
		mark[v] = valid;
		if (w != t) w = ST(w);
		if (w != t && mark[w] == valid) return w;
		mark[w] = valid;
	    }
	return v;
    }
    private Edge augment(Edge x)
    { 
	int v = x.v(), w = x.w(); int r = lca(v, w);
	int d = x.capRto(w);
	for (int u = w; u != r; u = ST(u))
	    if (st[u].capRto(ST(u)) < d)
		d = st[u].capRto(ST(u));
	for (int u = v; u != r; u = ST(u))
	    if (st[u].capRto(u) < d)
		d = st[u].capRto(u);
	x.addflowRto(w, d); Edge e = x;
	for (int u = w; u != r; u = ST(u))
	    { st[u].addflowRto(ST(u), d);
		if (st[u].capRto(ST(u)) == 0) e = st[u]; }
	for (int u = v; u != r; u = ST(u))
	    { st[u].addflowRto(u, d);
		if (st[u].capRto(u) == 0) e = st[u]; }
	return e;
    }
    private boolean onpath(int a, int b, int c)
    {
	for (int i = a; i != c; i = ST(i))
	    if (i == b) return true;
	return false;
    }
    private void reverse(int u, int x)
    { 
	Edge e = st[u];
	for (int i = ST(u); i != x; i = ST(i))
	    { Edge y = st[i]; st[i] = e; e = y; }
    }
    private void subs(Edge w, Edge y)
    { 
	int u = y.w(), v = y.v(), x = w.w();
	if (st[x] != w) x = w.v();
	int r = lca(u, v);
	if (onpath(u, x, r))
	    { reverse(u, x); st[u] = y; return; }
	if (onpath(v, x, r))
	    { reverse(v, x); st[v] = y; return; }
    }
    private int costR(Edge e, int v)
    { 
	int R = e.cost() + phiR(e.w()) - phiR(e.v());
	return e.from(v) ? R : -R; }
    private Edge besteligible()
    { 
	Edge x = null; int V = G.V();
	for (int v = 0, min = Edge.C*V; v < V; v++)
	    {
		AdjList A = G.getAdjList(v);
		for (Edge e = A.beg(); !A.end(); e = A.nxt())
		    if (e.capRto(e.other(v)) > 0)
			if (e.capRto(v) == 0)
			    if (costR(e, v) < min)
				{ x = e; min = costR(e, v); }
	    }
	return x;
    }
    NetworkMinCost(Network G, int s, int t)
    { 
	int V = G.V();
	this.G = G; this.s = s; this.t = t;
	st = new Edge[V];
	mark = new int[V]; phi = new int[V];
	for (int i = 0; i < V; i++) mark[i] = -1;
	Edge z = new Edge(s, t, Edge.M*V, Edge.C*V);
	G.insert(z); z.addflowRto(t, z.cap());
	dfsR(z, t);
	int old = 0;
	for (valid = 1; valid != old; ) {
	    old = valid;
	    for (int v = 0; v < G.V(); v++) {
		AdjList A = G.getAdjList(v);
		for (Edge e = A.beg(); !A.end(); e = A.nxt())
		    if (e.capRto(e.other(v)) > 0)
			if (e.capRto(v) == 0)
			    { subs(augment(e), e); valid++; }
	    }
	}
	G.remove(z);
    }
}


Let $n$ be the number of nodes, $m$ the number of links and $k$ the average degree.

$$2m=nk$$

Let $m'$ be the number of links and $k'$ the average degree in the reduced network, and $p$ the probability of a link.

$$m'=pm$$

According to the percolation theory, a giant connected component emerges at $k'=1$. Thus,

$$2m'=nk'$$
$$2pm=n$$
$$p=1/k$$

The value $p=1/k$ is drawn on the plots in the paper and corresponds to the intuition that every paper has at least one relevant reference or citation.

Assuming that every paper has at least one relevant reference and one relevant citation, the value of $p$ should be based on the condition $k_{out}=k_{in}=1$, where $k_{out}$ and $k_{in}$ are the average out-degree and in-degree. Since $k_{out}=k_{in}=k/2$ for any network, the corresponding probability should then be $p=2/k$.
## Problem Statement
Given:
- Bipartite Graph G(V = X ∪ Y, E)
- For every edge e 𝞊 E, Ce is the cost of the edge.
- |X| = |Y| = n
- There exist a perfect matching.
- Ce ≥ 0 ∀ e 𝞊 E.

Goal:
We need to find the minimum cost bipartite matching.


## Algorithm
```
Start with M equal to the empty set
Define p(x) = 0 for x 𝛜 X, and p(y) = min(e) into y ce for y 𝛜 Y
While M is not a perfect matching..,
	Compute shortest path distances d
P ← shortest alternating path using reduced cost & Dijkstra’s Algo.
Augment along P to produce a new matching M’.
foreach v ∈ X ∪ Y: p(v) ← p(v) + d(v).
Endwhile;
```

## Code Written By
- [Amit Srivastava](https://github.com/cynicphoenix)
- Gyan Prakash Singh

import java.util.*;

// Class to represent a node in the graph 
class Node implements Comparator<Node> { 
	public int node; 
    public int cost;

	public Node() { 
	} 

	public Node(int node, int cost) { 
		this.node = node; 
        this.cost = cost;
	} 

	@Override
	public int compare(Node node1, Node node2) { 
		if (node1.cost < node2.cost) 
			return -1; 
		if (node1.cost > node2.cost) 
			return 1; 
		return 0; 
	} 
}

// Class to represent edge for matching
class edge {
    int start;
    int end;

    public edge(int start, int end) {
        this.start = start;
        this.end = end;
    }
}

// Represent graph
class Graph {
    int V, P, Q, source, sink;
    List<List<Node>> adj;
    List<Integer> pSet, qSet;

    public Graph(int V, int P, int Q, int source, int sink) {
        this.V = V;
        this.P = P;
        this.Q = Q;
        this.source = source;
        this.sink = sink;
        adj = new ArrayList<List<Node>>();
        pSet = new ArrayList<Integer>();
        qSet = new ArrayList<Integer>();
    }
}


// Dijkstra Algo to find minimum path using min Heap: Modified version taken from Geeksforgeeks
class DPQ { 
    int dist[];
    int parent[];
	Set<Integer> settled; 
	PriorityQueue<Node> pq; 
    int V;
    int source;
    int sink; 
	List<List<Node>> adj; 

	public DPQ(int V, int source, int sink, List<List<Node>> adj) { 
        this.V = V;
        this.source = source;
        this.sink = sink;
        dist = new int[V];
        parent = new int[V]; 
		settled = new HashSet<Integer>(); 
        pq = new PriorityQueue<Node>(V, new Node()); 
        this.adj = adj;
	} 

	// Function for Dijkstra's Algorithm 
	public void dijkstra() {

		for (int i = 0; i < V; i++) 
			dist[i] = Integer.MAX_VALUE; 

		// Add source node to the priority queue 
		pq.add(new Node(source, 0)); 
        dist[source] = 0;
        parent[source] =- 1;
		while (settled.size() != V && !pq.isEmpty()) { 
			// remove the minimum distance node 
			// from the priority queue 
			int u = pq.remove().node; 
			// adding the node whose distance is 
			// finalized 
            settled.add(u); 
            if (u != sink)
			    e_Neighbours(u); 
		} 
	}

	// Function to process all the neighbours 
	// of the passed node 
	private void e_Neighbours(int u) { 
		int edgeDistance = -1; 
		int newDistance = -1; 

		// All the neighbors of v 
		for (int i = 0; i < adj.get(u).size(); i++) { 
			Node v = adj.get(u).get(i); 

			// If current node hasn't already been processed 
			if (!settled.contains(v.node)) { 
				edgeDistance = v.cost; 
				newDistance = dist[u] + edgeDistance; 

				// If new distance is cheaper in cost 
				if (newDistance < dist[v.node]) {
                    dist[v.node] = newDistance; 
                    parent[v.node] = u;
                }

				// Add the current node to the queue 
				pq.add(new Node(v.node, dist[v.node])); 
			} 
		} 
    
    }
} 


class sspa {
    public static void main(String[] args) {
        // Create Original Graph
        int V = 12, P = 5, Q = 5, source = 0, sink = 11;
        Graph originalGraph = new Graph(V, P, Q, source, sink);

		for (int i = 0; i < V; i++) { 
			List<Node> item = new ArrayList<Node>(); 
			originalGraph.adj.add(item); 
        }
        originalGraph.pSet.add(1);
        originalGraph.pSet.add(2);
        originalGraph.pSet.add(3);
        originalGraph.pSet.add(4);
        originalGraph.pSet.add(5);

        originalGraph.qSet.add(6);
        originalGraph.qSet.add(7);
        originalGraph.qSet.add(8);
        originalGraph.qSet.add(9);
        originalGraph.qSet.add(10);
    
        // Input from user
        originalGraph.adj.get(1).add(new Node(6, 1)); 
        originalGraph.adj.get(1).add(new Node(7, 2));
        originalGraph.adj.get(1).add(new Node(8, 2)); 
        originalGraph.adj.get(1).add(new Node(9, 2));
        originalGraph.adj.get(1).add(new Node(10, 2));

        originalGraph.adj.get(2).add(new Node(6, 2)); 
        originalGraph.adj.get(2).add(new Node(7, 2));
        originalGraph.adj.get(2).add(new Node(8, 2)); 
        originalGraph.adj.get(2).add(new Node(9, 1));
        originalGraph.adj.get(2).add(new Node(10, 2));

        
        originalGraph.adj.get(3).add(new Node(6, 2)); 
        originalGraph.adj.get(3).add(new Node(7, 2));
        originalGraph.adj.get(3).add(new Node(8, 2)); 
        originalGraph.adj.get(3).add(new Node(9, 2));
        originalGraph.adj.get(3).add(new Node(10, 1));

        
        originalGraph.adj.get(4).add(new Node(6, 2)); 
        originalGraph.adj.get(4).add(new Node(7, 1));
        originalGraph.adj.get(4).add(new Node(8, 2)); 
        originalGraph.adj.get(4).add(new Node(9, 2));
        originalGraph.adj.get(4).add(new Node(10, 2));

        
        originalGraph.adj.get(5).add(new Node(6, 2)); 
        originalGraph.adj.get(5).add(new Node(7, 2));
        originalGraph.adj.get(5).add(new Node(8, 1)); 
        originalGraph.adj.get(5).add(new Node(9, 2));
        originalGraph.adj.get(5).add(new Node(10, 2));
        
        // Add source sink
        originalGraph.adj.get(0).add(new Node(1, 0)); 
        originalGraph.adj.get(0).add(new Node(2, 0));
        originalGraph.adj.get(0).add(new Node(3, 0));
        originalGraph.adj.get(0).add(new Node(4, 0)); 
        originalGraph.adj.get(0).add(new Node(5, 0));

        originalGraph.adj.get(6).add(new Node(11, 0)); 
        originalGraph.adj.get(7).add(new Node(11, 0));
        originalGraph.adj.get(8).add(new Node(11, 0)); 
        originalGraph.adj.get(9).add(new Node(11, 0));
        originalGraph.adj.get(10).add(new Node(11, 0));

        List<edge> matching = new ArrayList<edge>();
        int[] potential = new int[V];

        // Initialize Potential
        for(int i = 0; i < V; i++) {
            if(i == source || i == sink)
                potential[i] = 0;
            else {
                if(originalGraph.pSet.contains(i))
                    potential[i] = 0;
                else {
                    int min = Integer.MAX_VALUE;
                    for(int j = 0; j < originalGraph.adj.size(); j++)
                        for(int k = 0; k < originalGraph.adj.get(j).size(); k++)
                            if(originalGraph.adj.get(j).get(k).node == i && originalGraph.adj.get(j).get(k).cost < min)
                                min = originalGraph.adj.get(j).get(k).cost;
                    potential[i] = min;
                }
            }
        }

        // Create a duplicate list which will be updated in each iteration
        List<List<Node>> adj = new ArrayList<List<Node>>();
        for (int i = 0; i < V; i++) { 
            List<Node> item = new ArrayList<Node>(); 
            adj.add(item);
            for(int j = 0; j < originalGraph.adj.get(i).size(); j++)
                adj.get(i).add(new Node(originalGraph.adj.get(i).get(j).node, originalGraph.adj.get(i).get(j).cost)); 
        }

        // Main Loop
        while(matching.size() < P) {
            // Make new graph with reduced edges
            for (int i = 0; i < adj.size(); i++) { 
                for(int j = 0; j < adj.get(i).size(); j++) {
                    int end = adj.get(i).get(j).node;
                    int p = i, q = end;
                    if(originalGraph.qSet.contains(i)) {
                        p = end; q = i;
                    }
                    int reduced_cost;
                    if(p == source || p == sink || q == source || q == sink)
                        reduced_cost = 0;
                    else {
                        int cost = 0;
                        for(int k = 0; k < originalGraph.adj.get(p).size(); k++)
                            if(originalGraph.adj.get(p).get(k).node == q)
                                cost = originalGraph.adj.get(p).get(k).cost;
                        reduced_cost = cost + potential[p] - potential[q];
                    }
                    adj.get(i).set(j, new Node(end, reduced_cost)); 
                }
            }

            // Apply Djikstra to find augmenting path
            DPQ dpq = new DPQ(V, source, sink, adj);
            dpq.dijkstra();

            // Find shortest path
            List <edge> path = new ArrayList<edge>();
            int next_node = sink;
            while(dpq.parent[next_node] != -1) {
                path.add(new edge(dpq.parent[next_node], next_node));
                next_node = dpq.parent[next_node];
            }

            // Reverse edges
            for(int i = 0; i < adj.size(); i++) {
                for(int j = 0; j < adj.get(i).size(); j++) {
                    for(int k = 0; k < path.size(); k++){
                        if(path.get(k).start == i && path.get(k).end == adj.get(i).get(j).node) {
                            int node = adj.get(i).get(j).node;
                            int cost = adj.get(i).get(j).cost;
                            adj.get(i).remove(j);
                            adj.get(node).add(new Node(i, cost));
                        }
                    }
                }
            }

            // Update matching
            for(int i = 0; i < path.size(); i++) {
                if(path.get(i).start == source || path.get(i).end == source || path.get(i).start == sink || path.get(i).end == sink)
                    continue;
                boolean flag = false;
                for(int j = 0; j < matching.size(); j++) {
                    if(matching.get(j).start == path.get(i).start && matching.get(j).end == path.get(i).end) {
                        matching.remove(j);
                        flag = true;
                    }
                    if(matching.get(j).start == path.get(i).end && matching.get(j).end == path.get(i).start) {
                        matching.remove(j);
                        flag = true;
                    }
                }
                if(flag == false)
                    matching.add(path.get(i));
            }

            // Update Potential
            for(int i = 0; i < V; i++)
                if(i != source && i != sink)
                    potential[i] = potential[i] + dpq.dist[i];
        }

        // Print matching result
        for(int i = 0; i < matching.size(); i++)
            System.out.println(matching.get(i).start + " " + matching.get(i).end);

    }
}

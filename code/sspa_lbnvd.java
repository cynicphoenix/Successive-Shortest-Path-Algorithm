import java.util.*;
import java.io.File;  // Import the File class
import java.io.FileNotFoundException;  // Import this class to handle errors
import java.util.Scanner; 

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

class ServiceCenter {
    long node;
    int capacity;
    int penalty;
    ServiceCenter (long node, int capacity, int penalty) {
        this.node = node;
        this.capacity = capacity;
        this.penalty = penalty;
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
        int totalServiceCenters = 0, totalDemandVertex = 0;
        
        List<Long> demandVertexList = new ArrayList<Long>();
        List<Long> serviceCenterList = new ArrayList<Long>();
        List<ServiceCenter> serviceCenterListData = new ArrayList<ServiceCenter>();        
        // List<datasetEdge> edgesListData = new ArrayList<datasetEdge>();
        List<Integer> serviceCenterPosition = new ArrayList<Integer>();       
        
        try {
            File ServiceCentersFile = new File("ServiceCenter.txt");
            Scanner myReader = new Scanner(ServiceCentersFile);
            while (myReader.hasNextLine()) {
                String data = myReader.nextLine();
                String[] tokens = data.split(",");
                ServiceCenter temp = new ServiceCenter(Long.parseLong(tokens[0]), Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]));
                serviceCenterListData.add(temp);
                serviceCenterList.add(Long.parseLong(tokens[0]));
                totalServiceCenters++;
            }
            myReader.close();
            System.out.println("Service Centers : "+totalServiceCenters);
        } catch (FileNotFoundException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }

        try {
            File NodesFile = new File("nodes.txt");
            Scanner myReader = new Scanner(NodesFile);
            int i = 0;
            while (myReader.hasNextLine()) {
                long temp = Long.parseLong(myReader.nextLine());
                if(!serviceCenterList.contains(temp)) {
                    demandVertexList.add(temp);
                    totalDemandVertex++;
                }
                else {
                    serviceCenterPosition.add(i);
                }
                i++;
            }
            myReader.close();
            System.out.println("Demand Vertex : " + totalDemandVertex);
        } catch (FileNotFoundException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }

        int[][] shortestDistance = new int[totalDemandVertex][totalServiceCenters];
        try {
            File CostMatrixFile = new File("CostMatrix.txt");
            Scanner myReader = new Scanner(CostMatrixFile);
            int i = 0, j = 0;
            while (myReader.hasNextLine()) {
                String data = myReader.nextLine();
                if(!serviceCenterPosition.contains(i)) {
                    String[] tokens = data.split(",");
                    for(int k = 0; k < totalServiceCenters; k++)
                        shortestDistance[j][k] = Integer.parseInt(tokens[k]);
                    j++;
                }
                i++;
            }
            myReader.close();
        } catch (FileNotFoundException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }

        int P = totalServiceCenters * totalDemandVertex;
        int Q = P;
        int V = P + Q + 2;
        int source = 0, sink = V - 1;

        // Service -> s1 : 1 - d, s2 : d+1 - 2d, ..., ss = .. - sd
        // Demand -> d1 : sd + 1, d2 : sd + 2, ..., dd = sd + d
        // Dummy -> sd + d + 1 ,.., sink - 1.

        Graph originalGraph = new Graph(V, P, Q, source, sink);
        for (int i = 0; i < V; i++) { 
			List<Node> item = new ArrayList<Node>(); 
			originalGraph.adj.add(item); 
        }

        // Add P side
        for(int  i = 1; i <= P; i++)
            originalGraph.pSet.add(i);

        // Add Q side
        for(int  i = P + 1; i < V - 1; i++)
            originalGraph.qSet.add(i);

        // Add source edges
        for(int i = 1; i <= P; i++)
            originalGraph.adj.get(0).add(new Node(i, 0));
   
        // Add sink edges
        for(int  i = P + 1; i < V - 1; i++)
            originalGraph.adj.get(i).add(new Node(sink, 0));
       

        // Add rest edges
        for(int demandVertex = P + 1; demandVertex <= P + totalDemandVertex; demandVertex++) {
            for(int serviceCenter = 0; serviceCenter < totalServiceCenters; serviceCenter++) {
                int distance = shortestDistance[demandVertex - P - 1][serviceCenter];
                int capacity = serviceCenterListData.get(serviceCenter).capacity;
                int penalty = serviceCenterListData.get(serviceCenter).penalty;
                for(int i = 1; i <= totalDemandVertex; i++) {
                    
                    
                    if (i <= capacity)
                        originalGraph.adj.get(serviceCenter * totalDemandVertex + i).add(new Node(demandVertex, distance));
                    else
                        originalGraph.adj.get(serviceCenter * totalDemandVertex + i).add(new Node(demandVertex, distance + penalty));
                }
            }
        }

        // Add dummy vertex edges
        for(int i = 1; i <= P; i++)
            for(int j = P + totalDemandVertex + 1; j < V - 1; j++) {
               
                originalGraph.adj.get(i).add(new Node(j, Integer.MAX_VALUE));
            }
        System.out.println("Graph Created");
        // Print Graph
        // for (int i = 0; i < V; i++) {
        //     for(int j = 0; j < originalGraph.adj.get(i).size(); j++) {
        //         System.out.println(i + " " + originalGraph.adj.get(i).get(j).node + " " + originalGraph.adj.get(i).get(j).cost);
        //     }
        // }
        
        List<edge> matching = new ArrayList<edge>();
        int[] potential = new int[V];
        int [][] costOriginal = new int [V][V];


        // Initialize Potential
        for(int i = 0; i < V; i++) {
            
            if(i == source || i == sink)
                potential[i] = 0;
            else {
                if(originalGraph.pSet.contains(i))
                    potential[i] = 0;
                else {
                    potential[i]=Integer.MAX_VALUE;
                }
            }
        }
        for(int j = 0; j < originalGraph.adj.size(); j++)
          {
            for(int k = 0; k < originalGraph.adj.get(j).size(); k++)
              {

                  int nd =originalGraph.adj.get(j).get(k).node;
                  int cst = originalGraph.adj.get(j).get(k).cost ;

                 /////  To store the original cost 
                  costOriginal[j][nd]=cst;

                  if(potential[nd]>cst)
                     potential[nd]=cst;
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
        System.out.println("Duplicate Graph Created");
        
        
        final long startTime = System.currentTimeMillis();
        // Main Loop
        while(matching.size() < P) {

            // Make new graph with reduced edges
            for (int i = 0; i < adj.size(); i++) { 
                for(int j = 0; j < adj.get(i).size(); j++) { 
                    int end = adj.get(i).get(j).node;
                    int p = i, q = end;
                    
                    if( p>P) {
                        p = end; q = i;
                    }
                    
                    int reduced_cost;
                    if(p == source || p == sink || q == source || q == sink)
                        reduced_cost = 0;
                    else {
                        int cost = 0;
                        // for(int k = 0; k < originalGraph.adj.get(p).size(); k++)
                        //     if(originalGraph.adj.get(p).get(k).node == q)
                               // cost =   originalGraph.adj.get(p).get(k).cost;
                        cost = costOriginal[p][q];
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
            // for(int i = 0; i < adj.size(); i++) {
            //     for(int j = 0; j < adj.get(i).size(); j++) {
            //         for(int k = 0; k < path.size(); k++){
            //             if(path.get(k).start == i && path.get(k).end == adj.get(i).get(j).node) {
            //                 int node = adj.get(i).get(j).node;
            //                 int cost = adj.get(i).get(j).cost;
            //                 adj.get(i).remove(j);adj.get(str).remove(j);
            //                 adj.get(node).add(new Node(i, cost));
            //             } 
            //         }
            //     }
            // }

            for(int k = 0; k < path.size(); k++){
                  int str = path.get(k).start ;
                  int ending = path.get(k).end ;
                for( int loop =0;loop<adj.get(str).size(); loop++)
                {
                    if( ending == adj.get(str).get(loop).node) {
                        int cost = adj.get(str).get(loop).cost;
                        adj.get(str).remove(loop);
                        adj.get(ending).add(new Node( str, cost));
                        break;
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
        final long endTime = System.currentTimeMillis();
        System.out.println("Total execution time: " + (endTime - startTime));
        System.out.println("Matching: ");

        // Print matching result
        for(int i = 0; i < matching.size(); i++) {
            if(matching.get(i).end <= P + totalDemandVertex) {
                long start = serviceCenterListData.get(matching.get(i).start/totalDemandVertex).node;
                long end = demandVertexList.get(matching.get(i).end - P - 1);
                System.out.println(start + " " + end);
            }
        }
    }
}
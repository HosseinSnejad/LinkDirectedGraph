
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Objects;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import javax.swing.text.html.HTMLDocument.Iterator;

public class LinkDirectedGraph {

	private static final int K_PART_SIZE = 4;
    private static final int K_MER_SIZE = 7;
    private static final int edgeSupport = 5;
    private static final int pathSupport = 5;
	public static char [] AminoAcid = new char [] {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'Z'};
	 public static class KPart {
	        public String amino;
	        int[] gaps;

	        public KPart(String amino, int[] gaps) {
	            this.amino = amino;
	            this.gaps = gaps;
	        }

	        @Override
	        public String toString() {
	            return "KPart{" + "amino='" + amino + '\'' +
	                    ", gaps=" + Arrays.toString(gaps) +
	                    '}';
	        }

	        @Override
	        public boolean equals(Object o) {
	            if (this == o) return true;
	            if (o == null || getClass() != o.getClass()) return false;
	            KPart kPart = (KPart) o;
	            return Objects.equals(amino, kPart.amino) &&
	                    Arrays.equals(gaps, kPart.gaps);
	        }

	        @Override
	        public int hashCode() {
	            int result = Objects.hash(amino);
	            result = 31 * result + Arrays.hashCode(gaps);
	            return result;
	        }
	    }
	 
	 public static class Edge{
		 int startVertex;
		 int endVertex;
		 int weight;
		 
		 public Edge() {
			 
		 }
		 public Edge(int start, int end) {
			 startVertex = start;
			 endVertex = end;
			 weight = 1;
		 }
		 public void addWeight(int start, int end) {
			 startVertex = start;
			 endVertex = end;
			 weight++;
		 }
		 public void addW() {
			 weight++;
		 }
		 public void addEdge(int start, int end) {
			 startVertex = start;
			 endVertex = end;
			 weight = 1;
		 }
		 public String toString() {
			 String s = "from " + startVertex + " to " + endVertex + " with weight = " + weight;
			 return s;
		 }
		 public boolean equals(Object o) {
			 return this.startVertex==((Edge) o).startVertex & this.endVertex==((Edge) o).endVertex;
		 }
		 
	 }
	 public static class DisjointUnionSets {
		    int[] rank, parent;
		    int n;
		  
		    // Constructor
		    public DisjointUnionSets(int n)
		    {
		        rank = new int[n];
		        parent = new int[n];
		        this.n = n;
		        makeSet();
		    }
		  
		    // Creates n sets with single item in each
		    void makeSet()
		    {
		        for (int i = 0; i < n; i++) {
		            // Initially, all elements are in
		            // their own set.
		            parent[i] = i;
		        }
		    }
		  
		    // Returns representative of x's set
		    int find(int x)
		    {
		        // Finds the representative of the set
		        // that x is an element of
		        if (parent[x] != x) {
		            // if x is not the parent of itself
		            // Then x is not the representative of
		            // his set,
		            parent[x] = find(parent[x]);
		  
		            // so we recursively call Find on its parent
		            // and move i's node directly under the
		            // representative of this set
		        }
		  
		        return parent[x];
		    }
		  
		    // Unites the set that includes x and the set
		    // that includes x
		    void union(int x, int y)
		    {
		        // Find representatives of two sets
		        int xRoot = find(x), yRoot = find(y);
		  
		        // Elements are in the same set, no need
		        // to unite anything.
		        if (xRoot == yRoot)
		            return;
		  
		        // If x's rank is less than y's rank
		        if (rank[xRoot] < rank[yRoot])
		  
		            // Then move x under y  so that depth
		            // of tree remains less
		            parent[xRoot] = yRoot;
		  
		        // Else if y's rank is less than x's rank
		        else if (rank[yRoot] < rank[xRoot])
		  
		            // Then move y under x so that depth of
		            // tree remains less
		            parent[yRoot] = xRoot;
		  
		        else // if ranks are the same
		        {
		            // Then move y under x (doesn't matter
		            // which one goes where)
		            parent[yRoot] = xRoot;
		  
		            // And increment the result tree's
		            // rank by 1
		            rank[xRoot] = rank[xRoot] + 1;
		        }
		    }
	 }
	 public static class Motif {
		 List<String> Peptides;
		 double [][] p;
		 double [][] IC;
		 double ic;
		 String kSub;
		 Motif(){
		 }
		 
	 }
	 public static void main(String[] args) throws FileNotFoundException {
		 
		 File file = new File(args[0]);
		System.out.println("Now Processing ... " + args[0]);
	     Scanner sc = new Scanner(file);
	     
		 List<String> Peptides = new ArrayList<String>();

	     while(sc.hasNext()){
			String input = sc.next();
			Peptides.add(input);
		}
	     
		HashMap <KPart, Integer> vertexMap = new HashMap<>(); //Map of vertices (KParts) and their labels (Integer)
	    LinkedHashMap <KPart, Integer> kpartMap = new LinkedHashMap<>() ;//ACCACCA // Map of KParts and their position in peptide
	    
	    Set<KPart> kpartSet = new HashSet();
	    
	    for(int i = 0; i < Peptides.size(); i++) {
	    	kpartSet.addAll(getKMerList(Peptides.get(i)));
	    }
	     
	    List<Integer> vertices = new ArrayList<>();
	    List<KPart> kpartList = new ArrayList<>(kpartSet);
	    
	    for (int i = 0; i < kpartList.size(); i++) {
	    	vertices.add(i);
	    }

	    for(int i = 0; i < kpartList.size(); i++) {
	    	vertexMap.putIfAbsent(kpartList.get(i), i);
	    }

	    
	    System.out.println("Size of vertexMap = " + vertexMap.size());
	    
	    Map<List<Integer>, Edge> edgeList = new HashMap<>();

	    System.out.println("SIZE of kpartList is = " + kpartSet.size());
	    
	    HashMap<KPart, List<String>> kMerMap = new HashMap();//Map of KMers and the peptides containing them
	    int c = 0;
	    Map<Integer, List<Integer>> incomming = new HashMap();
	    Map<Integer, List<Integer>> outgoing = new HashMap();
	    
	    for (String peptide:Peptides) {
	    	kpartMap = getKMer(peptide);
	    //	System.out.println("Done with kpartMap");
	    	c++;
	    	if (c % 10 == 0) {
	              System.out.print("\r" + c);
	         }
	    	// ‌Building Graph
	    	for (Entry<KPart, Integer> entry1: kpartMap.entrySet()) {
	    		if (kMerMap.containsKey(entry1.getKey())) {
	    			kMerMap.get(entry1.getKey()).add(peptide);
	    		}else {
	    			List<String> pepList = new ArrayList<>();
	    			pepList.add(peptide);
	    			kMerMap.put(entry1.getKey(), pepList);
	    		}
	    		
		    	int score1 = kpartMap.get(entry1.getKey());
		    	
	    		for (Entry<KPart, Integer> entry2: kpartMap.entrySet()) {
	    			int score2 = kpartMap.get(entry2.getKey());
	    			if (score1 == score2) continue;
	    			if(score1 < score2) {

	    				int i = vertexMap.get(entry1.getKey());
    					int j = vertexMap.get(entry2.getKey());
	    				Edge eee = new Edge(i, j);
    					if(outgoing.containsKey(i)) {
    						outgoing.get(i).add(j);
    					}else {
    						List<Integer> b = new ArrayList();
        					b.add(j);
        					outgoing.put(i, b);
    					}
    					if(incomming.containsKey(j)) {
    						incomming.get(j).add(i);
    					}else {
    						List<Integer> b = new ArrayList();
        					b.add(i);
        					incomming.put(j, b);
    					}
	    				
	    				if(!hasEdge(edgeList, vertexMap.get(entry1.getKey()), vertexMap.get(entry2.getKey()))){
	    					
	    					List <Integer> a = new ArrayList();
	    					a.add(i);a.add(j);
	    					edgeList.put(a, eee);
	    					
	    				}else {
	    				
	    					List <Integer> a = new ArrayList();
	    					a.add(i);a.add(j);
	    					edgeList.get(a).addW();
	    				}
				    }
	    		}
	    	}	
	    }
	    
	    // Deleting edges with weight less than 5
	    
	    List<Edge> list = new ArrayList();
	    Set<Integer> vertexSet = new HashSet();
	    System.out.println();
	   // System.out.println("Size of list before Deletion = " + list.size());
	    for(Edge e: edgeList.values()) {
	    	if(e.weight >= edgeSupport) {
	    		list.add(e);
	    		vertexSet.add(e.startVertex);
	    		vertexSet.add(e.endVertex);
	    	}
	    }
	    
	    

	    System.out.println("There are " + vertices.size() + " vertices in the graph." );	    
	    System.out.println("Size of EdgeList is = "+edgeList.size());
	    System.out.println("Size of EdgeList(after deletion of edges with lower weights) is = " + list.size());
	    Set<List<Integer>> path = findPath(vertexSet, incomming, outgoing, edgeList);
	    System.out.println("size of path list = "+ path.size());
//	    System.out.println("VertexMap.get(u) = " + vertexMap.values());
	    List<List<String>> peptidesforLogo = new ArrayList<>();
	    
	 // Writing to a file

	      File file1 = new File(args[1]);
	      PrintStream out = new PrintStream(file1);
	      out.println();
	      out.println();
	      out.println("There are " + vertices.size() + " vertices in the graph." );	 
	      out.println("There are " + vertexSet.size() + " notIsolated in the graph." );	
	      int isolatedVertices = vertices.size()-vertexSet.size();
	      out.println("There are " + isolatedVertices + " isolated Vertices in the graph." );	
		  out.println("Size of EdgeList(after deletion of edges with lower weights) is = "+list.size());
		  out.println("size of path list = "+ path.size());
		  out.println();
		  out.println("================================================================================================================");
		  List<Double> icList = new ArrayList();
		  
		  // Map of ic and related Motifs (Set of Peptides of path)
		  HashMap<Double, Set<String>> ICPepMap = new HashMap();
		  HashMap<Double, List<String>> ICPepListMap = new HashMap();
		  	    
		  HashMap<Character, Double> backgroundFrequency = backgroundFreq(Peptides);
//		  HashMap<Character, Double> backgroundFrequency = new HashMap();
//			
//			backgroundFrequency.put('A', 0.05);backgroundFrequency.put('R', 0.05);backgroundFrequency.put('N', 0.05);
//			backgroundFrequency.put('D', 0.05);backgroundFrequency.put('C', 0.05);backgroundFrequency.put('Q', 0.05);
//			backgroundFrequency.put('E', 0.05);backgroundFrequency.put('G', 0.05);backgroundFrequency.put('H', 0.05);
//			backgroundFrequency.put('I', 0.05);backgroundFrequency.put('L', 0.05);backgroundFrequency.put('K', 0.05);
//			backgroundFrequency.put('M', 0.05);backgroundFrequency.put('F', 0.05);backgroundFrequency.put('P', 0.05);
//			backgroundFrequency.put('S', 0.05);backgroundFrequency.put('T', 0.05);backgroundFrequency.put('W', 0.05);
//			backgroundFrequency.put('Y', 0.05);backgroundFrequency.put('V', 0.05);
	  	
		  // Map of ic and nodes of the path (uxv)
		  Map<Double, List<KPart>> icKPartNodesMap = new HashMap();
		  
		  //
		  List<List<String>> motifList = new ArrayList<>();
		  
		  // List of raw (unaligned) motifs of the whole graph (all paths)
		  List<List<String>> rawmotifList = new ArrayList<>();
		  Map <List<String> , String > motifKSubMap = new HashMap();
		  Map <Integer, List<String>> idPepListMap = new HashMap();
		  Map <Integer, String> idKSubMap = new HashMap();
		  int id = 0;
		  List<Motif> allMotifsList = new ArrayList(); 
		  
		  for (List<Integer> l:path) {
			  int u = l.get(0);
			  //System.out.println("u is  equal to = " + u);
			  int x = l.get(1);
			  //System.out.println("x is  equal to = " + x);
			  int v = l.get(2);
			  //System.out.println("v is  equal to = " + v);
			  List<String> alignedPeptides = new ArrayList();
			  List<String> l1Pep = new ArrayList();
			  List<String> l2Pep = new ArrayList();
			  List<String> l3Pep = new ArrayList();
	    	
			  l1Pep.addAll(kMerMap.get(kpartList.get(u)));
			  l2Pep.addAll(kMerMap.get(kpartList.get(x)));
			  l3Pep.addAll(kMerMap.get(kpartList.get(v)));
            
			  l1Pep.retainAll(l2Pep); // Peptides that contain nodes u and x
			  l2Pep.retainAll(l3Pep); // Peptides that contain nodes x and v
	    	
			  l1Pep.addAll(l2Pep); // Peptides of the whole path
			  Set<String> s1 = new HashSet(l1Pep);
			  List<String> ls1 = new ArrayList(s1);
			 // List<String> aligned = alignPeptides(l1, kpartList.get(u), kpartList.get(x), kpartList.get(v));
			  
			  // Aligning Peptides of the path
			  List<List<String>> rawList = new ArrayList();
			  for(int i = 0 ; i < ls1.size()-1; i++) {
				  rawList.add(makeRawAlignment(ls1.get(i), ls1.get(i+1)));
			  }
			  List<List<String>> finalList = finalAlignment(rawList);
			  List<String> aligned = alignList(finalList);
				
			  List<KPart> kpartNodes = Arrays.asList(kpartList.get(u), kpartList.get(x), kpartList.get(v));
			  
			  List <String> newAligned = seqAlignment(ls1, kpartNodes.get(0), kpartNodes.get(1), kpartNodes.get(2));
			  
			  List<String> bestAlignment = SeqAlgnmnt(ls1, kpartNodes.get(0), kpartNodes.get(1), kpartNodes.get(2));
			  List<String> bestAlignmentCleaned = motifCleaning(bestAlignment, backgroundFrequency); 
			  int kp1Count = 0, kp3Count = 0;
			 // System.out.println("  before 2 hmmm ... " + bestAlignmentCleaned);
			  for (String s:bestAlignmentCleaned) {	  
				  //System.out.println(" before hmmm ... " + s);
				  if (getKMer(removeDash(s)).containsKey(kpartNodes.get(0))) {
					  kp1Count++;
				  }
				  if (getKMer(removeDash(s)).containsKey(kpartNodes.get(2))) {
					  kp3Count++;
				  }
			  }
			  if (kp1Count >= edgeSupport && kp3Count >= edgeSupport) { 
				  int ed1,ed2 = 0;
				  for(String s:bestAlignmentCleaned) {
					  
				  }
			  
				  //System.out.println(".........H E R E .......");
//				  System.out.println(kpartNodes.get(0));
//				  System.out.println(kpartNodes.get(1));
//				  System.out.println(kpartNodes.get(2));
//				  System.out.println("    newAligned     START   ");
//				  for(String s: newAligned) {
//					  System.out.println(s);
//				  }	
//				  System.out.println("    newAligned     END   ");
//				  System.out.println();
				  Set<String> alignedSet = new HashSet(bestAlignmentCleaned);
				  List<String> ls = new ArrayList(alignedSet);
			  
			//	  Set<String> alignedSet = new HashSet(aligned);
//				  System.out.println("    aligned     START   ");
//				  for(int i = 0; i < aligned.size(); i++) {
//					  System.out.println(aligned.get(i));
//				  }
//				  System.out.println("    aligned     END   ");
//				  System.out.println();
//				  System.out.println("    bestAlignment     START   ");
//				  System.out.println(" Motif = " + getKSubSet(kpartNodes.get(0), kpartNodes.get(1), kpartNodes.get(2)));
//				  for(int i = 0; i < aligned.size(); i++) {
//					  System.out.println(bestAlignment.get(i));
//				  }
//				  System.out.println("    bestAlignment     END   ");
		//		  motifParamComputations(aligned, backgroundFrequency);
				  double[][] pwm = makeProfile(bestAlignment);
		//		  printMatrix(pwm);
				  //System.out.println();
				  int n = bestAlignment.get(0).length();
				  double[][]IC = informationContent(pwm, n, backgroundFrequency);
		//		  printMatrix(IC);  		
	    		
				  double ic = icOfMotif(IC);
				  icList.add(ic);
				 // if (alignedSet.size() >= pathSupport) {
					  ICPepMap.put(ic, alignedSet);
					  ICPepListMap.put(ic, ls);
					  icKPartNodesMap.putIfAbsent(ic, kpartNodes);
				  //}
				  
//				  System.out.println("information content = " + ic);
//				  System.out.println();
		//		  motifList.add(aligned);
				  rawmotifList.add(ls1);			  
				  motifKSubMap.put(ls1, getKSubSet(kpartList.get(u), kpartList.get(x), kpartList.get(v)));
				  idPepListMap.put(id, ls1);
				  idKSubMap.put(id, getKSubSet(kpartList.get(u), kpartList.get(x), kpartList.get(v)));
				  id++;
				  Motif m = new Motif();
				  m.Peptides = ls; m.ic = ic; m.kSub = getKSubSet(kpartList.get(u), kpartList.get(x), kpartList.get(v));
				  m.p = pwm; m.IC = IC;
				  allMotifsList.add(m);
			  }
		  }
//		  System.out.println();
//		  System.out.println("****** F I N I S H ******");
//		  System.out.println("****** F I N I S H ******");
//		  System.out.println("****** F I N I S H ******");
//		  System.out.println();
	    // Sorting ICPepMap based on ic
		Map<Double,Set<String>> sortedICPepMap = new TreeMap<Double, Set<String>>(Collections.reverseOrder());
		sortedICPepMap.putAll(ICPepMap);
		Map<Double,List<String>> sortedICPepListMap = new TreeMap<Double, List<String>>(Collections.reverseOrder());
		sortedICPepListMap.putAll(ICPepListMap);
		Map<Integer, Set<String>> sortedICPepMapRank = new HashMap<>();
		int rank = 1;
		for(Double d:sortedICPepMap.keySet()) {
			sortedICPepMapRank.put(rank, sortedICPepMap.get(d));
			rank++;
		}
		
		Map<Integer, List<String>> sortedICPepListMapRank = new HashMap<>();
		int rank2 = 1;
		for(Double d:sortedICPepListMap.keySet()) {
			sortedICPepListMapRank.put(rank2, sortedICPepListMap.get(d));
			rank2++;
		}

			int index = 1;
			for(Entry<Double, Set<String>> entry: sortedICPepMap.entrySet()) {
				List<String> temp = new ArrayList(entry.getValue());
				out.println();
				out.println(index + ". ");
				index ++;
				out.println();
				out.println("KParts(NODES) for this motif: "+icKPartNodesMap.get(entry.getKey()));
				out.println();
				out.println("   KSubSet:   ");
				out.println("   "+getKSubSet(icKPartNodesMap.get(entry.getKey()).get(0), icKPartNodesMap.get(entry.getKey()).get(1),icKPartNodesMap.get(entry.getKey()).get(2)));
				out.println();
				System.out.println(icKPartNodesMap.get(entry.getKey()));
				if (temp.size() == 0) continue;
				for(int i = 0; i < temp.size(); i++) {
					out.println(temp.get(i));
		    		System.out.println(temp.get(i));
		    	}
//				 System.out.println(" getKey() = "+entry.getKey());
//				 System.out.println(" getValue() = "+entry.getValue());
//				 System.out.println();
//				 System.out.println("  ICPepMap  ");
//				 System.out.println(ICPepMap.get(entry.getKey()));
//				 System.out.println("     S T O P ");
		         out.println();
		         double[][] pwm = makeProfile(temp);
		//            printMatrix(pwm);
		            System.out.println();
		            out.println();
		            out.println();
		            out.println("Position Probability Matrix: ");
		            out.println();
		            /*
		             * Start Printing PWM matrix in the output file
		             */
		            int i1 = 0;
		    		out.println("Pos   " +1+"     "+2+"     "+3+"     "+4+"     "+5+"     "+6+"     "+7+"     "+8+"     "+9+"     "+10+"    "+11);
		    		out.println("----------------------------------------------------------------------");
		    		
		    		for (double[] row : pwm) {
		    			out.print(AminoAcid[i1] + ": ");
		    			i1++;
		    			for (double element : row) {
		    				out.printf("%6.2f", element);
		    			}
		    			out.println();
		    		}
		    		//Finish Printing PWM matrix to the output file
		    		
		            int n = temp.get(0).length();
		            double[][]IC = informationContent(pwm, n, backgroundFrequency);
		//    		printMatrix(IC);
		    		out.println();
		    		out.println("Information Content Matrix: ");
		    		out.println();
		    		
		    		/*
		             * Start Printing IC matrix in the output file
		             */
		            int i2 = 0;
		    		out.println("Pos   " +1+"     "+2+"     "+3+"     "+4+"     "+5+"     "+6+"     "+7+"     "+8+"     "+9+"     "+10+"    "+11);
		    		out.println("----------------------------------------------------------------------");
		    		
		    		for (double[] row : IC) {
		    			out.print(AminoAcid[i2] + ": ");
		    			i2++;
		    			for (double element : row) {
		    				out.printf("%6.2f", element);
		    			}
		    			out.println();
		    		}
		    		//Finish Printing IC matrix to the output file
		    		
		    		
//		            double[][] IC = new double[20][aligned.get(0).length()];
		            out.println();
		            double ic = icOfMotif(IC);
//		            System.out.println("information content = " + ic);
//		            System.out.println();
		            out.println("information content = " + ic);
		            out.println();
		            out.println("================================================================================================================");
		            out.println("================================================================================================================");
		            
		    }

		
		// Set<String> pureMotifs = new HashSet();
		// Map<String, Set<String>> motifIdToMotif = new HashMap<>();

		//  pureMotifs = combineMotifs(rawmotifList, backgroundFrequency, motifIdToMotif);
		  List<List<String>> pureMotifsList = new ArrayList<>();
		 // pureMotifsList = purify(rawmotifList, backgroundFrequency);
		  
		  Map <Integer, String> newIdKSubMap1 = new HashMap();
			int id_new1 = 0;
			for (Entry<Double, List<String>> entry:sortedICPepListMap.entrySet()) {
				for (Motif m:allMotifsList) {
					if (m.Peptides == entry.getValue()) {
						newIdKSubMap1.putIfAbsent(id_new1, m.kSub);
						break;
					}
				}
				id_new1++;
			}
		  
		  pureMotifsList = purify6(sortedICPepListMapRank, sortedICPepMap, newIdKSubMap1, backgroundFrequency);
		 // Set<List<String>> pureMotifSet = new HashSet(pureMotifsList);
		  
		  
		// Writing purified motifs to a new file
		   File file2 = new File(args[2]);
	       PrintStream outt = new PrintStream(file2);
	       outt.println();
	       outt.println();
	       outt.println("edgeSupport = " + edgeSupport);
	       outt.println("pathSupport = " + pathSupport);
	       outt.println();
	       outt.println();
	    //   outt.println("There are " + pureMotifSet.size() + " pure motifs");
	       outt.println();    
	       outt.println();
	       int indexOutt = 1;
	       for(List<String> s:pureMotifsList) {
//	    	   	List<List<String>> rawList33 = new ArrayList();
//	   			for(int i = 0 ; i < s.size()-1; i++) {
//	   				rawList33.add(makeRawAlignment(s.get(i), s.get(i+1)));
//	   			}
//	   			List<List<String>> finalList33 = finalAlignment(rawList33);
//	   			List<String> aligned33 = alignList(finalList33);
//	   			Set<String> alignedSet33 = new HashSet(aligned33);
//	   			List<String> alignedList33 = new ArrayList(alignedSet33);	    	 
//	   			
//	    	    alignedList33.forEach(outt::println);
	    	   	outt.println();
	    	   	double[][] pwm = makeProfile(s);
	    	    int n = s.get(0).length();
	    	    double[][]IC = informationContent(pwm, n, backgroundFrequency);	  
	    	    double ic = icOfMotif(IC);
	    	    if (ic > 5) {
	    	    	List<List<String>> rawList = new ArrayList();
	  			  for(int i = 0 ; i < s.size()-1; i++) {
	  				  rawList.add(makeRawAlignment(s.get(i), s.get(i+1)));
	  			  }
	  			  List<List<String>> finalList = finalAlignment(rawList);
		    	   	outt.println(indexOutt+ ". ");
		    	   	outt.println();
		    	   	indexOutt++;
		    	   	for (List<String>lstr:finalList) {
		    	   		Set<String> set = new HashSet(lstr);
		    	   		set.forEach(outt::println);
		    	   	}
		    	    outt.println();
	//	    	    System.out.println("S");
	//	    	    for(String str:s) {
	//	    	    	System.out.println(str);
	//	    	    }
		    	    
		    		outt.println();
		    		outt.println();
		    		outt.println("Position Probability Matrix: ");
		    		outt.println();
		    		/*
		    		* Start Printing PWM matrix in the output file
		    		*/
		    		int i1 = 0;
		    		outt.println("Pos   " +1+"     "+2+"     "+3+"     "+4+"     "+5+"     "+6+"     "+7+"     "+8+"     "+9+"     "+10+"    "+11);
		    		outt.println("----------------------------------------------------------------------");
		    		    		
		    		for (double[] row : pwm) {
		    		    outt.print(AminoAcid[i1] + ": ");
		    		    i1++;
		    		    for (double element : row) {
		    		    	outt.printf("%6.2f", element);
		    		    }
		    		    outt.println();
		    		 }
		    		 //Finish Printing PWM matrix to the output file
		    		 
		    		 outt.println();
		    		 outt.println("Information Content Matrix: ");
		    		 outt.println();
		    		    		
		    		 /*
		    		  * Start Printing IC matrix in the output file
		    		  */
		    		 int i2 = 0;
		    		 outt.println("Pos   " +1+"     "+2+"     "+3+"     "+4+"     "+5+"     "+6+"     "+7+"     "+8+"     "+9+"     "+10+"    "+11);
		    		 outt.println("----------------------------------------------------------------------");
		    		    		
		    		 for (double[] row : IC) {
		    			 outt.print(AminoAcid[i2] + ": ");
		    			 i2++;
		    			 for (double element : row) {
		    				 outt.printf("%6.2f", element);
		    			 }
		    			 outt.println();
		    		 }
		    		  //Finish Printing IC matrix to the output file
		    		    		
		    		    		
	//	    		   double[][] IC = new double[20][aligned.get(0).length()];
		    		   outt.println();
		    		   outt.println("information content = " + ic);
		    		   outt.println();
	    	    }
	   			outt.println("================================================================================================================");
	   			outt.println("================================================================================================================");
	       }
	       
	    	   outt.println("================================================================================================================");
               outt.println("================================================================================================================");
	    /**
	     * TEST of LinkDirectedGraph on hand-generated data
	     
               
               
	       System.out.println("------#############-------------------------------------------");
	       System.out.println("------#############-------------------------------------------");
	       System.out.println("------#############-------------------------------------------");
	       System.out.println("------#############-------------------------------------------");
	       System.out.println("------#############-------------------------------------------");
	       System.out.println("------#############-------------------------------------------");
	       System.out.println("------#############---------HEEEERRREEEEEEEE----------------------------------");
	       System.out.println("------#############-------------------------------------------");
	       List<String> list1 = Arrays.asList("IGLAIAD", "IGLAIPV", "IGLAIPY", "IGLAILT", "SIGLAIP", 
	    		   "LIGLAIP", "PIGLAIR", "IGLAIPN", "PIGLAIT", "PIGLAIS", "AIGLAIP", "PIGLAIP", "PIGLAIL" );
	       
	       List<String> list2 = Arrays.asList("IGLAIPV", "LGLAIPD", "IGLAIPY", "SIGLAIP", "PIGLAIR", "IGLAIPN", 
	    		   "PIGLAIT", "IGLAIPE", "PIGLAIP", "SPIGLAI", "TPIGLAI", "YPIGLAI", "PIGLAIA", "FGLAIPD", "NGLAIPD" );
	       
	       List<String> list3 = Arrays.asList("IGLAIPH", "IGLAIPI", "IGLAIPN", "PIGLAIS", "AIGLAIP", "IGLAIPV", 
	    		   "IGLAILT", "IGLAIPY", "ACIGLAI", "PIGLAIT", "YPIGLAI", "GAIGLAI", "PIGLASL", "PIGLAIL" );
	       
	       List<String> list4 = Arrays.asList("AETVQSC", "AETVASC", "AETVESC", "AETVDSC" );
	       
	       List<String> list5 = Arrays.asList("PIGLAIP", "SPIGLAI", "TPIGLAI", "PIGLAIL", "YPIGLAI", "PIGLAIH", 
	    		   "PIGLAIR", "PIGLAIA", "PIGLAIT", "PIGLAIS" );
	       
	       List<String> list6 = Arrays.asList("PIGLANP", "PIELAIP", "PIALAIP", "PIGLASP", "PIGLAIP", "PIRLAIP", 
	    		   "TLAHPTV", "SLAHPTV", "PIGLAFP", "PIGLAMP", "LAHPTML", "LAQPTKL", "LAQPTML", "PIVLAIP", "XXXXXXX", "YYYYYYY", "OOOOOOO" ); //LALPTSY
	       
	       List<String> list7 = Arrays.asList("VSPTSRA", "SSTSMAR", "ATSRAAW", "QTSVAAS", "MTSWAAP", "QSNTSVA", "MTSVAAP", 
	    		   "NSTTSTA", "QSLTSLA", "MTSGAAP", "HSNTSVA", "ITSGAAP", "ETSLAAP", "FSITSCA" );
	       
	       List<String> list8 = Arrays.asList("AETVEGC", "AETVEIC", "AETVECC", "AETVESC", "AETVETC", "AETVERC" );
	       
	       List<String> list9 = Arrays.asList("IGLAIPV", "LGLAIPD", "IGLAIPY", "IGLAILT", 
	    		   "SIGLAIP", "PIGLAIR", "IGLAIPN", "PIGLAIT", "IGLAIPE", "PIGLAIP", "IGLAIPD", "PNGLAIP",
	    		   "IGLAIPI", "PIGLAIL", "IGLAIPH", "NGLAIPD", "PIGLAIH", "IGLAITD", "IGLAIAD", "FGLAIPD",
	    		   "PLGLAIP", "LIGLAIP", "PVGLAIP", "PIGLAIS", "AIGLAIP", "PIGLAIA", "IGLAISD");

	       List<String> list10 = Arrays.asList(
	    		   "IGLAIPV", "LGLAIPD", "IGLAIPY", "SIGLAIP", "PIGLAIR", "IGLAIPN", "PIGLAIT", "IGLAIPE",
	    		   "PIGLAIP", "SPIGLAI", "IGLAIPD", "PNGLAIP", "TPIGLAI", "IGLAIPI", "PIGLAIL", "IGLAIPH",
	    		   "NGLAIPD", "PIGLAIH", "FGLAIPD", "PLGLAIP", "LIGLAIP", "PVGLAIP", "PIGLAIS", "AIGLAIP",
	    		   "YPIGLAI", "PIGLAIA");	 
	       
	       
	       KPart kp1 = new KPart("IG-A", new int[] {0, 0, 0});
	       KPart kp2 = new KPart("GL-I", new int[] {0, 0, 0});
	       KPart kp3 = new KPart("AI-D", new int[] {0, 0, 0});
	       String kSub1 = getKSubSet(kp1, kp2, kp3);
	       
	       KPart kp4 = new KPart("PI-L", new int[] {0, 0, 0});
	       KPart kp5 = new KPart("GLAI", new int[] {0, 0, 0});
	       KPart kp6 = new KPart("L-IP", new int[] {0, 0, 0});
	       String kSub2 = getKSubSet(kp4, kp5, kp6);
	       
	       KPart kp7 = new KPart("PIGL", new int[] {0, 0, 0});
	       KPart kp8 = new KPart("IGLA", new int[] {0, 0, 0});
	       KPart kp9 = new KPart("G-AI", new int[] {0, 0, 0});
	       String kSub3 = getKSubSet(kp7, kp8, kp9);
	       
	       KPart kp10 = new KPart("A-TV", new int[] {0, 0, 0});
	       KPart kp11 = new KPart("TV-S", new int[] {0, 0, 0});
	       KPart kp12 = new KPart("V-SC", new int[] {0, 0, 0});
	       String kSub4 = getKSubSet(kp10, kp11, kp12);
	       
	       KPart kp13 = new KPart("PI-L", new int[] {0, 0, 0});
	       KPart kp14 = new KPart("GLAI", new int[] {0, 0, 0});
	       KPart kp15 = new KPart("AI-D", new int[] {0, 0, 0});
	       String kSub5 = getKSubSet(kp13, kp14, kp15);
	       
	       KPart kp16 = new KPart("PI-L", new int[] {0, 0, 0});
	       KPart kp17 = new KPart("LA-P", new int[] {0, 0, 0});
	       KPart kp18 = new KPart("A-PT", new int[] {0, 0, 0});
	       String kSub6 = getKSubSet(kp16, kp17, kp18);
	       
	       KPart kp19 = new KPart("S-TS", new int[] {0, 0, 0});
	       KPart kp20 = new KPart("TS-A", new int[] {0, 0, 0});
	       KPart kp21 = new KPart("S-AA", new int[] {0, 0, 0});
	       String kSub7 = getKSubSet(kp19, kp20, kp21);
	       
	       KPart kp22 = new KPart("AE-V", new int[] {0, 0, 0});
	       KPart kp23 = new KPart("E-VE", new int[] {0, 0, 0});
	       KPart kp24 = new KPart("VE-C", new int[] {0, 0, 0});
	       String kSub8 = getKSubSet(kp22, kp23, kp24);
	       
	       KPart kp25 = new KPart("I-LA", new int[] {0, 0, 0});
	       KPart kp26 = new KPart("G-AI", new int[] {0, 0, 0});
	       KPart kp27 = new KPart("LA-P", new int[] {0, 0, 0});
	       String kSub9 = getKSubSet(kp25, kp26, kp27);
	       
	       KPart kp28 = new KPart("PI-L", new int[] {0, 0, 0});
	       KPart kp29 = new KPart("GLAI", new int[] {0, 0, 0});
	       KPart kp30 = new KPart("L-IP", new int[] {0, 0, 0});
	       String kSub10 = getKSubSet(kp28, kp29, kp30);
	       
	       List<List<String>> mainList = new ArrayList();
	       mainList.add(list1); mainList.add(list2); mainList.add(list3); mainList.add(list4);
	       mainList.add(list5); mainList.add(list6); mainList.add(list7); mainList.add(list8);
	       mainList.add(list9); mainList.add(list10);
	       
	       Map<List<String>, String> motksubmap = new HashMap();
	       motksubmap.put(list1, kSub1); motksubmap.put(list2, kSub2); motksubmap.put(list3, kSub3);
	       motksubmap.put(list4, kSub4); motksubmap.put(list5, kSub5); motksubmap.put(list6, kSub6); motksubmap.put(list7, kSub7);
	       motksubmap.put(list8, kSub8); motksubmap.put(list9, kSub9); motksubmap.put(list9, kSub9);
	       
	       Map <Integer, String> idksubmap = new HashMap();
	       idksubmap.put(0, kSub1); idksubmap.put(1, kSub2); idksubmap.put(2, kSub3); 
	       idksubmap.put(3, kSub4); idksubmap.put(4, kSub5); idksubmap.put(5, kSub6); idksubmap.put(6, kSub7);
	       idksubmap.put(7, kSub8); idksubmap.put(8, kSub9); idksubmap.put(9, kSub10);
	       
	       Map<Integer, List<String>> idPeptides = new HashMap();
	       idPeptides.put(1, list1); idPeptides.put(2, list2); idPeptides.put(3, list3);
	       idPeptides.put(4, list4); idPeptides.put(5, list5); idPeptides.put(6, list6); idPeptides.put(7, list7);
	       idPeptides.put(8, list8); idPeptides.put(9, list9); idPeptides.put(10, list10);
	       
	       System.out.println("M0 = " + kSub1);
	       System.out.println("M2 = " + kSub3);
	       System.out.println("M5 = " + kSub6);
	       System.out.println("Motif Alignment betw 0 and 5 = " + kSubAlignmentCalc(kSub1,kSub6));
	       System.out.println("Motif Alignment betw 2 and 5 = " + kSubAlignmentCalc(kSub3,kSub6));
	       
	       System.out.println("   Size of Path = " + path.size());
	       System.out.println("   Size of rawmotifList = " + rawmotifList.size());
	       
	       HashMap<Double, Set<String>> icPepMap = new HashMap();
	       //HashMap<Double, Set<String>> sortedicPepMap = new HashMap();
	       //HashMap<Double, Set<String>> sortedicPepMapRank = new HashMap();
	       
	       List <String> newAligned0 = seqAlignment(list1, kp1, kp2, kp3); Set<String> alignedSet0 = new HashSet(newAligned0);
	       List <String> newAligned1 = seqAlignment(list2, kp4, kp5, kp6); Set<String> alignedSet1 = new HashSet(newAligned1);
	       List <String> newAligned2 = seqAlignment(list3, kp7, kp8, kp9); Set<String> alignedSet2 = new HashSet(newAligned2);
	       List <String> newAligned3 = seqAlignment(list4, kp10, kp11, kp12); Set<String> alignedSet3 = new HashSet(newAligned3);
	       List <String> newAligned4 = seqAlignment(list5, kp13, kp14, kp15); Set<String> alignedSet4 = new HashSet(newAligned4);
	       List <String> newAligned5 = seqAlignment(list6, kp18, kp17, kp18); Set<String> alignedSet5 = new HashSet(newAligned5);
	       List <String> newAligned6 = seqAlignment(list7, kp19, kp20, kp21); Set<String> alignedSet6 = new HashSet(newAligned6);
	       List <String> newAligned7 = seqAlignment(list8, kp22, kp23, kp24); Set<String> alignedSet7 = new HashSet(newAligned7);
	       List <String> newAligned8 = seqAlignment(list9, kp25, kp26, kp27); Set<String> alignedSet8 = new HashSet(newAligned8);
	       List <String> newAligned9 = seqAlignment(list10, kp28, kp29, kp30); Set<String> alignedSet9 = new HashSet(newAligned9);
	       
	       double[][] pwm0 = makeProfile(newAligned0);
	       double[][] pwm1 = makeProfile(newAligned1);
	       double[][] pwm2 = makeProfile(newAligned2);
	       double[][] pwm3 = makeProfile(newAligned3);
	       double[][] pwm4 = makeProfile(newAligned4);
	       double[][] pwm5 = makeProfile(newAligned5);
	       double[][] pwm6 = makeProfile(newAligned6);
	       double[][] pwm7 = makeProfile(newAligned7);
	       double[][] pwm8 = makeProfile(newAligned8);
	       double[][] pwm9 = makeProfile(newAligned9);
	       
	       int n0 = newAligned0.get(0).length();
	       double[][]IC0 = informationContent(pwm0, n0, backgroundFrequency);
	       int n1 = newAligned1.get(0).length();
	       double[][]IC1 = informationContent(pwm1, n1, backgroundFrequency);
	       int n2 = newAligned2.get(0).length();
	       double[][]IC2 = informationContent(pwm2, n2, backgroundFrequency);
	       int n3 = newAligned3.get(0).length();
	       double[][]IC3 = informationContent(pwm3, n3, backgroundFrequency);
	       int n4 = newAligned4.get(0).length();
	       double[][]IC4 = informationContent(pwm4, n4, backgroundFrequency);
	       int n5 = newAligned5.get(0).length();
	       double[][]IC5 = informationContent(pwm5, n5, backgroundFrequency);
	       int n6 = newAligned6.get(0).length();
	       double[][]IC6 = informationContent(pwm6, n6, backgroundFrequency);
	       int n7 = newAligned7.get(0).length();
	       double[][]IC7 = informationContent(pwm7, n7, backgroundFrequency);
	       int n8 = newAligned8.get(0).length();
	       double[][]IC8 = informationContent(pwm8, n8, backgroundFrequency);
	       int n9 = newAligned9.get(0).length();
	       double[][]IC9 = informationContent(pwm9, n9, backgroundFrequency);
	       
	       double i0 = icOfMotif(IC0); double i1 = icOfMotif(IC1); double i2 = icOfMotif(IC2); double i3 = icOfMotif(IC3); double i4 = icOfMotif(IC4);
	       double i5 = icOfMotif(IC5); double i6 = icOfMotif(IC6); double i7 = icOfMotif(IC7); double i8 = icOfMotif(IC8); double i9 = icOfMotif(IC9);

	       icPepMap.put(i0, alignedSet0); icPepMap.put(i1, alignedSet1); icPepMap.put(i2, alignedSet2);
	       icPepMap.put(i3, alignedSet3); icPepMap.put(i4, alignedSet4); icPepMap.put(i5, alignedSet5);
	       icPepMap.put(i6, alignedSet6); icPepMap.put(i7, alignedSet7); icPepMap.put(i8, alignedSet8); icPepMap.put(i9, alignedSet9);
	       
	       Map<Double, List<String>> icPepMapList = new HashMap();
	       icPepMapList.put(i0, newAligned0); icPepMapList.put(i1, newAligned1); icPepMapList.put(i2, newAligned2);
	       icPepMapList.put(i3, newAligned3); icPepMapList.put(i4, newAligned4); icPepMapList.put(i5, newAligned5);
	       icPepMapList.put(i6, newAligned6); icPepMapList.put(i7, newAligned7); icPepMapList.put(i8, newAligned8);
	       icPepMapList.put(i9, newAligned9);
	       
	       Motif m0 = new Motif(); m0.Peptides = newAligned0; m0.ic = i0; m0.kSub = kSub1;
	       Motif m1 = new Motif(); m1.Peptides = newAligned1; m1.ic = i1; m1.kSub = kSub2;
	       Motif m2 = new Motif(); m2.Peptides = newAligned2; m2.ic = i2; m2.kSub = kSub3;
	       Motif m3 = new Motif(); m3.Peptides = newAligned3; m3.ic = i3; m3.kSub = kSub4;
	       Motif m4 = new Motif(); m4.Peptides = newAligned4; m4.ic = i4; m4.kSub = kSub5;
	       Motif m5 = new Motif(); m5.Peptides = newAligned5; m5.ic = i5; m5.kSub = kSub6;
	       Motif m6 = new Motif(); m6.Peptides = newAligned6; m6.ic = i6; m6.kSub = kSub7;
	       Motif m7 = new Motif(); m7.Peptides = newAligned7; m7.ic = i7; m7.kSub = kSub8;
	       Motif m8 = new Motif(); m8.Peptides = newAligned8; m8.ic = i8; m8.kSub = kSub9;
	       Motif m9 = new Motif(); m9.Peptides = newAligned9; m9.ic = i9; m9.kSub = kSub10;
	       
	       List<Motif> motList = new ArrayList();
	       motList.add(m0);motList.add(m1);motList.add(m2);motList.add(m3);motList.add(m4);
	       motList.add(m5);motList.add(m6);motList.add(m7);motList.add(m8);motList.add(m9);
	       
	       
	       
	       Map<Double,Set<String>> sortedicPepMap = new TreeMap<Double, Set<String>>(Collections.reverseOrder());
	       sortedicPepMap.putAll(icPepMap);
			Map<Integer, Set<String>> sortedicPepMapRank = new HashMap<>();
			int r = 1;
			for(Double d:sortedicPepMap.keySet()) {
				sortedicPepMapRank.put(r, sortedicPepMap.get(d));
				r++;
			}
	       
		       Map<Double,List<String>> sortedicPepMapList = new TreeMap<Double, List<String>>(Collections.reverseOrder());
		       sortedicPepMapList.putAll(icPepMapList);
				Map<Integer, List<String>> sortedicPepMapListRank = new HashMap<>();
				int r2 = 1;
				for(Double d:sortedicPepMapList.keySet()) {
					sortedicPepMapListRank.put(r2, sortedicPepMapList.get(d));
					r2++;
				}
				
				//Motif m = new Motif();
				
				Map <Integer, String> newIdKSubMap = new HashMap();
				int id_new = 0;
				for (Entry<Double, List<String>> entry:sortedicPepMapList.entrySet()) {
					for (Motif m:motList) {
						if (m.Peptides == entry.getValue()) {
							newIdKSubMap.putIfAbsent(id_new, m.kSub);
							break;
						}
					}
					id_new++;
				}
				
			System.out.println("newIdKSubMap " + newIdKSubMap);
			List<List<String>> pureMotifsListExample = new ArrayList<>();
			pureMotifsListExample = purify6(sortedicPepMapListRank,  sortedicPepMap, newIdKSubMap, backgroundFrequency);
			
			int ind = 1;
			for(List<String> s:pureMotifsListExample) {
	    	   	System.out.println();
	    	   	System.out.println(ind+ ". ");
	    	   	System.out.println();
	    	   	ind++;
	    	    s.forEach(System.out::println);
	   			System.out.println("================================================================================================================");
	   			System.out.println("================================================================================================================");
	       }
			
	   
			System.out.println("kp16 = " + kp16);
   			System.out.println("kp17 = " + kp17);
   			System.out.println("kp18 = " + kp18);
   			System.out.println("kSub6 = " + kSub6);
   			System.out.println();
   			System.out.println("kp19 = " + kp19);
   			System.out.println("kp20 = " + kp20);
   			System.out.println("kp21 = " + kp21);
   			System.out.println("kSub7 = " + kSub7);
   			System.out.println();
   			System.out.println();
 	       */
              
               //System.out.println("Check here = " + kSubAlignment("AE-VE-C", "AETVE-C"));
               
               
               List<String>l0 = Arrays.asList("AIGLAIP","LGLAIPD","IGLAISD","IGLAITD","IGLAIPD","VGLAIPD","PIGFAIP","IGHAIPD",
            		   "ASIGLAI","IGVAIPD","NGLAIPD","IGLAIPT","HIGLAIP","PIGLAIA","PIGLAIP","AWIGLAI","IGLAIPY","TIGLAIP",
            		   "MGLAIPD","IGLAIPV", "IGLAILT", "PIGLAIL", "PIGHAIP" );
               List<String>l1 = Arrays.asList("PIWLAIP","IGLAIPV","IRLAIPD","AIGLAIP","TIGLAIP","FPITLAH","PIGLAIA","PIRLAIP"
            		   ,"HIGLAIP","IGLANPD","IGLATPD","GGPIGLA","YVPIGLA","ILLAPPR","PIGLANP","IELAIPD","IGLAIPY","PIELAIP",
            		   "IGLAIPD","PIVLAIP","IGLAMPD","PIGLAIL","IGLAIPT","PIGLAIP");
               List<String> l2 = Arrays.asList("LLPPVLQ", "LLPPLLF", "LFPPISI", "LLTPPQI", "LLAPPYW", "FFTLLPP", "TLWPPNS", 
            		   "LLLPPSH", "ILLAPPR", "TKLLCPP", "TLLLPPW", "LLRPPTP", "DLLMPPL", "LNPPPFS", "TLLPPTH", "LLPPISI", "DVVLLPP", 
            		   "HLLMPPL", "FLLPPIS", "LLPPKFT", "ISGLLPP", "FLNPPIS", "VLLPPSP", "VLYPPTS", "LTPPISS", "GTSLLPP", "TLLPPPM",
            		   "VTLLTPP", "VLFPPTS", "LLLPPPK", "NHLLPPL");
               List<String> l3 = Arrays.asList("ALRPPSP","LLPPVLQ","LLPPLLF","YLFPPTP","LLTPPQI","LLAPPYW","FFTLLPP","LLLPPSH",
            		   "ILLAPPR","LVPPHPV","TKLLCPP","TLLLPPW","QLQPPVP","LLRPPTP","DLLMPPL","DLSPPTP",
            		   "TLLPPTH","FLHPPMP","LLPPISI","DVVLLPP","HLLMPPL","VLFPPLP","FLLPPIS","WLRPPFP",
            		   "LLPPKFT","ISGLLPP","VLLPPSP","GTSLLPP","TLLPPPM","VTLLTPP","LLLPPPK","RLPPPPV",
            		   "NHLLPPL", "FLHPPLP", "SLYPPSP");
               List<String> l4 = Arrays.asList("AETFESC","AETVEIS","TETVESC","AETVESC","AETVETC","AETDESS","AETLESC","DETDESC",
            		   "DETVESC","AETVESY","PETVESC","AETVESS","AETDESC","AETVERG","AETVESI","AETVECC",
            		   "AETVERC","SETVESC","ETAESRL","AETVESF","AETVENC","AETVESG");
               List<String> l5 = Arrays.asList("ATTAFTL","AETLESC","WAETLYS","TAITLTT","YGTAVTL","WAVTLYS","TAATLNS","TAITLAT",
            		   "ATAATLP","ANTLPSR","TALTLRL","ARTLYSP","AETLVSY","RAETLYS","TAHTLPQ");
               List<String> l6 = Arrays.asList("TLMSLDN","RAETLYS","AETLESC","TTLWSLH","TLVSLLE","ARTLYSP","ANTLPSR","WAVTLYS",
            		   "AETLVSY","TAATLNS","TLMSLEN","WAETLYS","RTLHSLE","DATLWSL","EATLWSL","TLASLSE",
            		   "TLLSLIR");
               
               Map<Integer, String> idKSubTest = new HashMap();
               idKSubTest.put(0, "IG-AI-D"); idKSubTest.put(1, "PI-LA-P");idKSubTest.put(2, "LL-PP-S"); idKSubTest.put(3, "LL-PP-P");
               idKSubTest.put(4, "AETVES");idKSubTest.put(5, "TA-TL-S");idKSubTest.put(6, "A-TL-SL");
//               System.out.println("Test areSame Function");
//               System.out.println("BTW 8 and 18 = "+areSame(l0, l1, idKSubTest, 0, 1, backgroundFrequency));
//               System.out.println();
//               System.out.println("BTW 18 and 198 = "+areSame(l1, l2, idKSubTest, 1, 2, backgroundFrequency));
//               System.out.println();
//               System.out.println("BTW 18 and 206 = "+areSame(l1, l3, idKSubTest, 1, 3, backgroundFrequency));
//               System.out.println();
//               System.out.println("BTW 198 and 206 = "+areSame(l2, l3, idKSubTest, 2, 3, backgroundFrequency));
//               System.out.println();
//               System.out.println("BTW 4 and 22 = "+areSame(l4, l5, idKSubTest, 4, 5, backgroundFrequency));
//               System.out.println();
//               System.out.println("BTW 4 and 23 = "+areSame(l4, l6, idKSubTest, 4, 6, backgroundFrequency));
//               System.out.println();
	}
	 
	 public static boolean hasEdge(Map<List<Integer>, Edge> edglist, int i, int j) {
		return edglist.containsKey(Arrays.asList(i, j));
	 }
	
//	public static LinkedHashMap <KPart, Integer>  getKMer(String sequence) {
//		LinkedHashMap <KPart, Integer> retVal = new LinkedHashMap();
//		 for (int i = 0; i < 4; i++) {//ABCDEFG
//			 String s = sequence.substring(i, i+4);
//			 retVal.put(new KPart(s, new int[]{0,0,0}), i);
//			 String ss = s.replace(s.charAt(2), '-');
//			 retVal.put(new KPart(ss, new int[]{0,0,0}), i);
//			 retVal.put(new KPart(s.replace(s.charAt(1), '-'), new int[]{0,0,0}), i);
//			 retVal.put(new KPart(ss.replace(s.charAt(1), '-'), new int[]{0,0,0}), i);
//		 }
//		 return retVal;
//	 }
	
	public static LinkedHashMap <KPart, Integer>  getKMer(String sequence) {
		LinkedHashMap <KPart, Integer> retVal = new LinkedHashMap();
		 for (int i = 0; i < 4; i++) {//ABCDEFG
			 
			 StringBuilder sb = new StringBuilder();
             StringBuilder sb1 = new StringBuilder();
          //   StringBuilder sb2= new StringBuilder();
             String s = sequence.substring(i, i+4);
             //System.out.println(" hmmm ... " + sequence);
             sb.append(s.charAt(0)); sb.append('-'); sb.append(s.charAt(2));  sb.append(s.charAt(3));
             sb1.append(s.charAt(0)); sb1.append(s.charAt(1)); sb1.append('-');  sb1.append(s.charAt(3));
         //    sb2.append(s.charAt(0)); sb2.append('-'); sb2.append('-');  sb2.append(s.charAt(3));
			 
			 
			 retVal.put(new KPart(s, new int[]{0,0,0}), i);
			 retVal.put(new KPart(sb.toString(), new int[]{0,0,0}), i);
			 retVal.put(new KPart(sb1.toString(), new int[]{0,0,0}), i);
			// retVal.put(new KPart(sb2.toString(), new int[]{0,0,0}), i);
		 }
		 return retVal;
	 }
	
//	public static List<KPart>  getKMerList(String sequence) {
//		List<KPart> retVal = new ArrayList();
//		 for (int i = 0; i < 4; i++) {//ABCDEFG
//			 String s = sequence.substring(i, i+K_PART_SIZE);
//			 retVal.add(new KPart(s, new int[]{0,0,0}));
//			 String ss = s.replace(s.charAt(2), '-');
//			 retVal.add(new KPart(ss, new int[]{0,0,0}));
//			 retVal.add(new KPart(s.replace(s.charAt(1), '-'), new int[]{0,0,0}));
//			 retVal.add(new KPart(ss.replace(s.charAt(1), '-'), new int[]{0,0,0}));
//		 }
//		 return retVal;
//	 }
	public static Set<KPart>  getKMerList(String sequence) {
        Set<KPart> retVal = new HashSet();
         for (int i = 0; i < 4; i++) {//ABCDEFG
                 StringBuilder sb = new StringBuilder();
                 StringBuilder sb1 = new StringBuilder();
                // StringBuilder sb2= new StringBuilder();
                 String s = sequence.substring(i, i+4);
                 sb.append(s.charAt(0)); sb.append('-'); sb.append(s.charAt(2));  sb.append(s.charAt(3));
                 sb1.append(s.charAt(0)); sb1.append(s.charAt(1)); sb1.append('-');  sb1.append(s.charAt(3));
                // sb2.append(s.charAt(0)); sb2.append('-'); sb2.append('-');  sb2.append(s.charAt(3));

                 retVal.add(new KPart(s, new int[]{0,0,0}));
                 retVal.add(new KPart(sb.toString(), new int[]{0,0,0}));
                 retVal.add(new KPart(sb1.toString(), new int[]{0,0,0}));
               //  retVal.add(new KPart(sb2.toString(), new int[]{0,0,0}));

         }
         return retVal;
 }
	
	public static String conv2Str(KPart kpart) {
    	String str = "";
    	for (int i = 0; i < K_PART_SIZE; i++) {
    		str += kpart.amino.charAt(i);
    		if (i < 3) {
    			for (int j = 0; j < kpart.gaps[i]; j++) {
    				str += '-';
    			}
    		}
    	}
    	return str;
    }
	
	public static Set<List<Integer>> findPath(Set<Integer> vertices, Map<Integer, List<Integer>> in, Map<Integer, List<Integer>> out, Map<List<Integer>, Edge> edgeList){
		Set<List<Integer>> path = new HashSet<>();
		for(Integer u:vertices) {
			if(in.containsKey(u) && out.containsKey(u)) {
				for(int x :in.get(u)) {
					if(edgeList.get(Arrays.asList(x, u)).weight >= edgeSupport) {
						for (int v: out.get(u)) {
							if (edgeList.get(Arrays.asList(u, v)).weight >= edgeSupport) {
								path.add(Arrays.asList(x, u, v));
							}
						}
					}
				}
			}
			
			
		}
		return path;
	}
	
	public static List<String> alignPeptides(List<String> peps, KPart u, KPart x, KPart v){
		List <String> retVal = new ArrayList();
		LinkedHashMap<KPart, Integer> kpmap = new LinkedHashMap();
		int[] temp = new int [3];
		for(String s:peps) {
			StringBuilder sb = new StringBuilder();
			kpmap = getKMer(s);
			if(kpmap.containsKey(u) && kpmap.containsKey(x) && kpmap.containsKey(v)) {
				temp[0] = kpmap.get(u);
				temp[1] = kpmap.get(x);
				temp[2] = kpmap.get(v);
				break;
			}
		}
		for(String s:peps) {
			StringBuilder sb = new StringBuilder();
			kpmap = getKMer(s);
			if(kpmap.containsKey(u) && kpmap.containsKey(x)) {
				if (kpmap.get(u) < temp[0]) {
					sb.append('-');
					sb.append('-');
					sb.append('-');
					sb.append(s);
					sb.append('-');
					retVal.add(sb.toString());
				}else if(kpmap.get(u) == temp[0]) {
					sb.append('-');
					sb.append('-');
					sb.append(s);
					sb.append('-');
					sb.append('-');
					retVal.add(sb.toString());
				}else if(kpmap.get(u) - 1 == temp[0]){
					sb.append('-');
					sb.append(s);
					sb.append('-');
					sb.append('-');
					sb.append('-');
					retVal.add(sb.toString());
				}
				else {
					sb.append(s);
					sb.append('-');
					sb.append('-');
					sb.append('-');
					sb.append('-');
					retVal.add(sb.toString());
				}
			}else {
				if (kpmap.get(x) - temp[1] == -2) {
					sb.append('-');
					sb.append('-');
					sb.append('-');
					sb.append('-');
					sb.append(s);
					retVal.add(sb.toString());
				}else if(kpmap.get(x) - temp[1] == -1) {
					sb.append('-');
					sb.append('-');
					sb.append('-');
					sb.append(s);
					sb.append('-');
					retVal.add(sb.toString());
				}else if (kpmap.get(x) - temp[1] == 0){
					sb.append('-');
					sb.append('-');
					sb.append(s);
					sb.append('-');
					sb.append('-');
					retVal.add(sb.toString());
				}else if(kpmap.get(x) - temp[1] == 1) {
					sb.append('-');
					sb.append(s);
					sb.append('-');
					sb.append('-');
					sb.append('-');
					retVal.add(sb.toString());
				}else {
					sb.append(s);
					sb.append('-');
					sb.append('-');
					sb.append('-');
					sb.append('-');
					retVal.add(sb.toString());
				}
			}
			
		}
		
		
		
		return retVal;
	}
	
	public static double[][] makeProfile(List<String> alignedPeptides){
		int size = alignedPeptides.size();
		if(alignedPeptides.size()==0) {
			//System.out.println("alignedPeptides size = 0");
			double[][] urgentExit = new double [20][7];
			return urgentExit;
		}
		int seqLength = alignedPeptides.get(0).length();
		double[][] pwm= new double[20][seqLength];
		for (int i = 0; i < seqLength; i++) {
			double vA = 0; double vR = 0; double vN = 0; double vD = 0;
			double vC = 0; double vQ = 0; double vE = 0; double vG = 0;
			double vH = 0; double vI = 0; double vL = 0; double vK = 0;
			double vM = 0; double vF = 0; double vP = 0; double vS = 0;
			double vT = 0; double vW = 0; double vY = 0; double vV = 0;
            
            for (String s : alignedPeptides) {
                switch (s.charAt(i)) {
                    case 'A':
                        vA += 1;
                        break;
                    case 'R':
                        vR += 1;
                        break;
                    case 'N':
                        vN += 1;
                        break;
                    case 'D':
                        vD += 1;
                        break;
                    case 'C':
                        vC += 1;
                        break;
                    case 'Q':
                        vQ += 1;
                        break;
                    case 'E':
                        vE += 1;
                        break;
                    case 'G':
                        vG += 1;
                        break;
                    case 'H':
                        vH += 1;
                        break;
                    case 'I':
                        vI += 1;
                        break;
                    case 'L':
                        vL += 1;
                        break;
                    case 'K':
                        vK += 1;
                        break;
                    case 'M':
                        vM += 1;
                        break;
                    case 'F':
                        vF += 1;
                        break;
                    case 'P':
                        vP += 1;
                        break;
                    case 'S':
                        vS += 1;
                        break;
                    case 'T':
                        vT += 1;
                        break;
                    case 'W':
                        vW += 1;
                        break;
                    case 'Y':
                        vY += 1;
                        break;
                    case 'V':
                        vV += 1;
                        break;
                }
            }
         // We are using psuedo counts = 1 in calculating PPM (position probability matrix)
            pwm[0][i] = (vA+0.05)/(size+1); pwm[1][i] = (vR+0.05)/(size+1); pwm[2][i] = (vN+0.05)/(size+1); pwm[3][i] = (vD+0.05)/(size+1);
            pwm[4][i] = (vC+0.05)/(size+1); pwm[5][i] = (vQ+0.05)/(size+1); pwm[6][i] = (vE+0.05)/(size+1); pwm[7][i] = (vG+0.05)/(size+1);
            pwm[8][i] = (vH+0.05)/(size+1); pwm[9][i] = (vI+0.05)/(size+1); pwm[10][i] = (vL+0.05)/(size+1);pwm[11][i] = (vK+0.05)/(size+1);
            pwm[12][i] = (vM+0.05)/(size+1); pwm[13][i] = (vF+0.05)/(size+1); pwm[14][i] = (vP+0.05)/(size+1); pwm[15][i] = (vS+0.05)/(size+1);
            pwm[16][i] = (vT+0.05)/(size+1); pwm[17][i] = (vW+0.05)/(size+1); pwm[18][i] = (vY+0.05)/(size+1); pwm[19][i] = (vV+0.05)/(size+1);
		}
		return pwm;
	}
//	public static void printMatrix(double[][] matrix) {
//		int i = 0;
//		int j = 1;
//		System.out.println("Position   " + j + "    "+ j+1  + "    "+j+2+ "    "+j+3);
//		System.out.println("-------------------------------------------------------");
//		for (double[] row : matrix) {
//			System.out.print(AminoAcid[i] + ": ");
//			i++;
//			for (double element : row) {
//				System.out.printf("%5.2f", element);
//			}
//			System.out.println();
//		}
//	}
	
	public static void printMatrix(double[][] matrix) {
		int i = 0;
		int j = 1;
		System.out.println("Pos   " +1+"     "+2+"     "+3+"     "+4+"     "+5+"     "+6+"     "+7+"     "+8+"     "+9+"     "+10+"    "+11);
		System.out.println("----------------------------------------------------------------------");
		
		for (double[] row : matrix) {
			System.out.print(AminoAcid[i] + ": ");
			i++;
			for (double element : row) {
				System.out.printf("%6.2f", element);
			}
			System.out.println();
		}
	}
	public static HashMap <Character , Double> backgroundFreq(List<String> list){
		HashMap <Character , Double> map = new HashMap();
		double vA = 0; double vR = 0; double vN = 0; double vD = 0;
		double vC = 0; double vQ = 0; double vE = 0; double vG = 0;
		double vH = 0; double vI = 0; double vL = 0; double vK = 0;
		double vM = 0; double vF = 0; double vP = 0; double vS = 0;
		double vT = 0; double vW = 0; double vY = 0; double vV = 0;
		
		for(String s:list) {
			for(int i = 0; i < s.length(); i++) {
				switch (s.charAt(i)) {
                case 'A':
                    vA += 1;
                    break;
                case 'R':
                    vR += 1;
                    break;
                case 'N':
                    vN += 1;
                    break;
                case 'D':
                    vD += 1;
                    break;
                case 'C':
                    vC += 1;
                    break;
                case 'Q':
                    vQ += 1;
                    break;
                case 'E':
                    vE += 1;
                    break;
                case 'G':
                    vG += 1;
                    break;
                case 'H':
                    vH += 1;
                    break;
                case 'I':
                    vI += 1;
                    break;
                case 'L':
                    vL += 1;
                    break;
                case 'K':
                    vK += 1;
                    break;
                case 'M':
                    vM += 1;
                    break;
                case 'F':
                    vF += 1;
                    break;
                case 'P':
                    vP += 1;
                    break;
                case 'S':
                    vS += 1;
                    break;
                case 'T':
                    vT += 1;
                    break;
                case 'W':
                    vW += 1;
                    break;
                case 'Y':
                    vY += 1;
                    break;
                case 'V':
                    vV += 1;
                    break;
				}
			}
		}
		map.put('A', vA/(K_MER_SIZE*list.size()));map.put('R', vR/(K_MER_SIZE*list.size()));map.put('N', vN/(K_MER_SIZE*list.size()));map.put('D', vD/(K_MER_SIZE*list.size()));map.put('C', vC/(K_MER_SIZE*list.size()));
		map.put('Q', vQ/(K_MER_SIZE*list.size()));map.put('E', vE/(K_MER_SIZE*list.size()));map.put('G', vG/(K_MER_SIZE*list.size()));map.put('H', vH/(K_MER_SIZE*list.size()));map.put('I', vI/(K_MER_SIZE*list.size()));
		map.put('L', vL/(K_MER_SIZE*list.size()));map.put('K', vK/(K_MER_SIZE*list.size()));map.put('M', vM/(K_MER_SIZE*list.size()));map.put('F', vF/(K_MER_SIZE*list.size()));map.put('P', vP/(K_MER_SIZE*list.size()));
		map.put('S', vS/(K_MER_SIZE*list.size()));map.put('T', vT/(K_MER_SIZE*list.size()));map.put('W', vW/(K_MER_SIZE*list.size()));map.put('Y', vY/(K_MER_SIZE*list.size()));map.put('V', vV/(K_MER_SIZE*list.size()));
		
		return map;
	}
	public static double[][] informationContent(double[][] M, int n, HashMap <Character , Double> map){
		double[][] IC = new double[20][n];
		for(int i = 0; i < n; i++) {
			IC[0][i] = Math.max(M[0][i]*(Math.log(M[0][i]/map.get('A'))), 0);
			IC[1][i] = Math.max(M[1][i]*(Math.log(M[1][i]/map.get('R'))), 0);
			IC[2][i] = Math.max(M[2][i]*(Math.log(M[2][i]/map.get('N'))), 0);
			IC[3][i] = Math.max(M[3][i]*(Math.log(M[3][i]/map.get('D'))), 0);
			IC[4][i] = Math.max(M[4][i]*(Math.log(M[4][i]/map.get('C'))), 0);
			IC[5][i] = Math.max(M[5][i]*(Math.log(M[5][i]/map.get('Q'))), 0);
			IC[6][i] = Math.max(M[6][i]*(Math.log(M[6][i]/map.get('E'))), 0);
			IC[7][i] = Math.max(M[7][i]*(Math.log(M[7][i]/map.get('G'))), 0);
			IC[8][i] = Math.max(M[8][i]*(Math.log(M[8][i]/map.get('H'))), 0);
			IC[9][i] = Math.max(M[9][i]*(Math.log(M[9][i]/map.get('I'))), 0);
			IC[10][i] = Math.max(M[10][i]*(Math.log(M[10][i]/map.get('L'))), 0);
			IC[11][i] = Math.max(M[11][i]*(Math.log(M[11][i]/map.get('K'))), 0);
			IC[12][i] = Math.max(M[12][i]*(Math.log(M[12][i]/map.get('M'))), 0);
			IC[13][i] = Math.max(M[13][i]*(Math.log(M[13][i]/map.get('F'))), 0);
			IC[14][i] = Math.max(M[14][i]*(Math.log(M[14][i]/map.get('P'))), 0);
			IC[15][i] = Math.max(M[15][i]*(Math.log(M[15][i]/map.get('S'))), 0);
			IC[16][i] = Math.max(M[16][i]*(Math.log(M[16][i]/map.get('T'))), 0);
			IC[17][i] = Math.max(M[17][i]*(Math.log(M[17][i]/map.get('W'))), 0);
			IC[18][i] = Math.max(M[18][i]*(Math.log(M[18][i]/map.get('Y'))), 0);
			IC[19][i] = Math.max(M[19][i]*(Math.log(M[19][i]/map.get('V'))), 0);
		}
		return IC;
	}
	//This method considers 7 higher columns of IC matrix and returns summation of those 7 columns
	public static double icOfMotif (double[][] ICMatrix) {
		double ic = 0;
		int cols = ICMatrix[0].length; 
		int rows = ICMatrix.length;
		int i = 1;
		Map<Double, Integer> icOfColMap = new HashMap<>();
		Map<Integer, Double> icOfColMap2 = new HashMap<>();
		List<Double> list = new ArrayList();
		for(int ii = 0; ii < cols; ii++) {
			double icinCol = 0;
			for (int j = 0; j < rows; j++) {
				icinCol += ICMatrix[j][ii];
			}
			list.add(icinCol);
			icOfColMap.put(icinCol, i);
			icOfColMap2.put(i, icinCol);
			i++;
		}
		Collections.sort(list, Collections.reverseOrder());
//		Stream<Map.Entry<Integer, Double>> sorted =
//				icOfColMap2.entrySet().stream()
//			       .sorted(Collections.reverseOrder(Map.Entry.comparingByValue()));
//		for (double[] col : ICMatrix) {
//			i++;
//			double icinCol = 0;
//			
//			  for (double element: col) {
//				  icinCol += element;
//				//  ic += element;
//			  }
//			  icOfColMap.put(icinCol, i);
//		  }
//		Map<Double, Integer> sortedicOfColMap = new TreeMap<Double, Integer>(Collections.reverseOrder());
//		sortedicOfColMap.putAll(icOfColMap);
//		List<Double> list = new ArrayList(sortedicOfColMap.keySet());
//		List<Double> list2 = new ArrayList(sorted.);
//		if(list.size()<5) {
//			System.out.println("Size of list = " + list.size());
//			System.out.println(list);
//			System.out.println("pp = " + pp);
//			System.out.println(" size of pp = " + pp.size());
//			System.out.println(Arrays.deepToString(ICMatrix));
//		}
		//We consider 7 highest values of information content in columns.
//		for(int j = 0; j < 5; j++) {
//			ic += list.get(j);
//		}
		for(int j = 0; j < cols; j++) {
			ic += list.get(j);
		}
		return ic;
		
	}
	public static void motifParamComputations (List<String> aligned, HashMap<Character, Double> backgroundFrequency) {
		double[][] pwm = makeProfile(aligned);
	//	printMatrix(pwm);
		System.out.println();
		int n = aligned.get(0).length();
		double[][]IC = informationContent(pwm, n, backgroundFrequency);
	//	printMatrix(IC);
		System.out.println();
		System.out.println(icOfMotif(IC));
	}
	public static boolean areSameMotifs(List<String> l1, List<String> l2, HashMap<Character, Double> backgroundFrequency) {
					List<String> unionList = new ArrayList();
					List<String> s1 = new ArrayList(l1);
					List<String> s2 = new ArrayList(l2);
					HashSet<String> set1 = new HashSet(l1);
					HashSet<String> set1Copy = new HashSet(set1);
					HashSet<String> set2 = new HashSet(l2);
					HashSet<String> set2Copy = new HashSet(set2);
					List<String> ls1 = new ArrayList(set1);
					List<String> ls2 = new ArrayList(set2);
					s1.retainAll(s2);
					if(l1.size()==0 | l2.size()==0) {
						return false;
					}
					if(s1.size() == 0){
						return false;
					}
					List<String> set1Cleaned = motifCleaning(l1, backgroundFrequency);
					List<String> set2Cleaned = motifCleaning(l2, backgroundFrequency);
//					List<String> set1Raw = rawListFromAlignedList(set1Cleaned);
//					List<String> set2Raw = rawListFromAlignedList(set2Cleaned);
//
//					if(set1Raw.containsAll(set2Raw) || set2Raw.containsAll(set1Raw)) {
//						return true;
//					}
					unionList.addAll(l1);
					unionList.addAll(l2);
					Set<String> unionSet = new HashSet(unionList);
					List<String> union = new ArrayList(unionSet);					
					
					List<List<String>> rawListUnion = new ArrayList();
					for(int i = 0 ; i < union.size()-1; i++) {
						rawListUnion.add(makeRawAlignment(union.get(i), union.get(i+1)));
					}
					List<List<String>> finalListUnion = finalAlignment(rawListUnion);
					List<String> listUnion = alignList(finalListUnion);
					List<String> listUnionCleaned = motifCleaning(listUnion, backgroundFrequency);
					
//					System.out.println("   * LIST UNION * ");
//					for (int i = 0; i < listUnion.size(); i++) {
//						System.out.println(listUnion.get(i));
//					}
					
					double[][] pUnion = makeProfile(listUnionCleaned);
					int n = listUnionCleaned.get(0).length();
		            double[][]ICUnion = informationContent(pUnion, n, backgroundFrequency);
		            double icUnion = icOfMotif(ICUnion);
		            
		            List<List<String>> rawListl1 = new ArrayList();
					for(int i = 0 ; i < ls1.size()-1; i++) {
						rawListl1.add(makeRawAlignment(ls1.get(i), ls1.get(i+1)));
					}
					List<List<String>> finalListl1 = finalAlignment(rawListl1);
					List<String> listl1 = alignList(finalListl1);
					List<String> listl1Cleaned = motifCleaning(listl1, backgroundFrequency);
					List<String> raw1 = rawListFromAlignedList(listl1Cleaned);

					
		            double[][] p1 = makeProfile(listl1Cleaned);
		            int n1 = listl1Cleaned.get(0).length();
		            double[][]IC1 = informationContent(p1, n1, backgroundFrequency);
		            double ic1 = icOfMotif(IC1);
		            
		            List<List<String>> rawListl2 = new ArrayList();
					for(int i = 0 ; i < ls2.size()-1; i++) {
						rawListl2.add(makeRawAlignment(ls2.get(i), ls2.get(i+1)));
					}
					List<List<String>> finalListl2 = finalAlignment(rawListl2);
					List<String> listl2 = alignList(finalListl2);
					List<String> listl2Cleaned = motifCleaning(listl2, backgroundFrequency);
					
					List<String> raw2 = rawListFromAlignedList(listl2Cleaned);
					if(raw1.containsAll(raw2) || raw2.containsAll(raw1)) {
						return true;
					}
					
		            double[][] p2 = makeProfile(listl2Cleaned);
		            int n2 = listl2Cleaned.get(0).length();
		            double[][]IC2 = informationContent(p2, n2, backgroundFrequency);
		            double ic2 = icOfMotif(IC2);
//		            System.out.println("ic1 = " + ic1);
//		            System.out.println("ic2 = " + ic2);
//		            System.out.println("icUnion = " + icUnion);            
					if(icUnion >= ic1 || icUnion >= ic2) {
						return true;
					}else return false;
				}
	
	public static boolean areSame(List<String> l1, List<String> l2, Map <Integer, String> idKSubMap, int i, int j, 
			HashMap<Character, Double> backgroundFrequency) {
		
		if(idKSubMap.get(i).equals(idKSubMap.get(j))) {
//			System.out.println("     M1 = " + idKSubMap.get(i) + "  found ");
//			System.out.println("     M2 = " + idKSubMap.get(j) + "  found ");
			return true;
		}
		List<String> l1Copy = new ArrayList(l1);
		List<String> l2Copy = new ArrayList(l2);
		l1Copy.retainAll(l2Copy);
		 if(l1Copy.size()==0) return false; //jaye bahs dare. 2 with 5 is a counter example
		
		List<String> list1 = new ArrayList();
		for(String s:l1) {
			List<String> windowList = window(s, idKSubMap.get(i).length());
			list1.add(bestShift(idKSubMap.get(i), windowList));
		}
		Set<String> s1 = new HashSet(list1);
		List<String> ls1 = new ArrayList(s1);
//		System.out.println("+++  LS1 +++ ");
//		for(String s:ls1) {
//			System.out.println(s);
//		}
		List<String> ls1Cleaned = motifCleaning(ls1, backgroundFrequency);
//		System.out.println("+++  LS1  CLEANED +++ ");
//		for(String s:ls1Cleaned) {//inja
//			System.out.println(s);
//		}
		
		List<String> list2 = new ArrayList();
		for(String s:l2) {
			List<String> windowList = window(s, idKSubMap.get(j).length());
			list2.add(bestShift(idKSubMap.get(j), windowList));
		}
		Set<String> s2 = new HashSet(list2);
		List<String> ls2 = new ArrayList(s2);
		List<String> ls1Copy = new ArrayList(ls1);List<String> ls1Copyy = new ArrayList(ls1);
		List<String> ls2Copy = new ArrayList(ls2);List<String> ls2Copyy = new ArrayList(ls2);
		ls1Copy.retainAll(ls2Copy); ls2Copyy.retainAll(ls1Copyy);
		if(ls1Copy==ls1) return false;
		if(ls2Copyy==ls2) return false;
//		if(ls1Copy.size() > 10) return true;
//		System.out.println();
//		System.out.println("+++  LS2 +++ ");
//		for(String s:ls2) {
//			System.out.println(s);
//		}
		
		List<String> ls2Cleaned = motifCleaning(ls2, backgroundFrequency);
//		System.out.println("+++  LS2  CLEANED +++ ");
//		for(String s:ls2Cleaned) { // inja
//			System.out.println(s);
//		}
//		System.out.println("+++  LS2  CLEANED +++ ");
//		for(String s:ls2Cleaned) {
//			System.out.println(s);
//		}
//		System.out.println("idKSubMap.get(i) = "+idKSubMap.get(i));
//		System.out.println("idKSubMap.get(j) = "+idKSubMap.get(j));
//		System.out.println("motifKSubMap.get(l1) = " + motifKSubMap.get(l1));
//		System.out.println("motifKSubMap.get(l2) = " + motifKSubMap.get(l2));
//		System.out.println();

//		System.out.println("list1 = " );
//		System.out.println();
//		for(String s:ls1) {
//			System.out.println(s);
//		}
//		System.out.println();
//		System.out.println("list2 = " );
//		System.out.println();
//		for(String s:ls2) {
//			System.out.println(s);
//		}
		//System.out.println("list1 = " + list1);
		// You may apply motifCleaning on l1 and l2
		//List<String> set1Cleaned = motifCleaning(l1, backgroundFrequency);
		//List<String> set2Cleaned = motifCleaning(l2, backgroundFrequency);
		List<String> ls1Cleaned_Copy = new ArrayList(ls1Cleaned);
		List<String> ls2Cleaned_Copy = new ArrayList(ls2Cleaned);
		System.out.println("Size of ls1Cleaned = " + ls1Cleaned.size());
		System.out.println();
		System.out.println("Size of ls2Cleaned = " + ls2Cleaned.size());
		ls1Cleaned_Copy.retainAll(ls2Cleaned_Copy);
		System.out.println("Size of ls1Cleaned = " + ls1Cleaned_Copy.size());
		System.out.println();
		//if (ls1Cleaned_Copy.size() == 0) return false;
		List<String> listUnion = new ArrayList();
		System.out.println("M1 = " + idKSubMap.get(i));
		System.out.println("M2 = " + idKSubMap.get(j));
		String m3 = kSubAlignment(idKSubMap.get(i), idKSubMap.get(j));		
		System.out.println("M3 = " + m3);
		if (kSubAlignmentCalc(idKSubMap.get(i), idKSubMap.get(j)) >= 3 && ls1Cleaned_Copy.size() == 0) {
			return true;
		}else if (kSubAlignmentCalc(idKSubMap.get(i), idKSubMap.get(j)) < 3 && ls1Cleaned_Copy.size() == 0) {
			return false;
		}
		for(String s:l1) {
			List<String> windowList = window(s, m3.length());
			listUnion.add(bestShift2(m3, windowList));
		}
		for(String s:l2) {
			List<String> windowList = window(s, m3.length());
			listUnion.add(bestShift2(m3, windowList));
		}

		Set<String> setUnion = new HashSet(listUnion);
		List<String> unionList = new ArrayList(setUnion);
		List<String> unionListCleaned = motifCleaning (unionList, backgroundFrequency);
//		System.out.println("unionList = " );
//		System.out.println();
//		System.out.println("+++  unionListCleaned  CLEANED +++ ");
//
//		for(String s:unionListCleaned) {
//			System.out.println(s);
//		}
	//	List<String> listUnionCleaned = motifCleaning(listUnion, backgroundFrequency);

		double[][] pUnion = makeProfile(unionListCleaned);
		int n = unionList.get(0).length();
//		System.out.println("  n  ===  " + n);
//		System.out.println("  unionListCleaned  ===  " + unionListCleaned.size());
//		System.out.println();
//		System.out.println("  unionList  " + unionList);
//		System.out.println();
//		System.out.println("  unionListCleaned  " + unionListCleaned);
        double[][]ICUnion = informationContent(pUnion, n, backgroundFrequency);
        double icUnion = icOfMotif(ICUnion);
 	
        double[][] p1;
        if (ls1Cleaned.size() == 0) {
        	 p1 = makeProfile(ls1);
        }else {
        	 p1 = makeProfile(ls1Cleaned);
        }
     //   double[][] p1 = makeProfile(ls1Cleaned); // inja
        int n1 = list1.get(0).length();
        double[][]IC1 = informationContent(p1, n1, backgroundFrequency);
        double ic1 = icOfMotif(IC1);
        
        double[][] p2;
        if (ls2Cleaned.size() == 0) {
        	 p2 = makeProfile(ls2);
        }else {
        	 p2 = makeProfile(ls2Cleaned);
        }
        //double[][] p2 = makeProfile(ls2Cleaned); // inja
        int n2 = list2.get(0).length();
//        System.out.println("check the error1 = " + ls2Cleaned);
//        System.out.println(" n2 = " + n2);
//        System.out.println("check the error2 = " + ls2);
//        System.out.println("check the error3 = " + list2);
        double[][]IC2 = informationContent(p2, n2, backgroundFrequency);
        double ic2 = icOfMotif(IC2);
//        System.out.println("Reached this check point");
//        System.out.println("ic1 = " + ic1);
//        System.out.println("ic2 = " + ic2);
//        System.out.println("icUnion = " + icUnion);
        
		if(icUnion  >= Math.max(ic1, ic2)  ) {
			//if((icUnion + 1 > Math.max(ic1, ic2)) ||( icUnion > ic1 )|| (icUnion > ic2)) {
//			System.out.println("List 1 = ");
//			for(String s:ls1) {
//				System.out.println(s);
//			}
//			System.out.println();
//			System.out.println("List 2 = ");
//			for(String s:ls2) {
//				System.out.println(s);
//			}
//			System.out.println("!@#$%^&*!@#$%^&*!@#$%^&*");
			return true;
			}else if((icUnion < ic2) || ( icUnion < ic1)) { 
				if (kSubAlignmentCalc (idKSubMap.get(i), idKSubMap.get(j)) >= 3){
					return true;
			}return false;
		
//		}else if(icUnion > ic2 & icUnion < ic1) { //drop list 2
		
		}
		else return false;
	
		}
	
	public static boolean areSameMotifs2 (List<String> l1, List<String> l2, HashMap<Character, Double> backgroundFrequency) {
		List<List<String>> rawListl1 = new ArrayList();
		for(int i = 0 ; i < l1.size()-1; i++) {
			rawListl1.add(makeRawAlignment(l1.get(i), l1.get(i+1)));
		}
		List<List<String>> finalListl1 = finalAlignment(rawListl1);
		List<String> lis1 = alignList(finalListl1);
		Set<String> s1 = new HashSet(lis1);
		List<String> list1 = new ArrayList(s1);
		List<String> list1Cleaned = motifCleaning(list1, backgroundFrequency);
		List<String> raw1 = rawListFromAlignedList(list1Cleaned);
		
		List<List<String>> rawListl2 = new ArrayList();
		for(int i = 0 ; i < l2.size()-1; i++) {
			rawListl2.add(makeRawAlignment(l2.get(i), l2.get(i+1)));
		}
		List<List<String>> finalListl2 = finalAlignment(rawListl2);
		List<String> lis2 = alignList(finalListl2);
		Set<String> s2 = new HashSet(lis2);
		List<String> list2 = new ArrayList(s2);
		List<String> list2Cleaned = motifCleaning(list2, backgroundFrequency);
		List<String> raw2 = rawListFromAlignedList(list2Cleaned);
		if(raw1.containsAll(raw2) || raw2.containsAll(raw1)) return true;
		
		//icUnion calculations
		List<String> u = new ArrayList();
		u.addAll(raw1);
		u.addAll(raw2);
		Set<String> s = new HashSet(u);
		List<String> lisU = new ArrayList(s);
		
		List<List<String>> rawListUnion = new ArrayList();
		for(int i = 0 ; i < lisU.size()-1; i++) {
			rawListUnion.add(makeRawAlignment(lisU.get(i), lisU.get(i+1)));
		}
		List<List<String>> finalListlUnion = finalAlignment(rawListUnion);
		List<String> listU = alignList(finalListlUnion);
		Set<String> sU = new HashSet(listU);
		List<String> listUU = new ArrayList(sU);
		//List<String> listUCleaned = motifCleaning(listUU, backgroundFrequency);
		//double[][] pUnion = makeProfile(listUCleaned);
		double[][] pUnion = makeProfile(listUU);
	//	int n = listUCleaned.get(0).length();
		int n = listUU.get(0).length();
        double[][]ICUnion = informationContent(pUnion, n, backgroundFrequency);
        double icUnion = icOfMotif(ICUnion);
        
        //ic1 calculations		
		
        double[][] p1 = makeProfile(list1);
        int n1 = list1.get(0).length();
        double[][]IC1 = informationContent(p1, n1, backgroundFrequency);
        double ic1 = icOfMotif(IC1);
		
        //ic2 calculations		
		
        double[][] p2 = makeProfile(list2);
        int n2 = list2.get(0).length();
        double[][]IC2 = informationContent(p2, n2, backgroundFrequency);
        double ic2 = icOfMotif(IC2);

//        System.out.println("ic1 = " + ic1);
//        System.out.println("ic2 = " + ic2);
//        System.out.println("icUnion = " + icUnion);
        
		if(icUnion >= ic1 || icUnion >= ic2) {
			return true;
		}else 
			return false;
	}
	
	
	
	public static Set<String> combineMotifs(List<List<String>> inputMotifs, HashMap<Character, Double> backgroundFrequency, Map<String, Set<String>> motifIdToMotif) {
		if (inputMotifs.isEmpty()) {
			return new HashSet<>();
		}
		Map<Integer, String> motifIndexToMotifId = new HashMap<>();

		Set<String> allCombinedMotifs = new HashSet<>();
		
		for (int i = 0; i < inputMotifs.size(); i++) {
			motifIndexToMotifId.put(i, String.valueOf(i));
			motifIdToMotif.put(String.valueOf(i), new HashSet<>(inputMotifs.get(i)));
		}
//		System.out.println("Input motif size is = " + inputMotifs.size());
//		System.out.println("Printing . . ." + motifIndexToMotifId);
		Set<Integer> skip = new HashSet<>();
		
		for (int i = 0; i < inputMotifs.size(); i++) {
			if(skip.contains(i)) {
				continue;
			}
			for (int j = i + 1; j < inputMotifs.size(); j++) {
//				if(i==5) {
//					return allCombinedMotifs;
//				}
				
				String motifId1 = motifIndexToMotifId.get(i);
				String motifId2 = motifIndexToMotifId.get(j);
				Set<String> motif1 = motifIdToMotif.get(motifId1);
				Set<String> motif2 = motifIdToMotif.get(motifId2);
//				System.out.println("motifIndexToMotifIdMAP");
//				System.out.println(motifIndexToMotifId);
//				System.out.println("motifIndexToMotifIdMAP .... Finish");
				
			System.out.println("Size of motifID1 = "+ motifId1.length()+ " " + " ----i = " + i);
//				System.out.println("motif 1 = "+ motif1 );
//				System.out.println();
			System.out.println("Size of motifID2 = "+ motifId2.length()+ " " + " ----j = " + j);
//				System.out.println("motif 2 = "+ motif2 );
				
				if (areSameMotifs2(new ArrayList(motif1), new ArrayList(motif2), backgroundFrequency)) {
					System.out.println("Combining two Motif . . . ");
					Set<String> motifUnion = combineTwoMotifs(motifId1, motifId2, motif1, motif2, allCombinedMotifs);
					String unionId =  motifId1 +"_"+ motifId2;
					motifIndexToMotifId.put(i, unionId);
					motifIndexToMotifId.put(j, unionId);
					motifIdToMotif.put(unionId, motifUnion);
					skip.add(j);
				} else {
					allCombinedMotifs.add(motifId1);
					allCombinedMotifs.add(motifId2);
				}
			}
		}
		return allCombinedMotifs;
	}

	private static Set<String> combineTwoMotifs(String motifId1, String motifId2, Set<String> motif1, Set<String> motif2, 
														Set<String> allCombinedMotifs) {
		// combine the two motifs and put it in the map allCombinedMotifs
		String unionMotifId = motifId1+"_" + motifId2;
		Set<String> motifUnion = new HashSet<>();
		motifUnion.addAll(motif1);
		motifUnion.addAll(motif2);
		allCombinedMotifs.add(unionMotifId);

		// remove existing motifId1 and motifId2 from map allCombinedMotifs
		if (allCombinedMotifs.contains(motifId1)) {
			allCombinedMotifs.remove(motifId1);
		}
		if (allCombinedMotifs.contains(motifId2)) {
			allCombinedMotifs.remove(motifId2);
		}

		return motifUnion;
	}

	public static int calculateScore(String s1, String s2) {
		int score = 0;
		for (int i = 0; i < s1.length(); i++) {
			if (s1.charAt(i) == s2.charAt(i) & s1.charAt(i) != '-') {
				score += 7;
			}else if(s1.charAt(i) == s2.charAt(i) & s1.charAt(i) == '-') {
				score++;
			}
		}
		return score;
	}

	public static String preShift(String s) {
		StringBuilder sb = new StringBuilder();
		sb = sb.append("-");
		return sb.toString()+s;
	}
	
	public static String postShift(String s) {
		StringBuilder sb = new StringBuilder();
		sb = sb.append("-");
		return s + sb.toString();
	}
	
	public static List<String> makeRawAlignment(String s1, String s2) {
		int score1 = calculateScore(s1,s2);
		int score2 = calculateScore(preShift(s1), postShift(s2));
		int score3 = calculateScore(preShift(preShift(s1)), postShift(postShift(s2)));
		int score4 = calculateScore(preShift(preShift(preShift(s1))), postShift(postShift(postShift(s2))));
		int score5 = calculateScore(preShift(preShift(preShift(preShift(s1)))), postShift(postShift(postShift(postShift(s2)))));
		
		int score6 = calculateScore(preShift(s2), postShift(s1));
		int score7 = calculateScore(preShift(preShift(s2)), postShift(postShift(s1)));
		int score8 = calculateScore(preShift(preShift(preShift(s2))), postShift(postShift(postShift(s1))));
		int score9 = calculateScore(preShift(preShift(preShift(preShift(s2)))), postShift(postShift(postShift(postShift(s1)))));
		
		List<Integer> list = Arrays.asList(score1, score2, score3, score4, score5, score6, score7, score8, score9);
		Collections.sort(list, Collections.reverseOrder());
		List<String> alignedButNotequal = new ArrayList();
		if (list.get(0)==score1) {
			alignedButNotequal.add(s1); alignedButNotequal.add(s2);
		}else if(list.get(0) == score2) {
			alignedButNotequal.add(preShift(s1)); alignedButNotequal.add(postShift(s2));
		}else if(list.get(0) == score3) {
			alignedButNotequal.add(preShift(preShift(s1))); alignedButNotequal.add(postShift(postShift(s2))); 
		}else if(list.get(0) == score4) {
			alignedButNotequal.add(preShift(preShift(preShift(s1)))); alignedButNotequal.add(postShift(postShift(postShift(s2)))); 
		}else if(list.get(0) == score5) {
			alignedButNotequal.add(preShift(preShift(preShift(preShift(s1))))); alignedButNotequal.add(postShift(postShift(postShift(postShift(s2))))); 
		}else if(list.get(0) == score6) {
			alignedButNotequal.add(postShift(s1)); alignedButNotequal.add(preShift(s2));
		}else if(list.get(0) == score7) {
			alignedButNotequal.add(postShift(postShift(s1))); alignedButNotequal.add(preShift(preShift(s2))); 
		}else if(list.get(0) == score8) {
			alignedButNotequal.add(postShift(postShift(postShift(s1)))); alignedButNotequal.add(preShift(preShift(preShift(s2))));
		}else if(list.get(0) == score9) {
			alignedButNotequal.add(postShift(postShift(postShift(postShift(s1))))); alignedButNotequal.add(preShift(preShift(preShift(preShift(s2)))));
		}
		return alignedButNotequal;
	}
	
	public static List<List<String>> finalAlignment(List<List<String>> alignedButNotequal){
		List<List<String>> finalAlignedPairs = new ArrayList(alignedButNotequal);
		for (int i = 0; i < alignedButNotequal.size()-1; i++) {
			String str0 = finalAlignedPairs.get(i).get(0);
			String str1 = finalAlignedPairs.get(i).get(1);
			String str2 = finalAlignedPairs.get(i+1).get(0);
			String str3 = finalAlignedPairs.get(i+1).get(1);
			
			if(dashCountBeforeStart(str1) > dashCountBeforeStart(str2)) {	
				int m = dashCountBeforeStart(str1)-dashCountBeforeStart(str2);			
				int p1 = dashCountAfterEnd(str1)-dashCountAfterEnd(str2);
				int p2 = dashCountAfterEnd(str2)-dashCountAfterEnd(str1);		
				for(int j = 0; j < m; j++) {
					
						str2 = preShift(str2);
						str3 = preShift(str3);
					
				}
				if (p1 > 0) {
					for (int k = 0; k < p1; k++) {
						str2 = postShift(str2);
						str3 = postShift(str3);
					}
				}
				if(p2 > 0 &  p2 < 3) {
					for (int k = 0; k < p2; k++) {
						str0 = postShift(str0);
						str1 = postShift(str1);
					}
					int a = i-1;
					while(a >= 0) {//1. b < p2
						String back1 = finalAlignedPairs.get(a).get(0);
						String back2 = finalAlignedPairs.get(a).get(1);
						for(int b = 0; b < p2; b++) {
							back1 = postShift(back1);
							back2 = postShift(back2);
						}
						finalAlignedPairs.get(a).remove(0);
						finalAlignedPairs.get(a).add(back1);
						finalAlignedPairs.get(a).remove(0);
						finalAlignedPairs.get(a).add(back2);
						a--;
					}
				}
				finalAlignedPairs.get(i).remove(0);
				finalAlignedPairs.get(i).add(str0);
				finalAlignedPairs.get(i).remove(0);
				finalAlignedPairs.get(i).add(str1);
				finalAlignedPairs.get(i+1).remove(0);
				finalAlignedPairs.get(i+1).add(str2);
				finalAlignedPairs.get(i+1).remove(0);
				finalAlignedPairs.get(i+1).add(str3);
			}else if (dashCountBeforeStart(str1) < dashCountBeforeStart(str2)) {
				int m = dashCountBeforeStart(str2)-dashCountBeforeStart(str1);			
				int p1 = dashCountAfterEnd(str1)-dashCountAfterEnd(str2);
				int p2 = dashCountAfterEnd(str2)-dashCountAfterEnd(str1);			
				for(int j = 0; j < m; j++) {
					str1 = preShift(str1);	
					str0 = preShift(str0);
				}		
				int a = i-1;
				while(a >= 0) {//2. b < m
					String back1 = finalAlignedPairs.get(a).get(0);
					String back2 = finalAlignedPairs.get(a).get(1);
					for(int b = 0; b < m; b++) {
						back1 = preShift(back1);
						back2 = preShift(back2);
					}
					finalAlignedPairs.get(a).remove(0);
					finalAlignedPairs.get(a).add(back1);
					finalAlignedPairs.get(a).remove(0);
					finalAlignedPairs.get(a).add(back2);
					a--;
				}			
				if (p1 > 0) {
					for (int k = 0; k < p1; k++) {
						str2 = postShift(str2);
						str3 = postShift(str3);
					}
				}
				if(p2 > 0) {
					for (int k = 0; k < p2; k++) {
						str0 = postShift(str0);
						str1 = postShift(str1);
					}
				}
				finalAlignedPairs.get(i).remove(0);
				finalAlignedPairs.get(i).add(str0);
				finalAlignedPairs.get(i).remove(0);
				finalAlignedPairs.get(i).add(str1);
				finalAlignedPairs.get(i+1).remove(0);
				finalAlignedPairs.get(i+1).add(str2);
				finalAlignedPairs.get(i+1).remove(0);
				finalAlignedPairs.get(i+1).add(str3);
				}else if(dashCountAfterEnd(str1) > dashCountAfterEnd(str2)) {
					int m = dashCountAfterEnd(str1) - dashCountAfterEnd(str2);	
					int p1 = dashCountBeforeStart(str1)-dashCountBeforeStart(str2);
					int p2 = dashCountBeforeStart(str2)-dashCountBeforeStart(str1);
					for(int j = 0; j < m; j++) {
						str2 = postShift(str2);	
						str3 = postShift(str3);
					}
					if (p1 > 0) {
						for (int k = 0; k < p1; k++) {
							str2 = preShift(str2);
							str3 = preShift(str3);
						}
					}
					if(p2 > 0) {
						for (int k = 0; k < p2; k++) {
							str0 = preShift(str0);
							str1 = preShift(str1);
						}
					}		
					finalAlignedPairs.get(i).remove(0);
					finalAlignedPairs.get(i).add(str0);
					finalAlignedPairs.get(i).remove(0);
					finalAlignedPairs.get(i).add(str1);
					finalAlignedPairs.get(i+1).remove(0);
					finalAlignedPairs.get(i+1).add(str2);
					finalAlignedPairs.get(i+1).remove(0);
					finalAlignedPairs.get(i+1).add(str3);
					}else if(dashCountAfterEnd(str2) > dashCountAfterEnd(str1)){
						int m = dashCountAfterEnd(str2) - dashCountAfterEnd(str1);
						int p1 = dashCountBeforeStart(str1)-dashCountBeforeStart(str2);
						int p2 = dashCountBeforeStart(str2)-dashCountBeforeStart(str1);				
						int a = i-1;
						while(a >= 0) {//3. b < m
							String back1 = finalAlignedPairs.get(a).get(0);
							String back2 = finalAlignedPairs.get(a).get(1);
							for(int b = 0; b < m; b++) {
								back1 = postShift(back1);
								back2 = postShift(back2);
							}
							finalAlignedPairs.get(a).remove(0);
							finalAlignedPairs.get(a).add(back1);
							finalAlignedPairs.get(a).remove(0);
							finalAlignedPairs.get(a).add(back2);
							a--;
						}
					if (p1 > 0) {
						for (int k = 0; k < p1; k++) {
							str2 = preShift(str2);
							str3 = preShift(str3);
						}
					}
					if(p2 > 0) {
						for (int k = 0; k < p2; k++) {
							str0 = preShift(str0);
							str1 = preShift(str1);
						}
					}					
					for(int j = 0; j < m; j++) {
						str1 = postShift(str1);
						str0 = postShift(str0);
					}
					finalAlignedPairs.get(i).remove(0);
					finalAlignedPairs.get(i).add(str0);
					finalAlignedPairs.get(i).remove(0);
					finalAlignedPairs.get(i).add(str1);
					finalAlignedPairs.get(i+1).remove(0);
					finalAlignedPairs.get(i+1).add(str2);
					finalAlignedPairs.get(i+1).remove(0);
					finalAlignedPairs.get(i+1).add(str3);
				}
		}	
		return finalAlignedPairs;
	}
	
	public static int dashCountBeforeStart(String s) {
		int dashCount = 0;
		for (int i = 0; i < 4; i++) {
			if (s.charAt(i) == '-') {
				dashCount++;
			}
		}
		return dashCount;
	}
	
	public static int dashCountAfterEnd(String s) {
		int dashCount = 0;
		for (int i = 4; i < s.length(); i++) {
			if (s.charAt(i) == '-') {
				dashCount++;
			}
		}
		return dashCount;
	}
	public static List<String> alignList(List<List<String>> l){
		List<String> list = new ArrayList();
		for (int i = 0; i < l.size(); i++) {
			list.add(l.get(i).get(0));
		}
		list.add(l.get(l.size()-1).get(1));
		return list;
	}
	
	public static List<List<String>> purify(List<List<String>> rawMotifs, HashMap<Character, Double> backgroundFrequency){
		int i = 1;
		while(i < rawMotifs.size()) {
			List<String> curr = rawMotifs.get(i);
			int j = i-1;
			boolean isDeleted = false;
			while(j >= 0) {
				if (areSameMotifs(curr, rawMotifs.get(j), backgroundFrequency)) {				
					System.out.println(i + " with " + j + " ");
					System.out.println("First = " + curr.get(0));
					System.out.println("Second = " + rawMotifs.get(j).get(0));
		            System.out.println("++++++++++++++++++++++++");
					List<String> newList = new ArrayList();
					newList.addAll(rawMotifs.get(j));
					newList.addAll(curr);
					rawMotifs.set(j, newList);
					rawMotifs.remove(i);
					isDeleted = true;
					break;
				}				
				j--;
			}
			if(isDeleted == false) {
				i++;
			}
		}		
		return rawMotifs;
	}
	public static List<List<String>> purify2(List<List<String>> rawMotifs){int i = 1;
	while(i < rawMotifs.size()) {
		List<String> curr = rawMotifs.get(i);
		int j = i-1;
		boolean isDeleted = false;
		while(j >= 0) {
			if (newAreSameMotifs(curr, rawMotifs.get(j))) {				
//				System.out.println(i + " with " + j + " ");
//				System.out.println("First = " + curr.get(0));
//				System.out.println("Second = " + rawMotifs.get(j).get(0));
//	            System.out.println("++++++++++++++++++++++++");
				List<String> newList = new ArrayList();
				newList.addAll(rawMotifs.get(j));
				newList.addAll(curr);
				rawMotifs.set(j, newList);
				rawMotifs.remove(i);
				isDeleted = true;
				break;
			}				
			j--;
		}
		if(isDeleted == false) {
			i++;
		}
	}		
	return rawMotifs;
	}
	
	
//	public static List<List<String>> purify4 (List<List<String>> rawMotifs, Map<List<String>, String> motifKSubMap, 
//			Map<Integer, List<String>> idPepListMap, Map<Integer, String> idKSubMap, HashMap<Character, Double> backgroundFrequency){
//			
//			List< List< String > > retVal = new ArrayList();
//			List< List< Integer > > adjList = new ArrayList();
//			for(int l = 0; l < rawMotifs.size(); l++) {
//				adjList.add(new ArrayList());
//			}
//			Map <Integer, List<Integer>> adjMap = new HashMap();
//			Map <Integer, Integer> adjMapNew = new HashMap();
//			for (int i = 0; i < rawMotifs.size(); i++) {
//				for(int j = 0; j < rawMotifs.size(); j++) {
//					if (i==j) continue;
//					adjList.get(i).add(i);
//					if (areSame(rawMotifs.get(i), rawMotifs.get(j), motifKSubMap, idPepListMap, idKSubMap, i, j, backgroundFrequency)) {
//						System.out.println(i + " with " + j +" is " + " true");
//						adjList.get(i).add(j);
//						
//					}
//				}
//			}
//			
//			for(int i = 0; i < adjList.size(); i++) {
//				List<Integer> list = adjList.get(i);
//				boolean flag = false;
//				for(Integer j:list) {
//					if(adjMapNew.containsKey(j)) {
//						flag = true;
//						for(Integer k: list) {
//							adjMapNew.put(k, adjMapNew.get(j));
//						}
//						break;
//					}
//					adjMapNew.putIfAbsent(j,i);
//				}
//				
//			}
////			System.out.println("adjMapNew = " + adjMapNew);
////			System.out.println();
////			System.out.println();
////			System.out.println("adjList = " + adjList);
////			for (int k = 0; k < adjList.size(); k++) {
////				adjMap.put(k, adjList.get(k));
////			}
//			Set<Set <Integer>> set = new HashSet();
//			for(Entry<Integer, Integer> entry1: adjMapNew.entrySet()) {
//				Set<Integer> temp = new HashSet();
//				
//				for(Entry<Integer, Integer> entry2: adjMapNew.entrySet()) {
//					if (entry1.getKey()==entry2.getKey()) continue;
//					if(entry1.getValue() == entry2.getValue()) {
//						temp.add(entry1.getKey());
//						temp.add(entry2.getKey());
//					}else {
//						temp.add(entry1.getKey());
//					}
//					
//				}
//				set.add(temp);
//			}
////			for (Entry <Integer, List < Integer > > entry:adjMap.entrySet()) {
////				Set<Integer> set1 = new HashSet();
////				set1.add(entry.getKey()); set1.addAll(entry.getValue());
////				set.add(set1);
////			}
////			Set<Set <Integer>> pureSet = new HashSet();
////			for(Set<Integer> s1:set) {
////				boolean flag = false;
////				for (Set<Integer> s2:set) {
////					if (s1==s2)continue;
////					java.util.Iterator<Integer> it1 = s1.iterator();
////					while(it1.hasNext()) {
////						if(s2.contains(it1.next())) {
////							flag = true;
////						}
////					}
////					if (flag = false) pureSet.add(s1);
////				}
////			}
////			for(Set<Integer> s: set) {
////				System.out.println(s);
////			}
////			System.out.println("size of rawmotifs = " + rawMotifs.size());
////			System.out.println("size of adjList = " + adjList.size());
////			System.out.println("adjList.get(0) = "+ adjList.get(0));
//			for (Set<Integer> s: set) {
////				System.out.println("s = " + s);
//				java.util.Iterator<Integer> iter = s.iterator();
//				List<String> temp = new ArrayList();
//				while(iter.hasNext()) {		
//					temp.addAll(rawMotifs.get(iter.next()));
//				}
//				Set<String> tempSet = new HashSet(temp);
//				List<String> tempList = new ArrayList(tempSet);
//				retVal.add(tempList);
//			}
//			return retVal;
//		}
	
//	public static List<List<String>> purify5 (Map<Double,Set<String>> sortedICPepMap, Map<List<String>, String> motifKSubMap, 
//			Map<Integer, List<String>> idPepListMap, Map<Integer, String> idKSubMap, HashMap<Character, Double> backgroundFrequency){
//			
//			List< List< String > > retVal = new ArrayList();
//			List< List< Integer > > adjList = new ArrayList();
//			for(int l = 0; l < sortedICPepMap.size(); l++) {
//				adjList.add(new ArrayList());
//			}
//			Map <Integer, Set<Integer>> adjMap = new HashMap();
//			Map <Integer, Integer> adjMapNew = new HashMap();
//			Set<Set<String>> mainListSet = new HashSet();
//			for(Entry<Double,Set<String>> entry: sortedICPepMap.entrySet()) {
//				Set<String> s = entry.getValue();
//				if (s.size()==0) continue;
////				List<String> ls = new ArrayList(s);
//				List<String> lsMain = new ArrayList();
//				
//				
//				for(String ss:s) {
//					lsMain.add(removeDash(ss));
//				}
//				Set<String> lsMainSet = new HashSet(lsMain);
//				mainListSet.add(lsMainSet);
//				
//			}
//			List<List<String>> mainList = new ArrayList();
//			for(Set<String> s:mainListSet) {
//				List<String> ltemp = new ArrayList(s);			
//				mainList.add(ltemp);
//			}
//			
//			System.out.println("mainList = " + mainListSet);
//			System.out.println();
//			System.out.println("   size  of  mainList  =  " + mainListSet.size());
//			for (int i = 0; i < mainList.size(); i++) {
//				for(int j = 0; j < mainList.size(); j++) {
//					if (i==j) continue;
//					adjList.get(i).add(i);
//					if (areSame(mainList.get(i), mainList.get(j), motifKSubMap, idPepListMap, idKSubMap, i, j, backgroundFrequency)) {
//						System.out.println(i + " with " + j +" is " + " true");
//						adjList.get(i).add(j);
//						
//					}
//				}
//			}
//			
//			for(int i = 0; i < adjList.size(); i++) {
//				List<Integer> list = adjList.get(i);
//				boolean flag = false;
//				for(Integer j:list) {
//					if(adjMapNew.containsKey(j)) {
//						flag = true;
//						for(Integer k: list) {
//							adjMapNew.put(k, adjMapNew.get(j));
//						}
//						break;
//					}
//					adjMapNew.putIfAbsent(j,i);
//				}
//				
//			}
////			System.out.println("adjMapNew = " + adjMapNew);
////			System.out.println();
////			System.out.println();
////			System.out.println("adjList = " + adjList);
////			for (int k = 0; k < adjList.size(); k++) {
////				adjMap.put(k, adjList.get(k));
////			}
//			Set<Set <Integer>> set = new HashSet();
//			for(Entry<Integer, Integer> entry1: adjMapNew.entrySet()) {
//				Set<Integer> temp = new HashSet();
//				
//				for(Entry<Integer, Integer> entry2: adjMapNew.entrySet()) {
//					if (entry1.getKey()==entry2.getKey()) continue;
//					if(entry1.getValue() == entry2.getValue()) {
//						temp.add(entry1.getKey());
//						temp.add(entry2.getKey());
//					}else {
//						temp.add(entry1.getKey());
//					}
//					
//				}
//				set.add(temp);
//			}
////			for (Entry <Integer, List < Integer > > entry:adjMap.entrySet()) {
////				Set<Integer> set1 = new HashSet();
////				set1.add(entry.getKey()); set1.addAll(entry.getValue());
////				set.add(set1);
////			}
////			Set<Set <Integer>> pureSet = new HashSet();
////			for(Set<Integer> s1:set) {
////				boolean flag = false;
////				for (Set<Integer> s2:set) {
////					if (s1==s2)continue;
////					java.util.Iterator<Integer> it1 = s1.iterator();
////					while(it1.hasNext()) {
////						if(s2.contains(it1.next())) {
////							flag = true;
////						}
////					}
////					if (flag = false) pureSet.add(s1);
////				}
////			}
////			for(Set<Integer> s: set) {
////				System.out.println(s);
////			}
////			System.out.println("size of rawmotifs = " + rawMotifs.size());
////			System.out.println("size of adjList = " + adjList.size());
////			System.out.println("adjList.get(0) = "+ adjList.get(0));
//			for (Set<Integer> s: set) {
////				System.out.println("s = " + s);
//				java.util.Iterator<Integer> iter = s.iterator();
//				List<String> temp = new ArrayList();
//				while(iter.hasNext()) {		
//					temp.addAll(mainList.get(iter.next()));
//				}
//				Set<String> tempSet = new HashSet(temp);
//				List<String> tempList = new ArrayList(tempSet);
//				retVal.add(tempList);
//			}
//			return retVal;
//		}
	
	public static List<List<String>> purify6 (Map<Integer, List<String>> sortedICPepListMapRank, 
			 Map<Double,Set<String>> sortedICPepMap, Map<Integer, String> idKSubMap, 
			 HashMap<Character, Double> backgroundFrequency){
			
			List< List< String > > retVal = new ArrayList();
			List< List< Integer > > adjList = new ArrayList();
			Map <Integer, List <Integer>> adjListMap = new HashMap();
			
			Map<Integer, Double> map = new HashMap();
			int c = 0;
			for (Entry <Double, Set<String>> entry:sortedICPepMap.entrySet() ) {
				map.put(c, entry.getKey());
				c++;
			}
			Map <Integer, Set<Integer>> adjMap = new HashMap();
			Map <Integer, Integer> adjMapNew = new HashMap();
			Set<Set<String>> mainListSet = new HashSet();
			List<List<String>> mainListList = new ArrayList();
			HashMap <Integer, List<String>> newMap = new HashMap();
			HashMap <Integer, List<String>> newMapRank2Index = new HashMap();
			
			for(Entry<Integer,List<String>> entry: sortedICPepListMapRank.entrySet()) {
				List<String> s = entry.getValue();
				Integer u = entry.getKey();
				if (s.size() == 0) continue;
//				List<String> ls = new ArrayList(s);
				List<String> lsMain = new ArrayList();
			
				for(String ss:s) {
					lsMain.add(removeDash(ss));
				}
				newMap.put(u, lsMain);
				newMapRank2Index.put(u-1, lsMain);
				mainListList.add(lsMain);
				
			}
//			 List<List<String>> mainList = new ArrayList();
//			// List<List<String>> mainList = mainListSet.stream().collect(Collectors.toList());
//			for(Set<String> s:mainListSet) {
//				List<String> ltemp = s.stream().collect(Collectors.toList());
//				//List<String> ltemp = new ArrayList(s);			
//				mainList.add(ltemp);
//			}
			
			for (Entry<Integer, List<String>> entry1: newMapRank2Index.entrySet()) {   //28 DEC
				List <Integer> tempList = new ArrayList();                   //28 DEC
				for (Entry<Integer, List<String>> entry2: newMapRank2Index.entrySet()) {//28 DEC
					if (entry1.getKey()==entry2.getKey()) continue;
					if (areSame(entry1.getValue(), entry2.getValue(), idKSubMap, entry1.getKey(), entry2.getKey(), backgroundFrequency)) {
						tempList.add(entry2.getKey());					     //28 DEC	
					}
				}
				adjListMap.put(entry1.getKey(), tempList);                   //28 DEC
			}
			
//			System.out.println("  adjList = "  );
//			for (Entry<Integer, List<Integer>> entry: adjListMap.entrySet()) {
//				System.out.print(idKSubMap.get(entry.getKey()) + "  ");
//				System.out.println(entry);
//			}
//			for (List<Integer> l:adjListMap.values()) {// 28 DEC
//				System.out.println(l);
//			}
			int n_initial = adjListMap.size(); // 28 DEC
			DisjointUnionSets dus = new DisjointUnionSets(n_initial);
			for (Entry<Integer , List<Integer> > entry:adjListMap.entrySet()) {
				for (Integer j : entry.getValue()) {
					dus.union(entry.getKey(), j);
				}
			}
//			for(int i = 0; i < adjListMap.size(); i++) { // 28 DEC
//				for (Integer j : adjList.get(i)) {       // 28 DEC
//					dus.union(i, j);                     // 28 DEC
//				}
//			}
//			System.out.println("    Here dus ");
//			System.out.println("[0] = "+dus.parent[0]);
//			System.out.println("[1] = "+dus.parent[1]);
//			System.out.println("[2] = "+dus.parent[2]);
//			System.out.println("[3] = "+dus.parent[3]);
//			System.out.println("[4] = "+dus.parent[4]);
//			System.out.println("[5] = "+dus.parent[5]);
//			System.out.println("[6] = "+dus.parent[6]);
//			System.out.println("[7] = "+dus.parent[7]);
//			System.out.println("[8] = "+dus.parent[8]);
//			System.out.println();
//			System.out.println("0, " + idKSubMap.get(0));
//			System.out.println("1, " + idKSubMap.get(1));
//			System.out.println("2, " + idKSubMap.get(2));
//			System.out.println("3, " + idKSubMap.get(3));
//			System.out.println("4, " + idKSubMap.get(4));
//			System.out.println("5, " + idKSubMap.get(5));
//			System.out.println("6, " + idKSubMap.get(6));
//			System.out.println("7, " + idKSubMap.get(7));
//			System.out.println("8, " + idKSubMap.get(8));
//			System.out.println(Arrays.toString(dus.parent));
			Set<Integer> set = new HashSet();
			for (int i = 0; i < dus.parent.length; i++) {
				set.add(i);
			}
			List<Integer> list = new ArrayList(set);
			List<List <Integer> > finalList = new ArrayList();
			for (int i = 0; i < list.size(); i++) {
				Set<Integer> s = new HashSet();
				List<Integer> l = new ArrayList();
				for (int j = 0; j < dus.parent.length; j++) {
					if (dus.parent[j] == list.get(i)) {
						s.add(j);
					}
					 l = new ArrayList(s);
				}
				finalList.add(l);
			}
			for(int i = 0; i < finalList.size(); i++) {
				List<String> s = new ArrayList();
				for (int j: finalList.get(i)) {
					if (sortedICPepMap.get(map.get(j)) != null) {
						s.addAll(sortedICPepMap.get(map.get(j)));
					}
					
				}
				retVal.add(s);
			}
			List<List<String>> filtered = retVal.stream().filter(x -> !x.isEmpty()).collect(Collectors.toList());
//			System.out.println("   List of disjoint Sets  ");
//			for(int i = 0; i < finalList.size();i++) {
//				System.out.println(finalList.get(i));
//			}
	//		System.out.println("Check here for null = " + retVal.get(5));
			return filtered;
		}
	
	public static List<String> motifCleaning (List<String> list, HashMap<Character, Double> backgroundFrequency){
//		System.out.println("before cleaning ... here is begining of method motif cleaning");
//		System.out.println(list);
//		System.out.println("After cleaning ... here is begining of method motif cleaning");
		List<String> retVal = new ArrayList(list);
//		System.out.println("##########       Start Cleaning the List . . . . .   ########    ");
//		for(String s:retVal) {
//			System.out.println(s);
//		}
//		System.out.println("#################################################################");
//		System.out.println("#################################################################");
//		System.out.println("#################################################################");
		double[][] pwm = makeProfile(retVal);
		int n = retVal.get(0).length();
		double[][]IC = informationContent(pwm, n, backgroundFrequency);
		double ic = icOfMotif(IC);
		for (String s:list) {
			List<String> temp = new ArrayList(list);
			temp.remove(s);
			if(temp.size()==0) System.out.println("size of temp = 0");
			double[][] pwmTemp = makeProfile(temp);
			int nTemp = temp.get(0).length();
			double[][]ICTemp = informationContent(pwmTemp, nTemp, backgroundFrequency);
			double icTemp = icOfMotif(ICTemp);
			if (icTemp > ic) {
				retVal.remove(s);
			}
		}
		return retVal;
	}
	
	
	
	public static List<String> rawListFromAlignedList (List<String> list){
		List<String> retVal = new ArrayList();
		for(String s:list) {
			StringBuilder sb = new StringBuilder();
			for(int i = 0; i < s.length(); i++) {
				if(s.charAt(i)!='-') {
					sb.append(s.charAt(i));
				}
			}
			retVal.add(sb.toString());
		}
		return retVal;
	}
	public static List<String> seqAlignment(List<String> rawList, KPart u, KPart x, KPart v){
		List<String> retVal = new ArrayList();
		Map<KPart, Integer> kPartMap = new HashMap();
		int len = maxLenofAlignment(rawList);
		len = len - 7;
		for(String peptide:rawList) {
			kPartMap = getKMer(peptide);
			if (kPartMap.containsKey(u) & kPartMap.containsKey(x) & kPartMap.containsKey(v)) {
				if(kPartMap.get(x) == 1 & len == 2) {
					String s = peptide;
					s = preShift(s);
					s = preShift(s);
					retVal.add(s);
				}else if (kPartMap.get(x) == 2 & len == 2) {
					String s = peptide;
					s = preShift(s);
					s = postShift(s);
					retVal.add(s);
				}
				else if(kPartMap.get(x) == 2 & len == 3) {
					String s = peptide;
					s = preShift(s);
					s = postShift(s);
					s = postShift(s);
					retVal.add(s);
				}else if(kPartMap.get(x) == 0 & len == 3) {
					String s = peptide;
					s = preShift(s);
					s = preShift(s);
					s = postShift(s);
					retVal.add(s);
				}else if(kPartMap.get(x) == 1 & len == 3) {
					String s = peptide;
					s = preShift(s);
					s = preShift(s);
					s = postShift(s);
					retVal.add(s);
				}else if(kPartMap.get(x) == 2 & len == 1) {
					String s = peptide;
					s = preShift(s);
					retVal.add(s);
				}else if(kPartMap.get(x) == kPartMap.get(v) || kPartMap.get(u) == kPartMap.get(x)) {
					continue;
				}
				// ??????????
				else{retVal.add(peptide);}
			}else if (kPartMap.containsKey(u) & kPartMap.containsKey(x)) {
				if (kPartMap.get(x) == 3 & kPartMap.get(u) == 1) {
					String s = peptide;
					for (int i = 0; i < len; i++) {
						s = postShift(s);
					}
					retVal.add(s);
				}
				else if(kPartMap.get(x) == 3 & kPartMap.get(u) == 0 & len == 3) {
					String s = peptide;
					s= preShift(s);
					s= preShift(s);
					s= postShift(s);
					retVal.add(s);
				}
				else if (kPartMap.get(x) == 2 & len == 1) {
					String s = peptide;
					s = preShift(s);
					retVal.add(s);
				}else if (kPartMap.get(x) == 2 & len == 2) {
					String s = peptide;
					s = preShift(s);
					s = postShift(s);
					retVal.add(s);
				} else if (kPartMap.get(x) == 2 & len == 3) {
					String s = peptide;
					s = preShift(s);
					s = postShift(s);
					s = postShift(s);
					retVal.add(s);
				} else if (kPartMap.get(x) == 2 & len == 4) {
					// ????
					retVal.add(peptide);
				}else if (kPartMap.get(x) == 1 & len == 1) {
					// ???
					retVal.add(peptide);
				}else if (kPartMap.get(x) == 1 & len == 2) {
					String s = peptide;
					s = preShift(s);
					s = preShift(s);
					retVal.add(s);
				}else if (kPartMap.get(x) == 1 & len == 3) {
					String s = peptide;
					s = preShift(s);
					s = preShift(s);
					s = postShift(s);
					retVal.add(s);
				}
			}else if (kPartMap.containsKey(x) & kPartMap.containsKey(v)){
				if(kPartMap.get(v) == 3 & kPartMap.get(x) == 0) {
					String s = peptide;
					for (int i = 0; i < len; i++) {
						s = preShift(s);
					}
					retVal.add(s);
				}else if(kPartMap.get(v) == 3 & kPartMap.get(x) == 2 & len == 2) {
					String s = peptide;
					s = preShift(s);
					s = postShift(s);
					retVal.add(s);
				}
				else if (kPartMap.get(v) == 3 & kPartMap.get(x) == 2 & len == 3) {
					String s = peptide;
					s = preShift(s);
					s = postShift(s);
					s = postShift(s);
					retVal.add(s);
				}
				else if(kPartMap.get(v) == 2  & len == 3) {
					String s = peptide;
					s = preShift(s);
					s = preShift(s);
					s = postShift(s);
					retVal.add(s);
				}
				else if(kPartMap.get(v) == 2  & len == 2) {
					String s = peptide;
					s = preShift(s);
					s = preShift(s);
					retVal.add(s);
				}else if (kPartMap.get(v) == 2  & len == 1) {
					// ???
					retVal.add(peptide);
				}
				else if(kPartMap.get(v) == 1 & len == 3) {
					String s = peptide;
					s = preShift(s);
					s = preShift(s);
					s = preShift(s);
					retVal.add(s);
				}else if(kPartMap.get(v) == 1 & len == 2) {
					// ??????
					retVal.add(peptide);
				}else if(kPartMap.get(v) == 1 & len == 1) {
					// ?????
					retVal.add(peptide);
				}				
			}
		}		
		return retVal;
	}
	public static int maxLenofAlignment (List<String> list) {
		int retVal;
		List<Integer> values = new ArrayList();
		List<Integer> nums = new ArrayList();
		for (int i = 0; i < list.size(); i++) {
			for (int j = i+1; j < list.size()-1; j++) {
				int score1 = calculateScore(list.get(i),list.get(j));
				int score2 = calculateScore(preShift(list.get(i)), postShift(list.get(j)));
				int score3 = calculateScore(preShift(preShift(list.get(i))), postShift(postShift(list.get(j))));
				int score4 = calculateScore(preShift(preShift(preShift(list.get(i)))), postShift(postShift(postShift(list.get(j)))));
				int score5 = calculateScore(preShift(preShift(preShift(preShift(list.get(i))))), postShift(postShift(postShift(postShift(list.get(j))))));
				
				int score6 = calculateScore(preShift(list.get(j)), postShift(list.get(i)));
				int score7 = calculateScore(preShift(preShift(list.get(j))), postShift(postShift(list.get(i))));
				int score8 = calculateScore(preShift(preShift(preShift(list.get(j)))), postShift(postShift(postShift(list.get(i)))));
				int score9 = calculateScore(preShift(preShift(preShift(preShift(list.get(j))))), postShift(postShift(postShift(postShift(list.get(i))))));
				nums = Arrays.asList(score1, score2, score3, score4, score5, score6, score7, score8, score9);
				Collections.sort(nums, Collections.reverseOrder());
				if (nums.get(0) == score1) {
					values.add(7);
				}else if (nums.get(0) == score2 || nums.get(0) == score6) {
					values.add(8);
				}else if(nums.get(0) == score3 || nums.get(0) == score7) {
					values.add(9);
				}else if(nums.get(0) == score4 || nums.get(0) == score8) {
					values.add(10);
				}else if(nums.get(0) == score5 || nums.get(0) == score9) {
					values.add(11);
				}
				//values.add(nums.get(0));
			}
		}
		Collections.sort(values, Collections.reverseOrder());
		retVal = values.get(0);
		return retVal;
	}
	public static List<String> unionAlignment(List<String> alignedL1, List<String> alignedL2){
		List<String> retVal = new ArrayList();
		int len1 = alignedL1.get(0).length();
		int len2 = alignedL2.get(0).length();
		int diff;
		if(len1 > len2) {
			diff = len1 - len2;
		}else if(len2 > len1) {
			diff = len2 - len1;
		}else {
			diff = 0;
		}
		if (diff == 0) {
			
		}	
		return retVal;
	}
	public static boolean newAreSameMotifs(List<String> l1, List<String> l2) {
		Map<KPart, Integer> kPartMap1 = new HashMap();
		Map<KPart, Integer> kPartMap2 = new HashMap();
		Set<KPart> kPartSet1 = new HashSet();
		Set<KPart> kPartSet2 = new HashSet();
		int common = 0;
		for(String peptide:l1) {
			kPartMap1 = getKMer(peptide);
			kPartSet1.addAll(kPartMap1.keySet());
		}
		for(String peptide:l2) {
			kPartMap2 = getKMer(peptide);
			kPartSet2.addAll(kPartMap2.keySet());
		}
		for(KPart kp:kPartSet1) {
			if(kPartSet2.contains(kp)) {
				common++;
			}
		}
		System.out.println("common = " + common);
		System.out.println("kPartSet1.size() = " + kPartSet1.size());
		System.out.println("kPartSet2.size() = " + kPartSet2.size());
		System.out.println("(common*100/kPartSet1.size())  = " + (common*100/kPartSet1.size()) );
		System.out.println("(common*100/kPartSet2.size())  = " + (common*100/kPartSet2.size()) );
		if((common*100/kPartSet1.size()) >= 50 || (common*100/kPartSet2.size()) >= 50) {
			return true;
		}
		return false;
	}
	public static List<String> alignTest(List<String> list){
		List <String> retVal = new ArrayList();
		List <String> temp = new ArrayList();
		Map <List<String>, Integer> map = new HashMap();
		Map <List<String>, Integer> sortedmap = new HashMap();
		for(int i = 0; i < list.size()-1; i++) {
			for (int j = i+1; j < list.size(); j++) {
				String s1 = list.get(i); String s2 = list.get(j);
				int score1 = calculateScore(s1,s2);
				int score2 = calculateScore(preShift(s1), postShift(s2));
				int score3 = calculateScore(preShift(preShift(s1)), postShift(postShift(s2)));
				int score4 = calculateScore(preShift(preShift(preShift(s1))), postShift(postShift(postShift(s2))));
				int score5 = calculateScore(preShift(preShift(preShift(preShift(s1)))), postShift(postShift(postShift(postShift(s2)))));
				
				int score6 = calculateScore(preShift(s2), postShift(s1));
				int score7 = calculateScore(preShift(preShift(s2)), postShift(postShift(s1)));
				int score8 = calculateScore(preShift(preShift(preShift(s2))), postShift(postShift(postShift(s1))));
				int score9 = calculateScore(preShift(preShift(preShift(preShift(s2)))), postShift(postShift(postShift(postShift(s1)))));
				List<Integer> numList = Arrays.asList(score1, score2, score3, score4, score5, score6, score7, score8, score9);
				Collections.sort(numList, Collections.reverseOrder());
				List<String> l = Arrays.asList(s1,s2);
				map.put(l, numList.get(0));								
			}
		}
		sortedmap = sortByValue(map);
		List<Integer> values = (List<Integer>) sortedmap.values();
		int max = values.get(0);
		for (Entry<List<String>, Integer> entry:map.entrySet()) {
			if (entry.getValue()==max) {
				temp = entry.getKey();
				break;
			}
		}
		
		StringBuilder sb = new StringBuilder();
		int pre = (16-max)/2;
		int post = 16 - (pre + 7);
		for (int i = 0; i < pre; i++) {
			sb.append('-');
		}
		sb.append(temp.get(0));
		for (int i = 0; i < post; i++) {
			sb.append('-');
		}
		retVal.add(sb.toString());
		
		
		return retVal;
	}
	private static Map<List<String>, Integer> sortByValue(Map<List<String>, Integer> unsortMap) {
		//Map<KPart, Integer> sortedMap = new HashMap<KPart, Integer>();
        // 1. Convert Map to List of Map
        List<Map.Entry<List<String>, Integer>> list =
                new LinkedList<Map.Entry<List<String>, Integer>>(unsortMap.entrySet());

        // 2. Sort list with Collections.sort(), provide a custom Comparator
        //    Try switch the o1 o2 position for a different order
        Collections.sort(list, new Comparator<Map.Entry<List<String>, Integer>>() {
            public int compare(Map.Entry<List<String>, Integer> o1,
                               Map.Entry<List<String>, Integer> o2) {
                return (o2.getValue()).compareTo(o1.getValue());
            }
        });
        Collections.reverse(list);
     // 3. Loop the sorted list and put it into a new insertion order Map LinkedHashMap
        Map<List<String>, Integer> sortedMap = new LinkedHashMap<List<String>, Integer>();
        for (Map.Entry<List<String>, Integer> entry : list) {
            sortedMap.put(entry.getKey(), entry.getValue());
        }

        /*
        //classic iterator example
        for (Iterator<Map.Entry<String, Integer>> it = list.iterator(); it.hasNext(); ) {
            Map.Entry<String, Integer> entry = it.next();
            sortedMap.put(entry.getKey(), entry.getValue());
        }*/
        return sortedMap;
	}
	public static List<String> SeqAlgnmnt (List<String> list, KPart kp1, KPart kp2, KPart kp3){
//		System.out.println("Before SeqAlgnmnt Starts ");
//		System.out.println(list);
//		System.out.println("After SeqAlgnmnt Starts ");
		List<String> retVal = new ArrayList();
		String ksubset = getKSubSet(kp1, kp2, kp3);
		for (String s:list) {
			//if (s.length() >= ksubset.length()) {
				List<String> win = window(s, ksubset.length());
				//System.out.println(" win = " + win);
				String str = bestShift(ksubset, win);
				retVal.add(str);
			//}else {
				//List<String> win = window(ksubset, s.length());
				//String str = bestShift(s, win);
				//retVal.add(str);
				//}
			}
		return retVal;
	}
	public static String getKSubSet (KPart kp1, KPart kp2, KPart kp3) {
		String retVal = "";
		String temp1;
		String s1 = kp1.amino;String s2 = kp2.amino; String s3 = kp3.amino;
		int score1 = calculateScore(postShift(s1), preShift(s2));
		int score2 = calculateScore(postShift(postShift(s1)), preShift(preShift(s2)));
		int score3 = calculateScore(postShift(postShift(postShift(s1))), preShift(preShift(preShift(s2))));
		
		List<Integer> list1 = Arrays.asList(score1, score2, score3);
		Collections.sort(list1, Collections.reverseOrder());
		
		if (list1.get(0) == score1) {
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < postShift(s1).length(); i++) {
				if (postShift(s1).charAt(i) == preShift(s2).charAt(i)) {
					sb.append(postShift(s1).charAt(i));
				}else if (postShift(s1).charAt(i) == '-' & preShift(s2).charAt(i) != '-') {
					sb.append(preShift(s2).charAt(i));
				}else if (postShift(s1).charAt(i) != '-' & preShift(s2).charAt(i) == '-') {
					sb.append(postShift(s1).charAt(i));
				}
				else if (postShift(s1).charAt(i) == '-' & preShift(s2).charAt(i) == '-') {
					sb.append('-');
				}
			}
			temp1 = sb.toString();
		}else if (list1.get(0) == score2) {
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < postShift(postShift(s1)).length(); i++) {
				if (postShift(postShift(s1)).charAt(i) == preShift(preShift(s2)).charAt(i)) {
					sb.append(postShift(postShift(s1)).charAt(i));
				}else if (postShift(postShift(s1)).charAt(i) == '-' & preShift(preShift(s2)).charAt(i) != '-') {
					sb.append(preShift(preShift(s2)).charAt(i));
				}else if (postShift(postShift(s1)).charAt(i) != '-' & preShift(preShift(s2)).charAt(i) == '-') {
					sb.append(postShift(postShift(s1)).charAt(i));
				}
				else if (postShift(postShift(s1)).charAt(i) == '-' & preShift(preShift(s2)).charAt(i) == '-') {
					sb.append('-');
				}
			}
			temp1 = sb.toString();
		}else {

			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < postShift(postShift(postShift(s1))).length(); i++) {
				if (postShift(postShift(postShift(s1))).charAt(i) == preShift(preShift(preShift(s2))).charAt(i)) {
					sb.append(postShift(postShift(postShift(s1))).charAt(i));
				}else if (postShift(postShift(postShift(s1))).charAt(i) == '-' & preShift(preShift(preShift(s2))).charAt(i) != '-') {
					sb.append(preShift(preShift(preShift(s2))).charAt(i));
				}else if (postShift(postShift(postShift(s1))).charAt(i) != '-' & preShift(preShift(preShift(s2))).charAt(i) == '-') {
					sb.append(postShift(postShift(postShift(s1))).charAt(i));
				}
				else if (postShift(postShift(postShift(s1))).charAt(i) == '-' & preShift(preShift(preShift(s2))).charAt(i) == '-') {
					sb.append('-');
				}
			}
			temp1 = sb.toString();
		}
		//System.out.println("temp 1 = " + temp1);
		int len = temp1.length();
		String s3_1 = s3; String s3_2 = s3; String s3_3 = s3;
		for(int i = 0; i < len-3; i++) {
			s3_1 = preShift(s3_1);
		}
		for(int i = 0; i < len-2; i++) {
			s3_2 = preShift(s3_2);
		}
		for(int i = 0; i < len-1; i++) {
			s3_3 = preShift(s3_3);
		}
		int score4 = calculateScore(postShift(temp1), s3_1);
		int score5 = calculateScore(postShift(postShift(temp1)),s3_2);
		int score6 = calculateScore(postShift(postShift(postShift(temp1))), s3_3);
		
		List<Integer> list2 = Arrays.asList(score4, score5, score6);
		Collections.sort(list2, Collections.reverseOrder());
		
		if (list2.get(0) == score4) {
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < postShift(temp1).length(); i++) {
				if (postShift(temp1).charAt(i) == s3_1.charAt(i)) {
					sb.append(postShift(temp1).charAt(i));
				}else if (postShift(temp1).charAt(i) == '-' & s3_1.charAt(i) != '-') {
					sb.append(s3_1.charAt(i));
				}else if (postShift(temp1).charAt(i) != '-' & s3_1.charAt(i) == '-') {
					sb.append(postShift(temp1).charAt(i));
				}
				else if (postShift(temp1).charAt(i) == '-' & s3_1.charAt(i) == '-') {
					sb.append('-');
				}
			}
			retVal = sb.toString();
		} else if (list2.get(0) == score5) {
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < postShift(postShift(temp1)).length(); i++) {
				if (postShift(postShift(temp1)).charAt(i) == s3_2.charAt(i)) {
					sb.append(postShift(postShift(temp1)).charAt(i));
				}else if (postShift(postShift(temp1)).charAt(i) == '-' & s3_2.charAt(i) != '-') {
					sb.append(s3_2.charAt(i));
				}else if (postShift(postShift(temp1)).charAt(i) != '-' & s3_2.charAt(i) == '-') {
					sb.append(postShift(postShift(temp1)).charAt(i));
				}
				else if (postShift(postShift(temp1)).charAt(i) == '-' & s3_2.charAt(i) == '-') {
					sb.append('-');
				}
			}
			retVal = sb.toString();
		
		}else {
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < postShift(postShift(postShift(temp1))).length(); i++) {
				if (postShift(postShift(postShift(temp1))).charAt(i) == s3_3.charAt(i)) {
					sb.append(postShift(postShift(postShift(temp1))).charAt(i));
				}else if (postShift(postShift(postShift(temp1))).charAt(i) == '-' & s3_3.charAt(i) != '-') {
					sb.append(s3_3.charAt(i));
				}else if (postShift(postShift(postShift(temp1))).charAt(i) != '-' & s3_3.charAt(i) == '-') {
					sb.append(postShift(postShift(postShift(temp1))).charAt(i));
				}
				else if (postShift(postShift(postShift(temp1))).charAt(i) == '-' & s3_3.charAt(i) == '-') {
					sb.append('-');
				}
			}
			retVal = sb.toString();		
		}
		//System.out.println("End of getKSubSet method ==>  " + retVal);
		return retVal;
	}
	private static Map<String, Integer> sortByValueReversev1(Map<String, Integer> unsortMap) {
		//Map<KPart, Integer> sortedMap = new HashMap<KPart, Integer>();
        // 1. Convert Map to List of Map
        List<Map.Entry<String, Integer>> list =
                new LinkedList<Map.Entry<String, Integer>>(unsortMap.entrySet());

        // 2. Sort list with Collections.sort(), provide a custom Comparator
        //    Try switch the o1 o2 position for a different order
        Collections.sort(list, new Comparator<Map.Entry<String, Integer>>() {
            public int compare(Map.Entry<String, Integer> o1,
                               Map.Entry<String, Integer> o2) {
                return (o1.getValue()).compareTo(o2.getValue());
            }
        });
        Collections.reverse(list);
     // 3. Loop the sorted list and put it into a new insertion order Map LinkedHashMap
        Map<String, Integer> sortedMap = new LinkedHashMap<String, Integer>();
        for (Map.Entry<String, Integer> entry : list) {
            sortedMap.put(entry.getKey(), entry.getValue());
        }

        /*
        //classic iterator example
        for (Iterator<Map.Entry<String, Integer>> it = list.iterator(); it.hasNext(); ) {
            Map.Entry<String, Integer> entry = it.next();
            sortedMap.put(entry.getKey(), entry.getValue());
        }*/
        return sortedMap;
	}
	public static List<String> window(String s, int len){
		
		List<String> retVal = new ArrayList();
		
		int l = len - s.length();
		l = Math.max(4,l);
		for(int i = 0; i <= l; i++) {
			String sCopy = s;
			int k = i;
			while(k > 0) {
				sCopy = preShift(sCopy);
				k--;
			}
			
			for(int j = l-i; j > 0; j--) {																				
					sCopy = postShift(sCopy);
					
			}
			retVal.add(sCopy);
		}
		
		return retVal;
	}
	public static String bestShift (String kSubSet, List<String> windowList){
		 String retVal = "";
		 Map<String, Integer> map = new HashMap();
		 StringBuilder sb = new StringBuilder();
		 if (kSubSet.length() < 10) {
			 sb.append('-');sb.append('-');sb.append(kSubSet);
		 }else if (kSubSet.length() == 10){
			 sb.append('-');sb.append(kSubSet);
		 }else {
			 sb.append(kSubSet);
		 }
//		 System.out.println("windowList " + windowList);
//		 System.out.println("kSubSet " + kSubSet);
		 for (String s: windowList) {			
//			 System.out.println("s = ");
//			 System.out.println(s);
//			 System.out.println("sb.toString() = ");
//			 System.out.println(sb.toString());
			 int score = calculateScore(sb.toString(), s);
			 map.put(s, score);
		 }
		 Map<String, Integer> mapReversed = sortByValueReversev1(map);
		 List<Integer> value = mapReversed.values().stream().collect(Collectors.toList());
		 
		 int max = value.get(0);
		 String temp = "";
		 for(Entry<String, Integer> entry:mapReversed.entrySet()) {
			 if (entry.getValue() == max) {
				 temp = entry.getKey();
				 break;
			 }
		 }			 		 
		 return temp;		 
	 }
	public static String removeDash(String s) {
		
		StringBuilder sb = new StringBuilder();
		for(int i = 0; i < s.length(); i++) {
			if (s.charAt(i) != '-') {
				sb.append(s.charAt(i));
			}
		}
		return sb.toString();
	}
	
	public static String kSubAlignment (String ksub1, String ksub2) {
		String retVal = "";
		List<String> list = new ArrayList();
//		System.out.println("KSUB 1 = " + ksub1);
//		System.out.println("KSUB 2 = " + ksub2);
		if (ksub1.length() > ksub2.length()) {
			int len = ksub1.length() - ksub2.length();
			List<String> winList = window2 (ksub2, len);
//			System.out.println("First branch");
//			System.out.println("winList = ");
//			System.out.println(winList);
			list.add(bestShift2(ksub1, winList));
			list.add(ksub1);
		}else if (ksub2.length() > ksub1.length()) {
			int len = ksub2.length() - ksub1.length();
			List<String> winList = window2 (ksub1, len);
//			System.out.println("Second branch");
//			System.out.println("winList = ");
//			System.out.println(winList);
			list.add(bestShift2(ksub2, winList));
			list.add(ksub2);
		}else {
//			System.out.println("Third branch");
			list = makeRawAlignment(ksub1, ksub2);
		}
		StringBuilder sb = new StringBuilder();
		for(int i = 0; i < list.get(0).length(); i++) {
			if(list.get(0).charAt(i)==list.get(1).charAt(i)) {
				sb.append(list.get(0).charAt(i));
			}else if (list.get(0).charAt(i)=='-' & list.get(1).charAt(i)!='-') {
				sb.append(list.get(1).charAt(i));
			}else if (list.get(1).charAt(i)=='-' & list.get(0).charAt(i)!='-') {
				sb.append(list.get(0).charAt(i));
			}else {
				sb.append('-');
			}
		}
		retVal = sb.toString();
//		System.out.println("M3 = " + retVal);
		return retVal;
	}
	
	public static int kSubAlignmentCalc (String ksub1, String ksub2) {
		String retVal = "";
		List<String> list = new ArrayList();
//		System.out.println("KSUB 1 = " + ksub1);
//		System.out.println("KSUB 2 = " + ksub2);
		if (ksub1.length() > ksub2.length()) {
			int len = ksub1.length() - ksub2.length();
			List<String> winList = window2 (ksub2, len);
//			System.out.println("First branch");
//			System.out.println("winList = ");
//			System.out.println(winList);
			list.add(bestShift2(ksub1, winList));
			list.add(ksub1);
		}else if (ksub2.length() > ksub1.length()) {
			int len = ksub2.length() - ksub1.length();
			List<String> winList = window2 (ksub1, len);
//			System.out.println("Second branch");
//			System.out.println("winList = ");
//			System.out.println(winList);
			list.add(bestShift2(ksub2, winList));
			list.add(ksub2);
		}else {
//			System.out.println("Third branch");
			list = makeRawAlignment(ksub1, ksub2);
		}
		int score = 0;
		for(int i = 0; i < list.get(0).length(); i++) {
			if(list.get(0).charAt(i)==list.get(1).charAt(i)) {
				score++;
			}else if (list.get(0).charAt(i)=='-' & list.get(1).charAt(i)!='-') {
				continue;
			}else if (list.get(1).charAt(i)=='-' & list.get(0).charAt(i)!='-') {
				continue;
			}
		}
		
//		System.out.println("Score of Motif alginments = " + score);
		return score;
	}
	
public static List<String> window2 (String s, int len){
		
		List<String> retVal = new ArrayList();
		
		for(int i = 0; i <= len; i++) {
			String sCopy = s;
			int k = i;
			while(k > 0) {
				sCopy = preShift(sCopy);
				k--;
			}
			
			for(int j = len-i; j > 0; j--) {																				
					sCopy = postShift(sCopy);
					
			}
			retVal.add(sCopy);
		}
		
		return retVal;
	}
public static String bestShift2 (String kSubSet, List<String> windowList){
	 String retVal = "";
	 Map<String, Integer> map = new HashMap();
	 StringBuilder sb = new StringBuilder();
	 for (String s: windowList) {			
//		 System.out.println("s = ");
//		 System.out.println(s);
//		 System.out.println("sb.toString() = ");
//		 System.out.println(sb.toString());
		 int score = calculateScore(kSubSet, s);
		 map.put(s, score);
	 }
	 Map<String, Integer> mapReversed = sortByValueReversev1(map);
	 List<Integer> value = mapReversed.values().stream().collect(Collectors.toList());
	 
	 int max = value.get(0);
	 String temp = "";
	 for(Entry<String, Integer> entry:mapReversed.entrySet()) {
		 if (entry.getValue() == max) {
			 temp = entry.getKey();
			 break;
		 }
	 }			 		 
	 return temp;		 
}
}

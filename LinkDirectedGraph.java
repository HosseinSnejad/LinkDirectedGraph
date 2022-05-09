
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
			  
				  Set<String> alignedSet = new HashSet(bestAlignmentCleaned);
				  List<String> ls = new ArrayList(alignedSet);

				  double[][] pwm = makeProfile(bestAlignment);
				  int n = bestAlignment.get(0).length();
				  double[][]IC = informationContent(pwm, n, backgroundFrequency);	    		
				  double ic = icOfMotif(IC);
				  icList.add(ic);
				  ICPepMap.put(ic, alignedSet);
				  ICPepListMap.put(ic, ls);
				  icKPartNodesMap.putIfAbsent(ic, kpartNodes);

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

		         out.println();
		         double[][] pwm = makeProfile(temp);
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
		    		out.println();
		    		out.println("Information Content Matrix: ");
		    		out.println();
		    		
		    		
		            double ic = icOfMotif(IC);
		            out.println("information content = " + ic);
		            out.println();
		            out.println("================================================================================================================");
		            out.println("================================================================================================================");
		            
		    }
	}
	 
	 public static boolean hasEdge(Map<List<Integer>, Edge> edglist, int i, int j) {
		return edglist.containsKey(Arrays.asList(i, j));
	 }
	
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

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;



// PREPROCESSING AN INSTANCE BEFORE SOLVING THE CORRESPONDING PROBLEM

// Methods of this class answer to questions : Q2, Q3, Q4, Q5
// the principal methods of this class are :
// #) quickPreprocessing() : Q2
// #) removeIPdominated()  : Q3
// #) removeLPdominated()  : Q4
// #) testProcessing()    : Q5

class Preprocessing {
	
	
	// The following method is the preprocessing method we will use after.
	// It returns false if the problem has no solution, and true otherwise
	
	public static boolean preprocess(Instance Inst) {
		if (quickPreprocessing(Inst) == false) {
			System.out.println("Their is no solution for this problem");
			return false;
		}
		removeIPdominated(Inst);
		findLPdominated(Inst);
		return true;
	}
	
	// The following method answers to Q2
	// it returns false if the ILP problem corresponding to the instance Inst has no feasible solution.
	// And it removes triplets (k,m,n) that prevent, when chosen, any solution to be feasible
	
	public static boolean quickPreprocessing(Instance Inst) {
		// Computing the minimal total transmission power
		int minTotalP = 0;
		for (int n=0; n<Inst.N; n++)
			// Instance.Indices[n] is sorted in ascending order of p_kmn
			minTotalP += Inst.p(Inst.Indices[n]);
		// Checking the existence of at least one solution
		if (minTotalP > Inst.P)
			return false;
		// Removing the triplets preventing any solution to be feasible
		for (int n=0; n<Inst.N; n++) {
			DoubleLinked s = Inst.Indices[n].next;
			while (s!=null && Inst.p(s) <= Inst.P - minTotalP + Inst.p(Inst.Indices[n]))
				s = s.next;
			if (s != null)
				s = s.prev.next = null;
			}
		return true;
	}
	
	// The following method answers to Q3
	// It removes all the IP-dominated triplets and returns their number
	
	public static int removeIPdominated(Instance Instance) {
		int count = 0;
		for (int n=0; n<Instance.N; n++) {
			DoubleLinked i = Instance.Indices[n];
			DoubleLinked j = i.next;
			while (j != null) {
				if ( Instance.r(i) >= Instance.r(j) ) {
					j.remove();
					count++;
					j = j.next;
				}
				else {
					if (Instance.p(i) == Instance.p(j)) {
						i.remove();
						if (Instance.Indices[n]==i)
							Instance.Indices[n] = j;
						count++;
					}
					i = i.next;
					j = j.next;
				}
			}
		}
		return count;
	}
	
	// The two following methods answer to Q4.	
	// The following method finds all the LP-dominated triplets and returns their number
	// It sets their attribute lpDominated to true
	
	public static int findLPdominated(Instance Inst) {
		int count = 0;
		for (int n=0; n<Inst.N; n++) {
			DoubleLinked i = Inst.Indices[n];
			DoubleLinked j = i.next;
			if (j == null )
				continue;
			DoubleLinked k = j.next;
			while (k != null) {
				if ( (Inst.r(k)-Inst.r(j))*(Inst.p(j)-Inst.p(i)) >= (Inst.r(j)-Inst.r(i))*(Inst.p(k)-Inst.p(j))) {
					j.lpDominated = true;
					count++;
					if (i == Inst.Indices[n]) {
						j = k;
						k = k.nextLP();
					}
					else {
						j = i;
						i = i.prevLP();
					}
				}
				else {
					i = j;
					j = k;
					k = k.nextLP();
				}
			}
		}
		return count;
	}
	
	// The following method returns a copy of Instance.Indices from which all the LP-dominated triplets were removed
	// It does not modify the instance Instance
	// It assumes that we have already applied the method findLPdminated to Instance
	
	public static DoubleLinked[] removeLPdominated(Instance Inst) {
		DoubleLinked[] Indices = Inst.Indices;
		DoubleLinked[] IndicesLP = new DoubleLinked[Indices.length];
		for (int n=0; n<Indices.length; n++) {
			IndicesLP[n] = new DoubleLinked(Indices[n].k, Indices[n].m, Indices[n].n);
			DoubleLinked d = Indices[n].next;;
			DoubleLinked dLP = IndicesLP[n];
			while (d!=null) {
				if (d.lpDominated == false) {
					dLP.insertAfter(d.k, d.m, d.n);
					dLP = dLP.next;
				}
				d = d.next;
			}
		}
		return IndicesLP;
	}
	
	
	
	// The following method answers to Q5.
	// It shows how the size of the instance changes after applying each of the previous methods
	
	public static void testPreprocessing(String path) throws IOException {
		Instance Inst = new Instance(path);
		int numberTriplets = Inst.N*Inst.K*Inst.M;
		System.out.println("Initial number of triplets : " + numberTriplets);
		
		// Test quickPreprocessing
		System.out.println("After applying the method quickPreprocessing :");
		if (quickPreprocessing(Inst) == false) {
			System.out.println("The IP problem for this instance has no solution");
			return;
		}
		else {
			int count = DoubleLinked.size(Inst.Indices);
			System.out.println("Number of removed triplets   : " + (numberTriplets - count) );
			System.out.println("Number of remaining triplets : " + count);
			numberTriplets = count;
		}
		
		// Test removeIPdominated
		int countIP = removeIPdominated(Inst);
		numberTriplets -= countIP;
		System.out.println("After applying the method removeIPdominated");
		System.out.println("Number of removed triplets   : " + countIP );
		System.out.println("Number of remaining triplets : " + numberTriplets);
		
		// Test removeLPdominated
		int count = findLPdominated(Inst);
		DoubleLinked[] IndicesLP = removeLPdominated(Inst);
		int countRemainingLP = DoubleLinked.size(IndicesLP);
		System.out.println("After applying the method removeLPdominated");
		System.out.println("Number of removed triplets   : " + (numberTriplets - countRemainingLP) );
		System.out.println("Number of remaining triplets : " + (numberTriplets-count));
	}

	
	
}





//GREEDY ALGORITHM FOR THE LP PROBLEM

//Methods of this class answer to Q6 and Q7


class LinearP {

	// The following method answers to Q6
	
	public static Solution solve(Instance Inst) {
		Solution sol = new Solution(Inst);
		solveSubLP(sol, constructSorted(Inst), 0, Inst.P);
		return sol;
	}
	
	
	public static boolean solveSubLP(Solution sol, DoubleLinked sortedE, int n0, int leftBudgetPower) {
		Instance Inst = sol.Inst;
		int usedPower = 0;
		for (int n=n0; n<Inst.N; n++) {
			sol.x[n] = new Couple(Inst.Indices[n].k, Inst.Indices[n].m);
			usedPower += Inst.p(Inst.Indices[n]);
		}
		if (usedPower > leftBudgetPower)
			return false;
		DoubleLinked d = sortedE;
		while (d!= null && d.n < n0)
			d = d.next;
		while (d!=null && (usedPower < leftBudgetPower) ) {
			int n = d.n;
			if (usedPower + Inst.p(d) - Inst.p(sol.x[n].k, sol.x[n].m, n) <= leftBudgetPower) {
				usedPower = usedPower + Inst.p(d) - Inst.p(sol.x[n].k, sol.x[n].m, n);
				sol.x[d.n] = new Couple(d.k, d.m);
			} else {
				int p_kmn = Inst.p(d);
				float lambda = (float)(leftBudgetPower - usedPower)/(p_kmn);
				sol.solpart(d.k,d.m,n,lambda);
				usedPower += (int) (lambda*p_kmn);
			}
			do {
			d = d.next;
			} while (d!=null && d.n < n0);
		}
		return true;
	}
	
	
	public static DoubleLinked constructSorted(Instance Inst) {
		DoubleLinked[] IndicesLP = Preprocessing.removeLPdominated(Inst);
		int size = DoubleLinked.size(IndicesLP) - IndicesLP.length;
		DoubleLinked[] a = new DoubleLinked[size];
		int i = 0;
		for (int n=0; n<IndicesLP.length; n++) {
			DoubleLinked d = IndicesLP[n].next;
			while (d!=null) {
				a[i++] = d;
				d = d.next;
			}
		}
		return sortDecreasing(Inst, a);
	}
	
	public static DoubleLinked sortDecreasing(Instance Inst, DoubleLinked[] a) {
		mergeSortRec(Inst, a, new DoubleLinked[a.length], 0, a.length);
		a[0].prev = null;
		for (int i=0; i<a.length-1; i++) {
			a[i].next = a[i+1];
			a[i+1].prev = a[i];
		}
		a[a.length - 1].next = null;
		return a[0];
	}
	
	// The following methods are for merge-sort descending
	
	public static void mergeSortRec(Instance Inst, DoubleLinked[] a, DoubleLinked[] tmp, int left, int right) {
		if (left >= right - 1)
			return;
		int med = left +(right - left)/2;
		mergeSortRec(Inst, a, tmp, med, right);
		mergeSortRec(Inst, a, tmp, left, med);
		for (int i=left; i< right; i++)
			tmp[i] = a[i];
		mergeDescending(Inst, tmp, a, left, med, right);
	}
	public static void mergeDescending(Instance Inst, DoubleLinked[] a1, DoubleLinked[] a2, int left, int med, int right) {
		int i = left, j = med;
		for (int s=left; s<right; s++) {
			if (i<med && (j == right || e(Inst, a1[i]) > e(Inst, a1[j]) ))
				a2[s] = a1[i++];
			else
				a2[s] = a1[j++];
		}
	}
	
	// The following method defines the efficiency e.
	
	public static double e(Instance Inst, DoubleLinked d) {
		return (double)(Inst.r(d) - Inst.r(d.prev))/(Inst.p(d) - Inst.p(d.prev));
	}
	
	// The following method answers to Q7
	
	public static void testSolve(String path) throws IOException {
		Instance Inst = new Instance(path);
		Solution sol = null;
		double ti=0,tf=0,tm=0;
		for(int i = 0; i < 1000; i++) {
			ti = System.currentTimeMillis();
			if (Preprocessing.preprocess(Inst) == false)
				return;
			sol = solve(Inst);
			tf = System.currentTimeMillis();
			tm += tf - ti;
		}
		System.out.println(sol);
		System.out.println("Run time : " + (tf-ti) + " ms");
		System.out.println("Average run time : " + (tm/1000) + " ms");
	}
	
}



class Solution {
	Instance Inst;
	Couple[] x;
	int k, m, n;
	float lambda;
	
	
	public Solution(Instance Inst) {
		this.Inst = Inst;
		this.x = new Couple[Inst.N];
		this.lambda = 0;
	}
	
	public Solution(Instance Inst, Couple[] x, int k, int m, int n, float lambda) {
		this.Inst = Inst;
		this.x = x;
		this.k = k;
		this.m = m;
		this.n = n;
		this.lambda = lambda;
	}
	
	public void solpart(int k , int m, int n, float lambda) {
		this.k = k;
		this.m = m;
		this.n = n;
		this.lambda = lambda;
	}
	
	public Solution copy() {
		Couple[] xCopy = new Couple[x.length];
		for (int i=0; i<x.length; i++)
			xCopy[i] = x[i];
		return new Solution(Inst, xCopy, k, m, n, lambda);
	}
	
	// This method returns true if the solution found is a solution for the ILP problem
	public boolean isSolutionILP() {
		return ( lambda == 0 );
	}
	
	// This method returns the used power for a solution of the ILP or LP problem
	public float usedPower() {
		float usedP = 0;
		for (int i=0; i<x.length; i++)
			usedP += Inst.p(x[i].k, x[i].m, i);
		if (lambda != 0)
			usedP += lambda*(Inst.p(k,m,n) - Inst.p(x[n].k, x[n].m, n));
		return usedP;
	}
	
	// This method returns the total rate for a solution of the ILP or LP problem
	public float rate() {
		float dataR = 0;
		for (int i=0; i<x.length; i++)
			dataR += Inst.r(x[i].k, x[i].m, i);
		if (lambda != 0)
			dataR += lambda*(Inst.r(k,m,n) - Inst.r(x[n].k, x[n].m, n));
		return dataR;
	}
	
	// This method returns a string corresponding to the solution found
	// it shows the chosen triplets, the initial budget power, the total used power, and the total data rate
	
	public String toString() {
		String s = "";
		s += "\nThis is a solution for the ILP problem !";
		s += "\nThe optimal data rate    : "+(int)rate();
		s += "\nPower buget              : " + Inst.P;
		s += "\nThe total used power     : "+(int)usedPower();
		return s;
	}
	
	

	
}


//The principal method of this class is the constructor 
//This class transforms a text file to an object we can manipulate in our program.
//This object represents an instance of the LP or IP problem, for which we will try to find a solution.



class Instance {
	int N,K,M,P;	//  Number of channels, users, different powers for a (user,channel) and The total budget.(respectively)
	int[][][] p,r;  // matrixes of power and data rate respectively. (Dimensions are M,K,N) 
	DoubleLinked[] Indices;
	
	
	public Instance(int[][][] p, int[][][] r, int P) {
		this.p = p;
		this.r = r;
		this.K = p.length;
		this.M = p[0].length;
		this.N = p[0][0].length;
		this.P = P;
		constructIndices();
	}
	
	public Instance(String path) throws IOException{
		BufferedReader bf = new BufferedReader(new FileReader(path));
		// Initializing N, K, M and P
		this.N = (int)Float.parseFloat(bf.readLine());
		this.M = (int)Float.parseFloat(bf.readLine());
		this.K = (int)Float.parseFloat(bf.readLine());
		this.P = (int)Float.parseFloat(bf.readLine());
		// constructing p and r
		this.p = constructMatrix_kmn(bf);
		this.r = constructMatrix_kmn(bf);
		bf.close();
		// constructing Indices
		constructIndices();
	}
	
	// Construction of p and r
	
	public int[][][] constructMatrix_kmn(BufferedReader bf) throws IOException {
		int[][][] matrix = new int[K][M][N];
		for(int n=0; n<N; n++) {
			for(int k=0; k<K; k++) {
				String[] line = bf.readLine().trim().split("\\s+");
				for (int m=0; m<M; m++) 
					matrix[k][m][n] = (int)Float.parseFloat(line[m]);
			}
		}
		return matrix;
	}
	
	public int p(int k, int m, int n) {
		return this.p[k][m][n];
	}
	
	public int r(int k, int m, int n) {
		return this.r[k][m][n];
	}
	
	public int p(DoubleLinked triplet) {
		return this.p[triplet.k][triplet.m][triplet.n];
	}
	public int r(DoubleLinked triplet) {
		return this.r[triplet.k][triplet.m][triplet.n];
	}
	
	
	public void constructIndices() {
		Indices = new DoubleLinked[N];
		for (int n=0; n<N; n++) {
			DoubleLinked[] a = new DoubleLinked[K*M];
			int i =0;
			for (int k=0; k<K; k++)
				for (int m=0; m<M; m++)
					a[i++] = new DoubleLinked(k,m,n);
			Indices[n] = sortAscendingP(a);
		}
	}
	
	// We use the merge-sort algorithm to sort triplets (k,m,n) in ascending order of p_kmn
	// the following method assumes that initially a[i].prev == a[i].next == null for each i
	
	public DoubleLinked sortAscendingP(DoubleLinked[] a) {
		mergeSortRec(a, new DoubleLinked[a.length], 0, a.length);
		// Now that a is sorted, we will link its elements
		for (int i=0; i<a.length-1; i++) {
			a[i].next = a[i+1];
			a[i+1].prev = a[i];
		}
		return a[0];
	}
	public void mergeSortRec(DoubleLinked[] a, DoubleLinked[] tmp, int left, int right) {
		if (left >= right - 1)
			return;
		int med = left +(right - left)/2;
		mergeSortRec(a, tmp, med, right);
		mergeSortRec(a, tmp, left, med);
		for (int i=left; i< right; i++)
			tmp[i] = a[i];
		mergeAscendingP(tmp, a, left, med, right);
	}
	public void mergeAscendingP(DoubleLinked[] a1, DoubleLinked[] a2, int left, int med, int right) {
		int i = left, j = med;
		for (int s=left; s<right; s++) {
			if (i<med && (j == right || p(a1[i]) < p(a1[j]) ))
				a2[s] = a1[i++];
			else
				a2[s] = a1[j++];
		}
	}
	



}

//Objects of this class are Double linked lists of triplets.

class DoubleLinked {
	int k,m,n;
	boolean lpDominated; // this one is for marking the LP-Dominated triples
	DoubleLinked next,prev;
	
	public DoubleLinked(int k, int m, int n) {
		this.k = k;
		this.m = m;
		this.n = n;
		this.lpDominated = false;
		this.next = this.prev = null;
	}
	
	public DoubleLinked(int[][] L) {
		this.k = L[0][0];
		this.m = L[0][1];
		this.n = L[0][2];
		for(int i=L.length-1; i>0; i--) {
			this.insertAfter(L[i][0], L[i][1], L[i][2]);
		}
	}
	
	public void insertAfter(int k, int m, int n) {
		DoubleLinked d = new DoubleLinked(k,m,n);
		d.next = this.next;
		d.prev = this;
		this.next = d;
		if (d.next!=null) d.next.prev = d;
	}
	
	public void insertBefore(int k, int m, int n) {
		DoubleLinked d = new DoubleLinked(k,m,n);
		d.prev = this.prev;
		d.next = this;
		this.prev = d;
		if (d.prev!=null) d.prev.next = d;
	}
	
	public void remove() {
		if (this.prev != null) this.prev.next = this.next;
		if (this.next != null) this.next.prev = this.prev;
	}
	
	public int size() {
		int size = 0;
		DoubleLinked s = this;
		while (s != null) {
			s = s.next;
			size++;
		}
		return size;
	}
	
	public static int size(DoubleLinked[] L) {
		int s = 0;
		for(int n=0; n<L.length; n++)
			s += L[n].size();
		return s;
	}
	
	public String toString() {
		String s = "";
		System.out.print("[HEAD] --> ");
		DoubleLinked d = this;
		while (d != null) {
			s += "("+d.k+","+d.m+","+d.n + ") --> ";
			d = d.next;
		}
		s += "[TAIL] .";
		return s;
	}

	
	public DoubleLinked nextLP() {
		DoubleLinked s = this.next;
		while( s != null && s.lpDominated == true)
			s = s.next;
		return s;
	}
	
	public DoubleLinked prevLP() {
		DoubleLinked s = this.prev;
		while(s != null && s.lpDominated == true)
			s = s.prev;
		return s;
	}
	
}


class Couple {
	public int k,m;
	
	public Couple(int k, int m) {
		this.k = k;
		this.m = m;
	}
}

//DYNAMIC PROGRAMMING ALGORITHM FOR IP PROBLEM

//Methods of this class answer to Q8, Q9 
//The principal methods of this class are :
//	- solve()		    : Q8
//  - solveQ9()       	: Q9

class DynamicP {
	

	// The following method answers to Q8,

	
	public static Solution solve(Instance ins) {
		Solution sol = new Solution(ins);
		int[][] Debit = computeDebit(ins);
		int p = ins.P;
		for (int n=0; n<ins.N; n++) {
			DoubleLinked d = ins.Indices[n];
			while (d != null) {
				if (Debit[n][p] == ins.r(d) + Debit[n+1][p-ins.p(d)]) {
					sol.x[n] = new Couple(d.k, d.m);
					p = p-ins.p(d);
					d = null;
				} else
					d = d.next;
			}
		}
		return sol;
	}

	// The following method computes the table of maximum data rate for each sub-problem 
	
	public static int[][] computeDebit(Instance ins) {
		// Initializing Debit
		int[][] Debit = new int[ins.N+1][ins.P+1];
		// Debit[ins.N][p] = Debit[n][0] = 0 for all values of n and p
		for (int n=ins.N-1; n>=0; n--) {
			for (int p=0; p<ins.P+1; p++) {
				DoubleLinked d = ins.Indices[n];
				while (d != null && ins.p(d) <= p ) {
					if (Debit[n][p]< ins.r(d) + Debit[n+1][p-ins.p(d)])
						Debit[n][p] = ins.r(d) + Debit[n+1][p-ins.p(d)];
					d = d.next;
				}
			}
		}
		return Debit;
	}
	
	// The following method answers to Q9,
	
	public static Solution solveQ9(Instance ins) {
		int U = (int)DynamicP.solve(ins).rate();
		Solution sol = new Solution(ins);
		int[][] P = computeP(ins, U);
		int r = U;
		while (P[0][r]==-1)
			r--;
		System.out.println("Max value == " +r);
		for (int n=0; n<ins.N; n++) {
			DoubleLinked d = ins.Indices[n];
			while (d != null) {
				if (ins.r(d)<=r && P[n+1][r-ins.r(d)]!=-1 && P[n][r] == P[n+1][r-ins.r(d)] + ins.p(d)) {
					sol.x[n] = new Couple(d.k, d.m);
					r = r-ins.r(d);
					d = null;
				} else
					d = d.next;
			}
		}
		return sol;
	}

	// The following function computes the table of minimum possible power achieving a data rate r<=U
	
	public static int[][] computeP(Instance ins, int U){
		int[][] Z = new int[ins.N+1][U+1];
		for (int r=0; r<U+1; r++) {
			// When Z[n][r]==-1 then no solution is possible
			Z[ins.N-1][r] = Integer.MAX_VALUE;
			DoubleLinked d = ins.Indices[ins.N-1];
			while (d!= null && ins.r(d)<=r) {
				if (r==ins.r(d) && Z[ins.N-1][r] > ins.p(d))
					Z[ins.N-1][r] = ins.p(d);
				d = d.next;
			}
			if (Z[ins.N-1][r] == Integer.MAX_VALUE)
				Z[ins.N-1][r] = -1;
		}
		// We construct the rest of Z
		for (int n=ins.N-2; n>=0; n--) {
			for (int r=0; r<U+1; r++) {
				Z[n][r] = -1;
				int min = Integer.MAX_VALUE;
				DoubleLinked d = ins.Indices[n];
				while (d!=null) {
					if (ins.r(d)>r || Z[n+1][r-ins.r(d)]==-1) {
						d = d.next;
						continue;
					}
					if (min > Z[n+1][r-ins.r(d)] + ins.p(d))
						min = Z[n+1][r-ins.r(d)] + ins.p(d);
					d = d.next;
				}
				if (min != Integer.MAX_VALUE)
					Z[n][r] = min;
			}
		}
		return Z;
	}
	
	
	
	// This method tests the first DP approach
	public static void testSolve(String path) throws IOException {
		Instance ins = new Instance(path);
		if (Preprocessing.preprocess(ins) == false)
			return;
		double ti = System.currentTimeMillis();
		Solution sol = solve(ins);
		double tf = System.currentTimeMillis();
		System.out.println(sol);
		System.out.println("Run time : " + (tf-ti) + " ms");
	}
	// This method tests the second DP approach
	public static void testSolveQ9(String path) throws IOException {
		Instance ins = new Instance(path);
		if (Preprocessing.preprocess(ins) == false)
			return;
		double ti = System.currentTimeMillis();
		Solution sol = solveQ9(ins);
		double tf = System.currentTimeMillis();	
		System.out.println(sol);
		System.out.println("Run time : " + (tf-ti) + " ms");
	}

	
}


// Branch-and-Bound Algorithm -----------------
class BB {
	Solution bestSolution;
	int bestRate;
	DoubleLinked sortedE;
	
	public BB(Instance Inst){
		this.sortedE = LinearP.constructSorted(Inst);
		this.bestSolution = new Solution(Inst);
		LinearP.solveSubLP(bestSolution, sortedE, 0, Inst.P);
		bestSolution.solpart(0, 0, 0 , 0);
		this.bestRate = (int)bestSolution.rate();
	}
	
	public static Solution solve(Instance Inst) {
		BB bb = new BB(Inst);
		bb.solveRecBB(new Solution(Inst), 0, Inst.P);
		return bb.bestSolution;
	}
	 
	public boolean solveRecBB(Solution sol, int n0, int leftBudgetPower) {
		bestSolution.k++;
		if (n0 == sol.Inst.N) {
			int rate = (int)sol.rate();
			if (rate > bestRate) {
				bestRate = (int)rate;
				bestSolution = sol.copy();
			}
			return true;
		}
		if (!LinearP.solveSubLP(sol, sortedE, n0, leftBudgetPower))
			return false;
		float rate = sol.rate();
		if (rate < bestRate + 1)
			return false;
		if (sol.isSolutionILP()) {
			bestRate = (int)rate;
			bestSolution = sol.copy();
			return true;
		}
		// Recursion
		DoubleLinked d = sol.Inst.Indices[n0];
		while ( d!= null ) {
			sol.solpart(0,0,0,0);
			sol.x[n0] = new Couple(d.k, d.m);
			solveRecBB(sol, n0+1, leftBudgetPower-sol.Inst.p(d));
			d = d.next;
		}
		return true;
	}
	public static void testSolve(String path) throws IOException {
		Instance Inst = new Instance(path);
		double ti = System.currentTimeMillis();
		if (Preprocessing.preprocess(Inst) == false)
			return;
		Solution sol = solve(Inst);
		double tf = System.currentTimeMillis();
		System.out.println(sol);
		System.out.println("Run time : " + (tf-ti) + " ms");
		System.out.println("Number of explored nodes : "+sol.k);
	}
	
	
}
//ONLINE PROGRAMMING

//Methods of this class answer to Q12 and Q13
//the principal methods of this class are :
//- solve()	    : Q12
//- testSolve()   : Q13

//(usedChannels[n]==0) <==> (n is still available)
//eAverage is the expected slope we will use in our efficiencyFunction function
class OnlineP {
	Instance Inst;
	int[] usedChannels;
	int pMax, rMax, budgetPower;
	float eAverage;
	
	// 
	public OnlineP(int K, int M, int N, int pMax, int rMax, int budgetPower) {
		this.Inst = randomInstance(K, M, N, pMax, rMax, budgetPower);
		this.usedChannels = new int[Inst.N];
		this.pMax = pMax;
		this.rMax = rMax;
		this.budgetPower = budgetPower;
		this.eAverage = computeEAverage(pMax, rMax);
	}
	
	
	// This method generates a random instance of IP problem
	public static Instance randomInstance(int K, int M, int N, int pMax, int rMax, int budgetPower) {
		int[][][] p = randomMatrix(K, M, N, pMax);
		int[][][] r = randomMatrix(K, M, N, rMax);
		Instance Inst = new Instance(p, r, budgetPower);
		return Inst;
	}
	
	public static int[][][] randomMatrix(int K, int M, int N, int vMax){
		int[][][] matrix = new int[K][M][N];
		for (int k=0; k<K; k++)
			for (int m=0; m<M; m++)
				for (int n=0; n<N; n++)
					// random value between 1 and vMax both included
					matrix[k][m][n] = (int)(vMax*Math.random())+1;
		return matrix;
	}
	
	public static float computeEAverage(int pMax, int rMax) {
		float e = 0;
		for (int p=1; p<=pMax; p++)
			for (int r=1; r<=rMax; r++)
				e = e + (float)r/p;
		return e/(pMax*rMax);
	}
	
	// The following method answers to Q12
	
	public static Solution solve(int K, int M, int N, int pMax, int rMax, int budgetPower) {
		OnlineP op = new OnlineP(K, M, N, pMax, rMax, budgetPower);
		Solution sol = new Solution(op.Inst);
		sol.k = 1;
		if (op.completeSolution(sol)==false) {
			sol.k = 0;
		}
		return sol;
	}
	
	public boolean completeSolution(Solution sol) {
		int k = 0;
		int UCk = Inst.N;
		while (UCk>0 && k<Inst.K && budgetPower>0) {
			for (int n=0; n<Inst.N; n++) {
				if (usedChannels[n]==-1)
					continue;
				int mChoice = 0;
				float efficiencyChoice = efficiencyFunction(Inst.p(k,0,n), Inst.r(k,0,n));
				for (int m=1; m<Inst.M; m++) {
					float qual = efficiencyFunction(Inst.p(k,m,n), Inst.r(k,m,n));
					if ( qual > 0) {
						mChoice = m;
						efficiencyChoice = qual;
					}
				}
				if (Inst.p(k,mChoice,n)<=budgetPower && BetterThanAverage(efficiencyChoice, k, UCk)) {
					sol.x[n] = new Couple(k,mChoice);
					usedChannels[n] = -1;
					budgetPower -= Inst.p(k,mChoice,n);
				}
			}
			k ++;
		}
		boolean solIsComplete = true;
		for (int n=0; n< Inst.N; n++) {
			if (usedChannels[n] == 0) {
				sol.x[n] = new Couple(-1,-1);
				solIsComplete = false;
			}		
		}
		return solIsComplete;
	}
	public float efficiencyFunction(int p, int r) {
		if (p <= budgetPower)
			return r-p*eAverage;
		else
			return -pMax*rMax;
	}
	
	public boolean BetterThanAverage(float efficiencyChoice, int k, int UCk) {
		if (k==Inst.K-1)
			return true;
		int meanEfficiency = 0;
		for (int p=1; p<=budgetPower; p++)
			for (int r=1; r<=rMax; r++)
				meanEfficiency += efficiencyFunction(p,r);
		meanEfficiency /= pMax*rMax;
		if(efficiencyChoice < meanEfficiency)
			return false;
		return true;

	}
	

	public static float[] testSolve(int K, int M, int N, int pMax, int rMax, int budgetPower) {
		Solution sol = solve(K, M, N, pMax, rMax, budgetPower);
		if (sol.k == 0)
			return new float[] {0,0,0,0};
		Solution optimalSol = DynamicP.solve(sol.Inst);
		return new float[] {sol.rate()/optimalSol.rate(),sol.usedPower()/optimalSol.usedPower(),sol.usedPower(),optimalSol.usedPower()};
	}
	
public static void OnlineSolve(int N) {
		float rateRatio = 0;
		float powerRatio = 0;
		float usedPowerOnline = 0;
		float usedPowerOffline = 0;
		for (int i=0; i<N; i++) {
			float[] r = testSolve(10, 2, 4, 50, 100, 100);
			rateRatio += r[0];
			powerRatio += r[1];
			usedPowerOnline += r[2]; 
			usedPowerOffline += r[3];
		}
		usedPowerOnline /= N;
		usedPowerOffline /= N;
		rateRatio /= N;
		powerRatio /= N;
		System.out.println("TEST ONLINE PROGRAM OVER "+ N +" EXPERIENCES");
		System.out.println("The average power ratio : " + powerRatio);
		System.out.println("The average rate ratio : " + rateRatio);
		System.out.println("The average used power of the Online algorithm : " +usedPowerOnline);
		System.out.println("The average used power of the Offline algorithm : " +usedPowerOffline);
	}
	
}




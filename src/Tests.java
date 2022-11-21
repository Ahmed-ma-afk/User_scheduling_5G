import java.io.IOException;

public class Tests {

	public static void main(String[] args) throws IOException {		
		System.out.println("\n\n[--------------------------------USER SCHEDULING IN 5G---------------------------------------------------------]\n");

		System.out.println("\n\n[--------------------------------Solution of Questions 2, 3, 4 et 5 ---------------------------------------------------------]\n");
		System.out.println("Test file : test1.txt");
		Preprocessing.testPreprocessing("testfiles/test1.txt");
		System.out.println("\nTest file : test2.txt");
		Preprocessing.testPreprocessing("testfiles/test2.txt");
		System.out.println("\nTest file : test3.txt");
		Preprocessing.testPreprocessing("testfiles/test3.txt");
		System.out.println("\nTest file : test4.txt");
		Preprocessing.testPreprocessing("testfiles/test4.txt");
		System.out.println("\nTest file : test5.txt");
		Preprocessing.testPreprocessing("testfiles/test5.txt");
		System.out.println("\n[-----------------------------------------------------------------------------------------------------------------------------]\n\n");
		
		
		System.out.println("\n\n[------------------------------Solution of Questions 6 et 7 ----------------------------------------------------------------]\n");
		System.out.println("Test file : test1.txt");
		LinearP.testSolve("testfiles/test1.txt");
		System.out.println("\nTest file : test2.txt");
		LinearP.testSolve("testfiles/test2.txt");
		System.out.println("\nTest file : test3.txt");
		LinearP.testSolve("testfiles/test3.txt");
		System.out.println("\nTest file : test4.txt");
		LinearP.testSolve("testfiles/test4.txt");
		System.out.println("\nTest file : test5.txt");
		LinearP.testSolve("testfiles/test5.txt");
		System.out.println("\n[-----------------------------------------------------------------------------------------------------------------------------]\n\n");

		
		System.out.println("\n\n[----------------------------------Solution of Questions 8 et 9 ------------------------------------------------------------]\n");
		System.out.println("Test file : 'test1.txt' with First DP approach Q8");
		// First DP approach Q8
		DynamicP.testSolve("testfiles/test1.txt");
		System.out.println("\nTest file : 'test1.txt' with Second DP approach Q9");
		// Second DP approach Q9
		DynamicP.testSolveQ9("testfiles/test1.txt");
		
		System.out.println("\nTest file : 'test2.txt' with First DP approach Q8");
		// First DP approach Q8
		DynamicP.testSolve("testfiles/test2.txt");
		System.out.println("\nTest file : 'test2.txt' with Second DP approach Q9");
		// Second DP approach Q9
		DynamicP.testSolveQ9("testfiles/test2.txt");
		
		System.out.println("\nTest file : 'test3.txt' with First DP approach Q8");
		// First DP approach Q8
		DynamicP.testSolve("testfiles/test3.txt");
		System.out.println("\nTest file : 'test3.txt' with Second DP approach Q9");
		// Second DP approach Q9
		DynamicP.testSolveQ9("testfiles/test3.txt");
		
		System.out.println("\nTest file : 'test4.txt' with First DP approach Q8");
		// First DP approach Q8
		DynamicP.testSolve("testfiles/test4.txt");
		System.out.println("\nTest file : 'test4.txt' with Second DP approach Q9");
		// Second DP approach Q9
		DynamicP.testSolveQ9("testfiles/test4.txt");
		
		System.out.println("\nTest file : 'test5.txt' with First DP approach Q8");
		// First DP approach Q8
		DynamicP.testSolve("testfiles/test5.txt");
		System.out.println("\nTest file : 'test5.txt' with Second DP approach Q9");
		// Second DP approach Q9
		DynamicP.testSolveQ9("testfiles/test5.txt");
		System.out.println("\n[------------------------------------------------------------------------------------------------------------------------------]\n\n");

		
		System.out.println("\n\n[----------------------------------Solution of Questions 10 et 11 -----------------------------------------------------------]\n");
		System.out.println("Test file : test1.txt");
		BB.testSolve("testfiles/test1.txt");
		System.out.println("\nTest file : test2.txt");
		BB.testSolve("testfiles/test2.txt");
		System.out.println("\nTest file : test3.txt");
		BB.testSolve("testfiles/test3.txt");
		System.out.println("\nTest file : test4.txt");
		BB.testSolve("testfiles/test4.txt");
		System.out.println("\nTest file : test5.txt");
		BB.testSolve("testfiles/test5.txt");
		System.out.println("\n[-----------------------------------------------------------------------------------------------------------------------------]\n\n");
		
		System.out.println("\n\n[----------------------------------Solution of Questions 13 ----------------------------------------------------------------]\n");
		OnlineP.OnlineSolve(10000); // you can change the number of experiences
		System.out.println("\n[---------------------------------------------------------FIN-----------------------------------------------------------------]\n\n");

	}

}

 # User-scheduling-in-5G 

  
         THE PROBLEM 

In 5G, an antenna transmits data packets to smartphones (or users) through a
wireless medium, which is divided into a set of frequency channels.
The higher the power dedicated to a user, the higher the data rate it can experience. 
The exact dependence between power and data rate is however user and channel specific. 
With the same transmit power, a user close to the antenna will
enjoy for example a higher data rate than a user far away.
A wireless packet scheduler is thus responsible to allocated channels to users and to divide the total power
budget of the antenna among the available channels. 
The goal of this project is to design optimal packet schedulers in this context.  

The detailed subject and its tasks as well as test files are available on this page:  
https://marceaucoupechoux.wp.imt.fr/enseignement/english-inf421-pi/  




        PROCESSING
  
Working on an Multiple Choice Knapsack Problem (which is an integer linear problem) to design optimal packet schedulers for 5G   
We will first proceed to a data preprocessing
Then we will use different approaches to solve this ILP :
- A greedy algorithm
- A dynamic programming algorithm
- A branch an bound algorithm
- A stochastic online algorithm





         THE CODE 

The code is written in Java. Please read the following to run it :

- All classes are in the file "Code.java" including all the methods.
	Classes answering to the questions :

		The class Preprocessing answers to each of the questions 1,2,3,4 and 5.   

		The class LinearP answers to the questions 6, 7  

		The class DynamicP answers to the questions 8 and 9 and partially 11  

		The class BB answers to the questions 10 and partially 11  

		The class OnlineP answers to the questions 12 and 13  

- To see the results you can run the file "test.java", which contains all the tests for all classes. 



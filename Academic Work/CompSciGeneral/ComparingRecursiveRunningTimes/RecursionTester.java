/*
 * Timothy Rubio
 * 09-30-2018
 * Lab 03 - Comparing Recursive Running Times
 */

//Based on the two methods we wrote in class, it seems overwhelmingly apparent that the runtime of min2 is MUCH faster than min1. Both functions began giving me 
//a runtime of 0 milliseconds but min1's time grows very quickly as the elements approach our max boundary of 33 elements. This appears to be due to the use of overloading
//and recursion simultaneously in our min2 function.

import java.util.Random;

public class RecursionTester {
	// suppose you are given the following two recursive
	//   methods for finding the minimum of an array:
	// version 1:
	public static int min1(int[] array) 
	{
		return min1(array, array.length);
	}
	
	public static int min1(int[] array, int n) 
	{
		if(n <= 0) {
			// no minimum exists
			throw new IllegalArgumentException("array is empty");
		} else if(n == 1) {
			// base case: if there's a single element, that's it!
			return array[0];
		} else {
			// recursive case: it's either the last element, or the
			//                 minimum of the first n-1 elements
			if(min1(array, n-1) < array[n-1]) {
				return min1(array, n-1);
			} else {
				return array[n-1];
			}
		}
	}
	
	// version 2:
	public static int min2(int[] array) 
	{
		return min2(array, array.length);
	}
	
	public static int min2(int[] array, int n) 
	{
		if(n <= 0) {
			// no minimum exists
			throw new IllegalArgumentException("array is empty");
		} else if(n == 1) {
			// base case: if there's a single element, that's it!
			return array[0];
		} else {
			// recursive case: it's either the last element, or the
			//                 minimum of the first n-1 elements
			int tmp = min2(array, n-1);
			if(tmp < array[n-1]) {
				return tmp;
			} else {
				return array[n-1];
			}
		}
	}
	
	// main method
	public static void main(String[] args) 
	{
		
		
		final int SIZE = 33; //Size variable for array
		long startTime;
		long endTime;
		long totalTime;
		
		int[] array01 = new int[SIZE];
		int[] array02 = new int[SIZE];
		int[] array03 = new int[SIZE];
		
		//Pseudo-random number generator
		Random rand = new Random();
		rand.setSeed(System.currentTimeMillis()); //Using current time as seed
		System.out.println("Elements:		Runtime of min1:	Runtime of Min2:\n");
		//Filling arrays
		for (int i = 0; i < SIZE; i++)
		{
			
			array01[i] += i;
			startTime = System.currentTimeMillis();
			min1(array01);
			endTime = System.currentTimeMillis();
			totalTime = endTime - startTime;
			System.out.print( (i+1) + "			" + totalTime + "			");
			startTime = System.currentTimeMillis();
			min2(array01);
			endTime = System.currentTimeMillis();
			totalTime = endTime - startTime;
			System.out.println(totalTime + "\n");
		}

		int n = 0; //Variable for element counter since in descending order
		for (int i = SIZE-1; i >= 0; i--)
		{
			array02[i] += i;
			n++;
			startTime = System.currentTimeMillis();
			min1(array02);
			endTime = System.currentTimeMillis();
			totalTime = endTime - startTime;
			System.out.print( (n) + "			" + totalTime + "			");
			startTime = System.currentTimeMillis();
			min2(array02);
			endTime = System.currentTimeMillis();
			totalTime = endTime - startTime;
			System.out.println(totalTime + "\n");
		}
		
		for (int i = 0; i < SIZE; i++)
		{
			array03[i] = rand.nextInt(100);
			startTime = System.currentTimeMillis();
			min1(array03);
			endTime = System.currentTimeMillis();
			totalTime = endTime - startTime;
			System.out.print( (i+1) + "			" + totalTime + "			");
			startTime = System.currentTimeMillis();
			min2(array03);
			endTime = System.currentTimeMillis();
			totalTime = endTime - startTime;
			System.out.println(totalTime + "\n");
		}
	}
}

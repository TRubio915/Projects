/*
 * Timothy Rubio
 * 10-09-2018
 * Lab 06 - Iterators
 */

public class Driver
{

	public static void main(String[] args)
	{
		//creating lists
		DoublyLinkedList<Integer> list01 = new DoublyLinkedList<Integer>();
		DoublyLinkedList<Integer> list02 = new DoublyLinkedList<Integer>();
		DoublyLinkedList<Integer> list03 = new DoublyLinkedList<Integer>();
		
		//creating iterators
		TwoWayIterator<Integer> it01 = list01.iterator();
		TwoWayIterator<Integer> it02 = list02.iterator();
		TwoWayIterator<Integer> it03 = list03.iterator();
		
		//filling lists
		for (int i = 0; i <= 10; i++)
		{
			list01.addFirst(i);
		}
		for (int j = 0; j <= 100; j++)
		{
			list02.addLast(j);
		}
		for (int g = 0; g <= 1000; g++)
		{
			list03.addFirst(g);
		}
		
		//testing iterators
		while(it01.hasNext())
		{
			System.out.println(it01.next());
		}
		while(it01.hasPrevious())
		{
			System.out.println(it01.previous());
		}
		while(it02.hasNext())
		{
			System.out.println(it02.next());
		}
		while(it02.hasPrevious())
		{
			System.out.println(it02.previous());
		}
		while(it03.hasNext())
		{
			System.out.println(it03.next());
		}
		while(it03.hasPrevious())
		{
			System.out.println(it03.previous());
		}
		
		//For each loop
		for (Integer foo : list01) //For each loop
		{
			System.out.println(foo);
		}
		
	}
	
}

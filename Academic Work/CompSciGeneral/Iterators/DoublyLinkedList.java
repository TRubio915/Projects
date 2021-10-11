import java.util.Iterator;

public class DoublyLinkedList<E> implements Iterable<E>
{
	private int size;
	private Node<E> head;
	private Node<E> tail;
	
	public DoublyLinkedList()
	{
		size = 0;
		head = new Node( null, null, null ); 
		tail = new Node( null, null, head );
		head.setNext( tail ); //Connects head to tail
	}
	
	public int size() { return size; } //Returns size of list
	
	public boolean isEmpty() //Returns true or false based on if list is empty or not.
	{
		return size == 0;
	}
	
	public E first() //Returns data from first node in list
	{
		return head.getNext().getData();
	}
	
	public E last() //Returns data from last node in list
	{
		return tail.getPrev().getData();
	}
	
	public void addFirst( E value ) //Adds to beginning of list
	{
		addBetween( value, head, head.getNext());
	}
	
	public void addLast( E value ) //Adds to end of list
	{
		addBetween( value, tail.getPrev(), tail );
	}
	
	private void addBetween ( E value, Node<E> prev, Node<E> next ) //adds between next and previous
	{
		Node<E> tempNode = new Node<E>( value, next, prev );
		prev.setNext(tempNode);
		next.setPrev(tempNode);
		size++;
	}
	
	public E removeFirst() //Removes node after head
	{
		return remove( head.getNext() );
	}
	
	public E removeLast() //Removes node before tail
	{
		return remove( tail.getPrev() );
	}
	
	private E remove(Node<E> theNode) //Removes specified node
	{
		if ( isEmpty() )
		{
			return null;
		}
		else
		{
			theNode.getPrev().setNext( theNode.getNext() );
			theNode.getNext().setPrev( theNode.getPrev() );
			size--;
			return theNode.getData();
		}
	}
	
	public boolean search(E element) //Searches for passed in element
	{
		Node<E> tempNode = head.getNext();
		if(isEmpty())
		{
			return false;
		}
		else
		{
			while(tempNode.getNext() != null)
			{
				if(tempNode.getData() == element)
				{
					return true;
				}
			}
			return false;
		}
	}
	
	public String toString() //Displays as string
	{
		Node<E> tempNode = head.getNext();
		String theResult = "";
		while ( tempNode.getNext() != null)
		{
			theResult = theResult + tempNode.getData().toString();
			tempNode = tempNode.getNext();
		}
		return theResult;
	}
	
	public TwoWayIterator<E> iterator()
	{
		return new DoublyLLIterator();
	}
	

private static class Node <E> //Variable for node
{
	private E myData; //Data will be in form of passed in variable
	private Node<E> myNext; //<E> defines the type held by Node
	private Node<E> myPrev;
	
	public Node (E data, Node<E> next, Node<E> prev) //Constructor
	{
		myData = data;
		myNext = next;
		myPrev = prev;
	}
	
	public E getData() 
	{
		return myData;
	}
	
	public Node<E> getNext()
	{
		return myNext;
	}
	
	public void setNext(Node<E> theNode)
	{
		this.myNext = theNode;
	}
	
	public Node<E> getPrev()
	{
		return myPrev;
	}
	
	public void setPrev(Node<E> theNode)
	{
		this.myPrev = theNode;
	}
}


	
	private class DoublyLLIterator<E> implements TwoWayIterator<E>
	{
		private Node<E> current;
		private Node<E> prev;
		
		public DoublyLLIterator()
		{
			current = null;
			prev = null;
		}
		
		public boolean hasNext() 
		{
			return ((current == null && head != null) || (current.myNext != null));
		}
		
		// Implement next
		public E next() 
		{
			if(current == null && head != null) //Iterator just created but list not empty
			{
				current = (Node<E>) head;
				prev = null;
			}
			
			// Advance current and prev for next time
			// and return value of current
			prev = current;
			current = current.myNext;
			return current.myData;
		}

		public boolean hasPrevious()
		{
			return ((current == null && head != null) || (current.myPrev != null));
		}

		public E previous()
		{
			if(current == null && head != null) //Iterator just created but list not empty
			{
				current = (Node<E>) head;
				prev = null;
			}
			
			// Advance current and prev for next time
			// and return value of current
			prev = current;
			current = current.myPrev;
			return current.myData;
		}
	}
	
}

public class BinaryTree<E>
{
	private int mySize;
	private Node<E> myRoot;
	
	public BinaryTree()
	{
		mySize = 0;
		myRoot = null;
	}
	
	private Node<E> validateLocation(Location<E> loc)
	{
		if((loc == null) || !(loc instanceof Node))
		{
			throw new IllegalArgumentException("Not a Valid Location!");
		}
		return (Node<E>) loc;
	}
	
	public Location<E> addRoot(E val)
	{
		if(myRoot != null)
		{
			throw new IllegalStateException("Tree is Not Empty! Cannot Add Root!");
		}
		myRoot = new Node<>(val,null);
		mySize++;
		return myRoot;
	}
	
	public Location<E> addLeftChild(E val, Location<E> loc)
	{
		Node<E> n = validateLocation(loc);
		if(n.myLeftChild != null)
		{
			throw new IllegalStateException("Left Child Already Exists!");
		}
		n.myLeftChild = new Node<>(val, n);
		mySize++;
		return n.myLeftChild;
	}
	
	public Location<E> addRightChild(E val, Location<E> loc)
	{
		Node<E> n = validateLocation(loc);
		if(n.myRightChild != null)
		{
			throw new IllegalStateException("Right Child Already Exists!");
		}
		n.myRightChild = new Node<>(val, n);
		mySize++;
		return n.myRightChild;
	}
	
	public boolean isRoot(Location<E> loc)
	{
		Node<E> n = validateLocation(loc);
		return n == myRoot;
		//return n.parent == null;
	}
	
	public int size()
	{
		return mySize;
	}
	
	public boolean isEmpty()
	{
		return mySize == 0;
	}
	
	public boolean isLeaf(Location<E> loc)
	{
		Node<E> n = validateLocation(loc);
		return ((n.myLeftChild == null) && (n.myRightChild == null));
	}
	
	public boolean isInternal(Location<E> loc)
	{
		//No validation, isLeaf() handles this
		return (!isLeaf(loc));
	}
	
	public String toString()
	{
		StringBuilder sb = new StringBuilder();
		toString(myRoot, sb, 0);
		return sb.toString();
	}
	
	private void toString(Node<E> loc, StringBuilder sb, int level)
	{
		if(loc != null)
		{
			for(int i = 0; i < 2 * level; i++)
			{
				sb.append(" ");
			}
			sb.append(loc.myValue + "\n");
			toString(loc.myLeftChild, sb, level + 1);
			toString(loc.myRightChild, sb, level + 1);
		}
	}
	
	public static int evaluate(BinaryTree<String> theExprTree)
	{
		return evaluateHelper(theExprTree.myRoot);
	}
	
	private static int evaluateHelper(Node<String> theNode)
	{
		if(theNode == null)
		{
			throw new NullPointerException("Tree is Empty!");
		}
		if ((theNode.myLeftChild == null) && (theNode.myRightChild == null))
		{
			return Integer.parseInt(theNode.getValue());
		}
		switch(theNode.getValue())
		{
		case "+":
			return evaluateHelper(theNode.myLeftChild) + evaluateHelper(theNode.myRightChild);
			
		case "*":
			return evaluateHelper(theNode.myLeftChild) * evaluateHelper(theNode.myRightChild);
			
		case "/":
			return evaluateHelper(theNode.myLeftChild) / evaluateHelper(theNode.myRightChild);
			
		case "-":
			return evaluateHelper(theNode.myLeftChild) - evaluateHelper(theNode.myRightChild);
			
		default: 
			throw new IllegalStateException("Invalid Symbol");
		}
	}
	
	public static void printExpression(BinaryTree<String> theTree)
	{
		printExpressionHelper(theTree.myRoot);
	}
	
	public static void printExpressionHelper(Node<String> theNode)
	{
		if(theNode == null)
		{
			throw new NullPointerException("Tree is Empty!");
		}
		if(theNode.myLeftChild != null)
		{
			System.out.print("(");
			printExpressionHelper(theNode.myLeftChild);
		}
		System.out.print(theNode.getValue());
		if(theNode.myRightChild != null)
		{
			printExpressionHelper(theNode.myRightChild);
			System.out.print(")");
		}
	}
	
/*		
	public int add()
	{
		return addHelper(myRoot);
	}
	
	private int addHelper(Node<E> theNode)
	{
		if(theNode == null)
		{
			return 0;
		}
		else
		{
			return theNode.getValue() + addHelper(theNode.left) + addHelper(theNode.right);
		}
	}
*/	
	
	private static class Node<T> implements Location<T>
	{
		private T myValue;
		private Node<T> myParent;
		private Node<T> myLeftChild;
		private Node<T> myRightChild;
		
		public Node(T value, Node<T> parent)
		{
			myValue = value;
			myParent = parent;
			myLeftChild = null;
			myRightChild = null;
		}
		
		public T getValue()
		{
			return myValue;
		}		
		
		public String toString()
		{
			return myValue.toString();
		}
	}
}

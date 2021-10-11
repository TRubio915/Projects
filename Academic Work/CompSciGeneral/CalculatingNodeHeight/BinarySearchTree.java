import java.io.Serializable;

public class BinarySearchTree<E extends Comparable<E>> implements Serializable
{
	private transient int mySize;
	public transient int height = 0;
	private Node<E> myRoot;
	
	public BinarySearchTree()
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
	
	public Location<E> insert(E val)
	{
		if(this.isEmpty())
		{
			return addRoot(val);
		}
		else
		{
			return insertHelper(val, myRoot);
		}
	}
	
	private Location<E> insertHelper(E val, Location<E> loc)
	{
		Node<E> theNode = validateLocation(loc);
		if(val.compareTo(loc.getValue()) <= 0)
		{
			if(theNode.myLeftChild == null)
			{
				return addLeftChild(val, theNode);
			}
			else
			{
				return insertHelper(val, theNode.myLeftChild);
			}
		}
		else
		{
			if(theNode.myRightChild == null)
			{
				return addRightChild(val, theNode);
			}
			else
			{
				return insertHelper(val, theNode.myRightChild);
			}
		}
	}
	
	public boolean contains(E val)
	{
		if(this.isEmpty())
		{
			return false;
		}
		else
		{
			return containsHelper(val, myRoot);
		}
	}
	
	private boolean containsHelper(E val, Location<E> loc)
	{
		Node<E> theNode = validateLocation(loc);
		if(val.compareTo(loc.getValue()) == 0)
		{
			return true;
		}
		if(val.compareTo(loc.getValue()) <= 0)
		{
			if(theNode.myLeftChild == null)
			{
				return false;
			}
			else
			{
				return containsHelper(val, theNode.myLeftChild);
			}
		}
		else
		{
			if(theNode.myRightChild == null)
			{
				return false;
			}
			else
			{
				return containsHelper(val, theNode.myRightChild);
			}
		}
	}
	
	public boolean delete(E val)
	{
		if(this.isEmpty())
		{
			return false;
		}
		else
		{
			return deleteHelper(val, myRoot);
		}
	}
	
	private boolean deleteHelper(E val, Location<E> loc)
	{
		Node<E> theNode = validateLocation(loc);
		if(val.compareTo(loc.getValue()) == 0)
		{
			Node<E> parent = theNode.myParent;
			if(isLeaf(theNode))
			{
				disown(parent, theNode);
			}
			else
			{
				if((theNode.myLeftChild == null) || (theNode.myRightChild == null)) //One Child
				{
					if(theNode.myLeftChild != null) //One child on the Left
					{
						if(parent.myLeftChild == theNode)
						{
							parent.myLeftChild = theNode.myLeftChild;
							theNode.myLeftChild.myParent = parent;
						}
						else
						{
							parent.myRightChild = theNode.myLeftChild;
							theNode.myLeftChild.myParent = parent;
						}
					}
					else //One child on the right of theNode
					{
						if(parent.myLeftChild == theNode)
						{
							parent.myLeftChild = theNode.myRightChild;
							theNode.myRightChild.myParent = parent;
						}
						else
						{
							parent.myRightChild = theNode.myRightChild;
							theNode.myRightChild.myParent = parent;
						}
					}
				}
				else
				{
					//Two Children
					Node<E> theNodeToPromote = findLargestOnLeft(theNode.myLeftChild);
					deleteHelper(theNodeToPromote.getValue(), theNodeToPromote);
					theNode.myValue = theNodeToPromote.myValue;
				}
			}
			return true;
		}
		if(val.compareTo(loc.getValue()) <= 0)
		{
			if(theNode.myLeftChild == null)
			{
				return false;
			}
			else
			{
				return deleteHelper(val, theNode.myLeftChild);
			}
		}
		else
		{
			if(theNode.myRightChild == null)
			{
				return false;
			}
			else
			{
				return deleteHelper(val, theNode.myRightChild);
			}
		}
	}
	
	private Node<E> findLargestOnLeft(Node<E> theNode)
	{
		//Take one step to the left
		if(theNode.myRightChild != null)
		{
			return findLargestOnLeft(theNode.myRightChild);
		}
		else
		{
			return theNode;
		}
	}
	
	private void disown(Node<E> parent, Node<E> child)
	{
		if(parent.myLeftChild == child)
		{
			parent.myLeftChild = null;
		}
		else if(parent.myRightChild == child)
		{
			parent.myRightChild = null;
		}
		else
		{
			throw new IllegalStateException("Disown: You are NOT the father!");
		}
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
	
	public static int evaluate(BinarySearchTree<String> theExprTree)
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
	
	public static void printExpression(BinarySearchTree<String> theTree)
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
	
	//Lab 09 starts here
	
	//Compares heights for final determination
	public int calcNodeHeight(Location<E> loc)
	{
		int rightHeight = calcHeightRight(loc, 0);
		int leftHeight = calcHeightLeft(loc, 0);
		if (leftHeight > rightHeight)
		{
			return leftHeight;
		}
		else
		{
			return rightHeight;
		}
	}
	
	//Checks height based on right subtree
	private int calcHeightRight(Location<E> loc, int rHeight)
	{
		Node<E> theNode = validateLocation(loc);
		if(isLeaf(theNode))
		{
			return rHeight;
		}
		else
		{
			rHeight++;
			return calcHeightRight(theNode.myRightChild, rHeight);
		}
	}
	
	//Checks height based on left subtree
	private int calcHeightLeft(Location<E> loc, int lHeight)
	{
		Node<E> theNode = validateLocation(loc);
		if(isLeaf(theNode))
		{
			return lHeight;
		}
		else
		{
			lHeight++;
			return calcHeightLeft(theNode.myLeftChild, lHeight);
		}
	}
	
	//Compares heights for a balance check
	public boolean balanceCheck(Location<E> loc)
	{
		Node<E> theNode = validateLocation(loc);
		int lChildHeight = calcHeightLeft(theNode.myLeftChild, 0);
		int rChildHeight = calcHeightRight(theNode.myRightChild, 0);
		if((lChildHeight < rChildHeight - 1) || rChildHeight < lChildHeight - 1)
		{
			return false;
		}
		else
		{
			return true;
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
	
	private static class Node<T> implements Location<T>, Serializable
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

/*
 * Timothy Rubio
 * 10/18/2018
 * Lab 09 - Tree Traversal
 */

public class TestDriver
{

	public static void main(String[] args)
	{
		BinaryTree<String> theTree = new BinaryTree();
		
		Location<String> root = theTree.addRoot("*");
		Location<String> lChild = theTree.addLeftChild("2", root);
		Location<String> rChild = theTree.addRightChild("+", root);
		Location<String> rlChild = theTree.addLeftChild("3", rChild);
		Location<String> llChild = theTree.addRightChild("4", rChild);
		
		BinaryTree.printExpression(theTree);
		System.out.print(" = ");
		System.out.println(BinaryTree.evaluate(theTree));
	}
}

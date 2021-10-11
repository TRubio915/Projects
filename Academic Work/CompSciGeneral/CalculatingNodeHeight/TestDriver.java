import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

public class TestDriver
{

	public static void main(String[] args)
	{
		BinarySearchTree<Integer> theTree = new BinarySearchTree<Integer>();
		
		Location<Integer> root = theTree.addRoot(5);
		
		theTree.insert(3);
		theTree.insert(1);
		theTree.insert(4);
		theTree.insert(8);
		theTree.insert(15);
		theTree.insert(12);
		theTree.insert(18);
		
		System.out.println("Height of Root: " + theTree.calcNodeHeight(root));
		
		System.out.println("Tree Balanced?: " + theTree.balanceCheck(root));
		
		System.out.println(theTree);
		
		/*try
		{
		FileOutputStream fs = new FileOutputStream("Result.ser");
		ObjectOutputStream os = new ObjectOutputStream(fs);
		os.writeObject(theTree);
		os.close();
		}
		catch(Exception ex)
		{
			ex.printStackTrace();
		}
		
		try
		{
			FileInputStream fs = new FileInputStream("Result.ser");
			ObjectInputStream is = new ObjectInputStream(fs);
			BinarySearchTree otherTree = (BinarySearchTree)is.readObject();
			System.out.println(otherTree);
			is.close();
		}
		catch(Exception ex)
		{
			ex.printStackTrace();
		}*/
	}
}

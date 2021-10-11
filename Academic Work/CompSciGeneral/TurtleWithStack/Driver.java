/*
 * Timothy Rubio
 * 09/27/2018
 * Lab 04 - Turtle with Stack
 */

import java.awt.Point;
import java.util.Scanner;

public class Driver 
{
	public static void main(String[] args) 
	{
		// set up standard draw (creates an 800x800 canvas)
		StdDraw.setCanvasSize(800, 800);
		StdDraw.setXscale(0, 800);
		StdDraw.setYscale(0, 800);
		
		// starts a turtle in the center of the canvas
		StackTurtle st = new StackTurtle(new Point(400,400));
		
		// read in commands from the keyboard
		String command;
		Scanner keyboard = new Scanner(System.in);
		do 
		{
			// get next command from the keyboard
			command = keyboard.next();
			
			// process that command
			if (command.equalsIgnoreCase("["))
			{
				st.pushCurrentPosition();
			}
			else if (command.equalsIgnoreCase("]"))
			{
				st.popCurrentPosition();
			}
			else if(command.equalsIgnoreCase("f"))
			{
				st.moveForward(keyboard.nextDouble());
			} 
			else if(command.equalsIgnoreCase("r"))
			{
				st.rotate(keyboard.nextDouble());
			} 
			else 
			{
				System.out.println("Unknown command: " + command);
			}
		} while(!command.equalsIgnoreCase("q"));
		
		// close the keyboard scanner and exit the program
		keyboard.close();
		System.exit(0);
	}

}

import java.awt.Point;
import java.util.Stack;

public class StackTurtle {
	// private member variables
	private Point currentPosition;
	private double angle;
	private Stack<Point> turtleStack = new Stack();
	
	// constructor: creates a turtle at the specified position, pointing up
	public StackTurtle(Point startingPosition) {
		currentPosition = startingPosition;
		angle = 90;
	}
	
	public void moveForward(double distance) {
		// calculate the next position
		Point newPosition = new Point(
				(int)(currentPosition.x + distance*Math.cos(Math.toRadians(angle))),
				(int)(currentPosition.y + distance*Math.sin(Math.toRadians(angle))));
		
		// draw a line from the current position to the new position
		StdDraw.line(currentPosition.x, currentPosition.y, newPosition.x, newPosition.y);
		
		// update the current position
		currentPosition = newPosition;
	}
	
	public void rotate(double degrees) {
		angle += degrees;
	}
	
	public void pushCurrentPosition()
	{
		turtleStack.push(currentPosition);
	}
	
	public void popCurrentPosition()
	{
		currentPosition = turtleStack.pop();
	}
}

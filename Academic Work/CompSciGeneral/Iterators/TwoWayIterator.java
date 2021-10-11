import java.util.Iterator;

public interface TwoWayIterator<ValueType> extends Iterator<ValueType> {
	public boolean hasPrevious();
	public ValueType previous();
}

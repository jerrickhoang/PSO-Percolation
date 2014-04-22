import java.util.Comparator;

// this is used when sorting DataOutput objects in ascending function-value order in IntervalSummaryData

public class DataOutputFunctionValueComparator implements Comparator<DataOutput> {

	// smaller function value "comes before"
	public int compare(DataOutput d1, DataOutput d2) {
		if (d1.getFunctionValue() < d2.getFunctionValue())
			return -1;
		else if (d1.getFunctionValue() > d2.getFunctionValue())
			return 1;
		else
			return 0;
	}

	
}

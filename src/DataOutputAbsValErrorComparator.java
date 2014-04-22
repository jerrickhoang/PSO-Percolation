import java.util.Comparator;

// this is used when sorting DataOutput objects in ascending absolute value error order in IntervalSummaryData

public class DataOutputAbsValErrorComparator implements Comparator<DataOutput> {

	// smaller error "comes before"
	public int compare(DataOutput d1, DataOutput d2) {
		if (d1.getAbsValError() < d2.getAbsValError())
			return -1;
		else if (d1.getAbsValError() > d2.getAbsValError())
			return 1;
		else
			return 0;
	}

	
}

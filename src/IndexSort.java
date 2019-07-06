import java.util.Arrays;
import java.util.Comparator;

public class IndexSort {
	
	public Integer[] Ind(double[] a){	
		Integer[] b = new Integer[a.length];
		for (int i = 0 ; i < b.length ; i++){
			b[i] = i; 
		}
		Arrays.sort(b, new Comparator<Integer>() {
			public int compare(Integer o1, Integer o2) {
			 return Double.valueOf(a[o1]).compareTo(Double.valueOf(a[o2]));
			 }
		});
		
		return b;
	}

}

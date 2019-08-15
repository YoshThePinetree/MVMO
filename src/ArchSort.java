import java.util.Arrays;
import java.util.Comparator;

// Sort the Archive fitness and the variables
public class ArchSort {
	class IndexSort1 extends IndexSort{
	}
	
	IndexSort1 Isort = new IndexSort1();
	
	public void sort(Archive arc, double F, double[] Xn, int a, int d){	
		double [] F2 = new double [a+1];	// Past fitness + new fitness
		double [][] X2 = new double [a+1][d];	// Past variables + new fitness
		Integer [] F2ind = new Integer[a+1];
	 
		for(int m=0; m<a; m++){
			F2[m] = arc.fitness[m];
			for(int j=0; j<d; j++){
				X2[m][j] = arc.variable[m][j];
			}
		}
		F2[a]=F;
		for(int j=0; j<d; j++){
			X2[a][j] = Xn[j];
		}
	 
		F2ind=Isort.Ind(F2);	//Sorted index for the past fitness + new fitness
		for(int m=0; m<a; m++){	// Put all values into the Archive except the worst fitness
			arc.fitness[m] = F2[F2ind[m]];
			for(int j=0; j<d; j++){
				arc.variable[m][j] = X2[F2ind[m]][j];
			}
		}
	}

}

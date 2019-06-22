// Normalization class
public class NormVars {
	//Normalize each variable between 0 and 1
	public double[][] Normalize(double[][] X, double[] xlim){
		double[][] X1 = new double[X.length][X[0].length];
		
		for(int i=0; i<X.length; i++) {
			for(int j=0; j<X[0].length; j++) {
				X1[i][j]=(X[i][j] - xlim[2*j]) / (xlim[2*j+1] - xlim[2*j]);
			}
		}

		return X1;
	}
	
	public double[][] DeNormalize(double[][] X, double[] xlim){
		double[][] X1 = new double[X.length][X[0].length];
		
		for(int i=0; i<X.length; i++) {
			for(int j=0; j<X[0].length; j++) {
				X1[i][j]=X[i][j]*(xlim[2*j+1] - xlim[2*j]) + xlim[2*j];
			}
		}

		return X1;
	}
}

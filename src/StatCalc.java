// Obtain statistic values
public class StatCalc {
		
	// Calculation of mean of all elements in X
	public double Mean(double[] X){
		int len = X.length;
		double sum = 0;
		double y = 0;
		
		for(int i=0; i<len; i++) {
			sum = sum + X[i];
		}
		
		y = (double) sum/len;		
		return y;
	}
	
	// Calculation of variance of all elements in X
	public double Var(double[] X){
		int len = X.length;
		double summean = 0;
		double y = 0;
		
		StatCalc stat = new StatCalc();
		double mu = stat.Mean(X);
		
		for(int i=0; i<len; i++) {
			summean = summean + Math.pow((X[i] - mu),2);
		}
		
		y = (double) summean/len;		
		return y;
	}
	
}

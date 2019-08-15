// Objective function class
public class ObjFunc {	
	// Method to evaluate the objective function
	public static double EvalFunc(double[] x, int a, int d){
		double f = 0;
		double st1=0,st2=0;
		
		switch(a) {		
		case 1:		//Rosenbrock function constrained with a cubic and a line
			// round in the ranges
			if(x[0]<-1.5) {
				x[0] = -1.5;
			}
			if(x[0]>1.5) {
				x[0] = 1.5;
			}
			if(x[1]<-0.5) {
				x[1] = -0.5;
			}
			if(x[1]>2.5) {
				x[1] = 2.5;
			}
			
			st1 = Math.pow((x[0]-1),3) - x[1] + 1;
			st2 = x[0] + x[1] - 2;
		
			if(st1<=0 && st2<=0) {
				f = Math.pow((1-x[0]),2) + 100*Math.pow((x[1]-Math.pow(x[0], 2)),2);
			}else{
				f = 1000;
			}
			break;
			
		case 2:		// Rosenbrock function constrained to a disk
			// round in the ranges
			if(x[0]<-1.5) {
				x[0] = -1.5;
			}
			if(x[0]>1.5) {
				x[0] = 1.5;
			}
			if(x[1]<-1.5) {
				x[1] = -0.5;
			}
			if(x[1]>1.5) {
				x[1] = 1.5;
			}
			
			st1 = x[0]*x[0] + x[1]*x[1];
		
			if(st1<=2) {
				f = Math.pow((1-x[0]),2) + 100*Math.pow((x[1]-Math.pow(x[0], 2)),2);
			}else{
				f = 1000;
			}
			break;
					
		case 3:		// Mishra's Bird function - constrained
			if(x[0]<-10) {
				x[0] = -10;
			}
			if(x[0]>0) {
				x[0] = 0;
			}
			if(x[1]<-6.5) {
				x[1] = -6.5;
			}
			if(x[1]>0) {
				x[1] = 0;
			}
			
			st1 = Math.pow(x[0]+5,2) + Math.pow(x[1]+5,2);
			
			if(st1<25) {
				f = Math.sin(x[1])*Math.exp(Math.pow(1-Math.cos(x[0]),2)) + Math.cos(x[0])*Math.exp(Math.pow(1-Math.sin(x[1]),2)) + Math.pow(x[0]-x[1],2);
			}else{
				f = 1000;
			}
			break;
			
		case 4:		// Simionescu function
			// round in the ranges
			if(x[0]<-1.25) {
				x[0] = -1.25;
			}
			if(x[1]>1.25) {
				x[1] = 1.25;
			}
			
			double rt=1, rs=0.2;
			int n=8;
			double u=Math.pow( rt + rs*Math.cos(n*Math.atan(x[0]/x[1])), 2);
			st1 = x[0]*x[0] + x[1]*x[1];
			
			if(st1<=u) {
				f = 0.1*x[0]*x[1];
			}else{
				f = 1000;
			}
			break;
			
		case 5:		// 2-D Rastrigin function
			int A=10;
			f = A*d;
			
			for(int i=0; i<d; i++){ // round in the ranges and fitness calculation
				if(x[i]<-5.12) {
					x[i] = -5.12;
				}
				if(x[i]>5.12) {
					x[i] = 5.12;
				}
				f = f + (x[i]*x[i] - A*Math.cos(2*Math.PI*x[i]));
			}
		}
		return f;
	}
	
	// Method to display
	public static void FuncName(int x){
		switch(x) {
		case 1:
			System.out.println("Rosenbrock function constrained with a cubic and a line");
			break;
		case 2:
			System.out.println("Rosenbrock function constrained to a disk");
			break;
		case 3:
			System.out.println("Mishra's Bird function - constrained");
			break;
		case 4:
			System.out.println("Simionescu function");
			break;
		case 5:
			System.out.println("2-D Rastrigin function");
			break;
		default :
			System.out.println("Error");
			break;
		}
	}
	
	// Method to give variable limitation
	public static double[][] GetLimit(int x, int d){
		double[][] xylim = new double [2][d];
		
		switch(x) {
		case 1:
			xylim[0][0] = -1.5;
			xylim[1][0] = 1.5;
			xylim[0][1] = -0.5;
			xylim[1][1] = 2.5;

			break;
		case 2:
			xylim[0][0] = -1.5;
			xylim[1][0] = 1.5;
			xylim[0][1] = -0.5;
			xylim[1][1] = 1.5;
						
			break;
		case 3:
			xylim[0][0] = -10;
			xylim[1][0] = 0;
			xylim[0][1] = -6.5;
			xylim[1][1] = 0;
			
			break;
		case 4:
			xylim[0][0] = -1.25;
			xylim[1][0] = 1.25;
			xylim[0][1] = -1.25;
			xylim[1][1] = 1.25;
	
			break;
		case 5:
			for(int i=0; i<2; i++){
				for(int j=0; j<d; j++){
					if(i == 0) {
						xylim[i][j] = -5.12;					
					}else {
						xylim[i][j] = 5.12;
					}
				}
			}
			
			break;
		default :
			System.out.println("Error");
			break;
		}
		
		return xylim;
	}
	
}

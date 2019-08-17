// Objective function class
public class ObjFunc {	
	// Method to evaluate the objective function
	public double LevyW(double x) {
		double y=0;
		y = 1 + ((x-1)/4);
		
		return y;
	}
	
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
			
		case 5:		// Rastrigin function
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
			break;
		
		case 6:		// Rosenbrock function
			for(int i=0; i<d-1; i++){ // round in the ranges and fitness calculation
				f = f + (100*Math.pow(x[i+1] - (x[i]*x[i]),2)) + Math.pow(1 - x[i]*x[i],2); 
			}
			break;
			
		case 7:		// Styblinski-Tang function
			for(int i=0; i<d; i++){ // round in the ranges and fitness calculation
				if(x[i]<-5.0) {
					x[i] = -5.0;
				}
				if(x[i]>5.0) {
					x[i] = 5.0;
				}
				f = f + Math.pow(x[i],4) - 16*Math.pow(x[i],2) + 5*x[i]; 
			}
			f = f / 2;
			break;
			
		case 8:		// Ackley function
			double b1 = 20, b2 = 0.2, c = Math.PI;
			double f1 = 0, f2 = 0;
			
			for(int i=0; i<d; i++){ // round in the ranges and fitness calculation
				if(x[i]<-5.0) {
					x[i] = -5.0;
				}
				if(x[i]>5.0) {
					x[i] = 5.0;
				}
				f1 = f1 + x[i]*x[i];
				f2 = f2 + Math.pow(Math.cos(c*x[i]),2);
			}
			f = (-b1 * Math.exp(-b2 * Math.sqrt(f1/d))) - (Math.exp(f2/d)) + b1 + Math.exp(1);
			break;
			
		case 9:		// Levy function
			double f3 = 0;
			ObjFunc obj = new ObjFunc();
			
			for(int i=0; i<d-1; i++){ // round in the ranges and fitness calculation
				if(x[i]<-10.0) {
					x[i] = -10.0;
				}
				if(x[i]>10.0) {
					x[i] = 10.0;
				}
				f3 = f3 + (Math.pow(obj.LevyW(x[i]) - 1,2)) * (1 + 10*Math.pow(Math.sin(Math.PI * obj.LevyW(x[i]) + 1), 2));
			}
			f = Math.pow(Math.sin(Math.PI * obj.LevyW(x[0])), 2) + f3 + ((Math.pow(obj.LevyW(x[(x.length)-1]) - 1, 2)) * (1 + Math.pow(Math.sin(2*Math.PI*obj.LevyW(x[(x.length)-1])), 2)));
					
			break;
			
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
			System.out.println("Rastrigin function");
			break;
		case 6:
			System.out.println("Rosenbrock function");
			break;
		case 7:
			System.out.println("Styblinski-Tang function");
			break;
		case 8:
			System.out.println("Ackley function");
			break;
		case 9:
			System.out.println("Levy function");
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
			
		case 6:
			for(int i=0; i<2; i++){
				for(int j=0; j<d; j++){
					if(i == 0) {
						xylim[i][j] = -10;					
					}else {
						xylim[i][j] = 10;
					}
				}
			}
			break;
			
		case 7:
			for(int i=0; i<2; i++){
				for(int j=0; j<d; j++){
					if(i == 0) {
						xylim[i][j] = -5;					
					}else {
						xylim[i][j] = 5;
					}
				}
			}
			break;
			
		case 8:
			for(int i=0; i<2; i++){
				for(int j=0; j<d; j++){
					if(i == 0) {
						xylim[i][j] = -5;					
					}else {
						xylim[i][j] = 5;
					}
				}
			}
			break;
			
		case 9:
			for(int i=0; i<2; i++){
				for(int j=0; j<d; j++){
					if(i == 0) {
						xylim[i][j] = -10;					
					}else {
						xylim[i][j] = 10;
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

// Main function

public class MVMOmain {
	public static void main (String arg[]) {
		System.out.println("MVMO: Mean Variance Mapping Optimization");
		
		///////////////////////////
		// Initial configuration //
		///////////////////////////
		
		// #Benchmark number #The global minimum 
		// 1:Rosenbrock function constrained with a cubic and a line, f(1.0,1.0)=0
		// 2:Rosenbrock function constrained to a disk, f(1.0,1.0)=0
		// 3:Mishra's Bird function - constrained, f(-3.1302468,-1.5821422)=-106.7645367
		// 4:Simionescu function, f(+-0.84852813,-+0.84852813)=-0.072
		// 5:Rastrigin function, f(0, ..., 0)=0
		// 6:Rosenbrock function, f(1, ..., 1)=0
		// 7:Styblinski-Tang function, f(-2.903534, ..., -2.903534)=-39.16599*dim
		// 8:Ackley function, f(0, ..., 0)=0
		// 9:Levy function, f(1, ..., 1)=0
			
		int fnum=9;			// the function number to solve
		int dim=2;			// the number of dimension of the decision variable
		int trial=5;		// the number of trials with different random initial
		int ite=1000;		// the number of iterations for a trial
		int pop=100;		// the number of particles
		int rseed=1;		// random seed: MT method			
		int an=50;			// the archive size
		double gpi=0.7;		// the GP selection rate at first
		double gpf=0.2;		// the GP selection rate at last
		double sf=0.95;		// the initial scaling factor
		double af=2.0;		// the asymmetry factor
		int mi=15;			// the first mutation representatives
		int mf=6;			// the final mutation representatives
		
		
		System.out.printf("General Parameters:\n");
		System.out.printf("trial\titeration\tparticles\n");
		System.out.printf("%d\t%d\t\t%d\n",trial,ite,pop);
		System.out.printf("an\tgpi\tgpf\tsf\taf\tmi\tmf\n");
		System.out.printf("%d\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\n",an,gpi,gpf,sf,af,mi,mf);
		
		System.out.printf("Function:\n");
		ObjFunc.FuncName(fnum);
		
		// Answer check
		double[] XX = {1, 1};
		double ans=ObjFunc.EvalFunc(XX, fnum, dim);
		System.out.printf("%f\n",ans);
		
		/////////////////////////////////////////////////
		// Individual Generation & Initial Evaluation  //
		/////////////////////////////////////////////////
		double[][] Frec = new double [ite][trial];
		double[][][] X1rec = new double [ite][dim][trial];
		
		Sfmt rnd = new Sfmt(rseed);
		NormVars norm = new NormVars();
		
		double gp;
		int GP;
		int mn1, mn;

		for(int l=0;l<trial;l++){
			
			double[][] X = new double[pop][dim];	// the decision variable matrix
			double[][] Xn = new double[pop][dim];	// the NORMALIZED decision variable matrix
			double[] F = new double[pop];		// the objective function for the personal best
			double[] Gb = new double[dim];		// the global best vector
			double[][] xylim=ObjFunc.GetLimit(fnum,dim);
			
			double minfnc=1000;
			int minind=0;
			
			for(int i=0;i<pop;i++) {
				for(int j=0;j<dim;j++) {
					X[i][j]=rnd.NextUnif()*(Math.abs(xylim[0][j])+Math.abs(xylim[1][j]))-Math.abs(xylim[0][j]);
				}
				F[i]=ObjFunc.EvalFunc(X[i], fnum, dim);	// variable evaluation
//				System.out.println(F[i]);
				if(minfnc>F[i]) {
					minfnc=F[i];
					minind=i;
					for(int j=0; j<dim; j++) {
						Gb[j]=X[minind][j];
					}
				}
			}
			
			Xn=norm.Normalize(X, xylim);	// Normalization of the variables between 0 and 1

			// restoration of the individuals into the archives
			// archive exists for every individual
			Archive[] arc = new Archive[pop];
			Archive.DevArc(an,dim);
			for(int i=0; i<arc.length; i++){
				arc[i] = new Archive();
				arc[i].fitness[0] = F[i];
				for(int j=0; j<dim; j++) {
					arc[i].variable[0][j] = Xn[i][j];
					arc[i].mean[j] = Xn[i][j];
					arc[i].variance[j] = 1;
					arc[i].shape[j] = (-1) * Math.log(arc[i].variance[j]) * sf;
				}				
			}
			
			/////////////////////////////
			//  MVMO solution search   //
			/////////////////////////////
			IndexSort Isort = new IndexSort();
			Integer [] anum = new Integer[an];
			anum=Isort.Ind(F);
//			System.out.println(Arrays.toString(anum));
			
//			System.out.println(F[anum[1]]);	// refer by this manner
			
			for(int k=0; k<ite; k++) {		// MVMO loop start
//				System.out.printf("Iteration %d\n",k+1);
				gp=gpi + ((k+1)*(gpf-gpi)) / (ite);
		        GP=(int) Math.round(pop*gp);                   // A number of the Good Particles
		        int[] q = new int [GP];
		        for (int i=0; i<GP; i++) {
	        		q[i]=rnd.NextInt(GP);	// random index for choosing parents
		        }	        
		        
		        ////////////////////////////////////////////////
		        // Crossover in good and bad population group //
		        ////////////////////////////////////////////////
		        double[][] Xp = new double [pop][dim];
		        double[][] Xo = new double [pop][dim];
		        double[][] Xm = new double [pop][dim];
		        double[][] Xs1 = new double [pop][dim];
		        double[][] Xs2 = new double [pop][dim];
		        
		        for(int i=0; i<GP; i++) {	// loop for the GOOD solution
	        		for(int j=0; j<dim; j++){
	        			Xp[i][j] = arc[anum[i]].variable[0][j];
	        			Xm[i][j] = arc[anum[q[i]]].mean[j];
	        			if(Xp[i][j]<Xm[i][j]){
	        				Xs1[i][j] = (arc[anum[i]].shape[j]);
	        				Xs2[i][j] = (arc[anum[i]].shape[j]) * af;
	        			}
	        			else{
	        				Xs1[i][j] = (arc[anum[i]].shape[j]) * af;
	        				Xs2[i][j] = (arc[anum[i]].shape[j]);	        				
	        			}
	        			
	        		}
		        }
		        for(int i=GP; i<pop; i++) {	// loop for the BAD solution: Multi-parent crossover
		        	for(int j=0; j<dim; j++){
		        		int stop = 1;
		        		double beta = 0, l1= -1;
		        		while(stop==1) {
		        			beta = 2 * (rnd.NextUnif() - (0.5*(1-Math.pow((k+1)/ite,2))));
		        			l1 = Xp[q[rnd.NextInt(GP)]][j] + (beta * Xp[0][j] - Xp[GP-1][j]); 
		        			
		        			if(l1>=0 && l1<=1) {
		        			 break;	
		        			}
		        		}
		        		Xp[i][j]=l1;
		        		
		        		if(Xp[i][j]<Xm[i][j]){
	        				Xs1[i][j] = (arc[anum[i]].shape[j]);
	        				Xs2[i][j] = (arc[anum[i]].shape[j]) * af;
	        			}
	        			else{
	        				Xs1[i][j] = (arc[anum[i]].shape[j]) * af;
	        				Xs2[i][j] = (arc[anum[i]].shape[j]);	        				
	        			}
		        		Xm[i][j]=Xp[i][j];
		        		
		        	}
		        }
		        Xo=Xp;	// parents to offsprings
		        
		        ///////////////////////////////////// END OF CROSSOVER
		        
		        ///////////////////////////////////
		        // Mutation via mapping function //
		        ///////////////////////////////////
		        
		        mn1=(int) Math.round((mi) - (((k+1)*(mi-mf))/ite));
	        	mn=(int) Math.round(mf + (rnd.NextUnif()*(mf+mn1-mf)));		// The number of mutated particles
	        	
	        	// nonlinear functions
	        	double hx = 0;
	        	double h0 = 0;
	        	double h1 = 0;
	        	
	        	for(int i=0; i<mn; i++) {		// loop for the mutation
	        		int r1 = rnd.NextInt(pop);
	        		double q1 = rnd.NextUnif();
        			for(int j=0; j<dim; j++) {
        				hx = (Xm[r1][j]*(1-Math.exp((-1)*q1*Xs1[r1][j]))) + ((1-Xm[r1][j])*Math.exp((-1)*(1-q1)*Xs2[r1][j]));
        				h0 = ((1-Xm[r1][j])*Math.exp((-1)*Xs2[r1][j]));
        				h1 = (Xm[r1][j]*(1-Math.exp((-1)*Xs1[r1][j]))) + (1-Xm[r1][j]);
        				Xo[r1][j] = Math.abs(hx + ((1-h1+h0)*q1) - h0 );
        				// rounding the normalized variable
            			if(Xo[r1][j]>1) {
            				Xo[r1][j]=1;
            			}else if(Xo[r1][j]<0){
            				Xo[r1][j]=0;
            			}
        			}
	        	}
	    
	        	///////////////////////////////////// END OF MUTATION
	        	
		        ///////////////////////////
		        // Individual Evaluation //
		        ///////////////////////////
	        	
	        	X=norm.DeNormalize(Xo, xylim);	// De-Normalization of the variables
	        	
	        	for(int i=0;i<pop;i++) {
					F[i]=ObjFunc.EvalFunc(X[i], fnum, dim);	// variable evaluation
//					System.out.println(F[i]);
					if(minfnc>F[i]) {
						minfnc=F[i];
						minind=i;
						for(int j=0; j<dim; j++) {
							Gb[j]=X[minind][j];
						}
					}
				}
	        	
	        	Xn=norm.Normalize(X, xylim);	// Normalization of the variables
	        	
        		///////////////////////////////////// END OF EVALUATION
	        	
				//////////////////////////////////////////////
	        	// Individual Preservation into the Archive //
	        	//////////////////////////////////////////////
	        	StatCalc stat = new StatCalc();
	        	ArchSort asort = new ArchSort();
	        	
	        	// Preservation
				if(k<an-1) {	// if the archive is not full
					int len = k+2;
					for(int i=0; i<pop; i++) {
	        			arc[anum[i]].fitness[k+1] = F[i];
	        			
		        		for(int j=0; j<dim; j++){
		    				arc[anum[i]].variable[k+1][j] = Xn[i][j];
		    				
		    		        double[] X1 = new double [len];
		    		        for(int n=0; n<len; n++){
		    		        	X1[n] = arc[anum[i]].variable[n][j];
		    		        }				
		    				arc[anum[i]].mean[j] = stat.Mean(X1);
		    				arc[anum[i]].variance[j] = stat.Var(X1);
		    				arc[anum[i]].shape[j] = (-1) * Math.log(arc[i].variance[j]) * sf;	    				
		        		}
					}
				
				}else{		// if the archive is already full : select by sort
					for(int i=0; i<pop; i++) {	
						 asort.sort(arc[anum[i]], F[i], Xn[i], an, dim);
						 
						 for(int j=0; j<dim; j++){	// Statistical amount calculation for new archive
			    		        double[] X1 = new double [an];
			    		        for(int n=0; n<an; n++){
			    		        	X1[n] = arc[anum[i]].variable[n][j];
			    		        }				
			    				arc[anum[i]].mean[j] = stat.Mean(X1);
			    				arc[anum[i]].variance[j] = stat.Var(X1);
			    				arc[anum[i]].shape[j] = (-1) * Math.log(arc[i].variance[j]) * sf;	    				
		        		 }
					}
				}
				
	        	// Archive sort
				if(k<an-1) {	// if the archive is not full
					int len = k+1;
					for(int i=0; i<pop; i++) {	
	 					 asort.sort(arc[anum[i]], F[i], Xn[i], len, dim);
					}
				}
	        	
				double [] F1 = new double [pop];
				for(int i=0; i<pop; i++){	
					F1[i] = arc[i].fitness[0];
				}
				anum=Isort.Ind(F1);		// new sorted number
				Frec[k][l]=F1[anum[0]];
				
				double [][] X1 = new double [1][dim];
				double [][] X2 = new double [1][dim];
				for(int j=0; j<dim; j++){
					X1[0][j] = arc[anum[0]].variable[0][j];
				}
				X2 = norm.DeNormalize(X1, xylim);	// De-Normalization of the variables
				for(int j=0; j<dim; j++){
					X1rec[k][j][l]=X2[0][j];
				}
					
				/*
				System.out.printf("Iteration number:\t");
				System.out.printf("%d-",k+1);
				System.out.println(l+1);
				System.out.printf("The best fitness:\t");
				System.out.println(F1[anum[0]]);
				*/
			}
		}		
		        
		
		// write out the objective function value & search log
		String pathname = "C:\\result\\MVMO\\";
		String sufi= ".csv";
		String fnameF = "Fitness";
		
		WriteResult.Output(Frec, ite, trial, pathname + fnameF + sufi);
		
		double[][] Xrec = new double [ite][dim];
		for(int m=0; m<trial; m++){
			for(int i=0; i<ite; i++){
				for(int j=0; j<dim; j++){
					Xrec[i][j]=X1rec[i][j][m];
				}
			}
			String fnameX1 = "x" + (m+1);
			WriteResult.Output(Xrec, ite, dim, pathname + fnameX1 + sufi);
		}
		
		System.out.println(pathname + "test" + sufi); 
		
	}

}


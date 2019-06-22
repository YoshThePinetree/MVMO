import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

// Result writing class
public class WriteResult {
	// Method to display
		public static void Output(double X[][], int imax, int jmax, String filename){
			// write out the objective function value & search log
			try {
	            // create export
	            FileWriter fw = new FileWriter(filename, false); // enter your directory to create the file, true: add, false: overwrite
	            PrintWriter pw = new PrintWriter(new BufferedWriter(fw));

	            // the contents of the file
	            for(int i=0;i<imax;i++) {
					for(int j=0;j<jmax;j++) {
			            pw.print(X[i][j]);
			            if(j<jmax-1){
				            pw.print(",");
			            }else {
				            pw.println();		            	
			            }
					}
				}

	            // write it out
	            pw.close();

	        } catch (IOException ex) {
	            ex.printStackTrace();
	        }
		}

}

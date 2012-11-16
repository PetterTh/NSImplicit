
public class TDMAP {
	/**
	 * Solves the tri diagonal matrix equation using TDMA for problems with Periodic boundaries
	 */
	
	double[] x,bb,u,z;
	TDMA tdma;
	
	TDMAP (int n, double[] a, double[] b, double[] c, double[] v )
	{
	        /**
	         * n - number of equations
	         * a - sub-diagonal (means it is the diagonal below the main diagonal) -- indexed from 0..n-1
	         * b - the main diagonal
	         * c - sup-diagonal (means it is the diagonal above the main diagonal) -- indexed from 0..n-1
	         * v - right part
	         * x - the answer
	         */
		
        if (n<=2){
        	throw new IllegalArgumentException("n is too small for a cyclic problem");
        }
			x = new double[n];
			bb = new double[n];
			u = new double[n];
			z = new double[n];
			
	        double fact,gamma,beta,alpha;
	        beta=a[0];
	        alpha=c[n-1];
	        
	        gamma = -b[0]; 
	        
	        bb[0]=b[0]-gamma;
	        
	        bb[n-1]=b[n-1]-alpha*beta/gamma;
	        
	        for (int i = 1;i<n-1;i++){
	        	bb[i]=b[i];
	        }
	        
	        tdma = new TDMA(n,a,bb,c,v);
	        x=tdma.solve();
	        u[0]=gamma;
	        u[n-1]=alpha;
	        for (int i=1;i<n-1;i++) {
	        	u[i]=0.0;
	        }
	        tdma = new TDMA(n,a,bb,c,u);
	        z=tdma.solve();
	        fact=(x[0]+beta*x[n-1]/gamma)/(1.0+z[0]+beta*z[n-1]/gamma);
	        
	        for (int i=0;i<n;i++){
	        	x[i] -= fact*z[i];
	        }
	        
	        
	        
	}
	
	public double[] solve() {
		return x;
	}
	

}

import java.util.ArrayList;

public class NSSolverImplicit {

	final boolean periodic = true; // Used to tell if the boundaries should be
									// periodic
	final int nx, ny;
	final double dx, dy, V; // size of each cell
	double relax = 0.0; // Relaxation factor. Should not be used for transient
						// cases
	int nIt = 3; // Number of iterations with implicit solver
	boolean ExplicitSolver = false;
	boolean negativeDensity = false;

	private double[][] u; // velocity x
	private double[][] v; // velocity y
	private double[][] u_prev;
	private double[][] v_prev;
	private double[][] p; // pressure

	private double[][] e; // internal energy
	private double[][] e_prev;
	private double[][] r; // density
	private double[][] r_prev;
	private double[][] T; // Temperature in Kelvin
	private double[][] T_prev;

	private ArrayList<double[]> e_add = new ArrayList<double[]>(); // A list of
																	// containing
																	// points
																	// which
																	// should
																	// have its
																	// energy
																	// increased

	private boolean[][] Wall;

	private double my = 1.81e-5; // Viscosity

	public double Cv = 4.1796 * (100 * 100 * 100);// 0.001297 * (100 * 100 *
													// 100); // Heat Capacity at
													// constant
	// volume J/(cm3*K)
	public double R = 287.04; // Spesiffic gass constant Air

	public NSSolverImplicit(final int nx, final int ny, final double dx,
			final double dy) {

		this.nx = nx;
		this.ny = ny;
		this.dx = dx;
		this.dy = dy;
		this.V = dx * dy * 1;

		Cv = Cv * (dx * dy);
		zeroAll();
	}

	/*
	 * Used to copy a 2multidimentional array
	 */
	double[][] copy2DArray(double[][] src) {
		double[][] dest = new double[src.length][src[0].length];
		for (int i = 0; i < src.length; i++) {
			System.arraycopy(src[i], 0, dest[i], 0, src[0].length);
		}
		return dest;
	}

	private void zeroAll() {
		u = new double[nx][ny];
		v = new double[nx][ny];
		e = new double[nx][ny];
		r = new double[nx][ny];
		T = new double[nx][ny];
		p = new double[nx][ny];
		Wall = new boolean[nx][ny];

		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {

				Wall[i][j] = false;
				u[i][j] = 0;
				v[i][j] = 0;
				if (i < nx && j < ny) {
					p[i][j] = 101130;

					e[i][j] = 0;

					T[i][j] = 300;

					r[i][j] = p[i][j] / (R * T[i][j]);
				}
			}
		}
		u_prev = copy2DArray(u);
		v_prev = copy2DArray(v);

		e_prev = copy2DArray(e);
		T_prev = copy2DArray(T);
		r_prev = copy2DArray(r);
	}

	void addEnergy(int i, int j, double de) { // Add energy to the current cell

		double l[] = { i, j, de / r[i][j] };
		e_add.add(l);

	}

	void addPressure(int i, int j, double de) { // Add energy to the current
												// cell

		r_prev[i][j] += de / (R * T[i][j]);

	}

	double uS(int i, int j) {
		/*
		 * Calculate the velocity at the staggered grid position
		 */
		return (u[i][j] + u[i][j + 1] + u[i + 1][j] + u[i + 1][j + 1]) / 4;
	}

	double vS(int i, int j) {
		/*
		 * Calculate the velocity at the staggered grid position
		 */
		return (v[i][j] + v[i][j + 1] + v[i + 1][j] + v[i + 1][j + 1]) / 4;
	}

	private int mod(int x, int y) {
		int result = x % y;
		if (result < 0) {
			result += y;
		}
		return result;
	}

	void stepRho(double dt, boolean stepI) { // Calculate the new density
		// PApplet.arrayCopy(r, r_prev);
		TDMA tdma;
		int start = 1;
		if (periodic) {
			start = 0;
		}
		// for (int ni = 0; ni < nIt; ni++) {

		double[] a = new double[Math.max(nx - 1, ny - 1)];
		double[] b = new double[Math.max(nx - 1, ny - 1)];
		double[] c = new double[Math.max(nx - 1, ny - 1)];
		double[] s = new double[Math.max(nx - 1, ny - 1)];
		double[] x = new double[Math.max(nx - 1, ny - 1)];

		double uw, ue, us, un;
		double vw, ve, vs, vn;
		double rp;

		for (int i = 0; i < nx - start; i++) {
			int ineg = mod(i - 1, nx);
			int ipos = mod(i + 1, nx);

			for (int j = 0; j < ny - start; j++) {
				int jneg = mod(j - 1, ny);
				int jpos = mod(j + 1, ny);
				rp = r_prev[i][j];
				uw = (u[i][j] + u[i][jpos]) / 2;
				ue = (u[ipos][j] + u[ipos][jpos]) / 2;
				us = (u[i][jpos] + u[ipos][jpos]) / 2;
				un = (u[i][j] + u[ipos][j]) / 2;
				vw = (v[i][j] + v[i][jpos]) / 2;
				ve = (v[ipos][j] + v[ipos][jpos]) / 2;
				vs = (v[i][jpos] + v[ipos][jpos]) / 2;
				vn = (v[i][j] + v[ipos][j]) / 2;
				r[i][j] = r_prev[i][j] + rp * ((uw - ue) / dx + (vn - vs) / dy)
						* dt;
				if (r[i][j] < 0) {
					negativeDensity = true;
				}
			}
		}

	}

	// }

	void applyWall() { // Zero out the wall velocities
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (Wall[i][j]) {
					u_prev[i][j] = 0;
					v_prev[i][j] = 0;
					u[i][j] = 0;
					v[i][j] = 0;
				}
			}
		}
	}

	void stepVel(double dt, boolean calcU, boolean stepI) { // Calculate
															// the
															// new
		// velocities
		boolean calcV = !calcU;

		double rho;

		/*
		 * Calculate with TDMA - Implicit solver
		 */

		for (int ni = 0; ni < nIt; ni++) {

			stepT();
			stepRho(dt, stepI);
			stepP();
			stepE(dt);

			// applyWall();
			double[] a = new double[Math.max(nx, ny)];
			double[] b = new double[Math.max(nx, ny)];
			double[] c = new double[Math.max(nx, ny)];
			double[] s = new double[Math.max(nx, ny)];
			double[] x = new double[Math.max(nx, ny)];

			double Cn, Cs, Ce, Cw, C;

			int k = 1;
			int start = 1;
			if (periodic) {
				start = 0;
			}

			if (calcU) {

				// Calculate U TMDA in I dir
				if (stepI) {

					for (int j = start; j > (start - 1); j += k) {
						int jneg = mod(j - 1, ny);
						int jpos = mod(j + 1, ny);
						for (int i = start; i < nx - start; i++) {
							int ineg = mod(i - 1, nx);
							int ipos = mod(i + 1, nx);

							rho = (r[i][j] + r[i][jneg] + r[ineg][j] + r[ineg][jneg]) / 4;

							Cw = (-rho * u[i][j] / (2 * dx) - my / (dx * dx));
							C = (rho / (dt) + 2 * my / (dx * dx) + 2 * my
									/ (dy * dy));
							Ce = (rho * u[i][j] / (2 * dx) - my / (dx * dx));
							Cn = (-rho * v[i][j] / (2 * dy) - my / (dy * dy));
							Cs = (rho * v[i][j] / (2 * dy) - my / (dy * dy));

							a[i - start] = Ce;
							b[i - start] = C;
							c[i - start] = Cw;

							s[i - start] = rho
									* u_prev[i][j]
									/ dt
									- ((p[i][jneg] + p[i][j]) / 2 - (p[ineg][j] + p[ineg][jneg]) / 2)
									/ (dx)
									- (Cs * u[i][jpos] + Cn * u[i][jneg]);

						}
						if (periodic) {
							TDMAP tdma;
							tdma = new TDMAP(nx, a, b, c, s);
							x = tdma.solve();
						} else {
							TDMA tdma;
							tdma = new TDMA(nx - 2, a, b, c, s);
							x = tdma.solve();
						}

						for (int i = start; i < nx - start; i++) {
							u[i][j] = u[i][j] * relax + x[i - start]
									* (1 - relax);
						}
						if (j == ny - 1 - start) {
							k *= -1;
						}
					}
				} else {
					// Calculate U TMDA in J dir
					for (int i = start; i > (start - 1); i += k) {
						int ineg = mod(i - 1, nx);
						int ipos = mod(i + 1, nx);
						for (int j = start; j < ny - start; j++) {
							int jneg = mod(j - 1, ny);
							int jpos = mod(j + 1, ny);

							rho = (r[i][j] + r[i][jneg] + r[ineg][j] + r[ineg][jneg]) / 4;

							Cw = (-rho * u[i][j] / (2 * dx) - my / (dx * dx));
							C = (rho / dt + 2 * my / (dx * dx) + 2 * my
									/ (dy * dy));
							Ce = (rho * u[i][j] / (2 * dx) - my / (dx * dx));
							Cn = (-rho * v[i][j] / (2 * dy) - my / (dy * dy));
							Cs = (rho * v[i][j] / (2 * dy) - my / (dy * dy));

							a[j - start] = Cs;
							b[j - start] = C;
							c[j - start] = Cn;

							s[j - start] = rho
									* u_prev[i][j]
									/ dt
									- ((p[i][jneg] + p[i][j]) / 2 - (p[ineg][j] + p[ineg][jneg]) / 2)
									/ (dx)
									- (Ce * u[ipos][j] + Cw * u[ineg][j]);

						}

						if (periodic) {
							TDMAP tdma;
							tdma = new TDMAP(ny, a, b, c, s);
							x = tdma.solve();
						} else {
							TDMA tdma;
							tdma = new TDMA(ny - 2, a, b, c, s);
							x = tdma.solve();
						}
						for (int j = start; j < ny - start; j++) {
							u[i][j] = u[i][j] * relax + x[j - start]
									* (1 - relax);
						}

						if (i == nx - 1 - start) {
							k *= -1;
						}
					}
				}
			}
			if (calcV) {

				a = new double[Math.max(nx, ny)];
				b = new double[Math.max(nx, ny)];
				c = new double[Math.max(nx, ny)];
				s = new double[Math.max(nx, ny)];
				x = new double[Math.max(nx, ny)];

				// Calculate V TMDA in J dir
				if (!stepI) {
					for (int i = start; i > (start - 1); i += k) {
						int ineg = mod(i - 1, nx);
						int ipos = mod(i + 1, nx);

						for (int j = start; j < ny - start; j++) {
							int jneg = mod(j - 1, ny);
							int jpos = mod(j + 1, ny);

							rho = (r[i][j] + r[i][jneg] + r[ineg][j] + r[ineg][jneg]) / 4;

							Cw = (-rho * u[i][j] / (2 * dx) - my / (dx * dx));
							C = (rho / dt + 2 * my / (dx * dx) + 2 * my
									/ (dy * dy));
							Ce = (rho * u[i][j] / (2 * dx) - my / (dx * dx));
							Cn = (-rho * v[i][j] / (2 * dy) - my / (dy * dy));
							Cs = (rho * v[i][j] / (2 * dy) - my / (dy * dy));

							a[j - start] = Cs;
							b[j - start] = C;
							c[j - start] = Cn;

							s[j - start] = rho
									* v_prev[i][j]
									/ dt
									- ((p[i][j] + p[ineg][j]) / 2 - (p[i][jneg] + p[ineg][jneg]) / 2)
									/ (dy)
									- (Ce * v[ipos][j] + Cw * v[ineg][j]);

						}

						if (periodic) {
							TDMAP tdma;
							tdma = new TDMAP(ny, a, b, c, s);
							x = tdma.solve();
						} else {
							TDMA tdma;
							tdma = new TDMA(ny - 2, a, b, c, s);
							x = tdma.solve();
						}
						for (int j = start; j < ny - start; j++) {
							v[i][j] = v[i][j] * relax + x[j - start]
									* (1 - relax);

						}

						if (i == nx - 1 - start) {
							k *= -1;
						}
					}
				} else {
					// Calculate V TMDA in I dir
					for (int j = start; j > (start - 1); j += k) {
						// Positive j direction
						int jneg = mod(j - 1, ny);
						int jpos = mod(j + 1, ny);

						for (int i = start; i < nx - start; i++) {
							int ineg = mod(i - 1, nx);
							int ipos = mod(i + 1, nx);

							rho = (r[i][j] + r[i][jneg] + r[ineg][j] + r[ineg][jneg]) / 4;

							Cw = (-rho * u[i][j] / (2 * dx) - my / (dx * dx));
							C = (rho / dt + 2 * my / (dx * dx) + 2 * my
									/ (dy * dy));
							Ce = (rho * u[i][j] / (2 * dx) - my / (dx * dx));
							Cn = (-rho * v[i][j] / (2 * dy) - my / (dy * dy));
							Cs = (rho * v[i][j] / (2 * dy) - my / (dy * dy));

							a[i - start] = Ce;
							b[i - start] = C;
							c[i - start] = Cw;

							s[i - start] = rho
									* v_prev[i][j]
									/ dt
									- ((p[i][j] + p[ineg][j]) / 2 - (p[i][jneg] + p[ineg][jneg]) / 2)
									/ (dy)
									- (Cs * v[i][jpos] + Cn * v[i][jneg]);

						}

						if (periodic) {
							TDMAP tdma;
							tdma = new TDMAP(nx, a, b, c, s);
							x = tdma.solve();
						} else {
							TDMA tdma;
							tdma = new TDMA(nx - 2, a, b, c, s);
							x = tdma.solve();
						}
						for (int i = start; i < nx - start; i++) {
							v[i][j] = v[i][j] * relax + x[i - start]
									* (1 - relax);

						}

						if (j == ny - 1 - start) {
							k *= -1;
						}
					}
				}
			}

		}

	}

	void stepE(double dt) { // Calculate the new energy

		int start = 1;
		if (periodic) {
			start = 0;
		}
		double uw, ue, us, un;
		double vw, ve, vs, vn;
		double rp;

		for (int i = 0; i < nx - start; i++) {
			int ineg = mod(i - 1, nx);
			int ipos = mod(i + 1, nx);

			for (int j = 0; j < ny - start; j++) {
				int jneg = mod(j - 1, ny);
				int jpos = mod(j + 1, ny);
				rp = r_prev[i][j];
				uw = (u[i][j] + u[i][jpos]) / 2;
				ue = (u[ipos][j] + u[ipos][jpos]) / 2;
				us = (u[i][jpos] + u[ipos][jpos]) / 2;
				un = (u[i][j] + u[ipos][j]) / 2;
				vw = (v[i][j] + v[i][jpos]) / 2;
				ve = (v[ipos][j] + v[ipos][jpos]) / 2;
				vs = (v[i][jpos] + v[ipos][jpos]) / 2;
				vn = (v[i][j] + v[ipos][j]) / 2;
				e[i][j] = e_prev[i][j]
						+ ((uw * (e_prev[ineg][j] + p[ineg][j]) - ue
								* (e_prev[ipos][j] + p[ipos][j]))
								/ dx + (vn * (e_prev[i][jneg] + p[i][jneg]) - vs
								* (e_prev[i][jpos] + p[i][jpos]))
								/ dy) * dt;

			}
		}
		double[] item;
		for (int i = 0; i < e_add.size(); i++) {
			item = e_add.get(i);
			int ii = (int) item[0];
			int jj = (int) item[1];
			double de = (double) item[2] / 4;
			e[ii][jj] += de;
			// Add 1/4th of the energy as this is done four times
			// see step()

		}
	}

	void stepT() { // Calculate new Temperatuers
		double e1, e2;
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				e1 = (e[i][j]) / r[i][j] - 0.5
						* (Math.pow(u[i][j], 2) + Math.pow(v[i][j], 2));
				e2 = (e_prev[i][j])
						/ r_prev[i][j]
						- 0.5
						* (Math.pow(u_prev[i][j], 2) + Math
								.pow(v_prev[i][j], 2));
				T[i][j] += (e1 - e2) / (Cv);

			}
		}
	}

	void stepP() { // Calculate new Pressure from ideal gass law

		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				p[i][j] = p[i][j] * (relax) + r[i][j] * R * T[i][j]
						* (1 - relax);

			}
		}
	}

	public void step(double dt) { // Step the function one time step, dt

		stepVel(dt / 4, true, true);

		u_prev = copy2DArray(u);
		v_prev = copy2DArray(v);
		e_prev = copy2DArray(e);
		T_prev = copy2DArray(T);
		r_prev = copy2DArray(r);

		stepVel(dt / 4, true, false);

		u_prev = copy2DArray(u);
		v_prev = copy2DArray(v);
		e_prev = copy2DArray(e);
		T_prev = copy2DArray(T);
		r_prev = copy2DArray(r);
		stepVel(dt / 4, false, true);

		u_prev = copy2DArray(u);
		v_prev = copy2DArray(v);
		e_prev = copy2DArray(e);
		T_prev = copy2DArray(T);
		r_prev = copy2DArray(r);

		stepVel(dt / 4, false, false);

		u_prev = copy2DArray(u);
		v_prev = copy2DArray(v);
		e_prev = copy2DArray(e);
		T_prev = copy2DArray(T);
		r_prev = copy2DArray(r);

		// Reset the added energy
		e_add = new ArrayList<double[]>();

	}

	public double getP(int i, int j) {
		return p[i][j];
	}

	public double getE(int i, int j) {
		return e[i][j];
	}

	public double getT(int i, int j) {
		return T[i][j];
	}

	public double getU(int i, int j) {
		return u[i][j];
	}

	public double getV(int i, int j) {
		return v[i][j];
	}

	public double getR(int i, int j) {
		return r[i][j];
	}

	public void setWall(int i, int j) {
		Wall[i][j] = true;
	}

	public void setVel(int i, int j, double vel, boolean udir) {
		if (udir) {

			u[i][j] = vel;

			u_prev[i][j] = vel;

		} else {
			v_prev[i][j] = vel;

			v[i][j] = vel;
		}
	}

	/*
	 * Create a refined grid of the velocities with linear interpolation between
	 * each point
	 */
	public double[][] refineU(int nnx, int nny) {

		double[][] r = new double[nnx][nny];

		double dx = (double) nx / nnx;
		double dy = (double) ny / nny;

		double une, unw, usw, use; //NorthEast,NorthWest,SouthWest,SouthEast
		int ii, jj;

		for (int i = 0; i < nnx; i++) {
			for (int j = 0; j < nny; j++) {
				ii = (int) Math.floor(i * dx);
				jj = (int) Math.floor(j * dy);
				unw = u[ii][jj];
				une = u[ii + 1][jj];
				usw = u[ii][jj + 1];
				use = u[ii + 1][jj + 1];

				r[i][j] = (unw * (1 - (i * dx-ii)) + une * (i * dx-ii))
						* (1 - (j * dy-ii))
						+ (usw * (1 - (i * dx-ii)) + use * (i * dx-ii))
						* ((j * dy-jj));
			}
		}

		return r;

	}

	public double[][] refineV(int nnx, int nny) {

		double[][] r = new double[nnx][nny];

		double dx = (double) nx / nnx;
		double dy = (double) ny / nny;

		double vne, vnw, vsw, vse;
		int ii, jj;

		for (int i = 0; i < nnx; i++) {
			for (int j = 0; j < nny; j++) {
				ii = (int) Math.floor(i * dx);
				jj = (int) Math.floor(j * dy);
				vnw = v[ii][jj];
				vne = v[ii + 1][jj];
				vsw = v[ii][jj + 1];
				vse = v[ii + 1][jj + 1];

				r[i][j] = (vnw * (1 - (i * dx-ii)) + vne * (i * dx-ii))
						* (1 - (j * dy-ii))
						+ (vsw * (1 - (i * dx-ii)) + vse * (i * dx-ii))
						* ((j * dy-jj));
			}
		}

		return r;

	}

}

public class Gas {
	private double lambda, alpha, beta, nu, m;
	private double delvol, lvol, l;
	private double Vcr, Pcr, Tcr;
	private double akr, bkr, ckr, dkr, ekr, fkr, betavol;
	private double delta, eps, S;
	private double x0, y0, akrsh, ckrsh;
	private double g;
	private double tappr;

	public Gas(double alpha, double beta, double nu, double m) throws Exception {
		this.alpha = alpha;
		this.beta = beta;
		this.nu = nu;
		this.m = m;

		delvol = 3 * beta * nu / (1 - 2 * nu);
		l = 8. / 27.;
		Tcr = lambda = lvol = Math.pow(alpha * l / (beta + delvol), 1. / (m + 1.));

		Vcr = 3. * beta / (1. - 2. * nu);
		Pcr = lvol * (1. - 2. * nu) / (8. * beta * (1. + nu));

		recalcConst();

		calcTappr();
	}

	public void changeLambda(double newlambda) throws Exception {
		lambda = newlambda;
		recalcConst();
	}

	private void recalcConst() throws Exception {
		double lamp = Math.pow(lambda, m + 1.);

		akr = lamp;
		bkr = delvol * lamp - alpha / 2.;
		ckr = lamp * delvol * delvol + alpha * beta;
		dkr = -alpha * delvol / 2.;
		ekr = -alpha * beta * beta / 2.;
		betavol = beta * (beta + delvol);
		fkr = -alpha * delvol * betavol;

		delta = -lamp * alpha * alpha * (beta + delvol) * (beta + delvol) * (beta + delvol) * (beta + delvol) / 4.;
		eps = alpha * ((beta + delvol) * lamp - alpha / 4.);
		S = lamp * (1 + delvol * delvol) + alpha * beta;

		if (delta > 0)
			throw new Exception("Delta > 0. Error!");

		x0 = (bkr * ekr - ckr * dkr) / eps; // TODO use raw parameters
		y0 = (dkr * bkr - ekr * akr) / eps; // TODO use raw parameters

		akrsh = (S + Math.signum(akr - ckr) * Math.sqrt(S * S - 4. * eps)) / 2.;
		ckrsh = (S - Math.signum(akr - ckr) * Math.sqrt(S * S - 4. * eps)) / 2.;

		if (akr - ckr > 0) {
			g = (akr - ckr - Math.sqrt((akr - ckr) * (akr - ckr) + 4 * bkr * bkr)) / (2. * bkr);
		} else {
			g = (Math.sqrt((akr - ckr) * (akr - ckr) + 4 * bkr * bkr) - Math.abs(akr - ckr)) / (2. * bkr);
		}
	}

	private void calcTappr() throws Exception {
		if (lambda != Tcr) {
			throw new Exception("calcTappr should be called only if lambda == Tcr");
		}

		double xshcr = Vcr * Vcr - x0;
		double yshcr = 2 * Vcr - y0 + delvol;

		double k = Math.sqrt(Math.abs(delta) / eps) / (Math.sqrt(1 + g * g));

		double a1 = 2 * g * k / Math.sqrt(ckrsh);
		double a2 = Math.sqrt(4 * g * g * k * k / Math.abs(ckrsh) - 4 * (xshcr * xshcr - k * k / Math.abs(akrsh)));
		double a3 = 2 * (xshcr + k / Math.sqrt(akrsh));

		double t1 = (a1 + a2) / a3;
		double t2 = (a1 - a2) / a3;

		if (Math.abs(y(t1) - 2 * Vcr) < Math.abs(y(t2) - 2 * Vcr)) {
			tappr = t1;
		} else {
			tappr = t2;
		}
	}

	/**
	 * 
	 * @param t
	 *            - in case eps>0 means "t", in case eps<0 means "z"
	 * @return
	 * @throws Exception
	 */
	public double x(double t) throws Exception {
		double xsh;
		if (eps > 0) {
			xsh = Math.sqrt(Math.abs(delta) / eps) * (1. / Math.sqrt(1 + g * g))
					* ((1. - t * t) / Math.sqrt(akrsh) + 2. * t * g / Math.sqrt(ckrsh)) / (1. + t * t);
		} else {
			double fi = Math.atan(2 * bkr / (ckr - akr)) / 2.;
			double xshsh, yshsh;
			if (akr - ckr > 0) {
				xshsh = Math.sqrt(delta / (eps * akrsh)) * Math.sinh(t);
				yshsh = Math.sqrt(delta / (eps * Math.abs(ckrsh))) * Math.cosh(t);
			} else {
				xshsh = Math.sqrt(delta / (eps * Math.abs(akrsh))) * Math.cosh(t);
				yshsh = Math.sqrt(delta / (eps * ckrsh)) * Math.sinh(t);
			}
			xsh = xshsh * Math.cos(fi) + yshsh * Math.sin(fi);
		}
		return x0 + xsh;
	}

	/**
	 * 
	 * @param t
	 *            - in case eps>0 means "t", in case eps<0 means "z"
	 * @return
	 * @throws Exception
	 */
	public double y(double t) throws Exception {
		double ysh;
		if (eps > 0) {
			ysh = Math.sqrt(Math.abs(delta) / eps) * (1. / Math.sqrt(1 + g * g))
					* (-g * (1. - t * t) / Math.sqrt(akrsh) + 2. * t / Math.sqrt(ckrsh)) / (1. + t * t);
		} else {
			double fi = Math.atan(2 * bkr / (ckr - akr)) / 2.;
			double xshsh, yshsh;
			if (akr - ckr > 0) {
				xshsh = Math.sqrt(delta / (eps * akrsh)) * Math.sinh(t);
				yshsh = Math.sqrt(delta / (eps * Math.abs(ckrsh))) * Math.cosh(t);
			} else {
				xshsh = Math.sqrt(delta / (eps * Math.abs(akrsh))) * Math.cosh(t);
				yshsh = Math.sqrt(delta / (eps * ckrsh)) * Math.sinh(t);
			}
			ysh = -xshsh * Math.sin(fi) + yshsh * Math.cos(fi);
		}
		return ysh + y0 - delvol;
	}

	public double getP(double V) {
		return lambda / (V - beta) - alpha / (Math.pow(lambda, m) * (V + delvol) * (V + delvol));
	}

	/**
	 * Object-safe way to calculate some values.
	 * 
	 * @param V
	 * @param T
	 * @return
	 */
	public double getP(double V, double T) {
		return T / (V - beta) - alpha / (Math.pow(T, m) * (V + delvol) * (V + delvol));
	}

	/**
	 * 
	 * @param t
	 *            - in case eps>0 means "t", in case eps<0 means "z"
	 * @return
	 * @throws Exception
	 */
	public double Q(double t) throws Exception {
		double x = x(t);
		double y = y(t);

		double q = 4 * x / y / y;

		double a = 1. / (y * (1 - Math.sqrt(1 - q)) - 2 * beta);
		double b = 1. / (2 * y * Math.sqrt(1 - q));
		double c = (y * (1 - Math.sqrt(1 - q)) - 2 * beta);
		double d = (y * (1 + Math.sqrt(1 - q)) - 2 * beta);
		double e = alpha * y * Math.sqrt(1 - q)
				/ ((y * (1 - Math.sqrt(1 - q)) + 2 * delvol) * (x + delvol * (y - delvol)));

		return Math.pow(lambda, m + 1.) * (a + b * Math.log(c / d)) - e;
	}

	public double getTappr() {
		return tappr;
	}

	public double getTcr() {
		return Tcr;
	}

	public double getPcr() {
		return Pcr;
	}

	public double getVcr() {
		return Vcr;
	}

	public double getAlpha() {
		return alpha;
	}

	public double getBeta() {
		return beta;
	}

	public double getNu() {
		return nu;
	}

	public double getM() {
		return m;
	}

	public double getDelvol() {
		return delvol;
	}

	public double getEps() {
		return eps;
	}
}

import java.awt.Color;
import java.awt.Font;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYPointerAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.ui.TextAnchor;

public class Main {
	public static void main(String[] args) throws Exception {

		Gas nitrogen = new Gas(0.5, 0.9, 0, 0.36);
		handleGas(nitrogen);
	}

	private static void handleGas(Gas gas) throws Exception {
		double tbegin = 0.2;
		double tend = gas.getTcr();
		int N = 5000;
		int tN = 7;
		double vbegin = 0.1;
		double vend = 100;

		double[] vaxis = new double[N];
		double[][] funcGraph = new double[tN][N];
		for (int q = 0; q < N; q++) {
			vaxis[q] = vbegin + (vend - vbegin) / N * q;
			// fcrGraph[q] = f(vaxis[q], Tcr, d, l, a, b);
			for (int w = 0; w < tN; w++) {
				double t = tbegin + (tend - tbegin) / (tN - 1) * w;
				funcGraph[w][q] = gas.getP(vaxis[q], t);
			}
		}

		DefaultXYDataset data;
		ChartFrame frame;
		JFreeChart chart;

		data = new DefaultXYDataset();
		chart = ChartFactory.createXYLineChart(" ", "V", "P", data);
		chart.getLegend().setVisible(false);
		XYPlot plot = (XYPlot) chart.getPlot();
		XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();

		data.addSeries("X axis", new double[][] { vaxis, new double[N] });
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);

		data.addSeries("Critical point", new double[][] { { gas.getVcr() }, { gas.getPcr() } });
		renderer.setSeriesLinesVisible(data.getSeriesCount() - 1, false);
		renderer.setSeriesShapesVisible(data.getSeriesCount() - 1, true);
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);

		for (int w = 0; w < tN; w++) {
			double t = tbegin + (tend - tbegin) / (tN - 1) * w;
			data.addSeries("" + t, new double[][] { vaxis, funcGraph[w] });
			renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);
		}

		final XYPointerAnnotation pointer = new XYPointerAnnotation(
				"Critical point. V = " + gas.getVcr() + "; P = " + gas.getPcr() + "; T = " + gas.getTcr() + ".",
				gas.getVcr(), gas.getPcr(), -Math.PI / 4.0);
		pointer.setBaseRadius(45.0);
		pointer.setTipRadius(10.0);
		pointer.setFont(new Font("SansSerif", Font.PLAIN, 12));
		pointer.setPaint(Color.black);
		pointer.setTextAnchor(TextAnchor.HALF_ASCENT_LEFT);
		plot.addAnnotation(pointer);

		// marker.setLabel("Delta = "+d+"; Lamda = "+l+"; Alpha = "+a+"; Beta =
		// "+b+".");
		final XYTextAnnotation annotation = new XYTextAnnotation("Alpha = " + gas.getAlpha() + "; Beta = "
				+ gas.getBeta() + "; Nu = " + gas.getNu() + "; M = " + gas.getM() + ".", 12.5, -0.005);
		annotation.setFont(new Font("SansSerif", Font.PLAIN, 12));
		plot.addAnnotation(annotation);

		zone(data, renderer, gas);

		frame = createFrame(chart, "Первый случай", 0, 25, -0.01, 0.06);
	}

	private static void zone(DefaultXYDataset data, XYLineAndShapeRenderer renderer, Gas gas) throws Exception {
		int N = 5000;
		double h = gas.getTcr() / (N + N / 10);

		double lambda = gas.getTcr();
		gas.changeLambda(lambda);

		double y = 2 * gas.getVcr();
		double x = gas.getVcr() * gas.getVcr();
		double t = gas.getTappr();
		double bineps = 0.01;

		double[][] arrVbVc = new double[2][2 * N];
		for (int i = 0; i < N; i++) {

			lambda -= h;
			double preveps = gas.getEps();
			gas.changeLambda(lambda);

			double newt;
			if (gas.getEps() < 0) {
				// break;

				bineps = 2;
				try{
					newt = binhQSearch(t - bineps, t + bineps, gas);					
				}
				catch(Exception ex)
				{
					//TODO 
					double[][] qqq = new double[2][N];
					for(int q=0;q<N;q++)
					{
						qqq[0][q] = t - bineps + 2*q*bineps/N;
						qqq[1][q] = gas.Q(Math.tanh(qqq[0][q]/2.));
					}

					data.addSeries("Q", qqq);
					renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.red);
					break;
				}
				
//				newt = binQSearch(Math.tanh((t - bineps) / 2.), Math.tanh((t + bineps) / 2.), gas);
				System.out.println("" + newt);
				y = gas.y(Math.tanh(newt));
				x = gas.x(Math.tanh(newt));

//				newt = Math.log((1 + newt) / (1 - newt)); // 2 * atanh
			} else {
				newt = binQSearch(t - bineps, t + bineps, gas);
				y = gas.y(newt);
				x = gas.x(newt);
			}
			bineps = Math.abs(4 * (newt - t));
			t = newt;

			double q = 4 * x / (y * y);
			double Vc = y * (1 + Math.sqrt(1 - q)) / 2.;
			double Vb = y * (1 - Math.sqrt(1 - q)) / 2.;

			arrVbVc[0][N - i - 1] = Vb;
			arrVbVc[1][N - i - 1] = gas.getP(Vb, lambda);
			arrVbVc[0][N + i] = Vc;
			arrVbVc[1][N + i] = gas.getP(Vc, lambda);
		}
		data.addSeries("Critical zone", arrVbVc);
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.red);
	}

	private static double binQSearch(double t1, double t2, Gas gas) throws Exception {
		double t = (t1 + t2) / 2.;
//		if (Math.abs(t2 - t1) < 1e-8)
//			return t;

		double Q1 = gas.Q(t1);
		double Q = gas.Q(t);
		double Q2 = gas.Q(t2);

		if (Math.abs(Q) < 1e-12)
			return t;

		if (Q * Q1 < 0)
			return binQSearch(t1, t, gas);
		if (Q * Q2 < 0)
			return binQSearch(t, t2, gas);

//		if ((Q - Q1) * (Q2 - Q) > 0) {
//			double tn2 = Math.abs(Q1) < Math.abs(Q2) ? t1 : t2;
//			double tn1 = Math.abs(Q1) < Math.abs(Q2) ? t1 - (t2 - t1) : t2 - (t1 - t2);
//			return binQSearch(tn1, tn2, gas);
//		}

		throw new Exception("Binsearch error! t1 = " + t1 + "; Q1 = " + Q1 + "; t2 = " + t2 + "; Q2 = " + Q2 + "; t = "
				+ t + "; Q = " + Q);
	}
	private static double binhQSearch(double t1, double t2, Gas gas) throws Exception {
		double t = (t1 + t2) / 2.;
		if (Math.abs(t2 - t1) < 1e-8)
			return t;

		double Q1 = gas.Q(Math.tanh(t1/2.));
		double Q = gas.Q(Math.tanh(t/2.));
		double Q2 = gas.Q(Math.tanh(t2/2.));

		if (Math.abs(Q) < 1e-12)
			return t;

		if (Q * Q1 < 0)
			return binhQSearch(t1, t, gas);
		if (Q * Q2 < 0)
			return binhQSearch(t, t2, gas);

		if ((Q - Q1) * (Q2 - Q) > 0) {
			double tn2 = Math.abs(Q1) < Math.abs(Q2) ? t1 : t2;
			double tn1 = Math.abs(Q1) < Math.abs(Q2) ? t1 - (t2 - t1)/2. : t2 - (t1 - t2)/2.;
			return binhQSearch(tn1, tn2, gas);
		}

		throw new Exception("Binsearch error! t1 = " + t1 + "; Q1 = " + Q1 + "; t2 = " + t2 + "; Q2 = " + Q2 + "; t = "
				+ t + "; Q = " + Q);
	}

	private static ChartFrame createFrame(JFreeChart chart, String name, double xb, double xe, double yb, double ye) {
		XYPlot xyplot = (XYPlot) chart.getPlot();
		xyplot.setBackgroundPaint(Color.white);
		xyplot.setRangeGridlinePaint(Color.gray);
		xyplot.setDomainGridlinePaint(Color.gray);
		xyplot.getRangeAxis().setRange(yb, ye);
		xyplot.getDomainAxis().setRange(xb, xe);

		XYLineAndShapeRenderer lineandshaperenderer = (XYLineAndShapeRenderer) xyplot.getRenderer();
		lineandshaperenderer.setDrawOutlines(true);
		lineandshaperenderer.setUseFillPaint(true);
		lineandshaperenderer.setBaseFillPaint(Color.white);

		ChartFrame frame = new ChartFrame(name, chart);

		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
		return frame;
	}
}

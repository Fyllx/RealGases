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

public class BeattieBridgeman
{
	// static double vbegin = 0.1;
	// static double vend = 100;
	// static int N = 20000;
	//
	// static double tbegin;
	// static double tend;
	// static int tN = 30;// (int) (tend - tbegin);

	public static double f(double v, double t, double d, double l, double a, double b)
	{
		return (t - d / (v * v * v)) / v - (1 - b / v) * a / (v * v);
	}

	public static void main(String[] args)
	{
		firstCase();
		secondCase();
	}

	private static double P(double x, double y, double T, double a, double b, double d)
	{
		return T * x * x * x - a * b * (y * y - x) * x - d * (y * y * y - 2 * x * y);
	}

	private static double Q(double x, double y, double T, double a, double b, double d)
	{
		double q = 4 * x / (y * y);
		double Vc = y / 2 * (1 + Math.sqrt(1 - q));
		double Vb = y / 2 * (1 - Math.sqrt(1 - q));
		return (T * (1 + Vb * Math.log(Vb / Vc) / (Vc - Vb)) - (Vc - Vb) * (a / x + d / 3 * (2*Vc*Vc+y*y) / x / x / x - a * b / 2 * (Vc + y) / x / x));
	}

	private static double getX(double x, double y, double T, double a, double b, double d)
	{
		double eps = 10e-15;
		double t = x / y;
		double t1 = t;
		do
		{
			t = t1;
			double L = T * t * t * t + a * b * t - d;
			t1 = (d + Math.sqrt(d * d + a * (y + b) * y * L)) / (a * (y + b));
		} while (Math.abs(t1 - t) >= eps);
		return t1 * y;
	}
	private static double[] binarySearch(double y1, double x1, double Q1, double y2, double x2, double Q2, double T, double a, double b, double d)
	{
		double eps = 10e-12;
		double y = (y1+y2)/2;
		double x = getX(x1, y, T, a, b, d);
		double Q = Q(x, y, T, a, b, d);
		
		if(Math.abs(Q) < eps) return new double[] {x, y};
		
		if(Q1*Q < 0)
		{
			return binarySearch(y1, x1, Q1, y, x, Q, T, a, b, d);
		}
		else
		{
			if(Q*Q2 < 0)
			{
				return binarySearch(y, x, Q, y2, x2, Q2, T, a, b, d);
			}
		}
		System.out.println("������ � �������� ������:");
		System.out.println("T = "+T+";");
		System.out.println("y1 = "+y1+"; x1 = "+x1+"; Q1 = "+Q1);
		System.out.println("y2 = "+y2+"; x2 = "+x2+"; Q2 = "+Q2);
		System.out.println("_____________________________________________________");
		return null;
	}

	private static void firstZone(DefaultXYDataset data, XYLineAndShapeRenderer renderer, double Vcr, double Tcr, double a, double b, double d, double l)
	{
		int N = 5000;
		double h = Tcr / (N+N/10);

		double T = Tcr;
		double y = 2 * Vcr;
		double delta = 0.1;
		double x = Vcr * Vcr;

		double[][] arrVbVc = new double[2][2*N];
		for (int q = 0; q < N; q++)
		{
			T -= h;
			double y1 = y+10e-15;
			double x1 = getX(x, y1, T, a, b, d);
			double Q1 = Q(x1, y1, T, a, b, d);
			double y2 = y+delta;
			double x2 = getX(x, y2, T, a, b, d);
			double Q2 = Q(x2, y2, T, a, b, d);
			
			double[] buff = binarySearch(y1, x1, Q1, y2, x2, Q2, T, a, b, d);
			if(buff == null) break;
			x = buff[0];
			double yn = buff[1];
			delta = Math.abs(4*(yn-y));
			y = yn;
			
			double qq = 4 * x / (y * y);
			double Vc = y / 2 * (1 + Math.sqrt(1 - qq));
			double Vb = y / 2 * (1 - Math.sqrt(1 - qq));
			
			arrVbVc[0][N-q-1] = Vb;
			arrVbVc[1][N-q-1] = f(Vb, T, d, l, a, b);
			arrVbVc[0][N+q] = Vc;
			arrVbVc[1][N+q] = f(Vc, T, d, l, a, b);
		}
		data.addSeries("���� ����������� ���������", arrVbVc);
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.red);
	}
	private static void secondZone(DefaultXYDataset data, XYLineAndShapeRenderer renderer, double Vcr, double Tcr, double a, double b, double d, double l)
	{
		int N = 5000;
		double h = Tcr / (N+N/10);

		double T = Tcr;
		double y = 2 * Vcr;
		double delta = 0.1;
		double x = Vcr * Vcr;

		double[][] arrVbVc = new double[2][2*N];
		for (int q = 0; q < N; q++)
		{
			T -= h;
			double y1 = y+10e-15;
			double x1 = getX(x, y1, T, a, b, d);
			double Q1 = Q(x1, y1, T, a, b, d);
			double y2 = y+delta;
			double x2 = getX(x, y2, T, a, b, d);
			double Q2 = Q(x2, y2, T, a, b, d);
			
			double[] buff = binarySearch(y1, x1, Q1, y2, x2, Q2, T, a, b, d);
			if(buff == null) break;
			x = buff[0];
			double yn = buff[1];
			delta = Math.abs(4*(yn-y));
			y = yn;
			
			double qq = 4 * x / (y * y);
			double Vc = y / 2 * (1 + Math.sqrt(1 - qq));
			double Vb = y / 2 * (1 - Math.sqrt(1 - qq));
			
			arrVbVc[0][N-q-1] = Vb;
			arrVbVc[1][N-q-1] = f(Vb, T, d, l, a, b);
			arrVbVc[0][N+q] = Vc;
			arrVbVc[1][N+q] = f(Vc, T, d, l, a, b);
		}
		data.addSeries("���� ����������� ���������", arrVbVc);
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.red);
	}

	private static void firstCase()
	{
		double l = 0;
		double a = 1;
		double b = 0.5;
		double d = -1;

		double K = 1 + Math.sqrt(1 - 8 * d / (3 * a * b * b));

		double Vcr = 3 * b * K / 2;

		double Tcr = (3 * a * b - 8 * d / Vcr) / (Vcr * Vcr);

		double tbegin = 0.25;
		double tend = Tcr;
		int N = 20000;
		int tN = 7;
		double vbegin = 0.1;
		double vend = 100;

		double Pcr = f(Vcr, Tcr, d, l, a, b);

		double[] vaxis = new double[N];
		// double[] fcrGraph = new double[N];
		double[][] funcGraph = new double[tN][N];

		for (int q = 0; q < N; q++)
		{
			vaxis[q] = vbegin + (vend - vbegin) / N * q;
			// fcrGraph[q] = f(vaxis[q], Tcr, d, l, a, b);
			for (int w = 0; w < tN; w++)
			{
				double t = tbegin + (tend - tbegin) / (tN - 1) * w;
				funcGraph[w][q] = f(vaxis[q], t, d, l, a, b);
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

		data.addSeries("��� �", new double[][] { vaxis, new double[N] });
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);

		data.addSeries("����������� �����", new double[][] { { Vcr }, { Pcr } });
		renderer.setSeriesLinesVisible(data.getSeriesCount() - 1, false);
		renderer.setSeriesShapesVisible(data.getSeriesCount() - 1, true);
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);

		// data.addSeries("������ ��� �����������(-) ����������� "+Tcrm, new
		// double[][] { vaxis, fcrmGraph });
		for (int w = 0; w < tN; w++)
		{
			double t = tbegin + (tend - tbegin) / (tN - 1) * w;
			data.addSeries("" + t, new double[][] { vaxis, funcGraph[w] });
			renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);
		}

		final XYPointerAnnotation pointer = new XYPointerAnnotation("Критическая точка. V = " + Vcr + "; P = " + Pcr + "; T = " + Tcr + ".", Vcr, Pcr,
				-Math.PI / 4.0);
		pointer.setBaseRadius(45.0);
		pointer.setTipRadius(10.0);
		pointer.setFont(new Font("SansSerif", Font.PLAIN, 12));
		pointer.setPaint(Color.black);
		pointer.setTextAnchor(TextAnchor.HALF_ASCENT_LEFT);
		plot.addAnnotation(pointer);

		// marker.setLabel("Delta = "+d+"; Lamda = "+l+"; Alpha = "+a+"; Beta = "+b+".");
		final XYTextAnnotation annotation = new XYTextAnnotation("Delta = " + d + "; Lamda = " + l + "; Alpha = " + a + "; Beta = " + b + ".", 12.5, -0.005);
		annotation.setFont(new Font("SansSerif", Font.PLAIN, 12));
		plot.addAnnotation(annotation);

		firstZone(data, renderer, Vcr, Tcr, a, b, d, l);

		frame = createFrame(chart, "Первый случай", 0, 25, -0.01, 0.06);
	}

	private static void secondCase()
	{
		double l = 0;
		double a = 8;
		double b = 0.5;
		double d = 0.5;
		
//		double l = 0;
//		double a = 1;
//		double b = 0.5;
//		double d = a * b * b * 108. / 320.;

		double Kp = 1 + Math.sqrt(1 - 8 * d / (3 * a * b * b));
		double Km = 1 - Math.sqrt(1 - 8 * d / (3 * a * b * b));

		double Vcrp = 3 * b * Kp / 2;
		double Vcrm = 3 * b * Km / 2;

		double Tcrp = (3 * a * b - 8 * d / Vcrp) / (Vcrp * Vcrp);
		double Tcrm = (3 * a * b - 8 * d / Vcrm) / (Vcrm * Vcrm);

		double tbegin = Tcrm;
		double tend = Tcrp;
		int N = 20000;
		int tN = 14;
		double vbegin = 0.1;
		double vend = 100;

		double Pcrp = f(Vcrp, Tcrp, d, l, a, b);
		double Pcrm = f(Vcrm, Tcrm, d, l, a, b);

		double[] vaxis = new double[N];
		// double[] fcrmGraph = new double[N];
		double[][] funcGraph = new double[tN][N];

		for (int q = 0; q < N; q++)
		{
			vaxis[q] = vbegin + (vend - vbegin) / N * q;
			// fcrmGraph[q] = f(vaxis[q], Tcrm, d, l, a, b);
			for (int w = 0; w < tN; w++)
			{
				double t = tbegin + (tend - tbegin) / (tN - 1) * w;
				funcGraph[w][q] = f(vaxis[q], t, d, l, a, b);
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

		data.addSeries("��� �", new double[][] { vaxis, new double[N] });
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);

		data.addSeries("����������� ����� -", new double[][] { { Vcrm }, { Pcrm } });
		renderer.setSeriesLinesVisible(data.getSeriesCount() - 1, false);
		renderer.setSeriesShapesVisible(data.getSeriesCount() - 1, true);
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);

		data.addSeries("����������� ����� +", new double[][] { { Vcrp }, { Pcrp } });
		renderer.setSeriesLinesVisible(data.getSeriesCount() - 1, false);
		renderer.setSeriesShapesVisible(data.getSeriesCount() - 1, true);
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);

		// data.addSeries("������ ��� �����������(-) ����������� "+Tcrm, new
		// double[][] { vaxis, fcrmGraph });
		for (int w = 0; w < tN; w++)
		{
			double t = tbegin + (tend - tbegin) / (tN - 1) * w;
			data.addSeries("" + t, new double[][] { vaxis, funcGraph[w] });
			renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);
		}
		final XYPointerAnnotation pointerp = new XYPointerAnnotation("Критическа(+) точка. V = " + Vcrp + "; P = " + Pcrp + "; T = " + Tcrp + ".", Vcrp, Pcrp,
				-Math.PI / 4.0);
		pointerp.setBaseRadius(45.0);
		pointerp.setTipRadius(10.0);
		pointerp.setFont(new Font("SansSerif", Font.PLAIN, 12));
		pointerp.setPaint(Color.black);
		pointerp.setTextAnchor(TextAnchor.HALF_ASCENT_LEFT);
		plot.addAnnotation(pointerp);

		final XYPointerAnnotation pointerm = new XYPointerAnnotation("Критическая(-) точка. V = " + Vcrm + "; P = " + Pcrm + "; T = " + Tcrm + ".", Vcrm, Pcrm,
				Math.PI / 4.0);
		pointerm.setBaseRadius(45.0);
		pointerm.setTipRadius(10.0);
		pointerm.setFont(new Font("SansSerif", Font.PLAIN, 12));
		pointerm.setPaint(Color.black);
		pointerm.setTextAnchor(TextAnchor.HALF_ASCENT_LEFT);
		plot.addAnnotation(pointerm);

		final XYTextAnnotation annotation = new XYTextAnnotation("Delta = " + d + "; Lamda = " + l + "; Alpha = " + a + "; Beta = " + b + ".", 6, -0.05);
		annotation.setFont(new Font("SansSerif", Font.PLAIN, 12));
		plot.addAnnotation(annotation);

//		secondZone(data, renderer, Vcrp, Tcrp, a, b, d, l);
		
//		frame = createFrame(chart, "Второй случай", 0, 12, -0.1, 0.4);
		frame = createFrame(chart, "Второй случай", 0, 9, -0.1, 2.2);
	}

	private static ChartFrame createFrame(JFreeChart chart, String name, double xb, double xe, double yb, double ye)
	{
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

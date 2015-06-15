import java.math.RoundingMode;
import java.text.DecimalFormat;

public class FFT {

//	int n, m;

	// Lookup tables. Only need to recompute when size of FFT changes.
	double[] cos;
	double[] sin;

	double[] window;

	public FFT(int n) {
//		this.n = n;
//		this.m = (int) (Math.log(n) / Math.log(2));

		// Make sure n is a power of 2
		
		int m = (int) (Math.log(n) / Math.log(2));
		if (n != (1 << m))
			throw new RuntimeException("FFT length must be power of 2");

		// precompute tables
		cos = new double[n / 2];
		sin = new double[n / 2];

		// for(int i=0; i<n/4; i++) {
		// cos[i] = Math.cos(-2*Math.PI*i/n);
		// sin[n/4-i] = cos[i];
		// cos[n/2-i] = -cos[i];
		// sin[n/4+i] = cos[i];
		// cos[n/2+i] = -cos[i];
		// sin[n*3/4-i] = -cos[i];
		// cos[n-i] = cos[i];
		// sin[n*3/4+i] = -cos[i];
		// }

		for (int i = 0; i < n / 2; i++) {
			cos[i] = Math.cos(-2 * Math.PI * i / n);
			sin[i] = Math.sin(-2 * Math.PI * i / n);
		}

	}

	/***************************************************************
	 * 00089 * fft.c 00090 * Douglas L. Jones 00091 * University of Illinois at
	 * Urbana-Champaign 00092 * January 19, 1992 00093 *
	 * http://cnx.rice.edu/content/m12016/latest/ 00094 * 00095 * fft: in-place
	 * radix-2 DIT DFT of a complex input 00096 * 00097 * input: 00098 * n:
	 * length of FFT: must be a power of two 00099 * m: n = 2**m 00100 *
	 * input/output 00101 * x: double array of length n with real part of data
	 * 00102 * y: double array of length n with imag part of data 00103 * 00104
	 * * Permission to copy and use this program is granted 00105 * as long as
	 * this header is included. 00106
	 ****************************************************************/
	public void fft(double[] x, double[] y) {
		int i, j, k, n1, n2, a;
		double c, s, e, t1, t2;

		// Bit-reverse
		j = 0;
		int n = x.length;
		int m = (int) (Math.log(n) / Math.log(2));
		n2 = n / 2;
		for (i = 1; i < n - 1; i++) {
			n1 = n2;
			while (j >= n1) {
				j = j - n1;
				n1 = n1 / 2;
			}
			j = j + n1;

			if (i < j) {
				t1 = x[i];
				x[i] = x[j];
				x[j] = t1;
				t1 = y[i];
				y[i] = y[j];
				y[j] = t1;
			}
		}

		// FFT
		n1 = 0;
		n2 = 1;

		for (i = 0; i < m; i++) {
			n1 = n2;
			n2 = n2 + n2;
			a = 0;

			for (j = 0; j < n1; j++) {
				// c = cos[a];
				// s = sin[a];
				c = Math.cos(-2 * Math.PI * a / n);
				s = Math.sin(-2 * Math.PI * a / n);
				a += 1 << (m - i - 1);

				for (k = j; k < n; k = k + n2) {
					t1 = c * x[k + n1] - s * y[k + n1];
					t2 = s * x[k + n1] + c * y[k + n1];
					x[k + n1] = x[k] - t1;
					y[k + n1] = y[k] - t2;
					x[k] = x[k] + t1;
					y[k] = y[k] + t2;
				}
			}
		}
	}

	public void ifft(double[] x, double[] y) {

		for (int i = 0; i < y.length; i++) {

			y[i] = -y[i];

		}

		fft(x, y);

		for (int i = 0; i < y.length; i++) {

			y[i] = -y[i];

			x[i] = x[i] / 8;

		}

	}

	public void fft2(double[][] x, double[][] y) {

		double[] x1 = new double[x[0].length], y1 = new double[x[0].length];

		for (int m = 0; m < x.length; m++) {
			for (int n = 0; n < x[m].length; n++) {
				x1[n] = x[m][n];
				y1[n] = y[m][n];
				
				
			}
//			System.out.println("x: "+x1.length);
			fft(x1, y1);
			for (int n = 0; n < x[m].length; n++) {
				x[m][n] = x1[n];
				y[m][n] = y1[n];
			}
		}

		// transform second dimension
		double[] x2 = new double[x.length], y2 = new double[x.length];
		for (int n = 0; n < x[0].length; n++) {
			for (int m = 0; m < x.length; m++) {
				// v[m] = new Complex(GF[m][n]);
				x2[m] = x[m][n];
				y2[m] = y[m][n];

			}
//			System.out.println("x2: "+x2.length);
			fft(x2, y2);
			for (int m = 0; m < x.length; m++) {
				x[m][n] = x2[m];
				y[m][n] = y2[m];
			}
		}

	}

	// Test the FFT to make sure it's working
	public static void main(String[] args) {
		int N = 8;int M=4;

		FFT fft = new FFT(N);

		// double[] window = fft.getWindow();
		double[] re = new double[N];
		double[] im = new double[N];

		// Impulse
		re[0] = 1;
		im[0] = 0;
		for (int i = 1; i < N; i++)
			re[i] = im[i] = 0;
		beforeAfter(fft, re, im);

		System.out.println("----------------After IFFT--------------");
		fft.ifft(re, im);
		printReIm(re, im);

		// Nyquist
		for (int i = 0; i < N; i++) {
			re[i] = Math.pow(-1, i);
			im[i] = 0;
		}
		beforeAfter(fft, re, im);

		System.out.println("----------------After IFFT--------------");
		fft.ifft(re, im);
		printReIm(re, im);

		// Single sin
		for (int i = 0; i < N; i++) {
			re[i] = Math.cos(2 * Math.PI * i / N);
			im[i] = 0;
		}
		beforeAfter(fft, re, im);

		System.out.println("----------------After IFFT--------------");
		fft.ifft(re, im);
		printReIm(re, im);

		// Ramp
		for (int i = 0; i < N; i++) {
			re[i] = i;
			im[i] = 0;
		}
		beforeAfter(fft, re, im);

		System.out.println("----------------After IFFT--------------");
		fft.ifft(re, im);
		printReIm(re, im);

		System.out.println("----------------Before FFT2--------------");

		double[][] re2 = new double[N][M];
		double[][] im2 = new double[N][M];

		for (int i = 0; i < N; i++) {

			for (int j = 0; j < M; j++) {

				re2[i][j] = i;
				im2[i][j] = 0;

				System.out.print(" " + re2[i][j] + "+j" + im2[i][j]);
			}

			System.out.println();
		}

		fft.fft2(re2, im2);
		System.out.println("----------------After FFT2--------------");
		for (int i = 0; i < N; i++) {

			for (int j = 0; j < M; j++) {

				System.out.print(" " + re2[i][j] + "+j" + im2[i][j]);
			}

			System.out.println();
		}

		System.out.println("----------------After IFFT2--------------");
		fft.ifft2(re2, im2);
		for (int i = 0; i < N; i++) {

			for (int j = 0; j < M; j++) {
				DecimalFormat df = new DecimalFormat("##.##");
				df.setRoundingMode(RoundingMode.DOWN);
				System.out.print("    " + df.format(re2[i][j]) + "+j"
						+ df.format(im2[i][j]));
			}

			System.out.println();
		}

		long time = System.currentTimeMillis();
		double iter = 30000;
		for (int i = 0; i < iter; i++)
			fft.fft(re, im);
		time = System.currentTimeMillis() - time;
		System.out.println("Averaged " + (time / iter) + "ms per iteration");
	}

	public void ifft2(double[][] x, double[][] y) {

		double[] x1 = new double[x.length], y1 = new double[x.length];

		for (int m = 0; m < y.length; m++) {
			for (int n = 0; n < y[m].length; n++) {
				
				y[m][n] = -y[m][n];
			}
		}
		
		fft2(x, y);
		for (int m = 0; m < y.length; m++) {
			for (int n = 0; n < y[m].length; n++) {
				
				y[m][n] = -y[m][n]/(y.length*y[m].length);
				x[m][n] = x[m][n]/(y.length*y[m].length);
			}
		}
		// for (int m = 0; m < x.length; m++) {
		// for (int n = 0; n < x[m].length; n++) {
		// x1[n] = x[m][n];
		// y1[n] = y[m][n];
		// }
		// ifft(x1, y1);
		// for (int n = 0; n < x[m].length; n++) {
		// x[m][n] = x1[n];
		// y[m][n] = y1[n];
		// }
		// }
		//
		// // transform second dimension
		//
		// for (int n = 0; n < x[0].length; n++) {
		// for (int m = 0; m < x.length; m++) {
		// x1[m] = x[m][n];
		// y1[m] = y[m][n];
		// }
		// ifft(x1, y1);
		// for (int m = 0; m < x.length; m++) {
		// x[m][n] = x1[m];
		// y[m][n] = y1[m];
		// }
		// }
	}

	protected static void beforeAfter(FFT fft, double[] re, double[] im) {
		System.out.println("Before: ");
		printReIm(re, im);
		fft.fft(re, im);
		System.out.println("After: ");
		printReIm(re, im);
	}

	protected static void printReIm(double[] re, double[] im) {
		System.out.print("Re: [");
		for (int i = 0; i < re.length; i++)
			System.out.print(((int) (re[i] * 1000) / 1000.0) + " ");

		System.out.print("]\nIm: [");
		for (int i = 0; i < im.length; i++)
			System.out.print(((int) (im[i] * 1000) / 1000.0) + " ");

		System.out.println("]");
	}
}

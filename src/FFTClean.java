import java.math.RoundingMode;
import java.text.DecimalFormat;

public class FFTClean {

	public FFTClean(int n) {

		int m = (int) (Math.log(n) / Math.log(2));
		if (n != (1 << m))
			throw new RuntimeException("FFT length must be power of 2");

	}

	public void fft(double[] x, double[] y,Boolean is_row) {
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

		// System.out.println();
		// FFT
		n1 = 0;
		n2 = 1;

		for (i = 0; i < m; i++) {
			n1 = n2;
			n2 = n2 + n2;
			a = 0;

			for (j = 0; j < n1; j++) {

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

		fft(x, y,true);

		for (int i = 0; i < y.length; i++) {

			y[i] = -y[i];

			x[i] = x[i] / y.length;

		}

	}

	public void fft2(double[][] x, double[][] y) {

		double[] x1 = new double[x[0].length], y1 = new double[y[0].length];

		for (int m = 0; m < x.length; m++) {
			for (int n = 0; n < x[m].length; n++) {
				x1[n] = x[m][n];
				y1[n] = y[m][n];
			
			}
	
			fft(x1, y1,true);
			for (int n = 0; n < x[m].length; n++) {
				x[m][n] = x1[n];
				y[m][n] = y1[n];
			}
		}

		// transform second dimension
		double[] x2 = new double[x.length], y2 = new double[x.length];
		for (int n = 0; n < x[0].length; n++) {

			for (int m = 0; m < x.length; m++) {

				x2[m] = x[m][n];
				y2[m] = y[m][n];
			}

			fft(x2, y2,false);
			for (int m = 0; m < x.length; m++) {
				x[m][n] = x2[m];
				y[m][n] = y2[m];
			}
		}

	}

	// Test the FFT to make sure it's working
	public static void main(String[] args) {
		int N = 8;
		int M = 8;

		FFTClean fft = new FFTClean(N);

		System.out.println("----------------Before FFT2--------------");
		double[][] re2 = {
				{ 0.28235295, 0.36078432, 0.4117647, 0.3254902, 0.2627451,
						0.2509804, 0.24313726, 0.27058825 },

				{ 0.89411765, 0.88235295, 0.43137255, 0.25490198, 0.45882353,
						0.019607844, 0.0, 0.07058824 },

				{ 0.9254902, 0.9137255, 0.75686276, 0.88235295, 0.8980392,
						0.7176471, 0.03137255, 0.16862746 },
				{ 0.8862745, 0.5294118, 0.80784315, 0.9490196, 0.8509804,
						0.8509804, 0.7607843, 0.78431374 },

				{ 0.7411765, 0.8509804, 0.93333334, 0.85490197, 0.8745098,
						0.9137255, 0.93333334, 0.85490197 },

				{ 0.19215687, 0.20392157, 0.7137255, 0.74509805, 0.8156863,
						0.81960785, 0.81960785, 0.8784314 },

				{ 0.27058825, 0.050980393, 0.0, 0.050980393, 0.08627451,
						0.06666667, 0.6156863, 0.8039216 },

				{ 0.27058825, 0.007843138, 0.0, 0.019607844, 0.078431375,
						0.078431375, 0.60784316, 0.7372549 } };

		// double[][] re2 = new double[N][M];
		double[][] im2 = new double[N][M];
		//
		//
		//
		for (int i = 0; i < N; i++) {

			for (int j = 0; j < M; j++) {

				// re2[i][j] = i;
				// im2[i][j] = 0;

				// re2[i][j] = i + j;
				im2[i][j] = 0;

				System.out.print(" " + re2[i][j] + "+j" + im2[i][j]);
			}

			System.out.println();
		}

		// double[][] re2={{ 0.28235295,0.36},{0.28235295,0.36}};

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

	}

	public void image_filter(int D0) {
		int M = 9;
		int N = 9;

		double[][] distance = new double[M][N];

		for (int i = 1; i < M; i++) {
			for (int j = 1; j < N; j++) {
				int a = N / 2;
				int b = M / 2;
				if (i == 1 && j <= a)
					distance[i][j] = j - 1;
				else if (i == 1 && j > a)
					distance[i][j] = 1 + a - (j - a);
				else if (j == 1 && i <= b)
					distance[i][j] = i - 1;
				else if (j == 1 && i > b)
					distance[i][j] = 1 + b - (i - b);
				else {
					double w = distance[i][1];
					double x = distance[1][j];
					distance[i][j] = Math.sqrt(w * w + x * x);
				}

				double H = 1 - Math.exp(-(distance[i][j]) * (distance[i][j])
						/ (2 * (D0 * D0)));

			}

		}
	}

	public void ifft2(double[][] x, double[][] y) {

		for (int m = 0; m < y.length; m++) {
			for (int n = 0; n < y[m].length; n++) {

				y[m][n] = -y[m][n];
			}
		}

		fft2(x, y);
		for (int m = 0; m < y.length; m++) {
			for (int n = 0; n < y[m].length; n++) {

				y[m][n] = -y[m][n] / (y.length * y[m].length);
				x[m][n] = x[m][n] / (y.length * y[m].length);
			}
		}

	}

}

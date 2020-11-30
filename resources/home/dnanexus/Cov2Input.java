/*********************************************************************
 * Create input file for CONSERTING algorithm
 * @auth Xiang Chen
 * @date 08/11/2010
 *
 * Usage:
 *   Using a fixed size window:
 *	java Cov2Input -f window_size -i input -o output [-sd val]
 *	    window_size: size of the window
 *	    input: input filename
 *	    output: output filename
 *	    -sd val: (optional paramemter) if set, position with depth val * sd
 *		away from the mean will be removed 
 ********************************************************************/

import java.io.*;
import java.util.*;

public class Cov2Input {
    public void process(String[] args) {
	int mode = -1; //unset mode, 0: fixed size; 1: variable size
	int window_size = -1;
	int total_coverage = -1;
	double thesh = -1;
	String input = null, output = null;
	double thresh = -1;
	if (args.length %2 > 0) {
	    System.err.println("Incorrect number of parameters.");
	    System.exit(-1);
	} else {
	    for (int i = 0; i < args.length; i++) {
		if (args[i].compareToIgnoreCase("-f") == 0) {
		   if (mode < 0) {
			mode = 0;
			window_size = Integer.parseInt(args[++i]);
		   } else {
			System.err.println("Mode set as variable window already.  Ignored: " + args[i++] + " " + args[i]);
		   }
		} else if (args[i].compareToIgnoreCase("-s") == 0) {
		   if (mode < 0) {
			mode = 1;
			total_coverage = Integer.parseInt(args[++i]);
		   } else {
			System.err.println("Mode set as fixed window already.  Ignored: " + args[i++] + " " + args[i]);
		   }
		} else if (args[i].compareToIgnoreCase("-i") == 0) {
		    input = args[++i];
		} else if (args[i].compareToIgnoreCase("-o") == 0) {
		    output = args[++i];
		} else if (args[i].compareToIgnoreCase("-sd") == 0) {
		    thresh = Double.parseDouble(args[++i]);
		} else {
		   System.out.println("Parameter not recognized.  Quitting.");
		   System.exit(-1);
		}
	    }
	}
	if (mode < 0 || input == null || output == null) {
	    System.err.println("Missing parameters.  Quitting!");
	    System.exit(-1);
	}
	if (mode == 1) {
	    System.out.println("Not implemented yet.  Sorry.");
	    System.exit(1);
	}
	Summary[] summary = new Summary[1000];
	try {
	    BufferedReader br = new BufferedReader(new FileReader(input));
	    PrintWriter pw = new PrintWriter(new FileWriter(output));
	    String line = br.readLine();
	    if (line == null) {
	    	pw.println("Position\tMean\tLog\tSD\tOutliers");
		pw.close();
		br.close();
		System.exit(0);
	    }
	    StringTokenizer st = null;
	    boolean hg19 = (line.indexOf("chr") >= 0);
	    int pos = 0, count = 0, zeros = 0;
	    int[] cov = new int[window_size];
	    if (hg19) {
		st = new StringTokenizer(line, " =");
		for (int i = 0; i < 4; i++) st.nextToken();
		pos = Integer.parseInt(st.nextToken()) - 1;
		count = pos / window_size;
		pos = pos % window_size;
		if (count >= summary.length) {
			Summary[] tmp = new Summary[2 * count];
			System.arraycopy(summary, 0, tmp, 0, summary.length);
			summary = tmp;
		}
		for (int i = 0; i < count; i++) summary[i] = getSummary(cov, window_size, thresh);
		line = br.readLine();
		//System.out.println(line);
	    }
	    while (line != null) {
		if (hg19) {
		    //System.out.println(line);
		    if (line.indexOf("chr") >= 0) {
			st = new StringTokenizer(line, " =");
			for (int i = 0; i < 4; i++) st.nextToken();
			int pos2 = Integer.parseInt(st.nextToken()) - 1;
			int count2 = pos2 / window_size;
			if (count2 >= summary.length) {
			    Summary[] tmp = new Summary[2 * count2];
			    System.arraycopy(summary, 0, tmp, 0, summary.length);
			    summary = tmp;
			}
			if (count2 > count) {
			    for (int i = pos; i < window_size; i++) cov[i] = 0;
			    summary[count] = getSummary(cov, window_size, thresh);
			    for (int i = 0; i < window_size; i++) cov[i] = 0;
			}
			for (int i = count + 1; i < count2; i++) summary[i] = getSummary(cov, window_size, thresh);
			count = count2;
			pos = pos2 % window_size;
		    } else {
			cov[pos++] = Integer.parseInt(line);
		    }
		} else {
		    st = new StringTokenizer(line, " \t,");
		    st.nextToken();
		    try {
			cov[pos++] = Integer.parseInt(st.nextToken());
		    } catch (NoSuchElementException ne) {
			System.err.println(line);
			ne.printStackTrace();
			System.exit(-3);
		    }
		}
		if (pos == window_size) {
		    if (count == summary.length) {
			Summary[] tmp = new Summary[2 * summary.length];
			System.arraycopy(summary, 0, tmp, 0, summary.length);
			summary = tmp;
		    }
		    summary[count++] = getSummary(cov, window_size, thresh);
		    pos = 0;
		}
		line = br.readLine();
	    }
	    if (pos > 0) {
		if (count == summary.length) {
		    Summary[] tmp = new Summary[1 + summary.length];
		    System.arraycopy(summary, 0, tmp, 0, summary.length);
		    summary = tmp;
		}
		summary[count++] = getSummary(cov, pos, thresh);
		if (summary[count - 1].mean < 2e-3) zeros++;
	    }
	    br.close();
	    if (summary.length > count) {
		Summary[] tmp = new Summary[count];
		System.arraycopy(summary, 0, tmp, 0, count);
		summary = tmp;
	    }
	    double sum = 0;
	    for (int i = 0; i < count; i++) sum += summary[i].mean;
	    sum = Math.log(sum / (count - zeros));
	    pw.println("Position\tMean\tLog\tSD\tOutliers");
	    for (int i = 0; i < count - 1; i++) {
		pw.print((i * window_size + window_size / 2) + "\t" +  summary[i].mean + "\t");
		if (summary[i].mean < 2e-3) {
		    pw.print("NA");
		} else {
		    pw.print((Math.log(summary[i].mean) - sum) / Math.log(2));
		}
		pw.println("\t" + summary[i].sd + "\t" + summary[i].outliers);
	    }
	    pw.print(((count - 1) * window_size + pos / 2) + "\t" + summary[count - 1].mean + "\t");
	    if (summary[count - 1].mean < 2e-3) {
		pw.print("NA");
	    } else {
		pw.print((Math.log(summary[count - 1].mean) - sum) / Math.log(2));
	    }
	    pw.println("\t" + summary[count - 1].sd + "\t" + summary[count - 1].outliers);
	    pw.close();
	} catch (IOException e) {
	    e.printStackTrace();
	    System.exit(-2);
	} 
    }





    public Summary getSummary(int[] val, int num, double thresh) {
	double x = 0, xx = 0, sd = 0, max = val[0], min = val[0];
	int outliers = 0;
	for (int i = 0; i < num; i++) {
	    x += val[i];
	    xx += (val[i] * val[i]);
	    if (val[i] > max) max = val[i];
	    if (val[i] < min) min = val[i];
	}
	sd = (num > 1) ? Math.sqrt((xx - x * x / num) / (num - 1)) : 0;
	if (thresh > 0) {
	    double upper = x / num + thresh * sd, lower = x / num - thresh * sd;
	    if ((max > upper) || (min < lower)) {
		for (int i = 0; i < num; i++) {
		    if ((val[i] > upper) || (val[i] < lower)) {
			x -= val[i];
			xx -= (val[i] * val[i]);
			outliers++;
		    }
		}
		sd = (num - outliers > 1) ? Math.sqrt((xx - x * x / (num - outliers)) / (num - outliers - 1)) : 0;
	    }
	}
	Summary s = new Summary(x / (num - outliers), sd, outliers);
	return s;
    }
	
    public static void main(String[] args) {
	Cov2Input c = new Cov2Input();
	c.process(args);
    }


    class Summary {
	double mean, sd;
	int outliers;

	Summary(double means, double sds, int outlier) {
	    mean = means;
	    sd = sds;
	    outliers = outlier;
	}
    }
}

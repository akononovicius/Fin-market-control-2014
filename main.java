import java.io.*;
import java.math.*;
import java.math.BigDecimal;
import java.util.Random;
import java.text.SimpleDateFormat;
import java.util.Calendar;

public class main {
	public static void getCommandLine(String[] str) {
		if(!(str.length>0)) {
			System.out.println("Initiating with default settings...");
			return ;
		}
		int i=0;
		try {
			for(i=0;i<str.length;i++) {
				if(str[i].equalsIgnoreCase("--jobs") || str[i].equalsIgnoreCase("-j")) {
					i++;
					commonVariables.core=Integer.parseInt(str[i]);
					i++;
					commonVariables.realizations=Integer.parseInt(str[i]);
				} else if(str[i].equalsIgnoreCase("--realizationPoints") || str[i].equalsIgnoreCase("-rp")) {
					i++;
					int pts=Integer.parseInt(str[i]);
					pts=(int)Math.pow(2,Math.ceil(Math.log(pts)/Math.log(2)));
					commonVariables.realizationPoints=pts;
				} else if(str[i].equalsIgnoreCase("--outPoints") || str[i].equalsIgnoreCase("-op")) {
					i++;
					commonVariables.outputPoints=Integer.parseInt(str[i]);
				} else if(str[i].equalsIgnoreCase("-dt")) {
					i++;
					commonVariables.dt=Double.parseDouble(str[i]);
				} else if(str[i].equalsIgnoreCase("-dtSec")) {
					i++;
					commonVariables.dtSec=Double.parseDouble(str[i]);
				} else if(str[i].equalsIgnoreCase("--outputTemplate") || str[i].equalsIgnoreCase("-ot")) {
					i++;
					commonVariables.outputTemplate=str[i];
				} else if(str[i].equalsIgnoreCase("-m")) {
					i++;
					commonVariables.M=Double.parseDouble(str[i]);
				} else if(str[i].equalsIgnoreCase("-pars")) {
					i++;
					commonVariables.ecf=Double.parseDouble(str[i]);
					i++;
					commonVariables.efc=Double.parseDouble(str[i]);
					i++;
					commonVariables.ecc=Double.parseDouble(str[i]);
					i++;
					commonVariables.H=Double.parseDouble(str[i]);
				} else if(str[i].equalsIgnoreCase("-k")) {
					i++;
					commonVariables.k=Double.parseDouble(str[i]);
				} else if(str[i].equalsIgnoreCase("--outTime")) {
					commonVariables.outTime=true;
				} else if(str[i].equalsIgnoreCase("--symmetricImpact")) {
					commonVariables.impact=2;
				} else if(str[i].equalsIgnoreCase("--chartistImpact")) {
					commonVariables.impact=1;
				} else if(str[i].equalsIgnoreCase("--limits") || str[i].equalsIgnoreCase("-lim")) {
					i++;
					commonVariables.min=Double.parseDouble(str[i]);
					i++;
					commonVariables.max=Double.parseDouble(str[i]);
					commonVariables.useLimits=true;
				} else if(str[i].equalsIgnoreCase("--pdfLimits") || str[i].equalsIgnoreCase("-pdflim")) {
					i++;
					commonVariables.pdfMin=Double.parseDouble(str[i]);
					i++;
					commonVariables.pdfMax=Double.parseDouble(str[i]);
				} else {
					System.out.println("Error in input!");
					System.out.println("-> "+str[i]+" <-");
					System.exit(0);
				}
				if(!(i<str.length-1)) System.out.println("Commandline arguments were fully processed.");
			}
		} catch (Exception e) {
			System.out.println("Error in input!");
			System.out.println("-> "+str[i]+" <-");
			System.exit(0);
		}
	}
	public static void main(String [] args) {
		getCommandLine(args);
		launcher threadLauncher=new launcher();
	}
}

class launcher {
	private SimpleDateFormat sdf=new SimpleDateFormat("HH:mm:ss.SSS yyyy-MM-dd");
	private long startTime=0;//time (ms) when program was started
	private long timerPause=30000;//minimal pause (ms) between reports
	private long lastReport=0;//time (ms) of the last report
	private double[] pdf=new double[commonVariables.multiOutputPoints];//histogram
	private double[] spec=new double[commonVariables.realizationPoints];//spectra
	private int[] perc=null;//completion percentages
	private thread[] threads=null;//all threads
	private int completed=0;//how many threads are completed
	public launcher() {
		//which impact scenario is used?
		// adjust epsilon values accordingly
		if(commonVariables.impact==2) {
			commonVariables.ecf+=0.5*commonVariables.M;
			commonVariables.efc+=0.5*commonVariables.M;
			commonVariables.ecc+=0.5*commonVariables.M;
		} else if(commonVariables.impact==1) {
			commonVariables.efc+=commonVariables.M;
			commonVariables.ecc+=0.5*commonVariables.M;
		} else if(commonVariables.impact==0) {
			commonVariables.ecf+=commonVariables.M;
		}
		startTime=(Calendar.getInstance()).getTimeInMillis();
		System.out.println("Thread launcher has started!");
		int realCore=Math.max(commonVariables.realizations/commonVariables.core,1);
		perc=new int[commonVariables.core];
		threads=new thread[commonVariables.core];
		for(int i=0;i<commonVariables.core;i++) {
			perc[i]=0;
			threads[i]=new thread(commonVariables.ecf,commonVariables.efc,commonVariables.ecc,commonVariables.H,commonVariables.k,commonVariables.min,commonVariables.max,commonVariables.pdfMin,commonVariables.pdfMax,commonVariables.useLimits,commonVariables.dt,commonVariables.realizationPoints,realCore,commonVariables.multiOutputPoints,commonVariables.outTime,i,this);
			System.out.println("Thread nr. "+i+" was launched. Sent "+realCore+" jobs.");
		}
	}
	//thread should pass histogram array, spectra array and its number
	synchronized public void put(double[] p, double[] s, int nr) {
		completed++;
		for(int i=0;i<pdf.length;i++) pdf[i]+=p[i];
		for(int i=0;i<spec.length;i++) spec[i]+=s[i];
		perc[nr]=100;
		if(completed>=commonVariables.core) {
			System.out.println("All threads has reported their results. Preparing results for final output.");
			int first=-1;
			int last=-1;
			for(int i=0;i<pdf.length;i++) {
				if((pdf[i]>0)&&(first==-1)&&(i>0)) first=i;
				if(pdf[i]>0) last=i;
				pdf[i]/=((double)commonVariables.core);
			}
			if(first>0) first--;
			if(last<pdf.length-1) last++;
			for(int i=0;i<spec.length;i++) spec[i]/=((double)commonVariables.core);
			double xstep=(commonVariables.pdfMax-commonVariables.pdfMin)/((double)commonVariables.multiOutputPoints);
			outputarr(commonFunctions.pdfModification(pdf,commonVariables.pdfMin+first*xstep,commonVariables.pdfMin+last*xstep,commonVariables.pdfMin,xstep,commonVariables.outputPoints),commonVariables.outputTemplate+"dist",6,false);
			outputarr(commonFunctions.specModification(spec,commonVariables.dt/commonVariables.dtSec,commonVariables.outputPoints,false),commonVariables.outputTemplate+"spec",6,false);
			System.out.println("Done. Exiting.");
			System.exit(0);
		}
	}
	//function for a thread to report its progess (number, percentage)
	synchronized public void report(int nr, int vperc) {
		perc[nr]=vperc;
		Calendar cal=Calendar.getInstance();
		long nowMs=cal.getTimeInMillis();
		if(nowMs-lastReport>timerPause) {
			lastReport=nowMs;
			System.out.println(sdf.format(cal.getTime()));
			nowMs=cal.getTimeInMillis()-startTime;
			System.out.println("  Runtime: "+nowMs+" ms");
			int tmin=101;
			for(int i=0;i<perc.length;i++) {
				if(tmin>perc[i]) tmin=perc[i];
			}
			System.out.println("  Completion: "+tmin+"%");
			for(int i=0;i<perc.length;i++) {
				System.out.println("    ["+i+" - "+perc[i]+"%]");
			}
			long eta=(long)Math.floor((100.0/((double)tmin)-1.0)*nowMs);
			cal.setTimeInMillis(startTime+nowMs+eta);
			System.out.println("  ETA: "+sdf.format(cal.getTime())+" (+"+eta+" ms)");
			System.out.flush();
		}
	}
	//output 2D array
	public static void outputarr(double[][] arr, String name, int afterComma, boolean aboveZero) {
		BufferedWriter out=null;
		String fileName=name;
		try {
			File file=new File(fileName);
			file.createNewFile();
			out=new BufferedWriter(new FileWriter(fileName));
			for(int i=0;i<arr.length;i++) {
				if(((aboveZero)&&(arr[i][1]>0))||(!aboveZero)) {
					for(int j=0;j<arr[i].length-1;j++) out.write(round(arr[i][j],afterComma)+" ");
					out.write(round(arr[i][arr[i].length-1],afterComma)+"");
					out.newLine();
				}
			} 
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	//rounding
	public static double round(double d, int i) {
		int j = i;
		BigDecimal bigdecimal = new BigDecimal(d);
		bigdecimal = bigdecimal.setScale(j,bigdecimal.ROUND_HALF_UP);
		d = bigdecimal.doubleValue();
		return d;
	}
}

class thread implements Runnable {
	private Random gen=new Random();
	private double dt=1;//discretization time step
	private double pdfMin=0.01;//minimum value to include in pdf
	private double pdfMax=100;//maximum value to include in pdf
	private double[] x=null;//time series array
	private double[] pdf=null;//pdf array
	private double[] spec=null;//psd array
	private int realizationPoints=32768;//points in realization
	private int realizations=1;//number of realizations to obtain
	private doubleFundMood model=null;//the model itself
	private launcher parrent=null;//parrent object
	private int reportNext=328;//report progress after this number of new points was obtained
	private int percentageDone=1;//percentageDone-1 % already completed
	private int nr=1;//thread number
	private boolean outTime=true;//output time series?
	public thread() {}
	//properly initialize concurrent evaluation of the model
	public thread(double vecf, double vefc, double vecc, double vH, double vk, double vfmin, double vfmax, double vpmin, double vpmax, boolean vf, double vdt, int vrp, int vr, int vo, boolean outT, int vnr, launcher vpar) {
		model=new doubleFundMood(vecf,vefc,vecc,vH,vk,vfmin,vfmax,vf,vnr);
		pdfMin=vpmin;
		pdfMax=vpmax;
		dt=vdt;
		realizationPoints=vrp;
		realizations=vr;
		reportNext=(int)Math.floor(percentageDone*realizations*(((double)realizationPoints)/100.0));
		x=new double[realizationPoints];
		pdf=new double[vo];
		spec=new double[realizationPoints];
		for(int i=0;i<realizationPoints;i++) spec[i]=0;
		for(int i=0;i<vo;i++) pdf[i]=0;
		parrent=vpar;
		nr=vnr;
		outTime=outT;
		new Thread(this,"Thread nr. "+nr).start();
	}
	public void run() {
		double pdfStep=(pdfMax-pdfMin)/((double)pdf.length);
		for(int j=0;j<realizations;j++) {
			double[] pdft=new double[pdf.length];
			for(int i=0;i<pdf.length;i++) pdft[i]=0;
			int alreadyMade=j*realizationPoints;
			for(int i=0;i<realizationPoints;i++) {
				x[i]=Math.abs(model.step(dt));//obtain absolute relative log-price
				if((x[i]>=pdfMin)&&(x[i]<=pdfMax)) {//appropriately modify histogram
					int tmp=(int)Math.floor((x[i]-pdfMin)/pdfStep);
					if((tmp>=0)&&(tmp<pdft.length)) pdft[tmp]++;
				}
				if((i+alreadyMade)>=reportNext) {
					parrent.report(nr,percentageDone);
					percentageDone++;
					reportNext=(int)Math.floor(percentageDone*realizations*(((double)realizationPoints)/100.0));
				}
			}
			if(outTime) outputarr(x,commonVariables.outputTemplate+"c"+nr+".r"+j+".time",3);
			for(int i=0;i<pdf.length;i++) {
				pdft[i]/=((double)realizationPoints);
				pdf[i]+=pdft[i];
			}
			pdft=null;
			double[] spect=preformRealFFT(x);
			for(int i=0;i<=spec.length/2;i++) spec[i]+=spect[i];
		}
		if(realizations>1) {//normalize by the number of realizations
			for(int i=0;i<pdf.length;i++) pdf[i]/=((double)realizations);
			for(int i=0;i<=spec.length/2;i++) spec[i]/=((double)realizations);
		}
		parrent.put(pdf,spec,nr);
	}
	//output 1D array
	public static void outputarr(double[] arr, String name, int afterComma) {
		BufferedWriter out=null;
		String fileName=name;
		try {
			File file=new File(fileName);
			file.createNewFile();
			out=new BufferedWriter(new FileWriter(fileName));
			for(int i=0;i<arr.length;i++) {
				out.write(round(arr[i],afterComma)+"");
				out.newLine();
			} 
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	//rounding
	public static double round(double d, int i) {
		int j = i;
		BigDecimal bigdecimal = new BigDecimal(d);
		bigdecimal = bigdecimal.setScale(j,bigdecimal.ROUND_HALF_UP);
		d = bigdecimal.doubleValue();
		return d;
	}
	//fast fourier transform algorithm
	private static double matlog2=Math.log(2);
	public static double[] preformRealFFT(double[] arrx) {
		return preformRealFFT(arrx,true);
	}
	public static double[] preformRealFFT(double[] arrx, boolean subtractMean) {
		double[] rex=new double[arrx.length];
		if(subtractMean) {
			double mean=0;
			double tmean=0;
			for(int i=0;i<rex.length;i++) {
				tmean+=rex[i];
				if(tmean>Double.MAX_VALUE/100.0) {
					mean+=(tmean/((double)rex.length));
					tmean=0;
				}
			}
			if(tmean!=0) mean+=(tmean/((double)rex.length));
			for(int i=0;i<rex.length;i++) rex[i]-=mean;
		}
		System.arraycopy(arrx,0,rex,0,arrx.length);
		int n=rex.length;
		double[] imx=new double[n];
		for(int j=0;j<n;j++) {
			imx[j]=0;
		}
		int nm1=n-1;
		int nd2=(int)(n/2);
		int m=cint(Math.log(n)/matlog2);
		int j=nd2;
		for(int i=1;i<nm1;i++) {//bit reversal
			if(i<=j) {
				double tr=rex[j];
				double ti=imx[j];
				rex[j]=rex[i];
				imx[j]=imx[i];
				rex[i]=tr;
				imx[i]=ti;
			}
			int k=nd2;
			while(k<=j) {
				j-=k;
				k=(int)(k/2);
			}
			j+=k;
		}		
		for(int l=1;l<m+1;l++) {
			int le=cint(Math.pow(2,l));
			int le2=(int)(le/2);
			double ur=1;
			double ui=0;
			double tr=0;
			double ti=0;
			double sr=Math.cos(Math.PI/((double)le2));
			double si=-Math.sin(Math.PI/((double)le2));
			for(j=1;j<le2+1;j++) {
				int jm1=j-1;
				for(int i=jm1;i<=nm1;i+=le) {
					int ip=i+le2;
					tr=rex[ip]*ur-imx[ip]*ui;
					ti=imx[ip]*ur+rex[ip]*ui;
					rex[ip]=rex[i]-tr;
					imx[ip]=imx[i]-ti;
					rex[i]+=tr;
					imx[i]+=ti;
				}
				tr=ur;
				ur=tr*sr-ui*si;
				ui=tr*si+ui*sr;
			}
		}
		for(int i=0;i<=n/2;i++) {
			rex[i]=(rex[i]*rex[i]+imx[i]*imx[i]);
		}
		return rex;
	}
	private static int cint(double expr) {
		double dif=Math.abs((expr-(int)expr));
		if(dif==0.5) {
			if(((int)expr)%2==0) {
				return (int)expr;
			} else {
				return (int)expr+1;
			}
		} else if(dif>0.5) {
			return (int)expr+1;
		} else {
			return (int)expr;
		}
	}
}

//a class to store model variables
class commonVariables {
	public static String outputTemplate="rez.";//prefix of relative output path
	public static int core=1;//number of concurent threads to launch
	public static int realizations=1;//total number of realizations
	public static int realizationPoints=32768;//number of points in each realization
	public static int outputPoints=640;//number of points to output in final PDF and PSD
	public static int multiOutputPoints=1000000;//number of points used in the intermediate steps of calculating PDF
	public static double pdfMin=0.001;//minimal value included in PDF
	public static double pdfMax=1000;//maximal value included in PDF
	public static boolean useLimits=true;//enforce limiting of values?
	public static double min=0.0001;//minimal value allowed
	public static double max=0.9999;//maximal value allowed
	public static boolean outTime=false;//output time series?
	public static double dt=60*1e-8;//relative time step between consecutive points in realization
	public static double dtSec=1e-8;//how much relative time corresponds to second in real time
	public static double ecf=0.1;// c -> f individual transition rate
	public static double efc=3;// f -> c individual transition rate
	public static double ecc=3;// o -> p, p -> individual transition rate
	public static double H=300;//difference in time scale
	public static double k=0.1;//sprendinio tikslumas (maziau geriau)
	public static double M=0;//number of controlled agents
	public static int impact=0;//use fundamentalist impact (0), chartist impact (1) or chartist symmetric impact (2)
}

//a class with function to cleanup PSD and PDF
class commonFunctions {
	private static double matlog10=Math.log(10);
	//Faster lg calculation
	public static double LogBase10(double a) {
		return Math.log(a)/matlog10;
	}
	//linear spectra to logarithmic spectra
	public static double[][] specModification(double[] spec, double timeTick, int outPoints, boolean smoothen) {
		double normalization=commonFunctions.LogBase10(2*timeTick/spec.length);
		double scale=commonFunctions.LogBase10(spec.length)+commonFunctions.LogBase10(timeTick);
		double llim=0;
		double rlim=commonFunctions.LogBase10(spec.length/2.0);
		double lstep=(rlim-llim)/((double)outPoints);
		double clim=llim+lstep;
		int i=1;
		int inInterval=0;
		double[][] rez=new double[outPoints+1][2];
		int used=0;
		double total=0;
		double oldX=0;
		while(clim<=rlim) {
			while(commonFunctions.LogBase10(i)<clim) {
				total+=spec[i];
				i++;
				inInterval++;
			}
			if(total>0) {
				if(used==0) {
					oldX=Math.pow(10,clim-scale);
					rez[used][0]=commonFunctions.LogBase10(oldX/2.0);
					rez[used][1]=commonFunctions.LogBase10(total/((double)inInterval));
				} else {
					double newX=Math.pow(10,clim-scale);
					rez[used][0]=commonFunctions.LogBase10((newX+oldX)/2.0);
					rez[used][1]=commonFunctions.LogBase10(total/((double)inInterval));
					oldX=newX;
				}
				rez[used][1]+=normalization;
				used++;
			}
			inInterval=0;
			total=0;
			clim+=lstep;
		}
		double[][] rez2=new double[used][2];
		for(int ii=0;ii<used;ii++) {
			rez2[ii][0]=rez[ii][0];
			rez2[ii][1]=rez[ii][1];
		}
		if(smoothen) {
			/*movavg time window 3*/
			for(int ii=0;ii<used;ii++) {
				if((ii>0)&&(ii<used-1)) {
					rez2[ii][1]=(rez2[ii-1][1]+rez2[ii][1]+rez2[ii+1][1])/3.0;
				}
			}
		}
		return rez2;
	}
	//linear pdf to logarithmic pdf
	public static double[][] pdfModification(double[] pdf, double llim, double rlim, int outPoints) {
		return pdfModification(pdf,llim,rlim,1,1,outPoints);
	}
	public static double[][] pdfModification(double[] pdf, double llim, double rlim, double xlim, double xstep, int outPoints) {
		double[][] rez=new double[outPoints][2];
		for(int i=0;i<rez.length;i++) {
			rez[i][0]=0;
			rez[i][1]=0;
		}
		int wentThrough=0;
		double curlim=xlim;
		while(curlim<llim) {
			curlim+=xstep;
			wentThrough++;
		}
		llim=commonFunctions.LogBase10(llim);
		rlim=commonFunctions.LogBase10(rlim);
		double lstep=(rlim-llim)/((double)(outPoints-1));
		int used=0;
		while((llim<=rlim)&&(used<outPoints)) {
			double integral=0;
			llim+=lstep;
			while((commonFunctions.LogBase10(curlim)<llim)&&(wentThrough<pdf.length)) {
				curlim+=xstep;
				integral+=pdf[wentThrough];
				wentThrough++;
			}
			if(integral>0) {
				rez[used][0]=llim-0.5*lstep;
				if(used>0) rez[used][1]=commonFunctions.LogBase10(integral/(Math.pow(10,rez[used][0])-Math.pow(10,rez[used-1][0])));
				else rez[used][1]=commonFunctions.LogBase10(integral/(Math.pow(10,rez[used][0])-Math.pow(10,rez[used][0]-lstep)));
				used++;
			}
		}
		if(used<outPoints) {
			double[][] rez2=new double[used][2];
			for(int ii=0;ii<used;ii++) {
				rez2[ii][0]=rez[ii][0];
				rez2[ii][1]=rez[ii][1];
			}
			rez=new double[used][2];
			for(int ii=0;ii<used;ii++) {
				rez[ii][0]=rez2[ii][0];
				rez[ii][1]=rez2[ii][1];
			}
			rez2=null;
		}
		
		return rez;
	}
}

//implementation of numerical solution of the model itself
class doubleFundMood {
	private Random gen=new Random();
	public double lastX=gen.nextDouble();//last nf
	public double lastKsi=2.0*gen.nextDouble()-1.0;//last mood
	private double ecf=1;//individual transition rate c -> f
	private double efc=1;//individual transition rate f -> c
	private double ecc=1;//individual transition rate o -> p and p -> o
	private double H=1;//difference in time scales
	private boolean force=true;//enforce limits?
	private double forcedMin=0.001;//lower limit
	private double forcedMax=1.0-forcedMin;//upper limit
	private double kappa=0.1;//precision
	private double timestep=1;//constant part of the variable time step
	public doubleFundMood(double vecf, double vefc, double vecc, double vH, double vk, double vfmin, double vfmax, boolean vf, int nr) {
		gen=new Random(System.currentTimeMillis()+gen.nextInt(19+nr)+nr*13);
		ecf=vecf;
		efc=vefc;
		ecc=vecc;
		H=vH;
		kappa=vk;
		forcedMin=vfmin;
		forcedMax=vfmax;
		force=vf;
		timestep=kappa*kappa/(1.0+ecf+efc+(1.0+2.0*ecc)*H);
	}
	//discretization of the sollution at given intervals
	public double step(double dt) {
		double t=0;
		while(t<dt) {
			double innerDt=variableTimeStep(lastX,lastKsi);
			if(Double.isNaN(innerDt) || Double.isInfinite(innerDt)) innerDt=dt+dt-t-t;
			double whileDt=Math.min(dt-t,innerDt);
			double oldKsi=lastKsi;
			double oldX=lastX;
			lastKsi=solveSdeKsi(oldX,oldKsi,whileDt);
			lastX=solveSdeNf(oldX,oldKsi,whileDt);
			if(force) {
				lastKsi=Math.max(Math.min(lastKsi,forcedMax),-forcedMax);
				lastX=Math.max(Math.min(lastX,forcedMax),forcedMin);
			}
			t+=whileDt;
		}
		return lastKsi*(1.0-lastX)/lastX;
	}
	//desired timestep for given nf and ksi
	private double variableTimeStep(double nf, double ksi) {
		return timestep*tau(nf,ksi);
	}
	//trading activity scenario
	private double tau(double nf, double ksi) {
		double ret=Math.abs((1.0-nf)/nf*ksi);
		return 1.0/(1.0+0.5*ret*ret);
	}
	// numerical solution of stochastic differential equations
	private double driftNf(double nf, double ksi) {
		return (ecf*(1.0-nf)/tau(nf,ksi)-efc*nf);
	}
	private double diffusionNf(double nf, double ksi) {
		return Math.sqrt(2.0*nf*(1.0-nf)/tau(nf,ksi));
	}
	private double solveSdeNf(double nf, double ksi, double dt) {
		return nf+driftNf(nf,ksi)*dt+diffusionNf(nf,ksi)*Math.sqrt(dt)*gen.nextGaussian();
	}
	private double driftKsi(double nf, double ksi) {
		return -2.0*ksi*ecc*H/tau(nf,ksi);
	}
	private double diffusionKsi(double nf, double ksi) {
		return Math.sqrt(2.0*H*(1.0-ksi*ksi)/tau(nf,ksi));
	}
	private double solveSdeKsi(double nf, double ksi, double dt) {
		return ksi+driftKsi(nf,ksi)*dt+diffusionKsi(nf,ksi)*Math.sqrt(dt)*gen.nextGaussian();
	}
}
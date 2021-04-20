using System;
using Avl.Concerto.DataExplorer;
using System.IO;
using Lomont;
using MathNet.Numerics;
using System.Threading.Tasks;
using System.Globalization;
using System.Runtime.InteropServices;
using System.ComponentModel;
using System.Linq;

namespace NVHLibrary
{

    public struct Param_struct
    {
        public bool CDMTIME;
        public double delta_t;
        public double delta_f;
        public int windowtype;
        public int average;
        public double average_overlap;
        public int LQ;
        public int DC;
        public int y_axis;
        public int diagramtype;
        public int freq_weight;
        public int y_amplitude;
        public double f1;
        public double f2;
        public double n_octave;
        public bool mean;
        public double o1;
        public double o2;
        public double delta_o;
        public int tolerance;
        public string orderstring;
        public double dBref;
        public double y_unitconversion;
        public double dsrange1;
        public double dsrange2;
        public double fs;
        public int fromCA;
        public int toCA;
        public bool calcinverternoise;
        public double inverter_frequency;
        public int polepairs;

        public override bool Equals(object obj)
        {
            throw new NotImplementedException();
        }

        public override int GetHashCode()
        {
            throw new NotImplementedException();
        }

        public static bool operator ==(Param_struct left, Param_struct right)
        {
            return left.Equals(right);
        }

        public static bool operator !=(Param_struct left, Param_struct right)
        {
            return !(left == right);
        }
    }

    public struct XY_Data
    {
        public double[] xdata;
        public double[] ydata;

        public override bool Equals(object obj)
        {
            throw new NotImplementedException();
        }

        public override int GetHashCode()
        {
            throw new NotImplementedException();
        }

        public static bool operator ==(XY_Data left, XY_Data right)
        {
            return left.Equals(right);
        }

        public static bool operator !=(XY_Data left, XY_Data right)
        {
            return !(left == right);
        }
    }

    public struct XYZ_Data
    {
        public double[] xdata;
        public double[] ydata;
        public double[] zdata;
        public double[] xtime;
        public double[] freq_labels;
        public int[] dimensions;
        public double overalllevel;
        public int[] N_nonzero;
        public int rows;
        public string errorstring;
        public Param_struct paramset;
        public string xdescription;

        public override bool Equals(object obj)
        {
            throw new NotImplementedException();
        }

        public override int GetHashCode()
        {
            throw new NotImplementedException();
        }

        public static bool operator ==(XYZ_Data left, XYZ_Data right)
        {
            return left.Equals(right);
        }

        public static bool operator !=(XYZ_Data left, XYZ_Data right)
        {
            return !(left == right);
        }
    }

    internal class TonalityClass : IDisposable
    {
        [DllImport("Tonality.dll", CallingConvention = CallingConvention.Cdecl)]
        private static unsafe extern double Tonality(double* x0, double* y0, double fs, int blocksize, double energy, int wert, double* ton);
        // Wert: 2 wenn Diffusfeld, 1 wenn Freifeld

        internal static unsafe double Tonality_CallDLL(double[] fftreal, double[] fftimag, int N, double fs, double energy, int wert)
        {
            double result = 0;
            fixed (double* p1 = &fftreal[0])
            fixed (double* p2 = &fftimag[0])

            {

                double* p3 = &result;
                Tonality(p1, p2, fs, N, energy, wert, p3);
            }

            return result;
        }

        public void Dispose()
        {
            //Dispose();
            GC.Collect();
            GC.SuppressFinalize(this);
        }
    }

    public class NVHFunctions
    {
       
        internal struct NVHPackage
        {
            internal double[] rawdata_y;
            internal XY_Data n100;
            internal XY_Data enginespeed;
            internal double dt;
            internal double x0;
            internal CATimeStamps timestamps;
            internal Param_struct paramset; 
        }

        

        

        
        internal struct CycleData
        {
            internal double[] xdata;
            internal double[] ydata;
            internal int Count;
            internal int CycleCount;
        }
        
        

        

        internal struct CATimeStamps
        {
            public double[] t1;
            public double[] t2;
            public double[] cycles;
        }

        


        internal struct F_octave
        {
            public double[] f0;
            public double[] f1;
            public double[] f2;
        }

        

        

        

        public static double[] CalcTonality(double[][] audio, double fs, BackgroundWorker worker, double startpercent, double overallprogress)
        {
            if (audio == null) return null;
            int N = 8192;
            double overlap = 0.5;
            //int blocks = (int)((audio[0].Length - N) / (N * (1 - overlap))) + 1;
            int channels = audio.Length;

            double audio_length_seconds = audio[0].Length / fs;

            var fftr = new double[N];
            var ffti = new double[N];
            var oneblock = new double[N];

            LomontFFT Lomontclass = new LomontFFT()
            {
                A = 1
            };

            // Calculate Window Energy
            double energy = 0;
            for (int ii = 0; ii < N; ii++)
            {
                double t = (ii + 0.5) / N;
                double mul = (1 - Math.Cos(2 * Math.PI * t));
                energy += Math.Pow(0.5 * mul, 2);
            }

            var Ton2Ch = new double[channels][];
            var Ton = new XY_Data
            {
                xdata = DeltaArray(N / 2 / fs, (1 - overlap) * N / fs, audio_length_seconds - 0.1)
            };
            int blocks = Ton.xdata.Length;

            using (var TonClass = new TonalityClass())
            {
                for (int ch = 0; ch < channels; ch++)
                {
                    Ton2Ch[ch] = new double[blocks];
                    //Parallel.For(0, blocks, ii =>

                    startpercent += overallprogress * ch / channels;

                    for (int ii = 0; ii < blocks; ii++)
                    {
                        

                        for (int jj = 0; jj < N; jj++) oneblock[jj] = audio[ch][ii * (int)(N * (1 - overlap)) + jj];
                        oneblock = WindowFcn(oneblock, "hann");
                        Lomontclass.RealFFT(oneblock, true);
                        fftr[0] = oneblock[0];
                        fftr[N / 2] = oneblock[1];
                        for (int jj = 1; jj < N / 2; jj++)
                        {
                            fftr[jj] = oneblock[jj * 2];
                            fftr[jj + N / 2] = oneblock[N - jj * 2];
                            ffti[jj] = oneblock[jj * 2 + 1];
                            ffti[jj + N / 2] = -oneblock[N - jj * 2 + 1];
                        }

                        Ton2Ch[ch][ii] = TonalityClass.Tonality_CallDLL(fftr, ffti, N, fs, energy, 2);

                        worker.ReportProgress((int)(startpercent + overallprogress * ii / blocks / channels));
                    }

                    // Avoid Peaks in Tonality
                    double old = Ton2Ch[ch][0];
                    for (int i = 1; i < blocks - 1; i++)
                    {
                        if (Ton2Ch[ch][i] - old > 0.1) Ton2Ch[ch][i] = old + 0.1;
                        if (old - Ton2Ch[ch][i] > 0.1) Ton2Ch[ch][i] = old - 0.1;
                        old = Ton2Ch[ch][i];
                    }

                }
            }
              
            
            var Ton_mean = new double[Ton2Ch[0].Length];
            for (int ii = 0; ii < Ton_mean.Length; ii++)
            {
                double summe = 0;
                for (int ch = 0; ch < channels; ch++) summe += Ton2Ch[ch][ii];
                Ton_mean[ii] = summe / channels;
            }
            
            Ton.ydata = Interp1_linear(Ton.xdata, Ton_mean, DeltaArray(0.1, 0.025, audio_length_seconds - 0.1));
            Ton.xdata = DeltaArray(0.1, 0.025, audio_length_seconds - 0.1);
            ExtrapolateNaNs(Ton);
            int smoothingpoints = 5;
            var ton_smooth = Moving(Ton.ydata, 5);
            for (int ii = smoothingpoints / 2; ii < ton_smooth.Length - smoothingpoints / 2; ii++) Ton.ydata[ii] = ton_smooth[ii];

            return Ton.ydata;
            
        }


        internal void SaveSomeDoubles(string filename, params double[] doubles)
        {
            FileStream output;
            BinaryWriter binWrt;

            output = new FileStream(filename, FileMode.Create);
            binWrt = new BinaryWriter(output);

            binWrt.Write(doubles.Length);
            for (int ii = 0; ii < doubles.Length; ii++) binWrt.Write(doubles[ii]);
            binWrt.Close();
        }

        internal double[] LoadSomeDoubles(string filename)
        {

            BinaryReader br;
            double[] values;

            //reading from the file
            try
            {
                br = new BinaryReader(new FileStream(filename, FileMode.Open));
            }
            catch (IOException e)
            {
                Console.WriteLine(e.Message + "\n Cannot open file.");
                values = Array.Empty<double>();
                return values;
            }

            try
            {
                int length = br.ReadInt32();

                values = new double[length];
                for (int ii = 0; ii < values.Length; ii++) values[ii] = br.ReadDouble();


            }
            catch (IOException e)
            {
                Console.WriteLine(e.Message + "\n Cannot read from file.");
                values = Array.Empty<double>();
                return values;
            }
            br.Close();
            return values;
        }



        internal void SaveSomeIntegers(string filename, params int[] integers)
        {
            FileStream output;
            BinaryWriter binWrt;

            output = new FileStream(filename, FileMode.Create);
            binWrt = new BinaryWriter(output);

            binWrt.Write(integers.Length);
            for (int ii = 0; ii < integers.Length; ii++) binWrt.Write(integers[ii]);
            binWrt.Close();
        }

        internal int[] LoadSomeIntegers(string filename)
        {

            BinaryReader br;
            int[] values;

            //reading from the file
            try
            {
                br = new BinaryReader(new FileStream(filename, FileMode.Open));
            }
            catch (IOException e)
            {
                Console.WriteLine(e.Message + "\n Cannot open file.");
                values = Array.Empty<int>();
                return values;
            }

            try
            {
                int length = br.ReadInt32();

                values = new int[length];
                for (int ii = 0; ii < values.Length; ii++) values[ii] = br.ReadInt32();


            }
            catch (IOException e)
            {
                Console.WriteLine(e.Message + "\n Cannot read from file.");
                values = Array.Empty<int>();
                return values;
            }
            br.Close();
            return values;
        }



        internal void SaveXYData(string filename, XY_Data data)
        {
            FileStream output;
            BinaryWriter binWrt;

            output = new FileStream(filename, FileMode.Create);
            binWrt = new BinaryWriter(output);

            binWrt.Write(data.xdata.Length);
            for (int ii = 0; ii < data.xdata.Length; ii++) binWrt.Write(data.xdata[ii]);
            for (int ii = 0; ii < data.ydata.Length; ii++) binWrt.Write(data.ydata[ii]);
            binWrt.Close();
        }

        internal XY_Data LoadXYData(string filename)
        {

            BinaryReader br;
            XY_Data values;

            //reading from the file
            try
            {
                br = new BinaryReader(new FileStream(filename, FileMode.Open));
            }
            catch (IOException e)
            {
                Console.WriteLine(e.Message + "\n Cannot open file.");
                values = new XY_Data
                {
                    xdata = Array.Empty<double>(),
                    ydata = Array.Empty<double>()
                };
                return values;
            }

            try
            {
                int length = br.ReadInt32();

                values = new XY_Data
                {
                    xdata = new double[length],
                    ydata = new double[length]
                };
                for (int ii = 0; ii < length; ii++) values.xdata[ii] = br.ReadDouble();
                for (int ii = 0; ii < length; ii++) values.ydata[ii] = br.ReadDouble();


            }
            catch (IOException e)
            {
                Console.WriteLine(e.Message + "\n Cannot read from file.");
                values = new XY_Data
                {
                    xdata = Array.Empty<double>(),
                    ydata = Array.Empty<double>()
                };
                return values;
            }
            br.Close();
            return values;
        }


        internal void SaveCycleData(string filename, CycleData data)
        {
            FileStream output;
            BinaryWriter binWrt;

            output = new FileStream(filename, FileMode.Create);
            binWrt = new BinaryWriter(output);

            binWrt.Write(data.Count);
            binWrt.Write(data.CycleCount);
            for (int ii = 0; ii < data.Count; ii++) binWrt.Write(data.xdata[ii]);
            for (int ii = 0; ii < data.Count * data.CycleCount; ii++) binWrt.Write(data.ydata[ii]);
            binWrt.Close();
        }

        internal CycleData LoadCycleData(string filename)
        {

            BinaryReader br;
            CycleData values;

            //reading from the file
            try
            {
                br = new BinaryReader(new FileStream(filename, FileMode.Open));
            }
            catch (IOException e)
            {
                Console.WriteLine(e.Message + "\n Cannot open file.");
                values = new CycleData
                {
                    Count = 0,
                    CycleCount = 0,
                    xdata = Array.Empty<double>(),
                    ydata = Array.Empty<double>()
                };
                return values;
            }

            try
            {
                int count = br.ReadInt32();
                int cyclecount = br.ReadInt32();

                values = new CycleData
                {
                    Count = count,
                    CycleCount = cyclecount,  
                    xdata = new double[count],
                    ydata = new double[count * cyclecount]
                };
                for (int ii = 0; ii < count; ii++) values.xdata[ii] = br.ReadDouble();
                for (int ii = 0; ii < count * cyclecount; ii++) values.ydata[ii] = br.ReadDouble();


            }
            catch (IOException e)
            {
                Console.WriteLine(e.Message + "\n Cannot read from file.");
                values = new CycleData
                {
                    Count = 0,
                    CycleCount = 0,
                    xdata = Array.Empty<double>(),
                    ydata = Array.Empty<double>()
                };
                return values;
            }
            br.Close();
            return values;
        }


        internal static void SaveAPSInput(string filename, NVHPackage package)
        {
            FileStream output;
            BinaryWriter binWrt;

            var rawdata_y = package.rawdata_y;
            var n100 = package.n100;
            var enginespeed = package.enginespeed;
            var dt = package.dt;
            var x0 = package.x0;
            var timestamps = package.timestamps;
            var paramset = package.paramset;


            output = new FileStream(filename, FileMode.Create);
            binWrt = new BinaryWriter(output);

            // rawdata
            binWrt.Write(rawdata_y.Length);
            for (int i = 0; i < rawdata_y.Length; i++) binWrt.Write(BitConverter.GetBytes(rawdata_y[i]));

            // n100
            binWrt.Write(n100.xdata.Length);
            for (int i = 0; i < n100.xdata.Length; i++) binWrt.Write(BitConverter.GetBytes(n100.xdata[i]));
            binWrt.Write(n100.ydata.Length);
            for (int i = 0; i < n100.ydata.Length; i++) binWrt.Write(BitConverter.GetBytes(n100.ydata[i]));

            // enginespeed
            binWrt.Write(enginespeed.xdata.Length);
            for (int i = 0; i < enginespeed.xdata.Length; i++) binWrt.Write(BitConverter.GetBytes(enginespeed.xdata[i]));
            binWrt.Write(enginespeed.ydata.Length);
            for (int i = 0; i < enginespeed.ydata.Length; i++) binWrt.Write(BitConverter.GetBytes(enginespeed.ydata[i]));

            // dt
            binWrt.Write(BitConverter.GetBytes(dt));

            // x0
            binWrt.Write(BitConverter.GetBytes(x0));

            // timestamps
            binWrt.Write(timestamps.t1.Length);
            for (int i = 0; i < timestamps.t1.Length; i++) binWrt.Write(BitConverter.GetBytes(timestamps.t1[i]));
            binWrt.Write(timestamps.t2.Length);
            for (int i = 0; i < timestamps.t2.Length; i++) binWrt.Write(BitConverter.GetBytes(timestamps.t2[i]));
            binWrt.Write(timestamps.cycles.Length);
            for (int i = 0; i < timestamps.cycles.Length; i++) binWrt.Write(BitConverter.GetBytes(timestamps.cycles[i]));

            // params
            binWrt.Write(BitConverter.GetBytes(paramset.CDMTIME));
            binWrt.Write(BitConverter.GetBytes(paramset.delta_t));
            binWrt.Write(BitConverter.GetBytes(paramset.delta_f));
            binWrt.Write(BitConverter.GetBytes(paramset.windowtype));
            binWrt.Write(BitConverter.GetBytes(paramset.average));
            binWrt.Write(BitConverter.GetBytes(paramset.average_overlap));
            binWrt.Write(BitConverter.GetBytes(paramset.LQ));
            binWrt.Write(BitConverter.GetBytes(paramset.DC));
            binWrt.Write(BitConverter.GetBytes(paramset.y_axis));
            binWrt.Write(BitConverter.GetBytes(paramset.diagramtype));
            binWrt.Write(BitConverter.GetBytes(paramset.freq_weight));
            binWrt.Write(BitConverter.GetBytes(paramset.y_amplitude));
            binWrt.Write(BitConverter.GetBytes(paramset.f1));
            binWrt.Write(BitConverter.GetBytes(paramset.f2));
            binWrt.Write(BitConverter.GetBytes(paramset.n_octave));
            binWrt.Write(BitConverter.GetBytes(paramset.mean));
            binWrt.Write(BitConverter.GetBytes(paramset.o1));
            binWrt.Write(BitConverter.GetBytes(paramset.o2));
            binWrt.Write(BitConverter.GetBytes(paramset.delta_o));
            binWrt.Write(BitConverter.GetBytes(paramset.tolerance));
            binWrt.Write(BitConverter.GetBytes(paramset.dBref));
            binWrt.Write(BitConverter.GetBytes(paramset.y_unitconversion));
            binWrt.Write(BitConverter.GetBytes(paramset.dsrange1));
            binWrt.Write(BitConverter.GetBytes(paramset.dsrange2));
            binWrt.Write(BitConverter.GetBytes(paramset.fs));

            output.Close();

        }

        
        internal NVHPackage LoadAPSInput(string filename)
        {
             // Read full set of Data input from binary file (rawdata, leading quantity, params etc.)
            BinaryReader br;
            var result = new NVHPackage();

            //reading from the file
            try
            {
                br = new BinaryReader(new FileStream(filename, FileMode.Open));
            }
            catch (IOException e)
            {
                Console.WriteLine(e.Message + "\n Cannot open file.");
                result = new NVHPackage();
                return result;
            }

            try
            {
                // Rawdata
                byte[] bytes = br.ReadBytes(br.ReadInt32() * 8);
                result.rawdata_y = new double[bytes.Length / 8];
                for (int ii = 0; ii < result.rawdata_y.Length; ii++) result.rawdata_y[ii] = BitConverter.ToDouble(bytes, ii * 8);
                // n100
                bytes = br.ReadBytes(br.ReadInt32() * 8);
                result.n100.xdata = new double[bytes.Length / 8];
                for (int ii = 0; ii < result.n100.xdata.Length; ii++) result.n100.xdata[ii] = BitConverter.ToDouble(bytes, ii * 8);
                bytes = br.ReadBytes(br.ReadInt32() * 8);
                result.n100.ydata = new double[bytes.Length / 8];
                for (int ii = 0; ii < result.n100.ydata.Length; ii++) result.n100.ydata[ii] = BitConverter.ToDouble(bytes, ii * 8);
                // enginespeed
                bytes = br.ReadBytes(br.ReadInt32() * 8);
                result.enginespeed.xdata = new double[bytes.Length / 8];
                for (int ii = 0; ii < result.enginespeed.xdata.Length; ii++) result.enginespeed.xdata[ii] = BitConverter.ToDouble(bytes, ii * 8);
                bytes = br.ReadBytes(br.ReadInt32() * 8);
                result.enginespeed.ydata = new double[bytes.Length / 8];
                for (int ii = 0; ii < result.enginespeed.ydata.Length; ii++) result.enginespeed.ydata[ii] = BitConverter.ToDouble(bytes, ii * 8);
                // dt
                bytes = br.ReadBytes(8);
                result.dt = BitConverter.ToDouble(bytes, 0);
                // x0
                bytes = br.ReadBytes(8);
                result.x0 = BitConverter.ToDouble(bytes, 0);
                // enginespeed
                bytes = br.ReadBytes(br.ReadInt32() * 8);
                result.timestamps.t1 = new double[bytes.Length / 8];
                for (int ii = 0; ii < result.timestamps.t1.Length; ii++) result.timestamps.t1[ii] = BitConverter.ToDouble(bytes, ii * 8);
                bytes = br.ReadBytes(br.ReadInt32() * 8);
                result.timestamps.t2 = new double[bytes.Length / 8];
                for (int ii = 0; ii < result.timestamps.t2.Length; ii++) result.timestamps.t2[ii] = BitConverter.ToDouble(bytes, ii * 8);
                bytes = br.ReadBytes(br.ReadInt32() * 8);
                result.timestamps.cycles = new double[bytes.Length / 8];
                for (int ii = 0; ii < result.timestamps.cycles.Length; ii++) result.timestamps.cycles[ii] = BitConverter.ToDouble(bytes, ii * 8);
                // params
                bytes = br.ReadBytes(1);
                result.paramset.CDMTIME = BitConverter.ToBoolean(bytes, 0);
                bytes = br.ReadBytes(8);
                result.paramset.delta_t = BitConverter.ToDouble(bytes, 0);
                bytes = br.ReadBytes(8);
                result.paramset.delta_f = BitConverter.ToDouble(bytes, 0);
                bytes = br.ReadBytes(4);
                result.paramset.windowtype = BitConverter.ToInt32(bytes, 0);
                bytes = br.ReadBytes(4);
                result.paramset.average = BitConverter.ToInt32(bytes, 0);
                bytes = br.ReadBytes(8);
                result.paramset.average_overlap = BitConverter.ToDouble(bytes, 0);
                bytes = br.ReadBytes(4);
                result.paramset.LQ = BitConverter.ToInt32(bytes, 0);
                bytes = br.ReadBytes(4);
                result.paramset.DC = BitConverter.ToInt32(bytes, 0);
                bytes = br.ReadBytes(4);
                result.paramset.y_axis = BitConverter.ToInt32(bytes, 0);
                bytes = br.ReadBytes(4);
                result.paramset.diagramtype = BitConverter.ToInt32(bytes, 0);
                bytes = br.ReadBytes(4);
                result.paramset.freq_weight = BitConverter.ToInt32(bytes, 0);
                bytes = br.ReadBytes(4);
                result.paramset.y_amplitude = BitConverter.ToInt32(bytes, 0);
                bytes = br.ReadBytes(8);
                result.paramset.f1 = BitConverter.ToDouble(bytes, 0);
                bytes = br.ReadBytes(8);
                result.paramset.f2 = BitConverter.ToDouble(bytes, 0);
                bytes = br.ReadBytes(8);
                result.paramset.n_octave = BitConverter.ToDouble(bytes, 0);
                bytes = br.ReadBytes(1);
                result.paramset.mean = BitConverter.ToBoolean(bytes, 0);
                bytes = br.ReadBytes(8);
                result.paramset.o1 = BitConverter.ToDouble(bytes, 0);
                bytes = br.ReadBytes(8);
                result.paramset.o2 = BitConverter.ToDouble(bytes, 0);
                bytes = br.ReadBytes(8);
                result.paramset.delta_o = BitConverter.ToDouble(bytes, 0);
                bytes = br.ReadBytes(4);
                result.paramset.tolerance = BitConverter.ToInt32(bytes, 0);
                bytes = br.ReadBytes(8);
                result.paramset.dBref = BitConverter.ToDouble(bytes, 0);
                bytes = br.ReadBytes(8);
                result.paramset.y_unitconversion = BitConverter.ToDouble(bytes, 0);
                bytes = br.ReadBytes(8);
                result.paramset.dsrange1 = BitConverter.ToDouble(bytes, 0);
                bytes = br.ReadBytes(8);
                result.paramset.dsrange2 = BitConverter.ToDouble(bytes, 0);
                bytes = br.ReadBytes(8);
                result.paramset.fs = BitConverter.ToDouble(bytes, 0);


            }
            catch (IOException e)
            {
                Console.WriteLine(e.Message + "\n Cannot read from file.");
                result = new NVHPackage();
                return result;
            }
            br.Close();
            return result;


        }

        internal int[] GetCycleRange(double[] cyctime, double timerange1, double timerange2)
        {
            var cyclerange = new int[2];
            cyclerange[0] = -2;
            cyclerange[1] = -2;
            for (int ii = 0; ii < cyctime.Length; ii++)
            {
                if (cyctime[ii] > timerange1)
                {
                    cyclerange[0] = ii - 1;
                    break;
                }
            }
            for (int ii = 0; ii < cyctime.Length; ii++)
            {
                if (cyctime[ii] >= timerange2)
                {
                    cyclerange[1] = ii - 1;
                    break;
                }
            }

            if (cyclerange[0] == -1) cyclerange[0] = 0;
            if (cyclerange[1] == -1) cyclerange[1] = 0;
            if (cyclerange[0] == -2) cyclerange[0] = cyctime.Length - 1;
            if (cyclerange[1] == -2) cyclerange[1] = cyctime.Length - 1;
            return cyclerange;
        }


        private static double[] GetMdi(CycleData pcyl, double[] tdcoffset, int[] pinofftype, double bore, double conrod, double stroke, double pinoff, double pinoff2)
        {

            

            int cycleindex = pcyl.CycleCount - 1;
            int cycles = cycleindex + 1;
            int count = pcyl.Count;

            bore /= 1000;
            stroke /= 1000;
            conrod /= 1000;

            var mdi = new XY_Data
            {
                xdata = new double[pcyl.ydata.Length],
                ydata = new double[pcyl.ydata.Length]
            };

            var cyclesDS = DeltaArray(0, 1, cycleindex);

            int cylcount = pcyl.ydata.Length / pcyl.Count / pcyl.CycleCount;
            
            Parallel.For(0,cylcount, cyl =>
            {
                var phi = pcyl.xdata;

                var tdcoffset_cyl = tdcoffset[cyl];

                double pinoff_cyl = 0;
                if (pinofftype[cyl] == 1) pinoff_cyl = pinoff / 1000;
                else if (pinofftype[cyl] == -1) pinoff_cyl = -pinoff / 1000;
                else if (pinofftype[cyl] == 2) pinoff_cyl = pinoff2 / 1000;
                else if (pinofftype[cyl] == -2) pinoff_cyl = -pinoff2 / 1000;

                var alpha = new double[phi.Length];
                for (int ii = 0; ii < alpha.Length; ii++) alpha[ii] = Math.Asin((stroke / 2 * Math.Sin(phi[ii] * Math.PI / 180) - pinoff_cyl) / conrod);

                for (int c = 0; c < cycles; c++)
                {
                    var pcylc = new double[pcyl.Count];
                    for (int ii = 0; ii < pcyl.Count; ii++) pcylc[ii] = pcyl.ydata[ii + c * count + cyl * count * cycles] * 1000;
                    var fgas = new double[pcylc.Length];
                    for (int ii = 0; ii < fgas.Length; ii++) fgas[ii] = pcylc[ii] * (Math.Pow(bore, 2) * Math.PI / 4);

                    for (int ii = 0; ii < count; ii++) mdi.ydata[ii + cyl * count * cycles + c * count] = fgas[ii] / Math.Cos(alpha[ii]) * Math.Sin(phi[ii] * Math.PI / 180 + alpha[ii]) * stroke / 2;
                    for (int ii = 0; ii < count; ii++) mdi.xdata[ii + cyl * count * cycles + c * count] = phi[ii] + tdcoffset_cyl + 720 * c;

                }
            });
            
            
            // 

            var fromDS = new double[cylcount];
            for (int ii = 0; ii < cylcount; ii++) fromDS[ii] = mdi.xdata[ii * count * cycles];
            var toDS = new double[cylcount];
            for (int ii = 0; ii < cylcount; ii++) toDS[ii] = mdi.xdata[(ii + 1) * count * cycles - 1];
            double from = fromDS[0];
            for (int ii = 0; ii < fromDS.Length; ii++) if (fromDS[ii] > from) from = fromDS[ii];
            double to = toDS[0];
            for (int ii = 0; ii < toDS.Length; ii++) if (toDS[ii] < to) to = toDS[ii];

            var mdisum = new double[count * (cycles - 1)];

            var ii1 = new int[cylcount];
            // find index of from
            for (int cyl = 0; cyl < cylcount; cyl++)
            {
                for (int ii = 0; ii < count * cycles; ii++)
                {
                    if (mdi.xdata[ii + cyl * count * cycles] >= from)
                    {
                        ii1[cyl] = ii;
                        break;
                    }
                }
            }
            
            for (int cyl = 0; cyl < cylcount; cyl++)
            {
                for (int ii = 0; ii < mdisum.Length; ii++)
                {
                    mdisum[ii] += mdi.ydata[ii1[cyl] + cyl * cycles * count + ii];
                }
            }
            

            //var mdi = new XY_Data();
            return mdisum;
        }

        internal XY_Data GetHSI(CycleData pcyl, double[] tdcoffset, int[] pinofftype, double bore, double conrod, double stroke, double pinoff, double pinoff2, double[] cyctime)
        {
            var hsi = new XY_Data();

            if (cyctime[cyctime.Length - 1] - cyctime[0] < 1.101)
            {
                hsi.xdata = new double[cyctime.Length];
                hsi.ydata = new double[cyctime.Length];
                for (int ii = 0; ii < hsi.xdata.Length; ii++) hsi.xdata[ii] = cyctime[ii] - 0.5;
                return hsi;
            }
            
            var mdi = GetMdi(pcyl, tdcoffset, pinofftype, bore, conrod, stroke, pinoff, pinoff2);
            
            
            int fs_neu = 4096;
            

            var d_cyctime = new double[cyctime.Length - 1];
            for (int jj = 0; jj < d_cyctime.Length; jj++) d_cyctime[jj] = cyctime[jj + 1] - cyctime[jj];

            var seconds = new double[mdi.Length];
            for (int jj = 0; jj < cyctime.Length - 1; jj++)
            {
                seconds[jj * 720] = cyctime[jj];
                for (int kk = 1; kk < 720; kk++) seconds[jj * 720 + kk] = seconds[jj * 720 + kk - 1] + d_cyctime[jj] / 720;
            }

            /*
            // Low pass filter (Anti Aliasing)
            double maxorder = 16;
            int N_filter = 42;
            double fc1 = 2 * maxorder / 360;
            var b = MathNet.Filtering.FIR.FirCoefficients.LowPass(360, maxorder, N_filter / 2);
            b = WindowFcn(b, "hamming");
            var FilterObj = new MathNet.Filtering.FIR.OnlineFirFilter(b);
            var filtered_data_temp = FilterObj.ProcessSamples(rawdata_singleCAOcycle);
            var filtered_data = new double[filtered_data_temp.Length - N_filter / 2];
            Array.Resize(ref seconds, seconds.Length - N_filter / 2);
            for (int jj = 0; jj < filtered_data.Length; jj++) filtered_data[jj] = filtered_data_temp[jj + N_filter / 2];
            */


            var timedata_x = DeltaArray(seconds[0], 1.0 / fs_neu, seconds[seconds.Length - 1]);
            var timedata_y = Interp1_linear(seconds, mdi, timedata_x);


            double dt = timedata_x[1] - timedata_x[0];
            double x0 = timedata_x[0];

            var paramset = new Param_struct
            {
                diagramtype = 5,
                f1 = 3,
                f2 = 7,
                delta_f = 1,
                delta_t = 0.1,
                average = 1,
                average_overlap = 0,
                fs = 4096,
                LQ = 0,
                CDMTIME = false,
                windowtype = 1,
                y_amplitude = 1,
                y_axis = 0,
                y_unitconversion = 1,
                freq_weight = 0,
                dsrange1 = timedata_x[0],
                dsrange2 = timedata_x[timedata_x.Length - 1]
            };
            //var enginespeed = new XY_Data();
            var n100 = new XY_Data();
            var timestamps = new CATimeStamps();

            var result_xyz = CalcBPLevel(timedata_y, n100, dt, x0, timestamps, paramset);
            var result = new XY_Data
            {
                xdata = result_xyz.xdata,
                ydata = result_xyz.zdata
            };
            for (int ii = 0; ii < result.ydata.Length; ii++) result.ydata[ii] *= 12.5;

            return result;
        }

        internal static int[] Pulses2PulseIndices(bool[] pulses)
        {
            int L = 0;
            for (int ii = 0; ii < pulses.Length; ii++) if (pulses[ii]) L++;
            var result = new int[L];
            int jj = 0;
            for (int ii = 0; ii < pulses.Length; ii++)
            {
                if (pulses[ii])
                {
                    result[jj] = ii;
                    jj++;
                }
            }
            return result;
        }

        internal int[] RPM2PulseIndices(XY_Data rpm, int L_raw, int fs_raw)
        {
            var result = new int[L_raw];
            var samplediff = new double[rpm.ydata.Length];
            var x_samples = new int[rpm.xdata.Length];
            for (int ii = 0; ii < samplediff.Length; ii++)
            {
                samplediff[ii] = (double)fs_raw * 60 / rpm.ydata[ii];
                x_samples[ii] = (int)Math.Round((rpm.xdata[ii] - rpm.xdata[0]) * fs_raw, MidpointRounding.AwayFromZero);
            }
            result[0] = 0;
            int i = 0;
            int k = 1;
            int j;
            for (j = 1; j < x_samples.Length; j++)
            {
                while (i + Math.Round(samplediff[j - 1],MidpointRounding.AwayFromZero) - 2 < x_samples[j])
                {
                    result[k] = i + (int)Math.Round(samplediff[j - 1], MidpointRounding.AwayFromZero);
                    k++;
                    i += (int)Math.Round(samplediff[j - 1], MidpointRounding.AwayFromZero);
                }
            }
            j--;
            while (i + (int)Math.Round(samplediff[j],MidpointRounding.AwayFromZero) < result.Length)
            {
                result[k] = i + (int)Math.Round(samplediff[j],MidpointRounding.AwayFromZero);
                k++;
                i += (int)Math.Round(samplediff[j],MidpointRounding.AwayFromZero);
            }
            Array.Resize(ref result, k);

            return result;
        }

        internal static int[] RPMVel2PulseIndices(XY_Data rpm_vel, int L_raw, int fs_raw, bool IsVehicleSpeed)
        {
            var result = new int[L_raw];

            // Engine speed => Pulses = 60, Factor = 1
            // Vehicle speed => Pulses = 90, Factor = 16.667

            
            double rpmvelfactor;
            int pulses;
            double minpulsedistance;
            

            if (IsVehicleSpeed)
            {
                rpmvelfactor = 3.6;
                pulses = 90;
                minpulsedistance = 0.04;
                
            }
            else
            {
                rpmvelfactor = 60;
                pulses = 60;
                minpulsedistance = 0.0025;
                
            }
           
            // Time distance between two pulses
            var t_pulses = new double[rpm_vel.ydata.Length];
            for (int ii = 0; ii < t_pulses.Length; ii++) t_pulses[ii] = rpmvelfactor / rpm_vel.ydata[ii] / pulses;

            var t_minres = DeltaArray(0, 1 / (double)fs_raw, (L_raw - 1) / (double)fs_raw);
            var t_pulses_minres = Interp1_linear(rpm_vel.xdata, t_pulses, t_minres);
            for (int ii = 0; ii < t_pulses_minres.Length; ii++) if (t_pulses_minres[ii] > minpulsedistance) t_pulses_minres[ii] = Double.NaN;

            int L_indices = 0;
            //var pulsechannel = new bool[t_minres.Length];
            int jj = Array.FindIndex(t_pulses_minres, f => !double.IsNaN(f));
            double t_akt = t_minres[jj];
            while (jj < t_pulses_minres.Length)
            {
                result[L_indices] = jj;
                L_indices++;
                t_akt += t_pulses_minres[jj];
                if (double.IsNaN(t_akt))
                {
                    jj = Array.FindIndex(t_pulses_minres, jj, f => !double.IsNaN(f));
                    if (jj == -1) break;
                    t_akt = t_minres[jj];
                }
                else
                {
                    jj = (int)Math.Round(t_akt * fs_raw, MidpointRounding.AwayFromZero);
                }
            }

            Array.Resize(ref result, L_indices);
            return result;
        }

        internal static int[] FindPulses(double[] rawdata, int fs)
        {
            // Define trigger level
            
            var trigger = new double[rawdata.Length / fs + 1];
            double maxval, minval;
            int jj;
            for (int ii = 0; ii < trigger.Length - 1; ii++)
            {
                maxval = rawdata[ii * fs];
                minval = rawdata[ii * fs];
                for (jj = ii * fs; jj < (ii + 1) * fs; jj++)
                {
                    if (rawdata[jj] > maxval) maxval = rawdata[jj];
                    if (rawdata[jj] < minval) minval = rawdata[jj];
                }
                trigger[ii] = (maxval - minval) / 3;
            }
            trigger[trigger.Length - 1] = trigger[trigger.Length - 2];
            

            // Get Number of Trigger transgressions
            bool above;
            if (rawdata[0] > trigger[0]) above = true;
            else above = false;
            int L_Peaks = 0;
            for (int ii = 1; ii < rawdata.Length; ii++)
            {
                if (rawdata[ii] > trigger[ii / fs])
                {
                    if (!above) L_Peaks++;
                    above = true;
                }
                else above = false;
            }
            if (above) L_Peaks--;

            // Get positive and negative trigger transgression indices
            var positive = new int[L_Peaks];
            var negative = new int[L_Peaks];
            if (rawdata[0] > trigger[0]) above = true;
            else above = false;
            jj = 0;
            for (int ii = 1; ii < rawdata.Length; ii++)
            {
                if (rawdata[ii] > trigger[ii / fs])
                {
                    if (!above && jj < positive.Length) positive[jj] = ii;
                    above = true;
                }
                else
                {
                    if (above)
                    {
                        negative[jj] = ii;
                        if (positive[0] != 0) jj++;
                    }
                    above = false;
                }
            }

            // Get Peak indices
            var OT = new int[positive.Length];
            for (int ii = 0; ii < OT.Length; ii++)
            {
                OT[ii] = positive[ii];
                maxval = rawdata[positive[ii]];
                for (jj = positive[ii] + 1; jj <= negative[ii]; jj++)
                {
                    if (rawdata[jj] > maxval)
                    {
                        maxval = rawdata[jj];
                        OT[ii] = jj;
                    }
                }
            }
            
            var output = OT;
            return output;

        }
        

        internal static CATimeStamps GetCATimeStamps_CDMTIME(CycleData cdmtime_cyc, double fromCA, double toCA, double delta_t, int ms)
        {
            
            int fromCA_i = 0;
            int toCA_i = cdmtime_cyc.Count - 1;
            if (ms == 1) for (int ii = 0; ii < cdmtime_cyc.ydata.Length; ii++) cdmtime_cyc.ydata[ii] /= 1000;
            
            for (int ii = 0; ii < cdmtime_cyc.xdata.Length; ii++)
            {
                if (cdmtime_cyc.xdata[ii] >= fromCA)
                {
                    fromCA_i = ii;
                    break;
                }
            }
            for (int ii = cdmtime_cyc.xdata.Length - 1; ii >= 0; ii--)
            {
                if (cdmtime_cyc.xdata[ii] <= toCA)
                {
                    toCA_i = ii;
                    break;
                }
            }

            int L_lq = 0;
            int ii1 = fromCA_i;
            while (ii1 < cdmtime_cyc.CycleCount * cdmtime_cyc.Count)
            {
                ii1 += (int)(delta_t * cdmtime_cyc.Count);
                L_lq++;
            }
           
            var result = new CATimeStamps
            {
                t1 = new double[L_lq],
                t2 = new double[L_lq],
                cycles = new double[L_lq]

            };


            for (int ii = 0; ii < result.t1.Length; ii++)
            {
                result.t1[ii] = cdmtime_cyc.ydata[fromCA_i + (int)(ii * delta_t * cdmtime_cyc.Count)];
                result.t2[ii] = cdmtime_cyc.ydata[toCA_i + (int)(ii * delta_t * cdmtime_cyc.Count)];
                //result.cycles[ii] = ii * delta_t + 1 + cyclerange1;
            }

            return result;
        }

        public static double[] WindowFcn(double[] input, string windowtype)
        {
            if (Equals(windowtype, "hann")) for (int ii = 0; ii < input.Length; ii++) input[ii] *= 0.5 * (1 - Math.Cos(2 * Math.PI * ii / (input.Length - 1)));
            else if(Equals(windowtype, "hamming")) for (int ii = 0; ii < input.Length; ii++) input[ii] *= 0.54 - 0.46 * Math.Cos(2 * Math.PI * ii / (input.Length - 1));
            return input;
        }

        private static CycleData InterpolateUnevenCycleData(CycleData input, double minres)
        {
            
      
            var result = new CycleData
            {
                xdata = DeltaArray(-360, minres, 360 - minres),
                ydata = new double[(int)(720 / minres * input.CycleCount)],
                CycleCount = input.CycleCount,
                Count = input.Count
            };

            
            var xdata_long = VertcatCycleData(input);
            double[] result_xdata = DeltaArray(-360, minres, input.CycleCount * 720 - minres - 360);

            result.ydata = Interp1_linear(xdata_long, input.ydata, result_xdata);
            return result;
        }

        private static double[] VertcatCycleData(CycleData input)
        {
            
            var result = new double[input.Count * input.CycleCount];
            for (int ii = 0; ii < input.CycleCount; ii++)
            {
                for (int jj = 0; jj < input.Count; jj++)
                {
                    result[ii * input.Count + jj] = input.xdata[jj] + ii * 720;
                }
            }
            return result;
        }
        

        internal static double[] ResampleTMData(double[] rawdata,double x0, double dt, double fs_new)
        {
            
            double lastx = x0 + dt * (rawdata.Length - 1);
            double[] old_timeaxis = DeltaArray(x0, dt, lastx);
            
            double[] IP_input;

            if (fs_new < (1 / dt))
            {
                int N = 42;
                IP_input = new double[rawdata.Length - N / 2];
                double fc = fs_new / 2;
                
                double[] b = MathNet.Filtering.FIR.FirCoefficients.LowPass(fs_new, fc, N / 2);
                b = WindowFcn(b, "hamming");

                var FilterObj = new MathNet.Filtering.FIR.OnlineFirFilter(b);
                double[] filtered_data = FilterObj.ProcessSamples(rawdata);
                
                for (int jj = 0; jj < rawdata.Length - N / 2; jj++) IP_input[jj] = filtered_data[jj + N / 2];
                Array.Resize(ref old_timeaxis, old_timeaxis.Length - N / 2);
            }
            else
            IP_input = rawdata;
            

            double[] new_timeaxis = DeltaArray(x0, 1 / fs_new, old_timeaxis[old_timeaxis.Length - 1]);

            double[] result = Interp1_spline(old_timeaxis, IP_input, new_timeaxis);
            return result;
        }

        internal static double[] CDMTIME2CycRPM(CycleData cdmtime)
        {
            var cycrpm = new double[cdmtime.CycleCount];
            for (int ii = 0; ii < cycrpm.Length - 1; ii++) cycrpm[ii] = 2 * 1000 * 60 / (cdmtime.ydata[(ii + 1) * cdmtime.Count] - cdmtime.ydata[ii * cdmtime.Count]);
            cycrpm[cycrpm.Length - 1] = 2 * 1000 * 60 / ((cdmtime.ydata[cdmtime.ydata.Length - 1] - cdmtime.ydata[(cdmtime.CycleCount - 1) * cdmtime.Count]) * (1 + (cdmtime.xdata[0] + 720 - cdmtime.xdata[cdmtime.Count - 1]) / (cdmtime.xdata[cdmtime.Count - 1] - cdmtime.xdata[0])));
            return cycrpm;
        }

        internal static CycleData DT2CDMTIME(CycleData dT, double[] cyctime, int ms)
        {
            double time_factor = 1;
            if (ms == 0) time_factor = 1000;
            var dt_dT = new double[dT.Count];
            for (int ii = 1; ii < dt_dT.Length; ii++) dt_dT[ii - 1] = dT.xdata[ii] - dT.xdata[ii - 1];
            dt_dT[dt_dT.Length - 1] = dT.xdata[0] + 720 - dT.xdata[dT.Count - 1];
            var cdmtime = new CycleData
            {
                Count = dT.Count,
                CycleCount = dT.CycleCount,
                xdata = dT.xdata,
                ydata = new double[dT.ydata.Length]
            };
            //var cdmtime = new double[dT.ydata.Length];
            for (int ii = 0; ii < dT.CycleCount; ii++)
            {
                for (int jj = 0; jj < dT.Count; jj++)
                {
                    //if (!(ii == 0 && jj == 0)) cdmtime[ii * dT.Count + jj] = (cdmtime[ii * dT.Count + jj - 1] * 1000 + dt_dT[jj] * dT.ydata[ii * dT.Count + jj]) / 1000 + cyctime[ii];
                    if (jj == 0) cdmtime.ydata[ii * dT.Count] = cyctime[ii] * time_factor;
                    else cdmtime.ydata[ii * dT.Count + jj] = (cdmtime.ydata[ii * dT.Count + jj - 1] + dt_dT[jj] * dT.ydata[ii * dT.Count + jj] / 1000);
                }
            }
            return cdmtime;
        }
        
        internal static double[] DT2CycRPM(CycleData dT)
        {
            var result = new double[dT.CycleCount];
            double meanvalue;
            for (int ii = 0; ii < dT.CycleCount; ii++)
            {
                meanvalue = 0;
                for (int jj = 0; jj < dT.Count; jj++) meanvalue += dT.ydata[ii * dT.Count + jj];
                result[ii] = 60 / (360 * (meanvalue / dT.Count) / 1000 / 1000);
            }
            return result;
        }

        private static CycleData AntiAliasingCAFilter(CycleData rawdata, double[] cycrpm, double fs_target)
        {
            double fc = fs_target / 2.56;
            int N = 42;
            //
            double fs_order = 360 / (rawdata.xdata[1] - rawdata.xdata[0]);
            //double maxorder;
            //double fc1;
            //var b = new double[N + 1];
            var result = new CycleData
            {
                Count = rawdata.Count,
                CycleCount = rawdata.CycleCount,
                xdata = rawdata.xdata,
                ydata = new double[rawdata.ydata.Length]
            };

            //var filtered_data = new double[rawdata.Count + N + 1];
            //var window = new double[N + 1];
            //for (int ii = 0; ii < window.Length; ii++) window[ii] = 0.54 - 0.46 * Math.Cos(2 * Math.PI * ii / (window.Length - 1));

            // for (int i = 0; i < cycrpm.Length; i++)



            Parallel.For(0, cycrpm.Length, i =>
            {
                var b = new double[N + 1];
                var bufferarray = new double[rawdata.Count + N + 1];
                var filtered_data = new double[bufferarray.Length];
                if (i == cycrpm.Length - 1) for (int jj = 0; jj < rawdata.Count; jj++) bufferarray[jj] = rawdata.ydata[i * rawdata.Count + jj];
                else for (int jj = 0; jj < bufferarray.Length; jj++) bufferarray[jj] = rawdata.ydata[i * rawdata.Count + jj];
                double maxorder = fc / cycrpm[i] * 60;
                double fc1 = maxorder / (fs_order / 2);
                if (maxorder < fs_order / 2)
                {

                    b = MathNet.Filtering.FIR.FirCoefficients.LowPass(fs_order, maxorder, N / 2);
                    b = WindowFcn(b, "hamming");
                    var FilterObj = new MathNet.Filtering.FIR.OnlineFirFilter(b);
                    filtered_data = FilterObj.ProcessSamples(bufferarray);
                    /*if (i == cycrpm.Length - 1) for (int jj = 0; jj < rawdata.Count - N / 2; jj++) result.ydata[i * rawdata.Count + jj] = filtered_data[jj + N / 2];
                    else*/
                    for (int jj = 0; jj < rawdata.Count; jj++) result.ydata[i * rawdata.Count + jj] = filtered_data[jj + N / 2];
                }
                else
                {
                    for (int jj = 0; jj < rawdata.Count; jj++) result.ydata[i * rawdata.Count + jj] = bufferarray[jj];
                }
            }
            );
            //;

            return result;
        }

        internal static double[] ResampleCAData_dT(CycleData rawdata, CycleData dT, double fs_new)
        {
            // takes CA Data and returns TM data

            // possibly the there is still some potential for performance increase

            //var rawdata_long = new CycleData();
            //var dT_long = new CycleData();

            var raw_diff = Diff(rawdata.xdata);
            var dT_diff = Diff(dT.xdata);

            double minres_raw = Math.Round(Min(raw_diff), 2);
            double maxres_raw = Math.Round(Max(raw_diff), 2);
            double minres_dT = Math.Round(Min(dT_diff), 2);
            double maxres_dT = Math.Round(Max(dT_diff), 2);

            if (minres_raw != maxres_raw) rawdata = InterpolateUnevenCycleData(rawdata, minres_raw);
            if (minres_dT != maxres_dT) dT = InterpolateUnevenCycleData(dT, minres_dT);

            var cycrpm = DT2CycRPM(dT);
            rawdata = AntiAliasingCAFilter(rawdata, cycrpm, fs_new);

            double m;
            int duplicates = (int)(minres_dT / minres_raw);

            double[] rawdata_x_ip;
            double[] rawdata_new_xdata;
            if (minres_dT / minres_raw == duplicates)
            {
                rawdata_x_ip = new double[rawdata.CycleCount * rawdata.Count];
                for (int jj = 0; jj < duplicates; jj++) rawdata_x_ip[jj] = dT.ydata[0] * minres_raw * jj;
                for (int ii = 1; ii < dT.ydata.Length; ii++)
                {
                    m = (dT.ydata[ii] - dT.ydata[ii - 1]) / duplicates;
                    for (int jj = 0; jj < duplicates; jj++) rawdata_x_ip[ii * duplicates + jj] = rawdata_x_ip[ii * duplicates - 1] + (dT.ydata[ii] + m * (jj + 1)) * minres_raw * (jj + 1);
                }
                rawdata_new_xdata = rawdata_x_ip;
            }
            else
            {
                var rawdata_long_xdata = VertcatCycleData(rawdata);
                var dT_long_xdata = VertcatCycleData(dT);
                double dt_dT = dT_long_xdata[1] - dT_long_xdata[0];
                for (int ii = 0; ii < dT.ydata.Length; ii++) dT.ydata[ii] *= dt_dT;
                var microseconds = Cumsum(dT.ydata);
                rawdata_x_ip = Interp1_linear(dT_long_xdata, microseconds, rawdata_long_xdata);
                rawdata_new_xdata = RemoveNaN(rawdata_x_ip, rawdata_x_ip);
                rawdata.ydata = RemoveNaN(rawdata.ydata, rawdata_x_ip);
            }

            var new_timeaxis = DeltaArray(0, 1000000 / fs_new, rawdata_new_xdata[rawdata_new_xdata.Length - 1]);

            /*System.Diagnostics.Stopwatch tic = new System.Diagnostics.Stopwatch();
            tic.Start();*/

            var result = Interp1_spline(rawdata_new_xdata, rawdata.ydata, new_timeaxis);
            result = RemoveNaN(result, result);
            /*long calctime = tic.ElapsedMilliseconds;
            Console.WriteLine("Calculation Time: " + calctime);
            Console.ReadKey();*/

            return result;
        }

        internal static double[] ResampleCAData_CDMTIME(CycleData rawdata, CycleData cdmtime, double fs_new)
        {
            // takes CA Data and returns TM data

            // possibly the there is still some potential for performance increase

            //var rawdata_long = new CycleData();
            //var dT_long = new CycleData();
          
            var raw_diff = Diff(rawdata.xdata);
            var cdmtime_diff = Diff(cdmtime.xdata);

            double minres_raw = Math.Round(Min(raw_diff), 2);
            double maxres_raw = Math.Round(Max(raw_diff), 2);
            double minres_cdmtime = Math.Round(Min(cdmtime_diff), 2);
            double maxres_cdmtime = Math.Round(Max(cdmtime_diff), 2);

            if (minres_raw != maxres_raw) rawdata = InterpolateUnevenCycleData(rawdata, minres_raw);
            if (minres_cdmtime != maxres_cdmtime) cdmtime = InterpolateUnevenCycleData(cdmtime, minres_cdmtime);

            var cycrpm = CDMTIME2CycRPM(cdmtime);
            rawdata = AntiAliasingCAFilter(rawdata, cycrpm, fs_new);

            double m;
            int duplicates = (int)(minres_cdmtime / minres_raw);

            double[] rawdata_x_ip;
            double[] rawdata_new_xdata;
            if (minres_cdmtime / minres_raw == duplicates)
            {
                rawdata_x_ip = new double[rawdata.CycleCount * rawdata.Count];
                rawdata_x_ip[0] = cdmtime.ydata[0];
                //for (int jj = 0; jj < duplicates; jj++) rawdata_x_ip[jj] = cdmtime.ydata[0] * (jj + 1) / duplicates;
                for (int ii = 0; ii < cdmtime.ydata.Length; ii++)
                {
                    if (ii == cdmtime.ydata.Length - 1)
                    {
                        m = (cdmtime.ydata[ii] - cdmtime.ydata[ii - 1]) / duplicates;
                        for (int jj = 0; jj < duplicates - 1; jj++) rawdata_x_ip[ii * duplicates + (jj + 1)] = rawdata_x_ip[ii * duplicates] + m * (jj + 1);
                    }
                    else
                    {
                        m = (cdmtime.ydata[ii + 1] - cdmtime.ydata[ii]) / duplicates;
                        for (int jj = 0; jj < duplicates; jj++) rawdata_x_ip[ii * duplicates + (jj + 1)] = rawdata_x_ip[ii * duplicates] + m * (jj + 1);
                    }
                }
                rawdata_new_xdata = rawdata_x_ip;
            }
            else
            {
                var rawdata_long_xdata = VertcatCycleData(rawdata);
                var cdmtime_long_xdata = VertcatCycleData(cdmtime);
                //double dt_cdmtime = cdmtime_long_xdata[1] - cdmtime_long_xdata[0];
                //for (int ii = 0; ii < cdmtime.ydata.Length; ii++) cdmtime.ydata[ii] *= dt_cdmtime;
                //var microseconds = Cumsum(cdmtime.ydata);
                //rawdata_x_ip = Interp1_linear(cdmtime_long_xdata, microseconds, rawdata_long_xdata);
                rawdata_x_ip = Interp1_linear(cdmtime_long_xdata, cdmtime.ydata, rawdata_long_xdata);
                rawdata_new_xdata = RemoveNaN(rawdata_x_ip, rawdata_x_ip);
                rawdata.ydata = RemoveNaN(rawdata.ydata, rawdata_x_ip);
            }

            var new_timeaxis = DeltaArray(rawdata_new_xdata[0], 1000 / fs_new, rawdata_new_xdata[rawdata_new_xdata.Length - 1]);

            /*System.Diagnostics.Stopwatch tic = new System.Diagnostics.Stopwatch();
            tic.Start();*/

            var result = Interp1_spline(rawdata_new_xdata, rawdata.ydata, new_timeaxis);
            result = RemoveNaN(result, result);
            /*long calctime = tic.ElapsedMilliseconds;
            Console.WriteLine("Calculation Time: " + calctime);
            Console.ReadKey();*/

            return result;
        }

        /*private float[] Double2Float(double[] input)
        {
            var output = new float[input.Length];
            for (int ii = 0; ii < input.Length; ii++) output[ii] = (float)input[ii];
            return output;
        }*/

        private static int[] Cumsum(int[] input)
        {
            // Function similar to Matlabs Cumsum with one difference: result = [0;Cumsum_Matlab(input)];
            //                                                         result = result(1:end-1);
            var output = new int[input.Length];
            //output[0] = input[0];
            output[0] = 0;
            for (int ii = 1; ii < input.Length; ii++) output[ii] = input[ii - 1] + output[ii - 1];
            return output;
        }

        private static double[] Cumsum(double[] input)
        {
            // Function similar to Matlabs Cumsum with one difference: result = [0;Cumsum_Matlab(input)];
            //                                                         result = result(1:end-1);
            var output = new double[input.Length];
            //output[0] = input[0];
            output[0] = 0;
            for (int ii = 1; ii < input.Length; ii++) output[ii] = input[ii - 1] + output[ii - 1];
            return output;
        }

        public static double[] DeltaArray(double start, double delta, double end)
        {
            if (end < start) return null;
            double[] output = new double[(int)((end - start) / delta) + 1];
            for (int ii = 0; ii < output.Length; ii++) output[ii] = start + ii * delta;
            return output;
        }

        private static double[] Diff(double[] input)
        {
            double[] output = new double[input.Length - 1];
            for (int ii = 0; ii < output.Length; ii++) output[ii] = input[ii + 1] - input[ii];
            return output;
        }
        
        public static double Min(double[] input)
        {
            double output = input[0];
            for (int ii = 1; ii < input.Length; ii++) if (input[ii] < output) output = input[ii];
            return output;
        }

        public static double Max(double[] input)
        {
            double output = input[0];
            for (int ii = 1; ii < input.Length; ii++) if (input[ii] > output) output = input[ii];
            return output;
        }

        private static double OrderLevelCorrection(double[] z, int Index_Ord)
        {
            double z1 = z[Index_Ord - 1];
            double z2 = z[Index_Ord];
            double z3 = z[Index_Ord + 1];

            z1 = 20 * Math.Log10(z1 / 0.00002);
            z2 = 20 * Math.Log10(z2 / 0.00002);
            z3 = 20 * Math.Log10(z3 / 0.00002);

            double c = (z3 - 2 * z2 + z1) / 2;
            double x = (z1 - z2 - c) / 2 / c;
            if (x < -1) x = -1;
            if (x > 1) x = 1;
            double r = z2 + (z2 - z1 + c) * x + c * x * x;
            if (r < z2) r = z2;
            double Pegel_Ord = 0.00002 * Math.Pow(10, (r / 20));
            return Pegel_Ord;
        }

       private static double DOC_OrderAnalysis(double[] y, double[] z, double speed, int powersport)
       {
            // powersport Parameter: 0: Power, 1: Sport
            if (speed <= 0) return 0;
            var df = y[1] - y[0];
            double f_max = y[y.Length - 1];
            double Ord_Aufloesung = 0.5;
            double erste_Ord = Ord_Aufloesung;
            double hoechste_Ord = Math.Floor(600000 / speed); //%{ fg = 10000 HZ; hoechste Ord. := fg / (n / 60) }
            int Abstand = 1;
            double Rausch_Offs = 1.4125; //%{ Offset der zum Rauschen addiert wird(1.4125 entspricht 3dB) }
            int BB = 10;  //  %        { Bandbreite(%) in der die Ordnungsspitze gesucht wird}

            double Samples_Ordnungen = hoechste_Ord / Ord_Aufloesung * 2 + 1; //% { Länge des Ordnungs -}
            var Ord_Spectra = new double[(int)Math.Round(Samples_Ordnungen,MidpointRounding.AwayFromZero)]; //%   { spektrums }
            var Index_Data = new double[Ord_Spectra.Length];

            // Berechnung der Ordnungspegel
            // 10 % Suchbereich (hohe Frequenzen
            if (speed > df / 0.002 && speed < f_max * 6 / 7)
            {
                double f_Ord_unten = speed * (100 - BB / 2) / 6000;     // Frequenz der 1. Ordnung - 5 %
                double f_Ord_oben = speed * (100 + BB / 2) / 6000;      // Frequenz der 1. Ordnung + 5 %
                double delta_f = (f_Ord_oben - f_Ord_unten) / 2;        // Halbe Bandbreite in Hz
                int o = 2;
                double Ordnung = erste_Ord;

                while (Ordnung <= hoechste_Ord)
                {
                    f_Ord_unten = Ordnung * speed / 60 - delta_f;  
                    f_Ord_oben = Ordnung * speed / 60 + delta_f;
                    int Index_f_Ord_unten = (int)Math.Ceiling(f_Ord_unten / df);
                    int Index_f_Ord_oben = (int)Math.Floor(f_Ord_oben / df);
                    Ordnung += Ord_Aufloesung;

                    // Maximum im Toleranzbereich finden... daraus den Index der Ordnung Index_Ord bestimmen
                    int i_temp = Index_f_Ord_unten;
                    int Index_Ord = i_temp;
                    double Pegel = z[i_temp];
                    while (i_temp < Index_f_Ord_oben)
                    {
                        if (Pegel < z[i_temp + 1])
                        {
                            Pegel = z[i_temp + 1];
                            Index_Ord = i_temp + 1;
                        }
                        i_temp++;
                    }
                    // Pegel korrigieren... ??
                    var Pegel_Ord = OrderLevelCorrection(z, Index_Ord);
                    Ord_Spectra[o] = Pegel_Ord;
                    Index_Data[o] = Index_Ord;
                    o += 2;
                }
                
            }
            // Ohne Suchbereich (tiefe Frequenzen)
            else if ((240 < speed / df) && (speed <= df / 0.002))
            {
                int o = 2;
                double Ordnung = erste_Ord;
                while (Ordnung <= hoechste_Ord)
                {
                    double f_Ord = Ordnung * speed / 60;
                    Ordnung += Ord_Aufloesung;
                    int Index_Ord = (int)Math.Round(f_Ord / df, MidpointRounding.AwayFromZero);
                    var Pegel_Ord = OrderLevelCorrection(z, Index_Ord);
                    Ord_Spectra[o] = Pegel_Ord;
                    Index_Data[o] = Index_Ord;
                    o += 2;
                }
            }
            
            // Berechnung Rauschpegel

            int Index_fgu = (int)(9 / df);
            double Rauschen = erste_Ord - Ord_Aufloesung / 2; // { ergibt 0,25 }
            while (Rauschen <= hoechste_Ord)
            {
                int Index_f_Ord_unten = (int)Math.Round(Rauschen / Ord_Aufloesung * 2 - 1, MidpointRounding.AwayFromZero);
                int Index_Spektrum_u = (int)Math.Round(Index_Data[Index_f_Ord_unten], MidpointRounding.AwayFromZero);
                int I_Rauschen_u = Index_Spektrum_u;
                
                int Index_f_Ord_oben = Index_f_Ord_unten + 2;
                int Index_Spektrum_o = (int)Math.Round(Index_Data[Index_f_Ord_oben], MidpointRounding.AwayFromZero);
                int I_Rauschen_o = Index_Spektrum_o;

                int Anz_Linien = Index_Spektrum_o - Index_Spektrum_u - 1;

            //%{ Anzahl der Linien die sich zwischen 2 Ordnungen befinden }

                int Index_Rauschen = Index_f_Ord_unten + 1;
                double Rauschwert = 0;
                int Teiler = 0;

                if (Anz_Linien == 0)
                {
                    return double.NaN;
                }
                else if (Anz_Linien == 1) Ord_Spectra[Index_Rauschen] = z[Index_Spektrum_u + 1] * Rausch_Offs; //%{ Offsetkorrektur + 3 dB }
                else if (Anz_Linien == 2)
                {
                    Rauschwert = Math.Pow(z[Index_Spektrum_u + 1], 2);
                    Rauschwert += Math.Pow(z[Index_Spektrum_u + 2], 2);
                    Rauschwert = Math.Sqrt(Rauschwert / 2) * Rausch_Offs; //%{ Offsetkorrektur + 3 dB }
                    Ord_Spectra[Index_Rauschen] = Rauschwert;
                }
                else if (Anz_Linien >= 3)
                {
                    I_Rauschen_u += Abstand;
                    I_Rauschen_o -= Abstand;
                    // %{ Suchen des Tiefpunktes von der unteren Ord.ausgehend }
                    while ((I_Rauschen_o < z.Length - 1) && (z[I_Rauschen_u] > z[I_Rauschen_u + 1])) I_Rauschen_u ++;
                    // %{ Suchen des Tiefpunktes von der oberen Ord.ausgehend }
                    while ((I_Rauschen_o > 0) && (z[I_Rauschen_o] > z[I_Rauschen_o - 1])) I_Rauschen_o --;
                    // %   { wenn der Suchalgorithmus innerhalb der 2 Ordnungen nicht zum stehen kommt }
                    if (I_Rauschen_u >= Index_Spektrum_o) I_Rauschen_u = Index_Spektrum_o - 1;
                    if (I_Rauschen_o <= Index_Spektrum_u) I_Rauschen_o = Index_Spektrum_u + 1;
                    // %{ Beim ersten Rauschpegel werden die untersten 2 Freq.- Linien weggelassen }
                    if ((I_Rauschen_u == 1) && (I_Rauschen_o == 1))
                    {
                        I_Rauschen_u = Index_fgu;
                        I_Rauschen_o = Index_fgu;
                    }
                    // %{ Quadratischer Mittelwert der Rauschpegeln }
                    int j = I_Rauschen_u;
                    while (j <= I_Rauschen_o)
                    {
                        Rauschwert += Math.Pow(z[j], 2);
                        Teiler++;
                        j++;
                    }
                    Rauschwert = Math.Sqrt(Rauschwert / Teiler) * Rausch_Offs; //%{ Offsetkorrektur + 3 dB};

                    Ord_Spectra[Index_Rauschen] = Rauschwert;
                    j = I_Rauschen_u;
                    while (j <= I_Rauschen_o) j++;
                }
        
                Rauschen += Ord_Aufloesung;
            }
              //      Glättung der Rauschpegel 
            
            int Laenge_Ord_Sp = Ord_Spectra.Length; //%{ = 280 bei hoechste Ord. = 70}
            //%{ Der erste Rauschpegel wird mit dem nächsten arithmetisch gemittelt }
            double Mittelwert = (Math.Log10(Ord_Spectra[1] * Ord_Spectra[3])) / 2;
            Ord_Spectra[1] = Math.Pow(10, Mittelwert);
            //   Bildung des arithm.Mittelw.des Rauschpegels m.seinen benachbarten Rauschpegeln }
            int i = 3;
            while (i <= Laenge_Ord_Sp - 3)
            {
                Mittelwert = (Math.Log10(Ord_Spectra[i] * Ord_Spectra[i + 2] * Ord_Spectra[i - 2])) / 3;
                Ord_Spectra[i] = Math.Pow(10, Mittelwert);
                i += 2;
            }
            //     % { Der letzte Rauschpegel wird mit dem vorletzen gemittelt }
            Mittelwert = (Math.Log10(Ord_Spectra[Laenge_Ord_Sp - 2] * Ord_Spectra[Laenge_Ord_Sp - 4])) / 2;
            Ord_Spectra[Laenge_Ord_Sp - 2] = Math.Pow(10, Mittelwert);

            // Grenzfrequenzen und Bewertungen für den Bereich DOCPower
            var DOCweighting = new XY_Data
            {
                xdata = new double[4] { 60, 120, 500, 750},
                ydata = new double[4] { 0, 1, 1, 0}
            };
            if (powersport == 1)
            {
                DOCweighting.ydata[2] = 2;
                double k_u = 0.0625;
                double d_u = 300;
                double k_o = 0.175;
                double d_o = 600;
                double Anstieg_u = 80;
                double Anstieg_o = 80;

                DOCweighting.xdata[1] = speed * k_u + d_u;
                DOCweighting.xdata[2] = speed * k_o + d_o;
                DOCweighting.xdata[0] = DOCweighting.xdata[1] * (100 - Anstieg_u) / 100;
                DOCweighting.xdata[3] = DOCweighting.xdata[2] * (Anstieg_o / 100 + 1);
            }

            var Ord_Sp_W = new double[Ord_Spectra.Length];
            double delta_o = 0.25;
            // Umrechnung von Frequenz auf Index der Ordnung
            var Index_Ord_array = new int[DOCweighting.xdata.Length];
            for (int ii = 0; ii < Index_Ord_array.Length; ii++) Index_Ord_array[ii] = (int)Math.Round((DOCweighting.xdata[ii] * 60 / speed) / delta_o, MidpointRounding.AwayFromZero);
            // Berechnung von Steigung und Offset
            var k = new double[Index_Ord_array.Length - 1];
            var d = new double[k.Length];
            for (int ii = 0; ii < k.Length; ii++)
            {
                k[ii] = (DOCweighting.ydata[ii + 1] - DOCweighting.ydata[ii]) / (Index_Ord_array[ii + 1] - Index_Ord_array[ii]);
                d[ii] = DOCweighting.ydata[ii] - k[ii] * Index_Ord_array[ii];
            }
            // Bewertung des Ordnungsspektrums
            for (int ii = Index_Ord_array[0]; ii < Ord_Spectra.Length; ii++)
            {
                if (ii >= Index_Ord_array[0] && ii < Index_Ord_array[1]) Ord_Sp_W[ii] = Ord_Spectra[ii] * (k[0] * ii + d[0]);
                else if (ii >= Index_Ord_array[1] && ii < Index_Ord_array[2]) Ord_Sp_W[ii] = Ord_Spectra[ii] * (k[1] * ii + d[1]);
                else if (ii >= Index_Ord_array[2] && ii < Index_Ord_array[3]) Ord_Sp_W[ii] = Ord_Spectra[ii] * (k[2] * ii + d[2]);
            }

            // Auswertung des bewerteten Ordnungsspektrums
            double summe = 0;
            i = 2;
            double Rausch_Pegel;
            double Ord_Pegel;
            while (i < Ord_Sp_W.Length - 1)
            {
                if (Ord_Sp_W[i] != 0)
                {
                    if (Ord_Sp_W[i - 1] > Ord_Sp_W[i + 1]) Rausch_Pegel = Math.Pow(Ord_Sp_W[i - 1], 2);
                    else Rausch_Pegel = Math.Pow(Ord_Sp_W[i + 1], 2);
                    Ord_Pegel = Math.Pow(Ord_Sp_W[i], 2);
                    Ord_Pegel -= Rausch_Pegel;
                    if (Ord_Pegel < 0) Ord_Pegel = 0;
                    summe += Ord_Pegel;
                }
                i += 2;
            }
            
            return Math.Sqrt(summe);
        }


        private static double[] ABC_Weight(double[] f, int weight_curve)
        {
            //double df;
            double[] result;
            int ii;

            /*df = fs / N;
            f = new double[N / 2 + 1];*/
            result = new double[f.Length];

            //for (ii = 0; ii < f.Length; ii++) f[ii] = ii * df;
            //f[0] = 0.0001;
            //if (weight_curve == 1) for (ii = 0; ii < result.Length; ii++) result[ii] = Math.Pow(10, 0.1) * Math.Pow(12200, 2) * Math.Pow(f[ii], 4) / ((Math.Pow(f[ii], 2) + Math.Pow(20.6, 2)) * (Math.Pow(f[ii], 2) + Math.Pow(12200, 2)) * Math.Sqrt(Math.Pow(f[ii], 2) + Math.Pow(107.7, 2)) * Math.Sqrt(Math.Pow(f[ii], 2) + Math.Pow(737.9, 2)));
            if (weight_curve == 1) for (ii = 0; ii < result.Length; ii++) result[ii] = 187374174.3 * Math.Pow(f[ii], 4) / ((Math.Pow(f[ii], 2) + Math.Pow(20.6, 2)) * (Math.Pow(f[ii], 2) + Math.Pow(12200, 2)) * Math.Sqrt(Math.Pow(f[ii], 2) + Math.Pow(107.7, 2)) * Math.Sqrt(Math.Pow(f[ii], 2) + Math.Pow(737.9, 2)));
            
            else if (weight_curve == 2) for (ii = 0; ii < result.Length; ii++) result[ii] = Math.Pow(10, (0.0085)) * Math.Pow(12200, 2) * Math.Pow(f[ii], 3) / ((Math.Pow(f[ii], 2) + Math.Pow(20.6, 2)) * (Math.Pow(f[ii], 2) + Math.Pow(12200, 2)) * Math.Sqrt(Math.Pow(f[ii], 2) + Math.Pow(158.5, 2)));
            else if (weight_curve == 3) for (ii = 0; ii < result.Length; ii++) result[ii] = Math.Pow(10, (0.003)) * Math.Pow(12200, 2) * Math.Pow(f[ii], 2) / ((Math.Pow(f[ii], 2) + Math.Pow(20.6, 2)) * (Math.Pow(f[ii], 2) + Math.Pow(12200, 2)));
            return result;
        }

        private static F_octave CalcCornerFreqs(double n, double fs, double N)
        {
            F_octave result;
            /*double G = Math.Pow(10, (double)3 / 10);
            double fr = 1000;
            var f = DeltaArray(0, fs / N, fs / 2.56);
            int x = 30;
            int bw = 1;
            double fm_temp, f1_temp, f2_temp;
            int f1i, f2i;
            // check direction to 0 Hz
            while (bw > 0)
            {
                x--;
                if (n % 2 != 0) fm_temp = Math.Pow(G, (x - 30) / n) * fr;
                else fm_temp = Math.Pow(G, (2 * x - 59) / (2 * n)) * fr;
                f1_temp = Math.Pow(G, (-1 / (2 * n))) * fm_temp;
                f2_temp = Math.Pow(G, (1 / (2 * n))) * fm_temp;
                
                f2i = f.Length - 1;
                for (int ii = 0; ii < f.Length; ii++)
                {
                    if (f[ii] < f2_temp) f2i = ii;
                    else break;
                }
                f1i = 0;
                for (int ii = f.Length - 1; ii >= 0; ii--)
                {
                    if (f[ii] >= f1_temp) f1i = ii;
                    else break;
                }
                bw = f2i - f1i + 1;
            }
            int x1 = x + 1;

            // check direction to fs/2
            x = 30;
            bw = 1;
            f2_temp = 0;
            while (f2_temp <= f[f.Length - 1])
            {
                x++;
                if (n % 2 != 0) fm_temp = Math.Pow(G, (x - 30) / n) * fr;
                else fm_temp = Math.Pow(G, (2 * x - 59) / (2 * n)) * fr;
                f2_temp = Math.Pow(G, (1 / (2 * n))) * fm_temp;
            }
            int x2 = x - 1;

            var f1 = new double[x2 - x1 + 1];
            var f2 = new double[x2 - x1 + 1];
            var fm = new double[x2 - x1 + 1];

            for (x = x1; x <= x2; x++)
            {
                if (n % 2 != 0) fm[x - x1] = Math.Pow(G, (x - 30) / n) * fr;
                else fm[x - x1] = Math.Pow(G, (2 * x - 59) / (2 * n)) * fr;
                f1[x - x1] = Math.Pow(G, (-1 / (2 * n))) * fm[x - x1];
                f2[x - x1] = Math.Pow(G, (1 / (2 * n))) * fm[x - x1];
            }

            result = new F_octave
            {
                f0 = fm,
                f1 = f1,
                f2 = f2
            };*/
            double fstart;
            int ii, L;

            ii = 1;
            result = new F_octave();

            if (n == 1) fstart = 1000 / 128.0 * 16384 / N;     // 128 = 2^(14/2)
            else if (n == 3) fstart = 1000 / Math.Pow(2, (34 / 6.0)) * 16384 / N;
            else fstart = (n / 6) * 1000 / Math.Pow(2, (58 / 12.0)) * 16384 / N;
            while (fstart * Math.Pow(2, ((ii - 1) / n)) * (Math.Pow(2, (1 / (2 * n)))) < (fs / 2.56)) ii++;

            L = ii - 1;
            result.f0 = new double[L];
            result.f1 = new double[L];
            result.f2 = new double[L];

            for (ii = 0; ii < L; ii++)
            {
                result.f0[ii] = fstart * Math.Pow(2, (ii / n));
                
            }

            var normfreqs_terz = new double[29] { 20, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000 };
            var normfreqs_okt = new double[10] { 31.25, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000 };

            double[] normfreqs;
            if (n == 3)
            {
                normfreqs = normfreqs_terz.Where(f => f > result.f0[0] && f < result.f0[result.f0.Length - 1]).ToArray();
            }
            else
            {
                normfreqs = normfreqs_okt.Where(f => f > result.f0[0] && f < result.f0[result.f0.Length - 1]).ToArray();
            }
            for (ii = 0; ii < normfreqs.Length; ii++)
            {
                double minabs = Double.PositiveInfinity;
                int index = 0;
                for (int kk = 0; kk < result.f0.Length; kk++)
                {
                    if (Math.Abs(result.f0[kk] - normfreqs[ii]) < minabs)
                    {
                        minabs = Math.Abs(result.f0[kk] - normfreqs[ii]);
                        index = kk;
                    }
                }
                result.f0[index] = normfreqs[ii];
            }
            for (ii = 0; ii < L; ii++)
            {
                result.f1[ii] = result.f0[ii] / Math.Pow(2, (1 / (2 * n)));
                result.f2[ii] = result.f0[ii] * Math.Pow(2, (1 / (2 * n)));
            }

            return result;
        }

        internal static double[] ParseOrdersFromString(string orderstring)
        {
            //// Variable Definitions

            int L,                  // Number of orders
                NO,                 // Number of orders in one instruction
                ii;                 // Running variable

            double[] orderarray,    // Result array with every order
                orderarray_new;     // Result array with every order (singularity!)

            // Define number format
            NumberFormatInfo provider = new NumberFormatInfo
            {
                NumberDecimalSeparator = "."
            };
            

            // Split the input string into instructions (e.g. "2:4,6,8:0.5:9" results in "2:4" and "6" and "8:0.5:9")
            string[] splitstring = orderstring.Split(',', ';');

            // Find out how many orders are needed
            double[] length = new double[splitstring.Length];
            L = 0;
            for (ii = 0; ii < splitstring.Length; ii++)
            {
                string[] singlestringarray = splitstring[ii].Split(':');
                length[ii] = singlestringarray.Length;
                if (length[ii] == 1) L++;
                else if (length[ii] == 2) L += Convert.ToInt16(singlestringarray[1]) - Convert.ToInt16(singlestringarray[0]) + 1;
                else if (length[ii] == 3)
                {
                    L += (Int16)((Convert.ToDouble(singlestringarray[2], provider) - Convert.ToDouble(singlestringarray[0], provider)) / Convert.ToDouble(singlestringarray[1], provider) + 1);
                }
            }

            // Determine every order value
            orderarray = new double[L];
            L = 0;
            for (ii = 0; ii < splitstring.Length; ii++)
            {
                string[] singlestringarray = splitstring[ii].Split(':');
                length[ii] = singlestringarray.Length;
                if (length[ii] == 1)
                {
                    orderarray[L] = Convert.ToDouble(singlestringarray[0], provider);
                    L++;
                }
                else if (length[ii] == 2)
                {
                    NO = Convert.ToInt16(singlestringarray[1]) - Convert.ToInt16(singlestringarray[0]) + 1;
                    for (int jj = 0; jj < NO; jj++) orderarray[L + jj] = Convert.ToDouble(singlestringarray[0], provider) + (double)jj;
                    L += NO;
                }
                else if (length[ii] == 3)
                {
                    NO = (Int16)((Convert.ToDouble(singlestringarray[2], provider) - Convert.ToDouble(singlestringarray[0], provider)) / Convert.ToDouble(singlestringarray[1], provider) + 1);
                    for (int jj = 0; jj < NO; jj++) orderarray[L + jj] = Convert.ToDouble(singlestringarray[0], provider) + (double)jj * Convert.ToDouble(singlestringarray[1], provider);
                    L += NO;
                }
            }

            // Sort and check for singularity
            Array.Sort(orderarray);
            L = 1;
            for (ii = 1; ii < orderarray.Length; ii++) if (orderarray[ii] > orderarray[ii - 1]) L++;
            if (L < orderarray.Length)
            {
                orderarray_new = new double[L];
                L = 1;
                orderarray_new[0] = orderarray[0];
                for (ii = 1; ii < orderarray.Length; ii++) if (orderarray[ii] > orderarray[ii - 1])
                    {
                        orderarray_new[L] = orderarray[ii];
                        L++;
                    }
                orderarray = orderarray_new;
            }
            return orderarray;
        }

        internal static XYZ_Data AverageDiagramOverTime(XYZ_Data result, double[] aps_plaindata_xdata)
        {
            // Averaging over time if 2D Diagram was chosen

            int rows_out = result.zdata.Length / aps_plaindata_xdata.Length;
            double[] average_zdata = new double[rows_out];
            double[] average_ydata = new double[rows_out];
            double[] average_xdata = new double[1];
            average_xdata[0] = (aps_plaindata_xdata[aps_plaindata_xdata.Length - 1] - aps_plaindata_xdata[0]) / 2 + aps_plaindata_xdata[0];
            Parallel.For(0, rows_out, j =>
            {
                int i;
                double onefreq_alltime = 0;
                for (i = 0; i < aps_plaindata_xdata.Length; i++) onefreq_alltime += result.zdata[i * rows_out + j];
                average_zdata[j] = onefreq_alltime / (double)aps_plaindata_xdata.Length;
                average_ydata[j] = result.ydata[j];
            }
            );
            result.zdata = average_zdata;
            result.ydata = average_ydata;
            result.xdata = average_xdata;
            result.dimensions[0] = 1;
            return result;
        }

        public static XY_Data Interp1(double[] x, double[] y, double[] xq)
        {
            // Linear interpolation function. Does not extrapolate. Query points beyond the x-y-grid are omitted => yout could be shorter than x and y.
            // returns XY structure and works for both increasing and decreasing x values

            double direction = x[x.Length - 1] - x[0];

            double[] yq = new double[xq.Length];
            double[] xq1 = new double[xq.Length];
            var spline = Interpolate.Linear(x, y);
            int jj = 0;
            if (direction >= 0)
            {
                for (int ii = 0; ii < xq.Length; ii++)
                {
                    if (xq[ii] >= x[0] && xq[ii] <= x[x.Length - 1])
                    {
                        yq[jj] = spline.Interpolate(xq[ii]);
                        xq1[jj] = xq[ii];
                        jj++;
                    }
                }
            }
            else
            {
                for (int ii = 0; ii < xq.Length; ii++)
                {
                    if (xq[ii] <= x[0] && xq[ii] >= x[x.Length - 1])
                    {
                        yq[jj] = spline.Interpolate(xq[ii]);
                        xq1[jj] = xq[ii];
                        jj++;
                    }
                }
            }
                
            XY_Data output = new XY_Data
            {
                xdata = new double[jj],
                ydata = new double[jj]
            };
            for (int ii = 0; ii < jj; ii++)
            {
                output.xdata[ii] = xq1[ii];
                output.ydata[ii] = yq[ii];
            }
            return output;
        }

        private static double[] Interp1_spline(double[] x, double[] y, double[] xq)
        {
            // multi core spline interpolation

            double direction = x[x.Length - 1] - x[0];

            double[] yq = new double[xq.Length];
            double[] xq1 = new double[xq.Length];

            var spline = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkimaSorted(x,y);
            //var spline = MathNet.Numerics.Interpolation.LinearSpline.InterpolateSorted(x, y);

            //var spline = Interpolate.CubicSplineRobust(x, y);

            int L_singlepackage = xq.Length / 16;
            //int L_singlepackage = xq.Length;
            //int i = 0;
            /*System.Diagnostics.Stopwatch tic = new System.Diagnostics.Stopwatch();
            tic.Start();*/
            //for (int i = 0; i < 16; i++)
            Parallel.For(0, 16, i =>
            {
                for (int ii = 0; ii < L_singlepackage; ii++)
                {
                    if (xq[i * L_singlepackage + ii] >= x[0] && xq[i * L_singlepackage + ii] <= x[x.Length - 1])
                    {
                        yq[i * L_singlepackage + ii] = spline.Interpolate(xq[i * L_singlepackage + ii]);

                    }
                    else yq[i * L_singlepackage + ii] = Double.NaN;
                }
            }
            );
            for (int ii = L_singlepackage * 16; ii < xq.Length; ii++)
            {
                if (xq[ii] >= x[0] && xq[ii] <= x[x.Length - 1])
                {
                    yq[ii] = spline.Interpolate(xq[ii]);

                }
                else yq[ii] = Double.NaN;
            }
            /*long calctime = tic.ElapsedMilliseconds;
            Console.WriteLine("Calculation Time: " + calctime);
            Console.ReadKey();*/
            /*
                for (int ii = 0; ii < xq.Length; ii++)
                {
                    if (xq[ii] >= x[0] && xq[ii] <= x[x.Length - 1])
                    {
                        yq[ii] = spline.Interpolate(xq[ii]);
                        
                    }
                    else yq[ii] = Double.NaN;
                }
            */

            return yq;
        }


        public static double[] Moving(double[] x, int filterlength)
        {
            if (x == null || filterlength == 0) return null;
            // Moving Average Filter. Works only for odd filterlength. When even filterlength is provided, input is returned unchanged.
            if (filterlength % 2 != 1) return x;
            var result = new double[x.Length];
            double sum;

            for (int ii = filterlength / 2; ii < x.Length - filterlength / 2; ii++)
            {
                sum = 0;
                for (int jj = 0; jj < filterlength; jj++)
                {
                    sum += x[ii - filterlength / 2 + jj];
                }
                result[ii] = sum / filterlength;
            }
            for (int ii = 0; ii < filterlength / 2; ii++)
            {
                result[ii] = result[filterlength / 2];
                result[result.Length - ii - 1] = result[result.Length - filterlength / 2];
            }
            return result;

        }

        public static double[] REO(double[][] audio, double fs, BackgroundWorker worker = null, double startprogress = 0, double overallprogress = 0)
        {
            if (audio == null) return null;
            int N_REO = 8192;
            double overlap = 0.5;
            double percent = 0;
            worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
            double audio_length_seconds = audio[0].Length / fs;

            var params_REO = new Param_struct
            {
                delta_f = fs / N_REO,
                delta_t = N_REO / fs * (1 - overlap),
                average = 1,
                average_overlap = 0,
                fs = fs,
                windowtype = 1,
                DC = 1,
                y_amplitude = 0,
                y_axis = 0,
                dBref = 0.00002,
                freq_weight = 0,
                dsrange1 = 0,
                dsrange2 = audio_length_seconds,
                LQ = 0,
                y_unitconversion = 1
            };

            int channels = audio.Length;

            var REO2ch = new double[channels][];
            var REO = new XY_Data();

            var aps_y = Array.Empty<double>();

            for (int ch = 0; ch < channels; ch++)
            {
                var aps = CalcAPSData(audio[ch], new XY_Data(), 1 / fs, 0, new CATimeStamps(), params_REO);
                percent = (ch + 1) / channels * 50;
                worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
                if (ch == 0)
                {
                    aps_y = new double[aps.dimensions[1]];
                    Array.Copy(aps.ydata, aps_y, aps.dimensions[1]);
                    REO.xdata = new double[aps.dimensions[0]];
                    for (int ii = 0; ii < aps.dimensions[0]; ii++) REO.xdata[ii] = aps.xdata[ii * aps.dimensions[1]];
                }
                var aps_z = new double[aps.dimensions[1]];
                REO2ch[ch] = new double[aps.dimensions[0]];

                for (int ii = 0; ii < aps.dimensions[0]; ii++)
                {
                    Array.Copy(aps.zdata, ii * aps.dimensions[1], aps_z, 0, aps.dimensions[1]);
                    int Linien_Win = (int)Math.Floor(1000 / aps_y[1]) + 1;
                    var spec_w = new double[Linien_Win];
                    Array.Copy(aps_z, spec_w, Linien_Win);
                    for (int jj = 0; jj < spec_w.Length; jj++) spec_w[jj] = 20 * Math.Log10(spec_w[jj] / 0.00002);
                    spec_w = WindowFcn(spec_w, "hann");
                    double singleresult = 0;
                    for (int jj = 0; jj < spec_w.Length; jj++) singleresult += spec_w[jj] * aps_y[1]/18750;
                    if (singleresult > 0) singleresult = Math.Log(singleresult);
                    REO2ch[ch][ii] = singleresult;
                    percent = 50 + (double)ii / aps.dimensions[0] * 45;
                    worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
                }
            }
            var REO_mean = new double[REO2ch[0].Length];
            for (int ii = 0; ii < REO_mean.Length; ii++)
            {
                double summe = 0;
                for (int ch = 0; ch < channels; ch++) summe += REO2ch[ch][ii];
                REO_mean[ii] = summe / channels;
            }
            REO.ydata = Interp1_linear(REO.xdata, REO_mean, DeltaArray(0.1, 0.025, audio_length_seconds - 0.1));
            REO.xdata = DeltaArray(0.1, 0.025, audio_length_seconds - 0.1);
            ExtrapolateNaNs(REO);
            percent = 100;
            worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
            return REO.ydata;
        }

        private static XY_Data GetTerzSpectrumFromAPS(double[] z, double delta_f)
        {
            var TerzFreq = new double[44] {1,1.25,1.6,2,2.5,3.15,4,5,6.3,8,10,12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,
                                            1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000};
            int n = 3;
            double red = 1.23;  // Energy Correction Hann
            int i = 0;
            double v1 = Math.Pow(10,0.05 * 3 / n) - Math.Pow(10,-0.05 * 3 / n);
            double fp = Math.Pow(10,(-1) * 3.0 / n / 10);
            double fu = Math.Pow(10,(-0.5) * 3 / n / 10);
            double fu1 = fu - 0.1 * fp * v1;
            double fu2 = fu + 0.1 * fp * v1;
            double fm, fo, fo1, fo2, sum, v2;
            int ku, ko, maxi, z1;
            var result = new XY_Data
            {
                xdata = TerzFreq,
                ydata = new double[TerzFreq.Length]
            };

            while (i < TerzFreq.Length)
            {
                fm = Math.Pow(10, i * 3.0 / n / 10);
                fo = Math.Pow(10,(i + 0.5) * 3 / n / 10);
                fo1 = fo - 0.1 * fm * v1;
                fo2 = fo + 0.1 * fm * v1;
                ku = (int)Math.Round(fu1 / delta_f - 0.5, MidpointRounding.AwayFromZero);
                ko = (int)Math.Round(fo2 / delta_f + 0.5, MidpointRounding.AwayFromZero);
                sum = 0;
                if (ku < 0) ku = 0;
                if (ko > z.Length - 1) ko = z.Length - 1;
                z1 = 0;
                for (int j = ku; j <= ko; j++)
                {
                    if (j * delta_f >= fu1 && j * delta_f <= fu2)
                    {
                        v2 = 0.5 * (1 - Math.Cos(Math.PI / (fu2 - fu1) * (j * delta_f - fu1)));
                        sum += Math.Pow(z[j],2) * v2;
                        z1++;
                    }
                    if (j * delta_f >= fu2 && j * delta_f <= fo1)
                    {
                        sum += Math.Pow(z[j],2);
                        z1++;
                    }
                    if (j * delta_f >= fo1 && j * delta_f <= fo2)
                    {
                        v2 = 0.5 * (1 + Math.Cos(Math.PI / (fo2 - fo1) * (j * delta_f - fo1)));
                        sum += Math.Pow(z[j] * v2, 2);
                        z1++;
                    }
                }
                //if (z1 != 0) maxi = i;
                result.ydata[i] = Math.Sqrt(sum) / red;
                i++;
                fu1 = fo1;
                fu2 = fo2;
            }
            
            maxi = 0;
            for (i = 0; i < result.ydata.Length; i++) if (result.ydata[i] == 0) maxi = i + 1;
            for (i = 0; i < result.ydata.Length - maxi; i++)
            {
                result.xdata[i] = result.xdata[i + maxi];
                result.ydata[i] = result.ydata[i + maxi];
            }
            Array.Resize(ref result.xdata, result.xdata.Length - maxi);
            Array.Resize(ref result.ydata, result.ydata.Length - maxi);
            return result;
        }

        public static double Add5(double input)
        {
            return input + 5;
        }

        public static double[] Evenness(double[][] audio, double fs, BackgroundWorker worker = null, double startprogress = 0, double overallprogress = 100)
        {
            if (audio == null) return null;
            double overlap = 0.5;
            double N_DOC = 16384;
            double percent = 0;
            worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
            double audio_length_seconds = audio[0].Length / fs;
            
            var params_Evenness = new Param_struct
            {
                delta_f = fs / N_DOC,
                delta_t = N_DOC / fs * (1 - overlap),
                average = 1,
                average_overlap = 0,
                fs = fs,
                windowtype = 1,
                DC = 1,
                y_amplitude = 0,
                y_axis = 0,
                dBref = 0.00002,
                freq_weight = 0,
                dsrange1 = 0,
                dsrange2 = audio_length_seconds,
                LQ = 0,
                y_unitconversion = 1
            };

            int channels = audio.Length;

            var Even2Ch = new double[channels][];
            var Evenness = new XY_Data();


            //var aps_y = new double[0];

            for (int ch = 0; ch < channels; ch++)
            {
                var aps = CalcAPSData_VOICE(audio[ch], new XY_Data(), 1 / fs, 0, new CATimeStamps(), params_Evenness);
                percent = (ch + 1) / channels * 50;
                worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
                if (ch == 0)
                {
                    /*aps_y = new double[aps.dimensions[1]];
                    Array.Copy(aps.ydata, aps_y, aps.dimensions[1]);*/
                    Evenness.xdata = new double[aps.dimensions[0]];
                    for (int ii = 0; ii < aps.dimensions[0]; ii++) Evenness.xdata[ii] = aps.xdata[ii * aps.dimensions[1]];
                }
                var aps_z = new double[aps.dimensions[1]];
                Even2Ch[ch] = new double[aps.dimensions[0]];

                for (int ii = 0; ii < aps.dimensions[0]; ii++)
                {
                    // Calculate 1/3 Octave Spectra from APS

                    Array.Copy(aps.zdata, ii * aps.dimensions[1], aps_z, 0, aps.dimensions[1]);
                    var Terz = GetTerzSpectrumFromAPS(aps_z, aps.ydata[1]);

                    double sum = 0;
                    int z = 0;
                    for (int k = 0; k < Terz.ydata.Length; k++)
                    {
                        if (Terz.xdata[k] >= 20 && Terz.xdata[k] <= 9000)
                        {
                            z++;
                            double y = 20 * Math.Log10(Terz.ydata[k] / 0.00002);
                            double L;
                            if (Terz.xdata[k] <= 20) L = 75;
                            else if (Terz.xdata[k] > 20 && Terz.xdata[k] <= 1000) L = 50 + (75 - 50) / (Math.Log10(20) - Math.Log10(1000)) * (Math.Log10(Terz.xdata[k]) - Math.Log10(1000));
                            else L = 20 + (50 - 20) / (Math.Log10(1000) - Math.Log10(6000)) * (Math.Log10(Terz.xdata[k]) - Math.Log10(6000));
                            sum += y - L;
                        }
                    }
                    double h = sum / z;
                    sum = 0;
                    z = 0;
                    for (int k = 0; k < Terz.ydata.Length; k++)
                    {
                        if (Terz.xdata[k] >= 20 && Terz.xdata[k] <= 9000)
                        {
                            z++;
                            double y = 20 * Math.Log10(Terz.ydata[k] / 0.00002);
                            double L;
                            if (Terz.xdata[k] <= 20) L = h + 75;
                            else if (Terz.xdata[k] > 20 && Terz.xdata[k] <= 1000) L = h + 50 + (75 - 50) / (Math.Log10(20) - Math.Log10(1000)) * (Math.Log10(Terz.xdata[k]) - Math.Log10(1000));
                            else L = h + 20 + (50 - 20) / (Math.Log10(1000) - Math.Log10(6000)) * (Math.Log10(Terz.xdata[k]) - Math.Log10(6000));
                            if (y > L) sum += y - L;
                        } 
                    }
                    if (z > 0) Even2Ch[ch][ii] = 5 - (sum / z);
                    else Even2Ch[ch][ii] = 0;
                    percent = 50 + (double)ii / aps.dimensions[0] * 45;
                    worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
                }
            }
            var Even_mean = new double[Even2Ch[0].Length];
            for (int ii = 0; ii < Even_mean.Length; ii++)
            {
                double summe = 0;
                for (int ch = 0; ch < channels; ch++) summe += Even2Ch[ch][ii];
                Even_mean[ii] = summe / channels;
            }
            Evenness.ydata = Interp1_linear(Evenness.xdata, Even_mean, DeltaArray(0.1, 0.025, audio_length_seconds - 0.1));
            Evenness.xdata = DeltaArray(0.1, 0.025, audio_length_seconds - 0.1);
            ExtrapolateNaNs(Evenness);
            percent = 100;
            worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
            return Evenness.ydata;
        }


        public static double[] ExtArticulationIndex(double[][] audio, double fs, BackgroundWorker worker = null, double startprogress = 0, double overallprogress = 100)
        {
            if (audio == null) return null;
            double overlap = 0.5;
            double N_DOC = 16384;
            double percent = 0;
            worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
            double audio_length_seconds = audio[0].Length / fs;

            // Definition ab Terz Nr. 23
            var cAI = new double[16] { 1, 2, 3.25, 4.25, 4.5, 5.25, 6.5, 7.25, 8.5, 11.5, 11, 9.5, 9, 7.75, 6.25, 2.5 };
            var cL = new double[16] { 34, 39, 41, 43, 45, 45, 45, 44, 42, 40, 37, 35, 33, 30, 26, 21 };

            var params_AI = new Param_struct
            {
                delta_f = fs / N_DOC,
                delta_t = N_DOC / fs * (1 - overlap),
                average = 1,
                average_overlap = 0,
                fs = fs,
                windowtype = 1,
                DC = 1,
                y_amplitude = 0,
                y_axis = 0,
                dBref = 0.00002,
                freq_weight = 0,
                dsrange1 = 0,
                dsrange2 = audio_length_seconds,
                LQ = 0,
                y_unitconversion = 1
            };

            int channels = audio.Length;

            var ArtI2Ch = new double[channels][];
            var ArtI = new XY_Data();


            //var aps_y = new double[0];

            for (int ch = 0; ch < channels; ch++)
            {
                var aps = CalcAPSData_VOICE(audio[ch], new XY_Data(), 1 / fs, 0, new CATimeStamps(), params_AI);
                percent = (ch + 1) / channels * 50;
                worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
                if (ch == 0)
                {
                    /*aps_y = new double[aps.dimensions[1]];
                    Array.Copy(aps.ydata, aps_y, aps.dimensions[1]);*/
                    ArtI.xdata = new double[aps.dimensions[0]];
                    for (int ii = 0; ii < aps.dimensions[0]; ii++) ArtI.xdata[ii] = aps.xdata[ii * aps.dimensions[1]];
                }
                var aps_z = new double[aps.dimensions[1]];
                ArtI2Ch[ch] = new double[aps.dimensions[0]];
               
                for (int ii = 0; ii < aps.dimensions[0]; ii++)
                {
                    // Calculate 1/3 Octave Spectra from APS
                   
                    Array.Copy(aps.zdata, ii * aps.dimensions[1], aps_z, 0, aps.dimensions[1]);
                    var Terz = GetTerzSpectrumFromAPS(aps_z, aps.ydata[1]);
                    int i = 0;
                    while (Math.Round(Terz.xdata[i],MidpointRounding.AwayFromZero) < 200 && i < Terz.xdata.Length - 1) i++;
                    double sum = 0;
                    int tnr = 23;
                    while (i < Terz.xdata.Length && Math.Round(Terz.xdata[i],MidpointRounding.AwayFromZero) <= 6300)
                    {
                        var res = cAI[tnr - 23] * (1 + cL[tnr - 23] / 30 - (20 * Math.Log10(Terz.ydata[i] / 0.00002)) / 30);
                        sum += res;
                        i++;
                        tnr++;
                    }
                    if (tnr != 39) ArtI2Ch[ch][ii] = Double.PositiveInfinity;
                    else ArtI2Ch[ch][ii] = sum;
                    percent = 50 + (double)ii / aps.dimensions[0] * 45;
                    worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
                }
            }
            var AI_mean = new double[ArtI2Ch[0].Length];
            for (int ii = 0; ii < AI_mean.Length; ii++)
            {
                double summe = 0;
                for (int ch = 0; ch < channels; ch++) summe += ArtI2Ch[ch][ii];
                AI_mean[ii] = summe / channels;
            }
            ArtI.ydata = Interp1_linear(ArtI.xdata, AI_mean, DeltaArray(0.1, 0.025, audio_length_seconds - 0.1));
            ArtI.xdata = DeltaArray(0.1, 0.025, audio_length_seconds - 0.1);
            ExtrapolateNaNs(ArtI);
            percent = 100;
            worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
            return ArtI.ydata;
        }

        public static double[] DOC (double[][] audio, XY_Data rpm, double fs, int powersport, BackgroundWorker worker = null, double startprogress = 0, double overallprogress = 100)
        {
            if (rpm.xdata == null || rpm.ydata == null || audio == null) return null;
            // powersport: 0 : DOC Power, 1: DOC Sport
            double overlap = 0.5;
            double N_DOC = 16384;
            double percent = 0;
            worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));

            double audio_length_seconds = audio[0].Length / fs;

            var params_DOC = new Param_struct
            {
                delta_f = fs / N_DOC,
                delta_t = N_DOC / fs * (1 - overlap),
                average = 1,
                average_overlap = 0,
                fs = fs,
                windowtype = 1,
                DC = 1,
                y_amplitude = 0,
                y_axis = 0,
                dBref = 0.00002,
                freq_weight = 1,
                dsrange1 = 0,
                dsrange2 = audio_length_seconds,
                LQ = 0,
                y_unitconversion = 1
            };

            int channels = audio.Length;

            var DOC2ch = new double[channels][];
            var DOCP = new XY_Data();

            
            var aps_y = Array.Empty<double>();

            for (int ch = 0; ch < channels; ch++)
            {
                var aps = CalcAPSData(audio[ch], new XY_Data(), 1 / fs, 0, new CATimeStamps(), params_DOC);
                percent = (ch + 1) / channels * 50;
                worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
                if (ch == 0)
                {
                    aps_y = new double[aps.dimensions[1]];
                    Array.Copy(aps.ydata, aps_y, aps.dimensions[1]);
                    DOCP.xdata = new double[aps.dimensions[0]];
                    for (int ii = 0; ii < aps.dimensions[0]; ii++) DOCP.xdata[ii] = aps.xdata[ii * aps.dimensions[1]];
                }
                var aps_z = new double[aps.dimensions[1]];
                DOC2ch[ch] = new double[aps.dimensions[0]];
                var t1 = DeltaArray(0, 1, aps.dimensions[0] - 1);
                var t2 = new double[t1.Length];
                for (int ii = 0; ii < aps.dimensions[0]; ii++)
                {
                    // Calculate average rpm of FFT block
                    t1[ii] = Math.Round(t1[ii] * N_DOC * (1 - overlap), MidpointRounding.AwayFromZero) / fs;
                    t2[ii] = t1[ii] + N_DOC / fs;
                    /*var i1 = Array.FindIndex(rpm.xdata, w => w >= t1[ii]);
                    var i2 = Array.FindLastIndex(rpm.xdata, w => w <= t2[ii]);
                    var interval = new double[i2 - i1 + 1];
                    Array.Copy(rpm.ydata, interval, interval.Length);
                    var speed = interval.Sum() / interval.Length;*/

                    int z = 0;
                    double sum = 0;
                    int i1 = 0;
                    while (i1 < rpm.xdata.Length - 2 && rpm.xdata[i1 + 1] < t1[ii]) i1++;
                    int i = i1;
                    while (i < rpm.xdata.Length && rpm.xdata[i] <= t2[ii])
                    {
                        z++;
                        sum += rpm.ydata[i];
                        i++;
                    }
                    double speed;
                    if (z > 2) speed = sum / z;
                    //else speed = rpm.ydata[i1] + (rpm.ydata[i1 + 1] - rpm.ydata[i1]) / (rpm.xdata[i1 + 1] - rpm.xdata[i1]) * ((t2[ii] + t1[ii]) / 2 - rpm.xdata[i1]);
                    else speed = GetYfromXY(rpm, (t2[ii] + t1[ii]) / 2);

                    // Do Special Order Analysis
                    Array.Copy(aps.zdata, ii * aps.dimensions[1], aps_z, 0, aps.dimensions[1]);
                    DOC2ch[ch][ii] = DOC_OrderAnalysis(aps_y, aps_z, speed, powersport);
                    percent = 50 + (double)ii / aps.dimensions[0] * 45;
                    worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
                }
            }

            var DOCP_mean = new double[DOC2ch[0].Length];
            for (int ii = 0; ii < DOCP_mean.Length; ii++)
            {
                double summe = 0;
                for (int ch = 0; ch < channels; ch++) summe += DOC2ch[ch][ii];
                DOCP_mean[ii] = summe / channels;
            }
            DOCP.ydata = Interp1_linear(DOCP.xdata, DOCP_mean, DeltaArray(0.1, 0.025, audio_length_seconds - 0.1));
            DOCP.xdata = DeltaArray(0.1, 0.025, audio_length_seconds - 0.1);
            ExtrapolateNaNs(DOCP);
            for (int ii = 0; ii < DOCP.ydata.Length; ii++) if (DOCP.ydata[ii] < 0.00002) DOCP.ydata[ii] = 0.00002;
            for (int ii = 0; ii < DOCP.ydata.Length; ii++) DOCP.ydata[ii] = 20 * Math.Log10(DOCP.ydata[ii] / 0.00002);
            percent = 100;
            worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
            return DOCP.ydata;
        }


        public static double[] Sportiness(double[] time, double[] DOCSport, double[] ESC, double[] Tonality)
        {
            

            if (time == null || DOCSport == null || Tonality == null || ESC == null) return null;
            var result = new XY_Data
            {
                xdata = time,
                ydata = new double[time.Length]
            };
            for (int ii = 0; ii < time.Length; ii++) result.ydata[ii] = 1.288 * DOCSport[ii] + 23.7 * Tonality[ii] + 0.007421 * ESC[ii] - 70.09;
            return result.ydata;
        }

        public static double[] Sportiness_old(double[] time, double[] DOCSport, double[] ESC, double[] Tonality)
        {
            if (time == null || DOCSport == null || Tonality == null || ESC == null) return null;
            var result = new XY_Data
            {
                xdata = time,
                ydata = new double[time.Length]
            };
            for (int ii = 0; ii < time.Length; ii++) result.ydata[ii] = 1.288 * DOCSport[ii] + 23.7 * Tonality[ii] + 0.007421 * ESC[ii] - 70.09;
            return result.ydata;
        }

        public static XYZ_Data APSVOICE(double[][] audio_mch, double fs, double delta_t, double delta_f, XY_Data rpm = new XY_Data())
        {
            

            int LQ = 0;
            
            if (rpm.xdata != null)
            {
                LQ = 1;
            }
            

            var audio = new double[audio_mch[0].Length];
            int channels = audio_mch.Length;
            if (channels > 1)
            {
                for (int ch = 0; ch < channels; ch++)
                {
                    for (int ii = 0; ii < audio.Length; ii++) audio[ii] += audio_mch[ch][ii] / channels;
                }
            }
            
            
            

            double audio_length_seconds = audio.Length / fs;

            var params_ALevel = new Param_struct
            {
                delta_f = delta_f,
                delta_t = delta_t,
                //delta_t = 0.25,
                average = 1,
                average_overlap = 0,
                fs = fs,
                windowtype = 1,
                DC = 0,
                y_amplitude = 0,
                y_axis = 1,
                dBref = 0.00002,
                freq_weight = 1,
                dsrange1 = 0,
                dsrange2 = audio_length_seconds - (1 / fs),
                LQ = LQ,
                y_unitconversion = 1
            };

            var aps = CalcAPSData(audio, rpm, 1 / fs, 0, new CATimeStamps(), params_ALevel);
            
           
            //SaveSomeDoubles("C:\\Users\\u12o24\\Documents\\test1.bin",aps.zdata);
            return aps;
        }

        public static double[] ABCLinearLevel(double[][] audio, double fs, int weighting, BackgroundWorker worker = null, double startprogress = 0, double overallprogress = 100)
        {
            if (audio == null) return null;
            
            double overlap = 0.5;
            double N_ALevel = 8192;
            double percent = 0;
            worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
            double audio_length_seconds = audio[0].Length / fs;

            var params_ALevel = new Param_struct
            {
                delta_f = fs / N_ALevel,
                delta_t = N_ALevel / fs * (1 - overlap),
                //delta_t = 0.25,
                average = 1,
                average_overlap = 0,
                fs = fs,
                windowtype = 1,
                DC = 0,
                y_amplitude = 0,
                y_axis = 1,
                dBref = 0.00002,
                freq_weight = weighting,
                dsrange1 = 0,
                dsrange2 = audio_length_seconds,
                LQ = 0,
                y_unitconversion = 1
            };

            int channels = audio.Length;

            var ALevel2ch = new double[channels][];
            var ALevel = new XY_Data();


            //var aps_y = new double[0];

            for (int ch = 0; ch < channels; ch++)
            {
                var aps = CalcOverallLevel(audio[ch], new XY_Data(), 1 / fs, 0, new CATimeStamps(), params_ALevel);
                if (ch == 0) ALevel.xdata = aps.xdata;
                ALevel2ch[ch] = aps.zdata;
                percent = (ch+1)/channels * 95;
                worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
            }

            var ALevel_mean = new double[ALevel2ch[0].Length];
            for (int ii = 0; ii < ALevel_mean.Length; ii++)
            {
                double summe = 0;
                for (int ch = 0; ch < channels; ch++) summe += ALevel2ch[ch][ii];
                ALevel_mean[ii] = summe / channels;
            }
            ALevel.ydata = Interp1_linear(ALevel.xdata, ALevel_mean, DeltaArray(0.1, 0.025, audio_length_seconds - 0.1));
            ALevel.xdata = DeltaArray(0.1, 0.025, audio_length_seconds - 0.1);
            ExtrapolateNaNs(ALevel);
            //for (int ii = 0; ii < ALevel.ydata.Length; ii++) ALevel.ydata[ii] = 20 * Math.Log10(ALevel.ydata[ii] / 0.00002);
            percent = 100;
            worker.ReportProgress((int)(startprogress + overallprogress * percent / 100));
            return ALevel.ydata;
        }

        public static double[] Luxury(double[] time, double[] ExtArtIndex, double[] Evenness, double[] Tonality)
        {
            if (time == null || ExtArtIndex == null || Evenness == null || Tonality == null) return null;
            var result = new XY_Data
            {
                xdata = time,
                ydata = new double[time.Length]
            };
            for (int ii = 0; ii < time.Length; ii++) result.ydata[ii] = 0.44 * ExtArtIndex[ii] + 2.825 * Evenness[ii] - 27.5 * Tonality[ii] - 14.125;
            return result.ydata;
        }

        public static double[] Luxury_old(double[] time, double[] ExtArtIndex, double[] Evenness, double[] Tonality)
        {
            if (time == null || ExtArtIndex == null || Evenness == null || Tonality == null) return null;
            var result = new XY_Data
            {
                xdata = time,
                ydata = new double[time.Length]
            };
            for (int ii = 0; ii < time.Length; ii++) result.ydata[ii] = 0.44 * ExtArtIndex[ii] + 2.825 * Evenness[ii] - 27.5 * Tonality[ii] - 14.125;
            return result.ydata;
        }

        public static double[] Powerfulness_new(double[] time, double[] EngineSpeed, double[] ALevel, double[] Evenness, double[] ExtAI)
        {
            if (time == null || EngineSpeed == null || ALevel == null || Evenness == null || ExtAI == null) return null;
            var result = new XY_Data
            {
                xdata = time,
                ydata = new double[time.Length]
            };
            for (int ii = 0; ii < time.Length; ii++) result.ydata[ii] = -116.1909 + 0.0149 * EngineSpeed[ii] + 2.9123 * ALevel[ii] - 9.9332 * Evenness[ii] - 0.6223 * ExtAI[ii];
            return result.ydata;
        }

        public static double[] Powerfulness(double[] time, double[] ESC, double[] DOCP, double[] REO, double[] Tonality)
        {
            if (time == null || ESC == null || DOCP == null || REO == null || Tonality == null) return null;
            var result = new XY_Data
            {
                xdata = time,
                ydata = new double[time.Length]
            };
            for (int ii = 0; ii < time.Length; ii++) result.ydata[ii] = - 30 + 0.5 * DOCP[ii] + 33 * Tonality[ii] + 0.007 * ESC[ii] + 70 * REO[ii];
            return result.ydata;
        }

        public static double[] Pleasantness(double[] time, double[] EngineSpeed, double[] ESC, double[] REO, double[] LinearLevel, double[] DOCS, double[] DOCP, double[] Tonality)
        {
            if (time == null || EngineSpeed == null || LinearLevel == null || ESC == null || REO == null || DOCS== null || DOCP == null || Tonality == null) return null;
            var result = new XY_Data
            {
                xdata = time,
                ydata = new double[time.Length]
            };
            var sport = Sportiness(time, DOCS, ESC, Tonality);
            var power = Powerfulness(time, ESC, DOCP, REO, Tonality);
            for (int ii = 0; ii < time.Length; ii++) result.ydata[ii] = 22.845 + 0.0029 * EngineSpeed[ii] - 0.1214 * LinearLevel[ii] - 0.1636 * sport[ii] - 0.7027 * power[ii];
            //for (int ii = 0; ii < time.Length; ii++) result.ydata[ii] = -75.145 + 0.0047 * EngineSpeed[ii] + 1.3908 * LinearLevel[ii] -0.0185*DOCS[ii]-0.3413*Tonality[ii] - 0.00010686*ESC[ii] - 1.2834 * power[ii];
            //for (int ii = 0; ii < time.Length; ii++) result.ydata[ii] = -36.643 + 0.0047 * EngineSpeed[ii] + 1.3908 * LinearLevel[ii] - 0.0185 * DOCS[ii] - 0.6417 * DOCP[ii] - 42.6935 * Tonality[ii] - 0.0091 * ESC[ii] - 89.838 * REO[ii];
            return result.ydata;
        }

        public static double[] VOICEParamEngineSpeed(double[] time, XY_Data rpm)
        {
            var rpm_ip = new XY_Data
            {
                xdata = time,
                ydata = Interp1_linear(rpm.xdata, rpm.ydata, time)
            };

            //rpm_ip = ExtrapolateNaNs(rpm_ip);
            for (int ii = 0; ii < rpm_ip.ydata.Length; ii++) if (double.IsNaN(rpm_ip.ydata[ii])) rpm_ip.ydata[ii] = 0;

            return rpm_ip.ydata;
        }

        public static double[] VOICEParamTime(double audio_length_seconds)
        {
            var result = DeltaArray(0.1, 0.025, audio_length_seconds - 0.1);
            return result;
        }

        public static double[] EngineSpeedChange(XY_Data rpm, double audio_length_seconds)
        {
            if (rpm.xdata == null || rpm.ydata == null) return null;
            double timestep = 0.1;
            double timeinterval = 0.5;
            var timeaxis = DeltaArray(0, timestep, rpm.xdata[rpm.xdata.Length - 1] - timeinterval);
            int L = timeaxis.Length;
            var ESC = new XY_Data
            {
                xdata = new double[L],
                ydata = new double[L]
            };
            double starttime = 0;
            double endtime = timeinterval;
            int z = 0;
            while (starttime < rpm.xdata[rpm.xdata.Length - 1] && endtime < rpm.xdata[rpm.xdata.Length - 1] && z < ESC.xdata.Length)
            {
                var startspeed = GetYfromXY(rpm, starttime);
                var endspeed = GetYfromXY(rpm, endtime);
                ESC.xdata[z] = (endtime + starttime) / 2;
                ESC.ydata[z] = (endspeed - startspeed) / timeinterval;
                z ++;
                starttime += timestep;
                endtime += timestep;
            }

            timestep = 0.025;
            int smoothingpoints = 11;
            var xq = DeltaArray(0.1, timestep, audio_length_seconds - 0.1);

            int nulls = 0;
            for (int ii = ESC.xdata.Length - 1; ii > 0; ii--) if (ESC.xdata[ii] == 0) nulls++;
            Array.Resize(ref ESC.xdata, ESC.xdata.Length - nulls);
            Array.Resize(ref ESC.ydata, ESC.xdata.Length);

            var ESCip = new XY_Data
            {
                xdata = xq,
                ydata = Interp1_linear(ESC.xdata, ESC.ydata, xq)
            };

            ESCip = ExtrapolateNaNs(ESCip);

            var ESCsmooth = Moving(ESCip.ydata, smoothingpoints);
            ESC.xdata = xq;
            ESC.ydata = ESCsmooth;

            for (int ii = 0; ii < smoothingpoints / 2; ii++) ESC.ydata[ii] = ESCip.ydata[ii];
            for (int ii = ESC.ydata.Length - smoothingpoints / 2; ii < ESC.ydata.Length; ii++) ESC.ydata[ii] = ESCip.ydata[ii];
            return ESC.ydata;
        }

        public static XY_Data ExtrapolateNaNs(XY_Data xy)
        {
            if (xy.xdata == null || xy.ydata == null) return new XY_Data();
            // Post processing function for Interp1_linear to replace the NaNs outside grid with extrapolatet values
            // find first NaN
            int ii = 0;
            while (Double.IsNaN(xy.ydata[ii]) && ii < xy.ydata.Length) ii++;
            int gridindex1 = ii;
            ii = xy.ydata.Length - 1;
            while (Double.IsNaN(xy.ydata[ii]) && ii >= 0) ii--;
            int gridindex2 = ii;
            var xy_red = new XY_Data
            {
                xdata = RemoveNaN(xy.xdata, xy.ydata),
                ydata = RemoveNaN(xy.ydata, xy.ydata)
            };
            
            for (ii = 0; ii < gridindex1; ii++) xy.ydata[ii] = GetYfromXY(xy_red, xy.xdata[ii]);
            for (ii = gridindex2 + 1; ii < xy.ydata.Length; ii++) xy.ydata[ii] = GetYfromXY(xy_red, xy.xdata[ii]);
            return xy;
        }


        public static XY_Data Interpolate_smoothed(double[] x, double[] y, double[] xq)
        {
            //var result = new double[xq.Length];
            var result = new XY_Data
            {
                xdata = new double[xq.Length],
                ydata = new double[xq.Length]
            };
            double xrange1;
            double xrange2;
            for (int ii = 0; ii < xq.Length; ii++)
            {
                if (ii == 0) xrange1 = xq[ii];
                else xrange1 = xq[ii] - (xq[ii] - xq[ii - 1]) / 2;
                if (ii == xq.Length - 1) xrange2 = xq[ii];
                else xrange2 = xq[ii] + (xq[ii + 1] - xq[ii]) / 2;
                var xind1 = Array.FindIndex(x, f => f >= xrange1);
                var xind2 = Array.FindLastIndex(x, f => f <= xrange2);
                if (xind1 == -1 || xind2 == -1 || xind2 < xind1) result.ydata[ii] = Double.NaN;
                else
                {
                    double summe = 0;
                    for (int jj = xind1; jj <= xind2; jj++) summe += y[jj];
                    result.ydata[ii] = summe / (xind2 - xind1 + 1);
                }
                
            }
            var xwoNaN = RemoveNaN(xq, result.ydata);
            var ywoNaN = RemoveNaN(result.ydata, result.ydata);

            var yip = Interp1_linear(xwoNaN, ywoNaN, xq);
            result.xdata = RemoveNaN(xq, yip);
            result.ydata = RemoveNaN(yip, yip);

            return result;
        }

        public static double[] Interp1_linear(double[] x, double[] y, double[] xq)
        {
            if (x == null || y == null || xq == null) return null;
            if (x.Length != y.Length || x.Length == 0 || xq.Length == 0) return null;
            double direction = x[x.Length - 1] - x[0];

            double[] yq = new double[xq.Length];
            double[] xq1 = new double[xq.Length];

            var spline = MathNet.Numerics.Interpolation.LinearSpline.InterpolateSorted(x,y);


            int L_singlepackage = xq.Length / 16;
            /*System.Diagnostics.Stopwatch tic = new System.Diagnostics.Stopwatch();
            tic.Start();*/
            //for (int i = 0; i < 16; i++)
            Parallel.For(0, 16, i =>
            {
                for (int ii = 0; ii < L_singlepackage; ii++)
                {
                    if (xq[i * L_singlepackage + ii] >= x[0] && xq[i * L_singlepackage + ii] <= x[x.Length - 1])
                    {
                        yq[i * L_singlepackage + ii] = spline.Interpolate(xq[i * L_singlepackage + ii]);

                    }
                    else yq[i * L_singlepackage + ii] = Double.NaN;
                }
            }
            );
            for (int ii = L_singlepackage * 16; ii < xq.Length; ii++)
            {
                if (xq[ii] >= x[0] && xq[ii] <= x[x.Length - 1])
                {
                    yq[ii] = spline.Interpolate(xq[ii]);

                }
                else yq[ii] = Double.NaN;
            }
            /*
            for (int ii = 0; ii < xq.Length; ii++)
                {
                    if (xq[ii] >= x[0] && xq[ii] <= x[x.Length - 1])
                    {
                        yq[ii] = spline.Interpolate(xq[ii]);
                    }
                    else yq[ii] = Double.NaN;
                }
            */

            return yq;
        }

        public static double[] RemoveNaN(double[] input, double[] nanref)
        {
            if (input == null || nanref == null) return null;
            if (input.Length != nanref.Length || input.Length == 0) return null;
            int jj = 0;
            double[] output = new double[input.Length];
            for (int ii = 0; ii < input.Length; ii++)
            {
                if (!Double.IsNaN(nanref[ii]))
                {
                    output[jj] = input[ii];
                    jj++;
                }
            }
            Array.Resize(ref output, jj);
            return output;
        }
        public static int[] RemoveNaN(int[] input, double[] nanref)
        {
            if (input == null || nanref == null) return null;
            if (input.Length != nanref.Length || input.Length == 0) return null;
            int jj = 0;
            int[] output = new int[input.Length];
            for (int ii = 0; ii < input.Length; ii++)
            {
                if (!Double.IsNaN(nanref[ii]))
                {
                    output[jj] = input[ii];
                    jj++;
                }
            }
            Array.Resize(ref output, jj);
            return output;
        }


        /* public int FindIndex(XY_Data xy, double x, int roundingmode)
         {
             // Calculates average of the interval between x1 and x2
             if (roundingmode == 0) Array.FindIndex()

         }*/

        private static double GetYfromXY(XY_Data xy, double x)
        {
            // Interpolates or Extrapolates a Y value to a given X value on a XY grid of data
            int i;
            if (x < xy.xdata[xy.xdata.Length - 1])
            {
                i = 0;
                while (x > xy.xdata[i]) i++;
                if (i == 0) i = 1;
            }
            else i = xy.xdata.Length - 1;

            double result = xy.ydata[i - 1] + (xy.ydata[i] - xy.ydata[i - 1]) / (xy.xdata[i] - xy.xdata[i - 1]) * (x - xy.xdata[i - 1]);
            

            return result;
        }

        private static int Nextpow2(double x)
        {
            int result = (int)Math.Round(Math.Pow(2, Math.Ceiling(Math.Log(x) / Math.Log(2))),MidpointRounding.AwayFromZero);
            return result;
        }

        public static XY_Data ProcessLQ_mono(XY_Data lq)
        {
            // Test if increasing or decreasing
            int gradient = 0;
            int iMin = 0;
            int iMax = 0;
            double minvalue = lq.ydata[0];
            double maxvalue = lq.ydata[0];

            for (int ii = 1; ii < lq.xdata.Length; ii++)
            {
                if (lq.ydata[ii] - lq.ydata[ii - 1] > 0) gradient++;
                else if (lq.ydata[ii] - lq.ydata[ii - 1] < 0) gradient--;
                if (lq.ydata[ii] > maxvalue)
                {
                    maxvalue = lq.ydata[ii];
                    iMax = ii;
                }
                if (lq.ydata[ii] < minvalue)
                {
                    minvalue = lq.ydata[ii];
                    iMin = ii;
                }
            }
            
            int jj = 0;
            double minmaxval;
            int istart;

            XY_Data lq_mono = new XY_Data
            {
                xdata = new double[lq.xdata.Length],
                ydata = new double[lq.ydata.Length]
            };
            lq_mono.xdata[0] = lq.xdata[0];
            lq_mono.ydata[0] = lq.ydata[0];

            if (gradient >= 0)
            {
                minmaxval = minvalue;
                istart = iMin;
                for (int ii = 0; ii < lq.ydata.Length; ii++)
                {
                    if (lq.ydata[ii] > minmaxval && ii > istart)
                    {
                        minmaxval = lq.ydata[ii];
                        lq_mono.xdata[jj] = lq.xdata[ii];
                        lq_mono.ydata[jj] = lq.ydata[ii];
                        jj++;
                    }
                }
            }
            else
            {
                minmaxval = maxvalue;
                istart = iMax;
                for (int ii = 0; ii < lq.ydata.Length; ii++)
                {
                    if (lq.ydata[ii] < minmaxval && ii > istart)
                    {
                        minmaxval = lq.ydata[ii];
                        lq_mono.xdata[jj] = lq.xdata[ii];
                        lq_mono.ydata[jj] = lq.ydata[ii];
                        jj++;
                    }
                }
            }
            lq.xdata = new double[(jj - 1)];
            lq.ydata = new double[(jj - 1)];
            for (int ii = 0; ii < lq.xdata.Length; ii++)
            {
                if (gradient < 0)
                {
                    lq.xdata[ii] = lq_mono.xdata[lq.xdata.Length - 1 - ii];
                    lq.ydata[ii] = lq_mono.ydata[lq.xdata.Length - 1 - ii];
                }
                    else
                {
                    lq.xdata[ii] = lq_mono.xdata[ii];
                    lq.ydata[ii] = lq_mono.ydata[ii];
                }
            }
            

            
            return lq;
        }

        public static XY_Data SortandInterpolateData_old(XY_Data xydata, double resolution)
        {
            // Reads a set of X and Y Data, sorts it monotonically and interpolates it to an equidistant resolution
            var xysorted = new XY_Data
            {
                xdata = new double[xydata.xdata.Length],
                ydata = new double[xydata.ydata.Length]
            };
            int gradient = 0;
            for (int j = 1; j < xydata.xdata.Length; j++)
            {
                if (xydata.xdata[j] - xydata.xdata[j - 1] > 0) gradient++;
                else if (xydata.xdata[j] - xydata.xdata[j - 1] < 0) gradient--;
            }
            int j2 = 0;
            double minmaxval = xydata.xdata[0];
            xysorted.xdata[0] = xydata.xdata[0];
            xysorted.ydata[0] = xydata.ydata[0];
            if (gradient >= 0)
            {
                for (int j = 1; j < xydata.xdata.Length;j++)
                {
                    if(xydata.xdata[j] > minmaxval)
                    {
                        minmaxval = xydata.xdata[j];
                        j2++;
                        xysorted.xdata[j2] = xydata.xdata[j];
                        xysorted.ydata[j2] = xydata.ydata[j];
                    }
                }
            }
            else
            {
                for (int j = 1; j < xydata.xdata.Length;j++)
                {
                    if (xydata.xdata[j] < minmaxval)
                    {
                        minmaxval = xydata.xdata[j];
                        j2++;
                        xysorted.xdata[j2] = xydata.xdata[j];
                        xysorted.ydata[j2] = xydata.ydata[j];
                    }
                }
            }
            Array.Resize(ref xysorted.xdata, j2 + 1);
            Array.Resize(ref xysorted.ydata, j2 + 1);

            if (gradient < 0) Array.Reverse(xysorted.xdata);
            if (gradient < 0) Array.Reverse(xysorted.ydata);

            double[] newxdata;
            double[] newydata;

            if (xysorted.ydata.Length > 1)
            {
                double firstval = Math.Ceiling(xysorted.xdata[0] / resolution) * resolution;

                newxdata = DeltaArray(firstval, resolution, xysorted.xdata[xysorted.xdata.Length - 1]);

                newydata = Interp1_linear(xysorted.xdata, xysorted.ydata, newxdata);
            }
            else
            {

                double firstval = Math.Ceiling(xysorted.xdata[0] / resolution) * resolution;
                newxdata = new double[1] { firstval };
                newydata = new double[1] { xysorted.ydata[0] };
            }

                    //var ipval = Math.Round((xysorted.xdata[xysorted.xdata.Length - 1] - xysorted.xdata[0]) / (xysorted.xdata.Length - 1)) + 1;
                    
            var result = new XY_Data
            {
                xdata = newxdata,
                ydata = newydata
            };
            return result;
            
            /*
        if Form_Main.Form_QueryControl.Monotonic1.Checked then begin
            //
          val(Form_Main.Form_QueryControl.Edit_Smooth.Text,ipval,code);
          if xpnr <> P_time then begin
            if ipval < 1 then begin
              ipval := 1;
              Form_Main.Form_QueryControl.Edit_Smooth.Text := '1';
            end
          end else if ipval < 0.1 then begin
            ipval := 0.1;
            Form_Main.Form_QueryControl.Edit_Smooth.Text := '0.1';
          end;

          if round((xysorted^.x[high(xysorted^.x)]-xysorted^.x[0]) / ipval) + 1.0 > high(xysorted^.x) then begin
            Form_Main.Form_QueryControl.Edit_Smooth.Text := inttostr(round((xysorted^.x[high(xysorted^.x)]-xysorted^.x[0])/high(xysorted^.x))+1);
            ipval := round((xysorted^.x[high(xysorted^.x)]-xysorted^.x[0])/high(xysorted^.x))+1;
          end;

          newxy(xyip);
          //test :=  high(xysorted^.x);
          xyip := InterExtraPolate(xysorted,xysorted^.x[0],xysorted^.x[high(xysorted^.x)],ipval);

          SetLength(Curve[z]^.x,high(xyip^.x)+1);
          SetLength(Curve[z]^.y,high(xyip^.y)+1);
            //

          for j := 0 to high(xyip^.x) do begin
            Curve[z]^.x[j] := xyip^.x[j];
            Curve[z]^.y[j] := xyip^.y[j];

          end;

        end;
             * */
        }

        public static XY_Data SortandInterpolateData(XY_Data xydata, double resolution)
        {
            // Reads a set of X and Y Data, sorts it monotonically and interpolates it to an equidistant resolution

            var xydata_red = new XY_Data();
            //double zeroval;
            int zerovalfound = 0;
            for (int ii = 0; ii < xydata.xdata.Length; ii++)
            {
                if (xydata.xdata[ii] == 0)
                {
                    zerovalfound = 1;
                    //zeroval = xydata.ydata[ii];
                    break;
                }
            }
            int newlength = 0;
            for (int ii = 0; ii < xydata.xdata.Length; ii++) if (xydata.xdata[ii] != 0) newlength++;

            xydata_red.xdata = new double[newlength + zerovalfound];
            xydata_red.ydata = new double[newlength + zerovalfound];
            int jj = 0;
            for (int ii = 0; ii < xydata.xdata.Length; ii++)
            {
                if (xydata.xdata[ii] != 0)
                {
                    xydata_red.xdata[jj] = xydata.xdata[ii];
                    xydata_red.ydata[jj] = xydata.ydata[ii];
                    jj++;
                }
            }


            var xysorted = new XY_Data
            {
                xdata = new double[xydata_red.xdata.Length],
                ydata = new double[xydata_red.ydata.Length]
            };
            Array.Copy(xydata_red.xdata, xysorted.xdata, xysorted.xdata.Length);
            Array.Copy(xydata_red.ydata, xysorted.ydata, xysorted.ydata.Length);
            Array.Sort(xysorted.xdata, xysorted.ydata);
            /*var xysorted = new XY_Data
            {
                xdata = new double[xydata.xdata.Length],
                ydata = new double[xydata.ydata.Length]
            };
            int gradient = 0;
            for (int j = 1; j < xydata.xdata.Length; j++)
            {
                if (xydata.xdata[j] - xydata.xdata[j - 1] > 0) gradient++;
                else if (xydata.xdata[j] - xydata.xdata[j - 1] < 0) gradient--;
            }
            int j2 = 0;
            double minmaxval = xydata.xdata[0];
            xysorted.xdata[0] = xydata.xdata[0];
            xysorted.ydata[0] = xydata.ydata[0];
            if (gradient >= 0)
            {
                for (int j = 1; j < xydata.xdata.Length; j++)
                {
                    if (xydata.xdata[j] > minmaxval)
                    {
                        minmaxval = xydata.xdata[j];
                        j2++;
                        xysorted.xdata[j2] = xydata.xdata[j];
                        xysorted.ydata[j2] = xydata.ydata[j];
                    }
                }
            }
            else
            {
                for (int j = 1; j < xydata.xdata.Length; j++)
                {
                    if (xydata.xdata[j] < minmaxval)
                    {
                        minmaxval = xydata.xdata[j];
                        j2++;
                        xysorted.xdata[j2] = xydata.xdata[j];
                        xysorted.ydata[j2] = xydata.ydata[j];
                    }
                }
            }
            Array.Resize(ref xysorted.xdata, j2 + 1);
            Array.Resize(ref xysorted.ydata, j2 + 1);

            if (gradient < 0) Array.Reverse(xysorted.xdata);
            if (gradient < 0) Array.Reverse(xysorted.ydata);
            */
            
            double[] newxdata;
            double[] newydata;

            if (xysorted.ydata.Length > 1)
            {
                double firstval = Math.Ceiling(xysorted.xdata[0] / resolution) * resolution;

                newxdata = DeltaArray(firstval, resolution, xysorted.xdata[xysorted.xdata.Length - 1]);


                //newydata = Interpolate_smoothed(xysorted.xdata, xysorted.ydata, newxdata);
                var newxydata = Interpolate_smoothed(xysorted.xdata, xysorted.ydata, newxdata);
                newxdata = newxydata.xdata;
                newydata = newxydata.ydata;
            }
            else
            {

                double firstval = Math.Ceiling(xysorted.xdata[0] / resolution) * resolution;
                newxdata = new double[1] { firstval };
                newydata = new double[1] { xysorted.ydata[0] };
            }

            //var ipval = Math.Round((xysorted.xdata[xysorted.xdata.Length - 1] - xysorted.xdata[0]) / (xysorted.xdata.Length - 1)) + 1;

            var result = new XY_Data
            {
                xdata = newxdata,
                ydata = newydata
            };
            return result;

            /*
        if Form_Main.Form_QueryControl.Monotonic1.Checked then begin
            //
          val(Form_Main.Form_QueryControl.Edit_Smooth.Text,ipval,code);
          if xpnr <> P_time then begin
            if ipval < 1 then begin
              ipval := 1;
              Form_Main.Form_QueryControl.Edit_Smooth.Text := '1';
            end
          end else if ipval < 0.1 then begin
            ipval := 0.1;
            Form_Main.Form_QueryControl.Edit_Smooth.Text := '0.1';
          end;

          if round((xysorted^.x[high(xysorted^.x)]-xysorted^.x[0]) / ipval) + 1.0 > high(xysorted^.x) then begin
            Form_Main.Form_QueryControl.Edit_Smooth.Text := inttostr(round((xysorted^.x[high(xysorted^.x)]-xysorted^.x[0])/high(xysorted^.x))+1);
            ipval := round((xysorted^.x[high(xysorted^.x)]-xysorted^.x[0])/high(xysorted^.x))+1;
          end;

          newxy(xyip);
          //test :=  high(xysorted^.x);
          xyip := InterExtraPolate(xysorted,xysorted^.x[0],xysorted^.x[high(xysorted^.x)],ipval);

          SetLength(Curve[z]^.x,high(xyip^.x)+1);
          SetLength(Curve[z]^.y,high(xyip^.y)+1);
            //

          for j := 0 to high(xyip^.x) do begin
            Curve[z]^.x[j] := xyip^.x[j];
            Curve[z]^.y[j] := xyip^.y[j];

          end;

        end;
             * */
        }


        private static XYZ_Data FFT_Plaindata(double[] rawdata_y, XY_Data n100, double dt, double x0, CATimeStamps timestamps, Param_struct paramset)
        {
            //// Variable Definitions

            XYZ_Data aps_plaindata;     // Result of this function with single x, y and z axis. Dim(z) = Dim(x) * Dim(y)

            int               // Samples with data within one block
                lq_count,               // Number of main blocks
                ii1,                    // index of first valid leading quantity value (0 if leading quantity is not longer than rawdata)
                //ii2,                    // index of last valid leading quantity value (=leading_quantity.Length - 1 if leading quantity is not longer than rawdata)
                ii,                     // running index variable representing time
                jj,                     // running index variable representing frequency
                ColCnt;                 // number of calculated FFTs (=lq_count when no averaging is applied)
                                        // number of frequency points


            int[] pointerarray,         // for averaging: if overlapping leads to some blocks being identical, this array shows the index of the (singular) calculation result block, so pointerarray.Length > ColCnt
                time_indices,           // rawdata index of the center points of the calculation blocks
                N_nonzero,                   // 0 if leading quantity point is slightly later than the interpolated rawdata point, 1 if leading quantity point is earlier
                index;                  // sort index needed for the calculation of pointerarray

            double x_start,             // first query point if leaidng quantity is time
                x_end,                  // last query point if leaidng quantity is time
                lastx,
                average_factor;         // a factor often needed for average calculation


            double[] leading_quantity,  // leading quantity :-)
                querypoints,            // timing scheme of sub blocks for averaging with overlap. Example: average = 3, average_overlap = 50 => {-0.5, 0, 0.5}
                blockduration,
                leading_quantity_new,   // temporary leading quantity -> is copied to leading_quantity after some calculation steps
                sorted_singular,        // sorted leading quantity excluding identical subblocks (case averaging with overlapping) -> is copied to leading_quantity after some calculation steps
                apsvector,              // FFT result vector with all sub blocks
                windowfcn;              // window function (e.g. Hanning or rectangle)

            try
            {
                int N = Nextpow2(Math.Round(1 / dt) / paramset.delta_f);



                lastx = x0 + dt * (rawdata_y.Length - 1);
                double[] leading_quantity_temp;

                if (paramset.LQ == 2)
                {
                    leading_quantity_temp = Array.Empty<double>();
                    blockduration = new double[timestamps.t2.Length];
                    leading_quantity = new double[timestamps.t2.Length];
                    N_nonzero = new int[timestamps.t2.Length];
                    for (ii = 0; ii < blockduration.Length; ii++)
                    {
                        blockduration[ii] = timestamps.t2[ii] - timestamps.t1[ii];
                        leading_quantity[ii] = blockduration[ii] / 2.0 + timestamps.t1[ii];
                        N_nonzero[ii] = (int)(blockduration[ii] / dt);
                    }
                    paramset.dsrange1 = timestamps.t1[0];
                    paramset.dsrange2 = timestamps.t2[timestamps.t2.Length - 1];

                }
                else
                {



                    if (paramset.LQ == 1)
                    {

                        n100 = ProcessLQ_mono(n100);
                        leading_quantity_temp = new double[(int)((n100.ydata[n100.ydata.Length - 1] - paramset.delta_t * (Math.Ceiling(n100.ydata[0] / paramset.delta_t))) / paramset.delta_t) + 1];
                        for (ii = 0; ii < leading_quantity_temp.Length; ii++) leading_quantity_temp[ii] = paramset.delta_t * (Math.Ceiling(n100.ydata[0] / paramset.delta_t)) + ii * paramset.delta_t;
                        XY_Data lq_IP = Interp1(n100.ydata, n100.xdata, leading_quantity_temp);
                        leading_quantity = lq_IP.ydata;
                    }
                    else
                    {
                        leading_quantity_temp = Array.Empty<double>();
                        x_start = (int)((x0 + N / 2 * dt) / paramset.delta_t) * paramset.delta_t;
                        lq_count = (int)((lastx - x_start) / paramset.delta_t) + 1;
                        x_end = (lq_count - 1) * paramset.delta_t + x_start;
                        leading_quantity = new double[lq_count];
                        for (ii = 0; ii < lq_count; ii++) leading_quantity[ii] = ii * paramset.delta_t + x_start;
                    };
                    N_nonzero = new int[leading_quantity.Length];
                    for (ii = 0; ii < leading_quantity.Length; ii++) N_nonzero[ii] = N;

                }



                // Define leading quantity values - this has always the same unit(i.e.seconds) as the raw data!
                // external leading quantity is always read monotonically


                // Cut leading quantity if it is longer than rawdata or if it is beyond dsrange

                int dsrange1_index = (int)((paramset.dsrange1 - x0) / dt);
                int dsrange2_index = (int)((paramset.dsrange2 - x0) / dt);
                int[] ti = new int[leading_quantity.Length];
                for (ii = 0; ii < ti.Length; ii++) ti[ii] = (int)((leading_quantity[ii] - x0) / dt);

                average_factor = (paramset.average - 1) * (1 - paramset.average_overlap / 100);
                var nanref = new double[ti.Length];
                for (ii = 0; ii < ti.Length; ii++)
                {
                    if (ti[ii] - ((1 + average_factor) * Math.Floor((double)N_nonzero[ii] / 2)) < dsrange1_index || ti[ii] + ((1 + average_factor) * Math.Floor((double)N_nonzero[ii] / 2)) > dsrange2_index)
                    {
                        nanref[ii] = double.NaN;
                    }
                }
                RemoveNaN(ti, nanref);

                /*Array.Sort(ti);

                ii1 = 0;
                

                while (ti[ii1] - ((1 + average_factor) * Math.Floor((double)N_nonzero[ii1] / 2)) < dsrange1_index)
                {
                    ii1++;
                    if (ii1 > leading_quantity.Length - 1) break;
                }
                ii2 = leading_quantity.Length - 1;


                while (ti[ii2] + ((1 + average_factor) * Math.Floor((double)N_nonzero[ii2] / 2)) > dsrange2_index)
                {
                    ii2--;
                    if (ii2 < 0) break;
                }
                */
                //if (ii2 < ii1)
                if (ti.Length < 1)
                {
                    leading_quantity = new double[1];
                    leading_quantity[0] = (paramset.dsrange2 - paramset.dsrange1) / 2 + paramset.dsrange1;
                    leading_quantity_temp = new double[1];
                    leading_quantity_temp[0] = n100.ydata[0];
                    double diff;
                    if (n100.xdata[0] < leading_quantity[0])
                    {
                        for (ii = 1; ii < n100.xdata.Length; ii++)
                        {
                            diff = Math.Abs(n100.xdata[ii - 1] - leading_quantity[0]);
                            if (n100.xdata[ii] >= leading_quantity[0])
                            {
                                if (Math.Abs(n100.xdata[ii] - leading_quantity[0]) <= diff) leading_quantity_temp[0] = n100.ydata[ii];
                                else leading_quantity_temp[0] = n100.ydata[ii - 1];
                                break;
                            }
                        }
                    }

                    ti = new int[1];
                    ti[0] = (dsrange2_index - dsrange1_index) / 2 + dsrange1_index;


                    if (ti[0] - ((1 + average_factor) * N_nonzero[0] / 2) < dsrange1_index || ti[0] + ((1 + average_factor) * N_nonzero[0] / 2) > dsrange2_index)
                    {
                        aps_plaindata = new XYZ_Data
                        {
                            errorstring = "Not enough data for calculation"
                        };
                        return aps_plaindata;
                    }
                    else
                    {
                        //ii1 = 0;
                        //ii2 = 0;

                    }
                }

                /*leading_quantity_new = new double[(ii2 - ii1) + 1];
                for (ii = ii1; ii <= ii2; ii++) leading_quantity_new[ii - ii1] = leading_quantity[ii];
                leading_quantity = leading_quantity_new;*/
                leading_quantity = RemoveNaN(leading_quantity, nanref);
                N_nonzero = RemoveNaN(N_nonzero, nanref);
                for (ii = 0; ii < N_nonzero.Length; ii++) if (N_nonzero[ii] > N) N_nonzero[ii] = N;

                // Create output structure
                aps_plaindata = new XYZ_Data();

               /* var N_nonzero_new = new int[(leading_quantity.Length)];
                for (ii = 0; ii < N_nonzero_new.Length; ii++)
                {
                    if (N_nonzero[ii + ii1] > N) N_nonzero_new[ii] = N;
                    else N_nonzero_new[ii] = N_nonzero[ii + ii1];
                }
                N_nonzero = N_nonzero_new;*/



                // Write plain x-Data


                if (paramset.LQ == 2)
                {
                    //aps_plaindata.xdata = new double[(ii2 - ii1) + 1];
                    //for (ii = ii1; ii <= ii2; ii++) aps_plaindata.xdata[ii - ii1] = timestamps.cycles[ii];
                    aps_plaindata.xdata = RemoveNaN(timestamps.cycles, nanref);


                }
                else if (paramset.LQ == 1)
                {
                    //aps_plaindata.xdata = new double[(ii2 - ii1) + 1];
                    //for (ii = ii1; ii <= ii2; ii++) aps_plaindata.xdata[ii - ii1] = leading_quantity_temp[ii];
                    aps_plaindata.xdata = RemoveNaN(leading_quantity_temp, nanref);
                }
                else
                {
                    aps_plaindata.xdata = leading_quantity;
                }
                aps_plaindata.xtime = leading_quantity;



                // Calculate query points for averaging blocks avoiding double calculation of identical blocks
                if (paramset.average > 1)
                {
                    querypoints = new double[(int)(1 + average_factor / (1 - paramset.average_overlap / 100))];
                    for (ii = 0; ii < querypoints.Length; ii++) querypoints[ii] = ii * (1 - paramset.average_overlap / 100) - 0.5 * average_factor;


                    leading_quantity_new = new double[paramset.average * leading_quantity.Length];
                    for (ii = 0; ii < leading_quantity.Length; ii++) for (int jk = 0; jk < paramset.average; jk++)
                        {

                            leading_quantity_new[(ii * paramset.average + jk)] = leading_quantity[ii] + querypoints[jk] * N_nonzero[ii] * dt;
                        };

                    index = new int[leading_quantity_new.Length];
                    for (ii = 0; ii < index.Length; ii++) index[ii] = ii;

                    Array.Sort(leading_quantity_new, index);


                    pointerarray = new int[leading_quantity_new.Length];
                    sorted_singular = new double[leading_quantity_new.Length];
                    jj = 0;
                    pointerarray[0] = 0;
                    sorted_singular[0] = leading_quantity_new[0];
                    for (ii = 1; ii < leading_quantity_new.Length; ii++)
                    {
                        if (leading_quantity_new[ii] > leading_quantity_new[ii - 1])
                        {
                            jj++;
                            sorted_singular[jj] = leading_quantity_new[ii];
                        }
                        pointerarray[index[ii]] = jj;
                    }
                    leading_quantity = new double[jj + 1];
                    for (ii = 0; ii < leading_quantity.Length; ii++) leading_quantity[ii] = sorted_singular[ii];
                }
                else pointerarray = new int[2];
                ColCnt = leading_quantity.Length;

                if (paramset.LQ != 2)
                {
                    int N_nz = N_nonzero[0];
                    N_nonzero = new int[ColCnt];
                    for (ii = 0; ii < ColCnt; ii++) N_nonzero[ii] = N_nz;
                }


                // Write plain y-data



                // Calculate Indices of leading quantity query points (center of FFT block) for easy access
                time_indices = new int[ColCnt];

                for (ii = 0; ii < ColCnt; ii++)
                {
                    time_indices[ii] = (int)((leading_quantity[ii] - x0) / dt + 1) + 1;
                }


                // Calculate Window Function
                if (paramset.LQ == 2)
                {
                    int Nsum = 0;
                    for (ii = 0; ii < N_nonzero.Length; ii++) Nsum += N_nonzero[ii];
                    windowfcn = new double[Nsum];
                    ii1 = 0;
                    if (paramset.windowtype == 1)
                    {
                        for (ii = 0; ii < ColCnt; ii++)
                        {
                            for (jj = ii1; jj < ii1 + N_nonzero[ii]; jj++) windowfcn[jj] = 0.5 * (1 - Math.Cos(2 * Math.PI * (jj - ii1) / (N_nonzero[ii] - 1)));
                            ii1 += N_nonzero[ii];
                        }
                    }
                    else for (ii = 0; ii < windowfcn.Length; ii++) windowfcn[ii] = 1;

                }
                else
                {
                    windowfcn = new double[N_nonzero[0]];
                    if (paramset.windowtype == 1) for (ii = 0; ii < N_nonzero[0]; ii++) windowfcn[ii] = 0.5 * (1 - Math.Cos(2 * Math.PI * ii / (double)(N_nonzero[0] - 1))); //RICHTIG
                    //if (paramset.windowtype == 1) for (ii = 0; ii < N_nonzero[0]; ii++) windowfcn[ii] = 0.5 * (1 - Math.Cos(2 * Math.PI * ii / (double)(N_nonzero[0])));    // FALSCH (wie in Delphi)
                    else for (ii = 0; ii < N_nonzero[0]; ii++) windowfcn[ii] = 1;
                }
                int rows = (int)(N / 2.56) + 1;
                aps_plaindata.ydata = new double[rows];
                for (ii = 0; ii < aps_plaindata.ydata.Length; ii++) aps_plaindata.ydata[ii] = ii / dt / N;

                apsvector = new double[ColCnt * rows];

                int[] ColStart;

                if (paramset.LQ == 2) ColStart = Cumsum(N_nonzero);
                else ColStart = new int[1] { 0 };



                Parallel.For(0, ColCnt, i =>

                {
                    int j;
                    LomontFFT Lomontclass = new LomontFFT()
                    {
                        A = 1
                    };
                    double[] one_spectrum;
                    one_spectrum = new double[N];
                    if (paramset.LQ == 2) for (j = 1; j <= N_nonzero[i]; j++) one_spectrum[j - 1] = rawdata_y[time_indices[i] - N_nonzero[i] / 2 + j - 2] * paramset.y_unitconversion * windowfcn[ColStart[i] + j - 1];
                    else for (j = 1; j <= N_nonzero[i]; j++) one_spectrum[j - 1] = rawdata_y[time_indices[i] - N_nonzero[i] / 2 + j - 2] * paramset.y_unitconversion * windowfcn[j - 1];

                    Lomontclass.RealFFT(one_spectrum, true);
                    apsvector[i * rows] = Math.Pow(one_spectrum[0], 2);
                    for (j = 1; j < rows; j++) apsvector[i * rows + j] = Math.Pow(one_spectrum[(j - 1) * 2 + 2], 2) + Math.Pow(one_spectrum[(j - 1) * 2 + 3], 2);
                }
                );

                // Average Sub-Blocks to main blocks

                if (paramset.average > 1)
                {
                    aps_plaindata.zdata = new double[rows * pointerarray.Length / paramset.average];

                    Parallel.For(0, pointerarray.Length / paramset.average, i =>

                    {
                        int j;
                        for (j = 0; j < paramset.average; j++)
                        {
                            int col = pointerarray[i * paramset.average + j];
                            for (int kk = 0; kk < rows; kk++) aps_plaindata.zdata[i * rows + kk] += apsvector[col * rows + kk];
                        }
                        for (int kk = 0; kk < rows; kk++) aps_plaindata.zdata[i * rows + kk] /= paramset.average;
                    }
                    );


                }
                else aps_plaindata.zdata = apsvector;
                aps_plaindata.N_nonzero = N_nonzero;
                aps_plaindata.rows = rows;
                aps_plaindata.errorstring = null;
                // End of aps_plaindata function... Result can be used for any diagram type
                return aps_plaindata;
            }
            catch (IndexOutOfRangeException e)
            {

                aps_plaindata = new XYZ_Data
                {

                    errorstring = e.ToString()

                };
                return aps_plaindata;
            }





        }


        private static XYZ_Data FFT_Plaindata_VOICE(double[] rawdata_y, XY_Data n100, double dt, double x0, CATimeStamps timestamps, Param_struct paramset)
        {
            //// Variable Definitions

            XYZ_Data aps_plaindata;     // Result of this function with single x, y and z axis. Dim(z) = Dim(x) * Dim(y)

            int               // Samples with data within one block
                lq_count,               // Number of main blocks
                ii1,                    // index of first valid leading quantity value (0 if leading quantity is not longer than rawdata)
                ii2,                    // index of last valid leading quantity value (=leading_quantity.Length - 1 if leading quantity is not longer than rawdata)
                ii,                     // running index variable representing time
                jj,                     // running index variable representing frequency
                ColCnt;                 // number of calculated FFTs (=lq_count when no averaging is applied)
                                       // number of frequency points
      

            int[] pointerarray,         // for averaging: if overlapping leads to some blocks being identical, this array shows the index of the (singular) calculation result block, so pointerarray.Length > ColCnt
                time_indices,           // rawdata index of the center points of the calculation blocks
                N_nonzero,                   // 0 if leading quantity point is slightly later than the interpolated rawdata point, 1 if leading quantity point is earlier
                index;                  // sort index needed for the calculation of pointerarray
            
            double x_start,             // first query point if leaidng quantity is time
                x_end,                  // last query point if leaidng quantity is time
                lastx,
                average_factor;         // a factor often needed for average calculation


            double[] leading_quantity,  // leading quantity :-)
                querypoints,            // timing scheme of sub blocks for averaging with overlap. Example: average = 3, average_overlap = 50 => {-0.5, 0, 0.5}
                blockduration,
                leading_quantity_new,   // temporary leading quantity -> is copied to leading_quantity after some calculation steps
                sorted_singular,        // sorted leading quantity excluding identical subblocks (case averaging with overlapping) -> is copied to leading_quantity after some calculation steps
                apsvector,              // FFT result vector with all sub blocks
                windowfcn;              // window function (e.g. Hanning or rectangle)
            
            try
            {
                int N = Nextpow2(Math.Round(1/dt) / paramset.delta_f);
                

                
                lastx = x0 + dt * (rawdata_y.Length - 1);
                double[] leading_quantity_temp;

                if (paramset.LQ == 2)
                {
                    leading_quantity_temp = Array.Empty<double>();
                    blockduration = new double[timestamps.t2.Length];
                    leading_quantity = new double[timestamps.t2.Length];
                    N_nonzero = new int[timestamps.t2.Length];
                    for (ii = 0; ii < blockduration.Length; ii++)
                    {
                        blockduration[ii] = timestamps.t2[ii] - timestamps.t1[ii];
                        leading_quantity[ii] = blockduration[ii] / 2.0 + timestamps.t1[ii];
                        N_nonzero[ii] = (int)(blockduration[ii] / dt);
                    }
                    paramset.dsrange1 = timestamps.t1[0];
                    paramset.dsrange2 = timestamps.t2[timestamps.t2.Length - 1];

                }
                else
                {
                    

                    
                    if (paramset.LQ == 1)
                    {

                        n100 = ProcessLQ_mono(n100);
                        leading_quantity_temp = new double[(int)((n100.ydata[n100.ydata.Length - 1] - paramset.delta_t * (Math.Ceiling(n100.ydata[0] / paramset.delta_t))) / paramset.delta_t) + 1];
                        for (ii = 0; ii < leading_quantity_temp.Length; ii++) leading_quantity_temp[ii] = paramset.delta_t * (Math.Ceiling(n100.ydata[0] / paramset.delta_t)) + ii * paramset.delta_t;
                        XY_Data lq_IP = Interp1(n100.ydata, n100.xdata, leading_quantity_temp);
                        leading_quantity = lq_IP.ydata;
                    }
                    else
                    {
                        leading_quantity_temp = Array.Empty<double>();
                        x_start = (int)((x0 + N / 2 * dt) / paramset.delta_t) * paramset.delta_t;
                        lq_count = (int)((lastx - x_start) / paramset.delta_t) + 1;
                        x_end = (lq_count - 1) * paramset.delta_t + x_start;
                        leading_quantity = new double[lq_count];
                        for (ii = 0; ii < lq_count; ii++) leading_quantity[ii] = ii * paramset.delta_t + x_start;
                    };
                    N_nonzero = new int[leading_quantity.Length];
                    for (ii = 0; ii < leading_quantity.Length; ii++) N_nonzero[ii] = N;

                }


                
                // Define leading quantity values - this has always the same unit(i.e.seconds) as the raw data!
                // external leading quantity is always read monotonically
                

                // Cut leading quantity if it is longer than rawdata or if it is beyond dsrange

                int dsrange1_index = (int)((paramset.dsrange1 - x0) / dt);
                int dsrange2_index = (int)((paramset.dsrange2 - x0) / dt);
                int[] ti = new int[leading_quantity.Length];
                for (ii = 0; ii < ti.Length; ii++) ti[ii] = (int)((leading_quantity[ii] - x0) / dt);

                ii1 = 0;
                average_factor = (paramset.average - 1) * (1 - paramset.average_overlap / 100);

                while (ti[ii1] - ((1 + average_factor) * Math.Floor((double)N_nonzero[ii1] / 2)) < dsrange1_index)
                {
                    ii1++;
                    if (ii1 > leading_quantity.Length - 1) break;
                }
                ii2 = leading_quantity.Length - 1;

                
                while (ti[ii2] + ((1 + average_factor) * Math.Floor((double)N_nonzero[ii2] / 2)) > dsrange2_index)
                {
                    ii2--;
                    if (ii2 < 0) break;
                }

                if (ii2 < ii1)
                {
                    leading_quantity = new double[1];
                    leading_quantity[0] = (paramset.dsrange2 - paramset.dsrange1) / 2 + paramset.dsrange1;
                    leading_quantity_temp = new double[1];
                    leading_quantity_temp[0] = n100.ydata[0];
                    double diff;
                    if (n100.xdata[0] < leading_quantity[0])
                    {
                        for (ii = 1; ii < n100.xdata.Length; ii++)
                        {
                            diff = Math.Abs(n100.xdata[ii - 1] - leading_quantity[0]);
                            if (n100.xdata[ii] >= leading_quantity[0])
                            {
                                if (Math.Abs(n100.xdata[ii] - leading_quantity[0]) <= diff) leading_quantity_temp[0] = n100.ydata[ii];
                                else leading_quantity_temp[0] = n100.ydata[ii - 1];
                                break;
                            }
                        }
                    }
                    
                    ti = new int[1];
                    ti[0] = (dsrange2_index - dsrange1_index) / 2 + dsrange1_index;

                    
                    if (ti[0] - ((1 + average_factor) * N_nonzero[0] / 2) < dsrange1_index || ti[0] + ((1 + average_factor) * N_nonzero[0] / 2) > dsrange2_index)
                    {
                        aps_plaindata = new XYZ_Data
                        {
                            errorstring = "Not enough data for calculation"
                        };
                        return aps_plaindata;
                    }
                    else
                    {
                        ii1 = 0;
                        ii2 = 0;
                    }
                }
                
                leading_quantity_new = new double[(ii2 - ii1) + 1];
                for (ii = ii1; ii <= ii2; ii++) leading_quantity_new[ii - ii1] = leading_quantity[ii];
                leading_quantity = leading_quantity_new;

                // Create output structure
                aps_plaindata = new XYZ_Data();

                var N_nonzero_new = new int[(leading_quantity.Length)];
                for (ii = 0; ii < N_nonzero_new.Length; ii++)
                {
                    if (N_nonzero[ii + ii1] > N) N_nonzero_new[ii] = N;
                    else N_nonzero_new[ii] = N_nonzero[ii + ii1];
                }
                N_nonzero = N_nonzero_new;



                // Write plain x-Data
                

                if (paramset.LQ == 2)
                {
                    aps_plaindata.xdata = new double[(ii2 - ii1) + 1];
                    for (ii = ii1; ii <= ii2; ii++) aps_plaindata.xdata[ii - ii1] = timestamps.cycles[ii];
                    

                }
                else if (paramset.LQ == 1)
                {
                    aps_plaindata.xdata = new double[(ii2 - ii1) + 1];
                    for (ii = ii1; ii <= ii2; ii++) aps_plaindata.xdata[ii - ii1] = leading_quantity_temp[ii];
                }
                else
                {
                    aps_plaindata.xdata = leading_quantity;
                }
                aps_plaindata.xtime = leading_quantity;

                

                // Calculate query points for averaging blocks avoiding double calculation of identical blocks
                if (paramset.average > 1)
                {
                    querypoints = new double[(int)(1 + average_factor / (1 - paramset.average_overlap / 100))];
                    for (ii = 0; ii < querypoints.Length; ii++) querypoints[ii] = ii * (1 - paramset.average_overlap / 100) - 0.5 * average_factor;

                  
                    leading_quantity_new = new double[paramset.average * leading_quantity.Length];
                    for (ii = 0; ii < leading_quantity.Length; ii++) for (int jk = 0; jk < paramset.average; jk++)
                        {
                  
                            leading_quantity_new[(ii * paramset.average + jk)] = leading_quantity[ii] + querypoints[jk] * N_nonzero[ii] * dt;
                        };

                    index = new int[leading_quantity_new.Length];
                    for (ii = 0; ii < index.Length; ii++) index[ii] = ii;

                    Array.Sort(leading_quantity_new, index);
                  

                    pointerarray = new int[leading_quantity_new.Length];
                    sorted_singular = new double[leading_quantity_new.Length];
                    jj = 0;
                    pointerarray[0] = 0;
                    sorted_singular[0] = leading_quantity_new[0];
                    for (ii = 1; ii < leading_quantity_new.Length; ii++)
                    {
                        if (leading_quantity_new[ii] > leading_quantity_new[ii - 1])
                        {
                            jj++;
                            sorted_singular[jj] = leading_quantity_new[ii];
                        }
                        pointerarray[index[ii]] = jj;
                    }
                    leading_quantity = new double[jj + 1];
                    for (ii = 0; ii < leading_quantity.Length; ii++) leading_quantity[ii] = sorted_singular[ii];
                }
                else pointerarray = new int[2];
                ColCnt = leading_quantity.Length;

                if (paramset.LQ != 2)
                {
                    int N_nz = N_nonzero[0];
                    N_nonzero = new int[ColCnt];
                    for (ii = 0; ii < ColCnt; ii++) N_nonzero[ii] = N_nz;
                }

                   
                // Write plain y-data
                


                // Calculate Indices of leading quantity query points (center of FFT block) for easy access
                time_indices = new int[ColCnt];
                
                for (ii = 0; ii < ColCnt; ii++)
                {
                    time_indices[ii] = (int)((leading_quantity[ii] - x0) / dt + 1) + 1;
                }
                

                // Calculate Window Function
                if (paramset.LQ == 2)
                {
                    int Nsum = 0;
                    for (ii = 0; ii < N_nonzero.Length; ii++) Nsum += N_nonzero[ii];
                    windowfcn = new double[Nsum];
                    ii1 = 0;
                    if (paramset.windowtype == 1)
                    {
                        for (ii = 0; ii < ColCnt; ii++)
                        {
                            for (jj = ii1; jj < ii1 + N_nonzero[ii]; jj++) windowfcn[jj] = 0.5 * (1 - Math.Cos(2 * Math.PI * (jj - ii1) / (N_nonzero[ii] - 1)));
                            ii1 += N_nonzero[ii];
                        }
                    }
                    else for (ii = 0; ii < windowfcn.Length; ii++) windowfcn[ii] = 1;

                }
                else
                {
                    windowfcn = new double[N_nonzero[0]];
                    if (paramset.windowtype == 1) for (ii = 0; ii < N_nonzero[0]; ii++) windowfcn[ii] = 0.5 * (1 - Math.Cos(2 * Math.PI * ii / (double)(N_nonzero[0] - 1))); //RICHTIG
                    //if (paramset.windowtype == 1) for (ii = 0; ii < N_nonzero[0]; ii++) windowfcn[ii] = 0.5 * (1 - Math.Cos(2 * Math.PI * ii / (double)(N_nonzero[0])));    // FALSCH (wie in Delphi)
                    else for (ii = 0; ii < N_nonzero[0]; ii++) windowfcn[ii] = 1;
                }
                int rows = (int)(N / 2);
                aps_plaindata.ydata = new double[rows];
                for (ii = 0; ii < aps_plaindata.ydata.Length; ii++) aps_plaindata.ydata[ii] = ii / dt / N;
                
                apsvector = new double[ColCnt * rows];

                int[] ColStart;

                if (paramset.LQ == 2) ColStart = Cumsum(N_nonzero);
                else ColStart = new int[1] { 0 };
                
                
                
                Parallel.For(0, ColCnt, i =>
                
                {
                    int j;
                    LomontFFT Lomontclass = new LomontFFT()
                    {
                        A = 1
                    };
                    double[] one_spectrum;
                    one_spectrum = new double[N];
                    if (paramset.LQ == 2) for (j = 1; j <= N_nonzero[i]; j++) one_spectrum[j - 1] = rawdata_y[time_indices[i] - N_nonzero[i] / 2 + j - 2] * paramset.y_unitconversion * windowfcn[ColStart[i] + j - 1];
                    else for (j = 1; j <= N_nonzero[i]; j++) one_spectrum[j - 1] = rawdata_y[time_indices[i] - N_nonzero[i] / 2 + j - 2] * paramset.y_unitconversion * windowfcn[j - 1];
                    
                    Lomontclass.RealFFT(one_spectrum, true);
                    apsvector[i * rows] = Math.Pow(one_spectrum[0], 2);
                    for (j = 1; j < rows; j++) apsvector[i * rows + j] = Math.Pow(one_spectrum[(j - 1) * 2 + 2], 2) + Math.Pow(one_spectrum[(j - 1) * 2 + 3], 2);
                }
                );
                
                // Average Sub-Blocks to main blocks

                if (paramset.average > 1)
                {
                    aps_plaindata.zdata = new double[rows * pointerarray.Length / paramset.average];
                    
                    Parallel.For(0, pointerarray.Length / paramset.average, i =>
                    
                    {
                        int j;
                        for (j = 0; j < paramset.average; j++)
                        {
                            int col = pointerarray[i * paramset.average + j];
                            for (int kk = 0; kk < rows; kk++) aps_plaindata.zdata[i * rows + kk] += apsvector[col * rows + kk];
                        }
                        for (int kk = 0; kk < rows; kk++) aps_plaindata.zdata[i * rows + kk] /= paramset.average;
                    }
                    );
                    

                }
                else aps_plaindata.zdata = apsvector;
                aps_plaindata.N_nonzero = N_nonzero;
                aps_plaindata.rows = rows;
                aps_plaindata.errorstring = null;
                // End of aps_plaindata function... Result can be used for any diagram type
                return aps_plaindata;
            }
            catch (IndexOutOfRangeException e)
            {
                
                aps_plaindata = new XYZ_Data
                {
                    
                    errorstring = e.ToString()
                    
                };
                return aps_plaindata;
            }



            

        }

        internal static int CheckDimension(int L, XY_Data n100, double dt, double x0, CATimeStamps timestamps, Param_struct paramset)
        {
            int dimensions;

           // XYZ_Data aps_plaindata;     // Result of this function with single x, y and z axis. Dim(z) = Dim(x) * Dim(y)

            int               // Samples with data within one block
                lq_count,               // Number of main blocks
                ii1,                    // index of first valid leading quantity value (0 if leading quantity is not longer than rawdata)
                ii2,                    // index of last valid leading quantity value (=leading_quantity.Length - 1 if leading quantity is not longer than rawdata)
                ii;                     // running index variable representing time
                


            int[] N_nonzero;                   // 0 if leading quantity point is slightly later than the interpolated rawdata point, 1 if leading quantity point is earlier
                

            double x_start,             // first query point if leaidng quantity is time
                //x_end,                  // last query point if leaidng quantity is time
                lastx,
                average_factor;         // a factor often needed for average calculation


            double[] leading_quantity,  // leading quantity :-)
                blockduration,
                leading_quantity_new;   // temporary leading quantity -> is copied to leading_quantity after some calculation steps


            int N = Nextpow2(1 / dt / paramset.delta_f);
            
            lastx = x0 + dt * (L - 1);
            double[] leading_quantity_temp;

            if (paramset.LQ == 2)
            {
                //leading_quantity_temp = new double[0];
                blockduration = new double[timestamps.t2.Length];
                leading_quantity = new double[timestamps.t2.Length];
                N_nonzero = new int[timestamps.t2.Length];
                for (ii = 0; ii < blockduration.Length; ii++)
                {
                    blockduration[ii] = timestamps.t2[ii] - timestamps.t1[ii];
                    leading_quantity[ii] = blockduration[ii] / 2.0 + timestamps.t1[ii];
                    N_nonzero[ii] = (int)(blockduration[ii] / dt);
                }
                paramset.dsrange1 = timestamps.t1[0];
                paramset.dsrange2 = timestamps.t2[timestamps.t2.Length - 1];

            }
            else
            {



                if (paramset.LQ == 1)
                {

                    n100 = ProcessLQ_mono(n100);
                    leading_quantity_temp = new double[(int)((n100.ydata[n100.ydata.Length - 1] - paramset.delta_t * (Math.Ceiling(n100.ydata[0] / paramset.delta_t))) / paramset.delta_t) + 1];
                    for (ii = 0; ii < leading_quantity_temp.Length; ii++) leading_quantity_temp[ii] = paramset.delta_t * (Math.Ceiling(n100.ydata[0] / paramset.delta_t)) + ii * paramset.delta_t;
                    XY_Data lq_IP = Interp1(n100.ydata, n100.xdata, leading_quantity_temp);
                    leading_quantity = lq_IP.ydata;
                }
                else
                {
                    //leading_quantity_temp = new double[0];
                    x_start = (int)((x0 + N / 2 * dt) / paramset.delta_t) * paramset.delta_t;
                    lq_count = (int)((lastx - x_start) / paramset.delta_t) + 1;
                    //x_end = (lq_count - 1) * paramset.delta_t + x_start;
                    leading_quantity = new double[lq_count];
                    for (ii = 0; ii < lq_count; ii++) leading_quantity[ii] = ii * paramset.delta_t + x_start;
                };
                N_nonzero = new int[leading_quantity.Length];
                for (ii = 0; ii < leading_quantity.Length; ii++) N_nonzero[ii] = N;

            }


            int dsrange1_index = (int)((paramset.dsrange1 - x0) / dt);
            int dsrange2_index = (int)((paramset.dsrange2 - x0) / dt);
            int[] ti = new int[leading_quantity.Length];
            for (ii = 0; ii < ti.Length; ii++) ti[ii] = (int)((leading_quantity[ii] - x0) / dt);

            ii1 = 0;
            average_factor = (paramset.average - 1) * (1 - paramset.average_overlap / 100);
            //while (paramset.dsrange1 - (leading_quantity[ii1] - ((1 + average_factor) * N_nonzero[ii1] / 2) * dt) > dt / 2)
            //while (paramset.dsrange1 - (leading_quantity[ii1] - ((1 + average_factor) * N_nonzero[ii1] / 2) * dt) > dt / 2)
            while (ti[ii1] - ((1 + average_factor) * Math.Floor((double)N_nonzero[ii1] / 2)) < dsrange1_index)
            {
                ii1++;
                if (ii1 > leading_quantity.Length - 1) break;
            }
            ii2 = leading_quantity.Length - 1;

            //while (leading_quantity[ii2] + ((1 + average_factor) * N_nonzero[ii2] / 2 - 1) * dt - paramset.dsrange2 > dt / 2)
            while (ti[ii2] + ((1 + average_factor) * Math.Floor((double)N_nonzero[ii2] / 2)) > dsrange2_index)
            {
                ii2--;
                if (ii2 < 0) break;
            }

            if (ii2 < ii1)
            {
                ti = new int[1];
                ti[0] = (dsrange2_index - dsrange1_index) / 2 + dsrange1_index;

                if (ti[0] - ((1 + average_factor) * N_nonzero[0] / 2) < dsrange1_index || ti[0] + ((1 + average_factor) * N_nonzero[0] / 2) > dsrange2_index) return 0;
                else
                {
                    //ii1 = 0;
                    //ii2 = 0;
                    if (paramset.diagramtype < 4) return 2;
                    else return 1;

                }
            }

            leading_quantity_new = new double[(ii2 - ii1) + 1];
            for (ii = ii1; ii <= ii2; ii++) leading_quantity_new[ii - ii1] = leading_quantity[ii];
            leading_quantity = leading_quantity_new;

            if (leading_quantity.Length == 1 || paramset.mean)
            {
                if (paramset.diagramtype < 4) dimensions = 2;
                else dimensions = 1;
            }
            else
            {
                if (paramset.diagramtype < 4) dimensions = 3;
                else dimensions = 2;
            }

            return dimensions;
        }

        internal static XYZ_Data CalcAPSData(double[] rawdata_y, XY_Data n100, double dt, double x0, CATimeStamps timestamps, Param_struct paramset)
        {
            //// Variable Definitions

            XYZ_Data aps_plaindata,         // X,Y,Z data from Plaindata Function
                result;                     // Result Data X, Y, Z

            int /*rows_in,                    // number of frequency points with DC
                rows_out,                   // number of frequency points with/without DC
                lq_count,                   // number of FFT blocks
                N_nonzero,*/                  // number of non zero samples within one FFT block (zero padding)
                wincorr,                    // window amplitude correction factor (e.g. 2 for Hann window)
                ii;                         // running variable for loops

            int[] N_nonzero;

            double rmsfactor, overalllevel;

            double[] freq_weight_curve,     // Frequency weights vector (factor values)
                aps_ydata,                  // output frequency data in concerto format (same length of each axis)
                aps_xdata,                  // output leading quantity data in concerto format (same length of each axis)
                aps_zdata;

            //double blockduration;

            

            // FFT Calculations
            aps_plaindata = FFT_Plaindata(rawdata_y, n100, dt, x0, timestamps, paramset);
           

            if (aps_plaindata.errorstring != null)
            {
                result = new XYZ_Data
                {
                    errorstring = aps_plaindata.errorstring
                };
                return result;
            }

            /*if (paramset.CABase)
            {
                blockduration = paramset.toCA - paramset.fromCA;
                N_nonzero = (int)Math.Round(blockduration / dt, MidpointRounding.AwayFromZero);
                cycfactor = 1.0 / 720.0;
            }
            else*/
            //{
            N_nonzero = aps_plaindata.N_nonzero;
            
            //}

            // Simple variable calculations
            /*rows_in = aps_plaindata.ydata.Length;
            rows_out = rows_in - 1 + paramset.DC;
            lq_count = aps_plaindata.xdata.Length;
            N_nonzero = paramset.N / paramset.zeropad;*/

            int time_blocks = aps_plaindata.xdata.Length;

            // Calculate Frequency Weights
            if (paramset.freq_weight == 0)
            {
                freq_weight_curve = new double[aps_plaindata.ydata.Length];
                for (ii = 0; ii < freq_weight_curve.Length; ii++) freq_weight_curve[ii] = 1;
            }
            else
            {
                freq_weight_curve = ABC_Weight(aps_plaindata.ydata, paramset.freq_weight);
                /*else
                {
                    freq_weight_curve = new double[aps_plaindata.ydata.Length];
                    for (ii = 0; ii < time_blocks; ii++)
                    {
                        double[] this_freqaxis = new double[aps_plaindata.rows];
                        for (int jj = 0; jj < aps_plaindata.rows; jj++) this_freqaxis[jj] = aps_plaindata.ydata[aps_plaindata.ColStart[ii] + jj];
                        double[] this_freq_weight_curve = ABC_Weight(this_freqaxis, paramset.freq_weight);
                        for (int jj = 0; jj < aps_plaindata.rows; jj++) freq_weight_curve[aps_plaindata.ColStart[ii] + jj] = this_freq_weight_curve[jj];
                    }
                }*/
            }

            

            if (paramset.windowtype == 1) wincorr = 2;
            else wincorr = 1;

            
            if (paramset.y_amplitude == 1)
                rmsfactor = Math.Sqrt(2);
            else if (paramset.y_amplitude == 2)
                rmsfactor = 2 * Math.Sqrt(2);
            else rmsfactor = 1;

            //int N = Nextpow2(Math.Round(1 / paramset.delta_f / dt));
            overalllevel = CalcOverallLevelSingleValue(aps_plaindata.zdata, aps_plaindata.ydata.Length, N_nonzero, freq_weight_curve, paramset.windowtype, aps_plaindata.xdata.Length, paramset.DC);
            overalllevel *= rmsfactor;
            if (paramset.y_axis == 1) overalllevel = 20 * Math.Log10(overalllevel / paramset.dBref);

            // Calculate output Data
            //int[] rows_out = new int[aps_plaindata.rows.Length];
            int rows_out;
            //int[] ColStart_out = new int[aps_plaindata.ColStart.Length];
            //ColStart_out = aps_plaindata.ColStart;
            rows_out = aps_plaindata.rows - 1 + paramset.DC;
            /*if (paramset.DC == 0)
            {
                rows_out = aps_plaindata.rows - 1;
                //for (ii = 0; ii < ColStart_out.Length; ii++) ColStart_out[ii] = aps_plaindata.ColStart[ii] - ii;
            }
            else
            {
                rows_out = aps_plaindata.rows;
                //ColStart_out = aps_plaindata.ColStart;
            }*/
            aps_zdata = new double[rows_out * time_blocks];
            aps_ydata = new double[aps_zdata.Length];
            aps_xdata = new double[aps_zdata.Length];

            Parallel.For(0, aps_plaindata.xdata.Length, i =>
            //for (int ii = 0; ii < time_blocks; ii++)
            {
                int j;
                //int k;
                /*if (paramset.LQ == 2)*/
                //k = i;
                //else k = 0;
                if (paramset.DC == 1)
                {
                    if (paramset.y_axis == 1) aps_zdata[i * rows_out] = rmsfactor * 20 * Math.Log10(wincorr * Math.Sqrt(aps_plaindata.zdata[i * aps_plaindata.rows] / 2.0) * freq_weight_curve[0] / N_nonzero[i] / paramset.dBref);
                    else aps_zdata[i * rows_out] = rmsfactor * wincorr * Math.Sqrt(aps_plaindata.zdata[i * aps_plaindata.rows] / 2.0) * freq_weight_curve[0] / N_nonzero[i];
                    aps_xdata[i * rows_out] = aps_plaindata.xdata[i];
                    aps_ydata[i * rows_out] = aps_plaindata.ydata[0];
                }
                if (paramset.y_axis == 1)
                {
                    for (j = 0 + paramset.DC; j < rows_out; j++) //paramset.DC = 0
                    {
                        aps_zdata[i * rows_out + j] = rmsfactor * 20 * Math.Log10(2 * wincorr * Math.Sqrt(aps_plaindata.zdata[i * aps_plaindata.rows + j + 1 - paramset.DC] / 2.0) * freq_weight_curve[j + 1 - paramset.DC] / N_nonzero[i] / paramset.dBref);
                        aps_xdata[i * rows_out + j] = aps_plaindata.xdata[i];
                        aps_ydata[i * rows_out + j] = aps_plaindata.ydata[j + 1 - paramset.DC];
                    }
                }
                else
                {
                    for (j = 0 + paramset.DC; j < rows_out; j++) //paramset.DC = 0
                    {
                        aps_zdata[i * rows_out + j] = rmsfactor * 2 * wincorr * Math.Sqrt(aps_plaindata.zdata[i * aps_plaindata.rows + j + 1 - paramset.DC] / 2.0) * freq_weight_curve[j + 1 - paramset.DC] / N_nonzero[i];
                        aps_xdata[i * rows_out + j] = aps_plaindata.xdata[i];
                        aps_ydata[i * rows_out + j] = aps_plaindata.ydata[j + 1 - paramset.DC];
                    }
                }
                    
            }
                //;
            );

            //if (paramset.CABase) result.dimensions = new int[2] { aps_plaindata.xdata.Length, -1 };
            /*else */
            result.dimensions = new int[2] { aps_plaindata.xdata.Length, rows_out };
            

            result.zdata = aps_zdata;
            result.ydata = aps_ydata;
            result.xdata = aps_xdata;
            result.xtime = aps_plaindata.xtime;
            result.freq_labels = new double[1];
            result.freq_labels[0] = 0;
            result.overalllevel = overalllevel;
            result.rows = rows_out;
            result.N_nonzero = aps_plaindata.N_nonzero;
            result.errorstring = null;
            result.paramset = paramset;
            result.xdescription = "";
            //result.ColStart = ColStart_out;

            if (paramset.mean)
            {
                result = AverageDiagramOverTime(result, aps_plaindata.xdata);
            }


            return result;
        }


        internal static XYZ_Data CalcAPSData_VOICE(double[] rawdata_y, XY_Data n100, double dt, double x0, CATimeStamps timestamps, Param_struct paramset)
        {
            //// Variable Definitions

            XYZ_Data aps_plaindata,         // X,Y,Z data from Plaindata Function
                result;                     // Result Data X, Y, Z

            int /*rows_in,                    // number of frequency points with DC
                rows_out,                   // number of frequency points with/without DC
                lq_count,                   // number of FFT blocks
                N_nonzero,*/                  // number of non zero samples within one FFT block (zero padding)
                wincorr,                    // window amplitude correction factor (e.g. 2 for Hann window)
                ii;                         // running variable for loops

            int[] N_nonzero;

            double rmsfactor, overalllevel;

            double[] freq_weight_curve,     // Frequency weights vector (factor values)
                aps_ydata,                  // output frequency data in concerto format (same length of each axis)
                aps_xdata,                  // output leading quantity data in concerto format (same length of each axis)
                aps_zdata;

            //double blockduration;



            // FFT Calculations
            aps_plaindata = FFT_Plaindata_VOICE(rawdata_y, n100, dt, x0, timestamps, paramset);


            if (aps_plaindata.errorstring != null)
            {
                result = new XYZ_Data
                {
                    errorstring = aps_plaindata.errorstring
                };
                return result;
            }

            /*if (paramset.CABase)
            {
                blockduration = paramset.toCA - paramset.fromCA;
                N_nonzero = (int)Math.Round(blockduration / dt, MidpointRounding.AwayFromZero);
                cycfactor = 1.0 / 720.0;
            }
            else*/
            //{
            N_nonzero = aps_plaindata.N_nonzero;

            //}

            // Simple variable calculations
            /*rows_in = aps_plaindata.ydata.Length;
            rows_out = rows_in - 1 + paramset.DC;
            lq_count = aps_plaindata.xdata.Length;
            N_nonzero = paramset.N / paramset.zeropad;*/

            int time_blocks = aps_plaindata.xdata.Length;

            // Calculate Frequency Weights
            if (paramset.freq_weight == 0)
            {
                freq_weight_curve = new double[aps_plaindata.ydata.Length];
                for (ii = 0; ii < freq_weight_curve.Length; ii++) freq_weight_curve[ii] = 1;
            }
            else
            {
                freq_weight_curve = ABC_Weight(aps_plaindata.ydata, paramset.freq_weight);
                /*else
                {
                    freq_weight_curve = new double[aps_plaindata.ydata.Length];
                    for (ii = 0; ii < time_blocks; ii++)
                    {
                        double[] this_freqaxis = new double[aps_plaindata.rows];
                        for (int jj = 0; jj < aps_plaindata.rows; jj++) this_freqaxis[jj] = aps_plaindata.ydata[aps_plaindata.ColStart[ii] + jj];
                        double[] this_freq_weight_curve = ABC_Weight(this_freqaxis, paramset.freq_weight);
                        for (int jj = 0; jj < aps_plaindata.rows; jj++) freq_weight_curve[aps_plaindata.ColStart[ii] + jj] = this_freq_weight_curve[jj];
                    }
                }*/
            }



            if (paramset.windowtype == 1) wincorr = 2;
            else wincorr = 1;


            if (paramset.y_amplitude == 1)
                rmsfactor = Math.Sqrt(2);
            else if (paramset.y_amplitude == 2)
                rmsfactor = 2 * Math.Sqrt(2);
            else rmsfactor = 1;

            //int N = Nextpow2(Math.Round(1 / paramset.delta_f / dt));
            overalllevel = CalcOverallLevelSingleValue(aps_plaindata.zdata, aps_plaindata.ydata.Length, N_nonzero, freq_weight_curve, paramset.windowtype, aps_plaindata.xdata.Length, paramset.DC);
            overalllevel *= rmsfactor;
            if (paramset.y_axis == 1) overalllevel = 20 * Math.Log10(overalllevel / paramset.dBref);

            // Calculate output Data
            //int[] rows_out = new int[aps_plaindata.rows.Length];
            int rows_out;
            //int[] ColStart_out = new int[aps_plaindata.ColStart.Length];
            //ColStart_out = aps_plaindata.ColStart;
            rows_out = aps_plaindata.rows - 1 + paramset.DC;
            /*if (paramset.DC == 0)
            {
                rows_out = aps_plaindata.rows - 1;
                //for (ii = 0; ii < ColStart_out.Length; ii++) ColStart_out[ii] = aps_plaindata.ColStart[ii] - ii;
            }
            else
            {
                rows_out = aps_plaindata.rows;
                //ColStart_out = aps_plaindata.ColStart;
            }*/
            aps_zdata = new double[rows_out * time_blocks];
            aps_ydata = new double[aps_zdata.Length];
            aps_xdata = new double[aps_zdata.Length];

            Parallel.For(0, aps_plaindata.xdata.Length, i =>
            //for (int ii = 0; ii < time_blocks; ii++)
            {
                int j;
                //int k;
                /*if (paramset.LQ == 2)*/
                //k = i;
                //else k = 0;
                if (paramset.DC == 1)
                {
                    if (paramset.y_axis == 1) aps_zdata[i * rows_out] = rmsfactor * 20 * Math.Log10(wincorr * Math.Sqrt(aps_plaindata.zdata[i * aps_plaindata.rows] / 2.0) * freq_weight_curve[0] / N_nonzero[i] / paramset.dBref);
                    else aps_zdata[i * rows_out] = rmsfactor * wincorr * Math.Sqrt(aps_plaindata.zdata[i * aps_plaindata.rows] / 2.0) * freq_weight_curve[0] / N_nonzero[i];
                    aps_xdata[i * rows_out] = aps_plaindata.xdata[i];
                    aps_ydata[i * rows_out] = aps_plaindata.ydata[0];
                }
                if (paramset.y_axis == 1)
                {
                    for (j = 0 + paramset.DC; j < rows_out; j++) //paramset.DC = 0
                    {
                        aps_zdata[i * rows_out + j] = rmsfactor * 20 * Math.Log10(2 * wincorr * Math.Sqrt(aps_plaindata.zdata[i * aps_plaindata.rows + j + 1 - paramset.DC] / 2.0) * freq_weight_curve[j + 1 - paramset.DC] / N_nonzero[i] / paramset.dBref);
                        aps_xdata[i * rows_out + j] = aps_plaindata.xdata[i];
                        aps_ydata[i * rows_out + j] = aps_plaindata.ydata[j + 1 - paramset.DC];
                    }
                }
                else
                {
                    for (j = 0 + paramset.DC; j < rows_out; j++) //paramset.DC = 0
                    {
                        aps_zdata[i * rows_out + j] = rmsfactor * 2 * wincorr * Math.Sqrt(aps_plaindata.zdata[i * aps_plaindata.rows + j + 1 - paramset.DC] / 2.0) * freq_weight_curve[j + 1 - paramset.DC] / N_nonzero[i];
                        aps_xdata[i * rows_out + j] = aps_plaindata.xdata[i];
                        aps_ydata[i * rows_out + j] = aps_plaindata.ydata[j + 1 - paramset.DC];
                    }
                }

            }
            //;
            );

            //if (paramset.CABase) result.dimensions = new int[2] { aps_plaindata.xdata.Length, -1 };
            /*else */
            result.dimensions = new int[2] { aps_plaindata.xdata.Length, rows_out };


            result.zdata = aps_zdata;
            result.ydata = aps_ydata;
            result.xdata = aps_xdata;
            result.xtime = aps_plaindata.xtime;
            result.freq_labels = new double[1];
            result.freq_labels[0] = 0;
            result.overalllevel = overalllevel;
            result.rows = rows_out;
            result.N_nonzero = aps_plaindata.N_nonzero;
            result.errorstring = null;
            result.paramset = paramset;
            result.xdescription = "";
            //result.ColStart = ColStart_out;

            if (paramset.mean)
            {
                result = AverageDiagramOverTime(result, aps_plaindata.xdata);
            }


            return result;
        }

        private static double CalcOverallLevelSingleValue(double[] aps_plaindata_zdata, int rows_in, int[] N_nonzero, double[] freq_weight_curve, int windowtype, int time_blocks, int DC)
        {
            double[] aps_data_zdata = new double[time_blocks];
            double energycorr, OL;
            if (windowtype == 1) energycorr = Math.Sqrt(8 / 3);
            else energycorr = 1;

            Parallel.For(0, time_blocks, i =>
            //for (int i = 0; i < time_blocks; i++) 
            {
                double summe = 0;
                summe += 0.5 * DC * aps_plaindata_zdata[i * rows_in] * freq_weight_curve[0];
                for (int j = 1; j < rows_in; j++) summe += 2 * aps_plaindata_zdata[i * rows_in + j] * freq_weight_curve[j];
                aps_data_zdata[i] = Math.Sqrt(summe) / N_nonzero[i] * energycorr;
            }
            );
            OL = 0;
            for (int ii = 0; ii < aps_data_zdata.Length; ii++) OL += aps_data_zdata[ii];
            OL /= (double)aps_data_zdata.Length;
            return OL;
        }

        internal static XYZ_Data CalcNtelOctaveSpectra(double[] rawdata_y, XY_Data n100, double dt, double x0, CATimeStamps timestamps, Param_struct paramset)
        {
            //// Variable Definitions

            XYZ_Data aps_plaindata,     // X,Y,Z data from Plaindata Function
                result;                 // Result Data X, Y, Z

            F_octave f_octave;          // center frequencies for each ntel band

            int rows_out,               // number of band pass filters
                rows_in,                // number of frequency points (input data)
                ii;                     // general purpose running variable
            //kk;                     // running variable used for band pass summation

            int[] f1Pt,                 // f1 index of each band pass filter
                f2Pt;                   // f2 index of each band pass filter

            double df,                  // frequency resolution
                wincorr,                // window amplitude correction factor (e.g. 2 for Hann window)
                energycorrection,       // window energy correction factor (e.g. Math.Sqrt(8 / 3) / 2 for Hann window)
                rmsfactor,
                overalllevel;
            //BP_sum;                 // sum of one band pass

            //double[] onespec,           // single output spectrum of one leading quantity query point
            //aps_oneblock,           // single FFT input block
            double[] freq_weight_curve;      // Frequency weights vector (factor values)


            // FFT calculations

            aps_plaindata = FFT_Plaindata(rawdata_y, n100, dt, x0, timestamps, paramset);

            // Simple variable calculations

            rows_in = aps_plaindata.ydata.Length;
            int N = (int)(1 / dt / paramset.delta_f);
            int[] N_nonzero = aps_plaindata.N_nonzero;

            df = 1 / dt / N;

            if (paramset.windowtype == 1)
            {
                wincorr = 2;
                energycorrection = Math.Sqrt(8 / 3.0) / 2.0;
            }
            else
            {
                wincorr = 1;
                energycorrection = 1;
            }

            // Frequency weighting

            if (paramset.freq_weight == 0)
            {
                freq_weight_curve = new double[rows_in];
                for (ii = 0; ii < rows_in; ii++) freq_weight_curve[ii] = 1;
            }
            else freq_weight_curve = ABC_Weight(aps_plaindata.ydata, paramset.freq_weight);


            overalllevel = CalcOverallLevelSingleValue(aps_plaindata.zdata, aps_plaindata.ydata.Length, N_nonzero, freq_weight_curve, paramset.windowtype, aps_plaindata.xdata.Length, paramset.DC);
            
            // Calculation of band pass frequencies and indices

            f_octave = CalcCornerFreqs(paramset.n_octave, 1 / dt, (double)N);
            rows_out = f_octave.f0.Length;
            var f1Pt_tmp = new int[rows_out];
            var f2Pt_tmp = new int[rows_out];
            int k_start = rows_out;
            for (ii = rows_out - 1; ii >= 0; ii--)
            {
                f1Pt_tmp[ii] = (int)(f_octave.f1[ii] / df);
                f2Pt_tmp[ii] = (int)Math.Ceiling(f_octave.f2[ii] / df) - 1;
                if (f2Pt_tmp[ii] - f1Pt_tmp[ii] >= 0) k_start = ii;
                else break;
            }
            rows_out -= k_start;

            f1Pt = new int[rows_out];
            f2Pt = new int[rows_out];
            var f0 = new double[rows_out];
            Array.Copy(f1Pt_tmp, k_start, f1Pt, 0, rows_out);
            Array.Copy(f2Pt_tmp, k_start, f2Pt, 0, rows_out);
            Array.Copy(f_octave.f0, k_start, f0, 0, rows_out);

            // Memory allocation output data

            if (paramset.mean) result.zdata = new double[rows_out * aps_plaindata.xdata.Length];
            else result.zdata = new double[rows_out * aps_plaindata.xdata.Length];

            result.ydata = new double[aps_plaindata.xdata.Length * rows_out];
            result.xdata = new double[aps_plaindata.xdata.Length * rows_out];


            if (paramset.y_amplitude == 1)
                rmsfactor = Math.Sqrt(2);
            else if (paramset.y_amplitude == 2)
                rmsfactor = 2 * Math.Sqrt(2);
            else rmsfactor = 1;

            overalllevel *= rmsfactor;
            if (paramset.y_axis == 1) overalllevel = 20 * Math.Log10(overalllevel / paramset.dBref);
            // Output data calculation

            Parallel.For(0, aps_plaindata.xdata.Length, i =>
            {
                int j, k;
                double BP_sum;
                double[] onespec = new double[rows_out];
                double[] aps_oneblock = new double[rows_in];
                aps_oneblock[0] = Math.Sqrt(wincorr * aps_plaindata.zdata[i * rows_in]) * freq_weight_curve[0] / N_nonzero[i];
                for (j = 1; j < rows_in; j++) aps_oneblock[j] = Math.Sqrt(wincorr * aps_plaindata.zdata[i * rows_in + j]) * freq_weight_curve[j] * 2 / N_nonzero[i];
                if (paramset.y_axis == 1)
                {
                    for (j = 0; j < rows_out; j++)
                    {
                        BP_sum = 0;
                        for (k = f1Pt[j]; k <= f2Pt[j]; k++) BP_sum += Math.Pow(aps_oneblock[k], 2);
                        onespec[j] = Math.Sqrt(BP_sum);
                        if (f2Pt[j] > f1Pt[j]) onespec[j] *= energycorrection;
                        result.zdata[i * rows_out + j] = rmsfactor * 20 * Math.Log10(onespec[j] / paramset.dBref);
                        result.ydata[i * rows_out + j] = j + 1;
                        result.xdata[i * rows_out + j] = aps_plaindata.xdata[i];
                    }
                }
                else
                {
                    for (j = 0; j < rows_out; j++)
                    {
                        BP_sum = 0;
                        for (k = f1Pt[j]; k <= f2Pt[j]; k++) BP_sum += Math.Pow(aps_oneblock[k], 2);
                        onespec[j] = Math.Sqrt(BP_sum);
                        if (f2Pt[j] > f1Pt[j]) onespec[j] *= energycorrection;
                        result.zdata[i * rows_out + j] = rmsfactor * onespec[j];
                        result.ydata[i * rows_out + j] = j + 1;
                        result.xdata[i * rows_out + j] = aps_plaindata.xdata[i];
                    }
                }

            }
            );


            result.freq_labels = f0;
            for (int ij = 0; ij < result.freq_labels.Length; ij++)
            {
                if (Math.Round(f0[ij], MidpointRounding.AwayFromZero) > f0[ij]) result.freq_labels[ij] = Math.Round(f0[ij], 2, MidpointRounding.AwayFromZero);
                else result.freq_labels[ij] = Math.Round(f0[ij], 3, MidpointRounding.AwayFromZero);
            }
            result.dimensions = new int[2] { aps_plaindata.xdata.Length, rows_out };
            result.overalllevel = overalllevel;
            result.errorstring = null;
            //result.ColStart = new int[] { 0 };
            result.N_nonzero = new int[] { 0 };
            result.rows = 0;
            result.xtime = aps_plaindata.xtime;
            result.paramset = paramset;
            result.xdescription = "";
            // Averaging over time if 2D Diagram was chosen

            if (paramset.mean)
            {
                result = AverageDiagramOverTime(result, aps_plaindata.xdata);
            }



            return result;
        }

        internal static XYZ_Data CalcOrderSpectra(double[] rawdata_y, XY_Data n100, XY_Data enginespeed, double dt, double x0, CATimeStamps timestamps, Param_struct paramset)
        {


            XYZ_Data aps_plaindata,     // X,Y,Z data from Plaindata Function
                result;                 // Result Data X, Y, Z

            int rows_out,               // Number of orders to be calculated
                rows_in,                // Number of frequency points (plain data)
                ii;                     // Running variable
                                     // running variable used for band pass summation

            double wincorr,             // window amplitude correction factor (e.g. 2 for Hann window)
                energycorrection,       // window energy correction factor (e.g. Math.Sqrt(8 / 3) / 2 for Hann window)
                rmsfactor,
                overalllevel;
            //BP_sum;                 // sum of one band pass

            double[] oaxis,             // Array of the specific orders for calculation
                //aps_oneblock,           // single FFT input block
                freq_weight_curve;      // Frequency weights vector (factor values)
            //onespec;                // single output spectrum of one leading quantity query point

            //// Test parameters for inverter noise
            double f_t = paramset.inverter_frequency;
            int polepairs = paramset.polepairs;
            bool calcinverternoise = paramset.calcinverternoise;
            ////


            aps_plaindata = FFT_Plaindata(rawdata_y, n100, dt, x0, timestamps, paramset);
            int N = (int)(1 / dt / paramset.delta_f);
            int[] N_nonzero = aps_plaindata.N_nonzero;
            // Engine speed interpolation            
            XY_Data enginespeed_IP;
            /*if (paramset.LQ)*/ enginespeed_IP = Interp1(enginespeed.xdata, enginespeed.ydata, aps_plaindata.xtime);
            //else enginespeed_IP = Interp1(enginespeed.xdata, enginespeed.ydata, aps_plaindata.xdata);

            double[] apsz_neu = new double[enginespeed_IP.xdata.Length * aps_plaindata.ydata.Length];
            double[] aps_plaindata_xdata_neu = new double[enginespeed_IP.xdata.Length];
            int jj = 0;
            for (ii = 0; ii < aps_plaindata.xtime.Length; ii++)
            {
                if (aps_plaindata.xtime[ii] >= enginespeed_IP.xdata[0] && aps_plaindata.xtime[ii] <= enginespeed_IP.xdata[enginespeed_IP.xdata.Length - 1])
                {
                    aps_plaindata_xdata_neu[jj] = aps_plaindata.xdata[ii];
                    for (int kk = 0; kk < aps_plaindata.ydata.Length; kk++)apsz_neu[jj * aps_plaindata.ydata.Length + kk] = aps_plaindata.zdata[ii * aps_plaindata.ydata.Length + kk];
                    jj++;
                }
            }
            aps_plaindata.xdata = aps_plaindata_xdata_neu;
            aps_plaindata.zdata = apsz_neu;
            /*
            double[] apsz_neu = new double[enginespeed_IP.xdata.Length * aps_plaindata.ydata.Length];
            int jj = 0;
            for (ii = 0; ii < aps_plaindata.zdata.Length; ii++)
            {
                if (aps_plaindata.xdata[(int)(ii / aps_plaindata.ydata.Length)] >= enginespeed_IP.xdata[0])
                {
                    apsz_neu[jj] = aps_plaindata.zdata[ii];
                    jj++;
                }
            }*/
            //if (paramset.LQ) aps_plaindata.xdata = enginespeed_IP.ydata;
            //else aps_plaindata.xdata = enginespeed_IP.xdata;

            rows_in = aps_plaindata.ydata.Length;

            if (paramset.windowtype == 1)
            {
                wincorr = 2;
                energycorrection = Math.Sqrt(8 / 3.0) / 2.0;
            }
            else
            {
                wincorr = 1;
                energycorrection = 1;
            }

            // Frequency weighting

            if (paramset.freq_weight == 0)
            {
                freq_weight_curve = new double[rows_in];
                for (ii = 0; ii < rows_in; ii++) freq_weight_curve[ii] = 1;
            }
            else freq_weight_curve = ABC_Weight(aps_plaindata.ydata, paramset.freq_weight);

            overalllevel = CalcOverallLevelSingleValue(aps_plaindata.zdata, aps_plaindata.ydata.Length, N_nonzero, freq_weight_curve, paramset.windowtype, aps_plaindata_xdata_neu.Length, paramset.DC);
            

            // Order Parsing

            if (string.IsNullOrEmpty(paramset.orderstring))
            {
                rows_out = (int)((paramset.o2 - paramset.o1) / paramset.delta_o) + 1;
                oaxis = new double[rows_out];
                for (ii = 0; ii < rows_out; ii++) oaxis[ii] = ii * paramset.delta_o + paramset.o1;
            }
            else
            {
                oaxis = ParseOrdersFromString(paramset.orderstring);
                rows_out = oaxis.Length;
            }

            if (calcinverternoise)
            {
                rows_out *= 2;
                var oaxis_d_mirror = new double[rows_out];
                Array.Copy(oaxis, 0, oaxis_d_mirror, oaxis.Length, oaxis.Length);
                for (int ij = 1; ij <= oaxis.Length; ij++) oaxis_d_mirror[oaxis.Length - ij] = -oaxis[ij - 1];
                oaxis = oaxis_d_mirror;
            }


            // Memory allocation output data

            if (paramset.mean) result.zdata = new double[rows_out * aps_plaindata.xdata.Length];
            else result.zdata = new double[rows_out * aps_plaindata.xdata.Length];

            result.ydata = new double[aps_plaindata.xdata.Length * rows_out];
            result.xdata = new double[aps_plaindata.xdata.Length * rows_out];

            if (paramset.y_amplitude == 1)
                rmsfactor = Math.Sqrt(2);
            else if (paramset.y_amplitude == 2)
                rmsfactor = 2 * Math.Sqrt(2);
            else rmsfactor = 1;

            overalllevel *= rmsfactor;
            if (paramset.y_axis == 1) overalllevel = 20 * Math.Log10(overalllevel / paramset.dBref);

            
            
            
            // Output data calculation
            //for (int i = 0; i < aps_plaindata.xdata.Length; i++ )

            Parallel.For(0, aps_plaindata.xdata.Length, i =>
            {
                int j, k;
                double BP_sum;
                double[] onespec = new double[rows_out];
                double[] aps_oneblock = new double[rows_in];
                int[] f0Pt = new int[rows_out];
                int[] f1Pt = new int[rows_out];
                int[] f2Pt = new int[rows_out];
                for (j = 0; j < rows_out; j++)
                {
                    double frequencyoforder = enginespeed_IP.ydata[i] * oaxis[j] / 60.0;
                    if (calcinverternoise) frequencyoforder = frequencyoforder * polepairs +  f_t;
                    f0Pt[j] = (int)Math.Round(frequencyoforder / aps_plaindata.ydata[1], MidpointRounding.AwayFromZero) - 1;
                    f1Pt[j] = f0Pt[j] - paramset.tolerance;
                    f2Pt[j] = f0Pt[j] + paramset.tolerance;
                    if (f1Pt[j] < 0) f1Pt[j] = 0;
                    if (f2Pt[j] >= rows_in) f2Pt[j] = rows_in - 1;
                }

                aps_oneblock[0] = Math.Sqrt(wincorr * aps_plaindata.zdata[i * rows_in]) * freq_weight_curve[0] / N_nonzero[i];
                for (j = 1; j < rows_in; j++) aps_oneblock[j] = Math.Sqrt(wincorr * aps_plaindata.zdata[i * rows_in + j]) * freq_weight_curve[j] * 2 / N_nonzero[i];
                if (paramset.y_axis == 1)
                {
                    for (j = 0; j < rows_out; j++)
                    {
                        BP_sum = 0;
                        for (k = f1Pt[j]; k <= f2Pt[j]; k++) BP_sum += Math.Pow(aps_oneblock[k], 2);
                        onespec[j] = Math.Sqrt(BP_sum);
                        if (f2Pt[j] > f1Pt[j]) onespec[j] *= energycorrection;
                        result.zdata[i * rows_out + j] = rmsfactor * 20 * Math.Log10(onespec[j] / paramset.dBref);
                        result.ydata[i * rows_out + j] = j + 1;
                        result.xdata[i * rows_out + j] = aps_plaindata.xdata[i];
                    }
                }
                else
                {
                    for (j = 0; j < rows_out; j++)
                    {
                        BP_sum = 0;
                        for (k = f1Pt[j]; k <= f2Pt[j]; k++) BP_sum += Math.Pow(aps_oneblock[k], 2);
                        onespec[j] = Math.Sqrt(BP_sum);
                        if (f2Pt[j] > f1Pt[j]) onespec[j] *= energycorrection;
                        result.zdata[i * rows_out + j] = rmsfactor * onespec[j];
                        result.ydata[i * rows_out + j] = j + 1;
                        result.xdata[i * rows_out + j] = aps_plaindata.xdata[i];
                    }
                }
            }
            );
            
            // Averaging over time if 2D Diagram was chosen

            result.freq_labels = oaxis;
            result.dimensions = new int[2] { aps_plaindata.xdata.Length, rows_out };
            result.overalllevel = overalllevel;
            result.errorstring = null;
            //result.ColStart = new int[] { 0 };
            result.N_nonzero = new int[] { 0 };
            result.rows = 0;
            result.xtime = aps_plaindata.xtime;
            result.paramset = paramset;
            result.xdescription = "";
            if (paramset.mean)
            {
                result = AverageDiagramOverTime(result, aps_plaindata.xdata);
            }

            return result;

        }

        internal static XYZ_Data CalcOverallLevel(double[] rawdata_y, XY_Data n100, double dt, double x0, CATimeStamps timestamps, Param_struct paramset)
        {
            XYZ_Data result,
                aps_plaindata;

            int rows_in,
                ii;

            double rmsfactor;

            double[] freq_weight_curve;

            aps_plaindata = FFT_Plaindata(rawdata_y, n100, dt, x0, timestamps, paramset);

            int N = (int)(1 / dt / paramset.delta_f);
            int[] N_nonzero = aps_plaindata.N_nonzero;
            rows_in = aps_plaindata.ydata.Length;

            // Frequency weighting

            if (paramset.freq_weight == 0)
            {
                freq_weight_curve = new double[rows_in];
                for (ii = 0; ii < rows_in; ii++) freq_weight_curve[ii] = 1;
            }
            else freq_weight_curve = ABC_Weight(aps_plaindata.ydata, paramset.freq_weight);

            result.zdata = new double[aps_plaindata.xdata.Length];

            if (paramset.y_amplitude == 1)
                rmsfactor = Math.Sqrt(2);
            else if (paramset.y_amplitude == 2)
                rmsfactor = 2 * Math.Sqrt(2);
            else rmsfactor = 1;

            Parallel.For(0, aps_plaindata.xdata.Length, i =>
            //for (int i = 0; i < aps_plaindata.xdata.Length; i++)
            {
                int j;
                double levelsum = 0;
                //for (j = 1; j < rows_in; j++) levelsum += 2 * aps_plaindata.zdata[i * rows_in + j] * freq_weight_curve[j];
                for (j = 1; j < rows_in; j++) levelsum += 2 * aps_plaindata.zdata[i * rows_in + j] * Math.Pow(freq_weight_curve[j],2);
                levelsum += 0.5 * paramset.DC * aps_plaindata.zdata[i * rows_in] * freq_weight_curve[0];
                result.zdata[i] = rmsfactor * Math.Sqrt(levelsum) / N_nonzero[i];
                if (paramset.windowtype == 1) result.zdata[i] *= Math.Sqrt(8 / 3.0);
                if (paramset.y_axis == 1) result.zdata[i] = 20 * Math.Log10(result.zdata[i] / paramset.dBref);
            }
            );
            //;
            result.xdata = aps_plaindata.xdata;
            result.ydata = new double[1];
            result.ydata[0] = 0;
            result.freq_labels = new double[1];
            result.freq_labels[0] = 0;
            result.dimensions = new int[2] { aps_plaindata.xdata.Length, 1 };
            result.overalllevel = 0;
            for (ii = 0; ii < aps_plaindata.xdata.Length; ii++) result.overalllevel += result.zdata[ii];
            result.overalllevel /= (double)aps_plaindata.xdata.Length;
            

            result.errorstring = null;
            //result.ColStart = new int[] { 0 };
            result.N_nonzero = new int[] { 0 };
            result.rows = 0;
            result.xtime = aps_plaindata.xtime;
            result.paramset = paramset;
            result.xdescription = "";
            return result;
        }

        internal static XYZ_Data CalcBPLevel(double[] rawdata_y, XY_Data n100, double dt, double x0, CATimeStamps timestamps, Param_struct paramset)
        {
            XYZ_Data result,
                aps_plaindata;

            int rows_in,
                ii,
                f1Pt,
                f2Pt;

            double rmsfactor, overalllevel;

            double[] freq_weight_curve;

            aps_plaindata = FFT_Plaindata(rawdata_y, n100, dt, x0, timestamps, paramset);
            int N = (int)(1 / dt / paramset.delta_f);
            int[] N_nonzero = aps_plaindata.N_nonzero;
            rows_in = aps_plaindata.ydata.Length;

            // Frequency weighting

            if (paramset.freq_weight == 0)
            {
                freq_weight_curve = new double[rows_in];
                for (ii = 0; ii < rows_in; ii++) freq_weight_curve[ii] = 1;
            }
            else freq_weight_curve = ABC_Weight(aps_plaindata.ydata, paramset.freq_weight);

            overalllevel = CalcOverallLevelSingleValue(aps_plaindata.zdata, aps_plaindata.ydata.Length, N_nonzero, freq_weight_curve, paramset.windowtype, aps_plaindata.xdata.Length, paramset.DC);
            

            result.zdata = new double[aps_plaindata.xdata.Length];

            f1Pt = (int)Math.Round(paramset.f1 / paramset.delta_f, MidpointRounding.AwayFromZero);
            f2Pt = (int)Math.Round(paramset.f2 / paramset.delta_f, MidpointRounding.AwayFromZero);

            if (paramset.y_amplitude == 1)
                rmsfactor = Math.Sqrt(2);
            else if (paramset.y_amplitude == 2)
                rmsfactor = 2 * Math.Sqrt(2);
            else rmsfactor = 1;
            overalllevel *= rmsfactor;
            if (paramset.y_axis == 1) overalllevel = 20 * Math.Log10(overalllevel / paramset.dBref);

            Parallel.For(0, aps_plaindata.xdata.Length, i =>
            //for (ii = 0; ii < aps_plaindata.xdata.Length; ii++)
            {
                int j;
                double levelsum = 0;
                for (j = f1Pt; j <= f2Pt; j++) levelsum += 2 * aps_plaindata.zdata[i * rows_in + j] * freq_weight_curve[j];
                result.zdata[i] = rmsfactor * Math.Sqrt(levelsum) / N_nonzero[i];
                if (paramset.windowtype == 1) result.zdata[i] *= Math.Sqrt(8 / 3.0);
                if (paramset.y_axis == 1) result.zdata[i] = 20 * Math.Log10(result.zdata[i] / paramset.dBref);
            }
            );
            result.xdata = aps_plaindata.xdata;
            result.ydata = new double[1];
            result.ydata[0] = 0;
            result.freq_labels = new double[1];
            result.freq_labels[0] = 0;
            result.dimensions = new int[2] { aps_plaindata.xdata.Length, 1 };
            result.overalllevel = overalllevel;

            result.errorstring = null;
            //result.ColStart = new int[] { 0 };
            result.N_nonzero = new int[] { 0 };
            result.rows = 0;
            result.xtime = aps_plaindata.xtime;
            result.paramset = paramset;
            result.xdescription = "";
            return result;
        }
    }
}

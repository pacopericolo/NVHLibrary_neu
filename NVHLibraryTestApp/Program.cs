using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static NVHFunctions.Concerto;
using System.Xml.Serialization;
using NVHFunctions;
using System.IO;
using System.Xml;
using System.Globalization;
using IronPython;
using IronPython.Hosting;
using Microsoft.Scripting.Hosting;


namespace NVHLibraryTestApp
{
    class Program
    {
        static void Main()
        {

            var bx = new float[12];
            var ax = new float[12];

            LiquidDSPClass.liquid_iirdes(LiquidDSPClass.liquid_iirdes_filtertype.LIQUID_IIRDES_ELLIP, LiquidDSPClass.liquid_iirdes_bandtype.LIQUID_IIRDES_HIGHPASS, LiquidDSPClass.liquid_iirdes_format.LIQUID_IIRDES_SOS, 8, (float)(50.0/50000), 0.25f, 0.1f, 60f, bx, ax);








            var filepath = Path.GetFullPath(@"C:\Users\u12o24\Documents\Concerto_NVH");

            var rawdata = new double[500001];
            var sos1 = new double[24];
            var zi_tmp = new double[8];
            using (var br = new BinaryReader(new FileStream(Path.Combine(filepath, "rawdata.bin"), FileMode.Open)))
            {
                for (int ii = 0; ii < rawdata.Length; ii++) rawdata[ii] =  br.ReadDouble();
            };
            /*using (var br = new BinaryReader(new FileStream(Path.Combine(filepath, "sos1.bin"), FileMode.Open)))
            {
                for (int ii = 0; ii < sos1.Length; ii++) sos1[ii] = br.ReadDouble();
            };
            using (var br = new BinaryReader(new FileStream(Path.Combine(filepath, "zi.bin"), FileMode.Open)))
            {
                for (int ii = 0; ii < zi_tmp.Length; ii++) zi_tmp[ii] = br.ReadDouble();
            };*/
            using (var br = new BinaryReader(new FileStream(Path.Combine(filepath, "ziliquid.bin"), FileMode.Open)))
            {
                for (int ii = 0; ii < zi_tmp.Length; ii++) zi_tmp[ii] = br.ReadDouble();
            };

            using (var bw = new BinaryWriter(new FileStream(Path.Combine(filepath, "Liquid.bin"), FileMode.Create)))
            {
                for (int ii = 0; ii < bx.Length; ii++) bw.Write(bx[ii]);
                for (int ii = 0; ii < ax.Length; ii++) bw.Write(ax[ii]);
            };


            int L = 4;

            var b = new double[L][];
            var a = new double[L][];
            var zi = new double[L][];
            
            for (int ii = 0; ii < L; ii++)
            {
                b[ii] = new double[3];
                a[ii] = new double[3];
                zi[ii] = new double[2];

                for (int jj = 0; jj < 3; jj++)
                {
                    //b[ii][jj] = sos1[jj * L + ii];
                    //a[ii][jj] = sos1[jj * L + L * 3 + ii];
                    b[ii][jj] = bx[ii * 3 + jj];
                    a[ii][jj] = ax[ii * 3 + jj];
                    if (jj < 2) zi[ii][jj] = zi_tmp[ii * 2 + jj];
                }
            }
            var filtered = rawdata;

            for (int ii = 0; ii < L; ii++)
            {
                filtered = LiquidDSPClass.FiltFilt(filtered, b[ii], a[ii], zi[ii]);
            }

            using (var bw = new BinaryWriter(new FileStream(Path.Combine(filepath, "filtered_liquid.bin"), FileMode.Create)))
            {
                for (int ii = 0; ii < filtered.Length; ii++) bw.Write(filtered[ii]);
            };

            string test1 = "C:\\Programs\\pd";
            Console.WriteLine(test1.Replace('\\','/'));

            Microsoft.Scripting.Hosting.ScriptEngine pythonEngine =
                IronPython.Hosting.Python.CreateEngine();

            // Print the default search paths
            System.Console.Out.WriteLine("Search paths:");
            ICollection<string> searchPaths = pythonEngine.GetSearchPaths();
            foreach (string path in searchPaths)
            {
                System.Console.Out.WriteLine(path);
            }
            System.Console.Out.WriteLine();
            var doublearray = new double[3] { 3.5, 2.2, 0.7 };

            // Now modify the search paths to include the directory from
            // which we execute the script
            searchPaths.Add("..\\..");
            pythonEngine.SetSearchPaths(searchPaths);

            var scope = pythonEngine.CreateScope();
            scope.SetVariable("externalString", "how do you do?");
            scope.SetVariable("doublearray", doublearray);
            Microsoft.Scripting.Hosting.ScriptSource pythonScript =
                //pythonEngine.CreateScriptSourceFromString("print 'Hello World!'");
                pythonEngine.CreateScriptSourceFromFile("..\\..\\HelloWorld.py");
            pythonScript.Execute(scope);

            Console.WriteLine(scope.GetVariable("helloWorldString"));
            var arrayfrompython = (double[])scope.GetVariable("doublearray");
            
            Console.WriteLine(arrayfrompython[2]);

            Console.Read();

            //Console.WriteLine()

            /*var x = new double[10] { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
            var y = new double[10] { 5000, 6000, 4000, 3500, 3000, 2500, 1500, 1700, 1000, 500 };

            var test = Concerto.ProcessLQ_mono(new XY_Data { xdata = x, ydata = y });*/

            /*var arr1 = new int[3] { 7, 4, 3 };
            var arr2 = new int[3] { 1, 2, 3 };
            Array.Sort(arr1, arr2);*/

            //var test1 = GeneralMath.Nextpow2(50000 / 1);

            //var outstr = "";
            //NVHFunctions.F_octave f0 = Concerto.CalcCornerFreqs(6, 65536, 2*65536);
            //for (int ii = 0; ii < f0.f0.Length; ii++) if (f0.f0[ii] <= 20000 && f0.f0[ii] >= 20) outstr += Math.Round(f0.f0[ii], 1).ToString(CultureInfo.InvariantCulture) + ", ";
            //Console.WriteLine(outstr);
            string filename = "C:\\Users\\u12o24\\Documents\\nvhdebugdata.xls";
            NVHPackage obj;
            using (var filestream = new FileStream(filename, FileMode.Open))
            {
                using (var xmlreader = XmlReader.Create(filestream))
                {
                    XmlSerializer serializer = new XmlSerializer(typeof(NVHPackage));
                    obj = (NVHPackage)serializer.Deserialize(xmlreader);
                    xmlreader.Close();
                    filestream.Close();
                }
            }
            //obj.Paramset.Delta_t = 0.1;
            //obj.Paramset.Delta_f = 4;
            //var test = CalcAPSData(obj.Rawdata.ToArray(), obj.Enginespeed, obj.X0, new CATimeStamps(), obj.Paramset);
            //obj.Paramset.Delta_t = 1;
            //var test = CalcOrderSpectra(obj.Rawdata.ToArray(), new XY_Data(), obj.Enginespeed, obj.X0, new CATimeStamps(), obj.Paramset);
            //var test = GetXAxis((int)obj.Rawdata[0], new XY_Data(), obj.X0, new CATimeStamps(), obj.Paramset);

            //var test = CalcOverallLevel(obj.Rawdata.ToArray(), new XY_Data(), obj.X0, new CATimeStamps(), obj.Paramset);
            //Console.WriteLine(test.dimensions[0]);
            //Console.WriteLine(test.dimensions[1]);
            
            /*
            var x = new double[11];
            var y = new double[11];
            for (int ii = 0; ii < x.Length; ii++)
            {
                x[ii] = ii;
                y[ii] = Math.Sin(ii);
                Console.WriteLine("x: " + x[ii] + "\t y: " + y[ii]);
            }
            var xq = new double[49];
            for (int ii = 0; ii < xq.Length; ii++) xq[ii] = 0.25 * ii;
            var yq = GeneralMath.Interp1_linear(x, y, xq, true);
            for (int ii = 0; ii < yq.Length; ii++) Console.WriteLine("xq: " + xq[ii] + "\t yq:" + yq[ii]);

            var aps = new NVHLibrary.APS();
            var bplevel = new NVHLibrary.BandPassLevel();
            */

            /*double[] test = new double[5] { 4, 5, 11, 3, -2 };
            var result = NVHFunctions.GeneralMath.Cumsum(test);
            for (int ii = 0; ii < result.Length; ii++) Console.WriteLine(result[ii]);*/
            Console.ReadKey();
        }
    }
}

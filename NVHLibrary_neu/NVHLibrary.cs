
using System;
//using Avl.Concerto.DataExplorer;
//using Avl.Concerto.DotNetPlugin;
using Avl.Concerto.DataExplorer;

using System.Reflection;

using System.IO;
using Lomont;
using System.Threading.Tasks;
using System.Globalization;
using static NVHFunctions.GeneralMath;
using static NVHFunctions.WaveStream;
using static NVHFunctions.Concerto;
using NVHFunctions;
using System.Linq;
using System.Management.Instrumentation;
using System.Collections.Generic;
using Pd_Wrapper_x64;

namespace NVHLibrary
{

    internal struct XYZ_Concerto_Data
    {
        internal IDataSet xdata;
        internal IDataSet ydata;
        internal IDataSet zdata;
        internal IDataSet xtime;
        internal IDataSet freq_labels;
        internal IDataSet dimensions;
        internal double overalllevel;
    }


    public class BandPassLevel : APS
    {
        public double F1 { get; set; }
    }


    public class APS_Ring
    {
        private static readonly Param_struct paramset = new Param_struct
        {
            Delta_f = 4,
            Delta_t = 0.05,
            Windowtype = 1,
            Average = 1,
            Average_overlap = 0,
            LQ = 0,
            N_octave = 1,
            DC = 0,
            Y_axis = 1,
            Y_amplitude = 0,
            Diagramtype = 1,
            Freq_weight = 0,
            F1 = 0,
            F2 = 0,
            O1 = 0,
            O2 = 1,
            Delta_o = 1,
            Tolerance = 1,
            Mean = 0,
            Orderstring = "",
            Y_unitconversion = 1,
            DBref = 1e-06,
            Dsrange1 = double.NegativeInfinity,
            Dsrange2 = double.PositiveInfinity,
            Fs = 50000,
            FromCA = 0,
            ToCA = 360,
            Inverter_frequency = 0,
            Polepairs = 4
        };

        private static readonly double[][] xdata = new double[4][];
        private static readonly double[][] ydata = new double[4][];
        private static readonly double[][] zdata = new double[4][];
        /*private static double[] xdata = new double[193 * 6400];
        private static double[] ydata = new double[193 * 6400];
        private static double[] zdata = new double[193 * 6400];*/
        private static int ringindex = 0;
        private static double last_x = 0;
        private static int LRing = 193;

        

        public void Reset(int L)
        {
            LRing = L;
            ringindex = 0;
            for (int ii = 0; ii < 1; ii++)
            {
                xdata[ii] = new double[LRing * 6400];
                ydata[ii] = new double[LRing * 6400];
                zdata[ii] = new double[LRing * 6400];
            }
            
            last_x = 0;
        }

        public double LastXValue()
        {
            return last_x;
        }

        public string Calculate(IDataSet newdata)
        {
            string result = "";
            
            paramset.Dsrange1 = double.NegativeInfinity;
            paramset.Dsrange2 = double.PositiveInfinity;
            try
            {
                var apsresult = CalcAPSData_Ring(newdata.RealValues, newdata.x[0] / 1000, paramset, xdata[0], ydata[0], zdata[0], ringindex, LRing);
                ringindex += apsresult.dimensions[0];
                while (ringindex > LRing - 1) ringindex -= LRing;
                if (ringindex == 0) last_x = apsresult.xdata[apsresult.xdata.Length - 1];
                else last_x = apsresult.xdata[ringindex * 6400 - 1];
                xdata[0] = apsresult.xdata;
                ydata[0] = apsresult.ydata;
                zdata[0] = apsresult.zdata;
            }
            catch (Exception e)
            {
                result = e.Message;
                result += "Ringindex: " + ringindex;
            }
            
            

            return result;
        }

        public IDataSet GetXData()
        {
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = xdata.Length,
                RealValues = xdata[0]
            };
            
            return result;
        }

        public IDataSet GetYData()
        {
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = ydata.Length,
                RealValues = ydata[0]
            };
            return result;
        }
        public IDataSet GetZData()
        {
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = zdata.Length,
                RealValues = zdata[0]
            };
            return result;
        }


    }

    public class APS
    {
        double[] buffer1;
        CycleData cadatainput;
        CycleData[] cadatainputarray;
        int yunit = 0;
        //static IDataSet[] datasetarray = new IDataSet[32];
        static readonly double[][] datasetarray = new double[32][];
        static readonly XYZ_Data[] resultarray = new XYZ_Data[32];
        static private readonly double[] inlet_x0 = new double[32] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN };
        static readonly int[] calcflag = new int[32];
        static private readonly double[] inlet_dt = new double[32];
        
        

        /// <summary>
        /// Allocates memory for CA data in case of Concerto formula.
        /// </summary>
        /// <param name="count">Count of raw data.</param>
        /// <param name="cyclecount">Cycle count of raw data.</param>
        /// <param name="arraymode">For fast forward calculation if LQ == 2 and Delta_t > 2. 0: no fast-forward calculation (Delta_t <= 2), 1: fast-forward calculation (Delta_t > 2)</param>
        /// <param name="lastcyclesingular">Defines if the last calculated cycle is the last cycle of the measurement. 0: No, 1: Yes</param>
        public void ReadCyclesInit(int count, int cyclecount, int arraymode, int lastcyclesingular)
        {
            buffer1 = new double[count];
            if (arraymode == 1)
            {
                cadatainputarray = new NVHFunctions.CycleData[cyclecount];
                for (int ii = 0; ii < cyclecount - 1; ii++)
                {
                    cadatainputarray[ii].Count = count;
                    cadatainputarray[ii].CycleCount = 2;
                    cadatainputarray[ii].ydata = new double[count * 2];
                }
                if (lastcyclesingular == 1)
                {
                    cadatainputarray[cadatainputarray.Length - 1].Count = count;
                    cadatainputarray[cadatainputarray.Length - 1].CycleCount = 1;
                    cadatainputarray[cadatainputarray.Length - 1].ydata = new double[count];
                }
                else
                {
                    cadatainputarray[cadatainputarray.Length - 1].Count = count;
                    cadatainputarray[cadatainputarray.Length - 1].CycleCount = 2;
                    cadatainputarray[cadatainputarray.Length - 1].ydata = new double[count * 2];
                }
            }
            else
            {
                cadatainput = new NVHFunctions.CycleData
                {
                    Count = count,
                    CycleCount = cyclecount,
                    xdata = new double[count],
                    ydata = new double[count * cyclecount]
                };
            }

        }

        /// <summary>
        /// Writes raw data of one cycle into cadatainput/cadatainputarray. Only for CA data from Concerto formula.
        /// </summary>
        /// <param name="data">Dataset of Raw data.</param>
        /// <param name="arraymode">For fast forward calculation if LQ == 2 and Delta_t > 2. 0: no fast-forward calculation (Delta_t <= 2), 1: fast-forward calculation (Delta_t > 2)</param>
        /// <param name="arrayindex">Index of the current cycle.</param>
        /// <param name="cycleindex">Defines if this is the current or the consecutive cycle (in each array position, 2 cycles are written to take filter response into account). 0: current cycle, 1: consecutive cycle.</param>
        /// <param name="cyclecount">Number of cycles of raw data. Has to be provided, as it's only written into first cycle. Workaround hack.</param>
        public void ReadOneCycle(IDataSet data, int arraymode, int arrayindex, int cycleindex, int cyclecount)
        {
            if (arraymode == 1)
            {
                cadatainputarray[arrayindex].xdata = data.x.RealValues;
                buffer1 = data.RealValues;
                if (cycleindex == 0) for (int ii = 0; ii < data.Count; ii++) cadatainputarray[arrayindex].ydata[ii] = buffer1[ii];
                else for (int ii = 0; ii < data.Count; ii++) cadatainputarray[arrayindex].ydata[ii + data.Count] = buffer1[ii];
            }
            else
            {
                if (cycleindex == 0)
                {
                    cadatainput.Count = data.Count;
                    cadatainput.CycleCount = cyclecount;
                    cadatainput.xdata = data.x.RealValues;
                }
                //int kk = cycleindex - cyclestart;
                buffer1 = data.RealValues;
                for (int jj = 0; jj < data.Count; jj++) cadatainput.ydata[cycleindex * data.Count + jj] = buffer1[jj];
            }


        }

        public void AddRawDatasets(IDataSet data1, IDataSet data2, IDataSet data3, IDataSet data4, IDataSet data5, IDataSet data6, IDataSet data7, IDataSet x0)
        {
            var x0r = x0.RealValues;
            IDataSet[] dataSets = new IDataSet[7];
            dataSets[0] = data1;
            dataSets[1] = data2;
            dataSets[2] = data3;
            dataSets[3] = data4;
            dataSets[4] = data5;
            dataSets[5] = data6;
            dataSets[6] = data7;

            //for (int i = 0; i < 7; i++) datasetarray[i] = new double[dataSets[i].Count];

            //Parallel.For(0, 7, i =>
            for (int i = 0; i < 7; i++)
             {
                 AddRawDataset(dataSets[i], i, x0r[i]);
                 
             }
            //);
        }

        public void AddRawDataset(IDataSet data, int identifier, double x0_ds)
        {
            if (data.Count == 0)
            {
                datasetarray[identifier] = null;
                calcflag[identifier] = 0;
            }
            else
            {
                double x0 = x0_ds;
                if (x0 != inlet_x0[identifier])
                {
                    inlet_x0[identifier] = x0_ds;
                    inlet_dt[identifier] = data.x[1] - data.x[0];
                    datasetarray[identifier] = data.RealValues;
                    //Array.Copy(data.RealValues, datasetarray[identifier], data.Count);
                    calcflag[identifier] = 0;
                }
                else calcflag[identifier] = 1;
            }
            data.Release();
            
            
        }

        public IDataSet XData(int identifier)
        {
            if (resultarray[identifier].xdata != null)
            {
                var xdata = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = resultarray[identifier].xdata,
                    Count = resultarray[identifier].xdata.Length
                };
                return xdata;
            }
            else return new IDataSet(DataSetType.Numeric)
            {
                Count = 0
            };
        }
        public IDataSet YData(int identifier)
        {
            if (resultarray[identifier].ydata != null)
            {
                var ydata = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = resultarray[identifier].ydata,
                    Count = resultarray[identifier].ydata.Length
                };
                return ydata;
            }
            else return new IDataSet(DataSetType.Numeric)
            {
                Count = 0
            };
        }
        public IDataSet ZData(int identifier)
        {
            if (resultarray[identifier].zdata != null)
            {
                var zdata = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = resultarray[identifier].zdata,
                    Count = resultarray[identifier].zdata.Length
                };
                return zdata;
            }
            else return new IDataSet(DataSetType.Numeric)
            {
                Count = 0
            };
        }
        public IDataSet XTime(int identifier)
        {
            if (resultarray[identifier].xtime != null)
            {
                var xtime = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = resultarray[identifier].xtime,
                    Count = resultarray[identifier].xtime.Length
                };
                return xtime;
            }
            else return new IDataSet(DataSetType.Numeric)
            {
                Count = 0
            };
        }
        public IDataSet Dimensions(int identifier)
        {
            if (resultarray[identifier].dimensions != null)
            {

                var dimensions = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = new double[2] { resultarray[identifier].dimensions[0], resultarray[identifier].dimensions[1] },
                    Count = 2
                };
                return dimensions;
            }
            else return new IDataSet(DataSetType.Numeric)
            {
                Count = 0
            };
        }

        public double Overalllevel(int identifier)
        {
            
            return resultarray[identifier].overalllevel;
        }
        /*public IDataSet Overalllevel(int identifier)
        {
            var ol = new double[1] { resultarray[identifier].overalllevel};
            return new IDataSet(DataSetType.Numeric)
            {
                Count = 1,
                RealValues = ol
            };
        }*/


        /*
        public IDataSet YData { get; private set; } = new IDataSet(DataSetType.Numeric);
        public IDataSet ZData { get; private set; } = new IDataSet(DataSetType.Numeric);
        public IDataSet XTime { get; private set; } = new IDataSet(DataSetType.Numeric);
        public IDataSet Dimensions { get; private set; } = new IDataSet(DataSetType.Numeric);
        public double Overalllevel { get; private set; } = 0;*/


        /// <summary>
        /// Used leading quantity instead of time if LQ == 1. Sidenote: Dsrange and Delta_t also have to be supplied in the appropriate units!.
        /// </summary>
        public IDataSet N100 { private get; set; } = new IDataSet(DataSetType.Numeric);
        /// <summary>
        /// In case of CA based analysis (LQ == 2), this signal provides the connection between crank angle and time. Can be dT (CDMTIME == 0) or CDMTIME (CDMTIME == 1).
        /// </summary>
        public IDataSet DT { private get; set; } = new IDataSet(DataSetType.Numeric);
        /// <summary>
        /// In case of CA based analysis (LQ == 2) and dT signal (CDMTIME == 0), this signal provides the starting time for each cycle.
        /// </summary>
        public IDataSet Cyctime { private get; set; } = new IDataSet(DataSetType.Numeric);
        /// <summary>
        /// Frequency resolution in Hz.
        /// </summary>
        public double Delta_f { get; set; } = 0;
        /// <summary>
        /// Time/leading quantity resolution. Needs to be provided in the same unit as dsrange and n100.
        /// </summary>
        public double Delta_t { get; set; } = 0;
        /// <summary>
        /// Window type. 0: rectangular, 1: Hann (default)
        /// </summary>
        /// 
        public int Ms { get; set; } = 1;

        public int Windowtype { get; set; } = 1;
        /// <summary>
        /// Defines if the CA time information signal is CDMTIME or dT. 0: CDMTIME (default), 1: dT
        /// </summary>
        public bool CDMTIME { get; set; } = true;
        /// <summary>
        /// Number of sub-blocks to be averaged to one main block. Default is 1 (i.e. no sub-block averaging).
        /// </summary>
        public int Average { get; set; } = 1;
        /// <summary>
        /// Overlap of the averaged sub-blocks. Irrelevant if Average is 1. Default is 0.
        /// </summary>
        public double Average_overlap { get; set; } = 0;
        /// <summary>
        /// Leading quantity. 0: Time (default), 1: other signal if provided as n100, 2: Cycles
        /// </summary>
        public int LQ { get; set; } = 0;
        /// <summary>
        /// Defines if DC should be included in the result. 0 (default): DC is not included, 1: DC is included
        /// </summary>
        public int DC { get; set; } = 0;
        /// <summary>
        /// Defines the values of the magnitude axis. 0: linear, 1 (default): logarithmic - dB
        /// </summary>
        public int Y_axis { get; set; } = 1;
        /// <summary>
        /// Defines the values of the magnitude axis. 0 (default): RMS, 1: Peak, 2: Peak-Peak
        /// </summary>
        public int Y_amplitude { get; set; } = 0;
        /// <summary>
        /// Frequency weighting. 0: linear, 1 (default): A-weighted, 2: B-weighted, 3: C-weighted
        /// </summary>
        public int Freq_weight { get; set; } = 1;
        /// <summary>
        /// Defines if the result should be averaged over the complete leading quantity. Result would be only one 2D spectrum. 0 (default): 3D-Diagram, 1: 2-D Diagram.
        /// </summary>
        public int Mean { get; set; } = 1;
        /// <summary>
        /// Lower leading quantity range for calculation. Rawdata below this value is omitted. Default is -Infinity.
        /// </summary>
        public double Dsrange1 { get; set; } = double.NegativeInfinity;
        /// <summary>
        /// Upper leading quantity range for calculation. Rawdata above this value is omitted. Default is Infinity.
        /// </summary>
        public double Dsrange2 { get; set; } = double.PositiveInfinity;

        public int Cyclerange1 { get; set; } = 0;
        public int Cyclerange2 { get; set; } = int.MaxValue;
        /// <summary>
        /// Desired Sampling frequency in Hz. Default is "unchanged".
        /// </summary>
        public double Fs { get; set; } = 0;
        /// <summary>
        /// Lower crank angle range in degrees for special CA analysis. Data below this CA is omitted and if data is too small for FFT size, zero padding is applied. Default is -360°.
        /// </summary>
        public int FromCA { get; set; } = -360;
        /// <summary>
        /// Upper crank angle range in degrees for special CA analysis. Data above this CA is omitted and if data is too small for FFT size, zero padding is applied. Default is 360°.
        /// </summary>
        public int ToCA { get; set; } = 360;
        /// <summary>
        /// Information if rawdata is a Concerto formula or not. 0 (default): rawdata is no formula. 1: rawdata is a formula (certain calculation acceleration features are deactivated).
        /// </summary>
        public bool RawdataIsFormula { get; set; } = false;


        public string Calculate()
        {


            for (int ii = 0; ii < calcflag.Length; ii++) if (calcflag[ii] == 1) return "no calculation.";





            //for (int ii = 0; ii < inletpassed.Length; ii++) if (inletpassed[ii] == 0) return "not all inlet passed yet.";

            double time_factor;
            if (Ms == 0) time_factor = 1;
            else time_factor = 0.001;
            if (Average == 0) Average = 1;
            //if (Delta_f == 0) Delta_f = (1 / time_factor / (Rawdata.x[1] - Rawdata.x[0])) / 16384;
            Delta_f = 4;




            /*if (Delta_t == 0)
            {
                if (LQ == 1) Delta_t = (N100.x[N100.Count - 1] - N100.x[0]) / 500;
                else if (LQ == 2) Delta_t = Math.Max(Rawdata.CycleCount / 500, 1);
                else Delta_t = (Rawdata.x[Rawdata.Count - 1] - Rawdata.x[0]) * time_factor / 500;
            }*/

            //Parallel.For(0, 32, i =>
            //try
            //{
            /*if (concat == 1)
            {
                int L = 0;
                var pos = new int[10];
                for (int i = 0; i < 10; i++)
                {
                    pos[i] = L;
                    L += datasetarray[i].Count;
                }
                var bigdata = new double[L];
                var xdata = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = new double[2] { datasetarray[0].x[0], datasetarray[0].x[1] },
                    Count = 2
                };

                Parallel.For(0, 10, i =>
                //for (int i = 0; i < 10; i++)
                {
                    var data = datasetarray[i].RealValues;
                    for (int j = 0; j < datasetarray[i].Count; j++) bigdata[pos[i] + j] = data[j];
                    
                }
                );
                
                for (int i = 0; i < 10; i++) datasetarray[i] = null;

                datasetarray[0] = new IDataSet(DataSetType.Numeric)
                {
                    Count = L,
                    RealValues = bigdata,
                    x = xdata
                };

            }*/

            Parallel.For(0, 32, i =>

            //    for (int i = 0; i < 32; i++)
                {
                    //IDataSet Rawdata;
                    double[] Rawdata;
                    if (datasetarray[i] != null)
                    {
                        Rawdata = datasetarray[i];

                        //if (Rawdata.Units.Equals("m/s^2", StringComparison.OrdinalIgnoreCase)) yunit = 1;
                        //else if (Rawdata.Units.Equals("g", StringComparison.OrdinalIgnoreCase)) yunit = 2;
                        yunit = 1;

                        /*if (Rawdata.Count == 0) return "No rawdata provided.";
                        if (LQ == 1 && N100.Count == 0) return "Leading quantity selected, but no LQ data provided.";
                        if (LQ == 2 && DT.Count == 0) return "CA based analysis selected, but no dT or CDMTIME signal provided.";
                        if (!CDMTIME && Cyctime.Count == 0 && LQ == 2) return "For CA based analyis and dT signal, Cyctime data is needed additionaly.";
                        if (RawdataIsFormula && Rawdata.CycleCount > 1)
                        {
                            if (LQ == 2 && Delta_t > 2 && cadatainputarray == null) return "No raw data found. For CA based data from Concerto formula, use the ReadCyclesInit and ReadOneCycle method.";
                            else if (cadatainput == null) return "No raw data found. For CA based data from Concerto formula, use the ReadCyclesInit and ReadOneCycle method.";
                        }*/


                        double dbref = 0.00002;
                        double unitconversion = 1;
                        if (yunit == 2) unitconversion = 9.81;
                        if (yunit > 0) dbref = 0.0000001;



                        //if (Cyclerange1 < 1) Cyclerange1 = 1;
                        //if (Cyclerange2 > Rawdata.CycleCount) Cyclerange2 = Rawdata.CycleCount;

                        double _Fs;
                        if (Fs == 0) _Fs = 1 / time_factor / inlet_dt[i];
                        else _Fs = 50000;

                        var paramset = new Param_struct
                        {
                            Average = Average,
                            Average_overlap = Average_overlap,
                            Calcinverternoise = false,
                            CDMTIME = CDMTIME,
                            Cyclerange1 = Cyclerange1 - 1,
                            Cyclerange2 = Cyclerange2 - 1,
                            DBref = dbref,
                            DC = DC,
                            Delta_f = Delta_f,
                            Delta_t = Delta_t,
                            Diagramtype = 1,
                            Dsrange1 = Dsrange1,
                            Dsrange2 = Dsrange2,
                            Freq_weight = Freq_weight,
                            FromCA = FromCA,
                            ToCA = ToCA,
                            Fs = _Fs,
                            LQ = LQ,
                            Mean = Mean,
                            Windowtype = Windowtype,
                            Y_amplitude = Y_amplitude,
                            Y_axis = Y_axis,
                            Y_unitconversion = unitconversion
                        };

                        //if (Mean == 0) paramset.Mean = false;

                        //bool formulacalc = false;

                        //try
                        //{
                        //var apsresult = ConcertoWrapper.CalculateDiagram(Rawdata, DT, Cyctime, N100, new IDataSet(DataSetType.Numeric), paramset, formulacalc, time_factor, cadatainput, cadatainputarray);
                        var apsresult = ConcertoWrapper.CalculateDiagram(Rawdata, inlet_x0[i]*time_factor,inlet_dt[i]*time_factor, new IDataSet(DataSetType.Numeric), new IDataSet(DataSetType.Numeric), new IDataSet(DataSetType.Numeric), new IDataSet(DataSetType.Numeric), paramset, time_factor);
                        resultarray[i] = apsresult;
                        //}
                        /*catch (Exception e)
                        {
                            return e.Message;
                        }*/
                        
                        /*XData = apsresult.xdata;
                        YData = apsresult.ydata;
                        ZData = apsresult.zdata;
                        XTime = apsresult.xtime;
                        Dimensions = apsresult.dimensions;
                        Overalllevel = apsresult.overalllevel;*/
                    }
                    //datasetarray[i] = null;
                    calcflag[i] = 1;
                }
            );
            //}


            //inletpassed = new int[10];

            return "Calculation Done!";
        }
    }


    public class DataIOFunctions
    {
        FileStream output;
        BinaryWriter binWrt;


        


        public int SaveData(IDataSet data, String filename)
        {
            /*
             * This version uses cycle indexing within C# (not supported for old Concerto versions)
             * 
            */
            
            output = new FileStream(filename, FileMode.Create);
            binWrt = new BinaryWriter(output);
            binWrt.Write(data.CycleCount);
            binWrt.Write(data.Count);
            
            var buffer = data.x.RealValues;
            for (int i = 0; i < data.Count; i++) binWrt.Write(buffer[i]);
            for (int jj = 0; jj < data.CycleCount; jj++)
            {
                data.CycleIndex = jj;
                buffer = data.RealValues;
                for (int i = 0; i < data.Count; i++) binWrt.Write(buffer[i]);
            }
            output.Close();
            
            return 0;

        }

        internal XY_Data ReadConcertoData(IDataSet data)
        {

            // Don't change this to xdata = data.x.RealValues!! (Doesn't work!)
            var doubledata = new XY_Data
            {
                xdata = new double[data.Count],
                ydata = data.RealValues
            };
            
            for (int ii = 0; ii < data.Count; ii++)
            {
                doubledata.xdata[ii] = data.x[ii];
            }

            return doubledata;
        }



        internal double[] ReadConcertoYData(IDataSet data)
        {
            double[] doubledata = data.RealValues;
            return doubledata;
        }

        internal CycleData ReadConcertoCAData(IDataSet data, int cycle_start, int cycle_end, int delta = 1)
        {

            int cyclecount = (cycle_end - cycle_start) / delta + 1;
            var result = new NVHFunctions.CycleData
            {
                xdata = data.x.RealValues,
                ydata = new double[data.Count * cyclecount],
                Count = data.Count,
                CycleCount = cyclecount
            };
            double[] onearray;

            int kk = 0;
            for (int ii = cycle_start; ii <= cycle_end; ii+= delta)
            {
                data.CycleIndex = ii;
                onearray = data.RealValues;
                for (int jj = 0; jj < data.Count; jj++) result.ydata[kk * data.Count + jj] = onearray[jj];
                kk++;
            }
            
            return result;
        }

        
        internal double[] ReadData(string filename)
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
                byte[] bytes = br.ReadBytes((int)length * 8 * 2);
                
                values = new double[bytes.Length / 8];
                for (int ii = 0; ii < values.Length; ii++) values[ii] = BitConverter.ToDouble(bytes, ii * 8);
                

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

        internal double[] ReadYData(string filename)
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
                byte[] bytes = br.ReadBytes((int)length * 8 * 2);

                values = new double[bytes.Length / 16];
                for (int ii = values.Length / 2; ii < values.Length; ii++) values[ii] = BitConverter.ToDouble(bytes, ii * 8);


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

    }


    public class WaveExport
    {

        private static double[][] audiochannels;
        private static int currentchannel = 0;
        private static int L = 0;
        private static int fs = 0;
        private static double[] calibvalues;
        private static string[] units;
        private static string[] channelnames;

        public string InitChannels(int nchannels, int samplerate)
        {
            audiochannels = new double[nchannels][];
            calibvalues = new double[nchannels];
            units = new string[nchannels];
            channelnames = new string[nchannels];
            currentchannel = 0;
            L = 0;
            fs = samplerate;
            return "Opened " + nchannels + " channels.";
        }

        public string AddChannel(IDataSet channel, string unit = "Pa", string name = "")
        {
            if (name == "") channelnames[currentchannel] = "CH" + currentchannel;
            else channelnames[currentchannel] = name;
            units[currentchannel] = unit;
            double maxval = Max(channel.RealValues);
            double minval = Min(channel.RealValues);
            double absmaxval = Math.Max(Math.Abs(maxval), Math.Abs(minval));
            double VperPa = 0.99 / absmaxval;
            calibvalues[currentchannel] = VperPa;
            if (currentchannel == 0) L = channel.Count;
            if (L != channel.Count) return "Error: Channel count doesn't fit.";
            audiochannels[currentchannel] = channel.RealValues;
            currentchannel++;
            return "Added Channel " + name + " with " + VperPa + "V/" + unit;
        }

        public string WriteAudioFile(string filename)
        {
            var nchannels = audiochannels.Length;

            using (var writer = new NAudio.Wave.WaveFileWriter(filename, new NAudio.Wave.WaveFormat(fs, 24, nchannels)))
            {
                var floatarray = new float[nchannels * L];
                for (int ii = 0; ii < L; ii++)
                {
                    for (int ch = 0; ch < nchannels; ch++)
                    {
                        floatarray[ii * nchannels + ch] = (float)(audiochannels[ch][ii] * calibvalues[ch]);
                    }
                }
                writer.WriteSamples(floatarray, 0, floatarray.Length);
            }
                
            
            using (var tw = new StreamWriter(new FileStream(Path.Combine(Path.GetDirectoryName(filename), Path.GetFileNameWithoutExtension(filename) + ".calib"), FileMode.Create)))
            {
                for (int ii = 0; ii < calibvalues.Length; ii++) tw.WriteLine(channelnames[ii] + "\t" + calibvalues[ii].ToString(CultureInfo.InvariantCulture) + "\tV/" + units[ii]);
            }

            /*    var result = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = calibvalues,
                    Count = calibvalues.Length
                };
            return result;*/
            currentchannel = 0;
            L = 0;
            fs = 0;
            calibvalues = new double[0];
            audiochannels = new double[0][];
            return "Successfully written " + nchannels + " channels to file " + filename;
        }


    }

    public class WaveImport
    {
        HeaderObject headerobject;
        DataObject dataobject;
        //WaveFunctions.WaveStream.WaveFileObject waveobject;
        bool dataread = false;
        bool headerread = false;
        
        string filename_prop;
        //double dbref_prop;
        readonly IDataSet dummy = new IDataSet(DataSetType.Numeric)
        {
            Count = 0
        };

        public void SetProperties(double calib_l, double calib_r, int speedpulses_l, int speedpulses_r, double speedfactor_l, double speedfactor_r, double dbref)
        {
            if (calib_l == -1) headerobject.calib_l = Math.Pow(2, headerobject.nbits - 1);
            else headerobject.calib_l = Math.Pow(2, 23) / (Math.Pow(10, calib_l / 20) * dbref * 2 * Math.Sqrt(2));
            if (calib_r == -1) headerobject.calib_r = Math.Pow(2, headerobject.nbits - 1);
            else headerobject.calib_r = Math.Pow(2, 23) / (Math.Pow(10, calib_r / 20) * dbref * 2 * Math.Sqrt(2));
            headerobject.speedpulses_l = speedpulses_l;
            headerobject.speedpulses_r = speedpulses_r;
            headerobject.speedfactor_l = speedfactor_l;
            headerobject.speedfactor_r = speedfactor_r;
        }

        /*public void CollectData(string filename, double calib_l, double calib_r, int speedpulses_l, int speedpulses_r, double speedfactor_l, double speedfactor_r, double dbref)
        {
            var wavestreamclass = new WaveFunctions.WaveStream();
            waveobject = wavestreamclass.ReadWaveFile(filename, calib_l: calib_l, calib_r: calib_r, speedpulses_l: speedpulses_l, speedpulses_r: speedpulses_r, speedfactor_l: speedfactor_l, speedfactor_r: speedfactor_r, dbref: dbref);
            waveobject.calib_l = 20 * Math.Log10(Math.Pow(2, waveobject.nbits - 1) / waveobject.calib_l / dbref / 2 / Math.Sqrt(2));
            waveobject.calib_r = 20 * Math.Log10(Math.Pow(2, waveobject.nbits - 1) / waveobject.calib_r / dbref / 2 / Math.Sqrt(2));
            datacollected = true;
        }*/

        public IDataSet GetLeftChannel()
        {
            if (!headerread) return dummy;
            if (!dataread)
            {
                
                dataobject = ReadWaveData(filename_prop, headerobject);
                dataread = true;
            }
            
            var result_double = new double[dataobject.audio_L.Length];
            var time = new double[result_double.Length];
            for (int ii = 0; ii < result_double.Length; ii++)
            {
                result_double[ii] = dataobject.audio_L[ii] / headerobject.calib_l * Math.Pow(2, 16) / Math.Pow(2, headerobject.nbits);
                time[ii] = ii / (double)headerobject.samplerate;
            }

            var result_x = new IDataSet(DataSetType.Numeric)
            {
                Count = time.Length,
                RealValues = time
            };

            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = dataobject.audio_L.Length,
                RealValues = result_double,
                x = result_x
            };
            return result;
        }

        public IDataSet GetRightChannel()
        {
            if (!headerread) return dummy;
            if (headerobject.nchannels < 2) return dummy;
            
            if (!dataread)
            {
                
                dataobject = ReadWaveData(filename_prop, headerobject);
                dataread = true;
            }

            var result_double = new double[dataobject.audio_R.Length];
            var time = new double[result_double.Length];
            for (int ii = 0; ii < result_double.Length; ii++)
            {
                result_double[ii] = dataobject.audio_R[ii] / headerobject.calib_r * Math.Pow(2, 16) / Math.Pow(2, headerobject.nbits);
                time[ii] = ii / (double)headerobject.samplerate;
            }

            var result_x = new IDataSet(DataSetType.Numeric)
            {
                Count = time.Length,
                RealValues = time
            };

            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = dataobject.audio_R.Length,
                RealValues = result_double,
                x = result_x
            };
            return result;
        }

        public IDataSet GetLeftSpeedChannel()
        {
            if (!headerread) return dummy;
            if (headerobject.speedpulses_l == -1) return dummy;

            if (!dataread)
            {
                dataobject = ReadWaveData(filename_prop, headerobject);
                dataread = true;
            }

            var result_x = new IDataSet(DataSetType.Numeric)
            {
                Count = dataobject.rpm_l.xdata.Length,
                RealValues = dataobject.rpm_l.xdata
            };

            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = dataobject.rpm_l.ydata.Length,
                RealValues = dataobject.rpm_l.ydata,
                x = result_x
            };
            return result;

        }


        public IDataSet GetRightSpeedChannel()
        {
            if (!headerread) return dummy;
            if (headerobject.speedpulses_r == -1) return dummy;

            if (!dataread)
            {
                dataobject = ReadWaveData(filename_prop, headerobject);
                dataread = true;
            }

            var result_x = new IDataSet(DataSetType.Numeric)
            {
                Count = dataobject.rpm_r.xdata.Length,
                RealValues = dataobject.rpm_r.xdata
            };

            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = dataobject.rpm_r.ydata.Length,
                RealValues = dataobject.rpm_r.ydata,
                x = result_x
            };
            return result;
            
        }

        public IDataSet ReadHeader(string filename, double dbref)
        {
            filename_prop = filename;
            //dbref_prop = dbref;
            headerobject = ReadWaveHeader(filename);
            //if (headerobject.calib_l == -1) headerobject.calib_l = Math.Pow(2, headerobject.nbits - 1);
            //if (headerobject.calib_r == -1) headerobject.calib_r = Math.Pow(2, headerobject.nbits - 1);
            if (headerobject.calib_l != -1) headerobject.calib_l = 20 * Math.Log10(Math.Pow(2, headerobject.nbits - 1) / headerobject.calib_l / dbref / 2 / Math.Sqrt(2));
            if (headerobject.calib_r != -1) headerobject.calib_r = 20 * Math.Log10(Math.Pow(2, headerobject.nbits - 1) / headerobject.calib_r / dbref / 2 / Math.Sqrt(2));
            headerread = true;

            //if (!datacollected) return dummy;
            var propnames = new string[15] {    "Calib_l",
                                                "Calib_r",
                                                "nbits",
                                                "nchannels",
                                                "Description_left",
                                                "Description_right",
                                                "quantity_l",
                                                "quantity_r",
                                                "samplerate",
                                                "speedpulses_l",
                                                "speedpulses_r",
                                                "speedfactor_l",
                                                "speedfactor_r",
                                                "unit_l",
                                                "unit_r" };
            var props = new string[15] {        headerobject.calib_l.ToString(),
                                                headerobject.calib_r.ToString(),
                                                headerobject.nbits.ToString(),
                                                headerobject.nchannels.ToString(),
                                                headerobject.description_l,
                                                headerobject.description_r,
                                                headerobject.quantity_l,
                                                headerobject.quantity_r,
                                                headerobject.samplerate.ToString(),
                                                headerobject.speedpulses_l.ToString(),
                                                headerobject.speedpulses_r.ToString(),
                                                headerobject.speedfactor_l.ToString(),
                                                headerobject.speedfactor_r.ToString(),
                                                headerobject.unit_l,
                                                headerobject.unit_r};

            var result_x = new IDataSet(DataSetType.String)
            {
                Count = propnames.Length,
                StringValues = propnames
            };

            var result = new IDataSet(DataSetType.String)
            {
                Count = props.Length,
                StringValues = props,
                x = result_x
            };
            return result;
        }

    }

    public class HSI
    {
        /*Test of the HSI algorithm as workaround if cycle addressing does not work... Works within Concerto, Problems in Indicom*/
        NVHFunctions.CycleData pcyl;
        int cylcount;
        double[] tdcoffset, cyctime, xorg, buffer;
        int[] pinofftype;
        double bore, conrod, stroke, pinoff, pinoff2;

        public void Init(IDataSet pcyl1, IDataSet cyctime_dec, IDataSet tdcoffset_dec, IDataSet pinofftype_dec, double bore_dec, double conrod_dec, double stroke_dec, double pinoff_dec, double pinoff2_dec)
        {
            cylcount = tdcoffset_dec.Count;

            tdcoffset = tdcoffset_dec.RealValues;
            var pinofftype_d = pinofftype_dec.RealValues;
            pinofftype = new int[pinofftype_d.Length];
            for (int ii = 0; ii < pinofftype.Length; ii++) pinofftype[ii] = (int)pinofftype_d[ii];
            cyctime = cyctime_dec.RealValues;

            pcyl = new NVHFunctions.CycleData
            {
                Count = 720,
                CycleCount = cyctime.Length,
                xdata = new double[720],
                ydata = new double[720 * cyctime.Length * cylcount]
            };

            xorg = pcyl1.x.RealValues;
            buffer = new double[pcyl1.Count];
            
            int kk = 0;
            for (int ii = 0; ii < pcyl1.Count; ii++)
            {
                if (Math.Round(xorg[ii], 3, MidpointRounding.AwayFromZero) - Math.Round(xorg[ii], MidpointRounding.AwayFromZero) == 0)
                {
                    pcyl.xdata[kk] = xorg[ii];
                    kk++;
                }
            }
            bore = bore_dec;
            conrod = conrod_dec;
            stroke = stroke_dec;
            pinoff = pinoff_dec;
            pinoff2 = pinoff2_dec;


        }

        public void CollectPCYLData(IDataSet pcyl1, IDataSet pcyl2, IDataSet pcyl3, IDataSet pcyl4, IDataSet pcyl5, IDataSet pcyl6, double CycleIndex)
        {

            var jarray = new double[cylcount][];
            jarray[0] = pcyl1.RealValues;
            if (cylcount > 1) jarray[1] = pcyl2.RealValues;
            if (cylcount > 2) jarray[2] = pcyl3.RealValues;
            if (cylcount > 3) jarray[3] = pcyl4.RealValues;
            if (cylcount > 4) jarray[4] = pcyl5.RealValues;
            if (cylcount > 5) jarray[5] = pcyl6.RealValues;

            
            for (int cyl = 0; cyl < cylcount; cyl++)
            {
                int kk = (int)CycleIndex * 720 + cyl * 720 * cyctime.Length;
                buffer = jarray[cyl];
                for (int ii = 0; ii < pcyl1.Count; ii++)
                {
                    if (Math.Round(xorg[ii], 3, MidpointRounding.AwayFromZero) - Math.Round(xorg[ii], MidpointRounding.AwayFromZero) == 0)
                    {
                        pcyl.ydata[kk] = buffer[ii];
                        kk++;
                    }
                }
                
            }
        }

        public IDataSet Calculate()
        {
            

            var result = GetHSI(pcyl, tdcoffset, pinofftype, bore, conrod, stroke, pinoff, pinoff2, cyctime);
            
            var HSIResult_x = new IDataSet(DataSetType.Numeric)
            {
                Count = result.xdata.Length,
                RealValues = result.xdata
            };

            var HSIResult = new IDataSet(DataSetType.Numeric)
            {
                Count = result.xdata.Length,
                x = HSIResult_x,
                RealValues = result.ydata
            };
            
            return HSIResult;
        }

    }

    public class OnlinePlayback_Concerto
    {
        //private static readonly System.Diagnostics.Stopwatch timer = new System.Diagnostics.Stopwatch();
        private static NAudio.Wave.BufferedWaveProvider bufferedWaveProvider;
        private static NAudio.Wave.WaveOut waveOut;
        private static readonly System.Timers.Timer timer = new System.Timers.Timer(500);
        private static readonly System.Diagnostics.Stopwatch stopwatch = new System.Diagnostics.Stopwatch();
        private static readonly int nbytes = 2;
        private static int samplerate = 50000;
        private static IDataSet playdata;
        private static int blockcount = 0;
        private static int buffersize = 500;
        private static bool dataisinms = true;
        private static double calib = 1;
        private static bool timer_initialized = false;
        private static double startpos_ms = 0;
        private static double playuntil = double.NaN;


        public int PlayStatus()
        {
            if (playdata == null) return -1;
            if (timer.Enabled) return 1;
            else if (!timer.Enabled && blockcount != 0) return 2;
            else return 0;
        }

        public int PositionInfoIsValid()
        {
            if (waveOut != null)
            {
                try
                {
                    var test = waveOut.GetPosition();
                    return 1;
                }
                catch (Exception e)
                {
                    return 0;
                }
            }
            else return 0;
        }

        public double GetPosition()
        {
            //var position = blockcount * buffersize;
            //var position = stopwatch.ElapsedMilliseconds + startpos_ms;
            long streampos;
            double result;
            if (waveOut != null)
            {
                try
                {
                    streampos = waveOut.GetPosition();
                    result = (double)streampos * 1000 / (samplerate * nbytes);
                }
                catch (Exception e)
                {
                    return -1;
                }
                
            }
            else result = -1;
            return result + startpos_ms;
        }

        public void SetDataSet(IDataSet dataset, int ms, double calibrationfactor)
        {
            playdata = dataset;
            if (ms == 1) dataisinms = true;
            else if (ms == 0) dataisinms = false;
            calib = calibrationfactor;
            startpos_ms = dataset.x[0];
            if (!dataisinms) startpos_ms *= 1000;
            

        }

        // Testfunktion, um das statische Verbleiben des Datensets zu testen. Wird nicht mehr benötigt.
        /*public int GetDataCount()
        {
            return playdata.Count;
        }*/


        private void DisposeAudio()
        {
            if (bufferedWaveProvider != null) bufferedWaveProvider.ClearBuffer();
            
            if (waveOut != null)
            {
                waveOut.Stop();
                waveOut.Dispose();
            }
            

        }
        /*public void Reset()
        {
            
            if (bufferedWaveProvider != null)
            {
                if (bufferedWaveProvider.BufferedBytes > 0) bufferedWaveProvider.ClearBuffer();
            }

        }*/

        private void InitAudio()
        {
            bufferedWaveProvider = new NAudio.Wave.BufferedWaveProvider(new NAudio.Wave.WaveFormat(samplerate, nbytes * 8, 1))
            {
                BufferLength = samplerate * nbytes * 10,
                DiscardOnBufferOverflow = false
            };

            waveOut = new NAudio.Wave.WaveOut();
            waveOut.Init(bufferedWaveProvider);
            waveOut.Play();
            
        }

        public void StartPlayback(double position, int buffersize_ms = 500, double playtime = 0, int playtimelimit = 0)
        {
            if (playdata == null) return;
            if (playtimelimit != 0) playuntil = position + playtime;
            else playuntil = double.NaN;
            buffersize = buffersize_ms;
            if (timer.Interval != buffersize) timer.Interval = buffersize;
            if (!timer.Enabled)
            {
                stopwatch.Start();
                startpos_ms = position;
                double x0_ms = playdata.x[0];
                
                if (!dataisinms)
                {
                    x0_ms *= 1000;
                    startpos_ms *= 1000;
                }

                blockcount = (int)((startpos_ms - x0_ms) / buffersize);
                //AddSamples();
                timer.Start();
                
            }
            if (!timer_initialized)
            {
                timer.Elapsed += Timer_Elapsed;
                timer_initialized = true;
            }
            
        }


        public void StopPlayback()
        {
            timer.Stop();
            blockcount = 0;
            stopwatch.Reset();
            DisposeAudio();
        }

        /*public void PausePlayback()
        {
            timer.Stop();
            stopwatch.Stop();
        }*/
        

        private void Timer_Elapsed(object sender, System.Timers.ElapsedEventArgs e)
        {
            AddSamples();
            if (bufferedWaveProvider.BufferedBytes == 0) StopPlayback();
        }

        public void AddSamples()
        {
            try
            {
                //string message = "";
                //var newaudio = data.RealValues;
                //timestamp_onlineplayback = data.x[data.Count - 1];
                
                int sr;
                if (dataisinms) sr = (int)Math.Round(1000 * (double)playdata.Count / (playdata.x[playdata.Count - 1] - playdata.x[0]));
                else sr = (int)Math.Round((double)playdata.Count / (playdata.x[playdata.Count - 1] - playdata.x[0]));
                if (sr != samplerate)
                {
                    if (waveOut != null) DisposeAudio();
                    samplerate = sr;
                    InitAudio();
                }

                //data.Release();

                //bool wpInitialized = false;

                if (bufferedWaveProvider == null)
                {
                    bufferedWaveProvider = new NAudio.Wave.BufferedWaveProvider(new NAudio.Wave.WaveFormat(samplerate, nbytes * 8, 1))
                    {
                        BufferLength = samplerate * nbytes * 10,
                        DiscardOnBufferOverflow = false
                    };

                    //  wpInitialized = true;
                }
                if (waveOut == null)
                {
                    waveOut = new NAudio.Wave.WaveOut();
                    //wpInitialized = true;
                }
                if (waveOut.PlaybackState != NAudio.Wave.PlaybackState.Playing)
                {
                    waveOut.Init(bufferedWaveProvider);
                    waveOut.Play();
                    //wpInitialized = true;
                }
                //if (wpInitialized) message += "Online Playback engine started.";
                int L_block = (int)((double)buffersize / 1000 * samplerate);
                if (blockcount == 0) L_block *= 5;
                //var buffer = new byte[newaudio.Length * nbytes];
                var buffer = new byte[(int)(L_block * nbytes)];
                bool endofstream = false;

                if (blockcount * L_block > playdata.Count - 1) return;
                for (int ii = 0; ii < L_block; ii++)
                {
                    int index_audio = blockcount * L_block + ii;
                    if (index_audio > playdata.Count - 1)
                    {
                        endofstream = true;
                        //blockcount = -1;
                        //timer.Stop();
                        break;
                    }
                    else if(!double.IsNaN(playuntil))
                    {
                        if (index_audio > (playuntil - playdata.x[0]) * (samplerate / 1000))
                        {
                            endofstream = true;
                            break;
                        }
                    }
                    var bytes = BitConverter.GetBytes((int)(playdata[index_audio] * calib * Math.Pow(2, nbytes * 8 - 1)));
                    for (int jj = 0; jj < nbytes; jj++) buffer[ii * nbytes + jj] = bytes[jj];
                }
                if (!endofstream) bufferedWaveProvider.AddSamples(buffer, 0, buffer.Length);
                //if (blockcount == -1) DisposeAudio();
                
                if (blockcount == 0)
                {
                    
                    blockcount += 5;
                }
                else blockcount++;
                
            }
            catch (Exception)
            {
                StopPlayback();
                //timer.Stop();
                //DisposeAudio();
            }
        }


        /*public void PlayTestTone()
        {
            if (50000 != samplerate)
            {
                DisposeAudio();
                samplerate = 50000;
                InitAudio();
            }

            //bool wpInitialized = false;
            if (bufferedWaveProvider == null)
            {
                bufferedWaveProvider = new NAudio.Wave.BufferedWaveProvider(new NAudio.Wave.WaveFormat(samplerate, nbytes * 8, 1))
                {
                    BufferLength = samplerate * nbytes * 10,
                    DiscardOnBufferOverflow = false
                };

                //wpInitialized = true;
            }
            if (waveOut == null)
            {
                waveOut = new NAudio.Wave.WaveOut();
                //wpInitialized = true;
            }
            if (waveOut.PlaybackState != NAudio.Wave.PlaybackState.Playing)
            {
                waveOut.Init(bufferedWaveProvider);
                waveOut.Play();
                //wpInitialized = true;
            }
            //if (wpInitialized) message += "Online Playback engine started.";

            var sig = SignalGeneration.SineWave(1, 50000, 440, 0.1);

            var buffer = new byte[sig.Length * nbytes];

            for (int ii = 0; ii < sig.Length; ii++)
            {
                var bytes = BitConverter.GetBytes((int)(sig[ii] * Math.Pow(2, nbytes * 8 - 1)));
                for (int jj = 0; jj < nbytes; jj++) buffer[ii * nbytes + jj] = bytes[jj];
            }
            bufferedWaveProvider.AddSamples(buffer, 0, buffer.Length);



        }

        public IDataSet GetTestTone()
        {
            var sig = SignalGeneration.SineWave(1, 50000, 10, 0.5);
            var sig_x = new double[sig.Length];
            for (int ii = 0; ii < sig_x.Length; ii++) sig_x[ii] = (double)ii / 50000;
            var x = new IDataSet(DataSetType.Numeric)
            {
                Count = sig_x.Length,
                RealValues = sig_x
            };
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = sig.Length,
                RealValues = sig,
                x = x
            };
            return result;
        }*/
    }
    
    public class OnlinePlayback
    {

        private static NAudio.Wave.BufferedWaveProvider bufferedWaveProvider;
        private static double timestamp_onlineplayback = 0;     // dies ist der Zeitstempel des Datensets. Für Online-Datensets jeweils der letzte Datenpunkt.
        private static int blockcount_onlineplayback = 0;       // so viele Blöcke wurden bereits in den Datenpuffer übergeben.
        private static readonly System.Diagnostics.Stopwatch timer = new System.Diagnostics.Stopwatch();
        private static NAudio.Wave.WaveOut waveOut;
        private readonly int nbytes = 2;
        private int samplerate = 50000;
        private static int dynamicrange = 120;
        

        public double Get_timestamp_concerto()
        {
            return timer.ElapsedMilliseconds;
        }


        public double Get_timestamp()
        {
            return timestamp_onlineplayback;
        }
        public double Get_blockcount()
        {
            return blockcount_onlineplayback;
        }


        public string AddSamples(IDataSet data)
        {
            try
            {
                string message = "";
                var newaudio = data.RealValues;
                timestamp_onlineplayback = data.x[data.Count - 1];

                var sr = (int)Math.Round(1000 / (data.x[1] - data.x[0]));
                if (sr != samplerate)
                {
                    DisposeAudio();
                    samplerate = sr;
                    InitAudio();
                }

                data.Release();

                bool wpInitialized = false;
                
                if (bufferedWaveProvider == null)
                {
                    bufferedWaveProvider = new NAudio.Wave.BufferedWaveProvider(new NAudio.Wave.WaveFormat(samplerate, nbytes * 8, 1))
                    {
                        BufferLength = samplerate * nbytes * 10,
                        DiscardOnBufferOverflow = false
                    };

                    wpInitialized = true;
                }
                if (waveOut == null)
                {
                    waveOut = new NAudio.Wave.WaveOut();
                    wpInitialized = true;
                }
                if (waveOut.PlaybackState != NAudio.Wave.PlaybackState.Playing)
                {
                    waveOut.Init(bufferedWaveProvider);
                    waveOut.Play();
                    wpInitialized = true;
                }
                if (wpInitialized) message += "Online Playback engine started.";

                var buffer = new byte[newaudio.Length * nbytes];

                for (int ii = 0; ii < newaudio.Length; ii++)
                {
                    var bytes = BitConverter.GetBytes((int)(newaudio[ii] * Math.Pow(2, nbytes * 8 - 1)));
                    for (int jj = 0; jj < nbytes; jj++) buffer[ii * nbytes + jj] = bytes[jj];
                }
                bufferedWaveProvider.AddSamples(buffer, 0, buffer.Length);
                blockcount_onlineplayback++;
                return message;
                
            }
            catch (Exception e)
            {
                data.Release();
                return e.Message;

            }
        }
        
        public void Reset()
        {
            timer.Restart();
            blockcount_onlineplayback = 0;
            timestamp_onlineplayback = 0;
            if (bufferedWaveProvider != null)
            {
                if (bufferedWaveProvider.BufferedBytes > 0) bufferedWaveProvider.ClearBuffer();
            }
        }

        public void SetDynamicRange(int dynrange)
        {
            dynamicrange = dynrange;
        }

        public int GetDynamicRange()
        {
            return dynamicrange;
        }

        public void InitAudio()
        {
            bufferedWaveProvider = new NAudio.Wave.BufferedWaveProvider(new NAudio.Wave.WaveFormat(samplerate, nbytes * 8, 1))
            {
                BufferLength = samplerate * nbytes * 10,
                DiscardOnBufferOverflow = false
            };
            waveOut = new NAudio.Wave.WaveOut();
            waveOut.Init(bufferedWaveProvider);
            waveOut.Play();
        }


        public void DisposeAudio()
        {
            waveOut.Stop();
            waveOut.Dispose();

        }
    }


    public class TrendRecorder
    {
        
    }
    

    public class ConcertoWrapper
    {
        private static readonly bool allowcommonxfeature = true;

        //XYZ_Concerto_Data apsresult;
        private static double[] apsresult_xdata;
        private static double[] apsresult_ydata;
        private static double[] apsresult_zdata;
        private static double[] apsresult_commonxdata;
        private static double[] apsresult_dimensions;
        private static double[] apsresult_xtime;
        private static double[] apsresult_freq_labels;
        private static double apsresult_overalllevel;

        double[] buffer1;
        private static NAudio.Wave.BufferedWaveProvider bufferedWaveProvider;
        NVHFunctions.CycleData cadatainput;
        NVHFunctions.CycleData[] cadatainputarray;
        string errormessage = "ok";
        private static long playbacklength_bytes;
        NAudio.Wave.IWaveProvider provider;
        readonly System.Diagnostics.Stopwatch timer = new System.Diagnostics.Stopwatch();
        NVHFunctions.CATimeStamps time_stamps;
        private static double timestamp_onlineplayback;
        private static int blockcount_onlineplayback;
        private static double aps_timestamp = 0;
        private static NAudio.Wave.WaveOut waveOut;
        
        //private static IDataSet[] xyzdata = new IDataSet[3];

        public void ResetTime()
        {
            aps_timestamp = 0;
        }

        /*public static double[] ResampleTime2AngleDomain(double[] audiodata, double x0, int fs, XY_Data rpm, string ConcertoWorkDir)
        {
            // Apply Anti Aliasing Filter to prepare for resampling

        }*/


        public static double[] LP10_butterworth_RPMFrequency(double[] audiodata, int fs, XY_Data rpm, double samplesperrev, string pdpatchlocation)
        {
            double rpmres = (rpm.xdata[rpm.xdata.Length - 1] - rpm.xdata[0]) / rpm.xdata.Length;
            var filterfreq = rpm.ydata;
            if (rpmres != 0.125)
            {
                var rpm_x = DeltaArray(0, 0.125, rpm.xdata[rpm.xdata.Length - 1] - rpm.xdata[0]);
                filterfreq = Interp1_linear(rpm.xdata, rpm.ydata, rpm_x, true);
            }

            double factor = samplesperrev / 2.5 / 60;
            for (int ii = 0; ii < filterfreq.Length; ii++) filterfreq[ii] *= factor;
            var filtered_data = PdWrapper.LaunchPd(pdpatchlocation.Replace('\\', '/'), audiodata, fs, filterfreq, 0.125);
            return filtered_data;

        }

        public IDataSet CalcHSI_v2(IDataSet pcyl1, IDataSet pcyl2, IDataSet pcyl3, IDataSet pcyl4, IDataSet pcyl5, IDataSet pcyl6, IDataSet cyctime_dec, IDataSet tdcoffset_dec, IDataSet pinofftype_dec, double bore, double conrod, double stroke, double pinoff, double pinoff2)
        {
            /* Final version given to HMC */

            // Instantiate dummy result

            //var IOClass = new DataIOFunctions();

            int cylcount = tdcoffset_dec.Count;

            var tdcoffset = tdcoffset_dec.RealValues;
            var pinofftype_d = pinofftype_dec.RealValues;
            var pinofftype = new int[pinofftype_d.Length];
            for (int ii = 0; ii < pinofftype.Length; ii++) pinofftype[ii] = (int)pinofftype_d[ii];
            var cyctime = cyctime_dec.RealValues;


            var jarray = new IDataSet[cylcount];
            jarray[0] = pcyl1;
            if (cylcount > 1) jarray[1] = pcyl2;
            if (cylcount > 2) jarray[2] = pcyl3;
            if (cylcount > 3) jarray[3] = pcyl4;
            if (cylcount > 4) jarray[4] = pcyl5;
            if (cylcount > 5) jarray[5] = pcyl6;

            var pcyl = new NVHFunctions.CycleData
            {
                Count = 720,
                CycleCount = cyctime.Length,
                xdata = new double[720],
                ydata = new double[720 * cyctime.Length * cylcount]
            };

            var xorg = pcyl1.x.RealValues;
            double[] buffer;

            int kk = 0;
            for (int cyl = 0; cyl < cylcount; cyl++)
            {
                for (int c = 0; c < pcyl.CycleCount; c++)
                {
                    jarray[cyl].CycleIndex = c;
                    buffer = jarray[cyl].RealValues;
                    for (int ii = 0; ii < pcyl1.Count; ii++)
                    {
                        if (Math.Round(xorg[ii], 3, MidpointRounding.AwayFromZero) - Math.Round(xorg[ii], MidpointRounding.AwayFromZero) == 0)
                        {
                            if (cyl == 0 && c == 0) pcyl.xdata[kk] = xorg[ii];
                            pcyl.ydata[kk] = buffer[ii];
                            kk++;
                        }
                    }
                }
            }


            pcyl1.Release();
            pcyl2.Release();
            pcyl3.Release();
            pcyl4.Release();
            pcyl5.Release();
            pcyl6.Release();

            var result = GetHSI(pcyl, tdcoffset, pinofftype, bore, conrod, stroke, pinoff, pinoff2, cyctime);

            var HSIResult_x = new IDataSet(DataSetType.Numeric)
            {
                Count = result.xdata.Length,
                RealValues = result.xdata
            };

            var HSIResult = new IDataSet(DataSetType.Numeric)
            {
                Count = result.xdata.Length,
                x = HSIResult_x,
                RealValues = result.ydata
            };


            return HSIResult;
        }

        

        internal static XYZ_Data CalculateDiagram(double[] rawdataY, double x0, double dt, IDataSet dT, IDataSet cyctime, IDataSet n100_dec, IDataSet enginespeed_dec, Param_struct paramset, double time_factor)
        {
            DataIOFunctions DataIOClass = new DataIOFunctions();
            //XYZ_Concerto_Data apsresult;

            int ms = 0;
            if (time_factor == 0.001) ms = 1;

            XY_Data enginespeed;
            XY_Data n100;
            var timestamps = new CATimeStamps
            {
                t1 = new double[1],
                t2 = new double[1],
                cycles = new double[1]

            };
            /*double[] rawdataY;

            double x0;
            double dt;*/


            CycleData cdmtime;
            

            // Calculate Time Stamps
            if (paramset.LQ == 2)
            {
                if (paramset.CDMTIME)
                {
                    cdmtime = DataIOClass.ReadConcertoCAData(dT, paramset.Cyclerange1, paramset.Cyclerange2, (int)paramset.Delta_t);
                }
                else
                {
                    var dT_cyc = DataIOClass.ReadConcertoCAData(dT, paramset.Cyclerange1, paramset.Cyclerange2, (int)paramset.Delta_t);
                    var cyctime_d = new double[dT_cyc.CycleCount];
                    int kk = 0;
                    for (int ii = paramset.Cyclerange1; ii <= paramset.Cyclerange2; ii += (int)paramset.Delta_t)
                    {
                        cyctime_d[kk] = cyctime.RealValues[ii];
                        kk++;
                    }
                    cdmtime = DT2CDMTIME(dT_cyc, cyctime_d, ms);
                }



                timestamps = GetCATimeStamps_CDMTIME(cdmtime, paramset.FromCA, paramset.ToCA, ms);
                for (int ii = 0; ii < timestamps.cycles.Length; ii++) timestamps.cycles[ii] = ii * (int)paramset.Delta_t + paramset.Cyclerange1 + 1;

            }



            /*if (rawdata_dec.CycleCount > 1)
            {
                // Resample CA Data
                // Do special "fast forward" calculation for cyclic data and delta > 2
                if (paramset.LQ == 2 && paramset.Delta_t > 2)
                {

                    // Collect data of cycles to calcuate into an array
                    enginespeed = DataIOClass.ReadConcertoData(enginespeed_dec);
                    n100_dec.Release();
                    enginespeed_dec.Release();
                    for (int ii = 0; ii < enginespeed.xdata.Length; ii++) enginespeed.xdata[ii] *= time_factor;

                    //// Prepare Arrays for cyclewise Calculation
                    CycleData[] rawdata_array;
                    if (!formulacalc) rawdata_array = ReadCADataIntoCycleDataArray(rawdata_dec, timestamps);
                    else rawdata_array = cadatainputarray;
                    var dT_array = ReadCADataIntoCycleDataArray(dT, timestamps);


                    //// Do calculation for each cycle to calculate
                    var cyctime_s = cyctime.RealValues;
                    for (int ii = 0; ii < cyctime_s.Length; ii++) cyctime_s[ii] *= time_factor;
                    cyctime.Release();
                    apsresult = CalculateCycleWise(rawdata_array, enginespeed, dT_array, cyctime_s, timestamps, paramset);
                    var apsresult_d = ConvertResultBack(apsresult);
                    return apsresult_d;
                }
                else
                // Resample whole CA dataset if delta < 3 or if leading quantity is not cycle.
                {
                    int add = 0;
                    if (paramset.Cyclerange2 < rawdata_dec.CycleCount - 1) add = 1;
                    CycleData raw_d;
                    if (!formulacalc) raw_d = DataIOClass.ReadConcertoCAData(rawdata_dec, paramset.Cyclerange1, paramset.Cyclerange2 + add);
                    else raw_d = cadatainput;
                    rawdata_dec.Release();
                    var dT_d = DataIOClass.ReadConcertoCAData(dT, paramset.Cyclerange1, paramset.Cyclerange2 + add);
                    dT.Release();

                    if (paramset.CDMTIME) rawdataY = ResampleCAData_CDMTIME(raw_d, dT_d, paramset.Fs);
                    else rawdataY = ResampleCAData_dT(raw_d, dT_d, paramset.Fs);

                    x0 = (cyctime.RealValues[paramset.Cyclerange1]) * time_factor;
                    dt = 1 / paramset.Fs;
                }
            }
            else
            {
                rawdataY = rawdata_dec.RealValues;
                x0 = rawdata_dec.x[0] * time_factor;
                dt = (rawdata_dec.x[1] - rawdata_dec.x[0]) * time_factor;
            }*/


            cyctime.Release();

            //int dsrange1_index = (int)((paramset.Dsrange1 - x0) / dt);
            //int dsrange2_index = (int)((paramset.Dsrange2 - x0) / dt);

            // Resample TM data if necessary
            if (Math.Round(paramset.Fs) != Math.Round(1 / dt))
            {
                /*if (paramset.Dsrange2 - paramset.Dsrange1 < 0.7 * (rawdata_dec.x[rawdata_dec.Count - 1] * time_factor - x0))
                {
                    int add = 0;
                    if (dsrange2_index + 23 < rawdata_dec.Count) add = 23;
                    var rawdataY_red = new double[dsrange2_index - dsrange1_index + 1 + add];
                    for (int ii = 0; ii < rawdataY_red.Length; ii++) rawdataY_red[ii] = rawdataY[dsrange1_index + ii];
                    x0 = paramset.Dsrange1;
                    rawdataY = ResampleTMData(rawdataY_red, x0, dt, paramset.Fs);
                }
                else*/ rawdataY = ResampleTMData(rawdataY, x0, dt, paramset.Fs);
                dt = 1 / paramset.Fs;
            }
            //rawdata_dec.Release();

            if (paramset.Dsrange1 < x0) paramset.Dsrange1 = x0;
            if (paramset.Dsrange2 > x0 + dt * (rawdataY.Length - 1)) paramset.Dsrange2 = x0 + dt * (rawdataY.Length - 1);

            enginespeed = DataIOClass.ReadConcertoData(enginespeed_dec);
            enginespeed_dec.Release();
            n100 = DataIOClass.ReadConcertoData(n100_dec);
            n100_dec.Release();
            for (int ii = 0; ii < enginespeed.xdata.Length; ii++) enginespeed.xdata[ii] *= time_factor;
            for (int ii = 0; ii < n100.xdata.Length; ii++) n100.xdata[ii] *= time_factor;



            XYZ_Data apsresult_double;
            switch (paramset.Diagramtype)
            {
                case 1:

                    /*var output = new NVHFunctions.NVHPackage
                    {
                        rawdata_y = rawdataY,
                        n100 = n100,
                        enginespeed = enginespeed,
                        dt = dt,
                        x0 = x0,
                        timestamps = timestamps,
                        paramset = paramset
                    };
                    APSClass.SaveAPSInput("C:/Users/u12o24/Documents/nvhpackage1.bin", output);
                    APSClass.SaveSomeDoubles("C:/Users/u12o24/Documents/10kHz.bin",rawdataY);*/


                    apsresult_double = CalcAPSData(rawdataY, n100, x0, timestamps, paramset);
                    break;
                case 2:
                    /*var output = new NVHFunctions.NVHPackage
                    {
                        rawdata_y = rawdataY,
                        n100 = n100,
                        enginespeed = enginespeed,
                        dt = dt,
                        x0 = x0,
                        timestamps = timestamps,
                        paramset = paramset
                    };
                    APSClass.SaveAPSInput("C:/Users/u12o24/Documents/nvhpackage1.bin", output);*/

                    apsresult_double = CalcNtelOctaveSpectra(rawdataY, n100, dt, x0, timestamps, paramset);
                    break;
                case 3:
                    apsresult_double = CalcOrderSpectra_FFT(rawdataY, n100, enginespeed, x0, timestamps, paramset);
                    break;
                case 4:
                    apsresult_double = CalcOrderSpectra_FFT(rawdataY, n100, enginespeed, x0, timestamps, paramset);
                    break;
                case 5:
                    apsresult_double = CalcOverallLevel(rawdataY, n100, x0, timestamps, paramset);
                    break;
                case 6:
                    apsresult_double = CalcBPLevel(rawdataY, n100, x0, timestamps, paramset);
                    break;
                default:
                    apsresult_double = CalcAPSData(rawdataY, n100, x0, timestamps, paramset);
                    break;
            }
            //if (apsresult_double.errorstring != null) errormessage = apsresult_double.errorstring;

            //apsresult = ConvertResult(apsresult_double);
            //rawdata_dec.Release();

            return apsresult_double;
        }


        /// <summary>
        /// Calculates crank angle based data cycle-wise. Only for CA Data and LQ = 2 (Cycle). Calculation is much faster, because time-costly Resampling is only done for the cycles which need to be calculated.
        /// </summary>
        /// <param name="rawdata_array">Array of CycleData. Includes only cycles which are to be calculated.</param>
        /// <param name="enginespeed">Engine speed for e.g. order spectrum.</param>
        /// <param name="dT_array">Array of dT or CDMTIME signal. Includes only cycles which are to be calculated.</param>
        /// <param name="cyctime">Start time of each cycle. Includes every cycle. Must be in seconds.</param>
        /// <param name="timestamps">Time stamps for CA analysis. Includes only cycles which are to be calculated.</param>
        /// <param name="paramset">Set of FFT parameters.</param>
        /// <returns>2D/3D diagram result in Concerto-style.</returns>
        internal static XYZ_Concerto_Data CalculateCycleWise(CycleData[] rawdata_array, XY_Data enginespeed, CycleData[] dT_array, double[] cyctime, CATimeStamps timestamps, Param_struct paramset)
        {
            var result_array = new XYZ_Data[timestamps.cycles.Length];
            double dt = 1 / paramset.Fs;
            Parallel.For(0, timestamps.cycles.Length, i =>
            //for (int i = 0; i < timestamps.cycles.Length; i++)
            {
                int actcycle = (int)timestamps.cycles[i];
                double[] rawdataonecycle;
                if (paramset.CDMTIME) rawdataonecycle = ResampleCAData_CDMTIME(rawdata_array[i], dT_array[i], paramset.Fs);
                else rawdataonecycle = ResampleCAData_dT(rawdata_array[i], dT_array[i], paramset.Fs);
                var paramstruct = paramset;
                double x0_ = cyctime[actcycle - 1];
                if (paramstruct.Dsrange1 < x0_) paramstruct.Dsrange1 = x0_;
                if (paramstruct.Dsrange2 > x0_ + dt * (rawdataonecycle.Length - 1)) paramstruct.Dsrange2 = x0_ + dt * (rawdataonecycle.Length - 1);
                var timestamps_ = new NVHFunctions.CATimeStamps
                {
                    t1 = new double[1] { timestamps.t1[i] },
                    t2 = new double[1] { timestamps.t2[i] },
                    cycles = new double[1] { actcycle }
                };

                switch (paramset.Diagramtype)
                {
                    case 1:
                        result_array[i] = CalcAPSData(rawdataonecycle, new XY_Data(), x0_, timestamps_, paramstruct);
                        break;
                    case 2:
                        result_array[i] = CalcNtelOctaveSpectra(rawdataonecycle, new XY_Data(), dt, x0_, timestamps_, paramstruct);
                        break;
                    case 3:
                        result_array[i] = CalcOrderSpectra_FFT(rawdataonecycle, new XY_Data(), enginespeed, x0_, timestamps_, paramstruct);
                        break;
                    case 4:
                        result_array[i] = CalcOrderSpectra_FFT(rawdataonecycle, new XY_Data(), enginespeed, x0_, timestamps_, paramstruct);
                        break;
                    case 5:
                        result_array[i] = CalcOverallLevel(rawdataonecycle, new XY_Data(), x0_, timestamps_, paramstruct);
                        break;
                    case 6:
                        result_array[i] = CalcBPLevel(rawdataonecycle, new XY_Data(), x0_, timestamps_, paramstruct);
                        break;
                    default:
                        result_array[i] = CalcAPSData(rawdataonecycle, new XY_Data(), x0_, timestamps_, paramstruct);
                        break;
                }
                //if (result_array[i].errorstring != null) errormessage = result_array[i].errorstring;
            }
            );
            
            return ConvertResult_array(result_array, paramset);
        }

        

        public void CalculateTimeStamps_CDMTIME(IDataSet CDMTIME, int cyclerange1, int cyclerange2, int delta_t, double fromCA, double toCA, int ms)
        {
            var IOFcnClass = new DataIOFunctions();


            CDMTIME.CycleIndex = cyclerange1 - 1;
            var cdmtime_cyc = IOFcnClass.ReadConcertoCAData(CDMTIME, cyclerange1 - 1, cyclerange2 - 1, delta_t);
            CDMTIME.Release();
            time_stamps = GetCATimeStamps_CDMTIME(cdmtime_cyc, fromCA, toCA, ms);
            for (int ii = 0; ii < time_stamps.cycles.Length; ii++) time_stamps.cycles[ii] = ii * delta_t + cyclerange1;
        }
        public void CalculateTimeStamps_dT(IDataSet dT, IDataSet cyctime, int cyclerange1, int cyclerange2, int delta_t, double fromCA, double toCA, int ms)
        {
            var IOFcnClass = new DataIOFunctions();


            var dT_cyc = IOFcnClass.ReadConcertoCAData(dT, cyclerange1 - 1, cyclerange2 - 1, delta_t);
            dT.Release();
            var cyctime_d = new double[dT_cyc.CycleCount];
            int kk = 0;
            for (int ii = cyclerange1 - 1; ii <= cyclerange2 - 1; ii += delta_t)
            {
                cyctime_d[kk] = cyctime.RealValues[ii];
                kk++;
            }
            cyctime.Release();

            var cdmtime = DT2CDMTIME(dT_cyc, cyctime_d, ms);
            time_stamps = GetCATimeStamps_CDMTIME(cdmtime, fromCA, toCA, ms);
            for (int ii = 0; ii < time_stamps.cycles.Length; ii++) time_stamps.cycles[ii] = ii * delta_t + cyclerange1;

        }

        public IDataSet CDMTIME2TMRPM(IDataSet CDMTIME, IDataSet cyctime, int cyclerange1, int cyclerange2)
        {
            var IOFcnClass = new DataIOFunctions();

            var cdmtime_cyc = IOFcnClass.ReadConcertoCAData(CDMTIME, cyclerange1 - 1, cyclerange2 - 1);
            CDMTIME.Release();

            var cycrpm = CDMTIME2CycRPM(cdmtime_cyc);

            var tmrpm = new IDataSet(DataSetType.Numeric)
            {
                Count = cycrpm.Length,
                x = cyctime,
                RealValues = cycrpm
            };
            cyctime.Release();
            return tmrpm;

        }

        public IDataSet GetTimeStampsFromTimeAxis(IDataSet timeaxis, double fs, double deltaf, int average, int average_overlap)
        {
            var timestamps = GetTimeStampsFromXAxis(timeaxis.RealValues, fs, deltaf, average, average_overlap);
            var result_d = new double[timestamps.t1.Length * 2];
            for (int ii = 0; ii < timestamps.t1.Length; ii++)
            {
                result_d[ii * 2] = timestamps.t1[ii];
                result_d[ii * 2 + 1] = timestamps.t2[ii];
            }
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = result_d.Length,
                RealValues = result_d
            };
            return result;
        }

        internal static XYZ_Data ConvertResultBack(XYZ_Concerto_Data result_Concerto)
        {
            XYZ_Data result = new XYZ_Data()
            {
                xdata = result_Concerto.xdata.RealValues,
                ydata = result_Concerto.ydata.RealValues,
                zdata = result_Concerto.zdata.RealValues,
                xtime = result_Concerto.xtime.RealValues,
                dimensions = new int[2] { (int)result_Concerto.dimensions.RealValues[0], (int)result_Concerto.dimensions.RealValues[1] },
                overalllevel = result_Concerto.overalllevel
            };
            return result;
        }

        internal static XYZ_Concerto_Data ConvertResult(XYZ_Data result)
        {
            XYZ_Concerto_Data result_Concerto;

            result_Concerto = new XYZ_Concerto_Data()
            {

                xdata = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = result.xdata,
                    Count = result.xdata.Length
                },
                ydata = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = result.ydata,
                    Count = result.ydata.Length
                },
                zdata = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = result.zdata,
                    Count = result.zdata.Length
                },
                xtime = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = result.xtime,
                    Count = result.xtime.Length
                },
                freq_labels = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = result.freq_labels,
                    Count = result.freq_labels.Length
                },
                dimensions = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = new double[2] { result.dimensions[0], result.dimensions[1] },
                    Count = 2
                },
                overalllevel = result.overalllevel
            };

            return result_Concerto;
        }

        internal static XYZ_Concerto_Data ConvertResult_array(XYZ_Data[] result, Param_struct paramset)
        {

            XYZ_Concerto_Data result_Concerto;
            var result_packed = new XYZ_Data
            {
                xdata = new double[result.Length * result[0].xdata.Length],
                ydata = new double[result.Length * result[0].ydata.Length],
                zdata = new double[result.Length * result[0].zdata.Length],
                xtime = new double[result.Length * result[0].xtime.Length],
                freq_labels = new double[result.Length * result[0].freq_labels.Length],
                overalllevel = 0,
                dimensions = new int[2]
            };


            for (int ii = 0; ii < result.Length; ii++)
            {
                for (int jj = 0; jj < result[0].xdata.Length; jj++)
                {
                    result_packed.xdata[ii * result[0].xdata.Length + jj] = result[ii].xdata[jj];
                };
                for (int jj = 0; jj < result[0].ydata.Length; jj++)
                {
                    result_packed.ydata[ii * result[0].ydata.Length + jj] = result[ii].ydata[jj];
                };

                for (int jj = 0; jj < result[0].zdata.Length; jj++)
                {
                    result_packed.zdata[ii * result[0].zdata.Length + jj] = result[ii].zdata[jj];
                }

                for (int jj = 0; jj < result[0].xtime.Length; jj++)
                {
                    result_packed.xtime[ii * result[0].xtime.Length + jj] = result[ii].xtime[jj];
                };

                result_packed.overalllevel += result[ii].overalllevel;
            }
            result_packed.overalllevel /= result.Length;

            for (int jj = 0; jj < result[0].freq_labels.Length; jj++) result_packed.freq_labels[jj] = result[0].freq_labels[jj];

            result_packed.dimensions[0] = result.Length;
            result_packed.dimensions[1] = result[0].dimensions[1];

            if (paramset.Mean == 1)
            {
                double[] aps_xdata = new double[result.Length];
                for (int ii = 0; ii < result.Length; ii++) aps_xdata[ii] = result[ii].xdata[0];
                result_packed = AverageDiagramOverTime(result_packed, aps_xdata);
            }

            result_Concerto = new XYZ_Concerto_Data()
            {

                xdata = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = result_packed.xdata,
                    Count = result_packed.xdata.Length
                },
                ydata = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = result_packed.ydata,
                    Count = result_packed.ydata.Length
                },
                zdata = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = result_packed.zdata,
                    Count = result_packed.zdata.Length
                },
                xtime = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = result_packed.xtime,
                    Count = result_packed.xtime.Length
                },
                freq_labels = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = result_packed.freq_labels,
                    Count = result_packed.freq_labels.Length
                },
                dimensions = new IDataSet(DataSetType.Numeric)
                {
                    RealValues = new double[2] { result_packed.dimensions[0], result_packed.dimensions[1] },
                    Count = 2
                },
                overalllevel = result_packed.overalllevel
            };
            /*for (int ii = 0; ii < result_Concerto.zdata.Count; ii++)
            {
                if (Double.IsNegativeInfinity(result_Concerto.zdata[ii])) result_Concerto.zdata.
            }*/


            return result_Concerto;
        }

        public IDataSet DT2TMRPM(IDataSet dT, IDataSet cyctime, int cyclerange1, int cyclerange2)
        {
            var IOFcnClass = new DataIOFunctions();

            var dT_cyc = IOFcnClass.ReadConcertoCAData(dT, cyclerange1 - 1, cyclerange2 - 1);
            dT.Release();
            var cycrpm = DT2CycRPM(dT_cyc);

            var tmrpm = new IDataSet(DataSetType.Numeric)
            {
                Count = cycrpm.Length,
                x = cyctime,
                RealValues = cycrpm
            };
            cyctime.Release();
            return tmrpm;

        }


        private double GetSampleRate(IDataSet data, int ms)
        {
            //var dt = Diff(data.x.RealValues);
            //double samplerate = Math.Round(1 / Mean(dt), MidpointRounding.AwayFromZero);
            double samplerate = Math.Round(data.Count / (data.x[data.Count - 1] - data.x[0]));
            if (ms == 1) samplerate *= 1000;
            return samplerate;
        }

        public IDataSet Integrate(IDataSet data, int ms)
        {
            var samplerate = GetSampleRate(data, ms);
            double t_res = 1 / samplerate;
            var resultdata = new double[data.Count];
            var input = data.RealValues;
            for (int ii = 1; ii < resultdata.Length;ii++)
            {
                double area2add = (input[ii] + input[ii - 1]) / 2 * t_res;
                resultdata[ii] = area2add + resultdata[ii - 1];
            }
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = resultdata.Length,
                RealValues = resultdata,
                x = data.x
            };
            return result;
        }




        public IDataSet FilterwithPD(IDataSet audiodata, IDataSet rpm, int ms, string workdir)
        {
            double orderres = 0.05;
            int N = 512;
            double delta_t = 0.1;
            double order = 2;

            string pdpatchlocation = workdir + "\\PdPatches\\ProcessRPM.pd";
            var audiodata_x = audiodata.x.RealValues;
            var rpm_x = rpm.x.RealValues;
            var rpm_y = rpm.RealValues;

            if (ms == 1)
            {
                for (int ii = 0; ii < rpm_x.Length; ii++) rpm_x[ii] /= 1000;
                for (int ii = 0; ii < audiodata_x.Length; ii++) audiodata_x[ii] /= 1000;
                //fs *= 1000;
            }

            /*
            ///////// Test section
            string rpmfilename = "C:\\Users\\u12o24\\Documents\\ASIAT\\ASIAT\\Pd-patches\\rpm.txt";
            
            var rowscols = new int[2] { 481, 1 };
            rpm_y = new double[rowscols[0]];

            using (StreamReader sr = new StreamReader(rpmfilename))
            {

                string line;
                long ii_line = 0;
                bool success = true;

                while (!sr.EndOfStream)
                {
                    line = sr.ReadLine();
                    success = double.TryParse(line, NumberStyles.Any, CultureInfo.InvariantCulture, out rpm_y[ii_line]);
                    if (!success) break;
                    ii_line++;
                }
            }
            rpm_x = DeltaArray(0, 0.125, 60);
            //////////////////
            */

            // Calculate Samples per revolution based on order analysis parameters
            double samplesperrev = orderres * N;


            double x0 = audiodata_x[0];
            int fs = (int)(audiodata_x.Length / (audiodata_x[audiodata_x.Length - 1] - audiodata_x[0]));

            // Apply AntiAliasingFilter to prepare for resampling
            var filtered = LP10_butterworth_RPMFrequency(audiodata.RealValues, fs, new XY_Data { xdata = rpm_x, ydata = rpm_y }, samplesperrev, pdpatchlocation);
            // Interpolate RPM time axis to Rawdata resolution
            var rpm_highres = Interp1_linear(rpm_x, rpm_y, audiodata_x, true);
            // Calculate degrees per sample
            for (int ii = 0; ii < rpm_highres.Length; ii++) rpm_highres[ii] *= 6 / (double)fs; // Rev. per Minute * 360 = Degrees per minute / 60 = Degrees per second / fs = Degrees per sample
            // Add up the degrees
            var degrees = Cumsum(rpm_highres);
            
            // Calculate resolution of the new angular base in degrees
            double degreeres = 360 / samplesperrev;
            // Calculate new angular queries
            var degreequeries = DeltaArray(0, degreeres, degrees[degrees.Length - 1]);
            // Interpolate time equivalents of the new angular queries
            var timeaxis_new = Interp1_linear(degrees, audiodata_x, degreequeries, true);
            // Interpolate raw data onto new angular-equidistant time queries
            var rawdata_ip = Interp1_spline(audiodata_x, filtered, timeaxis_new, true);

            // Calculate FFT center points in time
            var t_fft = DeltaArray(x0, delta_t, audiodata_x[audiodata_x.Length - 1]);
            // Convert them into degree queries
            var ang_fft = Interp1_linear(timeaxis_new, degreequeries, t_fft, true);
            // Search the time indices on the angular axis
            for (int ii = 0; ii < ang_fft.Length; ii++) ang_fft[ii] = Math.Round(ang_fft[ii] / degreeres);
            ///////

            /*var n100 = new XY_Data
            {
                xdata = ang_fft,
                ydata = t_fft
            };*/




            //////
            var ang_fft_struct = new CATimeStamps { cycles = ang_fft };
            var paramstruct = new Param_struct
            {
                Fs = (1 / degreeres),
                Delta_f = orderres / 360,
                Windowtype = 1,
                Average = 1,
                LQ = 0,
                Delta_t = delta_t,
                DC = 1,
                Y_axis = 1,
                Freq_weight = 0,
                Y_amplitude = 0,
                Mean = 0,
                DBref = 2e-05,
                Dsrange1 = 0,
                Dsrange2 = degreequeries[degreequeries.Length - 1],
                Y_unitconversion = 1
            };
            // Calculate Order APS
            var apsresult = CalcAPSData(rawdata_ip, new XY_Data(), 0, ang_fft_struct, paramstruct);

            apsresult.overalllevel = 20 * Math.Log10(apsresult.overalllevel / paramstruct.DBref);
            for (int ii = 0; ii < apsresult.zdata.Length; ii++)
            {
                if (apsresult.zdata[ii] <= 0) apsresult.zdata[ii] = -400;
                else apsresult.zdata[ii] = 20 * Math.Log10(apsresult.zdata[ii] / paramstruct.DBref);
            }


            // x: Degrees, y: Order, z: Magnitude
            //int nx = apsresult.dimensions[0];
            int ny = apsresult.dimensions[1];

            // Extract second order
            int i_order = (int)(order / orderres);
            var indexvector = DeltaArray(i_order, ny, apsresult.zdata.Length);
            var order2 = new double[indexvector.Length];
            var degreeresult = new double[indexvector.Length];
            for (int ii = 0; ii < order2.Length; ii++)
            {
                order2[ii] = apsresult.zdata[indexvector[ii]];
                degreeresult[ii] = apsresult.xdata[indexvector[ii]];
            }

            var timeresult = Interp1_linear(degreequeries, timeaxis_new, degreeresult, true);

            for (int ii = 0; ii < timeresult.Length; ii++) timeresult[ii] = t_fft[(int)Math.Round(timeresult[ii]/delta_t)];


            //



            var time_IDS = new IDataSet(DataSetType.Numeric)
            {
                Count = timeresult.Length,
                RealValues = timeresult
            };

            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = order2.Length,
                RealValues = order2,
                x = time_IDS
            };
            return result;
        }

        /*public IDataSet ConvertTM2CAData(IDataSet rawdata, IDataSet enginespeed, string workdir, double samplesperrevolution = 360)
        {
            DataIOFunctions DataIOClass = new DataIOFunctions();
            var rawdata_y = DataIOClass.ReadConcertoYData(rawdata);
            int fs = (int)Math.Round(rawdata.Count / (rawdata.x[rawdata.Count - 1] - rawdata.x[0]));
            var rpmxy = DataIOClass.ReadConcertoData(enginespeed);
            string pdpatchlocation = workdir + "\\PdPatches\\ProcessRPM.pd";
            var result = ConvertTM2CAData_double(rawdata_y, fs, rpmxy, pdpatchlocation, samplesperrevolution);
            var result_ids = 
        }*/


        private double[] ConvertTM2CAData_double(double[] inputdata, int fs_input, XY_Data enginespeed, string pdpatchlocation, double samplesperrev = 360)
        {
            /// Samples per degree ist the sample rate in the CA domain, so the time-equivalent is samples per second.
            /// 

            var audiodata_x = DeltaArray(0, 1 / (double)fs_input, (inputdata.Length - 1) / (double)fs_input);
            

            // Apply AntiAliasingFilter to prepare for resampling
            var filtered = LP10_butterworth_RPMFrequency(inputdata, fs_input, enginespeed, samplesperrev, pdpatchlocation);

            // Interpolate RPM time axis to Rawdata resolution
            var rpm_highres = Interp1_linear(enginespeed.xdata, enginespeed.ydata, audiodata_x, true);
            // Calculate degrees per sample
            for (int ii = 0; ii < rpm_highres.Length; ii++) rpm_highres[ii] *= 6 / (double)fs_input; // Rev. per Minute * 360 = Degrees per minute / 60 = Degrees per second / fs = Degrees per sample
            // Add up the degrees
            var degrees = Cumsum(rpm_highres);
            
            // Calculate resolution of the new angular base in degrees
            double degreeres = 360 / samplesperrev;
            // Calculate new angular queries
            var degreequeries = DeltaArray(0, degreeres, degrees[degrees.Length - 1]);
            // Interpolate time equivalents of the new angular queries
            var timeaxis_new = Interp1_linear(degrees, audiodata_x, degreequeries, true);
            // Interpolate raw data onto new angular-equidistant time queries
            var rawdata_ip = Interp1_spline(audiodata_x, filtered, timeaxis_new, true);
            return rawdata_ip;
        }


        public XYZ_Data CalcOrderSpectra(double[] rawdata_y, XY_Data n100, XY_Data enginespeed, double x0, CATimeStamps timestamps, Param_struct paramset)
        {
            double orderres = paramset.Delta_o;

            double Napprox = 2.56 * paramset.O2 / orderres;
            int N = (int)Math.Pow(2, Math.Ceiling(Math.Log(Napprox) / Math.Log(2)));

            orderres = 2.56 * paramset.O2 / N;

            int freqweight = paramset.Freq_weight;
            int LQ = paramset.LQ;

            double dsrange1 = paramset.Dsrange1;
            double dsrange2 = paramset.Dsrange2;

            string workdir = paramset.Workdir;
            double delta_t = paramset.Delta_t;
            int fs = (int)Math.Round(paramset.Fs);

            // Calculate Samples per revolution based on order analysis parameters
            double samplesperrev = orderres * N;
            string pdpatchlocation = workdir + "\\PdPatches\\ProcessRPM.pd";

            //var rawdata_ip = ConvertTM2CAData(rawdata_y, fs, enginespeed, pdpatchlocation, samplesperrev);

            
            
            var audiodata_x = DeltaArray(x0, 1 / (double)fs, x0 + (rawdata_y.Length - 1) / (double)fs);
            var rpm_x = enginespeed.xdata;
            var rpm_y = enginespeed.ydata;

            //int fs = (int)(audiodata_x.Length / (audiodata_x[audiodata_x.Length - 1] - audiodata_x[0]));
            

            // Apply AntiAliasingFilter to prepare for resampling
            var filtered = LP10_butterworth_RPMFrequency(rawdata_y, fs, new XY_Data { xdata = rpm_x, ydata = rpm_y }, samplesperrev, pdpatchlocation);
            // Interpolate RPM time axis to Rawdata resolution
            var rpm_highres = Interp1_linear(rpm_x, rpm_y, audiodata_x, true);
            // Calculate degrees per sample
            for (int ii = 0; ii < rpm_highres.Length; ii++) rpm_highres[ii] *= 6 / (double)fs; // Rev. per Minute * 360 = Degrees per minute / 60 = Degrees per second / fs = Degrees per sample
            // Add up the degrees
            var degrees = Cumsum(rpm_highres);
            
            // Calculate resolution of the new angular base in degrees
            double degreeres = 360 / samplesperrev;
            // Calculate new angular queries
            var degreequeries = DeltaArray(0, degreeres, degrees[degrees.Length - 1]);
            // Interpolate time equivalents of the new angular queries
            var timeaxis_new = Interp1_linear(degrees, audiodata_x, degreequeries, true);

            degreequeries = RemoveNaN(degreequeries, timeaxis_new);
            timeaxis_new = RemoveNaN(timeaxis_new, timeaxis_new);

            // Interpolate raw data onto new angular-equidistant time queries
            var rawdata_ip = Interp1_spline(audiodata_x, filtered, timeaxis_new, true);
            
            // Calculate FFT center points in time
            var t_fft = DeltaArray(x0, delta_t, audiodata_x[audiodata_x.Length - 1]);
            // Convert them into degree queries
            var ang_fft = Interp1_linear(timeaxis_new, degreequeries, t_fft, true);
            // Search the time indices on the angular axis
            for (int ii = 0; ii < ang_fft.Length; ii++) ang_fft[ii] = Math.Round(ang_fft[ii] / degreeres);
            CATimeStamps ang_fft_struct;
            //ang_fft_struct = new CATimeStamps { cycles = ang_fft };
            ///////


            if (LQ == 0)
            {
                var t1 = new double[ang_fft.Length];
                var t2 = new double[ang_fft.Length];
                var cycles = new double[ang_fft.Length];
                for (int ii = 0; ii < ang_fft.Length; ii++)
                {
                    if ((int)ang_fft[ii] - N / 2 + 1 < 0 || (int)ang_fft[ii] + N / 2 >= degreequeries.Length)
                    {
                        t1[ii] = double.NaN;
                        t2[ii] = double.NaN;
                        cycles[ii] = double.NaN;
                    }
                    else
                    {
                        t1[ii] = degreequeries[(int)ang_fft[ii] - N / 2 + 1];
                        t2[ii] = degreequeries[(int)ang_fft[ii] + N / 2];
                        cycles[ii] = timeaxis_new[(int)ang_fft[ii]];
                    }
                }
                t1 = RemoveNaN(t1, cycles);
                t2 = RemoveNaN(t2, cycles);
                cycles = RemoveNaN(cycles, cycles);
                paramset.LQ = 2;

                ang_fft_struct = new CATimeStamps { cycles = cycles, t1 = t1, t2 = t2 };
            }
            else if (LQ == 1)
            {
                var n100_xCA = Interp1_linear(timeaxis_new, degreequeries, n100.xdata, true);
                n100.xdata = n100_xCA;
                ang_fft_struct = new CATimeStamps();
            }
            else //(LQ == 2)
            {
                ang_fft_struct.t1 = Interp1_linear(timeaxis_new, degreequeries, timestamps.t1, true);
                ang_fft_struct.t2 = Interp1_linear(timeaxis_new, degreequeries, timestamps.t2, true);
                ang_fft_struct.cycles = timestamps.cycles;
            }



            //////
            
            
            //paramset.Dsrange1 = 0;
            //paramset.Dsrange2 = degreequeries[degreequeries.Length - 1];
            var dsrange = Interp1_linear(timeaxis_new, degreequeries, new double[2] { dsrange1, dsrange2 }, true);
            paramset.Dsrange1 = dsrange[0];
            paramset.Dsrange2 = dsrange[1];

            paramset.Delta_f = orderres / 360;
            paramset.Fs = 1 / degreeres;
            paramset.Freq_weight = 0;

            // Calculate Order APS
            var apsresult = CalcAPSData(rawdata_ip, n100, 0, ang_fft_struct, paramset);




            // x: Degrees, y: Order, z: Magnitude
            int nx = apsresult.dimensions[0];
            int ny = apsresult.dimensions[1];

            var orderaxis = DeltaArray(0, orderres, (ny - 1) * orderres);

            var degreeresult = apsresult.xtime;
            var timeresult = Interp1_linear(degreequeries, timeaxis_new, degreeresult, true);
            for (int ii = 0; ii < timeresult.Length; ii++) timeresult[ii] = t_fft[(int)Math.Round((timeresult[ii] - x0) / delta_t)];
            //apsresult.xtime = timeresult;
            //if (paramset.LQ == 0)
            {
                for (int ii = 0; ii < nx; ii++)
                {
                    for (int jj = 0; jj < ny; jj++)
                    {
                        //apsresult.xdata[ii * ny + jj] = apsresult.xtime[ii];
                        apsresult.ydata[ii * ny + jj] = orderaxis[jj];
                    }
                }
            }

            apsresult.freq_labels = orderaxis;

            double[] f_FR = new double[0];
            double[] mag_FR = new double[0];

            if (paramset.DC == 0)
            {
                string filename_fr = Path.Combine(paramset.Workdir, "PdPatches", "DCFilter_FR.bin");
                
                using (var br = new BinaryReader(new FileStream(filename_fr, FileMode.Open)))
                {
                    var N2_FR = br.ReadInt32();
                    f_FR = new double[N2_FR / 2];
                    mag_FR = new double[f_FR.Length];
                    for (int ii = 0; ii < f_FR.Length; ii++) f_FR[ii] = br.ReadDouble();
                    for (int ii = 0; ii < mag_FR.Length; ii++) mag_FR[ii] = br.ReadDouble();
                    br.Close();
                }
            }


            // A-weighting
            if (freqweight != 0 || paramset.DC == 0)
            {
                var orderasfrequencies = new double[ny];
                double[] weightvector;
                var rpm_fft = Interp1_linear(rpm_x, rpm_y, timeresult, true);
                for (int ii = 0; ii < nx; ii++)
                {
                    for (int jj = 0; jj < ny; jj++) orderasfrequencies[jj] = rpm_fft[ii] / 60 * orderaxis[jj];
                    weightvector = ABC_Weight(orderasfrequencies, freqweight);

                    // DC Filter
                    if (paramset.DC == 0)
                    {
                        var DCFilter = Interp1_linear(f_FR, mag_FR, orderasfrequencies, true);
                        for (int jj = 0; jj < weightvector.Length; jj++) weightvector[jj] *= DCFilter[jj];
                    }

                    for (int jj = 0; jj < ny; jj++) apsresult.zdata[ii * ny + jj] *= weightvector[jj];
                }
            }

            //


            return apsresult;
            //



        }


        public IDataSet IIR_DCFilter(IDataSet data, int ms, string workdir)
        {
            var samplerate = GetSampleRate(data, ms);
            


            string pdpatchlocation = workdir + "\\PdPatches\\DCFilter.pd";
            var filtered_data = PdWrapper.LaunchPd(pdpatchlocation.Replace('\\', '/'), data.RealValues, (int)samplerate);
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = filtered_data.Length,
                RealValues = filtered_data,
                x = data.x
            };

            return result;
            
        }


        public IDataSet FIR_Highpass(IDataSet data, int ms, double fc = 87, int filterorder = 1024)
        {
            if (fc == 0) fc = 10;
            if (filterorder == 0) filterorder = 1024;
            var samplerate = GetSampleRate(data, ms);
            var filtered = Concerto.FIR_Highpass(data.RealValues, samplerate, fc,filterorder);
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = filtered.Length,
                RealValues = filtered,
                x = data.x
            };
            
         /*   
            string pdpatchlocation = workdir + "\\PdPatches\\ProcessRPM.pd";
            var filtered_data = PdWrapper.LaunchPd(pdpatchlocation.Replace('\\', '/'), data.RealValues, (int)samplerate);
         */

            return result;
        }

        public IDataSet GetCARangefromt1t2(IDataSet dT, int dTtype, double cyctime_ms, int cycle, double t1, double t2)
        {
            var DataIOClass = new DataIOFunctions();

            NVHFunctions.CycleData cdmtime;
            if (dTtype == 1)
            {
                cdmtime = DataIOClass.ReadConcertoCAData(dT, cycle - 1, cycle - 1, 1);
                dT.Release();
            }
            else
            {
                var dT_cyc = DataIOClass.ReadConcertoCAData(dT, cycle - 1, cycle - 1, 1);
                var cyctime_d = new double[1];
                cyctime_d[0] = cyctime_ms;
                // cyctime in ms !!!!!
                cdmtime = DT2CDMTIME(dT_cyc, cyctime_d, 1);
                dT.Release();
            }

            int i1 = 0;
            int i2 = cdmtime.Count - 1;
            double diff1 = Math.Abs(cdmtime.ydata[0] - t1);
            double diff2 = Math.Abs(cdmtime.ydata[cdmtime.Count - 1] - t2);
            for (int ii = 0; ii < cdmtime.Count; ii++)
            {
                if (Math.Abs(cdmtime.ydata[ii] - t1) <= diff1)
                {
                    diff1 = Math.Abs(cdmtime.ydata[ii] - t1);
                    i1 = ii;
                }
                else break;
            }
            for (int ii = cdmtime.Count - 1; ii >= 0; ii--)
            {
                if (Math.Abs(cdmtime.ydata[ii] - t2) <= diff2)
                {
                    diff2 = Math.Abs(cdmtime.ydata[ii] - t2);
                    i2 = ii;
                }
                else break;
            }
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = 2,
                RealValues = new double[2] { cdmtime.xdata[i1], cdmtime.xdata[i2] }
            };
            return result;


        }


        public IDataSet GetYAxis(IDataSet paramset_Concerto)
        {
            
            Param_struct paramset = new Param_struct()
            {
                Delta_f = paramset_Concerto[0],
                Delta_t = paramset_Concerto[1],
                Windowtype = (int)paramset_Concerto[2],
                Average = (int)paramset_Concerto[4],
                Average_overlap = paramset_Concerto[5],
                LQ = (int)paramset_Concerto[6],
                N_octave = paramset_Concerto[7],
                DC = (int)paramset_Concerto[8],
                Y_axis = (int)paramset_Concerto[9],
                Y_amplitude = (int)paramset_Concerto[10],
                Diagramtype = (int)paramset_Concerto[11],
                Freq_weight = (int)paramset_Concerto[12],
                F1 = paramset_Concerto[13],
                F2 = paramset_Concerto[14],
                O1 = paramset_Concerto[15],
                O2 = paramset_Concerto[16],
                Delta_o = paramset_Concerto[17],
                Tolerance = (int)paramset_Concerto[18],
                Mean = (int)paramset_Concerto[19],
                Orderstring = "",
                CDMTIME = false,
                Y_unitconversion = 1,
                DBref = 0.00002,
                Dsrange1 = paramset_Concerto[24],
                Dsrange2 = paramset_Concerto[25],
                Fs = paramset_Concerto[26]
            };

            if (paramset_Concerto[23] > 0) paramset.DBref = 0.000001;
            if (paramset_Concerto[23] == 2) paramset.Y_unitconversion = 9.81;

            var yaxis = Concerto.GetYAxis(paramset);

            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = yaxis.Length,
                RealValues = yaxis
            };

            return result;

        }

        public IDataSet GetXAxis(IDataSet rawdata_dec, IDataSet n100_dec, IDataSet t1, IDataSet t2, IDataSet cycles, IDataSet paramset_Concerto)
        {
            DataIOFunctions DataIOClass = new DataIOFunctions();


            double time_factor = 1;


            Param_struct paramset = new Param_struct()
            {
                Delta_f = paramset_Concerto[0],
                Delta_t = paramset_Concerto[1],
                Windowtype = (int)paramset_Concerto[2],
                Average = (int)paramset_Concerto[4],
                Average_overlap = paramset_Concerto[5],
                LQ = (int)paramset_Concerto[6],
                N_octave = paramset_Concerto[7],
                DC = (int)paramset_Concerto[8],
                Y_axis = (int)paramset_Concerto[9],
                Y_amplitude = (int)paramset_Concerto[10],
                Diagramtype = (int)paramset_Concerto[11],
                Freq_weight = (int)paramset_Concerto[12],
                F1 = paramset_Concerto[13],
                F2 = paramset_Concerto[14],
                O1 = paramset_Concerto[15],
                O2 = paramset_Concerto[16],
                Delta_o = paramset_Concerto[17],
                Tolerance = (int)paramset_Concerto[18],
                Mean = (int)paramset_Concerto[19],
                Orderstring = "",
                CDMTIME = false,
                Y_unitconversion = 1,
                DBref = 0.00002,
                Dsrange1 = paramset_Concerto[24],
                Dsrange2 = paramset_Concerto[25],
                Fs = paramset_Concerto[26]
            };




            if (paramset_Concerto[22] == 1) time_factor = 0.001;
            if (paramset_Concerto[23] > 0) paramset.DBref = 0.000001;
            if (paramset_Concerto[23] == 2) paramset.Y_unitconversion = 9.81;

            XY_Data n100;


            double x0;

            int L;

            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = 0
            };

            try
            {
                if (rawdata_dec.CycleCount > 1)
                {
                    L = (int)((paramset.Dsrange2 - paramset.Dsrange1) * paramset.Fs);
                    x0 = paramset.Dsrange1;

                }
                else
                {
                    x0 = rawdata_dec.x[0] * time_factor;


                    if (paramset.Dsrange1 < x0) paramset.Dsrange1 = x0;
                    if (paramset.Dsrange2 > time_factor * rawdata_dec.x[rawdata_dec.Count - 1]) paramset.Dsrange2 = time_factor * rawdata_dec.x[rawdata_dec.Count - 1];
                    L = (int)((paramset.Dsrange2 - paramset.Dsrange1) * paramset.Fs);
                }



                n100 = DataIOClass.ReadConcertoData(n100_dec);
                for (int ii = 0; ii < n100.xdata.Length; ii++) n100.xdata[ii] *= time_factor;


                var timestamps = new NVHFunctions.CATimeStamps
                {
                    t1 = t1.RealValues,
                    t2 = t2.RealValues,
                    cycles = cycles.RealValues
                };


                var rawdataY = new double[1];
                rawdataY[0] = L;



                var xaxis = Concerto.GetXAxis(L, n100, x0, timestamps, paramset);

                result = new IDataSet(DataSetType.Numeric)
                {
                    Count = xaxis.Length,
                    RealValues = xaxis
                };
            }
            catch (Exception e)
            {
                var rawdataYtest = new double[1] { rawdata_dec.Count };
                var enginespeed = new XY_Data();
                var x0test = rawdata_dec.x[0] * time_factor;
                var dttest = (rawdata_dec.x[1] - rawdata_dec.x[0]) * time_factor;
                var docdir = Environment.SpecialFolder.MyDocuments;
                var xlspath = Path.Combine(Environment.GetFolderPath(docdir), "nvhdebugdata.xls");

                if (!alreadywritten) IO.SaveAPSInput(xlspath, rawdataYtest, enginespeed, dttest, x0test, paramset, e.Message);
                alreadywritten = true;
            }

            

            return result;
        }

        public int GetDimensions(IDataSet rawdata_dec, IDataSet n100_dec, IDataSet enginespeed_dec, IDataSet t1, IDataSet t2, IDataSet cycles, IDataSet paramset_Concerto)
        {
            DataIOFunctions DataIOClass = new DataIOFunctions();


            double time_factor = 1;


            Param_struct paramset = new Param_struct()
            {
                Delta_f = paramset_Concerto[0],
                Delta_t = paramset_Concerto[1],
                Windowtype = (int)paramset_Concerto[2],
                Average = (int)paramset_Concerto[4],
                Average_overlap = paramset_Concerto[5],
                LQ = (int)paramset_Concerto[6],
                N_octave = paramset_Concerto[7],
                DC = (int)paramset_Concerto[8],
                Y_axis = (int)paramset_Concerto[9],
                Y_amplitude = (int)paramset_Concerto[10],
                Diagramtype = (int)paramset_Concerto[11],
                Freq_weight = (int)paramset_Concerto[12],
                F1 = paramset_Concerto[13],
                F2 = paramset_Concerto[14],
                O1 = paramset_Concerto[15],
                O2 = paramset_Concerto[16],
                Delta_o = paramset_Concerto[17],
                Tolerance = (int)paramset_Concerto[18],
                Mean = (int)paramset_Concerto[19],
                Orderstring = "",
                CDMTIME = false,
                Y_unitconversion = 1,
                DBref = 0.00002,
                Dsrange1 = paramset_Concerto[24],
                Dsrange2 = paramset_Concerto[25],
                Fs = paramset_Concerto[26]
            };

            


            if (paramset_Concerto[22] == 1) time_factor = 0.001;
            if (paramset_Concerto[23] > 0) paramset.DBref = 0.000001;
            if (paramset_Concerto[23] == 2) paramset.Y_unitconversion = 9.81;

            XY_Data n100;
            XY_Data enginespeed;

            if (paramset.Diagramtype != 3 && paramset.Diagramtype != 4) enginespeed_dec = new IDataSet(DataSetType.Numeric)
            {
                Count = 0,
                RealValues = new double[0]
            };

            double x0;
            
            int L;

            if (rawdata_dec.CycleCount > 1)
            {
                L = (int)((paramset.Dsrange2 - paramset.Dsrange1) * paramset.Fs);
                x0 = paramset.Dsrange1;
                
            }
            else
            {
                x0 = rawdata_dec.x[0] * time_factor;
                

                if (paramset.Dsrange1 < x0) paramset.Dsrange1 = x0;
                if (paramset.Dsrange2 > time_factor * rawdata_dec.x[rawdata_dec.Count - 1]) paramset.Dsrange2 = time_factor * rawdata_dec.x[rawdata_dec.Count - 1];
                L = (int)((paramset.Dsrange2 - paramset.Dsrange1) * paramset.Fs);
            }



            n100 = DataIOClass.ReadConcertoData(n100_dec);
            for (int ii = 0; ii < n100.xdata.Length; ii++) n100.xdata[ii] *= time_factor;
            enginespeed = DataIOClass.ReadConcertoData(enginespeed_dec);
            for (int ii = 0; ii < enginespeed.xdata.Length; ii++) enginespeed.xdata[ii] *= time_factor;

            var timestamps = new NVHFunctions.CATimeStamps
            {
                t1 = t1.RealValues,
                t2 = t2.RealValues,
                cycles = cycles.RealValues
            };


            var rawdataY = new double[1];
            rawdataY[0] = L;
            


            int dimensions = CheckDimension(L, n100, enginespeed,x0, timestamps, paramset);

            return dimensions;
        }

        public string GetErrorMessage()
        {
            if (errormessage == null) return "ok";
            else return errormessage;
        }

        public IDataSet GetOctaveFreqs(double n, double fs, double N)
        {
            var freqs = CalcCornerFreqs(n, fs, N);
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = freqs.f0.Length,
                RealValues = freqs.f0
            };
            return result;
        }

        public IDataSet GetFrequencyLabels()
        {
            if (apsresult_freq_labels == null)
            {
                var result = new IDataSet(DataSetType.Numeric)
                {
                    Count = 1,
                    RealValues = new double[1] { 0 }
                };
                return result;
            }
            else
            {
                var result = new IDataSet(DataSetType.Numeric)
                {
                    Count = apsresult_freq_labels.Length,
                    RealValues = apsresult_freq_labels
                };
                return result;
            }
            
        }

        public IDataSet GetOrderAxis(string orderstring, int mirrornegativeorders = 0)
        {
            bool negativeorders = false;
            if (mirrornegativeorders == 1) negativeorders = true;
            var oaxis = new IDataSet(DataSetType.Numeric);


            double[] oaxis_d = ParseOrdersFromString(orderstring);
            if (negativeorders) oaxis.Count = oaxis_d.Length * 2;
            else oaxis.Count = oaxis_d.Length;
            if (negativeorders)
            {
                var oaxis_d_mirror = new double[oaxis.Count];
                Array.Copy(oaxis_d, 0, oaxis_d_mirror, oaxis_d.Length, oaxis_d.Length);
                for (int ii = 1; ii <= oaxis_d.Length; ii++) oaxis_d_mirror[oaxis_d.Length - ii] = -oaxis_d[ii - 1];
                oaxis_d = oaxis_d_mirror;
            }
            oaxis.RealValues = oaxis_d;


            return oaxis;
        }

        public IDataSet GetOrderAxis(double o1, double o2, double delta_o, int mirrornegativeorders = 0)
        {
            bool negativeorders = false;
            o1 = 0;
            if (mirrornegativeorders == 1) negativeorders = true;
            int rows_out = (int)((o2 - o1) / delta_o) + 1;
            IDataSet oaxis = new IDataSet(DataSetType.Numeric);
            if (negativeorders) oaxis.Count = rows_out * 2;
            else oaxis.Count = rows_out;

            var oaxis_d = new double[rows_out];
            for (int ii = 0; ii < rows_out; ii++) oaxis_d[ii] = ii * delta_o + o1;

            if (negativeorders)
            {
                var oaxis_d_mirror = new double[oaxis.Count];
                Array.Copy(oaxis_d, 0, oaxis_d_mirror, oaxis_d.Length, oaxis_d.Length);
                for (int ii = 1; ii <= oaxis_d.Length; ii++) oaxis_d_mirror[oaxis_d.Length - ii] = -oaxis_d[ii - 1];
                oaxis_d = oaxis_d_mirror;
            }

            oaxis.RealValues = oaxis_d;
            return oaxis;
        }

        public double GetOverallLevel()
        {
            return apsresult_overalllevel;
        }

        public double GetPlaybackPosition()
        {
            double result = 0;
            if (waveOut != null)
            {
                if (waveOut.PlaybackState != NAudio.Wave.PlaybackState.Stopped)
                {
                    long playposition_bytes = waveOut.GetPosition();
                    result = 100 * (double)playposition_bytes / playbacklength_bytes;
                }
            }
            return result;

        }

        public string GetVersion()
        {
            Version version = Assembly.GetExecutingAssembly().GetName().Version;
            DateTime buildDate = new DateTime(2000, 1, 1).AddDays(version.Build).AddSeconds(version.Revision * 2);
            string displayableVersion = $"{version} ({buildDate})";
            return displayableVersion;
        }


        public IDataSet GetXData()
        {
            IDataSet result;
            if (apsresult_xdata == null)
            {
                result = new IDataSet(DataSetType.Numeric)
                {
                    Count = 1,
                    RealValues = new double[1] { 0 }
                };
            }
            else
            {
                result = new IDataSet(DataSetType.Numeric)
                {
                    Count = apsresult_xdata.Length,
                    RealValues = apsresult_xdata
                };
                
            }

            if (allowcommonxfeature)
            {
                var xds = new IDataSet(DataSetType.Numeric)
                {
                    Count = apsresult_commonxdata.Length,
                    RealValues = apsresult_commonxdata
                };
                result.x = xds;
                xds.Release();
            }
            

            return result;
            /*if (xyzdata[0] == null)
            {
                var result = new IDataSet(DataSetType.Numeric)
                {
                    Count = 1,
                    RealValues = new double[1] { 0 }
                };
                return result;
            }
            else
            {
                return xyzdata[0];
            }*/

        }

        public IDataSet GetLastZData()
        {
            IDataSet result;
            if (apsresult_zdata == null)
            {
                result = new IDataSet(DataSetType.Numeric)
                {
                    Count = 1,
                    RealValues = new double[1] { 0 }
                };
            }
            else
            {
                result = new IDataSet(DataSetType.Numeric)
                {
                    Count = apsresult_zdata.Length,
                    RealValues = apsresult_zdata
                };

            }
            return result;
            /*
            if (xyzdata[2] == null)
            {
                var result = new IDataSet(DataSetType.Numeric)
                {
                    Count = 1,
                    RealValues = new double[1] { 0 }
                };
                return result;
            }
            else
            {
                return xyzdata[2];
            }*/
        }

        public IDataSet GetXTime()
        {
            IDataSet result;
            if (apsresult_xtime == null)
            {
                result = new IDataSet(DataSetType.Numeric)
                {
                    Count = 1,
                    RealValues = new double[1] { 0 }
                };
            }
            else
            {
                result = new IDataSet(DataSetType.Numeric)
                {
                    Count = apsresult_xtime.Length,
                    RealValues = apsresult_xtime
                };
                
            }
            return result;
        }

        public IDataSet GetXYDimensions()
        {
            IDataSet result;
            if (apsresult_dimensions == null)
            {
                result = new IDataSet(DataSetType.Numeric)
                {
                    Count = 1,
                    RealValues = new double[1] { 0 }
                };
            }
            else
            {
                result = new IDataSet(DataSetType.Numeric)
                {
                    Count = apsresult_dimensions.Length,
                    RealValues = apsresult_dimensions
                };

            }
            return result;
        }

        public double GetAPSTimestamp()
        {
            return aps_timestamp;
        }

        public IDataSet GetYData()
        {
            IDataSet result;
            if (apsresult_ydata == null)
            {
                result = new IDataSet(DataSetType.Numeric)
                {
                    Count = 1,
                    RealValues = new double[1] { 0 }
                };
            }
            else
            {
                result = new IDataSet(DataSetType.Numeric)
                {
                    Count = apsresult_ydata.Length,
                    RealValues = apsresult_ydata
                };

            }

            if (allowcommonxfeature)
            {
                var xds = new IDataSet(DataSetType.Numeric)
                {
                    Count = apsresult_commonxdata.Length,
                    RealValues = apsresult_commonxdata
                };
                result.x = xds;
                xds.Release();
            }

            

            return result;
            /*
            if (xyzdata[1] == null)
            {
                var result = new IDataSet(DataSetType.Numeric)
                {
                    Count = 1,
                    RealValues = new double[1] { 0 }
                };
                return result;
            }
            else
            {
                return xyzdata[1];
            }*/
        }

        public double GetNfromDeltaF(double fs, double deltaf)
        {
            if (deltaf <= 0) deltaf = fs / 16384;
            else if (fs / deltaf > 262144) deltaf = fs / 16384;

            int N = Nextpow2(Math.Round(fs) / deltaf);
            return N;
        }


        bool alreadywritten = false;
        public void GetZData_Debug(IDataSet rawdata_dec, IDataSet enginespeed_dec, IDataSet paramset_Concerto, string ostring, string errormessage)
        {
            if (alreadywritten) return;
            double time_factor = 1;
            
            Param_struct paramset = new Param_struct()
            {
                Delta_f = paramset_Concerto[0],        // neu 14.11. (vorher N)
                Delta_t = paramset_Concerto[1],
                Windowtype = (int)paramset_Concerto[2],
                Average = (int)paramset_Concerto[4],
                Average_overlap = paramset_Concerto[5],
                LQ = (int)paramset_Concerto[6],
                N_octave = paramset_Concerto[7],
                DC = (int)paramset_Concerto[8],
                Y_axis = (int)paramset_Concerto[9],
                Y_amplitude = (int)paramset_Concerto[10],
                Diagramtype = (int)paramset_Concerto[11],
                Freq_weight = (int)paramset_Concerto[12],
                F1 = paramset_Concerto[13],
                F2 = paramset_Concerto[14],
                O1 = paramset_Concerto[15],
                O2 = paramset_Concerto[16],
                Delta_o = paramset_Concerto[17],
                Tolerance = (int)paramset_Concerto[18],
                Mean = (int)paramset_Concerto[19],
                Orderstring = ostring,
                Y_unitconversion = 1,
                DBref = 0.00002,
                Dsrange1 = paramset_Concerto[24],
                Dsrange2 = paramset_Concerto[25],
                Fs = paramset_Concerto[26],
                FromCA = (int)paramset_Concerto[20],
                ToCA = (int)paramset_Concerto[21],
                Inverter_frequency = paramset_Concerto[31],
                Polepairs = (int)paramset_Concerto[32]
            };

            if (paramset_Concerto[30] == 1) paramset.Calcinverternoise = true;
            
            if (paramset_Concerto[3] == 1) paramset.CDMTIME = true;

            
            if (paramset_Concerto[22] == 1)
            {
                time_factor = 0.001;
            
            }
            if (paramset_Concerto[23] > 0) paramset.DBref = 0.000001;
            if (paramset_Concerto[23] == 2) paramset.Y_unitconversion = 9.81;
            


            var rawdataYtest = rawdata_dec.RealValues;
            DataIOFunctions DataIOClass = new DataIOFunctions();
            var enginespeed = DataIOClass.ReadConcertoData(enginespeed_dec);
            if (time_factor != 1) for (int ii = 0; ii < enginespeed.xdata.Length; ii++) enginespeed.xdata[ii] *= time_factor;
            var x0test = rawdata_dec.x[0] * time_factor;
            var dttest = (rawdata_dec.x[1] - rawdata_dec.x[0]) * time_factor;
            var docdir = Environment.SpecialFolder.MyDocuments;
            var xlspath = Path.Combine(Environment.GetFolderPath(docdir), "nvhdebugdata.xls");

            IO.SaveAPSInput(xlspath, rawdataYtest, enginespeed, dttest, x0test, paramset, errormessage);
            alreadywritten = true;
        }


        

        
        public IDataSet GetZData(IDataSet rawdata_dec, IDataSet dT, IDataSet cyctime, IDataSet n100_dec, IDataSet enginespeed_dec, IDataSet paramset_Concerto, string ostring, string workdir)
        {
            //Calculate using rawdata from Concerto and converting them into double values before calculation start

            try
            {
                aps_timestamp = rawdata_dec.x[rawdata_dec.Count - 1];
                DataIOFunctions DataIOClass = new DataIOFunctions();
                

                bool formulacalc = false;
                double time_factor = 1;

                Param_struct paramset = new Param_struct()
                {
                    Delta_f = paramset_Concerto[0],        // neu 14.11. (vorher N)
                    Delta_t = paramset_Concerto[1],
                    Windowtype = (int)paramset_Concerto[2],
                    Average = (int)paramset_Concerto[4],
                    Average_overlap = paramset_Concerto[5],
                    LQ = (int)paramset_Concerto[6],
                    N_octave = paramset_Concerto[7],
                    DC = (int)paramset_Concerto[8],
                    Y_axis = (int)paramset_Concerto[9],
                    Y_amplitude = (int)paramset_Concerto[10],
                    Diagramtype = (int)paramset_Concerto[11],
                    Freq_weight = (int)paramset_Concerto[12],
                    F1 = paramset_Concerto[13],
                    F2 = paramset_Concerto[14],
                    O1 = paramset_Concerto[15],
                    O2 = paramset_Concerto[16],
                    Delta_o = paramset_Concerto[17],
                    Tolerance = (int)paramset_Concerto[18],
                    Mean = (int)paramset_Concerto[19],
                    Orderstring = ostring,
                    Y_unitconversion = 1,
                    DBref = 0.00002,
                    Dsrange1 = paramset_Concerto[24],
                    Dsrange2 = paramset_Concerto[25],
                    Fs = paramset_Concerto[26],
                    FromCA = (int)paramset_Concerto[20],
                    ToCA = (int)paramset_Concerto[21],
                    Inverter_frequency = paramset_Concerto[31],
                    Polepairs = (int)paramset_Concerto[32],
                    Workdir = workdir
                };

                if (paramset.Diagramtype != 3 && paramset.Diagramtype != 4) enginespeed_dec = new IDataSet(DataSetType.Numeric)
                {
                    Count = 0,
                    RealValues = new double[0]
                };

                if (paramset_Concerto[30] == 1) paramset.Calcinverternoise = true;
                //if (paramset_Concerto[19] == 1) paramset.Mean = true;
                if (paramset_Concerto[3] == 1) paramset.CDMTIME = true;

                int ms = 0;
                if (paramset_Concerto[22] == 1)
                {
                    time_factor = 0.001;
                    ms = 1;
                }
                if (paramset_Concerto[23] > 0) paramset.DBref = 0.000001;
                if (paramset_Concerto[23] == 2) paramset.Y_unitconversion = 9.81;
                if (paramset_Concerto[29] == 1) formulacalc = true;



                XY_Data enginespeed;
                XY_Data n100;
                var timestamps = new CATimeStamps();
                /*{
                    t1 = new double[1],
                    t2 = new double[1],
                    cycles = new double[1]

                };*/
                double[] rawdataY = Array.Empty<double>();

                double x0 = 0;
                double dt = 0;


                CycleData cdmtime;
                int[] cyclerange = new int[2];

                // Calculate cycle range
                if (paramset.LQ == 2 || rawdata_dec.CycleCount > 1)
                {
                    cyclerange[0] = (int)paramset_Concerto[27] - 1;
                    cyclerange[1] = (int)paramset_Concerto[28] - 1;
                }

                

                // Calculate Time Stamps
                if (paramset.LQ == 2)
                {
                    if (paramset.CDMTIME)
                    {
                        cdmtime = DataIOClass.ReadConcertoCAData(dT, cyclerange[0], cyclerange[1], (int)paramset.Delta_t);
                    }
                    else
                    {
                        var dT_cyc = DataIOClass.ReadConcertoCAData(dT, cyclerange[0], cyclerange[1], (int)paramset.Delta_t);
                        var cyctime_d = new double[dT_cyc.CycleCount];
                        int kk = 0;
                        for (int ii = cyclerange[0]; ii <= cyclerange[1]; ii += (int)paramset.Delta_t)
                        {
                            cyctime_d[kk] = cyctime.RealValues[ii];
                            kk++;
                        }
                        cdmtime = DT2CDMTIME(dT_cyc, cyctime_d, ms);
                    }



                    timestamps = GetCATimeStamps_CDMTIME(cdmtime, paramset.FromCA, paramset.ToCA, ms);
                    for (int ii = 0; ii < timestamps.cycles.Length; ii++) timestamps.cycles[ii] = ii * (int)paramset.Delta_t + cyclerange[0] + 1;

                }

                XYZ_Concerto_Data apsresult;
                
                if (rawdata_dec.CycleCount > 1)
                {
                    // Resample CA Data
                    // Do special "fast forward" calculation for cyclic data and delta > 2
                    if (paramset.LQ == 2 && paramset.Delta_t > 2)
                    {

                        // Collect data of cycles to calcuate into an array
                        enginespeed = DataIOClass.ReadConcertoData(enginespeed_dec);
                       
                        for (int ii = 0; ii < enginespeed.xdata.Length; ii++) enginespeed.xdata[ii] *= time_factor;
                        
                        //// Prepare Arrays for cyclewise Calculation
                        CycleData[] rawdata_array;
                        if (!formulacalc) rawdata_array = ReadCADataIntoCycleDataArray(rawdata_dec, timestamps);
                        else rawdata_array = cadatainputarray;
                        var dT_array = ReadCADataIntoCycleDataArray(dT, timestamps);


                        //// Do calculation for each cycle to calculate
                        var cyctime_s = cyctime.RealValues;
                        for (int ii = 0; ii < cyctime_s.Length; ii++) cyctime_s[ii] *= time_factor;
                        
                        apsresult = CalculateCycleWise(rawdata_array, enginespeed, dT_array, cyctime_s, timestamps, paramset);

                        rawdata_dec.Release();
                        dT.Release();
                        cyctime.Release();
                        paramset_Concerto.Release();
                        n100_dec.Release();
                        enginespeed_dec.Release();

                        apsresult.dimensions.Release();
                        apsresult.xdata.Release();
                        apsresult.ydata.Release();
                        apsresult.freq_labels.Release();
                        apsresult.xtime.Release();

                        return apsresult.zdata;
                    }
                    else
                    // Resample whole CA dataset if delta < 3 or if leading quantity is not cycle.
                    {
                        int add = 0;
                        if (cyclerange[1] < rawdata_dec.CycleCount - 1) add = 1;
                        CycleData raw_d;
                        if (!formulacalc) raw_d = DataIOClass.ReadConcertoCAData(rawdata_dec, cyclerange[0], cyclerange[1] + add);
                        else raw_d = cadatainput;
                        
                        var dT_d = DataIOClass.ReadConcertoCAData(dT, cyclerange[0], cyclerange[1] + add);
                        

                        if (paramset.CDMTIME) rawdataY = ResampleCAData_CDMTIME(raw_d, dT_d, paramset.Fs);
                        else rawdataY = ResampleCAData_dT(raw_d, dT_d, paramset.Fs);

                        x0 = (cyctime.RealValues[cyclerange[0]]) * time_factor;
                        dt = 1 / paramset.Fs;
                    }
                }
                else
                {
                    rawdataY = rawdata_dec.RealValues;
                    x0 = rawdata_dec.x[0] * time_factor;
                    dt = (rawdata_dec.x[rawdata_dec.Count - 1] - rawdata_dec.x[0]) / rawdata_dec.Count * time_factor;
                }

               
                

                int dsrange1_index = (int)((paramset.Dsrange1 - x0) / dt);
                int dsrange2_index = (int)((paramset.Dsrange2 - x0) / dt);

                // Resample TM data if necessary
                if (Math.Round(paramset.Fs) != Math.Round(1 / dt))
                {
                    if (paramset.Dsrange2 - paramset.Dsrange1 < 0.7 * (rawdata_dec.x[rawdata_dec.Count - 1] * time_factor - x0))  ///???
                    {
                        int add = 0;
                        if (dsrange2_index + 23 < rawdata_dec.Count) add = 23;
                        var rawdataY_red = new double[dsrange2_index - dsrange1_index + 1 + add];
                        for (int ii = 0; ii < rawdataY_red.Length; ii++) rawdataY_red[ii] = rawdataY[dsrange1_index + ii];
                        x0 = paramset.Dsrange1;
                        rawdataY = ResampleTMData(rawdataY_red, x0, dt, paramset.Fs);
                    }
                    else rawdataY = ResampleTMData(rawdataY, x0, dt, paramset.Fs);
                    dt = 1 / paramset.Fs;
                }
                

                if (paramset.Dsrange1 < x0) paramset.Dsrange1 = x0;
                if (paramset.Dsrange2 > x0 + dt * (rawdataY.Length - 1)) paramset.Dsrange2 = x0 + dt * (rawdataY.Length - 1);

                enginespeed = DataIOClass.ReadConcertoData(enginespeed_dec);
                
                n100 = DataIOClass.ReadConcertoData(n100_dec);
                
                for (int ii = 0; ii < enginespeed.xdata.Length; ii++) enginespeed.xdata[ii] *= time_factor;
                for (int ii = 0; ii < n100.xdata.Length; ii++) n100.xdata[ii] *= time_factor;



                XYZ_Data apsresult_double;
                switch (paramset.Diagramtype)
                {
                    case 1:

                        
                        
                        
                        apsresult_double = CalcAPSData(rawdataY, n100, x0, timestamps, paramset);
                        break;
                    case 2:
                        

                        apsresult_double = CalcNtelOctaveSpectra(rawdataY, n100, dt, x0, timestamps, paramset);
                        break;
                    case 3:
                        if (!paramset.Calcinverternoise) apsresult_double = CalcOrderSpectra(rawdataY, n100, enginespeed, x0, timestamps, paramset);
                        else apsresult_double = CalcOrderSpectra_FFT(rawdataY, n100, enginespeed, x0, timestamps, paramset);
                        break;
                    case 4:
                        if (!paramset.Calcinverternoise) apsresult_double = CalcOrderSpectra(rawdataY, n100, enginespeed, x0, timestamps, paramset);
                        else apsresult_double = CalcOrderSpectra_FFT(rawdataY, n100, enginespeed, x0, timestamps, paramset);
                        break;
                    case 5:
                        apsresult_double = CalcOverallLevel(rawdataY, n100, x0, timestamps, paramset);
                        break;
                    case 6:
                        apsresult_double = CalcBPLevel(rawdataY, n100, x0, timestamps, paramset);
                        break;
                    default:
                        apsresult_double = CalcAPSData(rawdataY, n100, x0, timestamps, paramset);
                        break;
                }
                if (paramset.Y_axis == 1)
                {
                    apsresult_double.overalllevel = 20 * Math.Log10(apsresult_double.overalllevel / paramset.DBref);
                    for (int ii = 0; ii < apsresult_double.zdata.Length; ii++)
                    {
                        if (apsresult_double.zdata[ii] <= 0) apsresult_double.zdata[ii] = -400;
                        else apsresult_double.zdata[ii] = 20 * Math.Log10(apsresult_double.zdata[ii] / paramset.DBref);
                    }
                    
                }

                


                if (apsresult_double.errorstring != null) errormessage = apsresult_double.errorstring;

                apsresult = ConvertResult(apsresult_double);

                rawdata_dec.Release();
                dT.Release();
                cyctime.Release();
                paramset_Concerto.Release();
                n100_dec.Release();
                enginespeed_dec.Release();
                apsresult_zdata = apsresult.zdata.RealValues;
                apsresult_xdata = apsresult.xdata.RealValues;
                apsresult_ydata = apsresult.ydata.RealValues;
                apsresult_xtime = apsresult.xtime.RealValues;
                apsresult_freq_labels = apsresult.freq_labels.RealValues;
                apsresult_dimensions = apsresult.dimensions.RealValues;
                apsresult_overalllevel = apsresult.overalllevel;
                /*xyzdata[0] = apsresult.xdata;
                xyzdata[1] = apsresult.ydata;
                xyzdata[2] = apsresult.zdata;*/

                if (allowcommonxfeature)
                {
                    if (apsresult.dimensions[0] == 1 || apsresult.dimensions[1] == 1)
                    {
                        if (paramset.Mean == 1)
                        {
                            apsresult_commonxdata = apsresult.ydata.RealValues;
                        }
                        else
                        {
                            apsresult_commonxdata = apsresult.xdata.RealValues;
                        }

                    }
                    else
                    {
                        apsresult_commonxdata = new double[apsresult_zdata.Length];
                        for (int ii = 0; ii < apsresult_commonxdata.Length; ii++) apsresult_commonxdata[ii] = ii + 1;
                    }

                    var xds = new IDataSet(DataSetType.Numeric)
                    {
                        Count = apsresult_commonxdata.Length,
                        RealValues = apsresult_commonxdata
                    };
                    apsresult.zdata.x = xds;
                    xds.Release();
                }
                
                apsresult.dimensions.Release();
                apsresult.xdata.Release();
                apsresult.ydata.Release();
                apsresult.freq_labels.Release();
                apsresult.xtime.Release();
                

                return apsresult.zdata;

            }
            catch (Exception e)
            {
                errormessage = e.ToString();
                /*GetZData_Debug(rawdata_dec, enginespeed_dec,paramset_Concerto, ostring, errormessage);
                IDataSet errorresult = new IDataSet(DataSetType.Numeric)
                {
                    Count = 1
                };
                errorresult.x[0] = 0;
                errorresult.y[0] = 0;*/
                /*var dummyresult = new XYZ_Data
                {
                    dimensions = new int[2] { 0, 0 },
                    errorstring = errormessage,
                    freq_labels = new double[0],
                    overalllevel = 0,
                    rows = 0,
                    xtime = new double[0],
                    xdata = new double[0],
                    ydata = new double[0],
                    zdata = new double[0]
                };
                var apsresult = ConvertResult(dummyresult);*/

                var result = new IDataSet(DataSetType.Numeric)
                {
                    Count = 1,
                    RealValues = new double[1] { 0 }
                };

                rawdata_dec.Release();
                dT.Release();
                cyctime.Release();
                paramset_Concerto.Release();
                n100_dec.Release();
                enginespeed_dec.Release();
                return result;

            }

        }



        
        

       

        
        public IDataSet Get_cycles()
        {
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = time_stamps.cycles.Length,
                RealValues = time_stamps.cycles
            };
            return result;
        }
        public IDataSet Get_t1()
        {
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = time_stamps.t1.Length,
                RealValues = time_stamps.t1
            };
            return result;
        }
        public IDataSet Get_t2()
        {
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = time_stamps.t2.Length,
                RealValues = time_stamps.t2
            };
            return result;
        }

        public IDataSet Interp1(IDataSet xy, IDataSet xq)
        {
            var xd = xy.x.RealValues;
            var yd = xy.RealValues;
            var xqd = xq.RealValues;
            xy.Release();
            xq.Release();

            var yqd = Interp1_linear(xd, yd, xqd, false);
            yqd = RemoveNaN(yqd, yqd);
            var yq = new IDataSet(DataSetType.Numeric)
            {
                RealValues = yqd,
                Count = yqd.Length
            };
            return yq;
        }
        private void OnPlaybackStopped(object sender, EventArgs e)
        {
            waveOut.Dispose();
        }


        public IDataSet MovingAverage(IDataSet data, int filterlength)
        {
            var ydata = data.RealValues;
            var ydata_ma = GeneralMath.Moving(ydata, filterlength);
            if (ydata_ma == null) return data;
            Array.Resize(ref ydata_ma, ydata_ma.Length - filterlength / 2);
            var x_short = new double[ydata_ma.Length];
            Array.Copy(data.x.RealValues, 0, x_short, 0, x_short.Length);

            var ydata_ep = GeneralMath.Interp1_linear(x_short, ydata_ma, data.x.RealValues, false, ydata_ma[ydata_ma.Length - 1]);


            
            var result = new IDataSet(DataSetType.Numeric)
            {
                x = data.x,
                RealValues = ydata_ep
            };
            data.Release();
            return result;
        }

        public void PausePlayback()
        {
            if (waveOut != null) if (waveOut.PlaybackState == NAudio.Wave.PlaybackState.Playing) waveOut.Pause();
        }

        public int PlayData(IDataSet data)
        {
            if (NAudio.Wave.WaveOut.DeviceCount < 1) return -1;


            waveOut = new NAudio.Wave.WaveOut();
            byte[] bytearr = new byte[2];
            byte[] sound_bytes = new byte[2 * data.Count];



            int fs = (Int32)(1 / (data.x[2] - data.x[1]));

            double maxpascal = data.Max;
            double dB = 20 * Math.Log10(2 * Math.Sqrt(2) / 0.00002);
            while (20 * Math.Log10(maxpascal / 0.00002) > dB) dB += 10;
            double pascal = 0.00002 * Math.Pow(10, dB / 20);
            double[] dataarray = data.RealValues;
            data.Release();


            for (int ii = 0; ii < dataarray.Length; ii++)
            {
                bytearr = BitConverter.GetBytes((Int16)(dataarray[ii] * 32768 / pascal));
                sound_bytes[ii * 2] = bytearr[0];
                sound_bytes[ii * 2 + 1] = bytearr[1];
            }
            playbacklength_bytes = sound_bytes.LongLength;

            provider = new NAudio.Wave.RawSourceWaveStream(new System.IO.MemoryStream(sound_bytes), new NAudio.Wave.WaveFormat(fs, 16, 1));

            waveOut.Init(provider);
            waveOut.Play();
            waveOut.PlaybackStopped += OnPlaybackStopped;

            return 0;

        }

        public string PlayOnline(IDataSet data, int latency, int buffersize)
        {
            // data has to be provided with samplerate 50 kHz
            // x axis has to be in ms
            // signal has to be already calibrated so that the range is between -1 and 1
            try
            {

                //int latency = 1000;
                //int buffersize = 500;
                int samplerate = 50000;
                int nbytes = 2;


                double L_ms = data.x[data.Count - 1] - data.x[0];

                if (L_ms < latency)
                {
                    blockcount_onlineplayback = 0;
                    timestamp_onlineplayback = 0;
                    if (bufferedWaveProvider != null)
                    {
                        if (bufferedWaveProvider.BufferedBytes > 0) bufferedWaveProvider.ClearBuffer();
                    }
                }


                if (L_ms < latency) return "Not enough data.";
                else if (data.x[data.Count - 1] - timestamp_onlineplayback < buffersize)
                {
                    return "skipped";    // not enough data available or still enough data in buffer
                }


                // init BufferedWaveProvider and WaveOut
                bool wpInitialized = false;
                bool woInitialized = false;
                if (bufferedWaveProvider == null)
                {
                    bufferedWaveProvider = new NAudio.Wave.BufferedWaveProvider(new NAudio.Wave.WaveFormat(samplerate, nbytes * 8, 1))
                    {
                        BufferLength = samplerate * nbytes * 10,
                        DiscardOnBufferOverflow = false
                    };
                    wpInitialized = true;
                }
                if (waveOut == null)
                {
                    waveOut = new NAudio.Wave.WaveOut();
                    woInitialized = true;
                }
                if (waveOut.PlaybackState != NAudio.Wave.PlaybackState.Playing)
                {
                    waveOut.Init(bufferedWaveProvider);
                    waveOut.Play();
                }

                int L_samples = (int)((double)buffersize / 1000 * samplerate);
                double index1_time = blockcount_onlineplayback * (double)L_samples / samplerate * 1000;
                int index1 = Array.FindIndex(data.x.RealValues, f => f >= index1_time);
                

                double[] newaudio = new double[L_samples];
                Array.Copy(data.RealValues, index1, newaudio, 0, L_samples);
                //var newaudio = data.RealValues;
                timestamp_onlineplayback = data.x[data.Count - 1];
                data.Release();


                int L_bytes = L_samples * nbytes;
                byte[] buffer = new byte[L_bytes];

                string message = "";
                if (wpInitialized) message += "WaveProvider initialized.";
                if (woInitialized) message += "WaveOut initialized.";

                for (int ii = 0; ii < L_samples; ii++)
                //for (int ii = index1; ii < L_samples + index1; ii++)
                {
                    var bytes = BitConverter.GetBytes((int)(newaudio[ii] * Math.Pow(2, nbytes * 8 - 1)));
                    for (int jj = 0; jj < nbytes; jj++)
                    {
                        //try
                        //{
                        buffer[ii * nbytes + jj] = bytes[jj];
                        //}
                        /*catch (Exception e)
                        { 
                            exc = e.Message; 
                        }*/


                    }
                }
                bufferedWaveProvider.AddSamples(buffer, 0, L_bytes);
                message += "Added " + L_samples + " samples at " + index1_time + " ms";

                

                blockcount_onlineplayback++;
                

                message += "Blockcount = " + blockcount_onlineplayback.ToString();
                message += "Timestamp = " + timestamp_onlineplayback.ToString();
                return message;
            }
            catch (Exception e)
            {
                return e.Message;
            }

        }

        public void InitAudio()
        {
            

            int samplerate = 50000;
            int nbytes = 2;
            bufferedWaveProvider = new NAudio.Wave.BufferedWaveProvider(new NAudio.Wave.WaveFormat(samplerate, nbytes * 8, 1))
            {
                BufferLength = samplerate * nbytes * 10,
                DiscardOnBufferOverflow = false
            };
            waveOut = new NAudio.Wave.WaveOut();
            waveOut.Init(bufferedWaveProvider);
            waveOut.Play();
        }


        public void DisposeAudio()
        {
            waveOut.Stop();
            waveOut.Dispose();
            
        }




        /*private void CreateWaveProvider(double timewindow, int fs)
        {
            byte[] bufferarray = new byte[(int)(fs * timewindow) * 3];
            for (int ii = 0; ii < bufferarray.Length; ii++) bufferarray[0] = 0;
            provider = new NAudio.Wave.RawSourceWaveStream(new System.IO.MemoryStream(bufferarray,true), new NAudio.Wave.WaveFormat(fs, 16, 1));
        }*/


        public IDataSet PAKDZF2RPM(IDataSet rawpulses, int pulses)
        {
            //int L = (int)(t_end - t_start)*fs;

            var IOClass = new DataIOFunctions();
            var rawpulses_double = IOClass.ReadConcertoYData(rawpulses);
            
            //var pulseindices = new int[rawpulses_double.Length];
            //for (int ii = 0; ii < pulseindices.Length; ii++) pulseindices[ii] = (int)(rawpulses_double[ii] * fs);
            var rpm_x = new double[rawpulses_double.Length - 1];
            var rpm_y = new double[rawpulses_double.Length - 1];

            for (int ii = 0; ii < rpm_x.Length; ii++)
            {
                rpm_x[ii] = rawpulses_double[ii + 1];
                rpm_y[ii] = 60 / (double)pulses / (rawpulses_double[ii + 1] - rawpulses_double[ii]);
            }
            rpm_y = GeneralMath.Moving(rpm_y, pulses);
            var x_red = new double[rpm_x.Length - (int)Math.Ceiling((double)pulses / 2)];
            var y_red = new double[x_red.Length];
            Array.Copy(rpm_x, 0, x_red, 0, x_red.Length);
            Array.Copy(rpm_y, 0, y_red, 0, y_red.Length);

            rpm_y = Interp1_linear(x_red, y_red, rpm_x, true);
            rawpulses.Release();

            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = rpm_x.Length,
                RealValues = rpm_y,
                x = new IDataSet(DataSetType.Numeric)
                {
                    Units = "s"
                }
            };
            result.x.RealValues = rpm_x;

            return result;
        }

        public IDataSet RawPulses2RPM(IDataSet rawpulses, int ms)
        {
            double timeunit;
            if (ms == 1) timeunit = 1000;
            else timeunit = 1;


            var IOClass = new DataIOFunctions();
            var rawpulses_double = IOClass.ReadConcertoYData(rawpulses);
            int fs_rpm = (int)(timeunit / (rawpulses.x[1] - rawpulses.x[0]));
            var pulseindices = FindPulses(rawpulses_double, fs_rpm);
            var rpm_x = new double[pulseindices.Length - 1];
            var rpm_y = new double[pulseindices.Length - 1];

            for (int ii = 0; ii < rpm_x.Length; ii++)
            {
                rpm_x[ii] = rawpulses.x[pulseindices[ii + 1]] / timeunit;
                rpm_y[ii] = 60 * (double)fs_rpm / (pulseindices[ii + 1] - pulseindices[ii]);
            }
            rawpulses.Release();

            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = rpm_x.Length,
                RealValues = rpm_y,
                x = new IDataSet(DataSetType.Numeric)
            };
            result.x.RealValues = rpm_x;

            return result;
        }

        public void ReadCyclesInit(int count, int cyclecount, int arraymode, int lastcyclesingular)
        {
            buffer1 = new double[count];
            if (arraymode == 1)
            {
                cadatainputarray = new NVHFunctions.CycleData[cyclecount];
                for (int ii = 0; ii < cyclecount - 1; ii++)
                {
                    cadatainputarray[ii].Count = count;
                    cadatainputarray[ii].CycleCount = 2;
                    cadatainputarray[ii].ydata = new double[count * 2];
                }
                if (lastcyclesingular == 1)
                {
                    cadatainputarray[cadatainputarray.Length - 1].Count = count;
                    cadatainputarray[cadatainputarray.Length - 1].CycleCount = 1;
                    cadatainputarray[cadatainputarray.Length - 1].ydata = new double[count];
                }
                else
                {
                    cadatainputarray[cadatainputarray.Length - 1].Count = count;
                    cadatainputarray[cadatainputarray.Length - 1].CycleCount = 2;
                    cadatainputarray[cadatainputarray.Length - 1].ydata = new double[count * 2];
                }
            }
            else
            {
                cadatainput = new NVHFunctions.CycleData
                {
                    Count = count,
                    CycleCount = cyclecount,
                    xdata = new double[count],
                    ydata = new double[count * cyclecount]
                };
            }

        }



        public void ReadOneCycle(IDataSet data, int arraymode, int arrayindex, int cycleindex, int cyclecount)
        {
            if (arraymode == 1)
            {
                cadatainputarray[arrayindex].xdata = data.x.RealValues;
                buffer1 = data.RealValues;
                if (cycleindex == 0) for (int ii = 0; ii < data.Count; ii++) cadatainputarray[arrayindex].ydata[ii] = buffer1[ii];
                else for (int ii = 0; ii < data.Count; ii++) cadatainputarray[arrayindex].ydata[ii + data.Count] = buffer1[ii];
            }
            else
            {
                if (cycleindex == 0)
                {
                    cadatainput.Count = data.Count;
                    cadatainput.CycleCount = cyclecount;
                    cadatainput.xdata = data.x.RealValues;
                }
                //int kk = cycleindex - cyclestart;
                buffer1 = data.RealValues;
                for (int jj = 0; jj < data.Count; jj++) cadatainput.ydata[cycleindex * data.Count + jj] = buffer1[jj];
            }


        }


        internal static CycleData[] ReadCADataIntoCycleDataArray(IDataSet CAData, CATimeStamps timestamps)
        {
            CycleData[] CAData_array;
            CAData_array = new CycleData[timestamps.cycles.Length];
            double[] onecycle;
            for (int ii = 0; ii < timestamps.cycles.Length; ii++)
            {
                if (timestamps.cycles[ii] < CAData.CycleCount - 1)
                {
                    CAData_array[ii].CycleCount = 2;
                    CAData_array[ii].ydata = new double[CAData.Count * 2];
                    CAData.CycleIndex = (int)timestamps.cycles[ii];
                    onecycle = CAData.RealValues;
                    for (int jj = 0; jj < onecycle.Length; jj++) CAData_array[ii].ydata[jj] = onecycle[jj];
                    CAData.CycleIndex = (int)timestamps.cycles[ii] + 1;
                    onecycle = CAData.RealValues;
                    for (int jj = 0; jj < onecycle.Length; jj++) CAData_array[ii].ydata[jj + onecycle.Length] = onecycle[jj];
                }
                else
                {
                    CAData_array[ii].CycleCount = 1;
                    CAData_array[ii].ydata = new double[CAData.Count];
                    CAData.CycleIndex = (int)timestamps.cycles[ii];
                    CAData_array[ii].ydata = CAData.RealValues;
                }
                CAData_array[ii].xdata = CAData.x.RealValues;
                CAData_array[ii].Count = CAData.Count;
            }
            CAData.Release();
            return CAData_array;
        }

        public IDataSet ResampleCASignal_dT(IDataSet CAData, IDataSet dT, double offset, double fs_new, int cyclestart, int cycleend, int ms)
        {
            // obsolete - done within dll

            // Units: 
            // CAData(x): deg
            // dT(y): us/deg
            // offset: sec
            // fs_new: samples/sec
            // cyclestart: first cycle (first index is 1)
            // cycleend: last cycle (last index is DS.Count)
            // ms: 1: result in ms, 2: result in sec


            var IOClass = new DataIOFunctions();

            double time_factor;
            if (ms == 1) time_factor = 1000;
            else time_factor = 1;

            var raw_d = IOClass.ReadConcertoCAData(CAData, cyclestart - 1, cycleend - 1);
            var dT_d = IOClass.ReadConcertoCAData(dT, cyclestart - 1, cycleend - 1);

            CAData.Release();
            dT.Release();


            var ydata = ResampleCAData_dT(raw_d, dT_d, fs_new);
            var xdata = new double[ydata.Length];
            
            for (int ii = 0; ii < xdata.Length; ii++) xdata[ii] = (offset + ii / fs_new) * time_factor;
            var xdata_IDS = new IDataSet(DataSetType.Numeric)
            {
                Count = xdata.Length,
                RealValues = xdata
            };
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = ydata.Length,
                x = xdata_IDS,
                RealValues = ydata

            };


            return result;
        }

        /* private void ResampleTMData(IDataSet rawdata, double fs)
         {
             using (var br = new BinaryWriter(new FileStream("C:/Pd/test.bin", FileMode.Create)))
             {
                 var doubles = rawdata.RealValues;
                 double x0 = rawdata.x[0];
                 double x1 = rawdata.x[1];
                 double dt = x1 - x0;
                 for (int ii = 0; ii < doubles.Length; ii++) br.Write(doubles[ii]);
                 br.Write(x0);
                 br.Write(dt);
                 br.Write(fs);
                 br.Close();
             }

         }*/


        public IDataSet ResampleCASignal_CDMTIME(IDataSet CAData, IDataSet cdmtime, double fs_new, int cyclestart, int cycleend, int ms)
        {
            // obsolete - done within dll

            // Units: 
            // CAData(x): deg
            // cdmtime(y): ms
            // fs_new: samples/sec
            // cyclestart: first cycle (first index is 1)
            // cycleend: last cycle (last index is DS.Count)
            // ms: 1: result in ms, 2: result in sec


            var IOClass = new DataIOFunctions();
            
            var raw_d = IOClass.ReadConcertoCAData(CAData, cyclestart - 1, cycleend - 1);
            var cdmtime_d = IOClass.ReadConcertoCAData(cdmtime, cyclestart - 1, cycleend - 1);

            CAData.Release();
            cdmtime.Release();

            double time_factor;
            if (ms == 1) time_factor = 1000;
            else time_factor = 1;

            var ydata = ResampleCAData_CDMTIME(raw_d, cdmtime_d, fs_new);
            var xdata = new double[ydata.Length];
            for (int ii = 0; ii < xdata.Length; ii++) xdata[ii] = (cdmtime_d.ydata[0]/1000 + ii / fs_new) * time_factor;
            
            var xdata_IDS = new IDataSet(DataSetType.Numeric)
            {
                Count = xdata.Length,
                RealValues = xdata
            };
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = ydata.Length,
                x = xdata_IDS,
                RealValues = ydata
                
            };
            
            
            return result;
        }

        public void ResumePlayback()
        {
            if (waveOut != null) if (waveOut.PlaybackState == NAudio.Wave.PlaybackState.Paused) waveOut.Play();
        }

        public void StopPlayback()
        {
            if (waveOut != null) if (waveOut.PlaybackState == NAudio.Wave.PlaybackState.Paused || waveOut.PlaybackState == NAudio.Wave.PlaybackState.Playing)
                {
                    waveOut.Stop();
                    waveOut.Dispose();
                }

        }


        public double Tictimer(int startstop)
        {
            double returnvalue;

            if (startstop == 1)
            {
                timer.Start();
                returnvalue = 0;
            }
            else
            {
                returnvalue = timer.ElapsedMilliseconds;
                timer.Restart();
            }
            return returnvalue;
        }

        public double Return2()
        {
            /*var input_d = data.RealValues;
            for (int ii = 0; ii < input_d.Length; ii++) input_d[ii] *= 2;
            var result = new IDataSet(DataSetType.Numeric)
            {
                Count = input_d.Length,
                RealValues = input_d,
                x = data.x
            };
            return result;*/
            return 2.0;
        }
        public double Return4()
        {
            return 4.0;
        }

        public string WaveExport(string filename, IDataSet audio_l, IDataSet audio_r, IDataSet rpm, double ms, double dbref, double calib_l, double calib_r)
        {
            string errormessage = "Audio data export successful.";
            double timeunit;
            if (ms == 1) timeunit = 1000;
            else timeunit = 1;

            //var FcnClass = new NVHFunctions();
            var IOClass = new DataIOFunctions();
            

            var rawdata_l = IOClass.ReadConcertoYData(audio_l);
            int fs_audio = (int)(timeunit / (audio_l.x[1] - audio_l.x[0]));

            audio_l.Release();

            var rpm_xy = new XY_Data();
            double[] rawdata_r = null;

            if (audio_r.Count != 0) rawdata_r = IOClass.ReadConcertoYData(audio_r);
            audio_r.Release();
            if (rpm.Count != 0)
            {
                rpm_xy.xdata = rpm.x.RealValues;
                rpm_xy.ydata = rpm.RealValues;
            }
            rpm.Release();
            if (dbref == -1) dbref = 0.00002;

            if (rpm_xy.xdata != null) for (int ii = 0; ii < rpm_xy.xdata.Length; ii++) rpm_xy.xdata[ii] /= timeunit;

            WaveStream.WriteAVLWaveFile(rawdata_l, fs_audio, filename, audio_r: rawdata_r, rpm_l: rpm_xy, calib_l: calib_l, calib_r: calib_r, dbref: dbref);


            return errormessage;
        }





        




        


        

        



        
        
    }
	

}
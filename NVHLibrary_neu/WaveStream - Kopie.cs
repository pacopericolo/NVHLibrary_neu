using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Collections;
using NVHLibrary;
using System.ComponentModel;

namespace WaveFunctions
{

    public class WaveStream
    {
        /*public struct XYData
        {
            public double[] xdata;
            public double[] ydata;
        }*/

        public struct WaveFileObject
        {
            public int[] audio_L;
            public int[] audio_R;
            public int speedchannel;   // 1: left, 2: right, 3: both
            public int nchannels;
            public int samplerate;
            public int nbits;
            public double calib_l;     // Digits/Unit
            public double calib_r;     // Digits/Unit
            public int speedpulses_l;
            public int speedpulses_r;
            public double speedfactor_l;
            public double speedfactor_r;
            public string quantity_l;
            public string quantity_r;
            public string description_l;
            public string description_r;
            public string unit_l;
            public string unit_r;
            public int max_l;
            public int min_l;
            public int max_r;
            public int min_r;
            
            public double maxspeed_l;
            public double minspeed_l;
            public double maxspeed_r;
            public double minspeed_r;
            public NVHFunctions.XY_Data rpm_l;
            public NVHFunctions.XY_Data rpm_r;
        }

        public struct HeaderObject
        {
            public int speedchannel;   // 1: left, 2: right, 3: both
            public int nchannels;
            public int samplerate;
            public int nbits;
            public double calib_l;     // Digits/Unit
            public double calib_r;     // Digits/Unit
            public int speedpulses_l;
            public int speedpulses_r;
            public double speedfactor_l;
            public double speedfactor_r;
            public string quantity_l;
            public string quantity_r;
            public string description_l;
            public string description_r;
            public string unit_l;
            public string unit_r;
            public long L;
        }

        public struct DataObject
        {
            public int[] audio_L;
            public int[] audio_R;
            public int max_l;
            public int min_l;
            public int max_r;
            public int min_r;
            public double maxspeed_l;
            public double minspeed_l;
            public double maxspeed_r;
            public double minspeed_r;
            public NVHFunctions.XY_Data rpm_l;
            public NVHFunctions.XY_Data rpm_r;
        }

        public HeaderObject ReadWaveHeader(string filename)
        {
            var output = new HeaderObject
            {
                
                // set default values

                calib_l = -1,
                calib_r = -1,

                speedpulses_l = -1,
                speedpulses_r = -1,

                speedfactor_l = 1,
                speedfactor_r = 1,


                unit_l = "",
                unit_r = ""
            };

            BinaryReader br;

            //var FcnClass = new NVHFunctions();
            try
            {
                using (br = new BinaryReader(new FileStream(filename, FileMode.Open), Encoding.ASCII))
                {
                    // read Riff chunk

                    var riffchunk = ReadByteString(br, 4);
                    uint riffchunksize = br.ReadUInt32();   // Filelength in bytes - 8 ("RIFF" "WAVE" not included)
                    var wavechunk = ReadByteString(br, 4);

                    // read "fmt " chunk

                    var fmtchunk = ReadByteString(br, 4);    // = "fmt "
                    uint fmtchunksize = br.ReadUInt32();            // normally 16 bytes

                    int wFormatTag = br.ReadByte();                 // 1 for PCM
                    br.ReadByte();
                    output.nchannels = br.ReadByte();                  // 2 for stereo
                    br.ReadByte();

                    output.samplerate = (int)br.ReadUInt32();          // = sample rate
                    uint nAvgBytesPerSec = br.ReadUInt32();         // = nChannels * nSamplesPerSec * nBitsPerSample / 8

                    int nBlockAlign = br.ReadByte();                // = nChannels * nBitsPerSample / 8
                    br.ReadByte();
                    output.nbits = br.ReadByte();             // = Bit Depth
                    br.ReadByte();



                    // read AVL/Comfort Info

                    var avlchunk = ReadByteString(br, 4);
                    if (String.Equals(avlchunk, "AVLI"))
                    {
                        uint avlchunksize = br.ReadUInt32();
                        ReadByteString(br, 4);                              // 4 "fill" bytes (Delphi Field Alignment)             
                        output.calib_l = br.ReadDouble();                    // calibration value Digits/Unit
                        output.calib_r = br.ReadDouble();                    // calibration value Digits/Unit
                        output.speedchannel = br.ReadByte();               // 1: left channel, 2: right channel
                        ReadByteString(br, 3);                              // 3 "fill" bytes (Delphi Field Alignment)
                        if (output.speedchannel == 1)
                        {
                            output.speedpulses_l = (int)br.ReadUInt32();
                            output.speedfactor_l = br.ReadDouble();
                        }
                        else
                        {
                            output.speedpulses_r = (int)br.ReadUInt32();   // Speed Pulses
                            output.speedfactor_r = br.ReadDouble();    // Speed Factor
                        }

                        var datachunk = ReadByteString(br, 4);   // "data"
                        uint datachunksize = br.ReadUInt32();
                        output.L = datachunksize / output.nbits * 8 / output.nchannels;
                    }
                    else if (String.Equals(avlchunk, "COMI"))
                    {
                        uint avlchunksize = br.ReadUInt32();
                        output.calib_l = br.ReadDouble();                    // calibration value Digits/Unit
                        output.calib_r = br.ReadDouble();                    // calibration value Digits/Unit
                        var speedchannel_L = br.ReadByte();
                        output.speedpulses_l = (int)br.ReadUInt32();             // Speed Pulses
                        output.speedfactor_l = br.ReadDouble();
                        output.quantity_l = ReadByteString(br, 3);
                        output.description_l = ReadByteString(br, 32);
                        output.unit_l = ReadByteString(br, 8);

                        var speedchannel_R = br.ReadByte();
                        output.speedpulses_r = (int)br.ReadUInt32();             // Speed Pulses
                        output.speedfactor_r = br.ReadDouble();
                        output.quantity_r = ReadByteString(br, 3);
                        output.description_r = ReadByteString(br, 32);
                        output.unit_r = ReadByteString(br, 8);

                        var datachunk = ReadByteString(br, 4);   // "data"
                        uint datachunksize = br.ReadUInt32();
                        output.L = datachunksize / output.nbits * 8 / output.nchannels;
                    }
                    br.Close();
                }
            }
            catch
            {

            }
            
            

            return output;
        }

        public double HMS2UnitsPerPascal(double hms, int nbits)
        {
            double unitsperpascal = Math.Pow(2, (nbits - 1)) / (Math.Pow(10, hms / 20) * 0.00002 * 2 * Math.Sqrt(2));
            return unitsperpascal;
        }

        public double UnitsPerPascal2HMS(double unitsperpascal, int nbits)
        {
            double hms = 20 * Math.Log10(Math.Pow(2, nbits - 1) / unitsperpascal / 0.00002 / 2 / Math.Sqrt(2));
            return hms;
        }

        private bool[] EliminateDoublePulses(bool[] input)
        {
            for (int ii = 0; ii < input.Length; ii++)
            {
                if(input[ii] && ii > 0)
                {
                    if (input[ii - 1]) input[ii] = false;
                }
            }
            return input;
        }

        public int[] GetNativeIntValuesFromPascalValues(double[] audiodouble, double hmscalib, int nbits)
        {
            var audioint = new int[audiodouble.Length];

            double unitsperpascal = HMS2UnitsPerPascal(hmscalib, nbits);
            for (int ii = 0; ii < audioint.Length; ii++) audioint[ii] = (int)(audiodouble[ii] * unitsperpascal);
            return audioint;
        }

        public DataObject ReadWaveData(string filename, HeaderObject input, long fromSample = -1, long toSample = -1, BackgroundWorker worker = null, double startprogress = 0, double overallprogress = 100)
        {
            var output = new DataObject();
            
            BinaryReader br;
            

            using (br = new BinaryReader(new FileStream(filename, FileMode.Open)))
            {
                // read Riff chunk
                br.BaseStream.Position = 36;

                var avlchunk = ReadByteString(br, 4);
                uint avlchunksize;
                if (String.Equals(avlchunk, "AVLI") || String.Equals(avlchunk, "COMI"))
                {
                    avlchunksize = br.ReadUInt32();
                    br.BaseStream.Position = 40 + avlchunksize + 8;
                }


                // Read Data Chunk

                uint datachunksize = br.ReadUInt32();           // size of data chunk in bytes

                if (input.nchannels > 2 || input.nchannels < 1) return output;

                output.audio_L = new int[datachunksize / input.nbits * 8 / input.nchannels];
                if (input.nchannels == 2) output.audio_R = new int[datachunksize / input.nbits * 8 / input.nchannels];

                bool[] pulses_L, pulses_R;

                pulses_L = new bool[output.audio_L.Length];
                if (input.nchannels == 2) pulses_R = new bool[output.audio_R.Length];
                else pulses_R = new bool[0];

                output.max_l = 0;
                output.min_l = 0;
                output.max_r = 0;
                output.min_r = 0;

                output.maxspeed_l = 0;
                output.minspeed_l = 0;
                output.maxspeed_r = 0;
                output.minspeed_r = 0;

                for (int ii = 0; ii < output.audio_L.Length; ii++)
                {
                    if (input.nbits == 24) output.audio_L[ii] = Read24BitSample(br);
                    else output.audio_L[ii] = br.ReadInt16();
                    if (output.audio_L[ii] > output.max_l) output.max_l = output.audio_L[ii];
                    if (output.audio_L[ii] < output.min_l) output.min_l = output.audio_L[ii];
                    var bit = (output.audio_L[ii] & (1 << input.nbits - 16)) != 0;
                    if (bit) pulses_L[ii] = true;


                    if (input.nchannels == 2)
                    {
                        if (input.nbits == 24) output.audio_R[ii] = Read24BitSample(br);
                        else output.audio_R[ii] = br.ReadInt16();
                        if (output.audio_R[ii] > output.max_r) output.max_r = output.audio_R[ii];
                        if (output.audio_R[ii] < output.min_r) output.min_r = output.audio_R[ii];

                        bit = (output.audio_R[ii] & (1 << input.nbits - 16)) != 0;
                        if (bit) pulses_R[ii] = true;
                    }
                    if (worker != null)
                    {
                        if (ii % 4 == 0) worker.ReportProgress((int)(startprogress + overallprogress / 2 * ii * 4 / output.audio_L.Length));
                    }
                    
                }


                br.Close();

                if (input.speedpulses_l != -1)
                {
                    // int L = Array.FindAll(pulses_L, x => x == true).Length;
                    
                    int L = 0;
                    pulses_L = EliminateDoublePulses(pulses_L);
                    for (int ii = 0; ii < pulses_L.Length; ii++) if (pulses_L[ii]) L++;

                    var indices = new int[(int)Math.Ceiling(L / (double)input.speedpulses_l)];
                    int jj = 0;
                    int kk = 0;
                    
                    for (int ii = 0; ii < pulses_L.Length; ii++)
                    {
                        if (pulses_L[ii])
                        {
                            if (kk == input.speedpulses_l) kk = 0;
                            if (kk == 0)
                            {
                                indices[jj] = ii;
                                jj++;
                            }
                            kk++;
                        }
                    }
                    if (indices.Length > 0)
                    {

                        var distances = new int[indices.Length - 1];
                        for (int ii = 0; ii < distances.Length; ii++) distances[ii] = indices[ii + 1] - indices[ii];
                        int distanceisone = Array.FindAll(distances, x => x == 1).Length;
                        var indices_red = new int[indices.Length - 1 - distanceisone];
                        var distances_red = new int[indices_red.Length];
                        jj = 0;
                        for (int ii = 1; ii < indices.Length; ii++)
                        {
                            if (distances[ii - 1] != 1)
                            {
                                indices_red[jj] = indices[ii];
                                distances_red[jj] = distances[ii - 1];

                                jj++;
                            }
                        }


                        double unitconversion = 1;
                        if (input.unit_l.Length > 2) if (String.Equals("m/s", input.unit_l.Substring(0, 3))) unitconversion = 16.667;


                        output.rpm_l = new NVHFunctions.XY_Data
                        {
                            xdata = new double[indices_red.Length],
                            ydata = new double[indices_red.Length]

                        };

                        for (int ii = 0; ii < indices_red.Length; ii++)
                        {

                            output.rpm_l.xdata[ii] = indices_red[ii] / (double)input.samplerate - (double)distances_red[0]/2/input.samplerate;
                            output.rpm_l.ydata[ii] = input.samplerate * 60 / distances_red[ii] / input.speedfactor_l / unitconversion;
                            if (output.rpm_l.ydata[ii] > output.maxspeed_l) output.maxspeed_l = output.rpm_l.ydata[ii];
                            if (output.rpm_l.ydata[ii] < output.minspeed_l) output.minspeed_l = output.rpm_l.ydata[ii];
                        }
                    }
                    else
                    {

                        output.rpm_l.xdata = null;
                        output.rpm_l.ydata = null;
                    }

                }
                if (worker != null) worker.ReportProgress((int)(startprogress + overallprogress / 2 + overallprogress / 4));

                if (input.speedpulses_r != -1)
                {
                    // int L = Array.FindAll(pulses_L, x => x == true).Length;
                    int L = 0;
                    pulses_R = EliminateDoublePulses(pulses_R);
                    for (int ii = 0; ii < pulses_R.Length; ii++) if (pulses_R[ii]) L++;

                    var indices = new int[(int)Math.Ceiling(L / (double)input.speedpulses_r)];
                    int jj = 0;
                    int kk = 0;
                    
                    for (int ii = 0; ii < pulses_R.Length; ii++)
                    {
                        if (pulses_R[ii])
                        {
                            if (kk == input.speedpulses_r) kk = 0;
                            if (kk == 0)
                            {
                                indices[jj] = ii;
                                jj++;
                            }
                            kk++;
                        }
                    }
                    if (indices.Length > 0)
                    {

                        var distances = new int[indices.Length - 1];
                        for (int ii = 0; ii < distances.Length; ii++) distances[ii] = indices[ii + 1] - indices[ii];
                        int distanceisone = Array.FindAll(distances, x => x == 1).Length;
                        var indices_red = new int[indices.Length - 1 - distanceisone];
                        var distances_red = new int[indices_red.Length];
                        jj = 0;
                        for (int ii = 1; ii < indices.Length; ii++)
                        {
                            if (distances[ii - 1] != 1)
                            {
                                indices_red[jj] = indices[ii];
                                distances_red[jj] = distances[ii - 1];

                                jj++;
                            }
                        }


                        double unitconversion = 1;
                        if (input.unit_r.Length > 2) if (String.Equals("m/s", input.unit_r.Substring(0, 3))) unitconversion = 16.667;


                        output.rpm_r = new NVHFunctions.XY_Data
                        {
                            xdata = new double[indices_red.Length],
                            ydata = new double[indices_red.Length]

                        };

                        for (int ii = 0; ii < indices_red.Length; ii++)

                        {

                            output.rpm_r.xdata[ii] = indices_red[ii] / (double)input.samplerate;
                            output.rpm_r.ydata[ii] = input.samplerate * 60 / distances_red[ii] / input.speedfactor_r / unitconversion;
                            if (output.rpm_r.ydata[ii] > output.maxspeed_r) output.maxspeed_r = output.rpm_r.ydata[ii];
                            if (output.rpm_r.ydata[ii] < output.minspeed_r) output.minspeed_r = output.rpm_r.ydata[ii];
                        }
                    }
                    else
                    {
                        output.rpm_r.xdata = null;
                        output.rpm_r.ydata = null;
                    }

                }
                if (worker != null) worker.ReportProgress((int)(startprogress + overallprogress));
            }
            
            
            return output;
        }

        public bool[][] ReadRPMPulses(string filename, HeaderObject input, BackgroundWorker worker = null, double startprogress = 0, double overallprogress = 100)
        {
            var output = new bool[2][];

            BinaryReader br;


            using (br = new BinaryReader(new FileStream(filename, FileMode.Open)))
            {
                // read Riff chunk
                br.BaseStream.Position = 36;

                var avlchunk = ReadByteString(br, 4);
                uint avlchunksize;
                if (String.Equals(avlchunk, "AVLI") || String.Equals(avlchunk, "COMI"))
                {
                    avlchunksize = br.ReadUInt32();
                    br.BaseStream.Position = 40 + avlchunksize + 8;
                }


                // Read Data Chunk

                uint datachunksize = br.ReadUInt32();           // size of data chunk in bytes

                if (input.nchannels > 2 || input.nchannels < 1) return output;
                int[] audio_R = new int[0];
                var audio_L = new int[datachunksize / input.nbits * 8 / input.nchannels];
                if (input.nchannels == 2) audio_R = new int[datachunksize / input.nbits * 8 / input.nchannels];

                bool[] pulses_L, pulses_R;

                pulses_L = new bool[audio_L.Length];
                if (input.nchannels == 2) pulses_R = new bool[audio_R.Length];
                else pulses_R = new bool[0];
                
                for (int ii = 0; ii < audio_L.Length; ii++)
                {
                    if (input.nbits == 24) audio_L[ii] = Read24BitSample(br);
                    else audio_L[ii] = br.ReadInt16();
                    
                    
                    var bit = (audio_L[ii] & (1 << input.nbits - 16)) != 0;
                    if (bit) pulses_L[ii] = true;


                    if (input.nchannels == 2)
                    {
                        if (input.nbits == 24) audio_R[ii] = Read24BitSample(br);
                        else audio_R[ii] = br.ReadInt16();
                        
                        bit = (audio_R[ii] & (1 << input.nbits - 16)) != 0;
                        if (bit) pulses_R[ii] = true;
                    }
                    if (worker != null)
                    {
                        if (ii % 4 == 0) worker.ReportProgress((int)(startprogress + overallprogress / 2 * ii * 4 / audio_L.Length));
                    }

                }


                br.Close();


                var pulses = new bool[2][];
                if (input.speedpulses_l != -1)
                {
                    // int L = Array.FindAll(pulses_L, x => x == true).Length;

                   // int L = 0;
                    output[0] = EliminateDoublePulses(pulses_L);
                    

                }
                
                if (input.speedpulses_r != -1)
                {
                    // int L = Array.FindAll(pulses_L, x => x == true).Length;
                    
                    output[1] = EliminateDoublePulses(pulses_R);
                    

                }
                
            }


            return output;
        }


        public WaveFileObject ReadWaveFile(string filename, double calib_l = -1, double calib_r = -1, int speedpulses_l = -1, int speedpulses_r = -1, double speedfactor_l = -1, double speedfactor_r = -1, double dbref = 0.00002)
        {
            var output = new WaveFileObject();
            //string inputfile = "Runup_GolfVII.wav";
            //string inputfile = "Runup_TeslaX.wav";
            //string inputfile = "test.wav";
            BinaryReader br;
            bool calib_l_defined = false;
            bool calib_r_defined = false;
            if (calib_l != -1) calib_l_defined = true;
            if (calib_r != -1) calib_r_defined = true;
            //var FcnClass = new NVHFunctions();
            br = new BinaryReader(new FileStream(filename, FileMode.Open),Encoding.ASCII);
            

            // read Riff chunk

            var riffchunk = ReadByteString(br, 4);
            uint riffchunksize = br.ReadUInt32();   // Filelength in bytes - 8 ("RIFF" "WAVE" not included)
            var wavechunk = ReadByteString(br, 4);

            // read "fmt " chunk

            var fmtchunk = ReadByteString(br, 4);    // = "fmt "
            uint fmtchunksize = br.ReadUInt32();            // normally 16 bytes

            int wFormatTag = br.ReadByte();                 // 1 for PCM
            br.ReadByte();
            output.nchannels = br.ReadByte();                  // 2 for stereo
            br.ReadByte();

            output.samplerate = (int)br.ReadUInt32();          // = sample rate
            uint nAvgBytesPerSec = br.ReadUInt32();         // = nChannels * nSamplesPerSec * nBitsPerSample / 8

            int nBlockAlign = br.ReadByte();                // = nChannels * nBitsPerSample / 8
            br.ReadByte();
            output.nbits = br.ReadByte();             // = Bit Depth
            br.ReadByte();

            // set default values

            /*if (calib_l == -1) output.calib_l = Math.Pow(2, output.nbits - 1);
            else output.calib_l = Math.Pow(2, output.nbits - 1) / (Math.Pow(10, calib_l / 20) * dbref * 2 * Math.Sqrt(2));
            if (calib_r == -1) output.calib_r = Math.Pow(2, output.nbits - 1);
            else output.calib_r = Math.Pow(2, output.nbits - 1) / (Math.Pow(10, calib_r / 20) * dbref * 2 * Math.Sqrt(2));*/
            
            if (speedpulses_l == -1) output.speedpulses_l = 0;
            else output.speedpulses_l = speedpulses_l;
            if (speedpulses_r == -1) output.speedpulses_r = 0;
            else output.speedpulses_r = speedpulses_r;
            if (speedfactor_l == -1) output.speedfactor_l = 0;
            else output.speedfactor_l = speedfactor_l;
            if (speedfactor_r == -1) output.speedfactor_r = 0;
            else output.speedfactor_r = speedfactor_r;

            output.unit_l = "";
            output.unit_r = "";

            // read AVL/Comfort Info

            var avlchunk = ReadByteString(br, 4);
            if (String.Equals(avlchunk, "AVLI"))
            {
                uint avlchunksize = br.ReadUInt32();
                ReadByteString(br, 4);                              // 4 "fill" bytes (Delphi Field Alignment)             
                if (calib_l_defined)
                {
                    br.ReadDouble();
                    output.calib_l = Math.Pow(2, output.nbits - 1) / (Math.Pow(10, calib_l / 20) * dbref * 2 * Math.Sqrt(2));
                }
                else output.calib_l = br.ReadDouble();                    // calibration value in digits/Unit
                if (calib_r_defined)
                {
                    br.ReadDouble();
                    output.calib_r = Math.Pow(2, output.nbits - 1) / (Math.Pow(10, calib_r / 20) * dbref * 2 * Math.Sqrt(2));
                }
                else output.calib_r = br.ReadDouble();                    // calibration value in digits/Unit
                output.speedchannel = br.ReadByte();               // 1: left channel, 2: right channel
                ReadByteString(br, 3);                              // 3 "fill" bytes (Delphi Field Alignment)
                if (output.speedchannel == 1)       
                {
                    output.speedpulses_l = (int)br.ReadUInt32();
                    output.speedfactor_l = br.ReadDouble();
                }
                else
                {
                    output.speedpulses_r = (int)br.ReadUInt32();   // Speed Pulses
                    output.speedfactor_r = br.ReadDouble();    // Speed Factor
                }

                var datachunk = ReadByteString(br, 4);   // "data"
            }
            else if (String.Equals(avlchunk, "COMI"))
            {
                uint avlchunksize = br.ReadUInt32();
                if (calib_l_defined)
                {
                    br.ReadDouble();
                    output.calib_l = Math.Pow(2, output.nbits - 1) / (Math.Pow(10, calib_l / 20) * dbref * 2 * Math.Sqrt(2));
                }
                else
                {
                    output.calib_l = br.ReadDouble();                    // calibration value in digits/Unit
                    //output.calib_l = 20 * Math.Log10(Math.Pow(2, output.nbits - 1) / output.calib_l / dbref / 2 / Math.Sqrt(2));
                }

                if (calib_r_defined)
                {
                    br.ReadDouble();
                    output.calib_r = Math.Pow(2, output.nbits - 1) / (Math.Pow(10, calib_r / 20) * dbref * 2 * Math.Sqrt(2));
                }
                else
                {
                    output.calib_r = br.ReadDouble();                    // calibration value in digits/Unit
                    //output.calib_r = 20 * Math.Log10(Math.Pow(2, output.nbits - 1) / output.calib_r / dbref / 2 / Math.Sqrt(2));
                }
                var speedchannel_L = br.ReadByte();
                output.speedpulses_l = (int)br.ReadUInt32();             // Speed Pulses
                output.speedfactor_l = br.ReadDouble();
                output.quantity_l = ReadByteString(br, 3);
                output.description_l = ReadByteString(br, 32);
                output.unit_l = ReadByteString(br, 8);

                var speedchannel_R = br.ReadByte();
                output.speedpulses_r = (int)br.ReadUInt32();             // Speed Pulses
                output.speedfactor_r = br.ReadDouble();
                output.quantity_r = ReadByteString(br, 3);
                output.description_r = ReadByteString(br, 32);
                output.unit_r = ReadByteString(br, 8);

                var datachunk = ReadByteString(br, 4);   // "data"
            }

            // Read Data Chunk

            uint datachunksize = br.ReadUInt32();           // size of data chunk in bytes

            if (output.nchannels > 2 || output.nchannels < 1) return output;

            output.audio_L = new int[datachunksize / output.nbits * 8 / output.nchannels];
            if (output.nchannels == 2) output.audio_R = new int[datachunksize / output.nbits * 8 / output.nchannels];

            bool[] pulses_L, pulses_R;

            pulses_L = new bool[output.audio_L.Length];
            if (output.nchannels == 2) pulses_R = new bool[output.audio_R.Length];
            else pulses_R = new bool[0];

            output.max_l = 0;
            output.min_l = 0;
            output.max_r = 0;
            output.min_r = 0;

            output.maxspeed_l = 0;
            output.minspeed_l = 0;
            output.maxspeed_r = 0;
            output.minspeed_r = 0;

            for (int ii = 0; ii < output.audio_L.Length; ii++)
            {
                if (output.nbits == 24) output.audio_L[ii] = Read24BitSample(br);
                else if (output.nbits == 32) output.audio_L[ii] = Read32BitSample(br, wFormatTag);
                else output.audio_L[ii] = br.ReadInt16();
                if (output.audio_L[ii] > output.max_l) output.max_l = output.audio_L[ii];
                if (output.audio_L[ii] < output.min_l) output.min_l = output.audio_L[ii];
                var bit = (output.audio_L[ii] & (1 << output.nbits - 16)) != 0;
                if (bit) pulses_L[ii] = true;
                
                if (output.nchannels == 2)
                {
                    if (output.nbits == 24) output.audio_R[ii] = Read24BitSample(br);
                    else if (output.nbits == 32) output.audio_R[ii] = Read32BitSample(br, wFormatTag);
                    else output.audio_R[ii] = br.ReadInt16();
                    if (output.audio_R[ii] > output.max_r) output.max_r = output.audio_R[ii];
                    if (output.audio_R[ii] < output.min_r) output.min_r = output.audio_R[ii];

                    bit = (output.audio_R[ii] & (1 << output.nbits - 16)) != 0;
                    if (bit) pulses_R[ii] = true;
                }
            }


            br.Close();

            if (output.speedpulses_l > 0)
            {
                // int L = Array.FindAll(pulses_L, x => x == true).Length;
                int L = 0;
                for (int ii = 0; ii < pulses_L.Length; ii++) if (pulses_L[ii]) L++;

                var indices = new int[(int)Math.Ceiling(L / (double)output.speedpulses_l)];
                int jj = 0;
                int kk = 0;
                for (int ii = 0; ii < pulses_L.Length; ii++)
                {
                    if (pulses_L[ii])
                    {
                        if (kk == output.speedpulses_l) kk = 0;
                        if (kk == 0)
                        {
                            indices[jj] = ii;
                            jj++;
                        }
                        kk++;
                    }
                }
                if (indices.Length > 0)
                {

                    var distances = new int[indices.Length - 1];
                    for (int ii = 0; ii < distances.Length; ii++) distances[ii] = indices[ii + 1] - indices[ii];
                    int distanceisone = Array.FindAll(distances, x => x == 1).Length;
                    var indices_red = new int[indices.Length - 1 - distanceisone];
                    var distances_red = new int[indices_red.Length];
                    jj = 0;
                    for (int ii = 1; ii < indices.Length; ii++)
                    {
                        if (distances[ii - 1] != 1)
                        {
                            indices_red[jj] = indices[ii];
                            distances_red[jj] = distances[ii - 1];

                            jj++;
                        }
                    }


                    double unitconversion = 1;
                    if (output.unit_l.Length > 2) if (String.Equals("m/s", output.unit_l.Substring(0, 3))) unitconversion = 16.667;


                    output.rpm_l = new NVHFunctions.XY_Data
                    {
                        xdata = new double[indices_red.Length],
                        ydata = new double[indices_red.Length]

                    };

                    for (int ii = 0; ii < indices_red.Length; ii++)
                    {

                        output.rpm_l.xdata[ii] = indices_red[ii] / (double)output.samplerate;
                        output.rpm_l.ydata[ii] = output.samplerate * 60 / distances_red[ii] / output.speedfactor_l / unitconversion;
                        if (output.rpm_l.ydata[ii] > output.maxspeed_l) output.maxspeed_l = output.rpm_l.ydata[ii];
                        if (output.rpm_l.ydata[ii] < output.minspeed_l) output.minspeed_l = output.rpm_l.ydata[ii];
                    }
                }
                else
                {





                    output.rpm_l.xdata = new double[0];
                    output.rpm_l.ydata = new double[0];
                }
                
            }

            if (output.speedpulses_r > 0)
            {
                // int L = Array.FindAll(pulses_L, x => x == true).Length;
                int L = 0;
                for (int ii = 0; ii < pulses_R.Length; ii++) if (pulses_R[ii]) L++;

                var indices = new int[(int)Math.Ceiling(L / (double)output.speedpulses_r)];
                int jj = 0;
                int kk = 0;
                for (int ii = 0; ii < pulses_R.Length; ii++)
                {
                    if (pulses_R[ii])
                    {
                        if (kk == output.speedpulses_r) kk = 0;
                        if (kk == 0)
                        {
                            indices[jj] = ii;
                            jj++;
                        }
                        kk++;
                    }
                }
                if (indices.Length > 0)
                {

                    var distances = new int[indices.Length - 1];
                    for (int ii = 0; ii < distances.Length; ii++) distances[ii] = indices[ii + 1] - indices[ii];
                    int distanceisone = Array.FindAll(distances, x => x == 1).Length;
                    var indices_red = new int[indices.Length - 1 - distanceisone];
                    var distances_red = new int[indices_red.Length];
                    jj = 0;
                    for (int ii = 1; ii < indices.Length; ii++)
                    {
                        if (distances[ii - 1] != 1)
                        {
                            indices_red[jj] = indices[ii];
                            distances_red[jj] = distances[ii - 1];

                            jj++;
                        }
                    }


                    double unitconversion = 1;
                    if (output.unit_r.Length > 2) if (String.Equals("m/s", output.unit_r.Substring(0, 3))) unitconversion = 16.667;


                    output.rpm_r = new NVHFunctions.XY_Data
                    {
                        xdata = new double[indices_red.Length],
                        ydata = new double[indices_red.Length]

                    };

                    for (int ii = 0; ii < indices_red.Length; ii++)

                    {

                        output.rpm_r.xdata[ii] = indices_red[ii] / (double)output.samplerate;
                        output.rpm_r.ydata[ii] = output.samplerate * 60 / distances_red[ii] / output.speedfactor_l / unitconversion;
                        if (output.rpm_r.ydata[ii] > output.maxspeed_r) output.maxspeed_r = output.rpm_r.ydata[ii];
                        if (output.rpm_r.ydata[ii] < output.minspeed_r) output.minspeed_r = output.rpm_r.ydata[ii];
                    }
                }
                else
                {
                    output.rpm_r.xdata = new double[0];
                    output.rpm_r.ydata = new double[0];
                }

            }

            return output;
        }

        
        public int WriteAVLWaveFile(double[] audio_l, double fs, string filename, double[] audio_r = null, NVHFunctions.XY_Data rpm_l = default(NVHFunctions.XY_Data), NVHFunctions.XY_Data rpm_r = default(NVHFunctions.XY_Data), double calib_l = 0, double calib_r = 0, double dbref = 0.00002, bool fadeinout = false, bool[][] pulses = null, int pulses_l = 0, int pulses_r = 0, double factor_l = 0, double factor_r = 0)
        // this is currently the used version of this method!
        {
            if (audio_r == null) audio_r = new double[0];
            if (audio_l.Length == 0 && audio_r.Length == 0) return 1;
            if ((audio_l.Length != 0 && audio_r.Length != 0) && audio_l.Length != audio_r.Length) return 1;
            int samplerate = (int)Math.Round(fs);
            int nchannels = 0;
            if (audio_l.Length > 0) nchannels++;
            if (audio_r.Length > 0) nchannels++;


            // calculate calibration value

            double maxpascal_l = Math.Abs(audio_l[0]);
            double maxpascal_r;
            if (audio_r.Length > 0) maxpascal_r = Math.Abs(audio_r[0]);
            else maxpascal_r = 0;
            for (int ii = 0; ii < audio_l.Length; ii++)
            {
                if (Math.Abs(audio_l[ii]) > maxpascal_l) maxpascal_l = Math.Abs(audio_l[ii]);
                if (audio_r.Length > 0) if (Math.Abs(audio_r[ii]) > maxpascal_r) maxpascal_r = Math.Abs(audio_r[ii]);
            }
            if (calib_l == -1)
            {
                if (maxpascal_l == 0) calib_l = 2958795.24;
                else calib_l = 32000 / maxpascal_l * 256;
            }
            else calib_l = Math.Pow(2, 23) / (Math.Pow(10, calib_l / 20) * dbref * 2 * Math.Sqrt(2));
            if (calib_r == -1)
            {
                if (maxpascal_r == 0) calib_r = 2958795.24;
                else calib_r = 32000 / maxpascal_r * 256;
            }
            else calib_r = Math.Pow(2, 23) / (Math.Pow(10, calib_r / 20) * dbref * 2 * Math.Sqrt(2));

            var audio_l_int = new Int32[audio_l.Length];
            var audio_r_int = new Int32[audio_r.Length];


            for (int ii = 0; ii < audio_l.Length; ii++)
            {
                //audio_l_int[ii] = (Int32)(audio_l[ii] * calib_l * 256);
                //if (audio_r.Length > 0) audio_r_int[ii] = (Int32)(audio_r[ii] * calib_l * 256);
                audio_l[ii] *= calib_l * 256;
                if (audio_r.Length > 0) audio_r[ii] *= calib_r * 256;
            }

            if (fadeinout)
            {
                int N_fade = (int)(0.18 * samplerate);
                for (int ii = 0; ii < N_fade - 1; ii++)
                {
                    audio_l[ii] *= (double)ii / N_fade;
                    audio_l[audio_l.Length - 1 - ii] *= (double)ii / N_fade;
                    if (audio_r.Length > 0)
                    {
                        audio_r[ii] *= (double)ii / N_fade;
                        audio_r[audio_r.Length - 1 - ii] *= (double)ii / N_fade;
                    }
                }
            }

            for (int ii = 0; ii < audio_l.Length; ii++)
            {
                audio_l_int[ii] = (int)audio_l[ii];
                if (audio_r.Length > 0) audio_r_int[ii] = (int)audio_r[ii];
            }

            int speedpulses_l = -1;
            if (rpm_l.ydata != null) speedpulses_l = 60;
            double speedfactor_l = 1;

            int speedpulses_r = -1;
            if (rpm_r.ydata != null) speedpulses_r = 90;
            double speedfactor_r = 16.667;

            if (pulses_l != 0) speedpulses_l = pulses_l;
            if (pulses_r != 0) speedpulses_r = pulses_r;
            if (factor_l != 0) speedfactor_l = factor_l;
            if (factor_r != 0) speedfactor_r = factor_r;


            //BinaryWriter br;

            //var FcnClass = new NVHFunctions();
            using (BinaryWriter br = new BinaryWriter(new FileStream(filename, FileMode.Create)))
            {
                var riffchunk = new char[4] { 'R', 'I', 'F', 'F' };
                br.Write(riffchunk);
                uint riffchunksize = (UInt32)(3 * audio_l.Length + 3 * audio_r.Length + 36);
                br.Write(riffchunksize);
                var wavechunk = new char[4] { 'W', 'A', 'V', 'E' };
                br.Write(wavechunk);
                var fmtchunk = new char[4] { 'f', 'm', 't', ' ' };
                br.Write(fmtchunk);
                br.Write((UInt32)16);
                br.Write((byte)1);
                br.Write((byte)0);
                br.Write((byte)nchannels);
                br.Write((byte)0);
                br.Write((UInt32)samplerate);
                br.Write((UInt32)(nchannels * samplerate * 3));
                br.Write((byte)(nchannels * 3));
                br.Write((byte)0);
                br.Write((byte)24);
                br.Write((byte)0);

                var cominfo = new char[4] { 'C', 'O', 'M', 'I' };

                br.Write(cominfo);
                br.Write((UInt32)128);
                br.Write(calib_l);
                br.Write(calib_r);
                br.Write((byte)1);
                br.Write((UInt32)speedpulses_l);
                br.Write(speedfactor_l);
                var quantity_l = new char[3] { 'R', 'P', 'M' };
                br.Write(quantity_l);
                var description_l = new char[32] { 'E', 'n', 'g', 'i', 'n', 'e', ' ', 'S', 'p', 'e', 'e', 'd', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0' };
                br.Write(description_l);
                var unit_l = new char[8] { 'r', 'p', 'm', '\0', '\0', '\0', '\0', '\0' };
                br.Write(unit_l);

                br.Write((byte)2);
                br.Write((UInt32)speedpulses_r);
                br.Write(speedfactor_r);
                var quantity_r = new char[3] { 'V', 'e', 'l' };
                br.Write(quantity_r);
                var description_r = new char[32] { 'V', 'e', 'h', 'i', 'c', 'l', 'e', ' ', 'V', 'e', 'l', 'o', 'c', 'i', 't', 'y', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0' };
                br.Write(description_r);
                var unit_r = new char[8] { 'm', '/', 's', '\0', '\0', '\0', '\0', '\0' };
                br.Write(unit_r);

                var datachunk = new char[4] { 'd', 'a', 't', 'a' };
                br.Write(datachunk);
                br.Write((UInt32)(audio_l.Length * 3 + audio_r.Length * 3));

                // Prepare Pulse Signal
                var FcnClass = new NVHFunctions();
                var pulse_indices_l = new int[0];
                var pulse_indices_r = new int[0];
                if (rpm_l.ydata != null) pulse_indices_l = FcnClass.RPMVel2PulseIndices(rpm_l, audio_l.Length, samplerate,false);


                if (rpm_r.ydata != null) pulse_indices_r = FcnClass.RPMVel2PulseIndices(rpm_r, audio_r.Length, samplerate, true);
                if (pulses != null)
                {
                    pulse_indices_l = FcnClass.Pulses2PulseIndices(pulses[0]);
                    if (pulses.Length > 1)
                    {
                        if (pulses[1] != null) pulse_indices_r = FcnClass.Pulses2PulseIndices(pulses[1]);
                    }
                }

                if (pulse_indices_r.Length > 0 && audio_r_int.Length == 0) audio_r_int = new int[audio_l_int.Length];


                byte[] byteArray;
                int kk1 = 0;
                int kk2 = 0;

                for (int ii = 0; ii < audio_l.Length; ii++)     // Audio data is stored as 32 bit integer little endian. This means, byteArray[0] is omitted for 24 bit audio output. Therefore, the LSB is stored at byteArray[1].
                {
                    byteArray = BitConverter.GetBytes(audio_l_int[ii]);
                    if (pulse_indices_l.Length > 0)
                    {
                        if (pulse_indices_l[kk1] == ii)
                        {
                            byteArray[1] |= 1;
                            if (kk1 < pulse_indices_l.Length - 2) kk1++;
                        }
                        else byteArray[1] = (byte)(byteArray[1] & ~1);
                    }
                    for (int jj = 1; jj < 4; jj++) br.Write(byteArray[jj]);
                    if (audio_r_int.Length > 0)
                    {
                        byteArray = BitConverter.GetBytes(audio_r_int[ii]);
                        //kk = 0;
                        if (pulse_indices_r.Length > 0)
                        {
                            if (pulse_indices_r[kk2] == ii)
                            {
                                byteArray[1] |= 1;
                                if (kk2 < pulse_indices_r.Length - 2) kk2++;
                            }
                            else byteArray[1] = (byte)(byteArray[1] & ~1);
                        }
                        for (int jj = 1; jj < 4; jj++) br.Write(byteArray[jj]);
                    }
                    //else if (kk < pulse_indices_l.Length - 2) kk++;

                }
                br.Close();
            }
            //br = new BinaryWriter(new FileStream(filename, FileMode.Create));

            // write Riff chunk


            return 0;
        }


        public int WriteAVLWaveFile(int[] audio_l, double fs, double calib_l, double calib_r, string filename, int[] audio_r = null, NVHFunctions.XY_Data rpm_l = default(NVHFunctions.XY_Data), NVHFunctions.XY_Data rpm_r = default(NVHFunctions.XY_Data), double dbref = 0.00002, bool fadeinout = false, bool[][] pulses = null, int pulses_l = 0, int pulses_r = 0, double factor_l = 0, double factor_r = 0)
        {
            if (audio_r == null) audio_r = new int[0];
            if (audio_l.Length == 0 && audio_r.Length == 0) return 1;
            if ((audio_l.Length != 0 && audio_r.Length != 0) && audio_l.Length != audio_r.Length) return 1;
            int samplerate = (int)Math.Round(fs);
            int nchannels = 0;
            if (audio_l.Length > 0) nchannels++;
            if (audio_r.Length > 0) nchannels++;


            // calculate calibration value

            /*double maxpascal_l = Math.Abs(audio_l[0]);
            double maxpascal_r;
            if (audio_r.Length > 0) maxpascal_r = Math.Abs(audio_r[0]);
            else maxpascal_r = 0;
            for (int ii = 0; ii < audio_l.Length; ii++)
            {
                if (Math.Abs(audio_l[ii]) > maxpascal_l) maxpascal_l = Math.Abs(audio_l[ii]);
                if (audio_r.Length > 0) if (Math.Abs(audio_r[ii]) > maxpascal_r) maxpascal_r = Math.Abs(audio_r[ii]);
            }
            if (calib_l == -1)
            {
                if (maxpascal_l == 0) calib_l = 2958795.24;
                else calib_l = 32000 / maxpascal_l * 256;
            }
            else calib_l = Math.Pow(2, 23) / (Math.Pow(10, calib_l / 20) * dbref * 2 * Math.Sqrt(2));
            if (calib_r == -1)
            {
                if (maxpascal_r == 0) calib_r = 2958795.24;
                else calib_r = 32000 / maxpascal_r * 256;
            }
            else calib_r = Math.Pow(2, 23) / (Math.Pow(10, calib_r / 20) * dbref * 2 * Math.Sqrt(2));

            var audio_l_int = new Int32[audio_l.Length];
            var audio_r_int = new Int32[audio_r.Length];*/
            var audio_l_int = audio_l;
            var audio_r_int = audio_r;


           /* for (int ii = 0; ii < audio_l.Length; ii++)
            {
                //audio_l_int[ii] = (Int32)(audio_l[ii] * calib_l * 256);
                //if (audio_r.Length > 0) audio_r_int[ii] = (Int32)(audio_r[ii] * calib_l * 256);
                audio_l[ii] *= calib_l * 256;
                if (audio_r.Length > 0) audio_r[ii] *= calib_r * 256;
            }

            if (fadeinout)
            {
                int N_fade = (int)(0.18 * samplerate);
                for (int ii = 0; ii < N_fade - 1; ii++)
                {
                    audio_l[ii] *= (double)ii / N_fade;
                    audio_l[audio_l.Length - 1 - ii] *= (double)ii / N_fade;
                    if (audio_r.Length > 0)
                    {
                        audio_r[ii] *= (double)ii / N_fade;
                        audio_r[audio_r.Length - 1 - ii] *= (double)ii / N_fade;
                    }
                }
            }

            for (int ii = 0; ii < audio_l.Length; ii++)
            {
                audio_l_int[ii] = (int)audio_l[ii];
                if (audio_r.Length > 0) audio_r_int[ii] = (int)audio_r[ii];
            }*/

            int speedpulses_l = -1;
            if (rpm_l.xdata != null) speedpulses_l = 1;
            
            double speedfactor_l = 1;
            int speedpulses_r = -1;
            if (rpm_l.xdata != null) speedpulses_l = 1;
            
            double speedfactor_r = 1;
            if (pulses_l != 0) speedpulses_l = pulses_l;
            if (pulses_r != 0) speedpulses_r = pulses_r;
            if (factor_l != 0) speedfactor_l = factor_l;
            if (factor_r != 0) speedfactor_r = factor_r;


            //BinaryWriter br;

            //var FcnClass = new NVHFunctions();
            using (BinaryWriter br = new BinaryWriter(new FileStream(filename, FileMode.Create)))
            {
                var riffchunk = new char[4] { 'R', 'I', 'F', 'F' };
                br.Write(riffchunk);
                uint riffchunksize = (UInt32)(3 * audio_l.Length + 3 * audio_r.Length + 36);
                br.Write(riffchunksize);
                var wavechunk = new char[4] { 'W', 'A', 'V', 'E' };
                br.Write(wavechunk);
                var fmtchunk = new char[4] { 'f', 'm', 't', ' ' };
                br.Write(fmtchunk);
                br.Write((UInt32)16);
                br.Write((byte)1);
                br.Write((byte)0);
                br.Write((byte)nchannels);
                br.Write((byte)0);
                br.Write((UInt32)samplerate);
                br.Write((UInt32)(nchannels * samplerate * 3));
                br.Write((byte)(nchannels * 3));
                br.Write((byte)0);
                br.Write((byte)24);
                br.Write((byte)0);

                var cominfo = new char[4] { 'C', 'O', 'M', 'I' };

                br.Write(cominfo);
                br.Write((UInt32)128);
                br.Write(calib_l);
                br.Write(calib_r);
                br.Write((byte)1);
                br.Write((UInt32)speedpulses_l);
                br.Write(speedfactor_l);
                var quantity_l = new char[3] { 'R', 'P', 'M' };
                br.Write(quantity_l);
                var description_l = new char[32] { 'E', 'n', 'g', 'i', 'n', 'e', ' ', 'S', 'p', 'e', 'e', 'd', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0' };
                br.Write(description_l);
                var unit_l = new char[8] { 'r', 'p', 'm', '\0', '\0', '\0', '\0', '\0' };
                br.Write(unit_l);

                br.Write((byte)2);
                br.Write((UInt32)speedpulses_r);
                br.Write(speedfactor_r);
                var quantity_r = new char[3] { 'V', 'e', 'l' };
                br.Write(quantity_r);
                var description_r = new char[32] { 'V', 'e', 'h', 'i', 'c', 'l', 'e', ' ', 'V', 'e', 'l', 'o', 'c', 'i', 't', 'y', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0' };
                br.Write(description_r);
                var unit_r = new char[8] { 'm', '/', 's', '\0', '\0', '\0', '\0', '\0' };
                br.Write(unit_r);

                var datachunk = new char[4] { 'd', 'a', 't', 'a' };
                br.Write(datachunk);
                br.Write((UInt32)(audio_l.Length * 3 + audio_r.Length * 3));

                // Prepare Pulse Signal
                var FcnClass = new NVHFunctions();
                var pulse_indices_l = new int[0];
                var pulse_indices_r = new int[0];
                if (rpm_l.xdata != null) pulse_indices_l = FcnClass.RPMVel2PulseIndices(rpm_l, audio_l.Length, samplerate, false);
                

                if (rpm_r.xdata != null) pulse_indices_r = FcnClass.RPMVel2PulseIndices(rpm_r, audio_r.Length, samplerate, true);
                if (pulses != null)
                {
                    pulse_indices_l = FcnClass.Pulses2PulseIndices(pulses[0]);
                    if (pulses.Length > 1)
                    {
                        if (pulses[1] != null) pulse_indices_r = FcnClass.Pulses2PulseIndices(pulses[1]);
                    }
                }

                byte[] byteArray;
                int kk = 0;
                for (int ii = 0; ii < audio_l.Length; ii++)
                {
                    byteArray = BitConverter.GetBytes(audio_l_int[ii]);
                    if (rpm_l.xdata != null)
                    {
                        if (pulse_indices_l[kk] == ii)
                        {
                            byteArray[1] |= 1;
                            if (kk < pulse_indices_l.Length - 2) kk++;
                        }
                        else byteArray[1] = (byte)(byteArray[1] & ~1);
                    }
                    for (int jj = 1; jj < 4; jj++) br.Write(byteArray[jj]);
                    if (audio_r.Length > 0)
                    {
                        byteArray = BitConverter.GetBytes(audio_r_int[ii]);
                        
                        for (int jj = 1; jj < 4; jj++) br.Write(byteArray[jj]);
                    }
                }
                br.Close();
            }
                //br = new BinaryWriter(new FileStream(filename, FileMode.Create));

            // write Riff chunk

            
            return 0;
        }

        internal string ReadByteString(BinaryReader br, int bytes)
        {
            var chararray = new char[bytes];
            for (int ii = 0; ii < bytes; ii++) chararray[ii] = br.ReadChar();
            var outstr = new string(chararray);
            return outstr;
        }

        internal int Read24BitSample(BinaryReader br)
        {
            var sample4bytes = new byte[4];
            for (int ii = 1; ii < 4; ii++) sample4bytes[ii] = br.ReadByte();
            int sample = BitConverter.ToInt32(sample4bytes, 0);
            return sample;
        }

        internal int Read32BitSample(BinaryReader br, int format)
        {
            var sample4bytes = new byte[4];
            int sample;
            for (int ii = 0; ii < 4; ii++) sample4bytes[ii] = br.ReadByte();
            if (format == 3)
            {
                float sample_f = BitConverter.ToSingle(sample4bytes, 0);
                sample = (int)(sample_f * Math.Pow(2,31));
            }
            else sample = BitConverter.ToInt32(sample4bytes, 0);
            return sample;
        }
    }
    
}


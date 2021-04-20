using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NVHLibrary;
using csmatio.io;
using System.IO;
using System.Security.Cryptography;
using System.Reflection;


namespace NVHLibrary
{
    public class CRUISE_Class
    {
        public struct CRUISEData
        {
            public double[] EMSpd;
            public double[] ICESpd;
            public double[] EMTq;
            public double[] ICETq;
            public double[] VehVel;
            public double[] Time;
            public string error;
        }

        public struct PRSData
        {
            public float[][] matrix;
            public float[] enginestart;
            public int fs_enginestart;
            public int volmain;
            public int volzero;
            public int volneg;
            public string error;
        }

        


        public string DecryptStringFromBytes_Aes(byte[] cipherText, byte[] Key, byte[] IV)
        {
            // Check arguments.
            if (cipherText == null || cipherText.Length <= 0)
                throw new ArgumentNullException("cipherText");
            if (Key == null || Key.Length <= 0)
                throw new ArgumentNullException("Key");
            if (IV == null || IV.Length <= 0)
                throw new ArgumentNullException("IV");

            // Declare the string used to hold
            // the decrypted text.
            string plaintext = null;

            // Create an Aes object
            // with the specified key and IV.
            using (Aes aesAlg = Aes.Create())
            {
                aesAlg.Key = Key;
                aesAlg.IV = IV;

                // Create a decryptor to perform the stream transform.
                ICryptoTransform decryptor = aesAlg.CreateDecryptor(aesAlg.Key, aesAlg.IV);

                // Create the streams used for decryption.
                using (MemoryStream msDecrypt = new MemoryStream(cipherText))
                {
                    using (CryptoStream csDecrypt = new CryptoStream(msDecrypt, decryptor, CryptoStreamMode.Read))
                    {
                        using (StreamReader srDecrypt = new StreamReader(csDecrypt))
                        {

                            // Read the decrypted bytes from the decrypting stream
                            // and place them in a string.
                            plaintext = srDecrypt.ReadToEnd();
                        }
                    }
                }

            }

            return plaintext;

        }
        public PRSData Readprsfile (string filename_prs)
        {
            var result = new PRSData
            {
                error = "ok"
            };

            try
            {
                result.volmain = 0;
                result.volzero = 0;
                result.volneg = 0;

                result.matrix = new float[1536][];
                

                using (BinaryReader reader = new BinaryReader(File.Open(filename_prs, FileMode.Open)))
                {
                    result.volmain = reader.ReadInt32();
                    result.volzero = reader.ReadInt32();
                    result.volneg = reader.ReadInt32();

                    if (result.volmain > 10 || result.volmain < 0 || result.volzero > 10 || result.volzero < 0 || result.volneg > 10 || result.volneg < 0) throw new System.ArgumentOutOfRangeException();

                    for (int ii = 0; ii < 1536; ii++)
                    {
                        result.matrix[ii] = new float[129];
                        result.matrix[ii][0] = ii + 1;
                        for (int jj = 1; jj < 129; jj++) result.matrix[ii][jj] = reader.ReadSingle();
                    }
                    while (reader.BaseStream.Position != reader.BaseStream.Length)
                    {
                        result.enginestart = new float[reader.ReadInt32()];
                        result.fs_enginestart = reader.ReadInt32();
                        for (int ii = 0; ii < result.enginestart.Length; ii++) result.enginestart[ii] = reader.ReadSingle();
                    }


                }
                
            }
            catch
            {
                result.error = "Corrupt file or no read access!";
            }


            return result;
        }

        public CRUISEData ReadCRUISEData(string filename_mat, int emtimeconst, int icetimeconst)
        {
            double fs = 20;
            var result = new CRUISEData
            {
                error = "ok"
            };

            try
            {

                var NVHclass = new NVHFunctions();
                bool success = false;
                var mfr = new MatFileReader(filename_mat);
               
                if (!(mfr.Content.Count == 1 && mfr.Content.TryGetValue(Path.GetFileName(filename_mat).Substring(0, Path.GetFileName(filename_mat).Length - 4), out csmatio.types.MLArray test)))
                {
                    result.error = "Matlab File needs to contain exactly one structure with the same name as the file!";
                    return result;
                }
                var mainstruct = mfr.Content[Path.GetFileName(filename_mat).Substring(0, Path.GetFileName(filename_mat).Length - 4)] as csmatio.types.MLStructure;

                success = FindKey(mainstruct, "Time");
                if (!success)
                {
                    result.error = "Could not find \"Time\" Signal in .mat-File!";
                    return result;
                }

                var time_signal = mainstruct["Time"] as csmatio.types.MLStructure;
                if (time_signal == null)
                {
                    result.error = "Time Signal has no member \"Data\"";
                    return result;
                }
                var time_data = time_signal["Data"] as csmatio.types.MLDouble;

                success = FindKey(mainstruct, "EM_Tq");
                if (!success)
                {
                    result.error = "Could not find \"EM_Tq\" Signal in .mat-File!";
                    return result;
                }
                var emachine_torque = mainstruct["EM_Tq"] as csmatio.types.MLStructure;
                if (emachine_torque == null)
                {
                    result.error = "EM_Tq Signal has no member \"Data\"";
                    return result;
                }
                var emachine_torque_data = emachine_torque["Data"] as csmatio.types.MLDouble;

                success = FindKey(mainstruct, "EM_Spd");
                if (!success)
                {
                    result.error = "Could not find \"EM_Spd\" Signal in .mat-File!";
                    return result;
                }
                var emachine_speed = mainstruct["EM_Spd"] as csmatio.types.MLStructure;
                if (emachine_speed == null)
                {
                    result.error = "EM_Spd Signal has no member \"Data\"";
                    return result;
                }
                var emachine_speed_data = emachine_speed["Data"] as csmatio.types.MLDouble;

                success = FindKey(mainstruct, "ICE_Tq");
                if (!success)
                {
                    result.error = "Could not find \"ICE_Tq\" Signal in .mat-File!";
                    return result;
                }
                var ice_torque = mainstruct["ICE_Tq"] as csmatio.types.MLStructure;
                if (ice_torque == null)
                {
                    result.error = "ICE_Tq Signal has no member \"Data\"";
                    return result;
                }
                var ice_torque_data = ice_torque["Data"] as csmatio.types.MLDouble;

                success = FindKey(mainstruct, "ICE_Spd");
                if (!success)
                {
                    result.error = "Could not find \"ICE_Spd\" Signal in .mat-File!";
                    return result;
                }
                var ice_speed = mainstruct["ICE_Spd"] as csmatio.types.MLStructure;
                if (ice_speed == null)
                {
                    result.error = "ICE_Spd Signal has no member \"Data\"";
                    return result;
                }
                var ice_speed_data = ice_speed["Data"] as csmatio.types.MLDouble;

                success = FindKey(mainstruct, "Veh_Vel");
                if (!success)
                {
                    result.error = "Could not find \"Veh_Vel\" Signal in .mat-File!";
                    return result;
                }
                var veh_vel = mainstruct["Veh_Vel"] as csmatio.types.MLStructure;
                if (veh_vel == null)
                {
                    result.error = "Veh_Vel Signal has no member \"Data\"";
                    return result;
                }
                var veh_vel_data = veh_vel["Data"] as csmatio.types.MLDouble;
                

                var time_tmp = new double[time_data.Size];
                for (int ii = 0; ii < time_tmp.Length; ii++) time_tmp[ii] = time_data.Get(ii);
                result.Time = NVHclass.DeltaArray(0, 1 / fs, time_tmp[time_tmp.Length - 1]);

                var etorque_tmp = new double[emachine_torque_data.Size];
                for (int ii = 0; ii < etorque_tmp.Length; ii++) etorque_tmp[ii] = emachine_torque_data.Get(ii);
                var etorque = NVHclass.Interp1_linear(time_tmp, etorque_tmp, result.Time);
                etorque = NVHclass.RemoveNaN(etorque, etorque);

                var espeed_tmp = new double[emachine_speed_data.Size];
                for (int ii = 0; ii < espeed_tmp.Length; ii++) espeed_tmp[ii] = emachine_speed_data.Get(ii);
                var espeed = NVHclass.Interp1_linear(time_tmp, espeed_tmp, result.Time);
                espeed = NVHclass.RemoveNaN(espeed, espeed);

                var icetorque_tmp = new double[ice_torque_data.Size];
                for (int ii = 0; ii < icetorque_tmp.Length; ii++) icetorque_tmp[ii] = ice_torque_data.Get(ii);
                var icetorque = NVHclass.Interp1_linear(time_tmp, icetorque_tmp, result.Time);
                icetorque = NVHclass.RemoveNaN(icetorque, icetorque);

                var icespeed_tmp = new double[ice_speed_data.Size];
                for (int ii = 0; ii < icespeed_tmp.Length; ii++) icespeed_tmp[ii] = ice_speed_data.Get(ii);
                var icespeed = NVHclass.Interp1_linear(time_tmp, icespeed_tmp, result.Time);
                icespeed = NVHclass.RemoveNaN(icespeed, icespeed);

                var velocity_tmp = new double[veh_vel_data.Size];
                for (int ii = 0; ii < velocity_tmp.Length; ii++) velocity_tmp[ii] = veh_vel_data.Get(ii);
                var velocity = NVHclass.Interp1_linear(time_tmp, velocity_tmp, result.Time);
                result.VehVel = NVHclass.RemoveNaN(velocity, velocity);


                // Moving Average Filter
                //int ma = 33;
                result.EMTq = NVHclass.Moving(etorque, emtimeconst);
                result.EMSpd = NVHclass.Moving(espeed, emtimeconst);
                result.ICETq = NVHclass.Moving(icetorque, icetimeconst);
                result.ICESpd = NVHclass.Moving(icespeed, icetimeconst);
            }
            catch
            {
                result.error = "Corrupt file!";
            }
            return result;
        }

        private bool FindKey(csmatio.types.MLStructure mainstruct, string key)
        {
            bool success = false;
            for (int ii = 0; ii < mainstruct.Keys.Count; ii++)
            {
                if (mainstruct.Keys[ii].ToString().Equals(key))
                {
                    success = true;
                    break;
                }
            }
            return success;
        }
    }
}

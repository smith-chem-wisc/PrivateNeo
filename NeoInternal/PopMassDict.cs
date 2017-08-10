using System;
using System.Collections.Generic;
using System.Linq;

namespace NeoInternal
{
    class PopMassDict
    {
        public static int maxMissingConsecutivePeaks = 3;
        public static int decimalDigitsForFragmentMassRounding = 3;
        public static Dictionary<double, List<string>> massDict = new Dictionary<double, List<string>>();
        public static double[] keys;
        public static char[] AANames = new char[20] { 'G', 'A', 'S', 'P', 'V', 'T', 'L', 'I', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'C', 'Y', 'W' }; //20 common AA, ordered by mass assuming carbamido

        public void PopulateMassDictionary()
        {
            List<double> tempKeys = new List<double>();
            double maxMass = (maxMissingConsecutivePeaks + 1) * 186.079313 + (Constants.WATER_MONOISOTOPIC_MASS + 3); //If one sequence knocks it out of range, the next might not (VVVVW vs VVVVV where the first is out of range but the second is not) W-G prevents this     //peaks+1 for converting missed peaks into number of ambiguous aa
            int mer = (4) * (maxMissingConsecutivePeaks + 1); //maximum length allowed (4>W/G>3)  //peaks+1 for converting missed peaks into number of ambiguous aa
            //int mer = (2*(maxMissingConsecutivePeaks+1)); //only two aa should fit into one... work around to reduce computational demand
            int[] indexes = new int[mer];
            for (int i = 0; i < mer; i++)
            {
                indexes[i] = 0;
            }
            string seq = "";
            int length = 1;
            /*         List<char> firstAA = new List<char>();
                     foreach(char aa in AANames)
                     {
                         firstAA.Add(aa);
                     }
                     Parallel.ForEach(Partitioner.Creat(0,))*/
            while (indexes[0] < AANames.Count())
            {
                //get new seq
                seq = "";
                if (indexes[2] == 15)
                { }
                for (int n = 0; n < length; n++)
                {
                    seq += AANames[indexes[n]];
                }

                //if new seq is within range
                double fragMass = MassCalculator.MonoIsoptopicMass(seq);
                if (fragMass < maxMass)
                {
                    var rounded = Math.Round(fragMass, decimalDigitsForFragmentMassRounding);
                    List<string> value;
                    if (massDict.TryGetValue(rounded, out value))
                    {
                        if (!value.Contains(seq))
                            value.Add(seq);
                    }
                    else
                    {
                        massDict.Add(rounded, new List<string> { seq });
                        tempKeys.Add(rounded);
                    }


                    if (mer != length && fragMass + 57.0214 < maxMass) //if not last position
                    {
                        length++; //allow m to increase
                    }
                    else //don't increment length, we're happy right now!
                    {
                        if (indexes[length - 1] < AANames.Count() - 1) //if not last aa in possible aa
                        {
                            indexes[length - 1]++;
                        }
                        else //if it is, we need to go back a bit
                        {
                            indexes[length - 1] = 0;
                            indexes[length - 2]++; //could cause crashing with weird aa
                            length--;
                        }
                    }
                }
                else //if not in range, move back one
                {
                    indexes[length - 1] = 0;
                    indexes[length - 2]++; //could cause crashing with weird aa
                    length--;
                }

                //important catch to make sure going back doesn't result in an out of index exception
                for (int i = indexes.Count() - 1; i >= 0; i--)
                {
                    if (indexes[i] == AANames.Count())
                    {
                        if (i > 0)
                        {
                            length--;
                            indexes[i] = 0;
                            indexes[i - 1]++;
                        }
                    }
                }
            }
            tempKeys.Sort();
            keys = new double[tempKeys.Count()];
            for (int k = 0; k < tempKeys.Count(); k++)
            {
                keys[k] = tempKeys[k];
            }

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\Dictionary" + "3" + ".txt"))
            {
                foreach (double key in keys)
                {
                    List<string> value = new List<string>();
                    massDict.TryGetValue(key, out value);
                    bool satisfyRule = false;
                    foreach (string s in value)
                    {
                        if (s.Length <= maxMissingConsecutivePeaks + 1) { satisfyRule = true; }
                    }
                    if (satisfyRule)
                    {
                        string output = key.ToString() + '\t';
                        foreach (string s in value)
                        {
                            output += s + ';';
                        }
                        output = output.Substring(0, output.Length - 1);
                        file.WriteLine(output);
                    }
                }
            }
            //test
            /*      foreach(double key in keys)
                  {
                      List<string> value;
                      if (massDict.TryGetValue(key, out value))
                      {
                          //MessageBox.Show(key.ToString() + " with " + value.Count());
                          foreach(string v in value)
                          {
                             // MessageBox.Show(v);
                          }
                      }
                      else
                      {
                          //MessageBox.Show("No sequences found for " + key);
                      }
                  }*/
        }
    }
}

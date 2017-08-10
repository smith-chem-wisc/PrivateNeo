using System;
using System.Collections.Generic;
using System.Linq;

namespace NeoInternal
{
    public class InitialID : PSM
    {
        public string id { get; set; }
        private double[] peakHits;
        private Boolean target;
        private double qValue;
        public string seq { get; set; }
        public string score { get; set; }
        private List<PTM> ptms;

        public InitialID(string file, int scan, double expMass, string id, string seq, string peaks, string score) : base(file, scan, expMass)
        {
            this.ptms = new List<PTM>();
            this.seq = cleanSeq(seq);
            this.id = id;
            this.score = score;
            try
            {
                setPeakHits(peaks);
            }
            catch
            {
                double[] temp = new double[] { 0, 0 };
                this.peakHits = temp;
            }
        /*    catch
            {
                MessageBox.Show("woops");
                MessageBox.Show(file+" "+scan.ToString()+" "+seq+" "+peaks);
            }*/
        }

        public void setPeakHits(double[] peakHits)
        {
            this.peakHits = peakHits;
        }

        public void setPeakHits(string peakHits)
        {
            peakHits = peakHits.Replace("];[", ","); //used for EThcD to convert partition into addition
            peakHits = peakHits.Replace('"', '[');
            peakHits = peakHits.Replace("[", "");
            peakHits = peakHits.Replace("]", "");
            peakHits = peakHits.Replace(" ", "");
            peakHits = peakHits.Replace(";", "");
            string[] strArray = peakHits.Split(',').ToArray(); //split as strings
            List<double> dbList = new List<double>(); //double list to sort multi-ion input
            foreach(string str in strArray)
            {
                dbList.Add(Convert.ToDouble(str));
            }
            dbList.Sort();
            double[] dbArray = new double[strArray.Count()]; //new double array
            for (int i = 0; i < strArray.Count(); i++)
            {
                dbArray[i] = dbList[i];
            }
            this.peakHits = dbArray;
        }
        public void setTarget(Boolean target)
        {
            this.target = target;
        }
        public void setQValue(double qValue)
        {
            this.qValue = qValue;
        }


        public string cleanSeq(string seq)
        {
            bool ModificationOn = false;
            string ModificationName = "";
            int aaIndex = 0;
            string cleanedSeq = "";
            foreach (char amino_acid in seq)
            {
                if (amino_acid == ')') //only occurs at end of mod
                {
                    ModificationOn = false;
                    double modMass = MassCalculator.getPTMMass(ModificationName);
                    PTM ptm = new PTM(ModificationName, aaIndex, modMass);
                    this.ptms.Add(ptm);                                      
                }
                if (ModificationOn == true) //only occurs if "(" already found
                {
                    ModificationName += amino_acid;
                }
                if (amino_acid == '(') //start collecting PTM name
                {
                    ModificationOn = true;
                }
                if (ModificationOn == false && amino_acid != ')')
                {
                    cleanedSeq += amino_acid;
                }
            }
            return cleanedSeq;
        }
        public string getID()
        {
            return this.id;
        }
        public double[] getPeakHits()
        {
            return this.peakHits;
        }
        public Boolean getTarget()
        {
            return this.target;
        }
        public double getQValue()
        {
            return this.qValue;
        }

        public List<PTM> getPTMs()
        {
            return this.ptms;
        }
    }
}

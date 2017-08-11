using System.Collections.Generic;

namespace NeoInternal
{
    public class PSM
    {
        private string file;
        private int scan;
        private double expMass;
        private InitialID nInfo;
        private InitialID cInfo;
        private List<FusionCandidate> candidates;
        public List<Variant> variants { get; set; }
        public FusionCandidate.FusionType fusionType { get; set; }


        public PSM(string file, int scan, double expMass)
        {
            this.file = file;
            this.scan = scan;
            this.expMass = expMass;
            this.fusionType = FusionCandidate.FusionType.TS; //default
            this.candidates = new List<FusionCandidate>();
            this.variants = new List<Variant>();
        }

        public PSM(string file, int scan, double expMass, InitialID nInfo, InitialID cInfo)
        {
            this.file = file;
            this.scan = scan;
            this.expMass = expMass;
            this.nInfo = nInfo;
            this.cInfo = cInfo;
            this.candidates = new List<FusionCandidate>();
            this.variants = new List<Variant>();
        }

        public void setFile(string file)
        {
            this.file = file;
        }
        public void setScan(int scan)
        {
            this.scan = scan;
        }
        public void setExpMass(double expMass)
        {
            this.expMass = expMass;
        }
        public void setBInfo(InitialID info)
        {
            this.nInfo = info;
        }
        public void setYInfo(InitialID info)
        {
            this.cInfo = info;
        }
        public void addFusionCandidate(FusionCandidate candidate)
        {
            this.candidates.Add(candidate);
        }
        

        public string getFile()
        {
            return this.file;
        }
        public int getScan()
        {
            return this.scan;
        }
        public double getExpMass()
        {
            return this.expMass;
        }
        public InitialID getNInfo()
        {
            return this.nInfo;
        }
        public InitialID getCInfo()
        {
            return this.cInfo;
        }
        public List<FusionCandidate> getFusionCandidates()
        {
            return this.candidates;
        }
        //public void setPotentialSequences_B(string bSeq, int[] peaksHit)
        //{

        //}
    }
}

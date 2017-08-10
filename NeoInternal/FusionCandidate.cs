using NeoInternal;
using System;
using System.Collections.Generic;
using System.Linq;

namespace NeoInternal
{
    public class FusionCandidate
    {
        public string seq { get; set; }
        private List<int> junctionIndexes;
        private bool[] foundIons;
        public List<ParentInfo> parentInfo { get; set; } //Will have a ParentInfo object for every possible fragment* (*every fragment with length of 6 or more OR just the first length that was found (<6))
        public List<TranslatedParent> translatedParents { get; set; }
        public List<CisParent> cisParents { get; set; }
        public List<TransParent> transParents { get; set; }

        public enum FusionType { TL, NC, RC, TS } //ordered by priority //translated, normalCis, reverseCis, trans
        public FusionType fusionType { get; set; }
        //private List<FusionCandidate> fragSources;    

        public FusionCandidate(String seq)
        {
            this.seq = seq;
            this.junctionIndexes = new List<int>();
            //this.foundIons = new bool[seq.Length];
            fusionType = FusionType.TS; //default
            this.parentInfo = new List<ParentInfo>();
            this.translatedParents = new List<TranslatedParent>();
            this.cisParents = new List<CisParent>();
            this.transParents = new List<TransParent>();
        }

        public void addJunctionIndex(int index)
        {
            this.junctionIndexes.Add(index);
        }

        public void setFoundIons(bool[] foundIons)
        {
            this.foundIons = foundIons;
        }

        public void makeFoundIons()
        {
            this.foundIons= new bool[this.seq.Length]; //A|B|C|D|E|F|K| where the whole peptide peak is always placed arbitrarly at the c term
            for (int i = 0; i < this.foundIons.Count(); i++)
            {
                this.foundIons[i] = false;
            }
        }

        public void deepCopyFoundIons(FusionCandidate original)
        {
            this.foundIons = new bool[original.getFoundIons().Count()];
            for (int index = 0; index < this.foundIons.Count(); index++)
            {
                this.foundIons[index] = original.getFoundIons()[index];
            }
        }

        public List<int> getJunctionIndexes()
        {
            return this.junctionIndexes;
        }

        public bool[] getFoundIons()
        {
            return this.foundIons;
        }
    }
}

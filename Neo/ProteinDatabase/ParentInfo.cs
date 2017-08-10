using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Neo
{
    //used for parents of a given sequence fragment for FusionCandidate objects
    class ParentInfo
    {
        public List<TheoreticalProtein> theoreticalProteins { get; set; }
        public string fragFound { get; set; }
        public enum terminal { N, C };
        public terminal parentType { get; set; } //What terminus is the fragment from
        public ParentInfo(List<TheoreticalProtein> proteins, terminal parentType, string seqFound)
        {
            this.theoreticalProteins = proteins;
            this.parentType = parentType;
            this.fragFound = seqFound;
        }
    }
}

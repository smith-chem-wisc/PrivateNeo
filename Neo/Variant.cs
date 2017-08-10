using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Neo
{
    class Variant:TranslatedParent
    {
        public enum variantType { UM, SNP, PTM }; //unmodified, single nucleotide polymorphism, post translational modification
        public variantType varType { get; set; }
        public string pepSeq { get; set; }
        public Variant(string pepSeq, variantType type, string id, string seq, int start, int length) :base(id, seq, start, length)
        {
            this.pepSeq = pepSeq;
            this.varType = type;
        }
    }
}

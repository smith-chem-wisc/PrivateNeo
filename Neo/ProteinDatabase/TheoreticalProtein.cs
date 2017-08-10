using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Neo
{
    class TheoreticalProtein
    {
        private string id;
        private string seq;

        public TheoreticalProtein(string id, string seq)
        {
            this.id = id;
            this.seq = seq;
        }
        public void setID(string id)
        {
            this.id = id;
        }
        public void setSeq(string seq)
        {
            this.seq = seq;
        }

        public string getID()
        {
            return this.id;
        }
        public string getSeq()
        {
            return this.seq;
        }
    }
}

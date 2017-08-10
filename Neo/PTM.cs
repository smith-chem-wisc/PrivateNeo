using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Neo
{
    class PTM
    {
        private string name;
        private int index;
        private double mass;

        public PTM(string name, int index, double mass)
        {
            this.name = name;
            this.index = index;
            this.mass = mass;
        }
        public string getName()
        {
            return this.name;
        }
        public int getIndex()
        {
            return this.index;
        }
        public double getMass()
        {
            return this.mass;
        }
    }
}

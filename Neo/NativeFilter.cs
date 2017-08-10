using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.ComponentModel;
namespace Neo
{
    class NativeFilter
    {
        private Neo neo = null;
        private BackgroundWorker worker = null;

        public NativeFilter(Neo neo, BackgroundWorker worker)
        {
            this.neo = neo;
            this.worker = worker;
        }


    }
}

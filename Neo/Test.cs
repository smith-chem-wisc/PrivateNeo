using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace Neo
{
    class Test
    {
        public static void testSoftware(List<TheoreticalProtein> database)
        {
            string[] pepRead = File.ReadAllLines(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\Neo\tests\testPepListInput.txt");
            int nr = pepRead.Count();
            List<PSM> realPSMs = new List<PSM>();
            foreach(string s in pepRead)
            {
                string[] sArray = s.Split('\t');
                PSM tempPSM = new PSM(sArray[0], Convert.ToInt16(sArray[1]), Convert.ToDouble(sArray[2]));

            }
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\Neo\Results\testPepListOutput.txt"))
            {
              /*  AlternativeSequences altSeq = new AlternativeSequences(new Neo(), new System.ComponentModel.BackgroundWorker());
                altSeq.readProteins(database);
                foreach (string str in pepRead)
                {
                    if (altSeq.foundParent(str, ParentInfo.terminal.C,new FusionCandidate(str)))
                    {
                        file.WriteLine(str + '\t' + "1");
                    }
                    else
                    {
                        file.WriteLine(str + '\t' + "0");
                    }
                }*/

            }
        }
/*
        public static void readProteins(List<TheoreticalProtein> list)
        {
            localTheoreticals = list.Select(t => t.getSeq()).ToList();
        }

        public static bool foundParent(List<string> localTheoreticals, string frag)
        {
            List<string> matches = localTheoreticals.AsParallel().Where(x => x.Contains(frag)).ToList();
            if (matches.Count() > 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }*/
    }
}

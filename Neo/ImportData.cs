using NeoInternal;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Windows.Forms;

namespace Neo
{
    class ImportData
    {
        private BackgroundWorker worker = null;

        public ImportData(BackgroundWorker worker)
        {
            this.worker = worker;
        }

        public List<PSM> ImportPSMs(string nFileName, string cFileName)
        {
            List<PSM> MMOutput = new List<PSM>();
            //     try
            //   {
            string[] nInput = File.ReadAllLines(nFileName);
            string[] cInput = File.ReadAllLines(cFileName);

            string[] header = nInput[0].Split('\t').ToArray(); //assume both files have identical headers

            int fileNameIndex = -1;
            int scanNumberIndex = -1;
            int scanPrecursorMassIndex = -1;
            int proteinAccessionIndex = -1;
            int fullSequenceIndex = -1;
            int matchedIonsIndex = -1;
            int matchedIonCountsIndex = -1;
            int scoreIndex = -1;
            for (int col = 0; col < header.Length; col++)
            {
                if (header[col].Equals("File Name"))
                {
                    fileNameIndex = col;
                }
                else if (header[col].Equals("Scan Number"))
                {
                    scanNumberIndex = col;
                }
                else if (header[col].Equals("Precursor Mass"))
                {
                    scanPrecursorMassIndex = col;
                }
                else if (header[col].Equals("Protein Accession"))
                {
                    proteinAccessionIndex = col;
                }
                else if (header[col].Equals("Base Sequence")) //"FullSequence" should be used for the detection of FPs containing PTMs and for missed cleave/nonspecific peptides containing PTMs
                {
                    fullSequenceIndex = col;
                }
                else if (header[col].Equals("Matched Ion Masses"))
                {
                    matchedIonsIndex = col;
                }
                else if (header[col].Equals("Matched Ion Counts"))
                {
                    matchedIonCountsIndex = col;
                }
                else if (header[col].Equals("Score"))
                {
                    scoreIndex = col;
                }
            }
            List<InitialID> nAssignment = new List<InitialID>();
            List<InitialID> cAssignment = new List<InitialID>();

            for (int i = 1; i < nInput.Count(); i++)
            {
                string[] line = nInput[i].Split('\t').ToArray();
             //   if (Convert.ToDouble(line[matchedIonCountsIndex]) >= 4)
                {
                    InitialID id = new InitialID(line[fileNameIndex], Convert.ToInt32(line[scanNumberIndex]), Convert.ToDouble(line[scanPrecursorMassIndex]), line[proteinAccessionIndex], line[fullSequenceIndex], line[matchedIonsIndex], line[scoreIndex]);
                    nAssignment.Add(id);
                }
            }

            for (int i = 1; i < cInput.Count(); i++)
            {
                string[] line = cInput[i].Split('\t').ToArray();
             //   if (Convert.ToDouble(line[matchedIonCountsIndex]) >= 4)
                {
                    InitialID id = new InitialID(line[fileNameIndex], Convert.ToInt32(line[scanNumberIndex]), Convert.ToDouble(line[scanPrecursorMassIndex]), line[proteinAccessionIndex], line[fullSequenceIndex], line[matchedIonsIndex], line[scoreIndex]);
                    cAssignment.Add(id);
                }
            }
            //sort by scan number
            List<InitialID> nAssignmentSorted = nAssignment.OrderBy(o => o.getScan()).ToList();
            List<InitialID> cAssignmentSorted = cAssignment.OrderBy(o => o.getScan()).ToList();

            //remove scans not found in both files
            for (int i = 0; i < nAssignmentSorted.Count(); i++)
            {
                //    neo.SetProgressMax(bAssignmentSorted.Count());
                //  neo.IncrementProgress();
                this.worker.ReportProgress(Convert.ToInt16((Convert.ToDouble(i) / Convert.ToDouble(nAssignmentSorted.Count())) * 100));
                if (i < cAssignmentSorted.Count())
                {
                    if (nAssignmentSorted[i].getScan().Equals(cAssignmentSorted[i].getScan()))
                    {
                        PSM psm = new PSM(nAssignment[i].getFile(), nAssignmentSorted[i].getScan(), nAssignmentSorted[i].getExpMass(), nAssignmentSorted[i], cAssignmentSorted[i]);
                        MMOutput.Add(psm);
                        continue;
                    }
                    else if (nAssignmentSorted[i].getScan() < cAssignmentSorted[i].getScan()) //no information was found for the b scan using y ions, so remove it
                    {
                        nAssignmentSorted.Remove(nAssignmentSorted[i]);
                        i--;
                    }
                    else  //no information was found for the y scan using b ions, so remove it
                    {
                        cAssignmentSorted.Remove(cAssignmentSorted[i]);
                        i--;
                    }
                }
            }
            return MMOutput;
            /*     }
                 catch (System.IO.IOException e)
                 {
                     MessageBox.Show("Looks like one of your files is open. Please close it and try again.");
                     System.Environment.Exit(0);
                     return new List<PSM>();
                 }*/
        }

        public static string cleanSeq(string fullSeq)
        {
            string temp = fullSeq.Replace("[Common Fixed:Carbamidomethyl of C]", "");
            temp = temp.Replace("ptmlist:", "");
            temp = MassCalculator.RemoveNestedParentheses(temp, false);
            temp = temp.Replace("'", "");
            temp = temp.Replace("\"", "");
            return temp;
        }

        public List<TheoreticalProtein> ImportDatabase(string databaseFileName)
        {
            List<TheoreticalProtein> database = new List<TheoreticalProtein>();
            string[] FASTARead = File.ReadAllLines(databaseFileName);
            int nr = FASTARead.GetLength(0);
            string ProteinName = "";
            string FullSequence = "";
            for (int r = 0; r < nr; r++)
            {
                //  neo.SetProgressMax(nr);
                //neo.IncrementProgress();
                this.worker.ReportProgress(Convert.ToInt16((Convert.ToDouble(r) / Convert.ToDouble(nr)) * 100));
                string lr = FASTARead[r];
                if (lr.StartsWith(">"))
                {
                    if (FullSequence != "") //prevents first false add.
                    {
                        database.Add(new TheoreticalProtein(ProteinName, FullSequence)); //load FASTA entry into datatable
                    }
                    ProteinName = lr.Substring(1); //remove >
                    FullSequence = "";
                }
                else
                {
                    FullSequence += lr;
                }
            }
            for (int i = 0; i < database.Count; i++)
            {
                try //try to remove other info and just produce an accession number
                {
                    database[i].id = (database[i].id.Split('|').ToArray()[1]);
                }
                catch { MessageBox.Show("FASTA FAIL " + database[i].id); } //usually caused by missing >sp|
            }
            return database;
        }
    }
}

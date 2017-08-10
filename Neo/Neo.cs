using NeoInternal;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Windows.Forms;

namespace Neo
{
    public partial class Neo : Form
    {

        System.Windows.Forms.OpenFileDialog ofdNIons = new OpenFileDialog();
        System.Windows.Forms.OpenFileDialog ofdYCons = new OpenFileDialog();
        System.Windows.Forms.OpenFileDialog ofdFASTA = new OpenFileDialog();
        public static string nFileName;//= @"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\170523_Comp_Cutoff\2017-05-19-16-27-56_Comp-BY\Task1Search\04-29-13_B6_Frac5_4uL-Calibrated_allPSMs_OpenSearch.psmtsv";//= @"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\Neo\bOutput.txt";
        public static string cFileName;//= @"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\170523_Comp_Cutoff\2017-05-19-16-27-56_Comp-BY\Task2Search\04-29-13_B6_Frac5_4uL-Calibrated_allPSMs_OpenSearch.psmtsv";//= @"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\Neo\yOutput.txt";
        public static string databaseFileName;//= @"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\170523_Comp_Cutoff\Mixed\ClassicSearchOutputFASTA.txt";//= @"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\Neo\Frac5_Scan_ClassicSearchOutputFASTA.fasta";

        // public static List<PSM> MMOutput;

        public Neo()
        {
            InitializeComponent();
        }

        private void Neo_Load(object sender, EventArgs e)
        {
            MS1Tolerance.Text = "5";
            MS2Tolerance.Text = "0.01";
            MassCalculator.importMasses();
            // To report progress from the background worker we need to set this property
            backgroundWorker1.WorkerReportsProgress = true;
            // This event will be raised on the worker thread when the worker starts
            backgroundWorker1.DoWork += new DoWorkEventHandler(Engine);
            // This event will be raised when we call ReportProgress
            backgroundWorker1.ProgressChanged += new ProgressChangedEventHandler(backgroundWorker1_ProgressChanged);
        }

        private void Engine(object sender, EventArgs e)
        {
            this.Invoke(new MethodInvoker(delegate { StatustxtBox.Text = "Reading N and C files... (Task 0/6)"; }));
            ImportData import = new ImportData(backgroundWorker1);
            List<PSM> MMOutput=new List<PSM>();
            List<TheoreticalProtein> database = new List<TheoreticalProtein>();
            bool passedFileIO = true;
         //   try
            {
                MMOutput = import.ImportPSMs(nFileName, cFileName, out string error_message);
                MessageBox.Show(error_message);

                this.Invoke(new MethodInvoker(delegate { StatustxtBox.Text = "Reading database... (Task 1/6)"; }));
                ClearProgress();
                database = import.ImportDatabase(databaseFileName, out string error_message2);
                MessageBox.Show(error_message2);

                {
                    //Cool bit of code to identify number of random sequences of a specific length existing within a database. 5mer takes ~3 hours.
                    /*    int mer = 5;
                        int[] indexes = new int[mer];
                        for (int i = 0; i < mer; i++)
                        {
                            indexes[i] = 0;
                        }
                        int[] hits = new int[mer];
                        string seq = "";
                        int length = 1;

                        while (indexes[0] < AANames.Count())
                        {
                            //get new seq
                            seq = "";

                            for (int n = 0; n < length; n++)
                            {
                                seq += AANames[indexes[n]];
                            }

                            //if new seq is present in db
                            if (database.AsParallel().Where(x => x.seq.Contains(seq)).Any())
                            {
                                hits[length-1]++;
                                if (mer != length) //if not last position
                                {
                                    length++; //allow m to increase
                                }
                                else //don't increment length, we're happy right now!
                                {
                                    if (indexes[length-1] < AANames.Count() - 1) //if not last aa in possible aa
                                    {
                                        indexes[length - 1]++;
                                    }
                                    else //if it is, we need to go back a bit
                                    {
                                        indexes[length - 1] = 0;
                                        indexes[length - 2]++; //could cause crashing with weird aa
                                        length--;
                                    }
                                }
                            }
                            else //if not hit, move back one
                            {
                                indexes[length - 1]++;
                            }

                            for (int i = indexes.Count()-1; i >= 0; i--)
                            {
                                if (indexes[i] == AANames.Count())
                                {
                                    if (i > 0)
                                    {
                                        length--;
                                        indexes[i] = 0;
                                        indexes[i - 1]++;
                                    }
                                }
                            }
                        }
                        foreach (int i in hits)
                        {
                            MessageBox.Show(i.ToString());
                        }*/
                }


        }
        //    catch
          //  {
            //    MessageBox.Show("We encountered some issues while reading in your files. Common problems include: Column header names changed from MM v 0.0.108, an input file is currently open, or the absolute file path is too long. Please check your input and try again.");
              //  passedFileIO = false;
                //this.Invoke(new MethodInvoker(delegate { StatustxtBox.Text = ":("; }));              
            //}
            if (passedFileIO)
            {
                this.Invoke(new MethodInvoker(delegate { StatustxtBox.Text = "Splicing peptides... (Task 2/6)"; }));
                ClearProgress();
                SpliceFragments sf = new SpliceFragments(backgroundWorker1);
                List<PSM> candidates = sf.ExperimentalTheoreticalMatching(MMOutput, out string error_message);
                MessageBox.Show(error_message);
                { }
      //          MessageBox.Show(candidates.Count().ToString() + " putative fusion peptide PSMs were selected");

               // this.Invoke(new MethodInvoker(delegate { StatustxtBox.Text = "Reading Fragment Index... (Task 3/6)"; }));
                ClearProgress();
                AlternativeSequences altSeq = new AlternativeSequences(backgroundWorker1);

                this.Invoke(new MethodInvoker(delegate { StatustxtBox.Text = "Finding ambiquity... (Task 3/6)"; }));
                ClearProgress();
                altSeq.FindAmbiguity(candidates, database, out string error_message2);
                MessageBox.Show(error_message2);
                FalsePositives.generateDecoys = false;
                this.Invoke(new MethodInvoker(delegate { StatustxtBox.Text = "Recording variants... (Task 4/6)"; }));
                ClearProgress();
                FalsePositives rfp = new FalsePositives(backgroundWorker1);
                rfp.FindCommonFalsePositives(candidates, database, out string error_message3);
                MessageBox.Show(error_message3);

                //        MessageBox.Show("After filtering, "+candidates.Count().ToString() + " putative fusion peptides were generated");

                this.Invoke(new MethodInvoker(delegate { StatustxtBox.Text = "Exporting results... (Task 5/6)"; }));
                ClearProgress();
                string e3 = ExportData.ExportAll(candidates, databaseFileName);
                if (e3.Length > 0) MessageBox.Show(e3);

                this.Invoke(new MethodInvoker(delegate { StatustxtBox.Text = "Complete! (6/6)"; }));

                MessageBox.Show("Complete!");
            }
        }

        #region Private GUI

        private void bButton_Click(object sender, EventArgs e)
        {
            ofdNIons.Filter = "All Files (*.*)|*.*";
            ofdNIons.FilterIndex = 1;
            ofdNIons.Multiselect = false;

            if (ofdNIons.ShowDialog() == DialogResult.OK)
            {
                btxtBox.Text = ofdNIons.FileName;
            }
        }

        private void yButton_Click(object sender, EventArgs e)
        {
            ofdYCons.Filter = "All Files (*.*)|*.*";
            ofdYCons.FilterIndex = 1;
            ofdYCons.Multiselect = false;

            if (ofdYCons.ShowDialog() == DialogResult.OK)
            {
                ytxtBox.Text = ofdYCons.FileName;
            }
        }

        private void FASTAButton_Click(object sender, EventArgs e)
        {
            ofdFASTA.Filter = "All Files (*.*)|*.*";
            ofdFASTA.FilterIndex = 1;
            ofdFASTA.Multiselect = false;

            if (ofdFASTA.ShowDialog() == DialogResult.OK)
            {
                FASTAtxtBox.Text = ofdFASTA.FileName;
            }
        }

        private void StartButton_Click(object sender, EventArgs e)
        {
            if(nFileName!=null && cFileName != null && databaseFileName != null)
            {
                if (!backgroundWorker1.IsBusy)
                {
                    backgroundWorker1.RunWorkerAsync();
                }
                else
                {
                    MessageBox.Show("It's already processing!");
                }                
            }
            else
            {
                MessageBox.Show("We seem to be missing a file... Please try again!");
            }
        }

        private void btxtBox_TextChanged(object sender, EventArgs e)
        {
            nFileName = btxtBox.Text;
        }

        private void ytxtBox_TextChanged(object sender, EventArgs e)
        {
            cFileName = ytxtBox.Text;
        }

        private void FASTAtxtBox_TextChanged(object sender, EventArgs e)
        {
            databaseFileName = FASTAtxtBox.Text;
        }

        private void label4_Click(object sender, EventArgs e)
        {

        }

        private void StatustxtBox_TextChanged(object sender, EventArgs e)
        {

        }

        private void progressBar1_Click(object sender, EventArgs e)
        {

        }

        public void SetProgressMax(int max)
        {
            progressBar1.Minimum = 0;
            progressBar1.Maximum = max;
        }

        public void IncrementProgress()
        {
            this.Invoke(new MethodInvoker(delegate { progressBar1.Increment(1); }));
        }

        public void ClearProgress()
        {
            this.Invoke(new MethodInvoker(delegate { progressBar1.Value = 0; }));   
        }

        public void backgroundWorker1_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            this.Invoke(new MethodInvoker(delegate { progressBar1.Value = e.ProgressPercentage; }));
        }

        private void backgroundWorker1_DoWork(object sender, DoWorkEventArgs e)
        {

        }

        private void MS1Tolerance_TextChanged(object sender, EventArgs e)
        {
            try
            {
                FalsePositives.precursorMassTolerancePpm = Convert.ToDouble(MS1Tolerance.Text);
            }
            catch { MS1Tolerance.Text = FalsePositives.precursorMassTolerancePpm.ToString(); }
        }

        private void MS2Tolerance_TextChanged(object sender, EventArgs e)
        {
            try
            {
                AlternativeSequences.productMassToleranceDa = Convert.ToDouble(MS2Tolerance.Text);
            }
            catch { MS2Tolerance.Text = FalsePositives.precursorMassTolerancePpm.ToString(); }
        }

        #endregion Private GUI

        private void label1_Click(object sender, EventArgs e)
        {

        }

        private void bCheckBox_CheckedChanged(object sender, EventArgs e)
        {
            if(bCheckBox.Checked)
            {
                AlternativeSequences.ionsUsed.Add(IonType.b);
            }
            else
            {
                AlternativeSequences.ionsUsed.Remove(IonType.b);
            }
        }

        private void yCheckBox_CheckedChanged(object sender, EventArgs e)
        {
            if (yCheckBox.Checked)
            {
                AlternativeSequences.ionsUsed.Add(IonType.y);
            }
            else
            {
                AlternativeSequences.ionsUsed.Remove(IonType.y);
            }
        }

        private void cCheckBox_CheckedChanged(object sender, EventArgs e)
        {
            if (cCheckBox.Checked)
            {
                AlternativeSequences.ionsUsed.Add(IonType.c);
            }
            else
            {
                AlternativeSequences.ionsUsed.Remove(IonType.c);
            }
        }

        private void zdotCheckBox_CheckedChanged(object sender, EventArgs e)
        {
            if (zdotCheckBox.Checked)
            {
                AlternativeSequences.ionsUsed.Add(IonType.zdot);
            }
            else
            {
                AlternativeSequences.ionsUsed.Remove(IonType.zdot);
            }
        }

        private void checkBox1_CheckedChanged(object sender, EventArgs e)
        {
            if(checkBox1.Checked)
            {
                FalsePositives.generateDecoys = true;
            }
            else
            {
                FalsePositives.generateDecoys = false;
            }
        }
    }
}

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.ComponentModel;
using System.IO;


namespace Neo
{
    class AlternativeSequences
    {
        private Neo neo;
        private BackgroundWorker worker = null;
        public List<TheoreticalProtein> theoreticalProteins;
        //    public static HashSet<string> foundSequences = new HashSet<string>();
        public HashSet<string> notFoundSequences = new HashSet<string>();
        public Dictionary<string, List<TheoreticalProtein>> foundSequences = new Dictionary<string, List<TheoreticalProtein>>();
        private const int maxMissingConsecutivePeaks = 2;
        private const int maxNumPossibleSequences = 250;
        private const int decimalDigitsForFragmentMassRounding = 3;
        public Dictionary<double, string[]> massDict = new Dictionary<double, string[]>();
        public double[] keys;
        public List<List<string>> fragmentSeqs = new List<List<string>>();

        public AlternativeSequences(Neo neo, BackgroundWorker worker)
        {
            this.neo = neo;
            this.worker = worker;
        }

        public void FindAmiguity(List<PSM> candidates,List<TheoreticalProtein> theoreticalProteins)
        {
            ReadProteins(theoreticalProteins);
            ReadMassDictionary();
            for (int i = 0; i < candidates.Count(); i++) //must be mutable while iterating
            {
                PSM psm = candidates[i];
                psm.fusionType = FusionCandidate.FusionType.TS; //for some maddening reason, this is not arriving here as trans, but instead translated
                if (IsTooMessy(psm)) //having explosion of combinations when greater than 3 consequtive peaks producing tens of thousands of sequences ids, causes hanging
                {
                    candidates.Remove(psm);
                    i--;
                }
                else
                {
                    if (GeneratePossibleSequences(psm)) //return true if fewer than specified number of ambiguities
                    {
                        if (!PossibleCandidate(psm))
                        {
                            candidates.Remove(psm);
                            i--;
                        }
                    }
                    else
                    {
                        candidates.Remove(psm);
                        i--;
                    }
                    {
                        /*Replacer(psm);
                    if (Flipper(psm)) //this flips sequences and return true if too many ambiguous sequences are made
                    {
                        candidates.Remove(psm);
                        i--;
                    }
                    else
                    {
                        if (FindParents(psm))
                        {
                            candidates.Remove(psm);
                            i--;
                        }
                    }*/
                    }
                }
                this.worker.ReportProgress(Convert.ToInt16((Convert.ToDouble(i) / Convert.ToDouble(candidates.Count())) * 100));
            }
            // RemoveTranslatablePeptides(candidates);
        }

        public bool IsTooMessy(PSM psm) //return true if too messy for confident identification
        {
            List<string> baseSequences = new List<string>();
            int currentBestScore = 0;
            for (int index = 0; index < psm.getFusionCandidates().Count(); index++)
            {
                bool badID = false;
                FusionCandidate fc = psm.getFusionCandidates()[index];
                findIons(fc, psm);
                int consecutiveMissedCounter = 0;
                int totalHitCounter = 0;
                foreach (bool b in fc.getFoundIons())
                {
                    if (consecutiveMissedCounter > maxMissingConsecutivePeaks) //if too many permutations possible because of an unmapped region
                    {
                        badID = true;
                    }
                    else if (!b)
                    {
                        consecutiveMissedCounter++;
                    }
                    else
                    {
                        totalHitCounter++;
                        consecutiveMissedCounter = 0;
                    } //only care about consecutive
                }
                bool isRepeat = false;
                if (baseSequences.Contains(psm.getFusionCandidates()[index].getSeq()))
                { isRepeat = true; }

                if (totalHitCounter > currentBestScore && !badID)//the others were worse, so delete them
                {
                    for(int i=0; i<index; i=0)
                    {
                        psm.getFusionCandidates().Remove(psm.getFusionCandidates()[0]);
                        index--;
                    }
                    currentBestScore = totalHitCounter;
                    baseSequences = new List<string> { psm.getFusionCandidates()[index].getSeq() };
                }
                else if(totalHitCounter<currentBestScore|badID|isRepeat)
                {
                    psm.getFusionCandidates().Remove(psm.getFusionCandidates()[index]);
                    index--;
                }
            }
            //If there's anything left
            if (psm.getFusionCandidates().Count() > 0) //It wasn't too messy! Yay!
            {
                return false;
            }
            else //this might be a fusion peptide, but we won't get any valuable information from this spectra, so discard it
            {
                return true;
            }
        }

        public void ReadProteins(List<TheoreticalProtein> list)
        {
            theoreticalProteins = list;
        }

        //use ion hits to know where peaks have been found by morpheus and where there is ambiguity
        public static void findIons(FusionCandidate fusionCandidate, PSM psm)
        {
            double[] nPeaks = psm.getNInfo().getPeakHits(); //get peaks
            double[] cPeaks = psm.getCInfo().getPeakHits();
            fusionCandidate.makeFoundIons();
            string candSeq = fusionCandidate.getSeq();
            bool[] foundIons = fusionCandidate.getFoundIons();
            
            //find which aa have peaks
            for (int i = 0; i < foundIons.Count()-1; i++)
            {
                //B IONS//
                if (Neo.ionsUsed.Contains(IonType.b))
                {
                    double bTheoMass = MassCalculator.MonoIsoptopicMass(candSeq.Substring(0, 1 + i)) - Constants.WATER_MONOISOTOPIC_MASS; 
                    foreach (PTM ptm in psm.getNInfo().getPTMs())
                    {
                        if (ptm.getIndex() <= i)
                        {
                            bTheoMass += ptm.getMass();
                        }
                    }
                    foreach (double expPeak in nPeaks)
                    {
                        if (expPeak > bTheoMass - Neo.productMassToleranceDa && expPeak < bTheoMass + Neo.productMassToleranceDa)
                        {
                            foundIons[i] = true;
                        }
                    }
                }
                //Y IONS//
                if (Neo.ionsUsed.Contains(IonType.y))
                {
                    double yTheoMass = MassCalculator.MonoIsoptopicMass(candSeq.Substring(candSeq.Length - 1 - i, i + 1));
                    foreach (PTM ptm in psm.getCInfo().getPTMs())
                    {
                        if (ptm.getIndex() >= candSeq.Length - 2 - i)
                        {
                            yTheoMass += ptm.getMass();
                        }
                    }
                    foreach (double expPeak in cPeaks)
                    {
                        if (expPeak > yTheoMass - Neo.productMassToleranceDa && expPeak < yTheoMass + Neo.productMassToleranceDa)
                        {
                            foundIons[foundIons.Count() - 2 - i] = true;
                        }
                    }
                }
                //C IONS//
                if (Neo.ionsUsed.Contains(IonType.c))
                {
                    double cTheoMass = MassCalculator.MonoIsoptopicMass(candSeq.Substring(0, 1 + i)) - Constants.WATER_MONOISOTOPIC_MASS + Constants.nitrogenMonoisotopicMass + 3 * Constants.hydrogenMonoisotopicMass; 
                    foreach (PTM ptm in psm.getNInfo().getPTMs())
                    {
                        if (ptm.getIndex() <= i)
                        {
                            cTheoMass += ptm.getMass();
                        }
                    }
                    foreach (double expPeak in nPeaks)
                    {
                        if (expPeak > cTheoMass - Neo.productMassToleranceDa && expPeak < cTheoMass + Neo.productMassToleranceDa)
                        {
                            foundIons[i] = true;
                        }
                    }
                }
                //ZDOT IONS//
                if (Neo.ionsUsed.Contains(IonType.zdot))
                {
                    double zdotTheoMass = MassCalculator.MonoIsoptopicMass(candSeq.Substring(candSeq.Length - 1 - i, i + 1)) - Constants.nitrogenMonoisotopicMass - 2 * Constants.hydrogenMonoisotopicMass;
                    foreach (PTM ptm in psm.getCInfo().getPTMs())
                    {
                        if (ptm.getIndex() >= candSeq.Length - 2 - i)
                        {
                            zdotTheoMass += ptm.getMass();
                        }
                    }
                    foreach (double expPeak in cPeaks)
                    {
                        if (expPeak > zdotTheoMass - Neo.productMassToleranceDa && expPeak < zdotTheoMass + Neo.productMassToleranceDa)
                        {
                            foundIons[foundIons.Count() - 2 - i] = true;
                        }
                    }
                }
            }
            //foundIons[0] = true; //AspN always starts with a D
            foundIons[foundIons.Count() - 1] = true;//A|B|C|D|E|F|K| where the whole peptide peak is always placed arbitrarly at the c term
        }
        
        public void Replacer(PSM psm)
        {
            //operate under assumption that no peak provides a greater confidence for Q than AG, and a peak for AG prevents Q from matching.
            string[] oldAA = new string[] { "I", "L" };//, "Q", "AG", "GA" }; //Q is only transformed into AG, since GA will be formed during flippage
            string[] newAA = new string[] { "L", "I" };//, "AG", "Q", "Q" };
            //produce alternative sequences
            foreach (FusionCandidate fusionCandidate in psm.getFusionCandidates())
            {
                findIons(fusionCandidate, psm);
            }
            bool done = false;
            int globalIndex = 0;
            //Do substitutions specified above
            while (!done)
            {
                done = true; //assume we're done
                List<FusionCandidate> tempList = new List<FusionCandidate>();
                foreach (FusionCandidate fusionCandidate in psm.getFusionCandidates())
                {
                    if (fusionCandidate.getSeq().Length > globalIndex) //prevent crashing, use to tell when done
                    {
                        done = false; //at least one sequence is still as long as the globalIndex, so we're not done
                        for (int aa = 0; aa < oldAA.Count(); aa++)
                        {
                            ExchangeAA(psm, fusionCandidate, globalIndex, tempList, oldAA[aa], newAA[aa]);
                        }
                    }
                }
                List<FusionCandidate> full = psm.getFusionCandidates();
                List<string> foundSeq = new List<string>();
                foreach (FusionCandidate fusionCandidate in full)
                {
                    foundSeq.Add(fusionCandidate.getSeq());
                }
                foreach (FusionCandidate tempCandidate in tempList)
                {
                    if (!foundSeq.Contains(tempCandidate.getSeq()))
                    {
                        findIons(tempCandidate, psm);
                        psm.getFusionCandidates().Add(tempCandidate);;
                    }
                }
                globalIndex++;
            }
        }

        public void HeapsAlgorithm(char[] aa, int currentSize, int seqLength, List<string> combinations)
        {
            if(currentSize==1)
            {
                combinations.Add(new string(aa));
            }
            for (int i=0; i<currentSize; i++)
            {
                HeapsAlgorithm(aa, currentSize - 1, seqLength, combinations);

                if(currentSize%2==1)
                {
                    Swap(aa, 0, currentSize - 1);
                }
                else
                {
                    Swap(aa, i, currentSize - 1);
                }
            }
        }

        public void Swap(char[] aa, int i, int j)
        {
            char temp = aa[i];
            aa[i] = aa[j];
            aa[j] = temp;
        }

        public bool Flipper(PSM psm) //IonType ion, string seq, double[] peaks) //This method flips aa with no detected peaks between them
        {
            /*     MessageBox.Show("In Flipper");
                 string str1 = "";
                 foreach(bool b in psm.getFusionCandidates()[0].getFoundIons())
                 {
                     str1 += b;
                 }
                 MessageBox.Show(str1);*/
            bool done = false; //use to determine if we are outside the length index of all possible peptides
            int globalIndex = 0;
            //Do substitutions specified above
            //    using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\Neo\Results\WeirdFlippage.txt"))
            //  {
            while (!done)
            {
                //     file.WriteLine("---" + globalIndex.ToString() + "---" + globalIndex.ToString() + "---" + globalIndex.ToString() + "---" + globalIndex.ToString() + "---");
                done = true; //let's assume we're done and correct it later if we're not
                List<FusionCandidate> tempCandidates = new List<FusionCandidate>();
                // MessageBox.Show(psm.getFusionCandidates().Count().ToString() + " " + globalIndex + " " + psm.getFusionCandidates()[0].getFoundIons().Count());
                if (psm.getFusionCandidates().Count() > 1000) //if there are more than 1000 possible sequences, this is junk and we are not searching them all
                {
                    // MessageBox.Show("Hit for psm" + psm.getScan());
                    return true;
                }
                for (int j = 0; j < psm.getFusionCandidates().Count(); j++)
                {
                    FusionCandidate fusionCandidate = psm.getFusionCandidates()[j];
                    //      file.WriteLine(fusionCandidate.getSeq());
                    if (fusionCandidate.getFoundIons().Count() > globalIndex) //prevent crashing, use to tell when done
                    {
                        done = false; //We're not done, because at least one fusioncandidate sequence length is still greater than the global index
                        string fusionSeq = fusionCandidate.getSeq();
                        bool[] IonFound = fusionCandidate.getFoundIons();
                        if (IonFound[globalIndex])
                        {
                            int mostRecent = -1; //most recent Ion found prior to this one
                            for (int i = 0; i < globalIndex; i++)
                            {
                                if (IonFound[i])
                                {
                                    mostRecent = i; //save most recent hit, exclusive of the current index
                                }
                            }
                            if (mostRecent + 1 != globalIndex) //if there's a gap
                            {
                                string ambiguousSeq = fusionSeq.Substring(mostRecent + 1, globalIndex - mostRecent);
                                List<string> combinations = new List<string>();
                                HeapsAlgorithm(ambiguousSeq.ToCharArray(), globalIndex - mostRecent, globalIndex - mostRecent, combinations); //populates combinations
                                foreach (string str in combinations)
                                {
                                    string nTermSeq = fusionSeq.Substring(0, mostRecent + 1);
                                    string cTermSeq = fusionSeq.Substring(globalIndex + 1, fusionSeq.Length - globalIndex - 1);
                                    string novelSeq = nTermSeq + str + cTermSeq;
                                    FusionCandidate tempCandidate = new FusionCandidate(novelSeq);
                                    if (!novelSeq.Equals(fusionSeq))
                                    {
                                        //                file.WriteLine('\t' + "Added " + novelSeq);
                                        tempCandidate.deepCopyFoundIons(fusionCandidate);
                                        tempCandidates.Add(tempCandidate);
                                    }
                                    if (str.Contains("GA") || str.Contains("AG")) //if these were formed without a peak between them, it could be a Q instead
                                    { //sloppy encoding for things like A|B|C|ADG|E|F, assumes low messiness tolerance
                                        str.Replace("AG", "Q");
                                        str.Replace("GA", "Q");
                                        if (!novelSeq.Equals(fusionSeq)) //if not already found, add it
                                        {
                                            //                file.WriteLine('\t' + "Added " + novelSeq);
                                            tempCandidate.makeFoundIons(); //has a new found Ions because length is no longer the same
                                            findIons(tempCandidate, psm);
                                            tempCandidates.Add(tempCandidate);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                List<FusionCandidate> full = psm.getFusionCandidates();
                List<string> foundSeq = new List<string>(); //get list of all old FP sequences
                foreach (FusionCandidate fusionCandidate in full)
                {
                    foundSeq.Add(fusionCandidate.getSeq());
                }
                foreach (FusionCandidate fusionCandidate in tempCandidates)
                {
                    if (!foundSeq.Contains(fusionCandidate.getSeq())) //if new list contains a FP sequence not already found, add it.
                    {
                        foundSeq.Add(fusionCandidate.getSeq());
                        full.Add(fusionCandidate);
                    }
                }
                globalIndex++;
            }
            return false;
            // }
        }

        public void ExchangeAA( PSM psm, FusionCandidate fusionCandidate, int globalIndex, List<FusionCandidate> tempList, string oldAA, string newAA)
        {
            int fragLength = oldAA.Length;
            string fusionSeq = fusionCandidate.getSeq();
            if (globalIndex + fragLength <= fusionCandidate.getSeq().Length) //if won't crash
            {
                if (fusionSeq.Substring(globalIndex,fragLength).Equals(oldAA)) //if found
                {
                    bool existingPeaks = false;
                    for (int j = globalIndex; j < globalIndex + oldAA.Length-1; j++)
                    {
                        if (fusionCandidate.getFoundIons()[j])
                        {
                            existingPeaks = true;
                        }
                    }
                    if (!existingPeaks)//if no peaks dictating arrangment of amino acids
                    {
                        //substring is better
                        string novelSeq = fusionCandidate.getSeq();
                        string nTermSeq = novelSeq.Substring(0, globalIndex);
                        string cTermSeq = novelSeq.Substring(globalIndex + oldAA.Length, novelSeq.Length - (globalIndex + oldAA.Length));
                        novelSeq = nTermSeq + newAA + cTermSeq;
                        FusionCandidate tempCandidate = new FusionCandidate(novelSeq);
                        if (oldAA.Length != newAA.Length)
                        {
                            findIons(tempCandidate, psm);
                        }
                        else
                        {
                            tempCandidate.deepCopyFoundIons(fusionCandidate);
                        }
                        tempList.Add(tempCandidate);
                    }
                }
            }
        }     

        public void ReadMassDictionary()
        {
            List<double> tempKeys = new List<double>();

            using (StreamReader masses = new StreamReader(Path.Combine(Path.GetDirectoryName(Environment.GetCommandLineArgs()[0]), "Dictionary" + maxMissingConsecutivePeaks + ".txt"))) //file located in Morpheus folder
            {
                while (masses.Peek() != -1)
                {
                    string line = masses.ReadLine();
                    string[] fields = line.Split('\t');
                    double key = Convert.ToDouble(fields[0]);
                    string[] sequences = fields[1].Split(';');
                    massDict.Add(key, sequences);
                    tempKeys.Add(key);
                }
            }
            keys = new double[tempKeys.Count()];
            for(int i=0; i<tempKeys.Count(); i++)
            {
                keys[i] = tempKeys[i];
            }
        }

        public bool GeneratePossibleSequences(PSM psm) //returns false if over the specified number of sequences are generated
        {
            List<string> foundSeq = new List<string>(); //get list of all FP sequences
            foreach (FusionCandidate fusionCandidate in psm.getFusionCandidates())
            {
                findIons(fusionCandidate, psm); //populate the foundIons array
                foundSeq.Add(fusionCandidate.getSeq());
            }
            bool done = false;
            int globalIndex = 0;
            while (!done)
            {
                done = true; //let's assume we're done and correct it later if we're not
                if (psm.getFusionCandidates().Count() > maxNumPossibleSequences) //if there are more than a set number of possible sequences, this is junk and we are not searching them all
                {
                    return false;
                }
                for (int fc = 0; fc < psm.getFusionCandidates().Count(); fc++)
                {
                    FusionCandidate fusionCandidate = psm.getFusionCandidates()[fc];
                    if (fusionCandidate.getFoundIons().Count() > globalIndex) //prevent crashing, use to tell when done by hitting end of fc
                    {
                        List<FusionCandidate> tempCandidates = new List<FusionCandidate>(); //fill with possible sequences

                        done = false; //We're not done, because at least one fusion candidate sequence length is still greater than the global index
                        string fusionSeq = fusionCandidate.getSeq();
                        bool[] IonFound = fusionCandidate.getFoundIons();
                        if (IonFound[globalIndex]) //only look for ambiguity if a peak was found to provide the stop point.
                        {
                            int mostRecent = -1; //most recent Ion found prior to this one (start point)
                            for (int i = 0; i < globalIndex; i++) //identify start point
                            {
                                if (IonFound[i])
                                {
                                    mostRecent = i; //save most recent hit, exclusive of the current index
                                }
                            }


                            string ambiguousFrag = fusionSeq.Substring(mostRecent + 1, globalIndex - mostRecent);
                            double key = MassCalculator.MonoIsoptopicMass(ambiguousFrag);

                            List<string> combinations = new List<string>();

                            double closestPeak = double.NaN;
                            var ipos = Array.BinarySearch(keys, key);
                            if (ipos < 0)
                                ipos = ~ipos;

                            if (ipos > 0)
                            {
                                var downIpos = ipos - 1;
                                // Try down
                                while (downIpos >= 0)
                                {
                                    closestPeak = keys[downIpos];
                                    if (closestPeak > key - Neo.productMassToleranceDa && closestPeak < key + Neo.productMassToleranceDa)
                                    {
                                        string[] value;
                                        if (massDict.TryGetValue(closestPeak, out value))
                                        {
                                            foreach (string frag in value)
                                            {
                                                combinations.Add(frag);
                                            }
                                        }
                                    }
                                    else
                                        break;
                                    downIpos--;
                                }
                            }
                            if (ipos < keys.Length)
                            {
                                var upIpos = ipos;
                                // Try here and up
                                while (upIpos < keys.Length)
                                {
                                    closestPeak = keys[upIpos];
                                    if (closestPeak > key - Neo.productMassToleranceDa && closestPeak < key + Neo.productMassToleranceDa)
                                    {
                                        string[] value;
                                        if (massDict.TryGetValue(closestPeak, out value))
                                        {
                                            foreach (string frag in value)
                                            {
                                                combinations.Add(frag);
                                            }
                                        }
                                    }
                                    else
                                        break;
                                    upIpos++;
                                }
                            }

                            foreach (string str in combinations)
                            {
                                string nTermSeq = fusionSeq.Substring(0, mostRecent + 1);
                                string cTermSeq = fusionSeq.Substring(globalIndex + 1, fusionSeq.Length - globalIndex - 1);
                                string novelSeq = nTermSeq + str + cTermSeq;
                                FusionCandidate tempCandidate = new FusionCandidate(novelSeq);
                                tempCandidates.Add(tempCandidate);
                            }
                        }
                        foreach (FusionCandidate newfc in tempCandidates)
                        {
                            if (!foundSeq.Contains(newfc.getSeq())) //if new FP sequence, add it.
                            {
                                foundSeq.Add(newfc.getSeq());
                                findIons(newfc,psm);
                                psm.getFusionCandidates().Add(newfc);
                            }
                        }
                    }
                }
                globalIndex++;
                if(psm.getFusionCandidates().Count()>maxNumPossibleSequences)
                { 
                    return false;
                }
            }
            return true;
        }

        //returns true if a full fusion sequence could not be made or was found in the database, making it translated instead of a novel fusion.
        public bool PossibleCandidate(PSM psm) 
        {
            foundSequences = new Dictionary<string, List<TheoreticalProtein>>(); //used for finding longer fragments than those previously identified. Also populates ParentInfo
            notFoundSequences = new HashSet<string>(); //don't bother looking for these fragments, since we know they don't exist. Good for multiple homologous putative fusion peptide sequences

            //  Stopwatch sw = new Stopwatch();

            //conduct an initial search of each candidate's full sequence to identify any that are translated
            for (int i = 0; i < psm.getFusionCandidates().Count(); i++) //foreach fusion peptide sequence that could map to this scan
            {
                string novelSeq = psm.getFusionCandidates()[i].getSeq();
                if (foundParent(novelSeq, ParentInfo.terminal.C, psm.getFusionCandidates()[i], false)) //check really quick to see if the whole thing exists as is. If so, assign it as translated. Terminal C was arbitrarily chosen
                {
                    foreach (ParentInfo info in psm.getFusionCandidates()[i].parentInfo)
                    {
                        foreach (TheoreticalProtein protein in info.theoreticalProteins)
                        {
                            if (protein.getSeq().Contains(novelSeq)) //if translated
                            {
                                psm.getFusionCandidates()[i].translatedParents.Add(new TranslatedParent(protein.getID(), protein.getSeq(), protein.getSeq().IndexOf(psm.getFusionCandidates()[i].getSeq()), psm.getFusionCandidates()[i].getSeq().Length));
                            }
                        }
                    }
                    psm.getFusionCandidates()[i].fusionType = FusionCandidate.FusionType.TL;
                    psm.fusionType = psm.getFusionCandidates()[i].fusionType;
                    for (int j = 0; j < psm.getFusionCandidates().Count(); j++)
                    {
                        if (j != i)
                        {
                            psm.getFusionCandidates().Remove(psm.getFusionCandidates()[j]);
                            j--;
                            i--;
                        }
                    }
                    return true;
                }
            }
            for (int i=0; i<psm.getFusionCandidates().Count(); i++) //foreach fusion peptide sequence that could map to this scan
            {
                //sw.StartFindParents
                if (!isViable(psm.getFusionCandidates()[i])) //remove this fusion peptide sequence if the parent fragments cannot be found with the given database
                {
                    psm.getFusionCandidates().Remove(psm.getFusionCandidates()[i]);
                    i--;
                }
                else
                {
                    DetermineFusionCandidateType(psm.getFusionCandidates()[i]); //cis, trans?
                    if(psm.fusionType>psm.getFusionCandidates()[i].fusionType) //if more likely than previous types, change the psm type (golf scoring)
                    {
                        psm.fusionType = psm.getFusionCandidates()[i].fusionType;
                    }
                    if (psm.fusionType.Equals(FusionCandidate.FusionType.TL))//if there's a possible sequence that's present in the database, it is likely correct and is it is not worth it to identify parents of other sequences.
                    {
                        //remove all other candidates
                        for(int j=0; j<psm.getFusionCandidates().Count();j++)
                        {
                            if(j!=i)
                            {
                                psm.getFusionCandidates().Remove(psm.getFusionCandidates()[j]);
                                j--;
                                i--;
                            }
                        }
                        return true;
                    }
                }
         //       else if(psm.getFusionCandidates()[i].getJunctionIndexes().Count()> psm.getFusionCandidates()[i].getSeq().Length) //if the PSM could match to the database
           //     {
                    //end this loop, the psm is likely not a fusion and instead a translated peptide
             //       return true;
               // }
                /*  sw.Stop();
                  TimeSpan ts = sw.Elapsed;
                  string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
                      ts.Hours, ts.Minutes, ts.Seconds,
                      ts.Milliseconds / 10);
                  MessageBox.Show("RunTime " + elapsedTime+" progress: "+i.ToString()+"/"+ psm.getFusionCandidates().Count().ToString());
                  sw.Reset();*/
            }

            if(psm.getFusionCandidates().Count()==0) //if no candidates are left, we couldn't make the sequence with the database and we'll discard it.
            {
                return false;
            }

            return true;
        }

        public bool isViable(FusionCandidate tempCandidate) //returns if sequence could be made from one or two proteins in database and writes fusion type, parents, and junctions to fusionCandidate
        {
            //need to check that each index is viable
            string novelSeq = tempCandidate.getSeq();
     /*       if (foundParent(novelSeq, ParentInfo.terminal.C, tempCandidate, false)) //check really quick to see if the whole thing exists as is. If so, assign it as translated. Terminal C was arbitrarily chosen
            {
                tempCandidate.fusionType = FusionCandidate.FusionType.translated;
                return true;
            }
            else*/
            {
                //N//
                int nTermParentLength = novelSeq.Length-1;//length-1, because if the whole thing existed we wouldn't have made it to the else loop. Several edits are made to reflect this
                if (nTermParentLength > 6) //used to speed up search by finding an ideal starting point
                {
                    nTermParentLength = 6; //low point of random probability (5 or 7 may also be suitable)
                }

                bool done = false;
                bool foundFirstSearch = false;

                //First pass search
                string testFrag = novelSeq.Substring(0, nTermParentLength);
                if (foundParent(testFrag, ParentInfo.terminal.N, tempCandidate, foundFirstSearch)) //if found
                {
                    foundFirstSearch = true;
                    nTermParentLength++; 
                }
                else //if not found
                {
                    foundFirstSearch = false;
                    nTermParentLength--;
                }

                //All other passes
                while (nTermParentLength < novelSeq.Length && nTermParentLength > 0 && !done) //while in range and not done
                {
                    testFrag = novelSeq.Substring(0, nTermParentLength);
                    if (foundParent(testFrag, ParentInfo.terminal.N, tempCandidate, foundFirstSearch)) //if found
                    {
                        if (!foundFirstSearch) 
                        {
                            nTermParentLength--;
                            done = true;
                        }
                        nTermParentLength++;
                    }
                    else //if not found
                    {
                        if(foundFirstSearch)
                        {
                            done = true;
                        }
                        nTermParentLength--;
                    }
                }
                //if (!done) //if a match was not found for b (i.e. the whole sequence was found), we are done and it is not a viable fusion
                //{ return false; }
                /*       List<ParentInfo> nTermPI = new List<ParentInfo>();
                       foreach(ParentInfo PI in tempCandidate.parentInfo) //save these matches, as they will be overwritten by Y searches
                       {
                           nTermPI.Add(PI);
                       }*/

                //C//
                done = false; //reset tracker
                foundFirstSearch = false;
                int cTermParentLength = novelSeq.Length-1;
                if (cTermParentLength > 6) //used to speed up search by finding an ideal starting point
                {
                    cTermParentLength = 6; //low point of random probability
                }

                testFrag = novelSeq.Substring(novelSeq.Length - cTermParentLength, cTermParentLength);
                //First pass search
                if (foundParent(testFrag, ParentInfo.terminal.C, tempCandidate, foundFirstSearch)) //if found
                {
                    foundFirstSearch = true;
                    cTermParentLength++;
                }
                else //if not found
                {
                    foundFirstSearch = false;
                    cTermParentLength--;
                }

                while (cTermParentLength > 0 && cTermParentLength < novelSeq.Length && !done)
                {
                    testFrag = novelSeq.Substring(novelSeq.Length - cTermParentLength, cTermParentLength);
                    if (foundParent(testFrag, ParentInfo.terminal.C, tempCandidate, foundFirstSearch))
                    {
                        if (!foundFirstSearch)
                        {
                            cTermParentLength--;
                            done = true;
                        }
                        cTermParentLength++;
                    }
                    else
                    {
                        if(foundFirstSearch)
                        {
                            done=true;
                        }
                        cTermParentLength--;
                    }
                }
                {
                    /*       foreach (ParentInfo PI in nTermPI) //append B searches 
                               {
                                   tempCandidate.parentInfo.Add(PI);
                               }*/
                    //Add Result
                }
                if (cTermParentLength + nTermParentLength < novelSeq.Length) //if no overlap
                {
                    return false;
                }
                else
                {
                    for (int junction = tempCandidate.getSeq().Length-cTermParentLength-1; junction < nTermParentLength; junction++)
                    {
                        tempCandidate.addJunctionIndex(junction);
                    }
                    return true;
                }
            }
            {
                /* //need to check that each index is viable
                    string novelSeq = tempCandidate.getSeq();
                    //X//
                    int nTermParentLength = novelSeq.Length;
                    int nTermIndex = 0;
                    while (nTermParentLength > 1 && nTermIndex==0)
                    {
                        string testFrag = novelSeq.Substring(0, nTermParentLength);
                        if (foundParent(testFrag))
                        {
                            nTermIndex = nTermParentLength - 1;
                        }
                        nTermParentLength--;
                    }
                    //Y//
                    int cTermParentLength = novelSeq.Length;
                    int cTermIndex = novelSeq.Length-1; //aribtrary start at end where cterm parent contributes 1 aa
                    while (cTermParentLength > 0 && cTermIndex==novelSeq.Length-1)
                    {
                        string testFrag = novelSeq.Substring(novelSeq.Length - cTermParentLength, cTermParentLength);
                        if (foundParent(testFrag))
                        {
                            cTermIndex = novelSeq.Length-cTermParentLength-1;
                        }
                        cTermParentLength--;
                    }
                    //Add Result
                    if (cTermIndex > nTermIndex) //if no overlap
                    {
                        return false;
                    }
                    else
                    {
                        for (int junction = cTermIndex; junction <= nTermIndex; junction++)
                        {
                            tempCandidate.addJunctionIndex(junction);
                        }
                        return true;
                    }*/
            }
        }

        public bool foundParent(string frag, ParentInfo.terminal terminal, FusionCandidate candidate, bool foundFirstSearch)
        {
            //localTheoreticals.AsParallel().Where(x => x.Contains(frag)).ToList();
            if (notFoundSequences.Contains(frag)) //has the fragment been searched but not found before?
            {
                return false;
            }

            List<TheoreticalProtein> matches = new List<TheoreticalProtein>();
            if (foundSequences.TryGetValue(frag, out matches)) //has the fragment been searched AND found before?
            {
                candidate.parentInfo.Add(new ParentInfo(matches, terminal, frag));
                return true;
            }
     /*       if (foundSequences.Where(x => x.Contains(frag)).Any())
            {
                return true;
            }*/


            if (foundFirstSearch) //Has something smaller been found before? Well, then we can just search against those found sequences
            {
                string shorterFrag = "";
                if(terminal.Equals(ParentInfo.terminal.N))
                {
                    shorterFrag = frag.Substring(0, frag.Length - 1);
                }
                else //if C
                {
                    shorterFrag = frag.Substring(1, frag.Length - 1);
                }

                foreach (ParentInfo info in candidate.parentInfo)
                {
                    if (info.parentType.Equals(terminal) && info.fragFound.Equals(shorterFrag))
                    {
                        List<TheoreticalProtein> tempProtList = new List<TheoreticalProtein>();
                        info.theoreticalProteins.ForEach(protein => tempProtList.Add(protein));
                        matches = tempProtList.AsParallel().Where(x => x.getSeq().Contains(frag)).ToList();
                    }
                }
            }
            else //it hasn't been found before... we need to search against the whole database :(
            {
                matches = theoreticalProteins.AsParallel().Where(x => x.getSeq().Contains(frag)).ToList();
            }
            if (matches != null && matches.Count()>0)
            {
                foundSequences.Add(frag, matches);

                //if(matches.Count()<=10)
                //{
                //         candidate.parentInfo = new List<ParentInfo>(); //clear it
                //       foreach(TheoreticalProtein match in matches)
                //     {
                //       ParentInfo parentInfo = new ParentInfo(match.getID(), terminal, frag);
                //     candidate.parentInfo.Add(parentInfo);
                //}
                // }
                candidate.parentInfo.Add(new ParentInfo(matches, terminal, frag));
                return true;
            }
            else
            {
                notFoundSequences.Add(frag);
                return false;
            }
        }

        public void DetermineFusionCandidateType(FusionCandidate fusionCandidate)
        {
            if(!fusionCandidate.fusionType.Equals(FusionCandidate.FusionType.TL))
            {
                string sequence = fusionCandidate.getSeq();
                foreach(ParentInfo info in fusionCandidate.parentInfo)
                {
                    int foundLength = info.fragFound.Length;

                    string compFrag = "";// fusionCandidate.getSeq().Substring()
                    if(info.parentType.Equals(ParentInfo.terminal.N))
                    {
                        compFrag = fusionCandidate.getSeq().Substring(foundLength, fusionCandidate.getSeq().Length - foundLength);
                    }
                    else
                    {
                        compFrag = fusionCandidate.getSeq().Substring(0, fusionCandidate.getSeq().Length - foundLength);
                    }

                    foreach (TheoreticalProtein protein in info.theoreticalProteins)
                    {
                        //get the index(es) of where the found fragment is
                        string subProt = protein.getSeq();
                        List<int> originalIndexes = new List<int>();
                        int pastIndex = 0;
                        while (subProt.Contains(info.fragFound))
                        {
                            int newIndex = subProt.IndexOf(info.fragFound);
                            originalIndexes.Add(newIndex + pastIndex);
                            subProt = subProt.Substring(newIndex + 1, subProt.Length - newIndex - 1); //need to remove old match
                            pastIndex += newIndex + 1;
                        }
                        fusionCandidate.transParents.Add(new TransParent(protein.getID(), protein.getSeq(), originalIndexes, foundLength, info.parentType));

                        //get the index(es) of where the complimentary fragment is (if it's a cis fusion peptide)
                        subProt = protein.getSeq();
                        List<int> complementaryIndexes = new List<int>();
                        pastIndex = 0;
                        while (subProt.Contains(compFrag))
                        {
                            int newIndex = subProt.IndexOf(compFrag);
                            complementaryIndexes.Add(newIndex + pastIndex);
                            subProt = subProt.Substring(newIndex + 1, subProt.Length - newIndex - 1); //need to remove old match
                            pastIndex += newIndex + 1;
                        }
                        if(complementaryIndexes.Count()>0) //if it is not trans
                        {
                            //if it is cis
                            if (info.parentType.Equals(ParentInfo.terminal.N))
                            {
                                fusionCandidate.cisParents.Add(new CisParent(protein.getID(), protein.getSeq(), originalIndexes, foundLength, complementaryIndexes, compFrag.Length));
                            }
                            else
                            {
                                fusionCandidate.cisParents.Add(new CisParent(protein.getID(), protein.getSeq(), complementaryIndexes, compFrag.Length, originalIndexes, foundLength));
                            }
                        }
                    }
                }
            }
            else
            {
                string seq = fusionCandidate.getSeq();
                foreach (ParentInfo info in fusionCandidate.parentInfo)
                {
                    foreach (TheoreticalProtein protein in info.theoreticalProteins)
                    {
                        if (protein.getSeq().Contains(seq)) //if translated
                        {
                            fusionCandidate.translatedParents.Add(new TranslatedParent(protein.getID(), protein.getSeq(), protein.getSeq().IndexOf(fusionCandidate.getSeq()), fusionCandidate.getSeq().Length));
                        }
                    }
                }
            }
            foreach(CisParent cisParent in fusionCandidate.cisParents)
            {
                if(cisParent.cisType<fusionCandidate.fusionType)
                {
                    fusionCandidate.fusionType = cisParent.cisType;
                }
            }

        }

        public void RemoveTranslatablePeptides(List<PSM> psms)
        {
            for (int i = 0; i < psms.Count(); i++)
            {
                bool remove = false;
                foreach (FusionCandidate fc in psms[i].getFusionCandidates())
                {
                    if (fc.getJunctionIndexes().Count() > fc.getSeq().Length) //delete it
                    {
                        remove = true;
                    }
                }
                if (remove)
                {
                    psms.Remove(psms[i]);
                    i--;
                }
            }
        }
    }
}

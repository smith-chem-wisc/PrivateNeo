using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;

namespace NeoInternal
{
    public class AlternativeSequences
    {
        private BackgroundWorker worker = null;
        public List<TheoreticalProtein> theoreticalProteins;
        public HashSet<string> notFoundSequences = new HashSet<string>();
        public Dictionary<string, List<TheoreticalProtein>> foundSequences = new Dictionary<string, List<TheoreticalProtein>>();
        private const int maxMissingConsecutivePeaks = 2;
        private const int maxNumPossibleSequences = 250;
        private const int decimalDigitsForFragmentMassRounding = 3;
        public Dictionary<double, string[]> massDict = new Dictionary<double, string[]>();
        public double[] keys;
        public List<List<string>> fragmentSeqs = new List<List<string>>();
        public static double productMassToleranceDa; //(Da)
        public static List<IonType> ionsUsed = new List<IonType> { IonType.b, IonType.y };

        public AlternativeSequences( BackgroundWorker worker)
        {
            this.worker = worker;
        }

        public void FindAmbiguity(List<PSM> candidates, List<TheoreticalProtein> theoreticalProteins, out string error_message)
        {
            error_message = "";
            ReadProteins(theoreticalProteins);
            ReadMassDictionary();
            for (int i = 0; i < candidates.Count(); i++) //must be mutable while iterating
            {
                PSM psm = candidates[i];
                psm.fusionType = FusionCandidate.FusionType.TS; //for some maddening reason, this is not arriving here as trans, but instead translated
                if (IsTooMessy(psm, out string e)) //having explosion of combinations when greater than 3 consequtive peaks producing tens of thousands of sequences ids, causes hanging
                {
                    candidates.Remove(psm);
                    i--;
                    error_message += e;
                }
                else
                {
                    if (GeneratePossibleSequences(psm, out string error_message1) && !PossibleCandidate(psm)) //return true if fewer than specified number of ambiguities
                    {
                        candidates.Remove(psm);
                        i--;
                    }
                    error_message += error_message1;
                }
                this.worker.ReportProgress(Convert.ToInt16((Convert.ToDouble(i) / Convert.ToDouble(candidates.Count())) * 100));
            }
        }

        public bool IsTooMessy(PSM psm, out string error_message) //return true if too messy for confident identification
        {
            error_message = "";
            List<string> baseSequences = new List<string>();
            int currentBestScore = 0;
            for (int index = 0; index < psm.getFusionCandidates().Count(); index++)
            {
                bool badID = false;
                FusionCandidate fc = psm.getFusionCandidates()[index];
                findIons(fc, psm, out string error_message1);
                error_message += error_message1;
                int consecutiveMissedCounter = 0;
                int totalHitCounter = 0;
                foreach (bool b in fc.getFoundIons())
                {
                    if (consecutiveMissedCounter > maxMissingConsecutivePeaks) //if too many permutations possible because of an unmapped region
                        badID = true;
                    else if (!b)
                        consecutiveMissedCounter++;
                    else
                    {
                        totalHitCounter++;
                        consecutiveMissedCounter = 0;
                    } //only care about consecutive
                }
                bool isRepeat = false;
                if (baseSequences.Contains(psm.getFusionCandidates()[index].seq))
                    isRepeat = true;

                if (totalHitCounter > currentBestScore && !badID)//the others were worse, so delete them
                {
                    for(int i=0; i<index; i=0)
                    {
                        psm.getFusionCandidates().Remove(psm.getFusionCandidates()[0]);
                        index--;
                    }
                    currentBestScore = totalHitCounter;
                    baseSequences = new List<string> { psm.getFusionCandidates()[index].seq };
                }
                else if(totalHitCounter<currentBestScore|badID|isRepeat)
                {
                    psm.getFusionCandidates().Remove(psm.getFusionCandidates()[index]);
                    index--;
                }
            }
            //If there's anything left
            if (psm.getFusionCandidates().Count() > 0) //It wasn't too messy! Yay!
                return false;
            else //this might be a fusion peptide, but we won't get any valuable information from this spectra, so discard it
                return true;
        }

        public void ReadProteins(List<TheoreticalProtein> list)
        {
            theoreticalProteins = list;
        }

        //use ion hits to know where peaks have been found by morpheus and where there is ambiguity
        public static void findIons(FusionCandidate fusionCandidate, PSM psm, out string error_message)
        {
            error_message = "";
            double[] nPeaks = psm.getNInfo().getPeakHits(); //get peaks
            double[] cPeaks = psm.getCInfo().getPeakHits();
            fusionCandidate.makeFoundIons();
            string candSeq = fusionCandidate.seq;
            bool[] foundIons = fusionCandidate.getFoundIons();
            
            //find which aa have peaks
            for (int i = 0; i < foundIons.Count()-1; i++)
            {
                //B IONS//
                if (ionsUsed.Contains(IonType.b))
                {
                    double bTheoMass = MassCalculator.MonoIsoptopicMass(candSeq.Substring(0, 1 + i), out string error_message2) - Constants.WATER_MONOISOTOPIC_MASS;
                    error_message += error_message2;
                    foreach (PTM ptm in psm.getNInfo().getPTMs())
                    {
                        if (ptm.index <= i)
                            bTheoMass += ptm.mass;
                    }
                    foreach (double expPeak in nPeaks)
                    {
                        if (expPeak > bTheoMass - productMassToleranceDa && expPeak < bTheoMass + productMassToleranceDa)
                            foundIons[i] = true;
                    }
                }
                //Y IONS//
                if (ionsUsed.Contains(IonType.y))
                {
                    double yTheoMass = MassCalculator.MonoIsoptopicMass(candSeq.Substring(candSeq.Length - 1 - i, i + 1), out string error_message3);
                    error_message += error_message3;
                    foreach (PTM ptm in psm.getCInfo().getPTMs())
                    {
                        if (ptm.index >= candSeq.Length - 2 - i)
                            yTheoMass += ptm.mass;
                    }
                    foreach (double expPeak in cPeaks)
                    {
                        if (expPeak > yTheoMass - productMassToleranceDa && expPeak < yTheoMass + productMassToleranceDa)
                            foundIons[foundIons.Count() - 2 - i] = true;
                    }
                }
                //C IONS//
                if (ionsUsed.Contains(IonType.c))
                {
                    double cTheoMass = MassCalculator.MonoIsoptopicMass(candSeq.Substring(0, 1 + i), out string error_message4) - Constants.WATER_MONOISOTOPIC_MASS + Constants.nitrogenMonoisotopicMass + 3 * Constants.hydrogenMonoisotopicMass;
                    error_message += error_message4;
                    foreach (PTM ptm in psm.getNInfo().getPTMs())
                    {
                        if (ptm.index <= i)
                            cTheoMass += ptm.mass;
                    }
                    foreach (double expPeak in nPeaks)
                    {
                        if (expPeak > cTheoMass - productMassToleranceDa && expPeak < cTheoMass + productMassToleranceDa)
                            foundIons[i] = true;
                    }
                }
                //ZDOT IONS//
                if (ionsUsed.Contains(IonType.zdot))
                {
                    double zdotTheoMass = MassCalculator.MonoIsoptopicMass(candSeq.Substring(candSeq.Length - 1 - i, i + 1), out string error_message5) - Constants.nitrogenMonoisotopicMass - 2 * Constants.hydrogenMonoisotopicMass;
                    error_message += error_message5;
                    foreach (PTM ptm in psm.getCInfo().getPTMs())
                    {
                        if (ptm.index >= candSeq.Length - 2 - i)
                            zdotTheoMass += ptm.mass;
                    }
                    foreach (double expPeak in cPeaks)
                    {
                        if (expPeak > zdotTheoMass - productMassToleranceDa && expPeak < zdotTheoMass + productMassToleranceDa)
                            foundIons[foundIons.Count() - 2 - i] = true;
                    }
                }
            }
            //foundIons[0] = true; //AspN always starts with a D
            foundIons[foundIons.Count() - 1] = true;//A|B|C|D|E|F|K| where the whole peptide peak is always placed arbitrarly at the c term
        }

        public void ReadMassDictionary()
        {
            List<double> tempKeys = new List<double>();

            using (StreamReader masses = new StreamReader(Path.Combine(Environment.CurrentDirectory, "Dictionary" + maxMissingConsecutivePeaks + ".txt"))) //file located in Morpheus folder
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

        public bool GeneratePossibleSequences(PSM psm, out string error_message) //returns false if over the specified number of sequences are generated
        {
            error_message = "";
            List<string> foundSeq = new List<string>(); //get list of all FP sequences
            foreach (FusionCandidate fusionCandidate in psm.getFusionCandidates())
            {
                findIons(fusionCandidate, psm, out string error_message1); //populate the foundIons array
                error_message += error_message1;
                foundSeq.Add(fusionCandidate.seq);
            }
            bool done = false;
            int globalIndex = 0;
            while (!done)
            {
                done = true; //let's assume we're done and correct it later if we're not
                if (psm.getFusionCandidates().Count() > maxNumPossibleSequences) //if there are more than a set number of possible sequences, this is junk and we are not searching them all
                    return false;

                for (int fc = 0; fc < psm.getFusionCandidates().Count(); fc++)
                {
                    FusionCandidate fusionCandidate = psm.getFusionCandidates()[fc];
                    if (fusionCandidate.getFoundIons().Count() > globalIndex) //prevent crashing, use to tell when done by hitting end of fc
                    {
                        List<FusionCandidate> tempCandidates = new List<FusionCandidate>(); //fill with possible sequences

                        done = false; //We're not done, because at least one fusion candidate sequence length is still greater than the global index
                        string fusionSeq = fusionCandidate.seq;
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
                            double key = MassCalculator.MonoIsoptopicMass(ambiguousFrag, out string error_message2);
                            error_message += error_message2;

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
                                    if (closestPeak > key - productMassToleranceDa && closestPeak < key + productMassToleranceDa)
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
                                    if (closestPeak > key - productMassToleranceDa && closestPeak < key + productMassToleranceDa)
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
                            if (!foundSeq.Contains(newfc.seq)) //if new FP sequence, add it.
                            {
                                foundSeq.Add(newfc.seq);
                                findIons(newfc, psm, out string error_message3);
                                error_message += error_message3;
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

            //conduct an initial search of each candidate's full sequence to identify any that are translated
            for (int i = 0; i < psm.getFusionCandidates().Count(); i++) //foreach fusion peptide sequence that could map to this scan
            {
                string novelSeq = psm.getFusionCandidates()[i].seq;
                if (foundParent(novelSeq, ParentInfo.terminal.C, psm.getFusionCandidates()[i], false)) //check really quick to see if the whole thing exists as is. If so, assign it as translated. Terminal C was arbitrarily chosen
                {
                    foreach (ParentInfo info in psm.getFusionCandidates()[i].parentInfo)
                    {
                        foreach (TheoreticalProtein protein in info.theoreticalProteins)
                        {
                            if (protein.seq.Contains(novelSeq)) //if translated
                            {
                                psm.getFusionCandidates()[i].translatedParents.Add(new TranslatedParent(protein.id, protein.seq, protein.seq.IndexOf(psm.getFusionCandidates()[i].seq), psm.getFusionCandidates()[i].seq.Length));
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
            string novelSeq = tempCandidate.seq;
            //N//
            int nTermParentLength = novelSeq.Length - 1;//length-1, because if the whole thing existed we wouldn't have made it to the else loop. Several edits are made to reflect this
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
                    if (foundFirstSearch)
                    {
                        done = true;
                    }
                    nTermParentLength--;
                }
            }

            //C//
            done = false; //reset tracker
            foundFirstSearch = false;
            int cTermParentLength = novelSeq.Length - 1;
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
                    if (foundFirstSearch)
                    {
                        done = true;
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
                for (int junction = tempCandidate.seq.Length - cTermParentLength - 1; junction < nTermParentLength; junction++)
                {
                    tempCandidate.addJunctionIndex(junction);
                }
                return true;
            }
        }

        public bool foundParent(string frag, ParentInfo.terminal terminal, FusionCandidate candidate, bool foundFirstSearch)
        {
            //localTheoreticals.AsParallel().Where(x => x.Contains(frag)).ToList();
            if (notFoundSequences.Contains(frag)) //has the fragment been searched but not found before?
                return false;

            List<TheoreticalProtein> matches = new List<TheoreticalProtein>();
            if (foundSequences.TryGetValue(frag, out matches)) //has the fragment been searched AND found before?
            {
                candidate.parentInfo.Add(new ParentInfo(matches, terminal, frag));
                return true;
            }

            if (foundFirstSearch) //Has something smaller been found before? Well, then we can just search against those found sequences
            {
                string shorterFrag = terminal.Equals(ParentInfo.terminal.N) ? frag.Substring(0, frag.Length - 1) : frag.Substring(1, frag.Length - 1);

                foreach (ParentInfo info in candidate.parentInfo)
                {
                    if (info.parentType.Equals(terminal) && info.fragFound.Equals(shorterFrag))
                    {
                        List<TheoreticalProtein> tempProtList = new List<TheoreticalProtein>();
                        info.theoreticalProteins.ForEach(protein => tempProtList.Add(protein));
                        matches = tempProtList.AsParallel().Where(x => x.seq.Contains(frag)).ToList();
                    }
                }
            }
            else //it hasn't been found before... we need to search against the whole database :(
            {
                matches = theoreticalProteins.AsParallel().Where(x => x.seq.Contains(frag)).ToList();
            }
            if (matches != null && matches.Count()>0)
            {
                foundSequences.Add(frag, matches);
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
                string sequence = fusionCandidate.seq;
                foreach(ParentInfo info in fusionCandidate.parentInfo)
                {
                    int foundLength = info.fragFound.Length;

                    string compFrag = "";// fusionCandidate.seq.Substring()
                    if(info.parentType.Equals(ParentInfo.terminal.N))
                    {
                        compFrag = fusionCandidate.seq.Substring(foundLength, fusionCandidate.seq.Length - foundLength);
                    }
                    else
                    {
                        compFrag = fusionCandidate.seq.Substring(0, fusionCandidate.seq.Length - foundLength);
                    }

                    foreach (TheoreticalProtein protein in info.theoreticalProteins)
                    {
                        //get the index(es) of where the found fragment is
                        string subProt = protein.seq;
                        List<int> originalIndexes = new List<int>();
                        int pastIndex = 0;
                        while (subProt.Contains(info.fragFound))
                        {
                            int newIndex = subProt.IndexOf(info.fragFound);
                            originalIndexes.Add(newIndex + pastIndex);
                            subProt = subProt.Substring(newIndex + 1, subProt.Length - newIndex - 1); //need to remove old match
                            pastIndex += newIndex + 1;
                        }
                        fusionCandidate.transParents.Add(new TransParent(protein.id, protein.seq, originalIndexes, foundLength, info.parentType));

                        //get the index(es) of where the complimentary fragment is (if it's a cis fusion peptide)
                        subProt = protein.seq;
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
                                fusionCandidate.cisParents.Add(new CisParent(protein.id, protein.seq, originalIndexes, foundLength, complementaryIndexes, compFrag.Length));
                            else
                                fusionCandidate.cisParents.Add(new CisParent(protein.id, protein.seq, complementaryIndexes, compFrag.Length, originalIndexes, foundLength));
                        }
                    }
                }
            }
            else
            {
                string seq = fusionCandidate.seq;
                foreach (ParentInfo info in fusionCandidate.parentInfo)
                {
                    foreach (TheoreticalProtein protein in info.theoreticalProteins)
                    {
                        if (protein.seq.Contains(seq)) //if translated
                        {
                            fusionCandidate.translatedParents.Add(new TranslatedParent(protein.id, protein.seq, protein.seq.IndexOf(fusionCandidate.seq), fusionCandidate.seq.Length));
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
                    if (fc.getJunctionIndexes().Count() > fc.seq.Length) //delete it
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

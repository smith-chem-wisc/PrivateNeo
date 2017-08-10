using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using NeoInternal;

namespace Neo
{
    class SpliceFragments
    {
        public static string test = "";

        private Neo neo;
        private BackgroundWorker worker = null;
        public SpliceFragments(Neo neo, BackgroundWorker worker)
        {
            this.neo = neo;
            this.worker = worker;
        }

        public List<PSM> ExperimentalTheoreticalMatching(List<PSM> psms)
        {
            List<PSM> candidates = new List<PSM>();
            int counter = 0;
            foreach (PSM psm in psms)
            {
                this.worker.ReportProgress(Convert.ToInt16((Convert.ToDouble(counter) / Convert.ToDouble(psms.Count())) * 100));
                counter++;
                if (psm.getExpMass() < 20000) //some incorrect charge state assignments can yield huge calculated peptide masses causing StackOverflowExceptions.
                {
                    //preliminary filters can be removed if MassMatch calls to IonCrop are set to true.
                    string B = IonCrop(psm.getNInfo().seq, psm.getExpMass(), 0, IonType.b, true); //used as a preliminary filter to prevent longer searches from seq ID's that are larger than the precursor mass
                    string Y = IonCrop(psm.getCInfo().seq, psm.getExpMass(), 0, IonType.y, true); //used as a preliminary filter to prevent longer searches from seq ID's that are larger than the precursor mass
                    for (int y = 0; y < Y.Length - Neo.ionsUsedMassVer; y++) //foreach y aa removed
                    {
                        for (int b = 0; b < B.Length - Neo.ionsUsedMassVer; b++) //foreach b aa removed
                        {
                            MassMatch(B, Y, psm, b, y); //return a FusionCandidate
                        }
                    }
                    if (psm.getFusionCandidates().Count()>0) //match was found
                    {
                      //  MessageBox.Show(psm.getScan()+" "+psm.getFusionCandidates()[0].seq+" "+psm.getFusionCandidates().Count()+" "+psms.Count());
                        candidates.Add(psm);
                    }
                }
            }            
            return candidates;
        }

              //method was originally written recursively, but large peptides result in stackoverflow exceptions
        public void MassMatch(string B, string Y, PSM psm, int BIndex, int YIndex) //this is the workhorse of SpliceFragments
        {
            test = psm.getScan().ToString();
            double ExperimentalMass = psm.getExpMass();         
            string BFrag = IonCrop(B, ExperimentalMass, BIndex, IonType.b, false); //returns a B ion sequence that has a mass smaller than the experimental mass by cleaving C term AA
            //BIndex = B.Length - BFrag.Length; //added 11/8/16 Useful first pass to record how many AA have been cleaved from C term 
            string YFrag = IonCrop(Y, ExperimentalMass, YIndex, IonType.y, false); //returns a Y ion sequence that has a mass smaller than the experimental mass by cleaving N term AA
            //YIndex = Y.Length - YFrag.Length; //added 11/8/16 Useful first pass to record how many AA have been cleaved from N term
            double TheoreticalMass = MassCalculator.MonoIsoptopicMass(BFrag) + MassCalculator.MonoIsoptopicMass(YFrag) - Constants.WATER_MONOISOTOPIC_MASS + Neo.fixedModMass; //water added once in b and once in y

            //add PTM masses
            foreach (PTM ptm in psm.getNInfo().getPTMs())
            {
                if (ptm.index < BFrag.Length)
                {
                    TheoreticalMass += ptm.mass;
                }
            }
            foreach(PTM ptm in psm.getCInfo().getPTMs())
            {
                if (Y.Length-ptm.index < YFrag.Length)
                {
                    TheoreticalMass += ptm.mass;
                }
            }

            if (YFrag.Length < Neo.ionsUsedMassVer) //If the number of AA from the C-term peptide is less than desired amount, end recursion. 
            {
                //we're done
            }
            else if (BFrag.Length < Neo.ionsUsedMassVer) //If the number of AA from the N-term peptide is less than desired amount, start over loop and remove a single aa from the C-term
            {
             //   MassMatch(B, Y, psm, 0, YIndex+1);
            }
            //if match
            //bool elif = true; //"else if" where not a match==true
            else if (FalsePositives.generateDecoys)
            {
                //else if (((TheoreticalMass - Constants.PEPTIDE_N_TERMINAL_MONOISOTOPIC_MASS * 7) > (ExperimentalMass - 1 * Constants.PEPTIDE_N_TERMINAL_MONOISOTOPIC_MASS) && (TheoreticalMass - Constants.PEPTIDE_N_TERMINAL_MONOISOTOPIC_MASS * 7) < (ExperimentalMass + 1 * Constants.PEPTIDE_N_TERMINAL_MONOISOTOPIC_MASS)) | ((TheoreticalMass + Constants.PEPTIDE_N_TERMINAL_MONOISOTOPIC_MASS * 7) > (ExperimentalMass - 1 * Constants.PEPTIDE_N_TERMINAL_MONOISOTOPIC_MASS) && (TheoreticalMass + Constants.PEPTIDE_N_TERMINAL_MONOISOTOPIC_MASS * 7) < (ExperimentalMass + 1 * Constants.PEPTIDE_N_TERMINAL_MONOISOTOPIC_MASS)))//if match                                                                                                                                                                                                                                                                                                                        //if ((TheoreticalMass) > (ExperimentalMass+ i -PrecursorMassToleranceDa) && (TheoreticalMass) < (ExperimentalMass+i +PrecursorMassToleranceDa)) //if match
                if (((TheoreticalMass) > (ExperimentalMass - 9.5) && (TheoreticalMass) < (ExperimentalMass - 4.5)) | ((TheoreticalMass) > (ExperimentalMass + 5.5) && (TheoreticalMass) < (ExperimentalMass + 7.5)))//if match                          //if ((TheoreticalMass) > (ExperimentalMass+ i -PrecursorMassToleranceDa) && (TheoreticalMass) < (ExperimentalMass+i +PrecursorMassToleranceDa)) //if match                                                                                                                                                                                                                 //else if(Math.Abs(ExperimentalMass - TheoreticalMass)<40 && (ExperimentalMass - TheoreticalMass)-Math.Floor(ExperimentalMass-TheoreticalMass)>0.3 && (ExperimentalMass - TheoreticalMass) - Math.Floor(ExperimentalMass - TheoreticalMass) < 0.8)
                {
                   // elif = false;
                    bool previouslyFound = false;
                    foreach (FusionCandidate oldCandidate in psm.getFusionCandidates())
                    {
                        if ((BFrag + YFrag).Equals(oldCandidate.seq)) //see if that sequence was already recorded
                        {
                            previouslyFound = true;
                        }
                    }
                    if (!previouslyFound) //if fusion sequence was not previously assigned to this psm
                    {
                        FusionCandidate candidate = new FusionCandidate(BFrag + YFrag);
                        psm.addFusionCandidate(candidate);
                        //         MassMatch(B, Y, psm, BIndex + 1, YIndex);
                    }
                }
            }
            else
            {
                if ((TheoreticalMass) > (ExperimentalMass * (1 - FalsePositives.precursorMassTolerancePpm / 1000000)) && (TheoreticalMass) < (ExperimentalMass * (1 + FalsePositives.precursorMassTolerancePpm / 1000000))) //if match                          //if ((TheoreticalMass) > (ExperimentalMass+ i -PrecursorMassToleranceDa) && (TheoreticalMass) < (ExperimentalMass+i +PrecursorMassToleranceDa)) //if match
                {
                   // elif = false;
                    bool previouslyFound = false;
                    foreach (FusionCandidate oldCandidate in psm.getFusionCandidates())
                    {
                        if ((BFrag + YFrag).Equals(oldCandidate.seq)) //see if that sequence was already recorded
                        {
                            previouslyFound = true;
                        }
                    }
                    if (!previouslyFound) //if fusion sequence was not previously assigned to this psm
                    {
                        FusionCandidate candidate = new FusionCandidate(BFrag + YFrag);
                        psm.addFusionCandidate(candidate);
                        //         MassMatch(B, Y, psm, BIndex + 1, YIndex);
                    }
                }
            }
           // if(elif) //not a match
            {
                /*        if (TheoreticalMass < ExperimentalMass && BIndex == 0) //first pass, theo less than exp and can't take away more ions
                        {
                            //we're done
                        }
                        else
                        {
                            if (TheoreticalMass < ExperimentalMass) //if b out of ions, but y not, crop off a y and start again
                            {
                                BIndex = 0;
                                YIndex++;
                                MassMatch(B,Y, psm, BIndex, YIndex);
                            }
                            else
                            { //crop off a b ion
                                MassMatch(B, Y, psm, BIndex + 1, YIndex);
                            }
                        }*/
            }
        }

        //This method removes the number of amino acids specified by FragNumber from the respecitve terminus specified by ion of IonSequence
        //If checkToRemoveExtraAA is true, additional AA will be removed to achieve a theoretical mass less than the Experimental mass
        public string IonCrop(string IonSequence, double ExperimentalMass, int FragNumber, IonType ion, bool checkToRemoveExtraAA)
        {
            string IonFrag;
            if (ion == IonType.b)
            {
                IonFrag = IonSequence.Substring(0, (IonSequence.Length - FragNumber));
                if (IonFrag.Substring(IonSequence.Length - FragNumber - 1, 1) == ")") //if end of a PTM annotation
                {
                    while (IonFrag.Substring(IonSequence.Length - FragNumber - 1, 1) != "(")
                    {
                        FragNumber++;
                        IonFrag = IonSequence.Substring(0, (IonSequence.Length - FragNumber));
                    }
                    FragNumber++; //removes "("
                    FragNumber++; //removes the AA the PTM was attached to
                    IonFrag = IonSequence.Substring(0, (IonSequence.Length - FragNumber));
                }
            }
            else //Ion==Y
            {
                IonFrag = IonSequence.Substring((0 + FragNumber), (IonSequence.Length - FragNumber));
                if (IonFrag.Substring(0, 1) == "(") //if start of a PTM annotation
                {
                    while (IonFrag.Substring(0, 1) != ")")
                    {
                        FragNumber++;
                        IonFrag = IonSequence.Substring((0 + FragNumber), (IonSequence.Length - FragNumber));
                    }
                    FragNumber++; //removes ")"
                    IonFrag = IonSequence.Substring((0 + FragNumber), (IonSequence.Length - FragNumber));
                }
            }
            if (checkToRemoveExtraAA == false)
            {
                return IonFrag;
            }
            else {
                double IonMass = MassCalculator.MonoIsoptopicMass(IonFrag);
                if (IonMass < ExperimentalMass) //end if the mass of the fragment is lower than the experimental
                {
                    return IonFrag;
                }
                else //call the function again to remove another amino acid.
                {
                    FragNumber++;
                    return IonCrop(IonSequence, ExperimentalMass, FragNumber, ion, true);
                }
            }
        }

    }
}

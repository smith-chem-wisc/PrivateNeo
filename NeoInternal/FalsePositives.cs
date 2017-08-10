using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;

namespace NeoInternal
{
    public class FalsePositives
    {
        private BackgroundWorker worker = null;
        public static bool generateDecoys;
        public static readonly int ionsUsedDigFilter = 6;
        public static double precursorMassTolerancePpm; //(ppm)

        public FalsePositives(BackgroundWorker worker)
        {
            this.worker = worker;
        }        
        
        public void FindCommonFalsePositives(List<PSM> psms, List<TheoreticalProtein> database)
        {
            int i = 0;
            foreach(PSM psm in psms)
            {
                this.worker.ReportProgress(Convert.ToInt16((Convert.ToDouble(i) / Convert.ToDouble(psms.Count())) * 100));
                i++;
                foreach (FusionCandidate fusionCandidate in psm.getFusionCandidates())
                {
                    foreach(ParentInfo info in fusionCandidate.parentInfo)
                    {
                        if(info.fragFound.Length>=6)
                        {
                            foreach(TheoreticalProtein protein in info.theoreticalProteins)
                            {
                                string protSeq = protein.seq;
                                char[] candidateSeq = fusionCandidate.seq.ToCharArray();
                                int index = protSeq.IndexOf(info.fragFound);
                                int fragLength = info.fragFound.Length;
                                string possibleTranslatedSequence = protSeq.Substring(index, fragLength);
                                if(!possibleTranslatedSequence.Equals(fusionCandidate.seq)) //if not already found as translated
                                { 
                                    double fragMass = MassCalculator.MonoIsoptopicMass(possibleTranslatedSequence);
                                    double expMass = psm.getExpMass();
                                    while (fragMass < expMass + 187.079 - 57.021+1)
                                    {
                                        //Find SNPs
                                        if (candidateSeq.Count() == possibleTranslatedSequence.Length)
                                        {
                                            char[] possibleSeqArray = possibleTranslatedSequence.ToCharArray();
                                            if (IsSNP(candidateSeq, possibleSeqArray))
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence, Variant.variantType.SNP, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                                psm.variants.Add(new Variant(fusionCandidate.seq, Variant.variantType.SNP, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                        }

                                        if (generateDecoys)
                                        {
                                            //search for unmodified sequences

                                            if (((fragMass) > (expMass - 9.5) && (fragMass) < (expMass - 4.5)) | ((fragMass) > (expMass + 5.5) && (fragMass) < (expMass + 7.5))) //if match, add it!
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence, Variant.variantType.UM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            //fragMass, expMass
                                            //search for PTMs
                                            else if (((fragMass + .984016) > (expMass - 9.5) && (fragMass + .98402) < (expMass - 4.5)) | ((fragMass + .98402) > (expMass + 5.5) && (fragMass + .98402) < (expMass + 7.5))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Deamidation", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if (((fragMass + 15.99491463) > (expMass - 9.5) && (fragMass + 15.99491463) < (expMass - 4.5)) | ((fragMass + 15.99491463) > (expMass + 5.5) && (fragMass + 15.99491463) < (expMass + 7.5))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Oxidation", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if (((fragMass + 15.99491463 * 2) > (expMass - 9.5) && (fragMass + 15.99491463 * 2) < (expMass - 4.5)) | ((fragMass + 15.99491463 * 2) > (expMass + 5.5) && (fragMass + 15.99491463 * 2) < (expMass + 7.5))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+DiOxidation", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if (((fragMass + 15.99491463 * 3) > (expMass - 9.5) && (fragMass + 15.99491463 * 3) < (expMass - 4.5)) | ((fragMass + 15.99491463 * 3) > (expMass + 5.5) && (fragMass + 15.99491463 * 3) < (expMass + 7.5))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+TriOxidation", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            //else if ((fragMass + 14.01565) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + 14.01565) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            //{
                                            //    psm.variants.Add(new Variant(Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            //}                                            else if (((TheoreticalMass + .984016) > (ExperimentalMass - 9.5) && (TheoreticalMass + .98402) < (ExperimentalMass - 4.5)) | ((TheoreticalMass + .98402) > (ExperimentalMass + 5.5) && (TheoreticalMass) < (ExperimentalMass + 7.5))) //if match
                                            else if (((fragMass + 42.01056) > (expMass - 9.5) && (fragMass + 42.01056) < (expMass - 4.5)) | ((fragMass + 42.01056) > (expMass + 5.5) && (fragMass + 42.01056) < (expMass + 7.5))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Acetylation", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if (((fragMass + 79.96633) > (expMass - 9.5) && (fragMass + 79.96633) < (expMass - 4.5)) | ((fragMass + 79.96633) > (expMass + 5.5) && (fragMass + 79.96633) < (expMass + 7.5))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Phosphorylation", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                        }
                                        else
                                        {
                                            //search for unmodified sequences
                                            if ((fragMass) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match, add it!
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence, Variant.variantType.UM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }

                                            //search for PTMs
                                            else if ((fragMass + .984016) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + .98402) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Deamidation", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if ((fragMass + 15.99491463) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + 15.99491) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Oxidation", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if ((fragMass + 15.99491463 * 2) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + 15.99491 * 2) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Dioxidation", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if ((fragMass + 15.99491463 * 3) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + 15.99491 * 3) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Trioxidation", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            //else if ((fragMass + 14.01565) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + 14.01565) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            //{
                                            //    psm.variants.Add(new Variant(Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            //}
                                            else if ((fragMass + 42.01056) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + 42.01056) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Acetyl", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if ((fragMass + 79.96633) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + 79.96633) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Phospho", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if ((fragMass + 61.913495) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + 61.913495) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Zinc", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if ((fragMass + 37.955588) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + 37.955588) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Potassium", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if ((fragMass + 21.981944) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + 21.981944) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Sodium", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if ((fragMass + 57.021464) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + 57.021464) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Carbamidomethyl", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if ((fragMass + 79.956815) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + 79.956815) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Sulfonation", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if ((fragMass + 14.01565) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + 14.01565) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Methyl", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if ((fragMass + 28.0313) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + 28.0313) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+DiMethyl", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if ((fragMass + 42.04695) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass + 42.04695) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+TriMethyl", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                            else if ((fragMass - -17.026549) > (expMass * (1 - precursorMassTolerancePpm / 1000000)) && (fragMass - -17.026549) < (expMass * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                            {
                                                psm.variants.Add(new Variant(possibleTranslatedSequence + "+Ammonia loss", Variant.variantType.PTM, protein.id, protSeq, protSeq.IndexOf(info.fragFound), info.fragFound.Length));
                                            }
                                        }
                                        //-17, 22, 

                                        //update the sequence
                                        fragLength++;
                                        if (info.parentType.Equals(ParentInfo.terminal.N))
                                        {
                                            if (index + fragLength < protSeq.Length)
                                            {
                                                possibleTranslatedSequence = protSeq.Substring(index, fragLength);
                                                fragMass = MassCalculator.MonoIsoptopicMass(possibleTranslatedSequence);
                                            }
                                            else
                                            {
                                                fragMass = expMass + 187.079 - 57.021 + 1;
                                            }

                                        }
                                        else
                                        {
                                            index--;
                                            if (index >= 0)
                                            {
                                                possibleTranslatedSequence = protSeq.Substring(index, fragLength);
                                                fragMass = MassCalculator.MonoIsoptopicMass(possibleTranslatedSequence);
                                            }
                                            else
                                            {
                                                fragMass = expMass + 187.079 - 57.021 + 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                List<Variant> variants = psm.variants;
                for (int v = variants.Count() - 1; v >= 0; v--) 
                {
                    for (int v2 = 0; v2 < v; v2++) 
                    {
                        if(variants[v].pepSeq.Equals(variants[v2].pepSeq) && variants[v].id.Equals(variants[v2].id))
                        {
                            variants.Remove(variants[v]);
                            v2 = variants.Count();
                        }
                    }
                }
            }
        }

        public bool IsSNP(char[] candidateSeq, char[] possibleSeqArray)
        {
            int missedAA = 0;
            for (int i = 0; i < candidateSeq.Count(); i++)
            {
                if (!candidateSeq[i].Equals(possibleSeqArray[i]))
                {
                    missedAA++;
                }
            }
            if (missedAA == 1)
            {
                return true;
            }
            return false;
        }

        //compare the 6 first and last aa of each fusion candidate with database and determine if precursor mass can be achieved within 5 ppm. If it can, remove the psm from the list
        public void removeTranslatedPeptides(List<PSM> psms, List<TheoreticalProtein> database)
        {
            for (int i = 0; i < psms.Count(); i++)
            {

                this.worker.ReportProgress(Convert.ToInt16((Convert.ToDouble(i) / Convert.ToDouble(psms.Count())) * 100));
                bool removed = false;
                int fcIndex = 0;//FusionCandidate
                while (!removed && fcIndex < psms[i].getFusionCandidates().Count())
                {
                    FusionCandidate fc = psms[i].getFusionCandidates()[fcIndex];
                    string seq = fc.seq;
                    if (seq.Length >= ionsUsedDigFilter)
                    {
                        string Nterm = seq.Substring(0, ionsUsedDigFilter);
                        string Cterm = seq.Substring(seq.Length - ionsUsedDigFilter, ionsUsedDigFilter);


                        //N-TERMINUS SEARCHING
                        List<TheoreticalProtein> matches = database.AsParallel().Where(x => x.seq.Contains(Nterm)).ToList();
                        int protIndex = 0; //TheoreticalProtein
                        while (!removed && protIndex < matches.Count())
                        {
                            TheoreticalProtein prot = matches[protIndex];
                            double prodMass = 0;
                            //use below code to catch multiple appearances of a frag in a parent protein
                            List<int> indexes = new List<int>();
                            string subProt = prot.seq;
                            int pastIndex = 0;
                            while (subProt.Contains(Nterm))
                            {
                                int newIndex = subProt.IndexOf(Nterm);
                                indexes.Add(newIndex + pastIndex);
                                subProt = subProt.Substring(newIndex + 1, subProt.Length - newIndex - 1); //need to remove old match
                                pastIndex += newIndex + 1;
                            }
                            string protProd = "";
                            int sameProtIndex = 0;
                            while (!removed && sameProtIndex < indexes.Count())
                            {
                                Boolean hitEndOfProt = false;
                                int numAAused = ionsUsedDigFilter;
                                while (prodMass < psms[i].getExpMass() + 200 && hitEndOfProt == false)
                                {
                                    try
                                    {
                                        //   for (int PTMInclude = 0; PTMInclude < PTMIncludeMax; PTMInclude++) //loop twice, once without incorporating identified ptm masses and once with
                                        {
                                            //Obtain the Product mass
                                            //if (Ion == 0) //if B main
                                            {
                                                protProd = prot.seq.Substring(indexes[sameProtIndex], numAAused);
                                                prodMass = MassCalculator.MonoIsoptopicMass(protProd);
                                                //MessageBox.Show("1 "+ProtProd+" "+FASTARow[0].ToString() + " " + ProdMass.ToString());
                                                /*        if (PTMInclude == 1)
                                                        {
                                                            ProdMass += BPTMMass;
                                                        }*/
                                            }
                                            //MissedCleavage and NonSpecific Cleavage/autolysis catch
                                            if(generateDecoys)
                                            {
                                                if (((prodMass) > (psms[i].getExpMass() - 9.5) && (prodMass) < (psms[i].getExpMass() - 4.5)) | ((prodMass) > (psms[i].getExpMass() + 5.5) && (prodMass) < (psms[i].getExpMass() + 7.5))) //if match, add it!
                                                {
                                                    psms.Remove(psms[i]);
                                                    i--;
                                                    removed = true;
                                                }
                                            }
                                            else
                                            {
                                                if ((prodMass) > (psms[i].getExpMass() * (1 - precursorMassTolerancePpm / 1000000)) && (prodMass) < (psms[i].getExpMass() * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                                {
                                                    psms.Remove(psms[i]);
                                                    i--;
                                                    removed = true;
                                                }
                                            }
                                        }
                                    }
                                    catch //sloppy patch used for when hitting the end of a protein
                                    {
                                        hitEndOfProt = true;
                                    }
                                    numAAused++;
                                }
                                sameProtIndex++;
                            }
                            protIndex++;
                        }


                        //C-TERMINUS SEARCHING
                        matches = database.AsParallel().Where(x => x.seq.Contains(Cterm)).ToList();
                        protIndex = 0; //TheoreticalProtein
                        while (!removed && protIndex < matches.Count())
                        {
                            TheoreticalProtein prot = matches[protIndex];
                            double prodMass = 0;
                            //use below code to catch multiple appearances of a frag in a parent protein
                            List<int> indexes = new List<int>();
                            string subProt = prot.seq;
                            int pastIndex = 0;
                            while (subProt.Contains(Cterm))
                            {
                                int newIndex = subProt.IndexOf(Cterm);
                                indexes.Add(newIndex + pastIndex);
                                subProt = subProt.Substring(newIndex + 1, subProt.Length - newIndex - 1); //need to remove old match
                                pastIndex += newIndex + 1;
                            }
                            string protProd = "";
                            int sameProtIndex = 0;
                            while (!removed && sameProtIndex < indexes.Count())
                            {
                                Boolean hitEndOfProt = false;
                                int numAAused = ionsUsedDigFilter;
                                while (prodMass < psms[i].getExpMass() + 200 && hitEndOfProt == false)
                                {
                                    try
                                    {
                                        //   for (int PTMInclude = 0; PTMInclude < PTMIncludeMax; PTMInclude++) //loop twice, once without incorporating identified ptm masses and once with
                                        {
                                            //Obtain the Product mass
                                            //else //if y main
                                            {
                                                protProd = prot.seq.Substring(indexes[sameProtIndex] - numAAused + ionsUsedDigFilter, numAAused);
                                                //MessageBox.Show(BProt);
                                                prodMass = MassCalculator.MonoIsoptopicMass(protProd);
                                                //MessageBox.Show(Prot.IndexOf(bigFrag).ToString() + " " + AASearchLength + " "+bigFrag.Count().ToString());
                                                //MessageBox.Show("2 "+bigFrag+" "+ProtProd + " " + ProdMass.ToString());
                                                /*      if (PTMInclude == 1)
                                                      {
                                                          ProdMass += YPTMMass;
                                                      }*/
                                            }
                                            //MissedCleavage and NonSpecific Cleavage/autolysis catch
                                            if(generateDecoys)
                                            {
                                                if (((prodMass) > (psms[i].getExpMass() - 9.5) && (prodMass) < (psms[i].getExpMass() - 4.5)) | ((prodMass) > (psms[i].getExpMass() + 5.5) && (prodMass) < (psms[i].getExpMass() + 7.5))) //if match, add it!
                                                {
                                                    psms.Remove(psms[i]);
                                                    i--;
                                                    removed = true;
                                                }
                                            }
                                            else
                                            {
                                                if ((prodMass) > (psms[i].getExpMass() * (1 - precursorMassTolerancePpm / 1000000)) && (prodMass) < (psms[i].getExpMass() * (1 + precursorMassTolerancePpm / 1000000))) //if match
                                                {
                                                    psms.Remove(psms[i]);
                                                    i--;
                                                    removed = true;
                                                }
                                            }
                                        }
                                    }
                                    catch //sloppy patch used for when hitting the end of a protein
                                    {
                                        hitEndOfProt = true;
                                    }
                                    numAAused++;
                                }
                                sameProtIndex++;
                            }
                            protIndex++;
                        }
                    }
                    fcIndex++;
                }
            }
        }
    }
}

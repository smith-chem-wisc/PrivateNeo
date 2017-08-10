namespace Neo
{
    class FalsePositiveOld
    {

        /*        private void CompareDigestionProducts() //this should be broken up into multiple methods
        {
            DigestCandidates.Columns.Add("BProd", typeof(string));
            DigestCandidates.Columns.Add("YProd", typeof(string));
            DigestCandidates.Columns.Add("Type", typeof(string));
            //DataRow[] foundRows;
            string[] DigestProduct = { "Z", "Z", "Z" };
            //   string BProt = "";
            // string YProt = "";
            bool Done = false;
            double ExpMass;
            string BPep;
            string YPep;
            //string BProd;
            //string YProd;
            //double BProdMass;
            //double YProdMass;
            double BPTMMass;
            double YPTMMass;

            for (int i = 0; i < DigestCandidates.Rows.Count; i++) //for each spectra 
            { //if it's a candidate
                if (DigestCandidates.Rows[i][5].ToString() != "Z") //&& DigestCandidates.Rows[i][1].ToString().Substring(0,5)!="DECOY" && DigestCandidates.Rows[i][3].ToString().Substring(0, 5) != "DECOY")
                {
                    DigestProduct[0] = "Z";
                    DigestProduct[1] = "Z";
                    DigestProduct[2] = "Z";

                    // MessageBox.Show(DigestCandidates.Rows[i][5].ToString());
                    BPep = DigestCandidates.Rows[i][0].ToString(); //get sequence of BFrag
                    BPTMMass = MonoIsoptopicMass(BPep);//get mass of BFrag
                    BPep = Regex.Replace(BPep, "(\\(.*?\\))", "");//remove PTM annotations
                    BPTMMass = BPTMMass - MonoIsoptopicMass(BPep);//calculate mass of PTMs removed-this will be used as an added mass for database comparisons
                    DigestCandidates.Rows[i][0] = BPep;//Replace PTM annotated sequence with PTM free sequence for downstream FASTA searching

                    //Do same with Y ID
                    YPep = DigestCandidates.Rows[i][2].ToString();
                    YPTMMass = MonoIsoptopicMass(YPep);
                    YPep = Regex.Replace(YPep, "(\\(.*?\\))", "");
                    YPTMMass = YPTMMass - MonoIsoptopicMass(YPep);
                    DigestCandidates.Rows[i][2] = YPep;

                    string[] fraglr = DigestCandidates.Rows[i][5].ToString().Split('|').ToArray(); //separate multiple potential fusion junctions
                    string[] BandYMains = { fraglr[fraglr.Count() - 2], fraglr[0] }; //just need last and first. Use -2 because each entry ends with a "|", producing a new (blank) entry. Big B first, then big Y
                    for (int Ion = 0; Ion < BandYMains.Count(); Ion++) //Ion specifies frag, where 0 is B and 1 is Y for the majority. This is used so that ABCD_EFGHIJ and ABCDEFGH_IJ doesn't search fasta for -IJ, which saves time
                    {
                        //begin searching for additional peptide sources
                        ExpMass = Convert.ToDouble(DigestCandidates.Rows[i][4])-fixedModMass;
                        String Prot = "";
                        String ProtProd = "";
                        double ProdMass = 0.0;
                        //interested in finding if one of the original IDs corresponded to a matching mass
                        string bigFrag = "";

                       /* if (Ion == 0)
                        {
                            bigFrag = BPep.Substring(0, IonsUsedDigFilter);
                        }
                        else
                        {
                            bigFrag = YPep.Substring(YPep.Count() - IonsUsedDigFilter, IonsUsedDigFilter);
                        }*/
     /*   bigFrag = BandYMains[Ion].Split('_').ToArray()[Ion];



                        //search for missed cleavage/autolysis/sequence in db (where first four aa match an ID)
                        if (bigFrag.Count() >= IonsUsedDigFilter) //optimizes false negatives
                        {
                            int PTMIncludeMax = 2;
                            if (Ion == 0)
                            {
                                if (BPTMMass > -0.5 && BPTMMass< 0.5)
                                {
                                    PTMIncludeMax = 1;
                                }
}
                            else
                            {
                                if (YPTMMass > -0.5 && YPTMMass< 0.5)
                                {
                                    PTMIncludeMax = 1;
                                }
                            }
                            int rowNum = 0;
                            while (rowNum<FASTA.Rows.Count && Done == false) //foreach fasta entry (in the ion loop, so search b and y frags separately)
                            {
                                if (FASTA.Rows[rowNum][1].ToString().Contains(bigFrag)) //if it contains the frag of interest 
                                {
                                    ProdMass = 0;
                                    Prot = FASTA.Rows[rowNum][1].ToString();
int numAAUsed = IonsUsedDigFilter;
Boolean hitEndOfProt = false;
//use below code to catch multiple appearances of a frag in a parent protein
List<int> indexes = new List<int>();
string subProt = Prot;
int pastIndex = 0;
                                    while (subProt.Contains(bigFrag))
                                    {
                                        int newIndex = subProt.IndexOf(bigFrag);
indexes.Add(newIndex + pastIndex);
                                        subProt = subProt.Substring(newIndex + 1, subProt.Length - newIndex - 1); //need to remove old match
                                        pastIndex += newIndex + 1;
                                    }
                                    foreach (int index in indexes)
                                    {
                                        while (ProdMass<ExpMass + 200 && hitEndOfProt == false)
                                        {
                                            try
                                            {
                                                for (int PTMInclude = 0; PTMInclude<PTMIncludeMax; PTMInclude++) //loop twice, once without incorporating identified ptm masses and once with
                                                {
                                                    //Obtain the Product mass
                                                    if (Ion == 0) //if B main
                                                    {
                                                        ProtProd = Prot.Substring(index, numAAUsed);
                                                        ProdMass = MonoIsoptopicMass(ProtProd);
                                                        //MessageBox.Show("1 "+ProtProd+" "+FASTARow[0].ToString() + " " + ProdMass.ToString());
                                                        if (PTMInclude == 1)
                                                        {
                                                            ProdMass += BPTMMass;
                                                        }
                                                    }

                                                    else //if y main
                                                    {
                                                        ProtProd = Prot.Substring(index - numAAUsed, numAAUsed);
                                                        //MessageBox.Show(BProt);
                                                        ProdMass = MonoIsoptopicMass(ProtProd);
                                                        //MessageBox.Show(Prot.IndexOf(bigFrag).ToString() + " " + AASearchLength + " "+bigFrag.Count().ToString());
                                                        //MessageBox.Show("2 "+bigFrag+" "+ProtProd + " " + ProdMass.ToString());
                                                        if (PTMInclude == 1)
                                                        {
                                                            ProdMass += YPTMMass;
                                                        }
                                                    }

                                                    //MissedCleavage and NonSpecific Cleavage/autolysis catch
                                                    if ((ProdMass) > (ExpMass* (1 - PrecursorMassTolerancePpm / 1000000)) && (ProdMass) < (ExpMass* (1 + PrecursorMassTolerancePpm / 1000000))) //if match
                                                    {
                                                        DigestProduct[Ion] = ProtProd;
                                                        Done = true;
                                                    }
                                                }
                                            }
                                            catch //sloppy patch used for when hitting the end of a protein
                                            {
                                                hitEndOfProt = true;
                                            }
                                            numAAUsed++;
                                        }
                                        numAAUsed = IonsUsedDigFilter; //reset the counter
                                    }
                                }
                                rowNum++;
                            }
                        }*/
                        /*
                        //Search for PTMs that might shift the mass of a db sequence (where first four aa match an ID) to match the expMass
                        rowNum = 0;
                        while (rowNum < FASTA.Rows.Count && Done == false) //foreach fasta entry (in the ion loop, so search b and y frags separately)
                        {
                            if (FASTA.Rows[rowNum][1].ToString().Contains(bigFrag)) //if it contains the frag of interest 
                            {
                                ProdMass = 0;
                                Prot = FASTA.Rows[rowNum][1].ToString();
                                int numAAUsed = IonsUsedDigFilter;
                                Boolean hitEndOfProt = false;                                
                                //use below code to catch multiple appearances of a frag in a parent protein
                                List<int> indexes = new List<int>();
                                string subProt = Prot;
                                int pastIndex = 0;
                                while (subProt.Contains(bigFrag))
                                {
                                    int newIndex = subProt.IndexOf(bigFrag);
                                    indexes.Add(newIndex + pastIndex);
                                    subProt = subProt.Substring(newIndex + 1, subProt.Length - newIndex - 1); //need to remove old match
                                    pastIndex += newIndex + 1;
                                }
                                foreach (int index in indexes)
                                {
                                    while (ProdMass < ExpMass + 200 && hitEndOfProt == false)
                                    {
                                        try
                                        {
                                            for (int PTMInclude = 0; PTMInclude < 2; PTMInclude++) //loop twice, once without incorporating identified ptm masses and once with
                                            {
                                                //Obtain the Product mass
                                                if (Ion == 0) //if B main
                                                {
                                                    ProtProd = Prot.Substring(index, IonsUsedDigFilter + numAAUsed); //what if it appears multiple times in a protein? FIXME http://stackoverflow.com/questions/6865419/indexof-for-multiple-results
                                                    ProdMass = MonoIsoptopicMass(ProtProd);
                                                    //MessageBox.Show("1 "+ProtProd+" "+FASTARow[0].ToString() + " " + ProdMass.ToString());
                                                    if (PTMInclude == 1)
                                                    {
                                                        ProdMass += BPTMMass;
                                                    }
                                                    if (BPTMMass > -0.5 && BPTMMass < 0.5)
                                                    {
                                                        PTMInclude = 2;
                                                    }
                                                }
                                                else //if y main
                                                {
                                                    ProtProd = Prot.Substring(index - numAAUsed, IonsUsedDigFilter + numAAUsed);
                                                    //MessageBox.Show(BProt);
                                                    ProdMass = MonoIsoptopicMass(ProtProd);
                                                    //MessageBox.Show(Prot.IndexOf(bigFrag).ToString() + " " + AASearchLength + " "+bigFrag.Count().ToString());
                                                    //MessageBox.Show("2 "+bigFrag+" "+ProtProd + " " + ProdMass.ToString());
                                                    if (PTMInclude == 1)
                                                    {
                                                        ProdMass += YPTMMass;
                                                    }
                                                    if (YPTMMass > -0.5 && YPTMMass < 0.5)
                                                    {
                                                        PTMInclude = 2;
                                                    }
                                                }
                                                //PTMDISCOVERY     // Variants.Rows.Add(BProdMass - ExpMass, BProd, ExpMass, 'B'); //recording every comparison; not necessary but can be used to manually look for trends in dMass
                                                //MessageBox.Show(BProd + " " + BProdMass.ToString() + " " + ExpMass.ToString());
                                                if ((ProdMass + .984016) > (ExpMass * (1 - PrecursorMassTolerancePpm / 1000000)) && (ProdMass + .98402) < (ExpMass * (1 + PrecursorMassTolerancePpm / 1000000))) //if match
                                                {
                                                    DigestProduct[Ion] = ProtProd + " + Deamidation";
                                                    Done = true;
                                                }
                                                else if ((ProdMass + 15.99491463) > (ExpMass * (1 - PrecursorMassTolerancePpm / 1000000)) && (ProdMass + 15.99491) < (ExpMass * (1 + PrecursorMassTolerancePpm / 1000000))) //if match
                                                {
                                                    DigestProduct[Ion] = ProtProd + " + Oxidation";
                                                    Done = true;
                                                }
                                                else if ((ProdMass + 15.99491463 * 2) > (ExpMass * (1 - PrecursorMassTolerancePpm / 1000000)) && (ProdMass + 15.99491 * 2) < (ExpMass * (1 + PrecursorMassTolerancePpm / 1000000))) //if match
                                                {
                                                    DigestProduct[Ion] = ProtProd + " + DiOxidation";
                                                    Done = true;
                                                }
                                                else if ((ProdMass + 15.99491463 * 3) > (ExpMass * (1 - PrecursorMassTolerancePpm / 1000000)) && (ProdMass + 15.99491 * 3) < (ExpMass * (1 + PrecursorMassTolerancePpm / 1000000))) //if match
                                                {
                                                    DigestProduct[Ion] = ProtProd + " + TriOxidation";
                                                    Done = true;
                                                }
                                                else if ((ProdMass + 14.01565) > (ExpMass * (1 - PrecursorMassTolerancePpm / 1000000)) && (ProdMass + 14.01565) < (ExpMass * (1 + PrecursorMassTolerancePpm / 1000000))) //if match
                                                {
                                                    DigestProduct[Ion] = ProtProd + " + Methylation";
                                                    Done = true;
                                                }
                                                else if ((ProdMass + 42.01056) > (ExpMass * (1 - PrecursorMassTolerancePpm / 1000000)) && (ProdMass + 42.01056) < (ExpMass * (1 + PrecursorMassTolerancePpm / 1000000))) //if match
                                                {
                                                    DigestProduct[Ion] = ProtProd + " + Acetylation";
                                                    Done = true;
                                                }
                                                else if ((ProdMass + 79.96633) > (ExpMass * (1 - PrecursorMassTolerancePpm / 1000000)) && (ProdMass + 79.96633) < (ExpMass * (1 + PrecursorMassTolerancePpm / 1000000))) //if match
                                                {
                                                    DigestProduct[Ion] = ProtProd + " + Phosphorylation";
                                                    Done = true;
                                                }
                                                else if (ProtProd.Contains("Q") && (ProdMass - 17.02655) > (ExpMass * (1 - PrecursorMassTolerancePpm / 1000000)) && (ProdMass - 17.02655) < (ExpMass * (1 + PrecursorMassTolerancePpm / 1000000))) //if match)
                                                {
                                                    DigestProduct[Ion] = ProtProd + " + PyroGlutamic Acid";
                                                    Done = true;
                                                }
                                                else if (ProtProd.Contains("E") && (ProdMass - 18.01056) > (ExpMass * (1 - PrecursorMassTolerancePpm / 1000000)) && (ProdMass - 18.01056) < (ExpMass * (1 + PrecursorMassTolerancePpm / 1000000))) //if match)
                                                {
                                                    DigestProduct[Ion] = ProtProd + " + PyroGlutamic Acid";
                                                    Done = true;
                                                }
                                                numAAUsed++;
                                                if (Done) //no reason to continue the search if a viable digest has already been found
                                                {
                                                    PTMInclude = 2;
                                                    hitEndOfProt = true;
                                                }
                                            }
                                        }
                                        catch //sloppy patch used for when hitting the end of a protein
                                        {
                                            hitEndOfProt = true;
                                        }
                                    }
                                    numAAUsed = IonsUsedDigFilter; //reset the counter
                                }
                            }
                            rowNum++;
                        }

                        //SNPDISCOVERY 
                        rowNum = 0;
                        while (rowNum < FASTA.Rows.Count && Done == false) //foreach fasta entry (in the ion loop, so search b and y frags separately) or until a match is found
                        {
                            if (FASTA.Rows[rowNum][1].ToString().Contains(bigFrag)) //if it contains the frag of interest 
                            {
                                //MessageBox.Show(rowNum.ToString() + " " + FASTA.Rows.Count.ToString());
                                Prot = FASTA.Rows[rowNum][1].ToString(); //get protein seq
                                //int numAAUsed = 1;  //Number of additional AA
                                //use below code to catch multiple appearances of a frag in a parent protein
                                List<int> indexes = new List<int>();
                                string subProt = Prot;
                                int pastIndex = 0;
                                while (subProt.Contains(bigFrag))
                                {
                                    int newIndex = subProt.IndexOf(bigFrag);
                                    indexes.Add(newIndex + pastIndex);
                                    subProt = subProt.Substring(newIndex + 1, subProt.Length - newIndex - 1); //need to remove old match
                                    pastIndex += newIndex + 1;
                                }
                                foreach (int index in indexes)
                                {
                                    string[] similarSequenceLengths = new string[4];
                                    for (int PTMInclude = 0; PTMInclude < 2; PTMInclude++) //loop twice, once without incorporating identified ptm masses and once with
                                    {
                                        int numAA = 0;
                                        Boolean overExpMass = false;
                                        if (Ion == 0) //if B main
                                        {
                                            bigFrag = BPep.Substring(0, IonsUsedDigFilter);
                                            while (index + numAA + IonsUsedDigFilter < Prot.Length && overExpMass == false)
                                            {
                                                ProtProd = Prot.Substring(index, IonsUsedDigFilter + numAA); //what if it appears multiple times in a protein? FIXME http://stackoverflow.com/questions/6865419/indexof-for-multiple-results
                                                ProdMass = MonoIsoptopicMass(ProtProd);                                                               //additionProd = ProtProd.Substring(fragIndex+IonsUsedDigFilter,numAAUsed); //causing crash
                                                if (PTMInclude == 1)
                                                {
                                                    ProdMass += BPTMMass;
                                                }
                                                if (ProdMass > ExpMass)
                                                {
                                                    overExpMass = true;
                                                }
                                                numAA++;
                                            }
                                            similarSequenceLengths[0] = Prot.Substring(index, IonsUsedDigFilter + numAA - 2);
                                            similarSequenceLengths[1] = Prot.Substring(index, IonsUsedDigFilter + numAA - 1);
                                            similarSequenceLengths[2] = Prot.Substring(index, IonsUsedDigFilter + numAA);
                                            if (Prot.Length > IonsUsedDigFilter + numAA+1+index)
                                            {
                                                similarSequenceLengths[3] = Prot.Substring(index, IonsUsedDigFilter + numAA + 1);
                                            }
                                            else { similarSequenceLengths[3] = similarSequenceLengths[2]; }

                                            if (BPTMMass > -0.5 && BPTMMass < 0.5)
                                            {
                                                PTMInclude = 2;
                                            }
                                        }
                                        else //if y main
                                        {
                                            bigFrag = YPep.Substring(YPep.Count() - IonsUsedDigFilter, IonsUsedDigFilter);
                                            int fragIndex = Prot.IndexOf(bigFrag);
                                            while (Prot.Length - fragIndex - numAA - IonsUsedDigFilter > 0 && overExpMass == false && fragIndex - numAA > 0)
                                            {
                                                ProtProd = Prot.Substring(fragIndex - numAA, IonsUsedDigFilter + numAA); //what if it appears multiple times in a protein? FIXME http://stackoverflow.com/questions/6865419/indexof-for-multiple-results
                                                ProdMass = MonoIsoptopicMass(ProtProd);                                                               //additionProd = ProtProd.Substring(fragIndex+IonsUsedDigFilter,numAAUsed); //causing crash
                                                if (PTMInclude == 1)
                                                {
                                                    ProdMass += YPTMMass;
                                                }
                                                if (ProdMass > ExpMass)
                                                {
                                                    overExpMass = true;
                                                }
                                                numAA++;
                                            }
                                            similarSequenceLengths[0] = Prot.Substring(fragIndex - numAA, IonsUsedDigFilter + numAA - 2);
                                            similarSequenceLengths[1] = Prot.Substring(fragIndex - numAA, IonsUsedDigFilter + numAA - 1);
                                            similarSequenceLengths[2] = Prot.Substring(fragIndex - numAA, IonsUsedDigFilter + numAA);
                                            if (Prot.Length > IonsUsedDigFilter + numAA + 1 + index)
                                            {
                                                similarSequenceLengths[3] = Prot.Substring(fragIndex - numAA, IonsUsedDigFilter + numAA + 1);
                                            }
                                            else { similarSequenceLengths[3] = similarSequenceLengths[2]; }
                                            if (BPTMMass > -0.5 && BPTMMass < 0.5)
                                            {
                                                PTMInclude = 2;
                                            }
                                        }
                                        int nextAA = 0;
                                        foreach (string seq in similarSequenceLengths) //foreach of the similar sequences
                                        {
                                            while (Done == false && nextAA < seq.Length) //foreach of the chars in said sequence
                                            {
                                                foreach (char c in AANames) //list of 20 common aa
                                                {
                                                    double snpDeltaMass = MonoIsoptopicMass(c.ToString()) - MonoIsoptopicMass(seq[nextAA].ToString());
                                                    if ((ProdMass + snpDeltaMass) > (ExpMass * (1 - PrecursorMassTolerancePpm / 1000000)) && (ProdMass + snpDeltaMass) < (ExpMass * (1 + PrecursorMassTolerancePpm / 1000000))) //if match
                                                    {
                                                        DigestProduct[Ion] = seq + " + " + seq.Substring(nextAA, 1) + "->" + c;
                                                        Done = true;
                                                        PTMInclude = 2;
                                                    }
                                                }
                                                nextAA++;
                                            }
                                        }
                                    }
                                }
                            }
                            rowNum++;
                        }*/
           /*             Done = false;
                    }

                    //start parent search
                    if (DigestProduct[0] == "Z" && DigestProduct[1] == "Z") //if it's still a candidate
                    {
                        for (int Ion = 0; Ion<BandYMains.Count(); Ion++) //Ion specifies frag, where 0 is B and 1 is Y for the majority. This is used so that ABCD_EFGHIJ and ABCDEFGH_IJ doesn't search fasta for -IJ, which saves time
                        {
                            //Find other parent masses
                            string[] bothFrag = BandYMains[Ion].Split('_').ToArray();
string bigFrag = bothFrag[Ion]; //grab B frag if looking at main B fusion, or y if main y
string lilFrag = bothFrag[1 - Ion]; //1-1=0, 1-0=1. Flips the two.
bigFrag = Regex.Replace(bigFrag, "(\\(.*?\\))", "");
                            lilFrag = Regex.Replace(lilFrag, "(\\(.*?\\))", "");
                            string lilOrigID = ""; //what access num is associated
string bigOrigID = "";
                            if (Ion == 1) //want B ion (short) when Y is long, since it has greater chance of a parent match being found.
                            {
                                lilOrigID = DigestCandidates.Rows[i][1].ToString();
bigOrigID = DigestCandidates.Rows[i][3].ToString();
                            }
                            else
                            {
                                bigOrigID = DigestCandidates.Rows[i][1].ToString();
lilOrigID = DigestCandidates.Rows[i][3].ToString();
                            }
                            int counter = 0;


//Finds additional parent proteins that could yield the short fragment. Limits local FDR...
Boolean foundMatch = false;
                            while (counter<FASTA.Rows.Count && foundMatch == false)
                            {
                                if (FASTA.Rows[counter][1].ToString().Contains(lilFrag) == true && FASTA.Rows[counter][0].ToString() != lilOrigID)
                                {
                                    if (Ion == 1)
                                    {
                                        DigestCandidates.Rows[i][1] = DigestCandidates.Rows[i][1].ToString() + " OR Additional Targets. Ex: " + FASTA.Rows[counter][0].ToString();
                                    }
                                    else
                                    {
                                        DigestCandidates.Rows[i][3] = DigestCandidates.Rows[i][3].ToString() + " OR Additional Targets. Ex: " + FASTA.Rows[counter][0].ToString();
                                    }
                                    foundMatch = true; //end the loop or we'll be here for hours. Just need to know if a single target matches
                                }
                                counter++;
                            }

                            //3/18/17
                            //Determine if Cis or Trans Fusion Peptide
                            if (bigFrag.Length >= IonsUsedDigFilter)
                            {
                                foreach (DataRow row in FASTA.Rows)
                                {
                                    if (row[1].ToString().Contains(bigFrag)) //&& row[1].ToString().Contains(lilFrag)) //if protein from other half
                                    {
                                        string Prot = row[1].ToString();
                                        for (int c = 0; c<Prot.Count(); c++) //look for viable mass that could fulfill the other half
                                        {
                                            int fragLength = 1;
double targetMass = MonoIsoptopicMass(lilFrag);
string frag = Prot.Substring(c, fragLength);
double mass = MonoIsoptopicMass(frag);
                                            while (mass <= targetMass - 1 && c + fragLength<Prot.Count())
                                            {
                                                fragLength++;
                                                frag = Prot.Substring(c, fragLength);
                                                mass = MonoIsoptopicMass(frag);
                                            }
                                            if ((mass) > (targetMass* (1 - PrecursorMassTolerancePpm / 1000000)) && (mass) < (targetMass* (1 + PrecursorMassTolerancePpm / 1000000)))
                                            {
                                                DigestProduct[2] = "C";
                                                c = Prot.Count(); //end search
                                                Done = true;
                                            }
                                        }
                                    }
                                }
                            }
                            //if a replacement was not found within the same protein
                            if (DigestProduct[2] != "C")
                            {
                                DigestProduct[2] = "T";
                            }
                        }

                        //if interested, here are the peptides that are likely digestion products without any PTMs or modifications
                        DigestCandidates.Rows[i][6] = DigestProduct[0];
                        DigestCandidates.Rows[i][7] = DigestProduct[1];
                        DigestCandidates.Rows[i][8] = DigestProduct[2];
                        DigestProduct[0] = "Z";
                        DigestProduct[1] = "Z";
                        DigestProduct[2] = "Z";
                    }
                }
                else //if not a candidate, just fill with 'Z's
                {
                    //defined earlier: string[] DigestProduct = { "Z", "Z", "Z" };
                    DigestCandidates.Rows[i][6] = DigestProduct[0];
                    DigestCandidates.Rows[i][7] = DigestProduct[1];
                    DigestCandidates.Rows[i][8] = DigestProduct[2];
                }
            }
        }*/
    }
}

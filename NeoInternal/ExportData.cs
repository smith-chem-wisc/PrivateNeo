using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace NeoInternal
{
    public class ExportData
    {
        public static string folder;
        private readonly BackgroundWorker worker;

        public ExportData(BackgroundWorker worker)
        {
            this.worker = worker;
        }

        public string ExportAll(List<PSM> psms, string databaseFileName)
        {
            folder = DateTime.Now.ToString("yyyy-MM-dd_hh-mm-ss");
            string path = "";
            string[] temp=databaseFileName.Split('\\').ToArray();
            for(int i=0; i<temp.Count()-1; i++)
            {
                path += temp[i]+'\\';
            }
            Directory.CreateDirectory(path + folder);
            ExportCandidates(psms,path, out string e);
            if (e.Length > 0)
                return e;
            ExportFullFASTA(psms, databaseFileName, path);
            ExportFASTAAppendix(psms, databaseFileName, path);
            ExportFalsePositiveAppendix(psms, databaseFileName, path, out string e2);
            if (e2.Length > 0)
                return e2;
            ExportFilteredFusionPeptideAppendix(psms, databaseFileName, path);
            return "";
        }

        public void ExportCandidates(List<PSM> psms, string path, out string error_message)
        {
            string mutableError_message = "";
            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "ExportedFusionCandidates.txt"))
            {
                file.WriteLine("Scan" + '\t' + "ExperimentalMass" + '\t' + "OriginalNSequence" + '\t' + "OriginalNScore" + '\t' + "OriginalCSequence" + '\t' + "OriginalCScore" + '\t' + "SampleSequence" + '\t' + "Ambiguity" + '\t' + "ProbableType" + '\t' + "MostProbableSequenceJunctions" + '\t' + "MostProbableSequence(s)" + '\t' + "MostProbableParents" + '\t' + "AllPossibleSequenceJunctions" + '\t' + "AllPossibleSequence(s)" + '\t' + "AllPossibleParent(s)" + '\t' + "NumberOfPossibleSequences" + "PotentialFalsePositives");

                int progress = 0;
                Parallel.ForEach(psms, (psm) =>
                {
                   //printout the scan, the mass, the sequences with and without junctions, the number of potential sequences
                   string allPossibleSequences = "";
                    string mostProbableSequences = "";
                    int indexOfFirstProbableSequence = -1;
                    string mostProbableParents = "";
                    string allPossibleParents = "";
                    FusionCandidate.FusionType probableType = psm.fusionType;
                    for (int fc = 0; fc < psm.getFusionCandidates().Count(); fc++)
                    {
                        FusionCandidate fusionCandidate = psm.getFusionCandidates()[fc]; //need fc for indexOfFirstProbableSequence
                       char[] tempArray = fusionCandidate.seq.ToCharArray();
                        //if most probable, add to all and most
                        if (fusionCandidate.fusionType.Equals(probableType))
                        {
                            //record sequences
                            if (indexOfFirstProbableSequence < 0)
                            {
                                indexOfFirstProbableSequence = fc;
                            }
                            for (int i = 0; i < tempArray.Count(); i++)
                            {
                                mostProbableSequences += tempArray[i];
                                allPossibleSequences += tempArray[i];
                                foreach (int junction in fusionCandidate.getJunctionIndexes())
                                {
                                    if (junction == i)
                                    {
                                        mostProbableSequences += "-";
                                        allPossibleSequences += "-";
                                    }
                                }
                            }
                            mostProbableSequences += "|";
                            allPossibleSequences += "|";

                            //record parents
                            string tempParents = "";
                            switch (probableType)
                            {
                                case FusionCandidate.FusionType.TL:
                                    tempParents += GenerateParentOutput(fusionCandidate.translatedParents, new List<CisParent>(), new List<TransParent>());
                                    break;

                                case FusionCandidate.FusionType.NC:
                                case FusionCandidate.FusionType.RC:
                                    tempParents += GenerateParentOutput(new List<TranslatedParent>(), fusionCandidate.cisParents, new List<TransParent>());
                                    break;

                                default: //if trans
                                    tempParents += GenerateParentOutput(new List<TranslatedParent>(), new List<CisParent>(), fusionCandidate.transParents);
                                    break;
                            }
                            mostProbableParents += tempParents;
                            allPossibleParents += tempParents;
                        }
                        else //not most probable, only add it to allPossibleSequences
                        {
                            //record sequences
                            for (int i = 0; i < tempArray.Count(); i++)
                            {
                                allPossibleSequences += tempArray[i];
                                if (fusionCandidate.getJunctionIndexes().Contains(i))
                                {
                                    allPossibleSequences += "-";
                                }
                            }
                            allPossibleSequences += "|";
                            //record parents
                            allPossibleParents += GenerateParentOutput(fusionCandidate.translatedParents, fusionCandidate.cisParents, fusionCandidate.transParents);
                        }
                       /*        foreach(ParentInfo PI in fusionCandidate.parentInfo)
                               {
                                   parents += PI.accession + "_" + PI.parentType.ToString() + "_" + PI.seqFound + "|";
                               }*/
                    }

                    allPossibleSequences = allPossibleSequences.Substring(0, allPossibleSequences.Length - 1); //remove last "|"
                   mostProbableSequences = mostProbableSequences.Substring(0, mostProbableSequences.Length - 1); //remove last "|"

                   string ambiguity = "";
                    AlternativeSequences.findIons(psm.getFusionCandidates()[indexOfFirstProbableSequence], psm, out string e); //this should be carried over, but it's not...
                   lock (mutableError_message)
                    {
                        mutableError_message += e;
                    }
                    bool[] foundIons = psm.getFusionCandidates()[indexOfFirstProbableSequence].getFoundIons();
                    char[] firstSeq = psm.getFusionCandidates()[indexOfFirstProbableSequence].seq.ToCharArray();
                   //   if(foundIons.Count()==firstSeq.Count()) //prevent crashing if something went wrong
                   // {
                   bool ambiguous = false;
                    for (int i = 0; i < foundIons.Count(); i++)
                    {
                        if (foundIons[i]) //if found
                       {
                            ambiguity += firstSeq[i]; //add aa
                           if (ambiguous) //if it is part of an ambiguous sequence
                           {
                                ambiguity += ")";
                                ambiguous = false; //no longer ambiguous
                           }
                        }
                        else
                        {
                            if (!ambiguous)
                            {
                                ambiguous = true;
                                ambiguity += "(";
                            }
                            ambiguity += firstSeq[i];
                        }
                    }
                    string potentialFalsePositives = "";
                    foreach (Variant v in psm.variants)
                    {
                        potentialFalsePositives += v.id + "_" + v.start + "-" + (v.start + v.peptideLength - 1) + "(" + v.pepSeq + ")" + v.varType + "|";
                    }
                    if (potentialFalsePositives.Length > 0)
                    {
                        potentialFalsePositives = potentialFalsePositives.Substring(0, potentialFalsePositives.Length - 1); //remove last |
                   }
                   //workarounds for excel. Actual limit is 32767, but that doesn't seem to work
                   if (mostProbableParents.Length > 30000)
                    {
                        mostProbableParents = mostProbableParents.Substring(0, 30000);
                    }
                    if (allPossibleParents.Length > 30000)
                    {
                        allPossibleParents = allPossibleParents.Substring(0, 30000);
                    }
                   //  }
                   //file.WriteLine("Scan" +                     '\t' + "ExperimentalMass" +         '\t' + "OriginalBSequence" + '\t' + "OriginalBScore" +      '\t' + "OriginalYSequence" +    '\t' + "OriginalYScore" +       '\t' + "SampleSequence" +                                        '\t' + "Ambiguity" +       '\t' + "ProbableType" +     '\t' + "MostProbablySequenceJunctions"+ '\t' + "MostProbableSequence(s)" +            '\t' + "MostProbableParents" + '\t' + "AllPossibleSequenceJunctions" + '\t' + "AllPossibleSequence(s)" + '\t' + "AllPossibleParent(s)" + '\t' + "NumberOfPossibleSequences");
                   lock (file)
                    {
                        file.WriteLine(psm.getScan().ToString() + '\t' + psm.getExpMass().ToString() + '\t' + psm.getNInfo().seq + '\t' + psm.getNInfo().score + '\t' + psm.getCInfo().seq + '\t' + psm.getCInfo().score + '\t' + psm.getFusionCandidates()[indexOfFirstProbableSequence].seq + '\t' + ambiguity + '\t' + psm.fusionType.ToString() + '\t' + mostProbableSequences + '\t' + mostProbableSequences.Replace("-", "") + '\t' + mostProbableParents + '\t' + allPossibleSequences + '\t' + allPossibleSequences.Replace("-", "") + '\t' + allPossibleParents + '\t' + psm.getFusionCandidates().Count().ToString() + '\t' + potentialFalsePositives);
                        progress++;
                        this.worker.ReportProgress(progress / psms.Count()* 100);
                    }
                });
            }
            /*      using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\Neo\Results\" + folder + @"\"+folder+"testFlipper.txt"))
                  {
                      foreach (PSM psm in psms)
                      {
                          //printout the scan, the mass, the sequences with and without junctions, the number of potential sequences
                          foreach (FusionCandidate fusionCandidate in psm.getFusionCandidates())
                          {
                              file.WriteLine(fusionCandidate.seq);
                          }
                      }
                  }*/
            error_message = mutableError_message;
        }

        private static void ExportFullFASTA(List<PSM> psms, string databaseFileName, string path)
        {
            //@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\Neo\Results\"
            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "FullFusionDatabase.fasta"))
            {
                //copy database
                string[] FASTARead = File.ReadAllLines(databaseFileName);
                foreach (string s in FASTARead)
                {
                    file.WriteLine(s);
                }

                //Write FalsePositive Sources
                foreach (PSM psm in psms)
                {
                    int fusNum = 1;

                    //fusionPeptides = Regex.Replace(candidateRow[5].ToString(), "(\\(.*?\\))", "");//remove PTM annotations

                    //the following code is for FASTA output
                    string scan = psm.getScan().ToString();
                    foreach (Variant v in psm.variants)
                    {
                        file.WriteLine(">sp|" + scan + v.varType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + v.varType.ToString() + " peptide GN=FUS PE=1 SV=1");
                        string seq = v.pepSeq;
                        for (int i = 0; i < seq.Count(); i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Count())
                            {
                                file.WriteLine(seq.Substring(i, 60));
                            }
                            else
                            {
                                file.WriteLine(seq.Substring(i, seq.Length - i));
                            }
                        }
                        fusNum++;
                    }
                }

                //Write Fusion Peptide Candidates
                foreach (PSM psm in psms)
                {
                    int fusNum = 1;

                    //fusionPeptides = Regex.Replace(candidateRow[5].ToString(), "(\\(.*?\\))", "");//remove PTM annotations

                    //the following code is for FASTA output
                    string scan = psm.getScan().ToString();
                    foreach (FusionCandidate fusionCandidate in psm.getFusionCandidates())
                    {
                        file.WriteLine(">sp|" + scan + fusionCandidate.fusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.fusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                        string seq = fusionCandidate.seq;
                        for (int i = 0; i < seq.Count(); i += 60) //60 used as number of AAs per line in a FASTA file. It is unlikely to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Count())
                            {
                                file.WriteLine(seq.Substring(i, 60));
                            }
                            else
                            {
                                file.WriteLine(seq.Substring(i, seq.Length - i));
                            }
                        }
                        fusNum++;
                    }
                }
            }
        }

        private static void ExportFASTAAppendix(List<PSM> psms, string databaseFileName, string path)
        {
            //@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\Neo\Results\"
            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "FusionDatabaseAppendix.fasta"))
            {
                foreach (PSM psm in psms)
                {
                    int fusNum = 1;

                    //fusionPeptides = Regex.Replace(candidateRow[5].ToString(), "(\\(.*?\\))", "");//remove PTM annotations

                    //the following code is for FASTA output
                    string scan = psm.getScan().ToString();
                    foreach (FusionCandidate fusionCandidate in psm.getFusionCandidates())
                    {
                        file.WriteLine(">sp|" + scan + fusionCandidate.fusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.fusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                        string seq = fusionCandidate.seq;
                        for (int i = 0; i < seq.Count(); i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Count())
                            {
                                file.WriteLine(seq.Substring(i, 60));
                            }
                            else
                            {
                                file.WriteLine(seq.Substring(i, seq.Length - i));
                            }
                        }
                        fusNum++;
                    }
                }
            }
        }
        private static void ExportFalsePositiveAppendix(List<PSM> psms, string databaseFileName, string path, out string error_message)
        {
            error_message = "";
            //@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\Neo\Results\"
            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "FalsePositiveAppendix.fasta"))
            {
                foreach (PSM psm in psms)
                {
                    int fusNum = 1;

                    //fusionPeptides = Regex.Replace(candidateRow[5].ToString(), "(\\(.*?\\))", "");//remove PTM annotations

                    //the following code is for FASTA output
                    string scan = psm.getScan().ToString();
                    foreach (Variant v in psm.variants)
                    {
                        file.WriteLine(">sp|" + scan + v.varType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + v.varType.ToString() + " peptide GN=FUS PE=1 SV=1");
                        string seq = v.pepSeq;
                        for (int i = 0; i < seq.Count(); i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Count())
                            {
                                file.WriteLine(seq.Substring(i, 60));
                            }
                            else
                            {
                                file.WriteLine(seq.Substring(i, seq.Length - i));
                            }
                        }
                        fusNum++;
                    }
                }
            }

            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "FalsePositiveAppendix.xml"))
            {
                using (StreamReader header = new StreamReader(Path.Combine(Environment.CurrentDirectory, "xmlHeader.txt"))) //file located in Morpheus folder
                {
                    while (header.Peek() != -1)
                    {
                        string line = header.ReadLine();
                        file.WriteLine(line);
                    }
                }
                foreach (PSM psm in psms)
                {
                    int fusNum = 1;

                    //fusionPeptides = Regex.Replace(candidateRow[5].ToString(), "(\\(.*?\\))", "");//remove PTM annotations

                    //the following code is for XML output
                    string scan = psm.getScan().ToString();
                    foreach (Variant v in psm.variants)
                    {
                        file.WriteLine("<entry>");
                        file.WriteLine("<accession>" + scan + v.varType.ToString() + fusNum + "</accession>");
                        file.WriteLine("<name>" + scan + v.varType.ToString() + fusNum + "</name>");
                        file.WriteLine("<protein>");
                        file.WriteLine("<recommendedName>");
                        file.WriteLine("<fullName>" + scan + v.varType.ToString() + fusNum + "</fullName>");
                        file.WriteLine("</recommendedName>");
                        file.WriteLine("</protein>");
                        file.WriteLine("<gene>");
                        file.WriteLine("<name type=" + '"' + "primary" + '"' + ">Fus</name>");
                        file.WriteLine("</gene>");
                        file.WriteLine("<dbReference type=" + '"' + "Neo" + '"' + " id=" + '"' + scan + fusNum + '"' + ">");
                        file.WriteLine("<property type=" + '"' + "entry name" + '"' + " value=" + '"' + scan + fusNum + '"' + "/>");
                        file.WriteLine("</dbReference>");
                        file.WriteLine("<feature type=" + '"' + "chain" + '"' + ">");
                        file.WriteLine("<location>");
                        file.WriteLine("<begin position=" + '"' + "1" + '"' + " />");
                        string[] fullSequence = v.pepSeq.Split('+');
                        file.WriteLine("<end position=" + '"' + fullSequence[0].Count() + '"' + " />");
                        file.WriteLine("</location>");
                        file.WriteLine("</feature>");
                        //EnterPTMInfo
                        if (v.varType.Equals(Variant.variantType.PTM) && fullSequence.Count() > 1)
                        {
                            string ptmName = fullSequence[1];
                            char[] fullSeqArray = fullSequence[0].ToCharArray();
                            for (int position = 0; position < fullSeqArray.Count(); position++)
                            {
                                switch (ptmName)
                                {
                                    case "Deamidation":
                                        if (fullSeqArray[position].Equals('N') | fullSeqArray[position].Equals('Q'))
                                        {
                                            printPTMInfo(ptmName + " on " + fullSeqArray[position], position + 1, file);
                                        }
                                        break;
                                    case "Trioxidation":
                                        if (fullSeqArray[position].Equals('C'))
                                        {
                                            printPTMInfo(ptmName + " of " + fullSeqArray[position], position + 1, file);
                                        }
                                        goto case "Dioxidation";
                                    case "Dioxidation":
                                        if (fullSeqArray[position].Equals('P') | fullSeqArray[position].Equals('W'))
                                        {
                                            printPTMInfo("Dioxidation" + " on " + fullSeqArray[position], position + 1, file);
                                        }
                                        goto case "Oxidation";
                                    case "Oxidation":
                                        if (fullSeqArray[position].Equals('D') | fullSeqArray[position].Equals('E') | fullSeqArray[position].Equals('F') | fullSeqArray[position].Equals('H') | fullSeqArray[position].Equals('P') | fullSeqArray[position].Equals('V') | fullSeqArray[position].Equals('W') | fullSeqArray[position].Equals('Y'))
                                        {
                                            printPTMInfo("Oxidation" + " on " + fullSeqArray[position], position + 1, file);
                                        }
                                        break;
                                    case "Acetyl":
                                        if (fullSeqArray[position].Equals('K'))
                                        {
                                            printPTMInfo(ptmName + "lysine", position + 1, file);
                                        }
                                        else if (fullSeqArray[position].Equals('S'))
                                        {
                                            printPTMInfo(ptmName + "serine", position + 1, file);
                                        }
                                        else if (position == 1)
                                        {
                                            printPTMInfo(ptmName, position + 1, file);
                                        }
                                        break;
                                    case "Phospho":
                                        if (fullSeqArray[position].Equals('S') | fullSeqArray[position].Equals('T') | fullSeqArray[position].Equals('Y') | fullSeqArray[position].Equals('H') | fullSeqArray[position].Equals('R') | fullSeqArray[position].Equals('K'))
                                        {
                                            printPTMInfo(ptmName + " on " + fullSeqArray[position], position + 1, file);
                                        }
                                        break;
                                    case "Zinc":
                                        if (fullSeqArray[position].Equals('D') | fullSeqArray[position].Equals('E'))
                                        {
                                            printPTMInfo(ptmName + " on " + fullSeqArray[position], position + 1, file);
                                        }
                                        break;
                                    case "Potassium":
                                        if (fullSeqArray[position].Equals('D') | fullSeqArray[position].Equals('E'))
                                        {
                                            printPTMInfo(ptmName + " on " + fullSeqArray[position], position + 1, file);
                                        }
                                        break;
                                    case "Sodium":
                                        if (fullSeqArray[position].Equals('D') | fullSeqArray[position].Equals('E'))
                                        {
                                            printPTMInfo(ptmName + " on " + fullSeqArray[position], position + 1, file);
                                        }
                                        break;
                                    case "Carbamidomethyl":
                                        if (fullSeqArray[position].Equals('K'))
                                        {
                                            printPTMInfo(ptmName, position + 1, file);
                                        }
                                        break;
                                    case "TriMethyl":
                                        if (fullSeqArray[position].Equals('K'))
                                        {
                                            printPTMInfo(ptmName + "ation", position + 1, file);
                                        }
                                        goto case "DiMethyl";
                                    case "DiMethyl":
                                        if (fullSeqArray[position].Equals('A') | fullSeqArray[position].Equals('K') | fullSeqArray[position].Equals('N') | fullSeqArray[position].Equals('R'))
                                        {
                                            printPTMInfo("DiMethyl" + " on " + fullSeqArray[position], position + 1, file);
                                        }
                                        goto case "Methyl";
                                    case "Methyl":
                                        if (fullSeqArray[position].Equals('D') | fullSeqArray[position].Equals('E') | fullSeqArray[position].Equals('H') | fullSeqArray[position].Equals('K') | fullSeqArray[position].Equals('N') | fullSeqArray[position].Equals('Q') | fullSeqArray[position].Equals('R') | fullSeqArray[position].Equals('S') | fullSeqArray[position].Equals('V'))
                                        {
                                            printPTMInfo("Methyl" + " on " + fullSeqArray[position], position + 1, file);
                                        }
                                        break;
                                    case "Ammonia loss":
                                        printPTMInfo(ptmName, position, file);
                                        break;
                                    case "Sulfonation":
                                        if (fullSeqArray[position].Equals('S') | fullSeqArray[position].Equals('T') | fullSeqArray[position].Equals('Y'))
                                        {
                                            printPTMInfo(ptmName + " on " + fullSeqArray[position], position + 1, file);
                                        }
                                        break;
                                    default:
                                        error_message += "PTM " + ptmName + " was not found";
                                        break;
                                }
                            }
                        }
                        file.WriteLine("<sequence length=" + '"' + fullSequence[0].Length + '"' + ">" + fullSequence[0] + "</sequence>");
                        file.WriteLine("</entry>");
                        fusNum++;
                    }
                }
                file.WriteLine("</mzLibProteinDb>");
            }
        }

        private static void printPTMInfo(string ptmName, int position, StreamWriter file)
        {
            file.WriteLine("<feature type=" + '"' + "modified residue" + '"' + " description=" + '"' + ptmName + '"' + ">");
            file.WriteLine("<location>");
            file.WriteLine("<position position=" + '"' + position + '"' + " />");
            file.WriteLine("</location>");
            file.WriteLine("</feature>");
        }

        private static void ExportFilteredFusionPeptideAppendix(List<PSM> psms, string databaseFileName, string path)
        {
            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "FilteredFusionDatabaseAppendix.fasta"))
            {
                foreach (PSM psm in psms)
                {
                    int fusNum = 1;
                    if (!psm.fusionType.Equals(FusionCandidate.FusionType.TL) && psm.variants.Count()==0)
                    {
                        //fusionPeptides = Regex.Replace(candidateRow[5].ToString(), "(\\(.*?\\))", "");//remove PTM annotations

                        //the following code is for FASTA output
                        string scan = psm.getScan().ToString();
                        foreach (FusionCandidate fusionCandidate in psm.getFusionCandidates())
                        {
                            file.WriteLine(">sp|" + scan + fusionCandidate.fusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.fusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                            string seq = fusionCandidate.seq;
                            for (int i = 0; i < seq.Count(); i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                            {
                                if ((i + 60) < seq.Count())
                                {
                                    file.WriteLine(seq.Substring(i, 60));
                                }
                                else
                                {
                                    file.WriteLine(seq.Substring(i, seq.Length - i));
                                }
                            }
                            fusNum++;
                        }
                    }
                }
            }
        }

        private static string GenerateParentOutput(List<TranslatedParent> translatedParents, List<CisParent> cisParents, List<TransParent> transParents)
        {
            string output = "";
            foreach (TranslatedParent tlp in translatedParents)
            {
                output += tlp.id + "_" + tlp.start + "-" + (tlp.start + tlp.peptideLength-1) + "(" + tlp.seq.Substring(tlp.start, tlp.peptideLength) + ")" + "|";
            }
            foreach (CisParent cp in cisParents)
            {
                foreach (int ns in cp.nStart)
                {
                    foreach (int cs in cp.cStart)
                    {
                        output += cp.id + "_" + ns + "-" + (ns + cp.nLength-1) + "(" + cp.seq.Substring(ns, cp.nLength) + ")"
                            + "&" + cs + "-" + (cs + cp.cLength-1) + "(" + cp.seq.Substring(cs, cp.cLength) + ")" + "|";
                    }
                }
            }
            string nOutput = "";
            string cOutput = "";
            foreach (TransParent tp in transParents)
            {
                if(tp.terminal.Equals(ParentInfo.terminal.N))
                {
                    foreach (int ts in tp.start)
                    {
                        nOutput += tp.id + "_" + ts + "-" + (ts + tp.peptideLength-1) + "(" + tp.seq.Substring(ts, tp.peptideLength) + ")" + "|";
                    }
                }
                else
                {
                    foreach (int ts in tp.start)
                    {
                        cOutput += tp.id + "_" + ts + "-" + (ts + tp.peptideLength-1) + "(" + tp.seq.Substring(ts, tp.peptideLength) + ")" + "|";
                    }
                }
            }
            if (nOutput.Length > 0) //if there was a trans, append and delimit the two lists
            {
                output += nOutput + "&&|" + cOutput;
            }
            output = output.Substring(0, output.Length - 1); //remove last |

            return output;
        }
    }
}

﻿using System;
using System.Data;
using System.Globalization;
using System.IO;
using System.Linq;

namespace NeoInternal
{
    public class MassCalculator
    {
        public static double[] MONOISOTOPIC_AMINO_ACID_MASSES;
        public static double[] AVERAGE_AMINO_ACID_MASSES;
        public static DataTable ModificationsDT = new DataTable();

        public static void importMasses()
        {
            AminoAcidMasses();
            ModificationMasses();
        }

        private static void AminoAcidMasses()
        {
            MONOISOTOPIC_AMINO_ACID_MASSES = new double['Z' - 'A' + 1]; //makes array with 26 0's
            AVERAGE_AMINO_ACID_MASSES = new double['Z' - 'A' + 1];
            for (int i = 0; i < MONOISOTOPIC_AMINO_ACID_MASSES.Length; i++)
            {
                MONOISOTOPIC_AMINO_ACID_MASSES[i] = double.NaN;
                AVERAGE_AMINO_ACID_MASSES[i] = double.NaN;
            }

            using (StreamReader amino_acids = new StreamReader(Path.Combine(Environment.CurrentDirectory, "amino_acids.tsv"))) //file located in Morpheus folder
            {
                amino_acids.ReadLine();

                while (amino_acids.Peek() != -1)
                {
                    string line = amino_acids.ReadLine();
                    string[] fields = line.Split('\t');

                    char one_letter_code = char.Parse(fields[0]);
                    if (!char.IsUpper(one_letter_code))
                    {
                        throw new ArgumentOutOfRangeException("Invalid amino acid abbreviation: " + one_letter_code);
                    }
                    double monoisotopic_mass = double.Parse(fields[1], CultureInfo.InvariantCulture);
                    MONOISOTOPIC_AMINO_ACID_MASSES[one_letter_code - 'A'] = monoisotopic_mass;
                    double average_mass = double.Parse(fields[2], CultureInfo.InvariantCulture);
                    AVERAGE_AMINO_ACID_MASSES[one_letter_code - 'A'] = average_mass;
                }
            }
        }

        private static void ModificationMasses()
        {
            ModificationsDT.Columns.Add("Name", typeof(string));
            ModificationsDT.Columns.Add("MonoisotopicMass", typeof(double));
            using (StreamReader uniprot_mods = new StreamReader(Path.Combine(Environment.CurrentDirectory, "Mods.txt")))///"ptmlist.txt")))
            {
                string description = null;
                string feature_type = null;
                double monoisotopic_mass_shift = double.NaN;
                double average_mass_shift = double.NaN;
                while (uniprot_mods.Peek() != -1)
                {
                    string line = uniprot_mods.ReadLine();
                    if (line.Length >= 2)
                    {
                        switch (line.Substring(0, 2))
                        {
                            case "ID":
                                description = line.Substring(5);
                                description = RemoveNestedParentheses(description, true);
                                description = description.Replace("'", "");
                                break;
                            case "FT":
                                feature_type = line.Substring(5);
                                break;
                            case "TG":
                                if (feature_type == "MOD_RES")
                                {
                                    string amino_acid = line.Substring(5);
                                }
                                break;
                            case "PP":
                                if (feature_type == "MOD_RES")
                                {
                                }
                                break;
                            case "MM":
                                monoisotopic_mass_shift = double.Parse(line.Substring(5));
                                break;
                            case "MA":
                                average_mass_shift = double.Parse(line.Substring(5));
                                break;
                            case "DR":
                                if (line.Contains("PSI-MOD"))
                                {
                                }
                                break;
                            case "//":
                                ModificationsDT.Rows.Add(description, monoisotopic_mass_shift);
                                break;
                        }
                    }
                }
            }
        }

        public static double MonoIsoptopicMass(string baseSequence, out string error_message)
        {
            //troubleshooter++;
            //if(Convert.ToDouble(troubleshooter)==Math.Round(Convert.ToDouble(troubleshooter)/1000)*1000)
            //{
            //MessageBox.Show(troubleshooter.ToString() + " " + baseSequence);
            //}
            error_message = "";
            double monoisotopic_mass = Constants.WATER_MONOISOTOPIC_MASS;
            Boolean ModificationOn = false;
            string ModificationName = "";
            foreach (char amino_acid in baseSequence)
            {
                if (amino_acid == ')') //only occurs at end of mod
                {
                    ModificationOn = false;
                    //if (ModificationName == "oxidation of M") //annotated differently than uniprot
                    //if (ModificationName == "v:Oxidation of M anywhere") //annotated differently than uniprot
                    if (ModificationName == "v:Oxidation") //Added 3/6/17
                    {
                        monoisotopic_mass += 15.99491463;
                        ModificationName = "";
                    }
                    else
                    {
                        DataRow[] PTMRow = ModificationsDT.Select("Name = '" + ModificationName + "'");
                        try
                        {
                            monoisotopic_mass += Convert.ToDouble(PTMRow[0][1]);
                        }
                        catch
                        {
                            error_message += "PTM " + "'" + ModificationName + "'" + " could not be found";
                        }
                        ModificationName = "";
                    }
                }
                if (ModificationOn == true) //only occurs if "(" already found
                {
                    ModificationName += amino_acid;
                }
                if (amino_acid == '(') //start collecting PTM name
                {
                    ModificationOn = true;
                }
                if (ModificationOn == false && amino_acid != ')')
                {
                    monoisotopic_mass += GetMonoisotopicMass(amino_acid, baseSequence, out string error_message2); //something making it here after (
                    error_message += error_message2;
                }
            }
            return monoisotopic_mass;
        }

        public static double GetMonoisotopicMass(char aminoAcid, string seq, out string error_message)
        {
            error_message = "";
            try
            {
                return MONOISOTOPIC_AMINO_ACID_MASSES[aminoAcid - 'A'];
            }
            catch
            {
                error_message += "Invalid amino acid '" + aminoAcid + "' caught " + seq;
                return 123;
            }
        }

        public static double GetAverageMass(char aminoAcid)
        {
            return AVERAGE_AMINO_ACID_MASSES[aminoAcid - 'A'];
        }

        public static string RemoveNestedParentheses(string lr, bool parentheses)
        {
            //Revised 3/6/2017. MetaMorpheus now uses brackets instead of parentheses for ptms
            int nested = 0;
            char[] lrCharArray = lr.ToCharArray();
            for (int i = 0; i < lrCharArray.Count(); i++)
            {
                if (lrCharArray[i] == '[') { lrCharArray[i] = '('; } //Added 3/6/17
                else if (lrCharArray[i] == ']') { lrCharArray[i] = ')'; } //Added 3/6/17
                if (lrCharArray[i] == '(')
                {
                    if (parentheses == true)
                    {
                        lrCharArray[i] = '[';
                        nested++;
                    }
                    else
                    {
                        parentheses = true;
                    }
                }
                if (lrCharArray[i] == ')')
                {
                    if (nested > 0)
                    {
                        lrCharArray[i] = ']';
                        nested--;
                    }
                    else {
                        parentheses = false;
                    }
                }
            }
            string cleanedSequence = new string(lrCharArray);
            return cleanedSequence;
        }

        public static double getPTMMass(string ptmName, out string error_message)
        {
            error_message = "";
            string searchString = "Name ==" + ptmName;
            try
            {
                return Convert.ToDouble(ModificationsDT.Select(searchString)[0]);
            }
            catch
            {
                error_message += "Oops! We weren't able to find the PTM named " + ptmName + " in our Mods.txt file.";
                return 0;
            }
        }
    }
}

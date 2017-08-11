namespace NeoInternal
{
    public class TranslatedParent : TheoreticalProtein
    {
        public int start { get; set; }
        public int peptideLength { get; set; }
        public FusionCandidate.FusionType translated { get; }

        public TranslatedParent(string id, string seq, int start, int length):base(id, seq)
        {
            this.start = start;
            this.peptideLength = length;
        }
    }
}

namespace NeoInternal
{
    public class TheoreticalProtein
    {
        public string id { get; set; }
        public string seq { get; set; }

        public TheoreticalProtein(string id, string seq)
        {
            this.id = id;
            this.seq = seq;
        }
    }
}

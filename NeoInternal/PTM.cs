namespace NeoInternal
{
    public class PTM
    {
        public string name { get; set; }
        public int index { get; set; }
        public double mass { get; set; }

        public PTM(string name, int index, double mass)
        {
            this.name = name;
            this.index = index;
            this.mass = mass;
        }
    }
}

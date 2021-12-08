namespace Aligner3
{
    public class Alignment
    {
        public string Represented { get; }

        private int[] counts;
        
        public Alignment(string longest, int count)
        {
            Represented = longest;
            counts = new int[longest.Length];
            for (int i = 0; i < counts.Length; i++)
                counts[i] = count;
        }

        public void AddSequence(int firstPosition, int size, int count)
        {
            int maxsize = firstPosition + size;
            for (int i = firstPosition; i < maxsize; i++)
                counts[i] += count;
        }

        public int[] GetCounts()
        {
            return counts;
        }
    }
}
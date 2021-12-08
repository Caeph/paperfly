namespace SamplerEulerian
{
    public class SequenceSampler : BaseSampler
    {
        public SequenceSampler(int k, DeBruijnGraph deBruijnGraph) : base(deBruijnGraph, k)
        {
        }
    }
}
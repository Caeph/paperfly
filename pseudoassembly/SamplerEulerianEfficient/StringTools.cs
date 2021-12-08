using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SamplerEulerian
{
    public class StringTools
    {
        public static string DebruijnJoin(List<string> kmers, int k)
        {
            var builder = new StringBuilder(kmers.First());
            foreach (var kmer in kmers.Skip(1))
                builder.Append(kmer.Substring(k - 1));

            return builder.ToString();
        }
    }
}
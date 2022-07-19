using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SamplerEulerian
{
    public class StringTools
    {
        public static string DebruijnJoin(List<string> kmers, int k)
        {
            //assumption -- kmers overlap by k-1 chars
            var builder = new StringBuilder(kmers.First());
            foreach (var kmer in kmers.Skip(1))
                builder.Append(kmer.Substring(k - 1));

            return builder.ToString();
        }

        // public static string DebruijnJoinNoAssumption(List<string> kmers, int k)
        // {
        //     //no assumption on the overlap size
        //     var builder = new StringBuilder(kmers.First());
        //     foreach (var kmer in kmers.Skip(1))
        //     {
        //         string toadd;
        //         int largestOverlap = 0;
        //         // get the largest overlap
        //         for (int overlapsize = k-1; overlapsize > 0; overlapsize--)
        //         {
        //             bool fits = false;
        //             
        //             if (fits)
        //             {
        //                 largestOverlap = overlapsize;
        //                 break;
        //             }
        //         }
        //
        //         toadd = kmer.Substring(largestOverlap);
        //         builder.Append(toadd);
        //     }
        //
        //     return builder.ToString();
        // }
    }
}
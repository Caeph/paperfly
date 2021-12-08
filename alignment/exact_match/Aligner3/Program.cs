using System;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace Aligner3
{
    internal class Program
    {
        public static void Main(string[] args)
        {
            
            
            //ARGUMENT PROCESSING
            string csvDir = args[0];
            string outputFile = args[1];
            int k = int.Parse(args[2]);
            //float coef = float.Parse(args[3]);
            //float distanceThrPercentage = coef * 1f/(k - 1);

            var sw_overall = Stopwatch.StartNew();
            
            //do stuff ON EACH CSV FILE
            using (var writer = new StreamWriter(outputFile))
            {
                foreach (var file in Directory.EnumerateFiles(csvDir))
                {
                    var sw = Stopwatch.StartNew();
                    
                    var aligner = new AlignerEngine(file, k);
                    aligner.Align();

                    foreach (var alignment in aligner.GetAlignments())
                    {
                        var sq = alignment.Represented;
                        writer.Write(
                            sq
                        );
                        writer.Write(
                            ';'
                        );
                        writer.WriteLine(
                            string.Join(",",alignment.GetCounts().Select(x => x.ToString()))
                        );
                    }
                    
                    sw.Stop();
                    Console.WriteLine($"{file}={sw.Elapsed.ToString()}");
                }
            }

            sw_overall.Stop();
            Console.WriteLine($"Overall elapsed time = {sw_overall.Elapsed.ToString()}");
        }
    }
}
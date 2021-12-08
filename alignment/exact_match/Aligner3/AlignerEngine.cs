using System;
using System.Collections.Generic;
using System.Diagnostics.Eventing.Reader;
using System.IO;
using System.Linq;
using System.Xml.XPath;
using Aligner;

namespace Aligner3
{
    public class AlignerEngine
    {
        private int k;
        private string csvFile;

        private List<Alignment> against = new List<Alignment>();
        
        public AlignerEngine(string csvFile, int k)
        {
            this.csvFile = csvFile;
            this.k = k;
        }

        private (string, int) readLine(string line)
        {
            var i = line.Split(';');
            return (i[0], int.Parse(i[1]));
        }

        public List<Alignment> GetAlignments()
        {
            return against;
        }
        public void Align()
        {
            var assembled = File.ReadLines(csvFile).Select(s => readLine(s)).GroupBy(tup => tup.Item1.Length).OrderByDescending(g => g.Key);
            //find longest, set as targets

            var enumerator = assembled.GetEnumerator();
            if (!enumerator.MoveNext()) return;
            
            var longest = enumerator.Current;

            foreach (var item in longest)
            {
                against.Add(
                    new Alignment(item.Item1, item.Item2)
                    );
            }

            while (enumerator.MoveNext())
            {
                var sequenceCounts = enumerator.Current.GroupBy(x => x.Item1).ToDictionary(gr => gr.Key, gr => gr.Sum(x => x.Item2));
                var tree = new AhoCorasickTree(sequenceCounts.Keys);
                //align 

                var matchesInGroup = new Dictionary<string, List<(Alignment, int)>>();
                foreach (var target in against)
                {
                    var matches = tree.FindAll(target.Represented);
                    
                    foreach (var match in matches)
                    {
                        if (matchesInGroup.ContainsKey(match.Item1))
                        {
                            matchesInGroup[match.Item1].Add((target,match.Item2));
                        }
                        else
                        {
                            matchesInGroup.Add(match.Item1, new List<(Alignment, int)>() { (target,match.Item2) }); //alignment, offset
                        }
                    }
                }

                var aligned = matchesInGroup.Keys;
                foreach (var sq in aligned)
                {
                    var to = matchesInGroup[sq].OrderByDescending(x => x.Item1.Represented.Length).ToList();

                    if (to.Count == 1)
                    {
                        var target = to.First();
                        
                        target.Item1.AddSequence(target.Item2, sq.Length, sequenceCounts[sq]);
                        continue;
                    }
                    
                    //aligned
                    var ratios = new Dictionary<(Alignment, int), float>();
                    foreach (var target in to) //todo nerozdelovat rovnomerne, ale v pomeru -- rozhoduje misto, kam se alignuje
                    {
                        var alignment = target.Item1;
                        var offset = target.Item2;

                        var countsOnTarget = alignment.GetCounts();
                        float sum = 0f;
                        for (int i = offset; i < (offset+sq.Length); i++)
                            sum += countsOnTarget[i];

                        float mean = sum / sq.Length;
                        ratios.Add(target, mean);
                    }
                    
                    //sort by ratio (descending), add 
                    var sorted_ratios = ratios.OrderByDescending(kv => kv.Value).ToList();
                    float part = sequenceCounts[sq] / ratios.Values.Sum();

                    var used = sequenceCounts[sq];
                    
                    foreach (var kv in sorted_ratios)
                    {
                        var alignment = kv.Key.Item1;
                        var offset = kv.Key.Item2;
                        var ratio = kv.Value;

                        var count = (int)Math.Round(ratio * part);
                        if (count >= used)
                        {
                            alignment.AddSequence(offset, sq.Length, used);
                            break;
                        }
                        alignment.AddSequence(offset, sq.Length, count);
                        used -= count;
                    }
                }

                // unaligned
                var unaligned = sequenceCounts.Keys.Except(aligned);
                foreach (var sq in unaligned)
                {
                    against.Add(
                        new Alignment(sq, sequenceCounts[sq])
                        );
                }
            }
        }
    }
}
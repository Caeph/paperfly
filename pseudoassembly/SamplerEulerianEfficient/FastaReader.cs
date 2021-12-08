using System.Collections.Generic;
using System.IO;

namespace SamplerEulerian
{
    public class FastaEntry
    {
        public string Header { get; }
        public string Sequence { get; }
        public FastaEntry(string header, string sequence)
        {
            Header = header;
            Sequence = sequence;
        }
    }

    public class FastaReader
    {
        public static IEnumerable<FastaEntry> Read(string fastaname)
        {
            using (var fileStream = File.OpenRead(fastaname))
            using (StreamReader streamReader = new StreamReader(fileStream)) {
                string header;
                string seq;
                while ((header = streamReader.ReadLine()) != null)
                {
                    seq = streamReader.ReadLine();
                    var entry = new FastaEntry(header.Substring(1), seq);
                    yield return entry;
                }

                // Process line
            }
        }
    }
}
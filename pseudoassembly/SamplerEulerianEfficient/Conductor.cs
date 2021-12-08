using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace SamplerEulerian
{
    public class Conductor
    {
        public void Run(string ComponentsPath, 
            string samplerType,
            int reportStep, 
            string outputdirname,
            int k=-1
        )
        {
            var overallSw = Stopwatch.StartNew();
            
            var files = Directory.GetFiles(ComponentsPath);
            Directory.CreateDirectory(outputdirname);

            foreach (string component in files)
            {
                string outputName;
                if (component.Contains("singles"))
                {
                    outputName = outputdirname + "/" + component.Split('/').Last().Replace("singles", "x") + ".csv";
                    using (var writer = new StreamWriter(outputName))
                    {
                        foreach (var entry in FastaReader.Read(component))
                        {
                            string[] splits = entry.Header.Split(new char[]{' '}, StringSplitOptions.RemoveEmptyEntries);

                            string count = splits[4].Split(',').First();
                            writer.WriteLine($"{entry.Sequence};{count}");
                        }
                    }
                    continue;
                }
                
                Console.WriteLine(component);
                var graph = new DeBruijnGraph(component);
                var sampler = Program.SamplerFactory(samplerType, graph, k);
                outputName = outputdirname + "/" + component.Split('/').Last() + ".csv";
                sampler.Sample(outputName, reportStep);
                
                //i++;
            }
            
            overallSw.Stop();
            Console.WriteLine(
                $"{ComponentsPath}: total time elapsed={overallSw.Elapsed.ToString()}"
            );
        }
    }
}
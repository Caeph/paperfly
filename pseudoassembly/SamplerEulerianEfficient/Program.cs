using System;
using System.Collections.Generic;

/// <summary>
/// --k 60 --input_path /Users/kf/UOCHB/phd/eulerianSampler/components
/// </summary>

namespace SamplerEulerian
{
    public static class IEnumerableExt
    {
        // usage: someObject.SingleItemAsEnumerable();
        public static IEnumerable<T> SingleItemAsEnumerable<T>(this T item)
        {
            yield return item; 
        }
        
        public static TSource MaxBy<TSource, TProperty>(this IEnumerable<TSource> source,
            Func<TSource, TProperty> selector)
        {
            // check args        

            using (var iterator = source.GetEnumerator())
            {
                if (!iterator.MoveNext())            
                    throw new InvalidOperationException();

                var max = iterator.Current; 
                var maxValue = selector(max);
                var comparer = Comparer<TProperty>.Default;

                while (iterator.MoveNext())
                {
                    var current = iterator.Current;
                    var currentValue = selector(current);

                    if (comparer.Compare(currentValue, maxValue) > 0)
                    {
                        max = current;
                        maxValue = currentValue;
                    }
                }

                return max;
            }
        }
    }
    
    internal class Program
    {
        public static ISampler SamplerFactory(string samplerType, DeBruijnGraph graph, int k=0)
        {
            ISampler sampler;
            switch (samplerType)
            {
                case "sequence":
                    if (k <= 0)
                    {
                        Console.WriteLine("The *k* parameter is needed to use the sequence sampler, exiting...");
                        Environment.Exit(1);
                    }
                    sampler = new BaseSampler(graph, k);
                    break;
                default:
                    throw new NotImplementedException();
            }

            return sampler;
        }
        
        public static void Main(string[] args)
        {
            string inputDir=null;
            int samplingDepth=-1;
            string outputDir=null;
            int k=-1; //negative if undefined
            string samplerType="sequence";
            int report_step = 100;
            bool canonical=false;
            bool singleItem = false;
            Dictionary<string, Action<string>> paramAction = new Dictionary<string, Action<string>>
            {
                {"--input_path", (val) => { inputDir = val;}},
                {"--samplingDepth", s => {
                    if (!int.TryParse(s, out samplingDepth))
                    {
                        Console.WriteLine("It must be possible to parse sampling depth to an integer, exiting...");
                        Environment.Exit(1);
                    }
                }},
                {"--k", s => {
                    if (!int.TryParse(s, out k))
                    {
                        Console.WriteLine("It must be possible to parse k to an integer, exiting...");
                        Environment.Exit(1);
                    }
                }},
                {"--output_path", s => { outputDir = s;}},
                {"--report_step", s => { if (!int.TryParse(s, out report_step))
                {
                    Console.WriteLine("It must be possible to parse report_step to an integer, exiting...");
                    Environment.Exit(1);
                }}},
                {"--canonical", s => { canonical = bool.Parse(s); }},
                {"--single_item", s => { singleItem = bool.Parse(s);}}
            };
            
            if ((args.Length % 2 != 0) || (args.Length == 0))
            {
                Console.WriteLine("Incomplete list of arguments, exiting...");
                Environment.Exit(1);
            }
            
            for (int i = 0; i < args.Length; i += 2)
            {
                //EXPECTATION -- arguments go in pairs --arg value
                string param = args[i];
                string value = args[i + 1];

                if (!paramAction.ContainsKey(param))
                {
                    Console.WriteLine($"Unknown argument {param}, exiting...");
                    Environment.Exit(1);
                }
                paramAction[param](value);
            }

            Console.WriteLine("RUNNING WITH PARAMS:");
            Console.WriteLine(String.Join(";", args));

            if (singleItem)
            {
                var graph = new DeBruijnGraph(inputDir);
                
                ISampler sampler;
                sampler = SamplerFactory(samplerType, graph, k);
                sampler.Sample(outputDir, report_step);
            }
            else
            {
                Conductor conductor = new Conductor();
                conductor.Run(inputDir, samplerType, report_step, outputDir, k);   
            }
        }
    }
}
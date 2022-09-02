using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using QuickGraph;

namespace SamplerEulerian
{
    public class PathInGraph
    {
        public IEnumerable<int> PathVertices { get; private set; }
        public PathInGraph(IEnumerable<int> path)
        {
            PathVertices = path;
        }
    }

    public class ProcessedPathInGraph : PathInGraph
    {
        public int Count { get; private set; }
        public string Sequence { get; private set; }

        public ProcessedPathInGraph(IEnumerable<int> path, int count, string seq) : base(path)
        {
            Count = count;
            Sequence = seq;
        }
    }

    public class PathInCondensation
    {
        public PathInCondensation(IEnumerable<CondensationNode> path)
        {
            PathVertices = path;
        }
        public IEnumerable<CondensationNode> PathVertices { get; private set; }
    }

    interface ISampler
    {
        void Sample(string outputPath, int report_step);
    }
    
    public class BaseSampler : ISampler
    {
        protected DeBruijnGraph deBruijnGraph;
        //protected BidirectionalGraph<int, Edge<int>> graph;
        protected Condensation structures;
        protected int k;
        protected List<(int, string, int[])> isolatedVertices = new List<(int, string, int[])>();
        public BaseSampler(DeBruijnGraph deBruijnGraph, int k)
        {
            this.deBruijnGraph = deBruijnGraph;
            //graph = deBruijnGraph.Graph;
            this.k = k;
        }

        private Stopwatch sw;

        private void report(int reportStep, int current)
        {
            if (current % reportStep == 0)
            {
                sw.Stop();
                Console.WriteLine($"{deBruijnGraph.Title}:{reportStep.ToString()} steps done in {sw.Elapsed.ToString()}, graph has {deBruijnGraph.VertexCount.ToString()} vertices and {deBruijnGraph.EdgeCount.ToString()} edges.");
                sw = Stopwatch.StartNew();
            }
        }
        
        public void Sample(string outputPath, int report_step)
        {
            //prep: construct condensation
            Console.WriteLine($"{deBruijnGraph.Title} processing started, graph has {deBruijnGraph.VertexCount.ToString()} vertices and {deBruijnGraph.EdgeCount.ToString()} edges.");
            prep();
            Console.WriteLine($"{deBruijnGraph.Title} preparation done.");
            using (var writer = new StreamWriter(outputPath)) //no appending to previous run
            {}

            int i = 0;
            while (deBruijnGraph.Graph.EdgeCount > 0)
            {
                int bottleneckCount;
                var paths = getPathsFromIteration(out bottleneckCount);
                var processedPaths = processPaths(paths, bottleneckCount).ToList();
                write(processedPaths, outputPath);

                actualize(processedPaths, outputPath);
                report(report_step, i);

                i++;
            }
            processLeftovers(outputPath);
        }

        protected void writeIsolated(string outputPath)
        {
            using (var writer = File.AppendText(outputPath))
            {
                foreach (var isolated in isolatedVertices)
                {
                    foreach (var result in pathFromNode(isolated.Item1, isolated.Item2, isolated.Item3))
                    {
                        writer.WriteLine($"{result.Item1};{result.Item2.ToString()}");
                        //Console.WriteLine($"{result.Item1};{result.Item2.ToString()}");
                    }
                }
            }

            isolatedVertices = new List<(int, string, int[])>();
        }

        protected void processLeftovers(string outputPath)
        {
            using (var writer = File.AppendText(outputPath))
            {
                foreach (var isolated in isolatedVertices)
                {
                    foreach (var result in pathFromNode(isolated.Item1, isolated.Item2, isolated.Item3))
                    {
                        writer.WriteLine($"{result.Item1};{result.Item2.ToString()}");
                        //Console.WriteLine($"{result.Item1};{result.Item2.ToString()}");
                    }
                }

                foreach (var vertex in deBruijnGraph.Graph.Vertices)
                {
                    foreach (var result in pathFromNode(vertex, deBruijnGraph.Sequences[vertex], deBruijnGraph.AllCounts[vertex]))
                    {
                        writer.WriteLine($"{result.Item1};{result.Item2.ToString()}");
                        //Console.WriteLine($"{result.Item1};{result.Item2.ToString()}");
                    }
                }
            }
        }

        private List<(string, int[])> getFragments(List<int> updatedCounts, string originalSequence)
        {
            var result = new List<(string, int[])>();
            if (updatedCounts.All(v => v == 0))
                return result;

            int firstPresent = -1;
            int extensions = 0;
            for (int i = 0; i < updatedCounts.Count; i++)
            {
                if (updatedCounts[i] == 0)
                {
                    if (firstPresent >= 0)
                    {
                        var subsequence = originalSequence.Substring(firstPresent, k + extensions);
                        var subcounts = updatedCounts.GetRange(firstPresent, extensions+1).ToArray();
                        result.Add((subsequence, subcounts));
                        
                        firstPresent = -1;
                        extensions = 0;
                    }
                }
                else
                {
                    if (firstPresent < 0) firstPresent = i;
                    else extensions++;
                }
            }
            if (firstPresent >= 0)
            {
                var subsequence = originalSequence.Substring(firstPresent, k + extensions);
                var subcounts = updatedCounts.GetRange(firstPresent, extensions+1).ToArray();
                result.Add((subsequence, subcounts));
            }

            return result;
        }

        private IEnumerable<(string,int)> pathFromNode(int vertex, string sequence, int[] originalCounts)
        {
            var stack = new Stack<(string, int[])>();
            stack.Push((sequence, originalCounts));

            while (stack.Count > 0)
            {
                var current = stack.Pop();
                int bottleneck = current.Item2.Min();
                yield return (current.Item1, bottleneck);

                var newCounts = current.Item2.Select(c => c - bottleneck).ToList();
                var fragments = getFragments(newCounts, current.Item1);
                foreach (var fragment in fragments)
                {
                    stack.Push(fragment);
                }
            }
        }

        protected IEnumerable<ProcessedPathInGraph> processPaths(IEnumerable<PathInGraph> paths, int bottleneckCount)
        {
            var unprocessed = paths.Where(path =>
                path.PathVertices.Min(v => deBruijnGraph.MinCounts[v])  == bottleneckCount).ToList();
            //heuristic: those that have a more narrow bottleneckCount are less likely to be seen, so these paths are thrown away.
            //there is always at least one with this bottleneck

            int fairShare = bottleneckCount / unprocessed.Count;
            int leftovers = bottleneckCount % unprocessed.Count;
            
            //take first <leftovers> paths, leftovers < unprocessed.Count ALWAYS
            for (int i = 0; i < leftovers; i++)
            {
                var newPath = new ProcessedPathInGraph(unprocessed[i].PathVertices, fairShare + 1,
                    StringTools.DebruijnJoin(
                        unprocessed[i].PathVertices.Select(v => deBruijnGraph.Sequences[v]).ToList(), k
                    ));
                yield return newPath;
            }

            for (int i = leftovers; i < unprocessed.Count; i++)
            {
                var newPath = new ProcessedPathInGraph(unprocessed[i].PathVertices, fairShare,
                    StringTools.DebruijnJoin(
                        unprocessed[i].PathVertices.Select(v => deBruijnGraph.Sequences[v]).ToList(), k
                    ));
                yield return newPath;
            }
        }

        protected void write(IEnumerable<ProcessedPathInGraph> processedPaths, string outputPath)
        {
            using (var writer = File.AppendText(outputPath))
            {
                foreach (var path in processedPaths)
                {
                    
                    writer.WriteLine($"{path.Sequence};{path.Count.ToString()}");
                    //Console.WriteLine($"{path.Sequence};{path.Count.ToString()}");
                }
            }
        }

        protected void actualize(IEnumerable<ProcessedPathInGraph> processedPaths, string outputPath)
        {
            foreach (var path in processedPaths)
            {
                var count = path.Count;

                //if (path.PathVertices.Distinct().Count() != path.PathVertices.Count())
                //{ }

                foreach (var vertex in path.PathVertices)
                {
                   deBruijnGraph.UpdateCount(vertex, count); 
                }
            }
            
            //find filled nodes (can occur IN MIDDLE OF CONTRACTED NODE)
            var toUpdate = deBruijnGraph.Graph.Vertices.Where(v => deBruijnGraph.MinCounts[v] == 0).ToList();
            //deal with them -- remove/split vertex
            var split = new Dictionary<int, int>(); //original --> last, those in the middle are isolated --> added to isolatedVertices collection and removed
            var removedVertices = new List<int>();
            var removedEdges = new List<Edge<int>>();

            bool constructNew = false;

            foreach (var vertex in toUpdate)
            {
                var counts = deBruijnGraph.AllCounts[vertex];
                if (structures.IsBigComponent(vertex)) constructNew = true;
                //if totally filled --> remove vertex
                if (counts.All(c => c == 0))
                {
                    deBruijnGraph.RemoveVertex(vertex);
                    removedVertices.Add(vertex);
                    continue;
                }
                //otherwise --> shatter
                var splitResults = deBruijnGraph.ShatterVertex(vertex, k);

                removedEdges.AddRange(splitResults.Item1);

                var newVertices = splitResults.Item2;
                if (newVertices.Length > 1) //otherwise only a single one, no change in graph vertices (only in edges)
                {
                    for (int i = 1; i < newVertices.Length - 1; i++)
                    {
                        isolatedVertices.Add((newVertices[i], deBruijnGraph.Sequences[newVertices[i]], deBruijnGraph.AllCounts[newVertices[i]]));
                        deBruijnGraph.RemoveVertex(newVertices[i]);
                    }
                    split.Add(vertex, newVertices[newVertices.Length - 1]);
                }
            }

            foreach (var isolatedVertex in deBruijnGraph.Graph.Vertices.Where(v => (deBruijnGraph.Graph.OutDegree(v) == 0) && (deBruijnGraph.Graph.InDegree(v) == 0)).ToList())
            {
                isolatedVertices.Add((isolatedVertex, deBruijnGraph.Sequences[isolatedVertex], deBruijnGraph.AllCounts[isolatedVertex]));
                deBruijnGraph.RemoveVertex(isolatedVertex);
                removedVertices.Add(isolatedVertex);
            }
            writeIsolated(outputPath);
            
            //update condensation
            //removed vertices, edges -- if not in scc, remove, otherwise construct new
            //split vertices -- if not in scc, split, otherwise construct new
            
            if (constructNew)
            {
                //structures = Condensation.Construct(deBruijnGraph, k);
                structures = structures.UpdateReconstruct(removedVertices, removedEdges, split, deBruijnGraph, k);
            }
            else
            {
                structures.Update(removedVertices, removedEdges, split);
            }
        }
        
        protected void prep()
        {
            structures = Condensation.Construct(deBruijnGraph, k);
            sw = Stopwatch.StartNew();
        }
        
        protected IEnumerable<PathInGraph> getPathsFromIteration(out int bottleneckCount)
        {
            //find longest paths in condensation 
            if (structures.CondensedGraph.EdgeCount == 0)
            {
                int maxTraversal = -1;
                CondensationNode argmax = null;
                foreach (var node in structures.CondensedGraph.Vertices)
                {
                    var traversal = structures.GetTraversalCost(node);
                    if (traversal > maxTraversal)
                    {
                        argmax = node;
                        maxTraversal = traversal;
                    }
                }
                var pathInGraph = new PathInGraph(argmax.Vertices);
                bottleneckCount = argmax.Vertices.Min(v => deBruijnGraph.MinCounts[v]);
                return new [] { pathInGraph };
            }
                
            var paths = get_paths_from_DAG(out bottleneckCount);
            return paths;
        }
        
        protected IEnumerable<PathInGraph> get_paths_from_DAG(out int bottleneckCount) // is potentially exponential BUT that is managed by LINQ Take(n) in upper scope
        {
            var toposortDistancesPredecessors = toposort_distances();
            var distances = toposortDistancesPredecessors.Item1;
            var previous_on_longest = toposortDistancesPredecessors.Item2;
            int longest_path_lenghts = distances.Values.Max();
            
            //get bottleneck and bottleneck count B
            var longest = extractPaths(distances, structures.TopologicalSort, previous_on_longest, longest_path_lenghts).First(); //! lazy evaluation -> O(n)
            //var toReturn = expand_paths(longest, structures, samplingDepth).First();
            //get B paths passing through bottleneck
            var bottleneck = getBottleneck(longest);
            var bottleneckVertex = bottleneck.Item1;
            bottleneckCount = bottleneck.Item2;

            var trails = enumerateBottleneckedPaths(bottleneckVertex, distances,
                previous_on_longest, longest_path_lenghts).Take(bottleneckCount).ToList(); //max amount of paths to take, not necessarily all of them occur in the result
            //doubling paths???
            var toReturn = expandPaths(trails);
            return toReturn.ToList(); 
        }

        protected PathInGraph expandPath(PathInCondensation condensationPath)
        {
            var graphVertices = new List<int>();
            foreach (var condensationNode in condensationPath.PathVertices)
            {
                if (typeof(GateNode) == condensationNode.GetType())
                {
                    continue; //that is empty
                }
                var walkthrough = condensationNode.Vertices;
                graphVertices.AddRange(walkthrough);
            }

            return new PathInGraph(graphVertices);
        }

        protected IEnumerable<PathInGraph> expandPaths(IEnumerable<PathInCondensation> condensationPaths)
        {
            foreach (var condensationPath in condensationPaths)
            {
                yield return expandPath(condensationPath);
            }
        }

        protected (CondensationNode, int) getBottleneck(PathInCondensation path) //(vertex in condensation, bottleneck count)
        {
            var savedPath = path.PathVertices.Where(v => typeof(GateNode) != v.GetType()).ToList();

            int bottleneckWidth = int.MaxValue;
            CondensationNode componentWithBottleneck = null;
            
            foreach (var condNode in savedPath)
            {
                var lowestCount = condNode.Vertices.Min(v => deBruijnGraph.MinCounts[v]);
                if (lowestCount < bottleneckWidth)
                {
                    bottleneckWidth = lowestCount;
                    componentWithBottleneck = condNode;
                }
            }

            return (componentWithBottleneck, bottleneckWidth);
        }
        
        protected IEnumerable<PathInCondensation> enumerateBottleneckedPaths(CondensationNode bottleneckNode,
            Dictionary<CondensationNode,int> distances,
            Dictionary<CondensationNode, List<CondensationNode>> prevOnLongest, 
            int longestLen)
        {
            //get path to bottleneck
            var to_bottleneck = getPathsTo(bottleneckNode, prevOnLongest); //includes the bottleneck node, LAZY! -- will not calculate what it does not have to
            var previous = constructPrevious(bottleneckNode);
            //get path from bottleneck
            var from_bottleneck = getPathsFrom(bottleneckNode, previous.Item1,
                structures.TopologicalSort.FindAll(r => distances[r] == longestLen), previous.Item2);
            
            //all combinations
            foreach (var to_path in to_bottleneck)
            {
                foreach (var from_path in from_bottleneck)
                {
                    var path = to_path.Concat(from_path.Skip(1)).ToList();
                    yield return new PathInCondensation(path);
                }
            }
        }
        
        
        private IEnumerable<IEnumerable<CondensationNode>> getPathsTo(CondensationNode bottleneckNode, 
            Dictionary<CondensationNode, List<CondensationNode>> previous)
        {
            var stack = new Stack<List<CondensationNode>>();
            stack.Push(new List<CondensationNode>() {bottleneckNode});

            while (stack.Count > 0)
            {
                var to_extend = stack.Pop(); // entire path
                var last_in_path = to_extend[to_extend.Count - 1];
                var prevs = previous[last_in_path];

                if (prevs.Count == 0)
                {
                    to_extend.Reverse();
                    yield return to_extend;
                }
                else
                {
                    foreach (var prev in prevs)
                    {
                        List<CondensationNode> path = new List<CondensationNode>(to_extend); 
                        path.Add(prev);
                        stack.Push(path);
                    }
                }
            }
        }

        private (Dictionary<CondensationNode, List<CondensationNode>>, HashSet<CondensationNode>) constructPrevious(
            CondensationNode bottleneckNode)
        {
            var stack = new Stack<CondensationNode>();
            stack.Push(bottleneckNode);

            HashSet<CondensationNode> seen = new HashSet<CondensationNode>() {bottleneckNode};

            while (stack.Count > 0)
            {
                var node = stack.Pop();
                
                foreach (var nextNode in structures.CondensedGraph.OutEdges(node).Select(x => x.Target))
                {
                    if (!seen.Contains(nextNode))
                    {
                        seen.Add(nextNode);
                        stack.Push(nextNode);
                    }
                }
            }

            var shortToposort = structures.TopologicalSort.Where(x => seen.Contains(x)).ToList();
            var desc = toposort_distances(shortToposort, seen);
            
            return (desc.Item2, seen);
        }

        private IEnumerable<IEnumerable<CondensationNode>> getPathsFrom(CondensationNode bottleneckNode, 
            Dictionary<CondensationNode,List<CondensationNode>> previous, 
            IEnumerable<CondensationNode> finishers, HashSet<CondensationNode> seen)
        {
            // hledam nejdelsi cestu z bottlenecku do jednoho z finisheru

            var stack = new Stack<List<CondensationNode>>();
            foreach (var node in finishers)
            {
                if (previous.ContainsKey(node))
                    stack.Push(new List<CondensationNode>() {node});
            }

            while (stack.Count > 0)
            {
                var to_extend = stack.Pop(); // entire path
                var last_in_path = to_extend[to_extend.Count - 1];
                var prevs = previous[last_in_path];

                if (prevs.Count == 0)
                {
                    to_extend.Reverse();
                    yield return to_extend;
                }
                else
                {
                    foreach (var prev in prevs)
                    {
                        if (!previous.ContainsKey(prev))
                            continue;
                        List<CondensationNode> path = new List<CondensationNode>(to_extend);

                        path.Add(prev);
                        stack.Push(path);
                    }
                }
            }
        }
        
        protected IEnumerable<PathInCondensation> extractPaths(
            Dictionary<CondensationNode,int> dsts, 
            List<CondensationNode> toposort,
            Dictionary<CondensationNode, List<CondensationNode>> previous, int longest_path_len
        )
        {
            //int counter = 0;

            var finishers = toposort.FindAll(r => dsts[r] == longest_path_len
            );

            var stack = new Stack<List<CondensationNode>>();

            foreach (var comp in finishers)
                stack.Push(new List<CondensationNode>() {comp});

            while (stack.Count > 0)
            {
                var to_extend = stack.Pop(); // entire path
                var last_in_path = to_extend[to_extend.Count - 1];
                var prevs = previous[last_in_path];

                if (prevs.Count == 0)
                {
                    to_extend.Reverse();
                    yield return new PathInCondensation(to_extend);
                }
                else
                {
                    foreach (var prev in prevs)
                    {
                        List<CondensationNode> path = new List<CondensationNode>(to_extend); 
                        path.Add(prev);
                        stack.Push(path);
                    }
                }
                
            }
        }
        
        
        
        protected (Dictionary<CondensationNode, int>, Dictionary<CondensationNode, List<CondensationNode>>) 
            toposort_distances(List<CondensationNode> toposort = null, HashSet<CondensationNode> accessible = null)
        {
            Func<CondensationNode, IEnumerable<Edge<CondensationNode>>> inEdges;
            if (toposort == null)
            {
                toposort = structures.TopologicalSort;
            }

            if (accessible == null)
            {
                inEdges = structures.CondensedGraph.InEdges;
                accessible = new HashSet<CondensationNode>();
            }
            else
            {
                inEdges = i => structures.CondensedGraph.InEdges(i).Where(x => accessible.Contains(x.Source));
            }
            
            Dictionary<CondensationNode, int> distances = new Dictionary<CondensationNode, int>();
            foreach (var component in toposort)
                distances[component] = 0;
            Dictionary<CondensationNode, List<CondensationNode>> previous_on_longest = new Dictionary<CondensationNode, List<CondensationNode>>();
            
            for (int i = 0; i < toposort.Count; i++)
            {
                var component = toposort[i];
                int furthest_neighbor_dst = 0;
                var preds = new List<CondensationNode>();

                var a = structures.CondensedGraph.ContainsVertex(component);

                foreach (var pred_e in inEdges(component))
                {
                    var pred = pred_e.Source;
                    int dst = distances[pred] + structures.GetTraversalCost(component);
                    if (furthest_neighbor_dst == dst)
                        preds.Add(pred);
                    if (furthest_neighbor_dst < dst)
                    {
                        preds.Clear();
                        preds.Add(pred);
                        furthest_neighbor_dst = dst;
                    }
                }
                
                distances[component] = furthest_neighbor_dst;
                previous_on_longest[component] = preds.OrderByDescending(
                    x => structures.GetTraversalCost(x)
                    ).ToList();
            }

            return (distances, previous_on_longest);
        }
    }
}
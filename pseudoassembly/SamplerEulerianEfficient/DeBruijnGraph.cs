using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using QuickGraph;

namespace SamplerEulerian
{
    public class DeBruijnGraph
    {
        public int VertexCount
        {
            get
            { return Graph.VertexCount; }
        }
        
        public string Title { get; }
        
        public int EdgeCount
        {
            get
            { return Graph.EdgeCount; }
        }

        public bool IsDirected
        {
            get { return Graph.IsDirected; }
        }

        public void RemoveVertices(IEnumerable<int> vertices)
        {
            foreach (var item in vertices)
            {
                RemoveVertex(item);
            }
        }

        public (List<Edge<int>>, int[]) ShatterVertex(int originalVertex, int k)
        {
            List<Edge<int>> removedEdges = new List<Edge<int>>();
            
            List<int[]> allNewCounts = new List<int[]>();
            List<int> nextNodeCounts = new List<int>();

            var counts = AllCounts[originalVertex].ToList();
            
            //get sequence fragments
            var sequenceFragments = new List<string>();
            var originalSequence = Sequences[originalVertex];
            int firstPresent = -1;
            int extensions = 0;
            for (int i = 0; i < counts.Count; i++)
            {
                if (counts[i] == 0)
                {
                    if (firstPresent >= 0)
                    {
                        var subsequence = originalSequence.Substring(firstPresent, k + extensions);
                        sequenceFragments.Add(subsequence);
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
                sequenceFragments.Add(subsequence);
            }

            //remove zeros on edges
            //zeros on edge negate edges in/out!
            var zerosLeft = Enumerable.Range(0, counts.Count).TakeWhile(i => counts[i] == 0);
            if (zerosLeft.Any())
            {
                removedEdges.AddRange(Graph.InEdges(originalVertex));
                RemoveEdges(Graph.InEdges(originalVertex).ToList()); 
                var firstNonZero = zerosLeft.Max() + 1;
                counts = counts.GetRange(firstNonZero, counts.Count - firstNonZero);
            }
                
            var zerosRight = Enumerable.Range(0, counts.Count).Reverse().TakeWhile(i => counts[i] == 0);
            if (zerosRight.Any())
            {
                removedEdges.AddRange(Graph.OutEdges(originalVertex));
                RemoveEdges(Graph.OutEdges(originalVertex).ToList());
                var lastNonZero = zerosRight.Min();
                counts = counts.GetRange(0, lastNonZero);
            }
            
            for (int i = 0; i < counts.Count; i++)
            {
                if (counts[i] == 0)
                {
                    if (nextNodeCounts.Count > 0)
                    {
                        allNewCounts.Add(nextNodeCounts.ToArray());
                        nextNodeCounts = new List<int>();
                    }
                }
                else
                {
                    nextNodeCounts.Add(counts[i]);
                }
            }
            if (nextNodeCounts.Count > 0) allNewCounts.Add(nextNodeCounts.ToArray());

            var newVertices = new List<int>();
            
            var last = allNewCounts[allNewCounts.Count - 1];
            var lastVertex = addNewVertex();
            AllCounts[lastVertex] = last;
            MinCounts[lastVertex] = last.Min();
            Sequences[lastVertex] = sequenceFragments[sequenceFragments.Count - 1];
            var outEdges = Graph.OutEdges(originalVertex).ToList();
            foreach (var edge in outEdges)
            {
                Graph.AddEdge(new Edge<int>(lastVertex, edge.Target));
            }
            
            var first = allNewCounts[0];
            RemoveEdges(outEdges);
            AllCounts[originalVertex] = first;
            MinCounts[originalVertex] = first.Min();
            Sequences[originalVertex] = sequenceFragments[0]; 
            newVertices.Add(originalVertex);

            int j = 1;
            if (allNewCounts.Count > 1)
            {
                foreach (var shard in allNewCounts.GetRange(1, allNewCounts.Count - 2))
                {
                    var vertex = addNewVertex();
                    AllCounts[vertex] = shard;
                    MinCounts[vertex] = shard.Min();
                    Sequences[vertex] = sequenceFragments[j]; 
                    newVertices.Add(vertex);
                    j++;
                }
            }
            newVertices.Add(lastVertex);
            return (removedEdges, newVertices.ToArray());
        }

        public void RemoveVertex(int item)
        {
            Graph.RemoveVertex(item);
            MinCounts.Remove(item);
            AllCounts.Remove(item);
            Sequences.Remove(item);
        }

        private int addNewVertex()
        {
            int newV = maxVertexNum + 1;
            maxVertexNum++;
            
            if (newV == 10272)
            {}

            Graph.AddVertex(newV);
            return newV;
        }
        
        public void RemoveEdges(IEnumerable<Edge<int>> edges)
        {
            foreach (var item in edges)
            {
                Graph.RemoveEdge(item);
            }
        }

        public void UpdateCount(int vertex, int count)
        {
            if (MinCounts[vertex] < count)
                throw new InvalidOperationException("Count of observations must be non-negative.");
            
            MinCounts[vertex] -= count;
            AllCounts[vertex] = AllCounts[vertex].Select(c => c - count).ToArray();
        }
        
        public BidirectionalGraph<int, Edge<int>> Graph { get; }
        private int maxVertexNum;
        public Dictionary<int, int> MinCounts { get; }
        public Dictionary<int, string> Sequences { get; }

        public Dictionary<int, int[]> AllCounts { get; }
        public DeBruijnGraph(string fastafilename)
        {
            Title = fastafilename;
            MinCounts = new Dictionary<int, int>();
            Graph = new BidirectionalGraph<int, Edge<int>>();
            Sequences = new Dictionary<int, string>();
            AllCounts = new Dictionary<int, int[]>();

            var fasta_reader = FastaReader.Read(fastafilename);

            var edges = new List<Edge<int>>();
            
            foreach (var entry in fasta_reader)
            {
                int id;
                int minCount = int.MaxValue;

                string[] splits = entry.Header.Split(new char[]{' '}, StringSplitOptions.RemoveEmptyEntries);
                id = int.Parse(splits[0]);

                string[] counts = splits[4].Split(',');
                
                for (int i = 0; i < counts.Length; i++)
                {
                    int now_count = int.Parse(counts[i]);
                    if (minCount > now_count)
                        minCount = now_count;
                }
                MinCounts.Add(id, minCount);
                Sequences.Add(id, entry.Sequence);

                var int_counts = counts.Select(c => int.Parse(c)).ToArray();
                
                AllCounts.Add(id, int_counts);
                for (int i = 5; i < splits.Length; i++)
                {
                    string[] edge_spec = splits[i].Split(':');
                    int out_id = int.Parse($"{edge_spec[3]}{edge_spec[2]}");

                    if (id != out_id)
                    {
                        if (!edges.Any(e => e.Source == out_id && e.Target == id)) 
                            edges.Add(new Edge<int>(id, out_id));  //removes bidirectional edges
                    }
                }
            }

            edges = edges.Where(e => e.Source != e.Target).ToList();

            //foreach (var edge in edges)
            //{
            //    var cont = edges.Any(e => e.Source == edge.Target && e.Target == edge.Source);
            //}
            
            Graph.AddVerticesAndEdgeRange(edges); //must be connected in the beginning
            maxVertexNum = Graph.Vertices.Max();
        }
    }
}
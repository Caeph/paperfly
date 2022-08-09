using System;
using System.Collections.Generic;
using System.Linq;
using QuickGraph;
using QuickGraph.Algorithms.ShortestPath;
using QuickGraph.Algorithms.Search;

namespace SamplerEulerian
{
    interface ITraversable
    {
        List<int> GetPassage(int beginning, int ending);
        List<int> GetCyclicTraversal(int beginning, bool reverse);
    }
    
    public class NonEulerianComponent : ITraversable
    {
        
        public List<int> GetPassage(int beginning, int ending)
        {
            if (beginning == ending)
            {
                return new List<int>() {beginning};
            }
            
            var vertices = new List<int>();

            var dfs = new DepthFirstSearchAlgorithm<int, Edge<int>>(subgraph);
            var dict = new Dictionary<int, Edge<int>>(); //

            dfs.TreeEdge += e =>
            {
                dict[e.Target] = e;
            };
            dfs.DiscoverVertex += v =>
            {
                if (v == ending)
                    dfs.Abort();
            };
            dfs.Compute(beginning);

            var edge = dict[ending];
            vertices.Add(ending);
            while (edge.Source != beginning)
            {
                vertices.Add(edge.Source);
                edge = dict[edge.Source];
            }
            vertices.Add(beginning);
            vertices.Reverse();
            
            /*var counts = vertices.GroupBy(v => v).Select(g => (g.Key, g.Count())).Any(x => x.Item2 > 1);
            if (counts)
            { }
            */

            return vertices;
        }

        public List<int> GetCyclicTraversal(int beginning, bool reverse)
        {
            var currentCircuit = new List<Edge<int>>();

            BidirectionalGraph<int, Edge<int>> subgGraph_oriented;
            if (reverse)
            {
                subgGraph_oriented = new BidirectionalGraph<int, Edge<int>>();
                foreach (var e in subgraph.Edges)
                {
                    subgGraph_oriented.AddVerticesAndEdge(new Edge<int>(e.Target, e.Source));
                }
            }
            else
            {
                subgGraph_oriented = subgraph;
            }

            var first = subgGraph_oriented.OutEdges(beginning).First();

            var seen = new HashSet<Edge<int>>();
            var seen_vertices = new HashSet<int>();

            seen_vertices.Add(beginning);
            
            var edge = first;
            while ((edge != null) && (edge.Target != first.Source))
            {
                currentCircuit.Add(edge);

                seen.Add(edge);
                seen_vertices.Add(edge.Target);
                
                var nextSource = edge.Target;
                edge = subgGraph_oriented.OutEdges(nextSource).FirstOrDefault(
                    e => (!seen.Contains(e)) & (!seen_vertices.Contains(e.Target))
                );
            }

            if (edge != null)
            {
                currentCircuit.Add(edge);
            }

            var result = currentCircuit.Select(e => e.Target).ToList();

            if (reverse)
            {
                result.Reverse();
                result.Add(beginning);
            }
            else
            {
                result.Insert(0,beginning);
            }
            
            return result;
        }

        private BidirectionalGraph<int, Edge<int>> subgraph;
        private Dictionary<int, Dictionary<int, Edge<int>>> accessEdgesRoots = new Dictionary<int, Dictionary<int, Edge<int>>>();

        public NonEulerianComponent(BidirectionalGraph<int, Edge<int>> subgraph)
        {
            this.subgraph = subgraph;
        }
    }

    public class EulerianCycle : ITraversable
    {
        private IEnumerable<int> yieldVerticesWithRepeats(int beginning, Predicate<Edge<int>> additionalEndingCondition)
        {
            var firstEdge = access[beginning];
            var index = edges.IndexOf(firstEdge);
            Edge<int> currentEdge;
            if (index == edges.Count - 1)
            {
                currentEdge = edges[0];
                index = 0;
            }
            else
            {
                currentEdge = edges[index + 1];
                index++;
            }

            yield return beginning;
            while (currentEdge != firstEdge)
            {
                yield return currentEdge.Source;

                //if (currentEdge.Source == end)
                //    yield break;
                if (additionalEndingCondition(currentEdge))
                    yield break;

                if (index == edges.Count - 1)
                {
                    currentEdge = edges[0];
                    index = 0;
                }
                else
                {
                    currentEdge = edges[index + 1];
                    index++;
                }
            }
        }

        public List<int> GetCyclicTraversal(int beginning, bool reverse=false)
        {
            var walkWithRepeatVertices = yieldVerticesWithRepeats(beginning, e => e.Source == beginning).ToList();

            var result = makeUnique(walkWithRepeatVertices);

            /*var counts = result.GroupBy(v => v).Select(g => (g.Key, g.Count())).Any(x => x.Item2 > 1);
            if (counts)
            { }*/
            //result.Add(result.First());
            
            return result;
        }

        protected List<int> makeUnique(List<int> walkWithRepeatVertices)
        {
            var result = walkWithRepeatVertices;

            var currentOccurences = result.GroupBy(v => v).ToDictionary(g => g.Key, g => g.Count());
            var touchpoint = result.Select(v => (int?)v).FirstOrDefault(v => currentOccurences[(int)v] > 1);

            HashSet<int> seen = new HashSet<int>();

            while (touchpoint != null)
            {
                var firsttouch =
                    result.IndexOf(touchpoint.Value); //can be beginning, not ending -- if ending is reached, enumeration stops
                var lasttouch = result.LastIndexOf(touchpoint.Value);

                if (firsttouch == lasttouch)
                {
                    seen.Add(touchpoint.Value);
                    touchpoint = result.Select(v => (int?)v).FirstOrDefault(v =>
                        (currentOccurences[(int)v] > 1) && (!seen.Contains(v.Value)));
                    continue;
                }

                var newResult = result.GetRange(0, firsttouch); //take only previous
                newResult.AddRange(
                    result.GetRange(lasttouch,
                        result.Count - lasttouch) //take touchpoint and the rest
                );

                result = newResult;

                currentOccurences = result.GroupBy(v => v).ToDictionary(g => g.Key, g => g.Count());

                touchpoint = result.Select(v => (int?)v).FirstOrDefault(v =>
                (currentOccurences[(int)v] > 1) && (!seen.Contains(v.Value)));
            }

            return result;
        }

        public List<int> GetPassage(int beginning, int ending)
        {
            var walkWithRepeatVertices =
                yieldVerticesWithRepeats(beginning, e => e.Source == ending)
                    .ToList(); //yields vertices as they are visited on the eulerian cycle
            var result = makeUnique(walkWithRepeatVertices);

            /*var counts = result.GroupBy(v => v).Select(g => (g.Key, g.Count())).Any(x => x.Item2 > 1);
            if (counts)
            { }*/

            return result;
        }

        private List<Edge<int>> edges;
        private Dictionary<int, Edge<int>> access = new Dictionary<int, Edge<int>>();
        private Dictionary<int, int> occurencesOnEul = new Dictionary<int, int>();

        public EulerianCycle(List<Edge<int>> edges)
        {
            foreach (var edge in edges)
            {
                access[edge.Source] = edge; //a single entrypoint is kept
                if (!occurencesOnEul.ContainsKey(edge.Source))
                    occurencesOnEul[edge.Source] = 1;
                else
                    occurencesOnEul[edge.Source] += 1;
            }

            this.edges = edges;
        }
    }
}
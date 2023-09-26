using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.IO;
using System.Linq;
using QuickGraph;
using QuickGraph.Algorithms;
using QuickGraph.Algorithms.AssigmentProblem;
using QuickGraph.Algorithms.ConnectedComponents;
using QuickGraph.Algorithms.ShortestPath;
using QuickGraph.Algorithms.TopologicalSort;

namespace SamplerEulerian
{
    
    [DebuggerDisplay("ID={ID}, scc={IsSCC}")]
    public class CondensationNode
    {
        public int ID { get; private set; }
        public bool IsSCC { get; private set; }
        public List<int> Vertices { get; private set; }
        
        public string Sequence { get; set; }
        public CondensationNode(int id, bool isSCC, List<int> graphVertices)
        {
            ID = id;
            IsSCC = isSCC;
            Vertices = graphVertices;
        }
    }

    [DebuggerDisplay("GATE {direction}: ID={ID}, scc={IsSCC}, node={RepresentedGraphNode}")]
    public class GateNode : CondensationNode
    {
        public int RepresentedGraphNode { get; }

        private string direction;
        //walkthrough cost is 0
        public GateNode(int id, bool isSCC, int representedNode, string direction) : base(id, isSCC, new List<int>())
        {
            RepresentedGraphNode = representedNode;
            this.direction = direction;
        }
    }

    [DebuggerDisplay("WALKTHROUGH {ID}: {InGraphNode}->{OutGraphNode}")]
    public class WalkthroughNode : CondensationNode
    {
        public int InGraphNode { get; }
        public int OutGraphNode { get; }
        public WalkthroughNode(int id, bool isSCC, int inGraphNode, int outGraphNode, List<int> walkthroughVertices) : base(id, isSCC, walkthroughVertices)
        {
            InGraphNode = inGraphNode;
            OutGraphNode = outGraphNode;
        }
    }

    public class WholeCycleNode : CondensationNode
    {
        public int Entry { get; }

        public WholeCycleNode(int id, bool isSCC, int entryVertex, List<int> walkthroughVertices) : base(id, isSCC,
            walkthroughVertices)
        {
            Entry = entryVertex;
        }
    }
    
    public class Condensation
    {
        public HashSet<int> BigComponents = new HashSet<int>();
        public StronglyConnectedComponentsAlgorithm<int, Edge<int>> StronglyConnected { get; private set; }
        public Dictionary<int, int> GraphVertexToComponent { get; private set; }
        public Dictionary<int, int[]> ComponentToMembers { get; private set; }
        public List<CondensationNode> TopologicalSort { get; private set; }
        public BidirectionalGraph<CondensationNode, Edge<CondensationNode>> CondensedGraph { get; private set; }

        private DeBruijnGraph _deBruijnGraph;
        private int k;
        private Dictionary<int, CondensationNode> smallComponentNodes;
        Dictionary<int, Dictionary<int, CondensationNode>> bigDictIn;
        Dictionary<int, Dictionary<int, CondensationNode>> bigDictOut;
        private Condensation(
            StronglyConnectedComponentsAlgorithm<int, Edge<int>> stronglyConnected,
            HashSet<int> bigComponents,
            Dictionary<int, int> graphVertexToComponent,
            List<CondensationNode> topologicalSort, 
            Dictionary<int, int[]> componentToMembers, 
            BidirectionalGraph<CondensationNode,
                Edge<CondensationNode>> condensedGraph,
            DeBruijnGraph deBruijnGraph, int k,
            Dictionary<int, CondensationNode> smallComponentNodes,
            Dictionary<int, Dictionary<int, CondensationNode>> bigDictIn,
            Dictionary<int, Dictionary<int, CondensationNode>> bigDictOut
        )
        {
            BigComponents = bigComponents;
            StronglyConnected = stronglyConnected;
            GraphVertexToComponent = graphVertexToComponent;
            TopologicalSort = topologicalSort;
            ComponentToMembers = componentToMembers;
            CondensedGraph = condensedGraph;
            _deBruijnGraph = deBruijnGraph;
            this.k = k;
            this.smallComponentNodes = smallComponentNodes;
            componentCount = StronglyConnected.ComponentCount;
            this.bigDictIn = bigDictIn;
            this.bigDictOut = bigDictOut;
        }

        /*public BidirectionalGraph<int, Edge<int>> GetSubgraph(int component)
        {
            return StronglyConnected.Graphs[component];
        }*/

        public bool IsBigComponent(int vertexInGraph)
        {
            return BigComponents.Contains(GraphVertexToComponent[vertexInGraph]);
        }

        public Condensation UpdateReconstruct(List<int> removedVertices, List<Edge<int>> removedEdges,
            Dictionary<int, int> splitVertices, DeBruijnGraph dbgraph, int k)
        {
            var newCond = construct(dbgraph, k);
            return newCond;
        }

        public void Update(List<int> removedVertices, List<Edge<int>> removedEdges,
            Dictionary<int, int> splitVertices)
        {
            //split vertices
            foreach (var removedEdge in removedEdges)
            {
                CondensationNode source;
                if (!GraphVertexToComponent.ContainsKey(removedEdge.Source))
                {
                    //edge is from shattered vertex, remove the original edge
                    var sourceVertex = removedEdge.Source;
                    if (!splitVertices.ContainsValue(sourceVertex))
                        throw new InvalidOperationException();
                    
                    //get original source
                    var originalSourceVertex = splitVertices.First(kv => kv.Value == sourceVertex).Key;
                    source = smallComponentNodes[GraphVertexToComponent[originalSourceVertex]];
                }
                //no such problem for target -- in edges are not moved
                else
                {
                    //source = smallComponentNodes[GraphVertexToComponent[removedEdge.Source]];
                    var sourceComponent = GraphVertexToComponent[removedEdge.Source];
                    if (BigComponents.Contains(sourceComponent))
                    {
                        if (bigDictOut[sourceComponent].ContainsKey(removedEdge.Source))
                        {
                            source = bigDictOut[sourceComponent][removedEdge.Source];
                        }
                        else continue;
                    }
                    else
                    {
                        source = smallComponentNodes[sourceComponent];
                    }
                }
                
                var targetComponent = GraphVertexToComponent[removedEdge.Target];
                CondensationNode targetNode;
                if (!BigComponents.Contains(targetComponent))
                    targetNode = smallComponentNodes[targetComponent];
                else
                {
                    var x = CondensedGraph.OutEdges(source).FirstOrDefault(e => e.Target.ID == targetComponent);
                    if (x == null) continue;
                    targetNode = x.Target;
                }
                //none of these two is a bigComponent
                //find edge, remove edge
                var edge = CondensedGraph.OutEdges(source).FirstOrDefault(e => e.Target == targetNode); //there is exactly one
                CondensedGraph.RemoveEdge(edge);
            }
            
            //shattered vertices
            foreach (var shattered in splitVertices)
            {
                var originalComponent = GraphVertexToComponent[shattered.Key]; //no change there -- vertex has the same number
                
                //edges
                var originalComponentNode = smallComponentNodes[originalComponent];
                var outEdges = CondensedGraph.OutEdges(originalComponentNode).ToList();
                foreach (var edge in outEdges) CondensedGraph.RemoveEdge(edge);

                //add the end vertex
                var newVertex = shattered.Value;

                var newComponent = componentCount + 1;
                componentCount++;
                var newNode = new CondensationNode(newComponent, false, new List<int>() { newVertex });

                CondensedGraph.AddVertex(newNode);
                foreach (var oldEdge in outEdges)
                {
                    CondensedGraph.AddEdge(new Edge<CondensationNode>(newNode, oldEdge.Target));
                }
                smallComponentNodes.Add(newComponent, newNode);
                GraphVertexToComponent.Add(newVertex, newComponent);
                ComponentToMembers.Add(newComponent, new int[] {newVertex});

                var oldVertexIndex = TopologicalSort.IndexOf(originalComponentNode);
                TopologicalSort.Insert(oldVertexIndex+1, newNode);
                //TopologicalSort.Remove(originalComponentNode);
            }
            
            //removed vertices
            var components = removedVertices.Select(v => GraphVertexToComponent[v]).Distinct().ToList();
            foreach (var component in components)
            {
                if (BigComponents.Contains(component))
                    throw new InvalidOperationException("Cannot remove a composite component");
                
                var member = ComponentToMembers[component][0]; //single item component
                GraphVertexToComponent.Remove(member);
                ComponentToMembers.Remove(component);
            }


            var nodes = components.Select(component => smallComponentNodes[component]).ToHashSet();
            TopologicalSort.RemoveAll(node => nodes.Contains(node));
            foreach (var node in nodes) CondensedGraph.RemoveVertex(node);
        }

        private int componentCount;

        public int GetTraversalCost(CondensationNode node, bool isFirst=false)
        {
            if (typeof(GateNode) == node.GetType())
            {
                return 0; //nothing here
            }

            if (node.Sequence == null)
            {
                var sequence = StringTools.DebruijnJoin(node.Vertices.Select(v => _deBruijnGraph.Sequences[v]).ToList(), k);
                node.Sequence = sequence;
                return sequence.Length;
            }
            else
            {
                return node.Sequence.Length;
            }
        }

        private static bool isEulerian(BidirectionalGraph<int, Edge<int>> componentGraph)
        {
            var isEul = componentGraph.Vertices.Any(v => componentGraph.InDegree(v) != componentGraph.OutDegree(v));
            return !isEul;
        }

        private static bool isIsolated(int vertex, BidirectionalGraph<int, Edge<int>> componentGraph)
        {
            if ((componentGraph.OutDegree(vertex) == 0) && (componentGraph.InDegree(vertex) == 0))
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        private static List<Edge<int>> getCircuit(BidirectionalGraph<int, Edge<int>> componentGraphCopy, Edge<int> first = null)
        {
            var currentCircuit = new List<Edge<int>>();

            if (first == null)
                first = componentGraphCopy.Edges.First();

            var edge = first;
            while ((edge != null) && (edge.Target != first.Source))
            {
                currentCircuit.Add(edge);
                componentGraphCopy.RemoveEdge(edge);

                var nextSource = edge.Target;
                edge = componentGraphCopy.OutEdges(nextSource).FirstOrDefault();
            }

            if (edge != null)
            {
                currentCircuit.Add(edge);
                componentGraphCopy.RemoveEdge(edge);
            }
            
            //remove all isolated vertices AT ONCE at the end
            var isolated = componentGraphCopy.Vertices.Where(v => (componentGraphCopy.OutDegree(v) == 0) &&
                                              (componentGraphCopy.InDegree(v) == 0)).ToList();
            foreach (var vertex in isolated)
            {
                componentGraphCopy.RemoveVertex(vertex);
            }

            return currentCircuit;
        }

        private static EulerianCycle getExistingEulerian(BidirectionalGraph<int, Edge<int>> componentGraph)
        {
            //find eulerian (componentGraph is eulerian)
            //Hielholzer algorithm
            var componentGraphCopy = componentGraph.Clone();
            //var intersecting = componentGraph.Vertices.Where(v => componentGraph.InDegree(v) > 1).ToHashSet();

            var protoEulerian = getCircuit(componentGraphCopy); //THE COMPONENT GRAPH COPY IS AFFECTED BY THE PROCEDURE
            //the edges in currentCircuit were removed from the copy
            while (componentGraphCopy.EdgeCount != 0)
            {
                var helper = protoEulerian.First(e => (componentGraphCopy.ContainsVertex(e.Target)) && (componentGraphCopy.OutDegree(e.Target) >= 1));
                var vertex = helper.Target;
                var firstEdge = componentGraphCopy.OutEdges(vertex).First();
                
                var nextCircuit = getCircuit(componentGraphCopy, firstEdge);
                
                //not to the end, but behind first
                var firstIndex = protoEulerian.IndexOf(helper);
                var i = firstIndex+1;
                foreach (var e in nextCircuit)
                {
                    protoEulerian.Insert(i, e);
                    i++;
                }
            }
            
            return new EulerianCycle(protoEulerian); //only a single cycle
        }

        private static EulerianCycle getEulerianApproximation(BidirectionalGraph<int, Edge<int>> orig)
        {
           var componentGraph = orig.Clone();
            //solve chinese postman on componentGraph
            //no negative cycles (no negative weights), strongly connected
            //find unbalanced vertices
            var unbalanced =
                componentGraph.Vertices.Where(v => componentGraph.InDegree(v) != componentGraph.OutDegree(v)).ToList();
            var unbalanced_in = unbalanced.Where(v => componentGraph.InDegree(v) > componentGraph.OutDegree(v)).ToList();
            var unbalanced_out = unbalanced.Where(v => componentGraph.InDegree(v) < componentGraph.OutDegree(v)).ToList();

            if (unbalanced_in.Count != unbalanced_out.Count)
            {
                return null; //imbalanced -- no solution here
            }
            
            //find paths from IN to OUT
            var dijkstra = new DijkstraShortestPathAlgorithm<int, Edge<int>>(orig, e => 1d);

            //find minimal cost perfect pairing -- hungarian algorithm
            int[,] costs = new int[unbalanced_in.Count,unbalanced_out.Count]; //spojeni iteho IN s jtym OUT
            for (int i = 0; i < unbalanced_in.Count; i++)
            {
                var vi = unbalanced_in[i];
                dijkstra.Compute(vi);
                for (int j = 0; j < unbalanced_out.Count; j++)
                {
                    var vo = unbalanced_out[j];
                    double cost;

                    cost = dijkstra.Distances[vo];
                    costs[i, j] = (int)cost;
                    //if (floydWarshall.TryGetDistance(vi, vo, out cost))
                    //{
                    //    costs[i, j] = (int)cost;
                    //}
                }
            }
            var hungarian = new HungarianAlgorithm(costs);
            var assignments = hungarian.Run();

            var adjustedGraph = componentGraph.Clone();
            //double edges
            for (int i = 0; i < assignments.Length; i++)
            {
                int j = assignments[i];
                var from = unbalanced_in[i];
                var to = unbalanced_out[j];

                IEnumerable<Edge<int>> pathToAdjust;
                var trypath = orig.ShortestPathsDijkstra(e => 1d, from);
                trypath(to, out pathToAdjust);

                //floydWarshall.TryGetPath(from, to, out pathToAdjust);
                adjustedGraph.AddEdgeRange(pathToAdjust);
            }
            var eulerian = getExistingEulerian(adjustedGraph);

            return eulerian;
        }

        private static ITraversable calculateEulerian(BidirectionalGraph<int, Edge<int>> componentSubgraph)
        {
            ITraversable eulerian;

            var test = new StronglyConnectedComponentsAlgorithm<int, Edge<int>>(componentSubgraph);
            test.Compute();

            if (isEulerian(componentSubgraph))
            {
                eulerian = getExistingEulerian(componentSubgraph);
            }
            else
            {
                eulerian = getEulerianApproximation(componentSubgraph);
            }

            if (eulerian == null)
            {
                //imbalanced ins and outs -- what now???
                eulerian = new NonEulerianComponent(componentSubgraph);
            }
            
            return eulerian;
        }

        private static List<int> eulerianTraversal(ITraversable eulerian, int from, int to)
        {
            return eulerian.GetPassage(from, to);
        }

        private static List<Edge<CondensationNode>> makeEdgesFromGateCombinations(int bigComponent, List<int> gatesIn, 
            List<int> gatesOut, Dictionary<int, Dictionary<int, CondensationNode>> bigDictIn, 
            Dictionary<int, Dictionary<int, CondensationNode>> bigDictOut, ITraversable eulerian
            )
        {
            var edges = new List<Edge<CondensationNode>>();
            foreach (var gateIn in gatesIn)
                {
                    foreach (var gateOut in gatesOut)
                    {
                        CondensationNode inNode;
                        CondensationNode outNode;
                        
                        if (bigDictIn[bigComponent] == null)
                        {
                            bigDictIn[bigComponent] = new Dictionary<int, CondensationNode>();
                            inNode = new GateNode(bigComponent, true, gateIn, "IN");
                            bigDictIn[bigComponent][gateIn] = inNode;
                        }
                        else
                        {
                            if (!bigDictIn[bigComponent].ContainsKey(gateIn))
                            {
                                inNode = new GateNode(bigComponent, true, gateIn, "IN");
                                bigDictIn[bigComponent][gateIn] = inNode;
                            }
                            else
                            {
                                inNode = bigDictIn[bigComponent][gateIn];
                            }
                        }
                        
                        if (bigDictOut[bigComponent] == null)
                        {
                            bigDictOut[bigComponent] = new Dictionary<int, CondensationNode>();
                            outNode = new GateNode(bigComponent, true, gateOut, "OUT");
                            bigDictOut[bigComponent][gateOut] = outNode;
                        }
                        else
                        {
                            if (!bigDictOut[bigComponent].ContainsKey(gateOut))
                            {
                                outNode = new GateNode(bigComponent, true, gateOut, "OUT");
                                bigDictOut[bigComponent][gateOut] = outNode;
                            }
                            else
                            {
                                outNode = bigDictOut[bigComponent][gateOut];
                            }
                        }
                        
                        var walkthrough = eulerianTraversal(eulerian, gateIn, gateOut); //unique for every gateIN->gateOUT
                        var walkthroughNode = new WalkthroughNode(bigComponent, true, gateIn, gateOut, walkthrough);
                        
                        edges.Add(new Edge<CondensationNode>(inNode, walkthroughNode));
                        edges.Add(new Edge<CondensationNode>(walkthroughNode, outNode));
                    }
                }

            return edges;
        }


        private static Condensation construct(DeBruijnGraph dbgraph, int k, StronglyConnectedComponentsAlgorithm<int, Edge<int>> stronglyConnected = null, Dictionary<int, ITraversable>cachedEulerians = null)
        {
            var graph = dbgraph.Graph;

            if (stronglyConnected == null)
            {
                stronglyConnected = new StronglyConnectedComponentsAlgorithm<int, Edge<int>>(graph);
                stronglyConnected.Compute();
            }

            // output from tarjan gives topological sort...
            var vertexToComponent = stronglyConnected.Components.ToDictionary(
                pair => pair.Key, pair => pair.Value
            );
            var componentToMembers = vertexToComponent.GroupBy(r => r.Value
            ).ToDictionary(t => t.Key, t => t.Select(r => r.Key).ToArray());

            var bigComponents = new HashSet<int>();

            if (cachedEulerians == null)
            {
                cachedEulerians = new Dictionary<int, ITraversable>();
            
                for (int component = 0; component < stronglyConnected.ComponentCount; component++)
                {
                    int size = componentToMembers[component].Length;
                    if (size > 1)
                    {
                        bigComponents.Add(component);

                        var eulerian = calculateEulerian(stronglyConnected.Graphs[component]);
                        cachedEulerians[component] = eulerian;
                    }
                }
            }
            
            var bigDictIn = new Dictionary<int, Dictionary<int, CondensationNode>>(); //component --> dict[ vertex entryID --> InGateNode ]
            var bigDictOut = new Dictionary<int, Dictionary<int, CondensationNode>>(); //component --> dict[ vertex entryID --> OutGateNode ]
            var smallDict = new Dictionary<int, CondensationNode>();

            var x = getBigPassEdges(
                graph, bigComponents, vertexToComponent, componentToMembers, cachedEulerians, stronglyConnected.ComponentCount, dbgraph, bigDictIn, bigDictOut, smallDict);
            var edges = x.Item1;
            var isolated = x.Item2.ToHashSet();
            
            //add bigcomponents on edges
            //add isolated vertices
            
            var condensedGraph = new BidirectionalGraph<CondensationNode, Edge<CondensationNode>>();
            condensedGraph.AddVerticesAndEdgeRange(edges);
            condensedGraph.AddVertexRange(isolated);

            List<CondensationNode> toposort = new List<CondensationNode>();
            foreach (var vertex in isolated)
            {
                toposort.Add(vertex); //do not affect toposort
            }

            var sources = edges.Select(e => e.Source
            ).OrderByDescending(s => s.ID).ToList();
            var help = sources.ToHashSet();
            var targets = edges.Select(e => e.Target).Where(v => !help.Contains(v)).ToList();
            toposort.AddRange(sources);
            toposort.AddRange(targets);

            var condensation = new Condensation(stronglyConnected, bigComponents, vertexToComponent,
                toposort, componentToMembers, condensedGraph, dbgraph, k, smallDict, bigDictIn, bigDictOut
            );
            
            return condensation;
        }

        private static (List<Edge<CondensationNode>>, List<CondensationNode>) getBigPassEdges(BidirectionalGraph<int, Edge<int>> graph, HashSet<int> bigComponents, 
            Dictionary<int,int> vertexToComponent, Dictionary<int,int[]> componentToMembers, Dictionary<int, ITraversable> cachedEulerians, int componentCount, DeBruijnGraph deBruijnGraph,
       Dictionary<int, Dictionary<int, CondensationNode>> bigDictIn,
       Dictionary<int, Dictionary<int, CondensationNode>> bigDictOut, //component --> dict[ vertex entryID --> OutGateNode ]
       Dictionary<int, CondensationNode> smallDict, int n=3
            )
        {
            var edges = new List<Edge<CondensationNode>>();
            var isolatedVertices = new List<CondensationNode>();

            var outGates = new Dictionary<int, List<int>>();
            //var inGates = new Dictionary<int, List<int>>();

            foreach (var bigComponent in bigComponents)
            {
                var incomingEdges = componentToMembers[bigComponent].SelectMany(v => graph.InEdges(v)).Where(e => vertexToComponent[e.Source] != bigComponent);
                var outgoingEdges = componentToMembers[bigComponent].SelectMany(v => graph.OutEdges(v)).Where(e => vertexToComponent[e.Target] != bigComponent);

                var gatesIn = incomingEdges.Select(e => e.Target).Distinct().Take(n).ToList();
                //inGates[bigComponent] = gatesIn;
                var gatesOut = outgoingEdges.Select(e => e.Source).Distinct().Take(n).ToList();
                outGates[bigComponent] = gatesOut;

                incomingEdges = null;
                outgoingEdges = null;

                bigDictIn[bigComponent] = null;
                bigDictOut[bigComponent] = null;

                //for each combination of a gate IN and gate OUT
                if (gatesIn.Count > 0 && gatesOut.Count > 0)
                {
                    var eulerian = cachedEulerians[bigComponent]; // precalculated
                    var newEdges = makeEdgesFromGateCombinations(bigComponent, gatesIn, gatesOut, bigDictIn, bigDictOut,
                        eulerian);
                    edges.AddRange(newEdges);
                }
                else
                {
                    if (gatesIn.Count > 0)
                    {
                        foreach (var gateIn in gatesIn)
                        {
                            CondensationNode inNode;
                            if (bigDictIn[bigComponent] == null)
                            {
                                bigDictIn[bigComponent] = new Dictionary<int, CondensationNode>();
                            }
                            inNode = new GateNode(bigComponent, true, gateIn, "IN");
                            bigDictIn[bigComponent][gateIn] = inNode;

                            var cycleNode = new WholeCycleNode(bigComponent, true, gateIn,
                                cachedEulerians[bigComponent].GetCyclicTraversal(gateIn, false));
                            edges.Add(new Edge<CondensationNode>(inNode, cycleNode));
                        }

                        continue; //gates out is empty
                    }

                    if (gatesOut.Count > 0)
                    {
                        foreach (var gateOut in gatesOut)
                        {
                            CondensationNode outNode;
                            if (bigDictOut[bigComponent] == null)
                            {
                                bigDictOut[bigComponent] = new Dictionary<int, CondensationNode>();
                            }
                            outNode = new GateNode(bigComponent, true, gateOut, "OUT");
                            bigDictOut[bigComponent][gateOut] = outNode;

                            var traversal = cachedEulerians[bigComponent].GetCyclicTraversal(gateOut, true);
                            
                            //traversal = traversal.GetRange(1, traversal.Count - 1);
                            
                            var cycleNode = new WholeCycleNode(bigComponent, true, gateOut,
                                traversal
                                );
                            
                            edges.Add(new Edge<CondensationNode>(cycleNode, outNode));
                        }

                        continue;
                    }
                    
                    //component is isolated
                    //select the node with lowest count -- set in on edge
                    var minCount = componentToMembers[bigComponent].Min(v => deBruijnGraph.MinCounts[v]);
                    var entry = componentToMembers[bigComponent].First(v => deBruijnGraph.MinCounts[v] == minCount);
                    //make cycle from the node
                    var pass = new WholeCycleNode(bigComponent, true, entry,
                        cachedEulerians[bigComponent].GetCyclicTraversal(entry, false));
                    //add to some collection??
                    isolatedVertices.Add(pass);
                }
            }
            
            foreach (var bigComponent in bigComponents)
            {
                var outgoingEdges = componentToMembers[bigComponent].SelectMany(v => graph.OutEdges(v)).Where(
                    e => (vertexToComponent[e.Target] != bigComponent) && (outGates[bigComponent].Contains(e.Source))
                    );
                foreach (var outgoingEdge in outgoingEdges)
                {
                    var fromNode = bigDictOut[bigComponent][outgoingEdge.Source];
                    var toComponent = vertexToComponent[outgoingEdge.Target];
                    CondensationNode toNode;
                    if (bigComponents.Contains(toComponent))
                    {
                        if (bigDictIn[toComponent].ContainsKey(outgoingEdge.Target))
                            toNode = bigDictIn[toComponent][outgoingEdge.Target];
                        else
                            continue;
                    }
                    else
                    {
                        if (!smallDict.ContainsKey(toComponent))
                        {
                            toNode = new CondensationNode(toComponent, false, componentToMembers[toComponent].ToList());
                            smallDict.Add(toComponent, toNode);
                        }
                        else
                        {
                            toNode = smallDict[toComponent];
                        }
                    }
                    edges.Add(new Edge<CondensationNode>(fromNode, toNode));
                }
            }
            
            foreach (var singleNodeComponent in Enumerable.Range(0, componentCount).Except(bigComponents))
            {
                int member = componentToMembers[singleNodeComponent][0];
                CondensationNode currentNode;
                
                if (!smallDict.ContainsKey(singleNodeComponent))
                {
                    currentNode = new CondensationNode(singleNodeComponent, false, new List<int>() {member});
                    smallDict.Add(singleNodeComponent, currentNode);
                }
                else
                {
                    currentNode = smallDict[singleNodeComponent];
                }
                
                if (!graph.OutEdges(member).Any() && !graph.InEdges(member).Any())
                    isolatedVertices.Add(currentNode);
                
                foreach (var outEdge in graph.OutEdges(member))
                {
                    var inGraphSuccessor = outEdge.Target;
                    var inCondSuccessor = vertexToComponent[inGraphSuccessor];
                    if (inCondSuccessor == singleNodeComponent)
                    {
                        //should not be of interest -- self-loop, skip...
                        continue;
                    }

                    if (bigComponents.Contains(inCondSuccessor))
                    {
                        //najit IN gatu tehle big componenty
                        if (bigDictIn[inCondSuccessor].ContainsKey(inGraphSuccessor))
                        {
                            var gate = bigDictIn[inCondSuccessor][inGraphSuccessor];
                            edges.Add(new Edge<CondensationNode>(currentNode, gate));
                        }
                        continue;
                    }
                    
                    //successor is a small component
                    CondensationNode successorNode;
                    if (smallDict.ContainsKey(inCondSuccessor))
                    {
                        //node was already added
                        successorNode = smallDict[inCondSuccessor];
                    }
                    else
                    {
                        successorNode = new CondensationNode(inCondSuccessor, false, componentToMembers[inCondSuccessor].ToList());
                        smallDict.Add(inCondSuccessor, successorNode);
                    }
                    edges.Add(new Edge<CondensationNode>(currentNode, successorNode));
                }
            }
            
            return (edges, isolatedVertices);
        }

        public static Condensation Construct(DeBruijnGraph dbgraph, int k)
        {
            var condensation = construct(dbgraph, k);
            return condensation;
        }
    }
}

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace Aligner
{
    //taken from https://gist.github.com/alexandrnikitin/e4176d6b472b39155a7e0e5d68264e65
    [DebuggerDisplay("Value = {Value}, TransitionCount = {_transitionsDictionary.Count}")]
    internal class AhoCorasickTreeNode
    {
        public char Value { get; private set; }
        public AhoCorasickTreeNode Failure { get; set; }

        private readonly List<string> _results;
        private readonly Dictionary<char, AhoCorasickTreeNode> _transitionsDictionary;
        private readonly AhoCorasickTreeNode _parent;

        public IEnumerable<string> Results { get { return _results; } }
        public AhoCorasickTreeNode ParentFailure { get { return _parent == null ? null : _parent.Failure; } }
        public IEnumerable<AhoCorasickTreeNode> Transitions { get { return _transitionsDictionary.Values; } }

        public AhoCorasickTreeNode() : this(null, ' ')
        {
        }

        private AhoCorasickTreeNode(AhoCorasickTreeNode parent, char value)
        {
            Value = value;
            _parent = parent;

            _results = new List<string>();
            _transitionsDictionary = new Dictionary<char, AhoCorasickTreeNode>();
        }

        public void AddResult(string result)
        {
            if (!_results.Contains(result))
            {
                _results.Add(result);
            }
        }

        public void AddResults(IEnumerable<string> results)
        {
            foreach (var result in results)
            {
                AddResult(result);
            }
        }

        public AhoCorasickTreeNode AddTransition(char c)
        {
            var node = new AhoCorasickTreeNode(this, c);
            _transitionsDictionary.Add(node.Value, node);

            return node;
        }

        public AhoCorasickTreeNode GetTransition(char c)
        {
            return _transitionsDictionary.ContainsKey(c)
                       ? _transitionsDictionary[c]
                       : null;
        }

        public bool ContainsTransition(char c)
        {
            return _transitionsDictionary.ContainsKey(c);
        }
    }
    
    public class AhoCorasickTree
    {
        internal AhoCorasickTreeNode Root { get; set; }

        public AhoCorasickTree(IEnumerable<string> keywords)
        {
            Root = new AhoCorasickTreeNode();

            if (keywords != null)
            {
                foreach (var p in keywords)
                {
                    addPatternToTree(p);
                }

                SetFailureNodes();
            }
        }

        public bool Contains(string text, out int startPosition)
        {
            var result = contains(text, out startPosition);
            return result;
        }

        private bool contains(string text, out int startPosition)
        {
            var pointer = Root;
            
            //foreach (var c in text)
            for (int i = 0; i < text.Length; i++)
            {
                var c = text[i];
                startPosition = i;
                var transition = GetTransition(c, ref pointer);

                if (transition != null)
                    pointer = transition;

                if (pointer.Results.Any())
                    return true;
            }
            startPosition = -1;
            return false;
        }

        public IEnumerable<Tuple<string,int>> FindAll(string text)
        {
            var pointer = Root;

            for (int i = 0; i < text.Length; i++)
            {
                var c = text[i];
                var transition = GetTransition(c, ref pointer);

                if (transition != null)
                    pointer = transition;

                foreach (var result in pointer.Results)
                    yield return new Tuple<string, int>(result, i - result.Length + 1);
            }
        }

        private AhoCorasickTreeNode GetTransition(char c, ref AhoCorasickTreeNode pointer)
        {
            AhoCorasickTreeNode transition = null;
            while (transition == null)
            {
                transition = pointer.GetTransition(c);

                if (pointer == Root)
                    break;

                if (transition == null)
                    pointer = pointer.Failure;
            }
            return transition;
        }

        private void SetFailureNodes()
        {
            var nodes = FailToRootNode();
            FailUsingBFS(nodes);
            Root.Failure = Root;
        }

        private void resetFailures()
        {
            var node = Root;
            var nextLevel = new Queue<AhoCorasickTreeNode>();
            foreach (var item in node.Transitions)
            {
                nextLevel.Enqueue(item);
            }
            
            while (nextLevel.Count > 0)
            {
                node.Failure = null;
                foreach (var item in node.Transitions)
                {
                    nextLevel.Enqueue(item);
                }

                node = nextLevel.Dequeue();
            }
        }

        public void AddPatternToTree(string pattern)
        {
            resetFailures();
            addPatternToTree(pattern);
            SetFailureNodes();
        }
        
        public void AddPatternRangeToTree(IEnumerable<string> patterns)
        {
            resetFailures();
            foreach (var pattern in patterns)
            {
                addPatternToTree(pattern);
            }
            SetFailureNodes();
        }

        private void addPatternToTree(string pattern)
        {
            var node = Root;
            foreach (var c in pattern)
            {
                node = node.GetTransition(c)
                       ?? node.AddTransition(c);
            }
            node.AddResult(pattern);
        }

        private List<AhoCorasickTreeNode> FailToRootNode()
        {
            var nodes = new List<AhoCorasickTreeNode>();
            foreach (var node in Root.Transitions)
            {
                node.Failure = Root;
                nodes.AddRange(node.Transitions);
            }
            return nodes;
        }

        private void FailUsingBFS(List<AhoCorasickTreeNode> nodes)
        {
            while (nodes.Count != 0)
            {
                var newNodes = new List<AhoCorasickTreeNode>();
                foreach (var node in nodes)
                {
                    var failure = node.ParentFailure;
                    var value = node.Value;

                    while (failure != null && !failure.ContainsTransition(value))
                    {
                        failure = failure.Failure;
                    }

                    if (failure == null)
                    {
                        node.Failure = Root;
                    }
                    else
                    {
                        node.Failure = failure.GetTransition(value);
                        node.AddResults(node.Failure.Results);
                    }

                    newNodes.AddRange(node.Transitions);
                }
                nodes = newNodes;
            }
        }
    }
}
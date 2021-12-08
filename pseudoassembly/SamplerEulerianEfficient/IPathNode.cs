using System;
using System.Collections.Generic;
using System.Linq;

namespace Sampler2
{
    public interface IPathNode
    {
        int PossibilitiesCount { get; }
        int RepresentedLength { get; }
        int[,] GetRepresentedAlignment();
        int[] GetNodes();
        int[] GetPosition(int index);
        IEnumerable<IEnumerable<int>> GetRepresentedPaths();
    }

    public class SingleItemNode : IPathNode
    {
        private int[,] repre;
        public SingleItemNode(int representedNode)
        {
            RepresentedNode = representedNode;
        }
        public int RepresentedNode { get; }
        public int PossibilitiesCount
        {
            get { return 1; }
        }
        
        public int RepresentedLength
        {
            get { return 1; }
        }

        public IEnumerable<IEnumerable<int>> GetRepresentedPaths()
        {
            return new List<IEnumerable<int>>() {IEnumerableExt.SingleItemAsEnumerable(RepresentedNode)};
        }

        public int[,] GetRepresentedAlignment() //possibilities X length
        {
            if (repre == null)
            {
                var val = new int[1, 1];
                val[0, 0] = RepresentedNode;
                repre = val;
                return val;
            }
            return repre;
        }

        public int[] GetNodes()
        {
            return new[] {RepresentedNode};
        }

        public int[] GetPosition(int index)
        {
            if ((index < 0) || (index >= RepresentedLength))
                throw new IndexOutOfRangeException();

            return new int[] {RepresentedNode};
        }
    }

    public class MultiPossibilityNode : IPathNode
    {
        private IEnumerable<IEnumerable<int>> represented;
        private int[,] representedBox;
        public MultiPossibilityNode(IEnumerable<IEnumerable<int>> represented, bool representsBig)
        {
            this.represented = represented;
            PossibilitiesCount = represented.Count();
            RepresentedLength = represented.First().Count();
            RepresentsBig = representsBig;

            representedBox = new int[PossibilitiesCount, RepresentedLength];
            int i = 0;
            foreach (var possiblePath in represented)
            {
                int j = 0;
                foreach (var item in possiblePath)
                {
                    representedBox[i, j] = item;
                    j++;
                }
                i++;
            }
        }

        public IEnumerable<IEnumerable<int>> GetRepresentedPaths()
        {
            return represented;
        }

        public int PossibilitiesCount { get; }
        public bool RepresentsBig { get; }
        public int RepresentedLength { get; }
        public int[,] GetRepresentedAlignment()
        {
            return representedBox;
        }

        public int[] GetNodes()
        {
            var hs = represented.SelectMany(p => p).ToHashSet().ToArray();
            return hs;
        }
        
        public int[] GetPosition(int index)
        {
            if ((index < 0) || (index >= RepresentedLength))
                throw new IndexOutOfRangeException();

            return Enumerable.Range(0, PossibilitiesCount).Select(i => representedBox[i, index]).ToArray();
        }
    }
}
//<----- Suffix Node ----->
class SuffixNode
{
  public:
    //Depth
    //Start, end specify the edge by which the node is connected to it's parent
    int depth, start;
    int *end;

    SuffixNode **children;  //Children of node
                            //Pointer to pointer as we don't know size yet since alphabet may be incredibly large
    SuffixNode *parent;     //Parent of node
    SuffixNode *suffixLink; //Suffix Link pointer

    int suffixIndex = -1;

    //Constructor
    SuffixNode(int newStart, int *newEnd, SuffixNode *root, int childrenSize)
    {
        start = newStart;
        end = newEnd;
        suffixLink = root;
        suffixIndex = -1;

        this->children = new SuffixNode *[childrenSize];
        for (int i = 0; i < childrenSize; i++)
        {
            this->children[i] = nullptr;
        }
    }
};

//<----- Suffix Tree ----->
template <typename ContainerType, typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0), typename MatchFnTy = std::function<bool(Ty, Ty)>>
class SuffixTree
{
  private:
    SuffixNode *lastNewNode = nullptr;
    SuffixNode *activeNode = nullptr;

    int activeEdge = -1;
    int activeLength = 0;

    int remainingSuffixCount = 0;
    int leafEnd = -1;
    int *rootEnd = nullptr;
    int *splitEnd = nullptr;
    int seqSize = -1;
    int size1 = -1;

    SuffixNode *root;

    ContainerType seq;

    MatchFnTy match;

  public:
    //Builds the suffix tree and returns the root node
    void buildTree(ContainerType &newSeq, int newSize1, MatchFnTy matchFn)
    {

        seq = newSeq;
        seqSize = seq.size();
        size1 = newSize1;

        int *rootEnd = new int;
        *rootEnd = -1;
        root = new SuffixNode(-1, rootEnd, nullptr, seqSize);
        SuffixNode *cn = root;
        root->suffixLink = root;
        SuffixNode *needsSuffixLink = nullptr;

        match = matchFn;

        activeNode = root;

        int lastRule = 0;

        //Begin construction
        for (int i = 0; i < seqSize; i++)
        {
            extendSuffixTree(i);
        }
        int labelHeight = 0;
        setSuffixIndexByDFS(root, labelHeight);
    }

    void setSuffixIndexByDFS(SuffixNode *n, int labelHeight)
    {
        if (n == nullptr)
        {
            return;
        }

        int leaf = 1;
        for (int i = 0; i < seqSize; i++)
        {
            if (n->children[i] != nullptr)
            {
                leaf = 0;
                setSuffixIndexByDFS(n->children[i], labelHeight + edgeLength(n->children[i]));
            }
        }

        if (leaf == 1)
        {
            for (int i = n->start; i <= *(n->end); i++)
            {
                if (seq[i] == '#')
                {
                    n->end = new int;
                    *(n->end) = i;
                }
            }
            n->suffixIndex = seqSize - labelHeight;
        }
    }

    int edgeLength(SuffixNode *n)
    {
        if (n == root)
        {
            return 0;
        }
        return *(n->end) - (n->start) + 1;
    }

    int walkDown(SuffixNode *currNode)
    {
        //If activeLength is greater than current edge length,
        //set next internal node as activeNode and adjust activeEdge
        //and activeLength accordingly
        if (activeLength >= edgeLength(currNode))
        {
            activeEdge += edgeLength(currNode);
            activeLength -= edgeLength(currNode);
            activeNode = currNode;
            return 1;
        }
        return 0;
    }

    void extendSuffixTree(int pos)
    {
        leafEnd = pos + 1;
        remainingSuffixCount++;
        lastNewNode = nullptr;

        //Add all suffixes (yet to be added) one by one
        while (remainingSuffixCount > 0)
        {
            if (activeLength == 0)
            {
                activeEdge = pos;
            }

            //std::cout << "Active edge: " << activeEdge << std::endl;
            //std::cout << "Active length: " << activeLength << std::endl;

            int i = 0;
            bool found = false;
            do
            {
                if (match(seq[activeEdge], seq[i]))
                {
                    activeEdge = i;
                    found = true;
                    break;
                }
                i++;
            } while (activeNode->children[activeEdge] != nullptr);

            if (!found)
            {
                activeEdge = pos;
            }

            //There is no outgoing edge starting with activeEdge from activeNode
            if (activeNode->children[activeEdge] == nullptr)
            {
                activeNode->children[activeEdge] = new SuffixNode(pos, &leafEnd, root, seqSize);

                if (lastNewNode != nullptr)
                {
                    lastNewNode->suffixLink = activeNode;
                    lastNewNode = nullptr;
                }
            }
            //There is an outgoing edge starting with activeEdge from activeNode
            else
            {

                //Get the next node at the end of edge starting with activeEdge
                SuffixNode *next = activeNode->children[activeEdge];
                if (walkDown(next))
                {
                    //Start from next node (the new activeNode)
                    continue;
                }

                //Extension rule 3
                //Current character being processed is already on the edge
                if (match(seq[next->start + activeEdge], seq[pos]))
                {

                    //If a newly created node waiting for it's
                    //suffix link to be set, then set suffix link
                    //to current active node
                    if (lastNewNode != nullptr && activeNode != root)
                    {
                        lastNewNode->suffixLink = activeNode;
                        lastNewNode = nullptr;
                    }

                    activeLength++;

                    //Stop all processing in this phase
                    //Move onto next phase
                    break;
                }

                //Extension rule 2
                //activePoint is in the middle of the edge being traversed
                //and current character being processed is not on the edge
                //So add a new node and a new leaf edge out of that node
                int *splitEnd = new int;
                *splitEnd = next->start + activeLength - 1;

                //New internal node
                SuffixNode *split = new SuffixNode(next->start, splitEnd, root, seqSize);
                activeNode->children[activeEdge] = split;

                //New leaf coming out of new internal node
                split->children[activeEdge] = new SuffixNode(pos, &leafEnd, root, seqSize);
                next->start += activeLength;
                split->children[activeEdge] = next;

                //We have a new internal node
                //If there is any internal node created in last extensions
                //which are still waiting for it's suffix link reset then reset
                if (lastNewNode != nullptr)
                {
                    //suffixLink of lastNewNode points to current newly created node
                    lastNewNode->suffixLink = split;
                }

                lastNewNode = split;
            }

            //One suffix got added in tree
            //Decrement the count of suffixes yet to be added
            remainingSuffixCount--;
            if (activeNode == root && activeLength > 0)
            {
                activeLength--;
                activeEdge = pos - remainingSuffixCount + 1;
            }
            else if (activeNode != root)
            {
                activeNode = activeNode->suffixLink;
            }
        }
    }

    void getLCS()
    {
        ContainerType ret;
        int maxHeight = 0;
        int substringStartIndex = 0;

        traverse(root, 0, &maxHeight, &substringStartIndex);

        std::cout << "printing the LCS: " << std::endl;
        for (int k = 0; k < maxHeight; k++)
        {
            std::cout << seq[substringStartIndex + k] << std::endl;
        }

        return;
    }

    int traverse(SuffixNode *n, int labelHeight, int *maxHeight, int *substringStartIndex)
    {
        if (n == nullptr)
        {
            return 0;
        }

        int i = 0;
        int ret = -1;

        if (n->suffixIndex < 0) //If it is an internal node
        {
            for (i = 0; i < seqSize; i++)
            {
                if (n->children[i] != nullptr)
                {
                    ret = traverse(n->children[i], labelHeight + edgeLength(n->children[i]), maxHeight, substringStartIndex);
                    if (n->suffixIndex == -1)
                    {
                        n->suffixIndex = ret;
                    }
                    else if ((n->suffixIndex == -2 && ret == -3) || (n->suffixIndex == -3 && ret == -2) || n->suffixIndex == -4)
                    {
                        n->suffixIndex = -4;
                        if (*maxHeight < labelHeight)
                        {
                            *maxHeight = labelHeight;
                            *substringStartIndex = *(n->end) - labelHeight + 1;
                        }
                    }
                }
            }
        }
        else if (n->suffixIndex > -1 && n->suffixIndex < size1)
            return -2;
        else if (n->suffixIndex >= size1)
            return -3;
        return n->suffixIndex;
    }
};

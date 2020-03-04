
//<----- Suffix Node ----->
class SuffixNode
{
  public:
    //Start, end specify the edge by which the node is connected to it's parent
    int start;
    int secondStart;
    int *end;

    SuffixNode **children;  //Children of node
                            //Pointer to pointer as we don't know size yet since alphabet may be incredibly large
    SuffixNode *suffixLink; //Suffix Link pointer

    SuffixNode *parent; //Parent of node
                        //Used for MUM finding

    int suffixIndex;

    bool leaf;

    //Constructor
    SuffixNode(int newStart, int *newEnd, SuffixNode *root, int childrenSize)
    {
        start = newStart;
        secondStart = newStart;
        end = newEnd;
        suffixLink = root;
        suffixIndex = -1;
        leaf = false;

        this->children = new SuffixNode *[childrenSize];
        for (int i = 0; i < childrenSize; i++)
        {
            this->children[i] = nullptr;
        }
    }

    SuffixNode *operator=(const SuffixNode &b)
    {
        this->children = b.children;
        this->start = b.start;
        this->end = b.end;
        this->suffixLink = b.suffixLink;
        this->suffixIndex = b.suffixIndex;
        this->secondStart = b.secondStart;
        this->leaf = b.leaf;
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
    int seqSize = -1; //Length of input string
    int size1 = 0;    //Length of 1st string

    SuffixNode *root;

    ContainerType seq;
    std::vector<Ty> alphabet;

    MatchFnTy match;

    using BaseType = SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>;

  public:
    //Builds the generalised suffix tree
    //Generalised suffix tree used for getting the LCS
    void buildTree(ContainerType &Seq1, ContainerType &Seq2, MatchFnTy matchFn)
    {
        //seq = Seq1 + "#" + Seq2 + "$";
        seq = Seq1;
        Ty hash = Blank;
        Ty dollar = Ty(0);
        seq.push_back(hash);

        for (int i = 0; i < Seq2.size(); i++)
        {
            seq.push_back(Seq2[i]);
        }

        seq.push_back(dollar);

        seqSize = seq.size();
        size1 = Seq1.size();

        int *rootEnd = new int;
        *rootEnd = -1;
        root = new SuffixNode(-1, rootEnd, nullptr, seqSize);
        root->suffixLink = root;

        match = matchFn;

        activeNode = root;

        //Begin construction
        for (int i = 0; i < seqSize; i++)
        {
            extendSuffixTree(i);
        }

        //print(root);
        int labelHeight = 0;
        setSuffixIndexByDFS(root, labelHeight);
    }

    //Builds the suffix tree
    void buildTree(ContainerType &Seq1, MatchFnTy matchFn)
    {
        seq = Seq1;
        Ty dollar = Ty(0);
        seq.push_back(dollar);

        seqSize = seq.size();
        size1 = Seq1.size();

        int *rootEnd = new int;
        *rootEnd = -1;
        root = new SuffixNode(-1, rootEnd, nullptr, seqSize);
        root->suffixLink = root;

        match = matchFn;

        activeNode = root;

        //Begin construction
        for (int i = 0; i < seqSize; i++)
        {
            extendSuffixTree(i);
        }

        //print(root);
        int labelHeight = 0;
        setSuffixIndexByDFS(root, labelHeight);
    }

    void print(SuffixNode *node)
    {
        for (int i = 0; i < alphabet.size(); i++)
        {
            if (node->children[i] != nullptr)
            {
                std::string subst = seq.substr(node->children[i]->start, *(node->children[i]->end) - node->children[i]->start + 1);
                std::cout << seq[node->children[i]->start] << std::endl;
                std::cout << subst << std::endl;
            }
        }
        std::cout << std::endl;
        for (int i = 0; i < alphabet.size(); i++)
        {
            if (node->children[i] != nullptr)
            {
                print(node->children[i]);
            }
        }
    }

    void setSuffixIndexByDFS(SuffixNode *n, int labelHeight)
    {
        if (n == nullptr)
        {
            return;
        }

        int leaf = 1;
        for (int i = 0; i < alphabet.size(); i++)
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
                if (match(seq[i], Blank))
                {
                    n->end = new int;
                    *(n->end) = i;
                }
            }
            n->suffixIndex = seqSize - labelHeight + 1;
            n->leaf = true;
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

            /*std::cout << "" << std::endl;
            std::cout << "" << std::endl;
            std::cout << "Remaining suffix count: " << remainingSuffixCount << std::endl;
            std::cout << "LeafEnd: " << leafEnd << std::endl;
            std::cout << "Reading character: " << seq[pos] << std::endl;*/
            if (activeLength == 0)
            {
                activeEdge = pos;
            }

            auto result1 = find(std::begin(alphabet), std::end(alphabet), seq[pos - activeLength]);

            if (result1 == std::end(alphabet))
            {
                //std::cout << "Can't find " << seq[pos - activeLength] << " in alphabet, adding now" << std::endl;
                alphabet.push_back(seq[pos - activeLength]);
                //std::cout << "Alphabet size: " << alphabet.size() << std::endl;
                activeEdge = alphabet.size() - 1;
            }
            else
            {
                //std::cout << "Found " << seq[pos - activeLength] << " in alphabet" << std::endl;
                // Get index of element from iterator
                activeEdge = std::distance(alphabet.begin(), result1);
            }

            //std::cout << "Active edge: " << activeEdge << std::endl;
            //std::cout << "Active length: " << activeLength << std::endl;

            //There is no outgoing edge starting with activeEdge from activeNode
            if (activeNode->children[activeEdge] == nullptr)
            {
                activeNode->children[activeEdge] = new SuffixNode(pos, &leafEnd, root, seqSize);

                if (lastNewNode != nullptr)
                {
                    //std::cout << "LAST NEW NODE SUFFIX LINK IS SET" << std::endl;
                    lastNewNode->suffixLink = activeNode;
                    lastNewNode = nullptr;
                }
            }
            //There is an outgoing edge starting with activeEdge from activeNode
            else
            {
                //std::cout << "Found " << alphabet[activeEdge] << " in activeNode" << std::endl;

                //std::cout << "there is an outgoing edge" << std::endl;
                //Get the next node at the end of edge starting with activeEdge
                SuffixNode *next = activeNode->children[activeEdge];
                //SuffixNode *next = new SuffixNode(activeEdge, &leafEnd, root, seqSize);
                //*next = *activeNode->children[activeEdge];
                //std::cout << "***** impotant ********" << std::endl;
                //std::cout << *next->end << " " << *activeNode->children[activeEdge]->end << std::endl;
                //std::cout << "active edge is when initialising next: " << activeEdge << std::endl;
                if (walkDown(next))
                {
                    //Start from next node (the new activeNode)
                    //std::cout << "(" << activeNode->start << "," << activeEdge << "," << activeLength << ")" << std::endl;
                    continue;
                }

                //std::cout << "no need for walkdown" << std::endl;
                //Extension rule 3
                //Current character being processed is already on the edge
                //std::cout << "next->start: " << next->start << " activeLength: " << activeLength << " pos: " << pos << std::endl;
                //std::cout << "Attempting to match: " << seq[next->start + activeLength] << " and " << seq[pos] << std::endl;
                if (match(seq[next->start + activeLength], seq[pos]))
                {
                    //std::cout << "Rule 3" << std::endl;
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
                    //std::cout << "(" << activeNode->start << "," << activeEdge << "," << activeLength << ")" << std::endl;

                    break;
                }

                //Extension rule 2
                //activePoint is in the middle of the edge being traversed
                //and current character being processed is not on the edge
                //So add a new node and a new leaf edge out of that node

                auto result1 = find(std::begin(alphabet), std::end(alphabet), seq[pos]);
                int leafEdge = activeEdge;

                if (result1 == std::end(alphabet))
                {
                    //std::cout << "Can't find " << seq[pos] << " in alphabet, adding now" << std::endl;
                    alphabet.push_back(seq[pos]);
                    //std::cout << "Alphabet size: " << alphabet.size() << std::endl;
                    //activeEdge = alphabet.size();
                    leafEdge = alphabet.size() - 1;
                }
                else
                {
                    //std::cout << "Found " << seq[pos] << " in alphabet" << std::endl;
                    // Get index of element from iterator
                    //activeEdge = std::distance(alphabet.begin(), result1);
                    leafEdge = std::distance(alphabet.begin(), result1);
                }
                //std::cout << "Rule 2" << std::endl;
                splitEnd = new int;
                *splitEnd = next->start + activeLength - 1;

                //New internal node
                //std::cout << "creating leaf out of new internal node when active edge is: " << activeEdge << std::endl;

                SuffixNode *split = new SuffixNode(next->start, splitEnd, root, seqSize);
                split->secondStart = pos;
                activeNode->children[activeEdge] = split;

                //New leaf coming out of new internal node
                split->children[leafEdge] = new SuffixNode(pos, &leafEnd, root, seqSize);

                //New edge coming out of internal node
                next->start += activeLength;
                //std::cout << "NEXT START IS: " << next->start << std::endl;
                auto result2 = find(std::begin(alphabet), std::end(alphabet), seq[next->start]);
                int continuingEdge = activeEdge;
                if (result2 == std::end(alphabet))
                {
                    alphabet.push_back(seq[next->start]);
                    continuingEdge = alphabet.size() - 1;
                }
                else
                {
                    // Get index of element from iterator
                    continuingEdge = std::distance(alphabet.begin(), result2);
                }
                //std::cout << "CONTINUING EDGE IS: " << continuingEdge << std::endl;
                split->children[continuingEdge] = new SuffixNode(pos, &leafEnd, root, seqSize);
                split->children[continuingEdge] = next;

                //We have a new internal node
                //If there is any internal node created in last extensions
                //which are still waiting for it's suffix link reset then reset
                if (lastNewNode != nullptr)
                {
                    //suffixLink of lastNewNode points to current newly created node
                    lastNewNode->suffixLink = split;
                }

                lastNewNode = split;
                //std::cout << "LAST NEW NODE IS SET AS SPLIT" << std::endl;
            }

            //One suffix got added in tree
            //Decrement the count of suffixes yet to be added
            remainingSuffixCount--;
            if (activeNode == root && activeLength > 0)
            {
                //std::cout << "active node detected as root" << std::endl;
                activeLength--;
                activeEdge = pos - remainingSuffixCount + 1;

                auto result1 = find(std::begin(alphabet), std::end(alphabet), seq[pos - activeLength]);

                if (result1 == std::end(alphabet))
                {
                    //std::cout << "Can't find " << seq[pos - activeLength] << " in alphabet, adding now" << std::endl;
                    alphabet.push_back(seq[pos - activeLength]);
                    //std::cout << "Alphabet size: " << alphabet.size() << std::endl;
                    activeEdge = alphabet.size() - 1;
                }
                else
                {
                    //std::cout << "Found " << seq[pos - activeLength] << " in alphabet" << std::endl;
                    // Get index of element from iterator
                    activeEdge = std::distance(alphabet.begin(), result1);
                }
            }
            else if (activeNode != root)
            {
                //std::cout << "active node not detected as root" << std::endl;
                activeNode = activeNode->suffixLink;
            }

            //std::cout << "(" << activeNode->start << "," << activeEdge << "," << activeLength << ")" << std::endl;
        }
    }

    //<----- Funcs used for LCS ----->
    AlignedSequence<Ty, Blank> getLCS(ContainerType &Seq1, ContainerType &Seq2, MatchFnTy matchFn)
    {

        buildTree(Seq1, Seq2, matchFn);
        ContainerType ret;
        int maxHeight = 0;
        int substringStartIndex = 0;
        int secondIndex = 0;

        traverse(root, 0, &maxHeight, &substringStartIndex, &secondIndex);

        /*ArrayView<ContainerType> seqSub(seq);
        seqSub.sliceWindow(substringStartIndex, substringStartIndex + maxHeight);
        auto it = seqSub.begin();

        while (it != seqSub.end())
        {
            std::cout << *it;
            it++;
        }
        std::cout << std::endl;

        std::cout << "FOUND ON SEQ 2 INDEX: " << secondIndex << std::endl;*/
        AlignedSequence<Ty, Blank> Result;

        for (int i = 0; i < maxHeight; i++)
        {
            Result.Data.push_back(
                typename BaseType::EntryType(Seq1[substringStartIndex + i], Seq2[secondIndex + i], true));
        }

        return Result;
    }

    void deleteTree()
    {
        freeSuffixTreeByPostOrder(root);
        alphabet.clear();
    }

    void freeSuffixTreeByPostOrder(SuffixNode *n)
    {
        if (n == nullptr)
            return;

        for (int i = 0; i < alphabet.size(); i++)
        {
            if (n->children[i] != nullptr)
            {
                freeSuffixTreeByPostOrder(n->children[i]);
            }
        }
        if (n->suffixIndex == -1)
            delete n->end;
        delete n;
    }

    int traverse(SuffixNode *n, int labelHeight, int *maxHeight, int *substringStartIndex, int *secondIndex)
    {
        if (n == nullptr)
        {
            return 0;
        }

        int i = 0;
        int ret = -1;

        if (n->suffixIndex < 0) //If it is an internal node
        {
            for (i = 0; i < alphabet.size(); i++)
            {
                if (n->children[i] != nullptr)
                {
                    ret = traverse(n->children[i], labelHeight + edgeLength(n->children[i]), maxHeight, substringStartIndex, secondIndex);
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
                            *secondIndex = n->secondStart - labelHeight - size1 - 1;
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

    //<----- Funcs used for pattern matching ----->
    int traverseEdge(ContainerType pattern, int idx, int start, int end)
    {
        int k = 0;
        //Traverse the edge with character by character matching
        for (k = start; k <= end && idx < pattern.size(); k++, idx++)
        {
            if (!match(seq[k], pattern[idx]))
                return -1; // mo match
        }

        if (idx == pattern.size())
            return 1; // match
        return 0;     // more characters yet to match
    }

    int doTraversalToCountLeaf(SuffixNode *n)
    {
        if (n == nullptr)
            return 0;
        if (n->suffixIndex > -1)
        {
            std::cout << "Found at position: " << n->suffixIndex << std::endl;
            return 1;
        }
        int count = 0;
        int i = 0;
        for (i = 0; i < alphabet.size(); i++)
        {
            if (n->children[i] != nullptr)
            {
                count += doTraversalToCountLeaf(n->children[i]);
            }
        }
        return count;
    }

    int countLeaf(SuffixNode *n)
    {
        if (n == nullptr)
            return 0;
        return doTraversalToCountLeaf(n);
    }

    int doTraversal(SuffixNode *n, ContainerType pattern, int idx)
    {
        if (n == nullptr)
        {
            return -1; // no built tree so return no match
        }
        int res = -1;
        //If node n is not root node, then traverse edge
        //from node n's parent to node n.
        if (n->start != -1)
        {
            res = traverseEdge(pattern, idx, n->start, *(n->end));
            if (res == -1) //no match
                return -1;
            if (res == 1) //match
            {
                if (n->suffixIndex > -1)
                    std::cout << "susbtring count: 1 at position: " << n->suffixIndex << std::endl;
                else
                    std::cout << "substring count: " << countLeaf(n) << std::endl;
                return 1;
            }
        }
        //Get the character index to search
        idx = idx + edgeLength(n);
        //If there is an edge from node n going out
        //with current character str[idx], traverse that edge
        auto result1 = find(std::begin(alphabet), std::end(alphabet), pattern[idx]);
        int index = idx;
        if (result1 != std::end(alphabet))
        {
            index = std::distance(alphabet.begin(), result1);
        }

        if (n->children[index] != nullptr)
            return doTraversal(n->children[index], pattern, idx);
        else
            return -1; // no match
    }

    void checkForSubString(ContainerType pattern)
    {
        int res = doTraversal(root, pattern, 0);
        if (res == 1)
            std::cout << "Pattern is a susbtring" << std::endl;
        else
            std::cout << "Pattern is NOT a substring" << std::endl;
    }

    //<----- Funcs for getting MUMs ----->
    std::vector<MUM> getMUMs(ContainerType &Seq1, ContainerType &Seq2, MatchFnTy matchFn)
    {
        //Build tree
        buildTree(Seq1, Seq2, matchFn);

        //Mark nodes
        std::vector<SuffixNode *> markedNodes;
        std::vector<int> mumStarts;
        std::vector<bool> isSuffix;
        unsigned int length = 0;
        markNodes(root, markedNodes, root, mumStarts, length);

        //Unmark nodes with suffix links to marked nodes
        for (int i = 0; i < markedNodes.size(); i++)
        {
            isSuffix.push_back(false);
            for (int j = 0; j < markedNodes.size(); j++)
            {
                //Skip itself
                if (i != j)
                {
                    //Remove nodes with suffix links to marked nodes
                    if (markedNodes[i] == markedNodes[j]->suffixLink)
                    {
                        isSuffix[i] = true;
                        break;
                    }
                }
            }
        }

        std::vector<SuffixNode *> mumNodes;
        std::vector<int> newMumStarts;
        std::vector<MUM> ret;
        for (int i = 0; i < markedNodes.size(); i++)
        {
            if (!isSuffix[i])
            {
                mumNodes.push_back(markedNodes[i]);
                newMumStarts.push_back(mumStarts[i]);
                MUM mum = MUM(mumStarts[i], (markedNodes[i]->secondStart - (*(markedNodes[i]->end) - mumStarts[i] + 1) - size1 - 1), (*(markedNodes[i]->end) - mumStarts[i] + 1));
                ret.push_back(mum);
            }
        }

        //If no MUM available then simply return the longest common substring
        if (ret.size() == 0)
        {
            int maxHeight = 0;
            int substringStartIndex = 0;
            int secondIndex = 0;

            traverse(root, 0, &maxHeight, &substringStartIndex, &secondIndex);

            MUM mum = MUM(substringStartIndex, secondIndex, maxHeight);
            ret.push_back(mum);
        }

        return ret;

        /*markedNodes.clear();
        markedNodes = mumNodes;
        mumStarts.clear();
        mumStarts = newMumStarts;

        for (int i = 0; i < markedNodes.size(); i++)
        {
            std::string subst = seq.substr(mumStarts[i], *(markedNodes[i]->end) - mumStarts[i] + 1);
            std::cout << subst << std::endl;
            std::cout << "Found at second position: " << markedNodes[i]->secondStart - (*(markedNodes[i]->end) - mumStarts[i] + 1) - size1 - 1 << std::endl;
        }*/
    }

    //Recursively loop over all nodes and mark them if they have two leaf children
    void markNodes(SuffixNode *n, std::vector<SuffixNode *> &markedNodes, SuffixNode *startNode, std::vector<int> &mumStarts, unsigned int length)
    {
        if (n == nullptr)
        {
            return;
        }

        int noChildren = 0;
        bool isSeq1 = false;
        bool isSeq2 = false;
        for (int i = 0; i < alphabet.size(); i++)
        {
            if (n->children[i] != nullptr)
            {
                if (n->children[i]->leaf)
                {
                    noChildren++;

                    if (n->children[i]->suffixIndex >= 0 && n->children[i]->suffixIndex < size1)
                    {
                        isSeq1 = true;
                    }
                    else if (n->children[i]->suffixIndex >= size1)
                    {
                        isSeq2 = true;
                    }
                }
                else
                {
                    if (n == root)
                    {
                        markNodes(n->children[i], markedNodes, n->children[i], mumStarts, edgeLength(n->children[i]));
                    }
                    else
                    {
                        markNodes(n->children[i], markedNodes, startNode, mumStarts, length + edgeLength(n->children[i]));
                    }
                    noChildren = 3;
                }
            }
        }

        //A POTENTIAL MUM
        if (noChildren == 2 && isSeq1 && isSeq2)
        {
            markedNodes.push_back(n);
            mumStarts.push_back(n->start + edgeLength(n) - length);
        }
    }

    template <class InputIt, class T>
    constexpr InputIt find(InputIt first, InputIt last, const T &value)
    {
        for (; first != last; ++first)
        {
            if (match(*first, value))
            {
                return first;
            }
        }

        return last;
    }
};

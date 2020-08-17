#include <iostream>
#include <fstream>
#include <memory>
//<----- Suffix Node ----->
class SuffixNode
{
public:
    // Start, end specify the edge by which the node is connected to it's parent
    int start;
    int secondStart;
    std::shared_ptr<int> end;

    SuffixNode **children;  // Children of node
                            // Pointer to pointer as we don't know size yet since
                            // alphabet may be incredibly large
    SuffixNode *suffixLink; // Suffix Link pointer

    int suffixIndex;

    bool leaf;

    // Constructor
    SuffixNode(int newStart, std::shared_ptr<int> newEnd, SuffixNode *root, int childrenSize)
    {
        start = newStart;
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

    // Destructor
    ~SuffixNode() { delete[] children; }

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
    std::shared_ptr<int> leafEnd = nullptr;
    std::shared_ptr<int> splitEnd = nullptr;
    std::shared_ptr<int> rootEnd = nullptr;
    int seqSize = -1; // Length of input string
    int size1 = 0;    // Length of 1st string

    SuffixNode *root;

    std::vector<Ty> seq;

    std::vector<Ty> alphabet;

    MatchFnTy match;

    using BaseType = SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>;

    bool seenHash = false;

    bool seenDollar = false;

    bool general = true;

public:
    SuffixTree() {}

    // Builds the generalised suffix tree
    // Generalised suffix tree used for getting the LCS
    void buildTree(ContainerType &Seq1, ContainerType &Seq2, MatchFnTy matchFn)
    {
        // seq = Seq1 + "#" + Seq2 + "$";
        // seq = Seq1;
        // basic file operations
        match = matchFn;
        Ty hash = Blank;
        Ty dollar = Blank;

        // Push Seq1
        for (int i = 0; i < Seq1.size(); i++)
        {
            seq.push_back(Seq1[i]);
        }
        
        // Push Hash
        seq.push_back(hash);

        // Push Seq2
        for (int i = 0; i < Seq2.size(); i++)
        {
            seq.push_back(Seq2[i]);
        }

        // Push Dollar
        seq.push_back(dollar);

        seqSize = seq.size();
        size1 = Seq1.size() + 1;

        int *a = new int;
        std::shared_ptr<int> newPtr(a);
        rootEnd = newPtr;
        *rootEnd = -1;

        int *b = new int;
        std::shared_ptr<int> newPtrb(b);
        leafEnd = newPtrb;
        *leafEnd = -1;
        root = new SuffixNode(-1, rootEnd, nullptr, seqSize);
        root->suffixLink = root;

        activeNode = root;

        general = true;

        for (int i = 0; i < seqSize; i++)
        {
            extendSuffixTree(i);
        }

        int labelHeight = 0;
        setSuffixIndexByDFS(root, labelHeight);
        //print(root);
    }

    // Builds the suffix tree
    void buildTree(ContainerType &Seq1, MatchFnTy matchFn)
    {
        // seq = Seq1;
        for (int i = 0; i < Seq1.size(); i++)
        {
            seq.push_back(Seq1[i]);
        }
        Ty dollar = Blank;
        seq.push_back(dollar);

        seqSize = seq.size();
        size1 = seqSize;

        int *a = new int;
        std::shared_ptr<int> newPtr(a);
        rootEnd = newPtr;
        *rootEnd = -1;

        int *b = new int;
        std::shared_ptr<int> newPtrb(b);
        leafEnd = newPtrb;
        *leafEnd = -1;

        root = new SuffixNode(-1, rootEnd, nullptr, seqSize);
        root->suffixLink = root;

        match = matchFn;

        activeNode = root;

        general = false;

        // Begin construction
        for (int i = 0; i < seqSize; i++)
        {
            extendSuffixTree(i);
        }

        int labelHeight = 0;
        setSuffixIndexByDFS(root, labelHeight);

        //print(root);
    }

    void print(SuffixNode *node)
    {

        for (int i = 0; i < alphabet.size(); i++)
        {
            if (node->children[i] != nullptr)
            {
                
                typename std::vector<Ty>::const_iterator first = seq.begin() + node->children[i]->start;
                typename std::vector<Ty>::const_iterator last = seq.begin() + *(node->children[i]->end)+1;
                std::vector<Ty> sub(first, last);

                for (int i = 0; i < sub.size(); i++)
                {
                    if (sub[i] == NULL)
                    {
                        std::cout << "$";
                    }
                    else 
                    {
                        std::cout << sub[i];
                    }
                }

                std::cout << " " << node->children[i]->suffixIndex;

                std::cout << std::endl;

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

            /*for (int i = n->start; i <= *(n->end); i++)
            {
                //if (myMatch(seq[i], nullptr, i))
                if (myMatch(seq[i], Blank, i))
                {
                    //newfile.open(newPath, std::ios_base::app);
                    //newfile << "This nonsense\n";
                    //newfile.close();
                    //n->end = std::shared_ptr<int>(new int);
                    //*(n->end) = i;
                    //newfile.open(newPath, std::ios_base::app);
                    //newfile << "Not this nonsense\n";
                    //newfile.close();
                }
            }*/
            n->suffixIndex = seqSize - labelHeight;
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
        // If activeLength is greater than current edge length,
        // set next internal node as activeNode and adjust activeEdge
        // and activeLength accordingly
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
        if (!leafEnd)
        {
            int *a = new int;
            std::shared_ptr<int> newPtr(a);
            leafEnd = newPtr;
        }
        *leafEnd = pos;
        remainingSuffixCount++;
        lastNewNode = nullptr;

        // Add all suffixes (yet to be added) one by one
        while (remainingSuffixCount > 0)
        {

            if (activeLength == 0)
            {
                activeEdge = pos;
            }

            // Find the active edge
            auto result1 = find(std::begin(alphabet), std::end(alphabet), seq[pos - activeLength], pos);
            if (result1 == std::end(alphabet))
            {
                alphabet.push_back(seq[pos - activeLength]);
                activeEdge = alphabet.size() - 1;
            }
            else
            {
                // Get index of element from iterator
                activeEdge = std::distance(alphabet.begin(), result1);
            }

            // There is no outgoing edge starting with activeEdge from activeNode
            if (activeNode->children[activeEdge] == nullptr)
            {
                std::shared_ptr<int> sharedLeaf = leafEnd;
                activeNode->children[activeEdge] = new SuffixNode(pos, sharedLeaf, root, seqSize);

                if (lastNewNode != nullptr)
                {
                    lastNewNode->suffixLink = activeNode;
                    lastNewNode = nullptr;
                }

            }
            // There is an outgoing edge starting with activeEdge from activeNode
            else
            {
                // Get the next node at the end of edge starting with activeEdge
                SuffixNode *next = activeNode->children[activeEdge];

                if (walkDown(next))
                {
                    // Start from next node (the new activeNode)
                    continue;
                }

                // Extension rule 3
                // Current character being processed is already on the edge
                if (myMatch(seq[next->start + activeLength], seq[pos], pos))
                {
                    // If a newly created node waiting for it's
                    // suffix link to be set, then set suffix link
                    // to current activeLength
                    if (lastNewNode != nullptr && activeNode != root)
                    {
                        lastNewNode->suffixLink = activeNode;
                        lastNewNode = nullptr;
                    }

                    activeLength++;

                    // Stop all processing in this phase
                    // Move onto next phase
                    break;
                }

                // Extension rule 2
                // activePoint is in the middle of the edge being traversed
                // and current character being processed is not on the edge
                // So add a new node and a new leaf edge out of that node
                auto result1 = find(std::begin(alphabet), std::end(alphabet), seq[pos], pos);
                int leafEdge = activeEdge;
                if (result1 == std::end(alphabet))
                {
                    alphabet.push_back(seq[pos]);
                    leafEdge = alphabet.size() - 1;
                }
                else
                {
                    // Get index of element from iterator
                    leafEdge = std::distance(alphabet.begin(), result1);
                }

                std::shared_ptr<int> newTing(new int);
                splitEnd = newTing;
                *splitEnd = next->start + activeLength - 1;

                // New internal node from active node
                SuffixNode *split = new SuffixNode(next->start, splitEnd, root, seqSize);
                split->secondStart = pos;
                activeNode->children[activeEdge] = split;


                // Create leaves of new internal node

                // New leaf coming out of new internal node
                std::shared_ptr<int> sharedLeaf = leafEnd;
                split->children[leafEdge] = new SuffixNode(pos, sharedLeaf, root, seqSize);


                // New edge coming out of internal node
                next->start += activeLength;
                auto result2 = find(std::begin(alphabet), std::end(alphabet), seq[next->start], next->start);
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

                split->children[continuingEdge] = new SuffixNode(pos, sharedLeaf, root, seqSize);
                delete split->children[continuingEdge];
                split->children[continuingEdge] = next;

                // We have a new internal node
                // If there is any internal node created in last extensions
                // which are still waiting for it's suffix link reset then reset
                if (lastNewNode != nullptr)
                {
                    // suffixLink of lastNewNode points to current newly created node
                    lastNewNode->suffixLink = split;
                }

                lastNewNode = split;
            }

            // One suffix got added in tree
            // Decrement the count of suffixes yet to be added
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

    void deleteTree()
    {
        freeSuffixTreeByPostOrder(root);

        alphabet.clear();
        seq.clear();
    }

    void freeSuffixTreeByPostOrder(SuffixNode* n)
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

        delete n;
    }

    //<----- Funcs used for LCS ----->
    AlignedSequence<Ty, Blank> getLCS(ContainerType &Seq1, ContainerType &Seq2, MatchFnTy matchFn)
    {

        buildTree(Seq1, Seq2, matchFn);
        ContainerType ret;
        int maxHeight = 0;
        int substringStartIndex = 0;
        int secondIndex = 0;
        int suffixIndexX = 0;
        int suffixIndexY = 0;


        traverse(root, 0, &maxHeight, &substringStartIndex, &secondIndex, &suffixIndexX, &suffixIndexY);

        AlignedSequence<Ty, Blank> Result;

        for (int i = 0; i < maxHeight; i++)
        {
            Result.Data.push_back(
                typename BaseType::EntryType(Seq1[substringStartIndex + i], Seq2[secondIndex + i], true)
            );
        }

        return Result;
    }

    int traverse(SuffixNode *n, int labelHeight, int *maxHeight, int *substringStartIndex, int *secondIndex, int *suffixIndexX, int *suffixIndexY)
    {
        if (n == nullptr)
        {
            return 0;
        }

        int i = 0;
        int ret = -1;

        if (n->suffixIndex < 0) // If it is an internal node
        {
            for (i = 0; i < alphabet.size(); i++)
            {
                if (n->children[i] != nullptr)
                {
                    ret = traverse(n->children[i], labelHeight + edgeLength(n->children[i]), maxHeight, substringStartIndex, secondIndex, suffixIndexX, suffixIndexY);

                    if (n->suffixIndex == -1)
                    {
                        n->suffixIndex = ret;
                    }
                    else if ((n->suffixIndex == -2 && ret == -3) || (n->suffixIndex == -3 && ret == -2) || n->suffixIndex == -4)
                    {
                        n->suffixIndex = -4; // Mark node as XY
                        if (*maxHeight < labelHeight)
                        {
                            *maxHeight = labelHeight;
                            *substringStartIndex = *(n->end) - labelHeight + 1;
                            *secondIndex = *suffixIndexY - size1;
                        }
                    }
                }
            }
        }
        else if (n->suffixIndex > -1 && n->suffixIndex < size1)
        {
            *suffixIndexX = n->suffixIndex;
            return -2; // Mark node as X
        }
        else if (n->suffixIndex >= size1)
        {
            *suffixIndexY = n->suffixIndex;
            return -3; // Mark node as Y
        }
        return n->suffixIndex;
    }

    //<----- Funcs used for pattern matching ----->
    template <typename ArrayType>
    int traverseEdge(ArrayType pattern, int idx, int start, int end)
    {
        int k = 0;
        // Traverse the edge with character by character matching
        for (k = start; k <= end && idx < pattern.size(); k++, idx++)
        {

            if (!myMatch(seq[k], pattern[idx], 0))
                return -1; // no match
        }

        if (idx == pattern.size())
        {
            return 1; // match
        }
        return 0; // more characters yet to match
    }

    int doTraversalToCountLeaf(SuffixNode *n, std::vector<candidateWord> &candidates)
    {
        if (n == nullptr)
            return 0;
        if (n->suffixIndex > -1)
        {
            candidateWord cWord;
            cWord.indexSeq1 = n->suffixIndex;
            candidates.push_back(cWord);
            return 1;
        }
        int count = 0;
        int i = 0;
        for (i = 0; i < alphabet.size(); i++)
        {
            if (n->children[i] != nullptr)
            {
                count += doTraversalToCountLeaf(n->children[i], candidates);
            }
        }
        return count;
    }

    int countLeaf(SuffixNode *n, std::vector<candidateWord> &candidates)
    {
        if (n == nullptr)
            return 0;
        return doTraversalToCountLeaf(n, candidates);
    }

    template <typename ArrayType>
    void doTraversal(SuffixNode *n, ArrayType pattern, int idx, std::vector<candidateWord> &candidates)
    {
        if (n == nullptr)
        {
            return; // no built tree so return no match
        }
        int res = -1;
        // If node n is not root node, then traverse edge
        // from node n's parent to node n.
        if (n->start != -1)
        {
            res = traverseEdge(pattern, idx, n->start, *(n->end));
            if (res == -1) // no match
                return;
            if (res == 1) // match
            {
                if (n->suffixIndex > -1)
                {
                    candidateWord cWord;
                    cWord.indexSeq1 = n->suffixIndex;
                    candidates.push_back(cWord);
                }
                else
                {
                    countLeaf(n, candidates);
                }

                return;
            }
        }

        // Get the character index to search
        idx = idx + edgeLength(n);

        // If there is an edge from node n going out
        // with current character str[idx], traverse that edge
        auto result1 = find(std::begin(alphabet), std::end(alphabet), pattern[idx], *(n->end));
        int index = idx;
        if (result1 != std::end(alphabet))
        {
            index = std::distance(alphabet.begin(), result1);
        }
        if (n->children[index] != nullptr)
            doTraversal(n->children[index], pattern, idx, candidates);
        else
            return; // no match
    }

    template <typename ArrayType>
    std::vector<candidateWord> getCandidates(ArrayType pattern)
    {
        std::vector<candidateWord> candidates;
        doTraversal(root, pattern, 0, candidates);
        return candidates;
    }

    //<----- Funcs for getting MUMs ----->
    std::vector<MUM> getMUMs(ContainerType &Seq1, ContainerType &Seq2, MatchFnTy matchFn)
    {

        // Build tree
        buildTree(Seq1, Seq2, matchFn);

        // Mark nodes
        std::vector<SuffixNode *> markedNodes;
        std::vector<int> mumStarts;
        std::vector<bool> isSuffix;
        unsigned int length = 0;
        markNodes(root, markedNodes, root, mumStarts, length);

        // Unmark nodes with suffix links to marked nodes
        for (int i = 0; i < markedNodes.size(); i++)
        {
            isSuffix.push_back(false);
            for (int j = 0; j < markedNodes.size(); j++)
            {
                // Skip itself
                if (i != j)
                {
                    // Remove nodes with suffix links to marked nodes
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
                MUM mum = MUM(mumStarts[i], markedNodes[i]->secondStart, (*(markedNodes[i]->end) - mumStarts[i]) + 1);
                ret.push_back(mum);
            }
        }

        // If no MUM available then simply return the longest common substring
        if (ret.size() == 0)
        {
            int maxHeight = 0;
            int substringStartIndex = 0;
            int secondIndex = 0;
            int suffixIndexX = 0;
            int suffixIndexY = 0;
            traverse(root, 0, &maxHeight, &substringStartIndex, &secondIndex, &suffixIndexX, &suffixIndexY);
            MUM mum = MUM(substringStartIndex, secondIndex, maxHeight);
            ret.push_back(mum);
        }

        // If any are bad then remove
        for (int i = 0; i < ret.size(); i++)
        {
            if (ret[i].indexSeq1 < 0 || ret[i].indexSeq2 < 0 || ret[i].indexSeq1 + ret[i].length > Seq1.size() || ret[i].indexSeq2 + ret[i].length > Seq2.size())
            {
                std::cout << "BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD" << std::endl;

                /*newfile << "idx seq1: " << ret[i].indexSeq1 << "\n";
        newfile << "idx seq2: " << ret[i].indexSeq2 << "\n";
        newfile << "mum length: " << ret[i].length << "\n";
        newfile << "seq1 size: " << Seq1.size() << "\n";
        newfile << "seq2 size: " << Seq2.size() << "\n";*/
                ret.erase(ret.begin() + i);
                // newfile << "REMOVED\n";
                i--;
            }
        }

        return ret;
    }

    // Recursively loop over all nodes and mark them if they have two leaf
    // children
    void markNodes(SuffixNode *n, std::vector<SuffixNode *> &markedNodes, SuffixNode *startNode, std::vector<int> &mumStarts, unsigned int length)
    {
        if (n == nullptr)
        {
            return;
        }

        int noChildren = 0;
        bool isSeq1 = false;
        bool isSeq2 = false;
        int secStart = 0;

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
                        secStart = n->children[i]->suffixIndex - size1;
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

        // A POTENTIAL MUM
        if (noChildren == 2 && isSeq1 && isSeq2)
        {
            n->secondStart = secStart;
            markedNodes.push_back(n);
            mumStarts.push_back(n->start + edgeLength(n) - length);
        }
    }

    bool myMatch(Ty a, Ty b, int pos)
    {
        if (a == Blank || b == Blank || a == Ty(0) || b == Ty(0) || a == NULL || b == NULL)
        {
            if (a == Blank && b == Blank)
            {
                if (pos == size1 - 1)
                {
                    return true;
                }
            }
            return false;
        }
        else
        {
            return match(a, b);
        }
    }

    template <class InputIt, class T>
    constexpr InputIt find(InputIt first, InputIt last, const T &value, int pos)
    {
        if (value == Blank)
        {
            if (general)
            {
                if (!seenHash)
                {
                    seenHash = true;
                    return last;
                }
                else
                {
                    if (pos < seq.size() - 1)
                    {
                        for (; first != last; ++first)
                        {
                            if (myMatch(*first, value, pos))
                            {
                                return first;
                            }
                        }
                        //This is just the index of the sequence
                        return first + size1 - 1;
                    }
                    else
                    {
                        if (!seenDollar)
                        {
                            seenDollar = true;
                            return last;
                        }
                        else
                        {
                            return last - 1;
                        }
                    }
                }
            }
            else
            {
                if (!seenDollar)
                {
                    seenDollar = true;
                    return last;
                }
                else
                {
                    return last - 1;
                }
            }
        }

        for (; first != last; ++first)
        {
            if (myMatch(*first, value, pos))
            {
                return first;
            }
        }

        return last;
    }
};
